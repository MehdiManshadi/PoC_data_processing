import gzip
import hashlib
import re
import shutil
import tempfile
import time
from pathlib import Path, PurePosixPath
from zipfile import ZipFile

import pandas as pd
import requests
from remotezip import RemoteZip
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from dotenv import load_dotenv
import os


# ============================================================
# SETTINGS
# ============================================================

PRIDE_ACCESSIONS = (
    pd.read_excel(
        "/Users/mehman/Projects/PoC_data_processing/Human_lfq_proteingroups_report.xlsx",
        sheet_name=1,  # Second sheet
        usecols=[0],   # First column
        dtype=str,
    )
    .iloc[:, 0]
    .dropna()
    .str.strip()
    .loc[lambda values: values.ne("")]
    .drop_duplicates()
    .tolist()
)
'''
PRIDE_ACCESSIONS = [
    "PXD010489"
]
'''
# Load the variables from the .env file
load_dotenv()

# Read the strings from environment and wrap them in Path()
OUTPUT_DIR = Path(os.getenv("OUTPUT_DIR"))

HYPERGLYCEMIA_FEATURES_FILE = Path(os.getenv("HYPERGLYCEMIA_FEATURES_FILE"))

MITOCHONDRIAL_MYOPATHY_FEATURES_FILE = Path(os.getenv("MITOCHONDRIAL_MYOPATHY_FEATURES_FILE"))

SAVE_PROCESSED_TABLES = True
FEATURE_PRESENCE_THRESHOLD = 0.10
MAX_FULL_ZIP_DOWNLOAD_GB = 20
REQUEST_TIMEOUT = 300

PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
UNIPROT_API = "https://rest.uniprot.org"

GENE_ALIASES = {
    "gene names", "gene name", "gene", "genes",
    "gene symbol", "gene symbols", "genesymbol", "genesymbols",
}
FASTA_ALIASES = {"fasta header", "fasta headers"}
MAJORITY_ALIASES = {"majority protein id", "majority protein ids"}

LFQ_RE = re.compile(r"^LFQ\s+intensity(?:\s|$)", re.I)
UNIPROT_RE = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)"
)


# ============================================================
# BASIC HELPERS
# ============================================================

def session_with_retries():
    retry = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods={"GET", "HEAD", "POST"},
    )
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


def norm(x):
    return re.sub(r"[^a-z0-9]+", " ", str(x).casefold()).strip()


def base(path):
    return PurePosixPath(str(path).replace("\\", "/")).name


def api_files(payload):
    if isinstance(payload, list):
        return payload
    if not isinstance(payload, dict):
        return []
    for key in ("files", "content", "results"):
        if isinstance(payload.get(key), list):
            return payload[key]
    embedded = payload.get("_embedded", {})
    if isinstance(embedded, dict):
        for value in embedded.values():
            if isinstance(value, list):
                return value
    return []


def get_pride_files(accession, session):
    """Read ALL PRIDE file pages; important for datasets with >20 files."""
    files, seen = [], set()

    for page in range(1000):
        r = session.get(
            f"{PRIDE_API}/projects/{accession}/files",
            params={"page": page},
            timeout=REQUEST_TIMEOUT,
        )
        r.raise_for_status()
        batch = api_files(r.json())

        if not batch:
            break

        added = 0
        for record in batch:
            name = str(record.get("fileName") or "").strip()
            key = str(record.get("accession") or name)
            if name and key not in seen:
                seen.add(key)
                files.append(record)
                added += 1

        if added == 0:
            break

    return files


def download_url(record):
    for loc in record.get("publicFileLocations") or []:
        value = str(loc.get("value") or "").strip()
        if value.startswith(("http://", "https://")):
            return value
        if value.startswith("ftp://ftp.pride.ebi.ac.uk/"):
            return value.replace(
                "ftp://ftp.pride.ebi.ac.uk/",
                "https://ftp.pride.ebi.ac.uk/",
                1,
            )
    raise ValueError(f"No HTTP/FTP URL for {record.get('fileName')}")


def file_size(record):
    try:
        return int(record.get("fileSizeBytes") or 0)
    except (TypeError, ValueError):
        return 0


def parse_header(line):
    return line.decode("utf-8-sig", errors="replace").rstrip("\r\n").split("\t")


def gene_sources(columns):
    normalized = {norm(c): c for c in columns}

    gene_col = next(
        (normalized[a] for a in GENE_ALIASES if a in normalized),
        None,
    )

    if gene_col is None:
        excluded = {"gene ontology", "gene ontology id", "gene ontology ids", "gene description"}
        gene_col = next(
            (
                c for c in columns
                if norm(c).startswith("gene") and norm(c) not in excluded
            ),
            None,
        )

    fasta_cols = [c for c in columns if norm(c) in FASTA_ALIASES]
    majority_col = next((c for c in columns if norm(c) in MAJORITY_ALIASES), None)
    lfq_cols = [c for c in columns if LFQ_RE.match(str(c).strip())]

    return gene_col, fasta_cols, majority_col, lfq_cols


def valid_header(header):
    gene_col, fasta_cols, majority_col, lfq_cols = gene_sources(header)
    return bool((gene_col or fasta_cols or majority_col) and lfq_cols), len(lfq_cols)


# ============================================================
# FIND / DOWNLOAD proteinGroups.txt
# ============================================================

def read_direct_header(url, session, gz=False):
    with session.get(url, stream=True, timeout=REQUEST_TIMEOUT) as r:
        r.raise_for_status()
        r.raw.decode_content = True
        if gz:
            with gzip.GzipFile(fileobj=r.raw) as f:
                return parse_header(f.readline(4 * 1024 * 1024))
        return parse_header(r.raw.readline(4 * 1024 * 1024))


def read_zip_header(archive, member):
    with archive.open(member) as f:
        if member.lower().endswith(".gz"):
            with gzip.GzipFile(fileobj=f) as g:
                return parse_header(g.readline(4 * 1024 * 1024))
        return parse_header(f.readline(4 * 1024 * 1024))


def copy_stream(source, destination):
    with open(destination, "wb") as out:
        shutil.copyfileobj(source, out, length=4 * 1024 * 1024)


def extract_zip_member(archive, member, output):
    with archive.open(member) as f:
        if member.lower().endswith(".gz"):
            with gzip.GzipFile(fileobj=f) as g:
                copy_stream(g, output)
        else:
            copy_stream(f, output)


def download_whole(url, output, session):
    with session.get(url, stream=True, timeout=REQUEST_TIMEOUT) as r:
        r.raise_for_status()
        with open(output, "wb") as out:
            for chunk in r.iter_content(4 * 1024 * 1024):
                if chunk:
                    out.write(chunk)


def is_proteingroups_file(path):
    filename = base(path)
    return bool(
        re.fullmatch(
            r"proteingroups(?:[._-].*)?\.txt(?:\.gz)?",
            filename,
            flags=re.I,
        )
    )


def protein_groups_address(url, member=None):
    """Return a stable, unique address for a direct file or a ZIP member."""
    if member is None:
        return url
    return f"{url}!/{str(member).lstrip('/')}"


def protein_groups_output_path(output_dir, accession, address, archive_name, member=None):
    """Create a readable, collision-safe local name for one proteinGroups table."""
    if member is None:
        label = "direct"
    else:
        archive_label = re.sub(r"\.zip$", "", base(archive_name), flags=re.I)
        member_parent = str(PurePosixPath(member).parent)
        label = archive_label if member_parent == "." else f"{archive_label}_{member_parent}"

    label = re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("._-") or "source"
    label = label[:80]
    digest = hashlib.sha256(address.encode("utf-8")).hexdigest()[:12]

    destination_dir = Path(output_dir) / accession / "proteinGroups_files"
    destination_dir.mkdir(parents=True, exist_ok=True)
    return destination_dir / f"{accession}__{label}__{digest}.txt"


def extracted_source(address, local_file, archive_name, member=None, error=None):
    return {
        "proteinGroups_address": address,
        "local_file": str(local_file) if local_file is not None else "",
        "archive_name": archive_name if member is not None else "",
        "archive_member": member or "",
        "download_error": str(error) if error else "",
    }


def extract_all_zip_members(
    archive,
    accession,
    output_dir,
    archive_name,
    url,
    seen_addresses,
    continue_on_error=False,
):
    """Extract every proteinGroups.txt(.gz) member from one open ZIP."""
    extracted = []
    newly_seen = set()

    for info in archive.infolist():
        if info.is_dir() or not is_proteingroups_file(info.filename):
            continue

        member = info.filename
        address = protein_groups_address(url, member)
        if address in seen_addresses or address in newly_seen:
            continue

        output_file = protein_groups_output_path(
            output_dir,
            accession,
            address,
            archive_name,
            member,
        )

        try:
            extract_zip_member(archive, member, output_file)
            extracted.append(
                extracted_source(address, output_file, archive_name, member)
            )
            newly_seen.add(address)
        except Exception as e:
            if not continue_on_error:
                raise
            extracted.append(
                extracted_source(address, None, archive_name, member, error=e)
            )
            newly_seen.add(address)

    # Update the shared set only after the archive pass finishes. If remote ZIP
    # extraction fails halfway, the full-download fallback can retry all members.
    seen_addresses.update(newly_seen)
    return extracted


def find_and_download_all_proteingroups(accession, output_dir, session):
    """
    Download every direct or ZIP-contained proteinGroups.txt(.gz) in a dataset.

    The returned list contains one item per unique source address. Files are not
    filtered here by LFQ or gene columns; each discovered table is passed to the
    unchanged downstream analysis, which reports any table it cannot analyse.
    """
    records = get_pride_files(accession, session)
    if not records:
        raise FileNotFoundError("PRIDE returned no files.")

    sources = []
    seen_addresses = set()

    # 1) Direct proteinGroups.txt / proteinGroups.txt.gz
    for record in records:
        record_name = str(record.get("fileName") or "")
        if not is_proteingroups_file(record_name):
            continue

        address = None
        try:
            url = download_url(record)
            address = protein_groups_address(url)
            if address in seen_addresses:
                continue

            output_file = protein_groups_output_path(
                output_dir,
                accession,
                address,
                record_name,
            )

            with session.get(url, stream=True, timeout=REQUEST_TIMEOUT) as r:
                r.raise_for_status()
                r.raw.decode_content = True
                if base(record_name).casefold().endswith(".gz"):
                    with gzip.GzipFile(fileobj=r.raw) as g:
                        copy_stream(g, output_file)
                else:
                    copy_stream(r.raw, output_file)

            sources.append(extracted_source(address, output_file, record_name))
            seen_addresses.add(address)
        except Exception as e:
            if address is not None and address not in seen_addresses:
                sources.append(extracted_source(address, None, record_name, error=e))
                seen_addresses.add(address)

    # 2) proteinGroups.txt inside ZIP archives
    for record in records:
        archive_name = str(record.get("fileName") or "")
        if not base(archive_name).casefold().endswith(".zip"):
            continue

        try:
            url = download_url(record)
        except Exception:
            continue

        try:
            with RemoteZip(
                url,
                session=session,
                timeout=REQUEST_TIMEOUT,
                initial_buffer_size=1024 * 1024,
                support_suffix_range=False,
            ) as z:
                archive_sources = extract_all_zip_members(
                    z,
                    accession,
                    output_dir,
                    archive_name,
                    url,
                    seen_addresses,
                )
            sources.extend(archive_sources)
            continue
        except Exception as remote_error:
            # Full ZIP fallback only when reasonably sized.
            size = file_size(record)
            if size and size > MAX_FULL_ZIP_DOWNLOAD_GB * 1024**3:
                print(
                    f"  Could not inspect {archive_name} remotely and its size "
                    f"exceeds MAX_FULL_ZIP_DOWNLOAD_GB: {remote_error}"
                )
                continue

            with tempfile.TemporaryDirectory() as td:
                local_zip = Path(td) / "archive.zip"
                try:
                    download_whole(url, local_zip, session)
                    with ZipFile(local_zip) as z:
                        sources.extend(
                            extract_all_zip_members(
                                z,
                                accession,
                                output_dir,
                                archive_name,
                                url,
                                seen_addresses,
                                continue_on_error=True,
                            )
                        )
                except Exception as e:
                    print(f"  Could not inspect ZIP {archive_name}: {e}")

    if not sources:
        raise FileNotFoundError("No proteinGroups.txt or proteinGroups.txt.gz found.")

    return sources


# ============================================================
# GENE EXTRACTION + UNIPROT FALLBACK
# ============================================================

def fasta_genes(value):
    if pd.isna(value):
        return []
    return re.findall(r"\bGN=([^\s;]+)", str(value), flags=re.I)


def uniprot_ids(value):
    if pd.isna(value):
        return []

    result = []
    for x in str(value).split(";"):
        x = x.strip()

        # Some Majority protein IDs cells contain empty entries.
        # Skip them before calling split()[0].
        if not x:
            continue

        if "|" in x:
            parts = x.split("|")
            if len(parts) >= 2:
                x = parts[1].strip()

        if not x:
            continue

        x = re.sub(r"-\d+$", "", x.split()[0]).upper()

        if UNIPROT_RE.fullmatch(x) and x not in result:
            result.append(x)

    return result


def map_uniprot(accessions, session, cache):
    todo = [x for x in dict.fromkeys(accessions) if x not in cache]
    if not todo:
        return cache

    for start in range(0, len(todo), 10000):
        batch = todo[start:start + 10000]

        r = session.post(
            f"{UNIPROT_API}/idmapping/run",
            data={"from": "UniProtKB_AC-ID", "to": "Gene_Name", "ids": ",".join(batch)},
            timeout=REQUEST_TIMEOUT,
        )
        r.raise_for_status()
        job = r.json()["jobId"]

        while True:
            s = session.get(
                f"{UNIPROT_API}/idmapping/status/{job}",
                timeout=REQUEST_TIMEOUT,
            )
            s.raise_for_status()
            status = s.json().get("jobStatus")
            if status in {"NEW", "RUNNING"}:
                time.sleep(3)
                continue
            if status == "FAILED":
                raise RuntimeError("UniProt mapping failed.")
            break

        r = session.get(
            f"{UNIPROT_API}/idmapping/stream/{job}",
            params={"format": "json"},
            timeout=REQUEST_TIMEOUT,
        )
        r.raise_for_status()

        for acc in batch:
            cache.setdefault(acc, [])

        for item in r.json().get("results", []):
            acc = str(item.get("from", "")).upper()
            mapped = item.get("to")

            if isinstance(mapped, dict):
                gene = mapped.get("value") or mapped.get("geneName") or mapped.get("id")
            else:
                gene = mapped

            if gene:
                gene = str(gene).strip()
                if gene and gene not in cache[acc]:
                    cache[acc].append(gene)

    return cache


def create_presence_dataframe(file_path, session, uniprot_cache):
    header = list(pd.read_csv(file_path, sep="\t", nrows=0).columns)
    gene_col, fasta_cols, majority_col, lfq_cols = gene_sources(header)

    if not lfq_cols:
        raise ValueError("No LFQ intensity columns found.")
    if not (gene_col or fasta_cols or majority_col):
        raise ValueError("No usable gene source found.")

    usecols = list(dict.fromkeys(
        ([gene_col] if gene_col else [])
        + fasta_cols
        + ([majority_col] if majority_col else [])
        + lfq_cols
    ))
    df = pd.read_csv(file_path, sep="\t", usecols=usecols, low_memory=False)

    genes_out, source_out = [], []
    unresolved_rows = {}
    all_uniprot = []

    for i, row in df.iterrows():
        genes, seen, sources = [], set(), []

        if gene_col and pd.notna(row[gene_col]):
            for g in str(row[gene_col]).split(";"):
                g = g.strip()
                if g and g.upper() not in seen:
                    seen.add(g.upper())
                    genes.append(g)
            if genes:
                sources.append(f"Gene column ({gene_col})")

        fasta_added = False
        for col in fasta_cols:
            for g in fasta_genes(row[col]):
                g = g.strip()
                if g and g.upper() not in seen:
                    seen.add(g.upper())
                    genes.append(g)
                    fasta_added = True
        if fasta_added:
            sources.append("FASTA GN=")

        if genes:
            genes_out.append(";".join(genes))
            source_out.append(" + ".join(sources))
        else:
            genes_out.append(pd.NA)
            source_out.append("Unresolved")
            if majority_col:
                ids = uniprot_ids(row[majority_col])
                if ids:
                    unresolved_rows[i] = ids
                    all_uniprot.extend(ids)

    df["Gene names"] = genes_out
    df["Gene matching source"] = source_out

    if unresolved_rows:
        try:
            map_uniprot(all_uniprot, session, uniprot_cache)
            for i, ids in unresolved_rows.items():
                genes = []
                for acc in ids:
                    for g in uniprot_cache.get(acc, []):
                        if g not in genes:
                            genes.append(g)
                if genes:
                    df.at[i, "Gene names"] = ";".join(genes)
                    df.at[i, "Gene matching source"] = "Majority protein IDs -> UniProt"
        except Exception as e:
            print(f"  UniProt fallback failed: {e}")

    lfq = df[lfq_cols].apply(pd.to_numeric, errors="coerce")
    present = lfq.notna() & lfq.ne(0)

    df["Presence_count"] = present.sum(axis=1)
    df["Total_LFQ_columns"] = len(lfq_cols)
    df["Presence"] = df["Presence_count"] / len(lfq_cols)
    df["Presence_percent"] = df["Presence"] * 100

    return df[
        ["Gene names", "Gene matching source"]
        + lfq_cols
        + ["Presence_count", "Total_LFQ_columns", "Presence", "Presence_percent"]
    ]


# ============================================================
# FEATURE MISSINGNESS
# ============================================================

def read_feature_genes(path):
    genes = []
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        for gene in re.split(r"[;,\t]+", line):
            gene = gene.strip()
            if gene and gene.casefold() not in {
                "gene", "genes", "gene name", "gene names",
                "feature", "features", "gene symbol", "gene symbols",
            }:
                genes.append(gene.upper())
    return list(dict.fromkeys(genes))


def feature_missingness(df, feature_genes):
    # Preserves the original script's matching logic:
    # any gene appearing in "Gene names" is considered found, provided that
    # its row is non-zero in at least 10% of all LFQ columns in this table.
    present = set(
        df.loc[df["Presence"] > FEATURE_PRESENCE_THRESHOLD, "Gene names"]
        .dropna()
        .str.split(";")
        .explode()
        .str.strip()
        .str.upper()
    )

    target = set(g.upper() for g in feature_genes)
    found = target & present
    missing = target - found

    return {
        "total_genes": len(target),
        "found_genes": len(found),
        "missing_genes": len(missing),
        "missingness_percent": 100 * len(missing) / len(target),
        "found_gene_names": ";".join(sorted(found)),
        "missing_gene_names": ";".join(sorted(missing)),
    }


# ============================================================
# RUN
# ============================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    hyper = read_feature_genes(HYPERGLYCEMIA_FEATURES_FILE)
    mito = read_feature_genes(MITOCHONDRIAL_MYOPATHY_FEATURES_FILE)

    session = session_with_retries()
    uniprot_cache = {}
    # Nested structure: dataframes[accession][proteinGroups_address] = dataframe
    # This preserves dataset-level access while allowing multiple tables per dataset.
    dataframes = {}
    results = []

    try:
        for n, accession in enumerate(dict.fromkeys(PRIDE_ACCESSIONS), 1):
            accession = accession.strip().upper()
            print(f"[{n}/{len(PRIDE_ACCESSIONS)}] {accession}")

            try:
                sources = find_and_download_all_proteingroups(
                    accession,
                    OUTPUT_DIR,
                    session,
                )
                dataframes.setdefault(accession, {})
                print(f"  Found {len(sources)} unique proteinGroups table(s).")

                for source_number, source in enumerate(sources, 1):
                    address = source["proteinGroups_address"]
                    raw_file = source["local_file"]

                    result = {
                        # Keep these as the first two report columns. Together
                        # they uniquely identify every analysed table.
                        "PRIDE_accession": accession,
                        "proteinGroups_address": address,
                        "downloaded_file": raw_file,
                        "source_archive": source["archive_name"],
                        "source_member": source["archive_member"],
                        "processed_file": "",
                        "analysis_status": "failed",
                        "analysis_error": source["download_error"],
                    }

                    if source["download_error"]:
                        results.append(result)
                        print(
                            f"  [{source_number}/{len(sources)}] FAILED to download "
                            f"{address}: {source['download_error']}"
                        )
                        continue

                    try:
                        df = create_presence_dataframe(
                            raw_file,
                            session,
                            uniprot_cache,
                        )
                        dataframes[accession][address] = df

                        if SAVE_PROCESSED_TABLES:
                            processed_dir = OUTPUT_DIR / accession / "processed_tables"
                            processed_dir.mkdir(parents=True, exist_ok=True)
                            processed_file = (
                                processed_dir
                                / f"{Path(raw_file).stem}_processed.txt"
                            )
                            df.to_csv(processed_file, sep="\t", index=False)
                            result["processed_file"] = str(processed_file)

                        h = feature_missingness(df, hyper)
                        m = feature_missingness(df, mito)

                        result.update({
                            "analysis_status": "success",
                            "analysis_error": "",

                            "hyperglycemia_total_genes": h["total_genes"],
                            "hyperglycemia_found_genes": h["found_genes"],
                            "hyperglycemia_missing_genes": h["missing_genes"],
                            "hyperglycemia_missingness_percent": h["missingness_percent"],
                            "hyperglycemia_found_gene_names": h["found_gene_names"],
                            "hyperglycemia_missing_gene_names": h["missing_gene_names"],

                            "mitochondrial_myopathy_total_genes": m["total_genes"],
                            "mitochondrial_myopathy_found_genes": m["found_genes"],
                            "mitochondrial_myopathy_missing_genes": m["missing_genes"],
                            "mitochondrial_myopathy_missingness_percent": m["missingness_percent"],
                            "mitochondrial_myopathy_found_gene_names": m["found_gene_names"],
                            "mitochondrial_myopathy_missing_gene_names": m["missing_gene_names"],
                        })
                        results.append(result)

                        print(
                            f"  [{source_number}/{len(sources)}] "
                            f"Hyperglycemia: {h['missingness_percent']:.2f}% | "
                            f"Mitochondrial: {m['missingness_percent']:.2f}%"
                        )

                    except Exception as e:
                        result["analysis_error"] = str(e)
                        results.append(result)
                        print(
                            f"  [{source_number}/{len(sources)}] FAILED analysis "
                            f"for {address}: {e}"
                        )

            except Exception as e:
                print(f"  FAILED: {e}")

    finally:
        session.close()

    report_columns = [
        "PRIDE_accession",
        "proteinGroups_address",
        "downloaded_file",
        "source_archive",
        "source_member",
        "processed_file",
        "analysis_status",
        "analysis_error",
        "hyperglycemia_total_genes",
        "hyperglycemia_found_genes",
        "hyperglycemia_missing_genes",
        "hyperglycemia_missingness_percent",
        "hyperglycemia_found_gene_names",
        "hyperglycemia_missing_gene_names",
        "mitochondrial_myopathy_total_genes",
        "mitochondrial_myopathy_found_genes",
        "mitochondrial_myopathy_missing_genes",
        "mitochondrial_myopathy_missingness_percent",
        "mitochondrial_myopathy_found_gene_names",
        "mitochondrial_myopathy_missing_gene_names",
    ]
    dataset_missingness_scores = pd.DataFrame(results).reindex(columns=report_columns)

    # Discovery already de-duplicates addresses. This check protects the report
    # invariant if PRIDE ever returns the same source through multiple records.
    duplicate_rows = dataset_missingness_scores.duplicated(
        subset=["PRIDE_accession", "proteinGroups_address"],
        keep=False,
    )
    if duplicate_rows.any():
        raise RuntimeError(
            "Duplicate dataset/proteinGroups address pairs found in the report."
        )

    dataset_missingness_scores.to_csv(
        OUTPUT_DIR / "dataset_missingness_scores.csv",
        index=False,
    )

    return dataframes, dataset_missingness_scores


if __name__ == "__main__":
    dataframes, dataset_missingness_scores = main()
