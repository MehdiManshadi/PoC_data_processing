"""Download PRIDE proteinGroups tables and report coverage without reference lists.

Set INPUT_FILE to the PRIDE finder's Excel or CSV report (column: accession).
The existing downloader discovers all proteinGroups tables for each accession.
Gene extraction and UniProt fallback are unchanged. A row passes only when its
LFQ values are non-zero and non-missing in strictly more than 10% of all LFQ
columns. Coverage counts unique gene names in passing rows, case-insensitively.

Writes dataset_coverage_report.xlsx with dataset_coverage and table_coverage
sheets. The dataset count is the union of qualifying genes across its tables;
failed tables have blank counts. A mix of successful and failed table analyses
is marked partial.
Both sheets retain input organisms, study year (from publication_date), paper
title/link and tissue information, plus the LFQ column names for each record.

Dependencies: pandas, openpyxl, requests, remotezip, python-dotenv.
"""

import gzip
import hashlib
import re
import shutil
import tempfile
import time
from pathlib import Path, PurePosixPath
from zipfile import ZipFile
from urllib.parse import unquote, urlsplit

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

load_dotenv()
DELETE_DOWNLOADED_FILES = False
INPUT_FILE = Path("/Users/mehman/Projects/PoC_data_processing/zebrafish_test.csv")
INPUT_SHEET = 0  # First Excel sheet; change to a sheet name or zero-based index.
OUTPUT_DIR = Path("/Users/mehman/Projects/PoC_data_processing/pp")

SAVE_PROCESSED_TABLES = True
PRESENCE_THRESHOLD = 0.10
MAX_FULL_ZIP_DOWNLOAD_GB = 20
REQUEST_TIMEOUT = 300

PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
UNIPROT_API = "https://rest.uniprot.org"

METADATA_COLUMNS = [
    "organisms", "study_year", "publication_url", "publication_title", "combined_tissues",
]

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

def read_input_report(path, sheet_name=0):
    path = Path(path)
    if path.suffix.casefold() == ".csv":
        report = pd.read_csv(path, dtype=str, encoding="utf-8-sig")
    else:
        report = pd.read_excel(path, sheet_name=sheet_name, dtype=str)
    columns = {str(column).strip().casefold(): column for column in report.columns}
    if "accession" not in columns:
        raise ValueError("Input must be a PRIDE finder report with an 'accession' column.")
    values = report[columns["accession"]].dropna().str.strip().str.upper()
    values = values[values.ne("")].drop_duplicates()
    invalid = values[~values.str.fullmatch(r"PXD\d+")].tolist()
    if invalid:
        raise ValueError(f"Invalid PRIDE accession(s): {', '.join(invalid)}")
    # The finder has one metadata row per accession; retain the first occurrence.
    report = report.rename(columns={value: key for key, value in columns.items()}).fillna("")
    report["accession"] = report["accession"].str.strip().str.upper()
    metadata = {}
    for _, row in report.drop_duplicates("accession").iterrows():
        accession = row["accession"]
        if not accession:
            continue
        year = re.search(r"\b(\d{4})\b", str(row.get("publication_date", "")))
        tissues = row.get("combined_tissues", "") or "; ".join(dict.fromkeys(
            value for value in (
                row.get("pride_structured_tissues", ""),
                row.get("omicsdi_structured_tissues", ""),
                row.get("sdrf_tissues", ""),
            ) if value
        ))
        metadata[accession] = {
            "organisms": row.get("organisms", ""),
            "study_year": year.group(1) if year else "",
            "publication_url": row.get("publication_url", ""),
            "publication_title": row.get("publication_title", ""),
            "combined_tissues": tissues,
        }
    return values.tolist(), metadata


def session_with_retries():
    retry = Retry(
        total=4,                    # initial request + up to 4 retries
        connect=4,
        read=4,
        status=4,
        backoff_factor=5,           # progressively longer waits
        status_forcelist=(403, 429, 500, 502, 503, 504),
        allowed_methods={"GET", "HEAD", "POST"},
        respect_retry_after_header=True,
        raise_on_status=False,
    )

    session = requests.Session()

    session.headers.update({
        "User-Agent": (
            "PoC-data-processing/1.0 "
            "(https://github.com/MehdiManshadi/PoC_data_processing)"
        )
    })

    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


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


# ============================================================
# FIND / DOWNLOAD proteinGroups.txt
# ============================================================


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


def protein_groups_output_path(output_dir, accession, address):
    """Keep original filename, but use a unique folder to prevent overwriting."""
    url, separator, member = address.partition("!/")

    if separator:
        # File inside ZIP
        filename = base(member)
    else:
        # Direct file
        filename = base(unquote(urlsplit(url).path))

    # We decompress .gz files
    if filename.lower().endswith(".gz"):
        filename = filename[:-3]

    # Unique folder for this specific source
    source_id = hashlib.sha256(address.encode("utf-8")).hexdigest()[:12]

    destination_dir = (
        Path(output_dir)
        / accession
        / "proteinGroups_files"
        / source_id
    )
    destination_dir.mkdir(parents=True, exist_ok=True)

    return destination_dir / filename


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
# COVERAGE
# ============================================================

def passing_genes(df):
    genes = set(
        df.loc[df["Presence"] > PRESENCE_THRESHOLD, "Gene names"]
        .dropna()
        .str.split(";")
        .explode()
        .str.strip()
        .str.upper()
    )
    genes.discard("")
    return genes


def main():
    accessions, metadata = read_input_report(INPUT_FILE, INPUT_SHEET)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    session = session_with_retries()
    uniprot_cache, dataframes = {}, {}
    dataset_rows, table_rows = [], []

    try:
        for number, accession in enumerate(accessions, 1):
            print(f"[{number}/{len(accessions)}] {accession}")
            dataframes[accession] = {}
            dataset_genes, rows = set(), []
            dataset_lfq_names = []
            try:
                sources = find_and_download_all_proteingroups(accession, OUTPUT_DIR, session)
            except Exception as error:
                # Keep a failed dataset in the report rather than omitting it.
                sources = [extracted_source("", None, "", error=error)]

            for source in sources:
                address = source["proteinGroups_address"]
                raw_file = source["local_file"]
                row = {
                    "accession": accession,
                    "protein_groups_url": address,
                    "total_measured_genes": pd.NA,
                    "downloaded_file": raw_file,
                    "processed_file": "",
                    "analysis_status": "failed",
                    "analysis_error": source["download_error"],
                    **metadata[accession],
                    "lfq_column_names": "",
                }
                if not source["download_error"]:
                    try:
                        header = pd.read_csv(raw_file, sep="\t", nrows=0).columns
                        lfq_names = [
                            LFQ_RE.sub("", str(c).strip()).strip()
                            for c in header
                            if LFQ_RE.match(str(c).strip())]
                        row["lfq_column_names"] = "; ".join(lfq_names)
                        dataset_lfq_names.extend(lfq_names)
                        df = create_presence_dataframe(raw_file, session, uniprot_cache)
                        genes = passing_genes(df)
                        dataframes[accession][address] = df
                        if SAVE_PROCESSED_TABLES:
                            processed_dir = OUTPUT_DIR / accession / "processed_tables"
                            processed_dir.mkdir(parents=True, exist_ok=True)
                            processed_file = processed_dir / f"{Path(raw_file).stem}_processed.txt"
                            df.to_csv(processed_file, sep="\t", index=False)
                            row["processed_file"] = str(processed_file)
                        dataset_genes.update(genes)
                        row.update({
                            "total_measured_genes": len(genes),
                            "analysis_status": "success",
                            "analysis_error": "",
                        })
                        print(f"  {Path(raw_file).name}: {len(genes)} genes passed >10%")
                    except Exception as error:
                        row["analysis_error"] = str(error)
                if row["analysis_status"] == "failed":
                    print(f"  FAILED {address}: {row['analysis_error']}")
                rows.append(row)

            table_rows.extend(rows)
            successful = sum(row["analysis_status"] == "success" for row in rows)
            status = "success" if successful == len(rows) else "partial"
            if not successful:
                status = "failed"
            dataset_rows.append({
                "accession": accession,
                "total_measured_genes": len(dataset_genes) if successful else pd.NA,
                "analysis_status": status,
                **metadata[accession],
                "lfq_column_names": "; ".join(dict.fromkeys(dataset_lfq_names)),
                "analysis_error": "; ".join(dict.fromkeys(
                    row["analysis_error"] for row in rows if row["analysis_error"]
                )),
            })
    finally:
        session.close()

    table_coverage = pd.DataFrame(table_rows, columns=[
            "accession", "protein_groups_url", "total_measured_genes",
            "downloaded_file", "processed_file", "analysis_status", "analysis_error",
            *METADATA_COLUMNS, "lfq_column_names",
        ])
    dataset_coverage = pd.DataFrame(dataset_rows, columns=[
        "accession", "total_measured_genes", "analysis_status", "analysis_error",
        *METADATA_COLUMNS, "lfq_column_names",
    ])
    
    if table_coverage.duplicated(["accession", "protein_groups_url"]).any():
        raise RuntimeError("Duplicate dataset/proteinGroups address pairs in the report.")
    for report in (dataset_coverage, table_coverage):
        report["total_measured_genes"] = report["total_measured_genes"].astype("Int64")
        report_file = OUTPUT_DIR / "dataset_coverage_report.xlsx"
    with pd.ExcelWriter(report_file, engine="openpyxl") as writer:
        table_coverage.to_excel(writer, sheet_name="table_coverage", index=False)
        dataset_coverage.to_excel(writer, sheet_name="dataset_coverage", index=False)

    print(f"Saved {report_file}")

    # Delete all downloaded/processed files after analysis
    if DELETE_DOWNLOADED_FILES:
        for accession in accessions:
            accession_dir = OUTPUT_DIR / accession
            if accession_dir.exists():
                shutil.rmtree(accession_dir)

        print("Deleted all downloaded files.")
        
    return dataframes, dataset_coverage

if __name__ == "__main__":
    dataframes, dataset_coverage = main()
