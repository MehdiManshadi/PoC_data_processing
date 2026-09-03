"""Reusable PRIDE LFQ download, gene matching, and group aggregation functions.

The public notebook keeps only inputs and compact reports. All implementation
details live here.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import tempfile
import time
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping, Sequence
from zipfile import ZipFile

import pandas as pd
import requests
from remotezip import RemoteZip
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
UNIPROT_API = "https://rest.uniprot.org"
DEFAULT_TIMEOUT = 300

GENE_ALIASES = {
    "gene names",
    "gene name",
    "gene",
    "genes",
    "gene symbol",
    "gene symbols",
    "genesymbol",
    "genesymbols",
}
FASTA_ALIASES = {"fasta header", "fasta headers"}
MAJORITY_ALIASES = {"majority protein id", "majority protein ids"}
IGNORED_FEATURE_HEADERS = {
    "gene",
    "genes",
    "gene name",
    "gene names",
    "feature",
    "features",
    "gene symbol",
    "gene symbols",
}

LFQ_RE = re.compile(r"^LFQ\s+intensity(?:\s|$)", re.IGNORECASE)
UNIPROT_RE = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)"
)
SUPPORTED_AGGREGATIONS = {"mean", "median", "max", "min", "sum"}


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------


def session_with_retries() -> requests.Session:
    retry = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods={"GET", "HEAD", "POST"},
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def normalized_text(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).casefold()).strip()


def base_name(path: object) -> str:
    return PurePosixPath(str(path).replace("\\", "/")).name


def _api_files(payload: object) -> list[dict]:
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


def get_pride_files(
    accession: str,
    session: requests.Session,
    timeout: int = DEFAULT_TIMEOUT,
) -> list[dict]:
    """Read all PRIDE file pages, including accessions with more than 20 files."""
    files: list[dict] = []
    seen: set[str] = set()
    for page in range(1000):
        response = session.get(
            f"{PRIDE_API}/projects/{accession}/files",
            params={"page": page},
            timeout=timeout,
        )
        response.raise_for_status()
        batch = _api_files(response.json())
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


def _download_url(record: Mapping[str, object]) -> str:
    for location in record.get("publicFileLocations") or []:
        value = str(location.get("value") or "").strip()
        if value.startswith(("http://", "https://")):
            return value
        if value.startswith("ftp://ftp.pride.ebi.ac.uk/"):
            return value.replace(
                "ftp://ftp.pride.ebi.ac.uk/",
                "https://ftp.pride.ebi.ac.uk/",
                1,
            )
    raise ValueError(f"No HTTP/FTP URL for {record.get('fileName')}")


def _file_size(record: Mapping[str, object]) -> int:
    try:
        return int(record.get("fileSizeBytes") or 0)
    except (TypeError, ValueError):
        return 0


def _parse_header(line: bytes) -> list[str]:
    return line.decode("utf-8-sig", errors="replace").rstrip("\r\n").split("\t")


def gene_sources(
    columns: Sequence[object],
) -> tuple[object | None, list[object], object | None, list[object]]:
    """Find gene, FASTA, UniProt-fallback, and LFQ columns."""
    normalized = {normalized_text(column): column for column in columns}
    gene_column = next(
        (normalized[alias] for alias in GENE_ALIASES if alias in normalized),
        None,
    )
    if gene_column is None:
        excluded = {
            "gene ontology",
            "gene ontology id",
            "gene ontology ids",
            "gene description",
        }
        gene_column = next(
            (
                column
                for column in columns
                if normalized_text(column).startswith("gene")
                and normalized_text(column) not in excluded
            ),
            None,
        )

    fasta_columns = [
        column for column in columns if normalized_text(column) in FASTA_ALIASES
    ]
    majority_column = next(
        (
            column
            for column in columns
            if normalized_text(column) in MAJORITY_ALIASES
        ),
        None,
    )
    lfq_columns = [
        column for column in columns if LFQ_RE.match(str(column).strip())
    ]
    return gene_column, fasta_columns, majority_column, lfq_columns


def _is_proteingroups_file(path: object) -> bool:
    return base_name(path).casefold() in {
        "proteingroups.txt",
        "proteingroups.txt.gz",
    }


def _read_direct_header(
    url: str,
    session: requests.Session,
    timeout: int,
    gzipped: bool,
) -> list[str]:
    with session.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        response.raw.decode_content = True
        if gzipped:
            with gzip.GzipFile(fileobj=response.raw) as handle:
                return _parse_header(handle.readline(4 * 1024 * 1024))
        return _parse_header(response.raw.readline(4 * 1024 * 1024))


def _read_zip_header(archive: ZipFile, member: str) -> list[str]:
    with archive.open(member) as handle:
        if member.casefold().endswith(".gz"):
            with gzip.GzipFile(fileobj=handle) as gz_handle:
                return _parse_header(gz_handle.readline(4 * 1024 * 1024))
        return _parse_header(handle.readline(4 * 1024 * 1024))


def _read_local_header(path: Path) -> list[str]:
    with path.open("rb") as handle:
        return _parse_header(handle.readline(4 * 1024 * 1024))


def _copy_stream(source, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".partial")
    with partial.open("wb") as output:
        shutil.copyfileobj(source, output, length=4 * 1024 * 1024)
    partial.replace(destination)


def _download_direct(
    url: str,
    destination: Path,
    session: requests.Session,
    timeout: int,
    gzipped: bool,
) -> None:
    with session.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        response.raw.decode_content = True
        if gzipped:
            with gzip.GzipFile(fileobj=response.raw) as handle:
                _copy_stream(handle, destination)
        else:
            _copy_stream(response.raw, destination)


def _download_whole(
    url: str,
    destination: Path,
    session: requests.Session,
    timeout: int,
) -> None:
    with session.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with destination.open("wb") as output:
            for chunk in response.iter_content(4 * 1024 * 1024):
                if chunk:
                    output.write(chunk)


def _extract_zip_member(archive: ZipFile, member: str, destination: Path) -> None:
    with archive.open(member) as handle:
        if member.casefold().endswith(".gz"):
            with gzip.GzipFile(fileobj=handle) as gz_handle:
                _copy_stream(gz_handle, destination)
        else:
            _copy_stream(handle, destination)


def _source_address(url: str, member: str | None = None) -> str:
    return url if member is None else f"{url}!/{member.lstrip('/')}"


def _output_path(
    download_dir: Path,
    accession: str,
    address: str,
    archive_name: str,
    member: str | None = None,
) -> Path:
    if member is None:
        label = "direct"
    else:
        archive_label = re.sub(r"\.zip$", "", base_name(archive_name), flags=re.I)
        member_parent = str(PurePosixPath(member).parent)
        label = (
            archive_label
            if member_parent == "."
            else f"{archive_label}_{member_parent}"
        )
    label = re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("._-") or "source"
    digest = hashlib.sha256(address.encode("utf-8")).hexdigest()[:12]
    folder = download_dir / accession / "proteinGroups_files"
    folder.mkdir(parents=True, exist_ok=True)
    return folder / f"{accession}__{label[:80]}__{digest}.txt"


def _describe_gene_sources(header: Sequence[object]) -> str:
    gene_column, fasta_columns, majority_column, _ = gene_sources(header)
    parts: list[str] = []
    if gene_column:
        parts.append(f"gene column: {gene_column}")
    if fasta_columns:
        parts.append("FASTA GN=")
    if majority_column:
        parts.append(f"UniProt fallback: {majority_column}")
    return " + ".join(parts) if parts else "none"


def _source_record(
    address: str,
    local_file: Path,
    archive_name: str,
    header: Sequence[object],
    member: str | None = None,
) -> dict:
    _, _, _, lfq_columns = gene_sources(header)
    return {
        "proteinGroups_address": address,
        "local_file": str(local_file),
        "archive_name": archive_name if member is not None else "",
        "archive_member": member or "",
        "gene_sources": _describe_gene_sources(header),
        "lfq_count": len(lfq_columns),
        "lfq_columns": list(lfq_columns),
    }


def _inspect_zip_for_lfq(
    archive: ZipFile,
    accession: str,
    archive_name: str,
    url: str,
    download_dir: Path,
    force_redownload: bool,
    source_records: list[dict],
    seen_addresses: set[str],
) -> None:
    for info in archive.infolist():
        if info.is_dir() or not _is_proteingroups_file(info.filename):
            continue
        member = info.filename
        address = _source_address(url, member)
        if address in seen_addresses:
            continue
        destination = _output_path(
            download_dir,
            accession,
            address,
            archive_name,
            member,
        )
        try:
            if destination.exists() and not force_redownload:
                header = _read_local_header(destination)
            else:
                header = _read_zip_header(archive, member)
            _, _, _, lfq_columns = gene_sources(header)
            if not lfq_columns:
                continue
            if force_redownload or not destination.exists():
                _extract_zip_member(archive, member, destination)
            source_records.append(
                _source_record(
                    address,
                    destination,
                    archive_name,
                    header,
                    member,
                )
            )
            seen_addresses.add(address)
        except Exception as exc:
            print(f"  Skipped ZIP member {member}: {exc}")


def download_lfq_tables(
    accession: str,
    download_dir: str | Path = "PRIDE_downloads",
    *,
    force_redownload: bool = False,
    timeout: int = DEFAULT_TIMEOUT,
    max_full_zip_download_gb: float = 5.0,
    allow_full_zip_download_when_size_unknown: bool = False,
    session: requests.Session | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Download every LFQ-containing proteinGroups table for one accession."""
    accession = accession.strip().upper()
    if not re.fullmatch(r"PXD\d+", accession):
        raise ValueError("accession must look like PXD044319")
    download_dir = Path(download_dir)
    session = session or session_with_retries()
    records = get_pride_files(accession, session, timeout)
    if not records:
        raise FileNotFoundError(f"PRIDE returned no files for {accession}")

    source_records: list[dict] = []
    seen_addresses: set[str] = set()

    for record in records:
        record_name = str(record.get("fileName") or "")
        if not _is_proteingroups_file(record_name):
            continue
        try:
            url = _download_url(record)
            address = _source_address(url)
            if address in seen_addresses:
                continue
            destination = _output_path(
                download_dir,
                accession,
                address,
                record_name,
            )
            gzipped = base_name(record_name).casefold().endswith(".gz")
            if destination.exists() and not force_redownload:
                header = _read_local_header(destination)
            else:
                header = _read_direct_header(
                    url,
                    session,
                    timeout,
                    gzipped,
                )
            _, _, _, lfq_columns = gene_sources(header)
            if not lfq_columns:
                continue
            if force_redownload or not destination.exists():
                _download_direct(
                    url,
                    destination,
                    session,
                    timeout,
                    gzipped,
                )
            source_records.append(
                _source_record(address, destination, record_name, header)
            )
            seen_addresses.add(address)
        except Exception as exc:
            print(f"  Skipped direct file {record_name}: {exc}")

    for record in records:
        archive_name = str(record.get("fileName") or "")
        if not base_name(archive_name).casefold().endswith(".zip"):
            continue
        try:
            url = _download_url(record)
        except Exception as exc:
            print(f"  Skipped ZIP {archive_name}: {exc}")
            continue

        remote_error: Exception | None = None
        try:
            with RemoteZip(
                url,
                session=session,
                timeout=timeout,
                initial_buffer_size=1024 * 1024,
                support_suffix_range=False,
            ) as archive:
                _inspect_zip_for_lfq(
                    archive,
                    accession,
                    archive_name,
                    url,
                    download_dir,
                    force_redownload,
                    source_records,
                    seen_addresses,
                )
            continue
        except Exception as exc:
            remote_error = exc

        size = _file_size(record)
        maximum = int(max_full_zip_download_gb * 1024**3)
        if size > maximum:
            print(
                f"  Skipped full download of {archive_name} "
                f"({size / 1024**3:.2f} GB) after remote inspection failed: "
                f"{remote_error}"
            )
            continue
        if not size and not allow_full_zip_download_when_size_unknown:
            print(
                f"  Skipped full download of {archive_name}; its size is unknown "
                f"and remote inspection failed: {remote_error}"
            )
            continue

        with tempfile.TemporaryDirectory() as temporary_directory:
            local_zip = Path(temporary_directory) / "archive.zip"
            try:
                _download_whole(url, local_zip, session, timeout)
                with ZipFile(local_zip) as archive:
                    _inspect_zip_for_lfq(
                        archive,
                        accession,
                        archive_name,
                        url,
                        download_dir,
                        force_redownload,
                        source_records,
                        seen_addresses,
                    )
            except Exception as exc:
                print(f"  Could not inspect ZIP {archive_name}: {exc}")

    if not source_records:
        raise FileNotFoundError(
            f"No proteinGroups.txt file with LFQ columns was found for {accession}"
        )

    for index, record in enumerate(source_records, start=1):
        record["table_id"] = f"table_{index}"
        record["accession"] = accession
    tables = pd.DataFrame(source_records)

    catalog_rows = [
        {
            "table_id": record["table_id"],
            "LFQ_column_name": column,
            "local_file": record["local_file"],
            "proteinGroups_address": record["proteinGroups_address"],
        }
        for record in source_records
        for column in record["lfq_columns"]
    ]
    return tables, pd.DataFrame(catalog_rows)


# ---------------------------------------------------------------------------
# Gene extraction and UniProt fallback
# ---------------------------------------------------------------------------


def fasta_genes(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return re.findall(r"\bGN=([^\s;]+)", str(value), flags=re.IGNORECASE)


def uniprot_ids(value: object) -> list[str]:
    if pd.isna(value):
        return []
    result: list[str] = []
    for item in str(value).split(";"):
        item = item.strip()
        if not item:
            continue
        if "|" in item:
            parts = item.split("|")
            if len(parts) >= 2:
                item = parts[1].strip()
        if not item:
            continue
        item = re.sub(r"-\d+$", "", item.split()[0]).upper()
        if UNIPROT_RE.fullmatch(item) and item not in result:
            result.append(item)
    return result


def load_uniprot_cache(path: str | Path) -> dict[str, list[str]]:
    path = Path(path)
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return {
            str(key).upper(): list(dict.fromkeys(map(str, value)))
            for key, value in data.items()
            if isinstance(value, list)
        }
    except Exception as exc:
        print(f"Could not read UniProt cache; starting empty: {exc}")
        return {}


def save_uniprot_cache(cache: Mapping[str, list[str]], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(cache, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def map_uniprot(
    accessions: Iterable[str],
    session: requests.Session,
    cache: dict[str, list[str]],
    timeout: int = DEFAULT_TIMEOUT,
) -> dict[str, list[str]]:
    pending = [
        accession
        for accession in dict.fromkeys(accessions)
        if accession not in cache
    ]
    for start in range(0, len(pending), 10000):
        batch = pending[start : start + 10000]
        response = session.post(
            f"{UNIPROT_API}/idmapping/run",
            data={
                "from": "UniProtKB_AC-ID",
                "to": "Gene_Name",
                "ids": ",".join(batch),
            },
            timeout=timeout,
        )
        response.raise_for_status()
        job = response.json()["jobId"]

        while True:
            status_response = session.get(
                f"{UNIPROT_API}/idmapping/status/{job}",
                timeout=timeout,
            )
            status_response.raise_for_status()
            status = status_response.json().get("jobStatus")
            if status in {"NEW", "RUNNING"}:
                time.sleep(3)
                continue
            if status == "FAILED":
                raise RuntimeError("UniProt mapping failed")
            break

        result_response = session.get(
            f"{UNIPROT_API}/idmapping/stream/{job}",
            params={"format": "json"},
            timeout=timeout,
        )
        result_response.raise_for_status()
        for accession in batch:
            cache.setdefault(accession, [])
        for item in result_response.json().get("results", []):
            accession = str(item.get("from", "")).upper()
            mapped = item.get("to")
            if isinstance(mapped, dict):
                gene = (
                    mapped.get("value")
                    or mapped.get("geneName")
                    or mapped.get("id")
                )
            else:
                gene = mapped
            if gene:
                gene = str(gene).strip()
                if gene and gene not in cache[accession]:
                    cache[accession].append(gene)
    return cache


def create_presence_dataframe(
    file_path: str | Path,
    session: requests.Session,
    uniprot_cache: dict[str, list[str]],
    *,
    timeout: int = DEFAULT_TIMEOUT,
    uniprot_cache_file: str | Path | None = None,
) -> pd.DataFrame:
    """Extract genes and calculate presence across every LFQ column in the table."""
    file_path = Path(file_path)
    header = list(pd.read_csv(file_path, sep="\t", nrows=0).columns)
    gene_column, fasta_columns, majority_column, lfq_columns = gene_sources(header)
    if not lfq_columns:
        raise ValueError(f"No LFQ intensity columns found in {file_path}")
    if not (gene_column or fasta_columns or majority_column):
        raise ValueError(f"No usable gene source found in {file_path}")

    use_columns = list(
        dict.fromkeys(
            ([gene_column] if gene_column else [])
            + fasta_columns
            + ([majority_column] if majority_column else [])
            + lfq_columns
        )
    )
    frame = pd.read_csv(
        file_path,
        sep="\t",
        usecols=use_columns,
        low_memory=False,
    )
    genes_output: list[object] = []
    sources_output: list[str] = []
    unresolved_rows: dict[int, list[str]] = {}
    all_uniprot: list[str] = []

    for row_index, row in frame.iterrows():
        genes: list[str] = []
        seen: set[str] = set()
        matching_sources: list[str] = []

        if gene_column and pd.notna(row[gene_column]):
            for gene in str(row[gene_column]).split(";"):
                gene = gene.strip()
                if gene and gene.upper() not in seen:
                    seen.add(gene.upper())
                    genes.append(gene)
            if genes:
                matching_sources.append(f"Gene column ({gene_column})")

        fasta_added = False
        for column in fasta_columns:
            for gene in fasta_genes(row[column]):
                gene = gene.strip()
                if gene and gene.upper() not in seen:
                    seen.add(gene.upper())
                    genes.append(gene)
                    fasta_added = True
        if fasta_added:
            matching_sources.append("FASTA GN=")

        if genes:
            genes_output.append(";".join(genes))
            sources_output.append(" + ".join(matching_sources))
        else:
            genes_output.append(pd.NA)
            sources_output.append("Unresolved")
            if majority_column:
                ids = uniprot_ids(row[majority_column])
                if ids:
                    unresolved_rows[row_index] = ids
                    all_uniprot.extend(ids)

    frame["Gene names"] = genes_output
    frame["Gene matching source"] = sources_output

    if unresolved_rows:
        try:
            map_uniprot(all_uniprot, session, uniprot_cache, timeout)
            if uniprot_cache_file is not None:
                save_uniprot_cache(uniprot_cache, uniprot_cache_file)
            for row_index, ids in unresolved_rows.items():
                genes: list[str] = []
                for accession in ids:
                    for gene in uniprot_cache.get(accession, []):
                        if gene not in genes:
                            genes.append(gene)
                if genes:
                    frame.at[row_index, "Gene names"] = ";".join(genes)
                    frame.at[row_index, "Gene matching source"] = (
                        "Majority protein IDs -> UniProt"
                    )
        except Exception as exc:
            print(f"UniProt fallback failed; unresolved rows remain blank: {exc}")

    numeric_lfq = frame[lfq_columns].apply(pd.to_numeric, errors="coerce")
    observed = numeric_lfq.notna() & numeric_lfq.ne(0)
    frame[lfq_columns] = numeric_lfq
    frame["Presence_count"] = observed.sum(axis=1)
    frame["Total_LFQ_columns"] = len(lfq_columns)
    frame["Presence"] = frame["Presence_count"] / len(lfq_columns)
    frame["Presence_percent"] = frame["Presence"] * 100

    return frame[
        ["Gene names", "Gene matching source"]
        + list(lfq_columns)
        + [
            "Presence_count",
            "Total_LFQ_columns",
            "Presence",
            "Presence_percent",
        ]
    ]


# ---------------------------------------------------------------------------
# Feature input, grouping, and reports
# ---------------------------------------------------------------------------


def read_feature_genes(path: str | Path) -> list[str]:
    """Read first-column genes from TXT/CSV/TSV/Excel and split delimiters."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"features_list file not found: {path.resolve()}")

    suffix = path.suffix.casefold()
    if suffix in {".xlsx", ".xls"}:
        raw_values = pd.read_excel(path, header=None, usecols=[0]).iloc[:, 0]
    elif suffix == ".csv":
        raw_values = pd.read_csv(path, header=None, usecols=[0]).iloc[:, 0]
    elif suffix == ".tsv":
        raw_values = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0],
        ).iloc[:, 0]
    else:
        raw_values = path.read_text(encoding="utf-8-sig").splitlines()

    genes: list[str] = []
    seen: set[str] = set()
    for raw_value in raw_values:
        if pd.isna(raw_value):
            continue
        text = str(raw_value).split("#", 1)[0].strip()
        if not text:
            continue
        for gene in re.split(r"[;,\t]+", text):
            gene = gene.strip()
            key = gene.upper()
            if (
                gene
                and gene.casefold() not in IGNORED_FEATURE_HEADERS
                and key not in seen
            ):
                genes.append(gene)
                seen.add(key)
    if not genes:
        raise ValueError(f"No feature genes found in {path}")
    return genes


def available_genes(frame: pd.DataFrame, mask=None) -> set[str]:
    if mask is None:
        values = frame["Gene names"]
    else:
        values = frame.loc[mask, "Gene names"]
    return set(
        values.dropna().str.split(";").explode().str.strip().str.upper()
    ) - {""}


def _validate_aggregation(method: str, label: str) -> None:
    if method not in SUPPORTED_AGGREGATIONS:
        raise ValueError(f"{label} must be one of {sorted(SUPPORTED_AGGREGATIONS)}")


def _aggregate_rows(
    frame: pd.DataFrame,
    columns: Sequence[str],
    method: str,
    zero_as_missing: bool,
) -> pd.Series:
    values = frame[list(columns)].apply(pd.to_numeric, errors="coerce")
    if zero_as_missing:
        values = values.mask(values.eq(0))
    if method == "sum":
        return values.sum(axis=1, min_count=1)
    return getattr(values, method)(axis=1, skipna=True)


def _aggregate_duplicate_genes(frame: pd.DataFrame, method: str) -> pd.Series:
    grouped = frame.groupby("Gene", sort=False)["representative"]
    if method == "sum":
        return grouped.sum(min_count=1)
    return grouped.agg(method)


def resolve_group_tables(
    groups: Mapping[str, Sequence[str]],
    group_tables: Mapping[str, str] | None,
    downloaded_tables: pd.DataFrame,
) -> dict[str, str]:
    if not isinstance(groups, Mapping) or not groups:
        raise ValueError("groups must be a non-empty dictionary")
    group_tables = dict(group_tables or {})
    table_lookup = downloaded_tables.set_index("table_id").to_dict("index")
    unknown_groups = set(group_tables) - set(groups)
    if unknown_groups:
        raise KeyError(f"group_tables contains unknown groups: {sorted(unknown_groups)}")

    resolved: dict[str, str] = {}
    for group_name, columns in groups.items():
        if not str(group_name).strip():
            raise ValueError("group names cannot be empty")
        if not isinstance(columns, (list, tuple)) or not columns:
            raise ValueError(f"Group {group_name!r} must contain LFQ columns")
        if len(columns) != len(set(columns)):
            raise ValueError(f"Group {group_name!r} contains duplicate columns")

        explicit_table = group_tables.get(group_name)
        if explicit_table is not None:
            if explicit_table not in table_lookup:
                raise KeyError(f"Unknown table_id {explicit_table!r}")
            missing = [
                column
                for column in columns
                if column not in table_lookup[explicit_table]["lfq_columns"]
            ]
            if missing:
                raise KeyError(
                    f"Group {group_name!r}: columns absent from "
                    f"{explicit_table}: {missing}"
                )
            resolved[group_name] = explicit_table
            continue

        candidates = [
            table_id
            for table_id, metadata in table_lookup.items()
            if set(columns).issubset(set(metadata["lfq_columns"]))
        ]
        if len(candidates) == 1:
            resolved[group_name] = candidates[0]
        elif not candidates:
            all_columns = {
                column
                for metadata in table_lookup.values()
                for column in metadata["lfq_columns"]
            }
            absent = [column for column in columns if column not in all_columns]
            if absent:
                raise KeyError(
                    f"Group {group_name!r}: LFQ columns not found: {absent}"
                )
            raise ValueError(
                f"Group {group_name!r} spans multiple proteinGroups tables"
            )
        else:
            raise ValueError(
                f"Group {group_name!r} matches multiple tables {candidates}; "
                "set group_tables"
            )
    return resolved


def _safe_label(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_")


def output_column_name(
    accession: str,
    group: str,
    tissue_or_cell_culture: str,
    disease_type: str,
) -> str:
    parts = [
        _safe_label(value)
        for value in (
            accession,
            group,
            tissue_or_cell_culture,
            disease_type,
        )
    ]
    if any(not part for part in parts):
        raise ValueError("Output-name components cannot be empty")
    return "_".join(parts)


def _write_report(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.casefold()
    if suffix == ".xlsx":
        frame.to_excel(path, index=False)
    elif suffix == ".csv":
        frame.to_csv(path, index=False)
    elif suffix == ".tsv":
        frame.to_csv(path, sep="\t", index=False)
    else:
        raise ValueError("report_file must end with .xlsx, .csv, or .tsv")


def _show_missing_as_na(frame: pd.DataFrame) -> pd.DataFrame:
    """Keep the Gene column unchanged and write missing report values as 'NA'."""
    report = frame.copy()
    for column in report.columns[1:]:
        report[column] = report[column].astype(object).where(
            report[column].notna(),
            "NA",
        )
    return report


def build_representative_report(
    *,
    accession: str,
    tissue_or_cell_culture: str,
    disease_type: str,
    features_list_file: str | Path,
    report_file: str | Path,
    groups: Mapping[str, Sequence[str]],
    downloaded_tables: pd.DataFrame,
    group_tables: Mapping[str, str] | None = None,
    presence_threshold: float = 0.10,
    group_aggregation: str = "mean",
    duplicate_gene_aggregation: str = "mean",
    treat_zero_as_missing: bool = False,
    allow_column_reuse: bool = False,
    download_dir: str | Path = "PRIDE_downloads",
    timeout: int = DEFAULT_TIMEOUT,
    session: requests.Session | None = None,
) -> dict[str, pd.DataFrame | Path]:
    """Build group representatives after strict whole-table presence filtering.

    Matching and numerical inclusion are intentionally reported separately:

    * annotation match: target symbol exists anywhere after all gene-mapping routes;
    * eligible: at least one matched protein row has Presence > threshold, where
      Presence is calculated over every LFQ column in its full proteinGroups table;
    * group value: eligible genes retain their numerical group representative,
      including 0; unmatched or ineligible genes remain NA.
    """
    accession = accession.strip().upper()
    if not 0 <= presence_threshold <= 1:
        raise ValueError("presence_threshold must be between 0 and 1")
    _validate_aggregation(group_aggregation, "group_aggregation")
    _validate_aggregation(
        duplicate_gene_aggregation,
        "duplicate_gene_aggregation",
    )

    resolved_tables = resolve_group_tables(
        groups,
        group_tables,
        downloaded_tables,
    )
    if not allow_column_reuse:
        used: dict[tuple[str, str], str] = {}
        for group_name, columns in groups.items():
            table_id = resolved_tables[group_name]
            for column in columns:
                key = (table_id, column)
                if key in used:
                    raise ValueError(
                        f"{column!r} from {table_id} is used in both "
                        f"{used[key]!r} and {group_name!r}"
                    )
                used[key] = group_name

    features = read_feature_genes(features_list_file)
    feature_keys = [gene.upper() for gene in features]
    target = set(feature_keys)
    table_lookup = downloaded_tables.set_index("table_id").to_dict("index")
    session = session or session_with_retries()
    cache_file = Path(download_dir) / "uniprot_gene_cache.json"
    uniprot_cache = load_uniprot_cache(cache_file)

    loaded_tables: dict[str, pd.DataFrame] = {}
    matching_rows: list[dict] = []
    feature_status_rows: list[dict] = []

    for table_id in dict.fromkeys(resolved_tables.values()):
        frame = create_presence_dataframe(
            table_lookup[table_id]["local_file"],
            session,
            uniprot_cache,
            timeout=timeout,
            uniprot_cache_file=cache_file,
        )
        loaded_tables[table_id] = frame

        annotation_genes = available_genes(frame)
        eligible_mask = frame["Presence"].gt(presence_threshold)
        eligible_genes = available_genes(frame, eligible_mask)
        annotation_matched = target & annotation_genes
        eligible_matched = target & eligible_genes

        matching_rows.append(
            {
                "table_id": table_id,
                "target_genes": len(target),
                "annotation_matched": len(annotation_matched),
                "annotation_missing": len(target - annotation_genes),
                "presence_rule": f"Presence > {presence_threshold:.2f}",
                "whole_table_LFQ_columns": int(
                    frame["Total_LFQ_columns"].iloc[0]
                ),
                "eligible_after_presence": len(eligible_matched),
                "matched_but_not_eligible": len(
                    annotation_matched - eligible_genes
                ),
            }
        )
        for feature in features:
            key = feature.upper()
            feature_status_rows.append(
                {
                    "table_id": table_id,
                    "Gene": feature,
                    "annotation_matched": key in annotation_genes,
                    "eligible_after_presence": key in eligible_genes,
                }
            )

    new_columns: dict[str, list[object]] = {"Gene": features}
    group_rows: list[dict] = []

    for group_name, columns in groups.items():
        table_id = resolved_tables[group_name]
        frame = loaded_tables[table_id]
        eligible_mask = frame["Presence"].gt(presence_threshold)
        representative = _aggregate_rows(
            frame,
            columns,
            group_aggregation,
            treat_zero_as_missing,
        )
        gene_values = pd.DataFrame(
            {
                "Gene": frame.loc[eligible_mask, "Gene names"].str.split(";"),
                "representative": representative.loc[eligible_mask],
            }
        ).explode("Gene")
        gene_values["Gene"] = (
            gene_values["Gene"].astype("string").str.strip().str.upper()
        )
        gene_values = gene_values[
            gene_values["Gene"].notna() & gene_values["Gene"].ne("")
        ]
        gene_map = _aggregate_duplicate_genes(
            gene_values,
            duplicate_gene_aggregation,
        )

        output_name = output_column_name(
            accession,
            group_name,
            tissue_or_cell_culture,
            disease_type,
        )
        if output_name in new_columns:
            raise ValueError(
                f"Two groups produce the same output name: {output_name}"
            )
        values = [gene_map.get(key, float("nan")) for key in feature_keys]
        new_columns[output_name] = values
        group_value_count = int(pd.Series(values).notna().sum())

        table_match = next(
            row for row in matching_rows if row["table_id"] == table_id
        )
        group_rows.append(
            {
                "group": group_name,
                "table_id": table_id,
                "LFQ_columns": len(columns),
                "aggregation": group_aggregation,
                "annotation_matched": table_match["annotation_matched"],
                "eligible_after_presence": table_match[
                    "eligible_after_presence"
                ],
                "group_values_written": group_value_count,
                "output_column": output_name,
            }
        )

    # A run represents the groups supplied now. Do not retain columns from an
    # older report: the output is Gene + exactly one column per current group.
    report = _show_missing_as_na(pd.DataFrame(new_columns))
    report_file = Path(report_file)
    _write_report(report, report_file)

    return {
        "report": report,
        "report_file": report_file.resolve(),
        "matching_summary": pd.DataFrame(matching_rows),
        "group_summary": pd.DataFrame(group_rows),
        "feature_status": pd.DataFrame(feature_status_rows),
    }


def _display(frame: pd.DataFrame) -> None:
    try:
        from IPython.display import display

        display(frame)
    except ImportError:
        print(frame.to_string(index=False))


def download_and_show(
    dataset: Mapping[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Notebook-facing wrapper: download tables and show the short LFQ report."""
    tables, catalog = download_lfq_tables(
        str(dataset["accession"]),
        dataset.get("download_dir", "PRIDE_downloads"),
        force_redownload=bool(dataset.get("force_redownload", False)),
        timeout=int(dataset.get("timeout", DEFAULT_TIMEOUT)),
        max_full_zip_download_gb=float(
            dataset.get("max_full_zip_download_gb", 5.0)
        ),
        allow_full_zip_download_when_size_unknown=bool(
            dataset.get("allow_full_zip_download_when_size_unknown", False)
        ),
    )
    print(f"Found {len(tables)} LFQ-containing proteinGroups table(s).")
    _display(
        tables[["table_id", "lfq_count", "gene_sources", "local_file"]]
    )
    print(f"LFQ columns ({len(catalog)} total):")
    with pd.option_context("display.max_rows", None, "display.max_colwidth", 120):
        _display(catalog[["table_id", "LFQ_column_name"]])
    return tables, catalog


def analyze_and_show(
    dataset: Mapping[str, object],
    groups: Mapping[str, Sequence[str]],
    settings: Mapping[str, object],
    downloaded_tables: pd.DataFrame,
) -> dict[str, pd.DataFrame | Path]:
    """Notebook-facing wrapper: run the analysis and show compact results."""
    result = build_representative_report(
        accession=str(dataset["accession"]),
        tissue_or_cell_culture=str(dataset["tissue_or_cell_culture"]),
        disease_type=str(dataset["disease_type"]),
        features_list_file=dataset["features_list_file"],
        report_file=dataset["report_file"],
        groups=groups,
        downloaded_tables=downloaded_tables,
        group_tables=settings.get("group_tables", {}),
        presence_threshold=float(settings.get("presence_threshold", 0.10)),
        group_aggregation=str(settings.get("group_aggregation", "mean")),
        duplicate_gene_aggregation=str(
            settings.get("duplicate_gene_aggregation", "mean")
        ),
        # Whole-table presence controls eligibility. Within an eligible row,
        # zero is a valid group value and must remain zero in the report.
        treat_zero_as_missing=False,
        allow_column_reuse=bool(settings.get("allow_column_reuse", False)),
        download_dir=dataset.get("download_dir", "PRIDE_downloads"),
        timeout=int(dataset.get("timeout", DEFAULT_TIMEOUT)),
    )
    print(f"Report saved to: {result['report_file']}")
    print("\nGene matching and whole-table presence:")
    _display(result["matching_summary"])
    print("\nGroup representatives:")
    _display(result["group_summary"])
    print("\nReport preview:")
    _display(result["report"].head(20))
    return result


__all__ = [
    "analyze_and_show",
    "available_genes",
    "build_representative_report",
    "create_presence_dataframe",
    "download_and_show",
    "download_lfq_tables",
    "fasta_genes",
    "gene_sources",
    "get_pride_files",
    "load_uniprot_cache",
    "map_uniprot",
    "output_column_name",
    "read_feature_genes",
    "resolve_group_tables",
    "save_uniprot_cache",
    "session_with_retries",
    "uniprot_ids",
]
