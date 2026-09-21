"""PRIDE LFQ download, numbered column selection, and row-wise aggregation.

Notebook workflow (the same three calls as before):
    tables, catalog = prepare_lfq_selection(DATASET)
    preview_lfq_range_groups(GROUPS, catalog)
    result = process_lfq_ranges_and_show(DATASET, GROUPS, SETTINGS, tables, catalog)

GROUPS uses inclusive, 1-based LFQ catalogue ranges, for example:
    GROUPS = {"control": [1, 3], "treated": [[4, 6], [10, 12]]}
    SETTINGS = {
        "group_aggregation": "median",  # mean, min, max, and sum also supported
        "presence_threshold": 0.10,
        "treat_zero_as_missing": False,
        "append_to_existing": False,  # True adds groups to an existing row report
        "report_file": "protein_rows_lfq_report.xlsx",
    }

Each output row is one original proteinGroups row, in source order. A row
qualifies when MORE THAN 10% of ALL LFQ columns in its source table contain
non-missing, nonzero values. Failed rows remain, with NA aggregate values.
The LFQ calculation is unchanged: aggregate the selected columns within that
row, ignore missing values, and include zeros unless configured otherwise.

The five annotation columns are copied as text from exact or unambiguous
header equivalents; unavailable fields stay empty. Nothing is extracted from
FASTA headers, mapped through UniProt, split by gene, or merged across rows.
Even identical rows remain separate. Source file and Source row (1-based data
row, excluding the header) identify rows when appending reports. A group that
spans tables is calculated separately for each table's rows.

Requires pandas; Excel output also requires openpyxl. PRIDE downloading
additionally requires requests and remotezip. Local processing makes no network
requests. The old feature-list workflow and duplicate-gene settings are removed;
unrelated keys in an existing SETTINGS dictionary are ignored.
"""

from __future__ import annotations

from numbers import Integral
import gzip
import hashlib
import re
import shutil
import tempfile
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING, Mapping, Sequence
from zipfile import ZipFile

import pandas as pd

if TYPE_CHECKING:
    import requests


PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
DEFAULT_TIMEOUT = 300
ANNOTATION_ALIASES = {
    "Protein IDs": ("protein ids", "protein id"),
    "Majority protein IDs": ("majority protein ids", "majority protein id"),
    "Protein names": ("protein names", "protein name"),
    "Gene names": ("gene names", "gene name", "gene symbols", "gene symbol"),
    "Fasta headers": ("fasta headers", "fasta header"),
    "Only identified by site": (
        "only identified by site",
        "onlyidentifiedbysite",
    ),
    "Reverse": ("reverse",),
    "Potential contaminant": (
        "potential contaminant",
        "potentialcontaminant"
    )
}
ANNOTATION_COLUMNS = list(ANNOTATION_ALIASES)
SOURCE_COLUMNS = ["Source file", "Source row"]
REPORT_COLUMNS = ANNOTATION_COLUMNS + SOURCE_COLUMNS
LFQ_RE = re.compile(r"^LFQ\s+intensity(?:\s|$)", re.IGNORECASE)
SUPPORTED_AGGREGATIONS = {"mean", "median", "max", "min", "sum"}

LFQ_SAMPLE_PREFIX_RE = re.compile(
    r"^LFQ\s+intensity\s+",
    re.IGNORECASE,
)
REPLICATE_SUFFIX_PATTERNS = [
    (
        re.compile(
            r"^(.*?)(?:[_ .-]?(rep(?:licate)?|ex|run)"
            r"[_-]?([A-Za-z0-9]+))$",
            re.IGNORECASE,
        ),
        "labelled suffix",
    ),
    (
        re.compile(r"^(.*?)[_ .-]([A-Za-z])$"),
        "letter suffix",
    ),
    (
        re.compile(r"^(.*?)[_ .-](\d+)$"),
        "numeric suffix",
    ),
    (
        re.compile(r"^(.*?)([A-Da-d])$"),
        "attached A-D suffix",
    ),
]


def session_with_retries() -> requests.Session:
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

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



def lfq_columns(columns: Sequence[object]) -> list[object]:
    """Identify LFQ columns using the original LFQ intensity header rule."""
    return [column for column in columns if LFQ_RE.match(str(column).strip())]


def annotation_columns(columns: Sequence[object]) -> dict[str, object | None]:
    """Use exact headers first, then unambiguous spelling equivalents only."""
    result = {}
    for target, aliases in ANNOTATION_ALIASES.items():
        exact = [column for column in columns if column == target]
        candidates = exact or [
            column for column in columns if normalized_text(column) in aliases
        ]
        result[target] = candidates[0] if len(candidates) == 1 else None
    return result


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



def _source_record(
    address: str,
    local_file: Path,
    archive_name: str,
    header: Sequence[object],
    member: str | None = None,
) -> dict:
    found_lfq = lfq_columns(header)
    return {
        "proteinGroups_address": address,
        "local_file": str(local_file),
        "archive_name": archive_name if member is not None else "",
        "archive_member": member or "",
        "lfq_count": len(found_lfq),
        "lfq_columns": list(found_lfq),
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
            found_lfq = lfq_columns(header)
            if not found_lfq:
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
    from remotezip import RemoteZip

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
            found_lfq = lfq_columns(header)
            if not found_lfq:
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

    source_records.sort(key=lambda record: record["proteinGroups_address"])
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



def create_presence_dataframe(file_path: str | Path) -> pd.DataFrame:
    """Copy source annotations and compute the unchanged whole-table presence."""
    options = dict(sep="\t", encoding="utf-8-sig", keep_default_na=False,
                   skip_blank_lines=False, low_memory=False)
    header = list(pd.read_csv(file_path, nrows=0, **options).columns)
    found_lfq = lfq_columns(header)
    if not found_lfq:
        raise ValueError(f"No LFQ intensity columns found in {file_path}")
    sources = annotation_columns(header)
    source_annotations = [column for column in sources.values() if column is not None]
    raw = pd.read_csv(
        file_path,
        usecols=source_annotations + found_lfq,
        dtype={column: str for column in source_annotations},
        **options,
    )
    frame = pd.DataFrame(index=raw.index)
    for target, column in sources.items():
        frame[target] = raw[column].fillna("") if column is not None else ""

    numeric_lfq = raw[found_lfq].apply(pd.to_numeric, errors="coerce")
    observed = numeric_lfq.notna() & numeric_lfq.ne(0)
    frame[found_lfq] = numeric_lfq
    frame["Presence_count"] = observed.sum(axis=1)
    frame["Total_LFQ_columns"] = len(found_lfq)
    frame["Presence"] = frame["Presence_count"] / len(found_lfq)
    frame["Presence_percent"] = frame["Presence"] * 100
    return frame


def _split_replicate_suffix(sample_name: str) -> tuple[str, str, str]:
    """Split a conservative final replicate suffix from one LFQ sample name."""
    for pattern, rule in REPLICATE_SUFFIX_PATTERNS:
        match = pattern.match(sample_name)
        if not match:
            continue
        stem = match.group(1).rstrip("_ .-")
        suffix = sample_name[len(match.group(1)) :].lstrip("_ .-")
        if stem:
            return stem, suffix, rule
    return sample_name, "", "no recognized suffix"



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



def _safe_label(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_")



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



def _display(frame: pd.DataFrame) -> None:
    try:
        from IPython.display import display

        display(frame)
    except ImportError:
        print(frame.to_string(index=frame.index.name is not None))



def _natural_key(value: object) -> tuple:
    """Natural ordering puts replicate 2 before replicate 10."""
    return tuple(
        (1, int(part)) if part.isdigit() else (0, part.casefold())
        for part in re.split(r"(\d+)", str(value))
    )



def _dataset_identity(dataset: Mapping[str, object]) -> tuple[str, str]:
    accession = str(dataset.get("accession", "")).strip().upper()
    dataset_id = str(dataset.get("ID", "")).strip()
    if not re.fullmatch(r"PXD\d+", accession):
        raise ValueError("Set DATASET['accession'] to a PRIDE accession such as PXD010393")
    if not dataset_id or dataset_id.casefold() in {"none", "your_id", "enter_id"}:
        raise ValueError("Set DATASET['ID'] to your dataset ID before running Cell 1")
    if not _safe_label(dataset_id):
        raise ValueError("Dataset ID must contain letters or numbers")
    return accession, dataset_id



def make_numbered_lfq_catalog(lfq_catalog: pd.DataFrame) -> pd.DataFrame:
    """Create one stable 1..N row numbering over LFQ columns from ALL tables.

    Similarity means shared replicate-stripped stem, then natural sample-name
    ordering. Ties use source address and table ID. Rows are never collapsed
    merely because two tables have the same LFQ header label.
    """
    required = {"table_id", "LFQ_column_name"}
    if not required.issubset(lfq_catalog.columns):
        raise ValueError(f"LFQ catalogue is missing {sorted(required - set(lfq_catalog.columns))}")
    if lfq_catalog.empty:
        raise ValueError("There are no LFQ columns to select")
    catalog = lfq_catalog.drop(columns=["row_number"], errors="ignore").copy()
    if catalog[list(required)].isna().any().any():
        raise ValueError("LFQ catalogue has a missing table ID or column label")
    for column in required:
        catalog[column] = catalog[column].astype(str)
    if catalog.duplicated(["table_id", "LFQ_column_name"]).any():
        raise ValueError("Duplicate LFQ labels within one table cannot be selected unambiguously")
    catalog["sample_name"] = catalog["LFQ_column_name"].map(
        lambda value: LFQ_SAMPLE_PREFIX_RE.sub("", value.strip()).strip()
    )
    catalog["similarity_stem"] = catalog["sample_name"].map(
        lambda sample: _split_replicate_suffix(sample)[0]
    )
    records = catalog.to_dict("records")
    records.sort(key=lambda record: (
        _natural_key(normalized_text(record["similarity_stem"])),
        _natural_key(record["sample_name"]),
        _natural_key(record.get("proteinGroups_address", "")),
        _natural_key(record["table_id"]),
        record["LFQ_column_name"],
    ))
    catalog = pd.DataFrame(records)
    catalog.insert(0, "row_number", range(1, len(catalog) + 1))
    return catalog



def prepare_lfq_selection(
    dataset: Mapping[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cell 1: download/cache tables and display all selectable LFQ rows."""
    accession, dataset_id = _dataset_identity(dataset)
    tables, raw_catalog = download_lfq_tables(
        accession,
        dataset.get("download_dir", "PRIDE_downloads"),
        force_redownload=bool(dataset.get("force_redownload", False)),
        timeout=int(dataset.get("timeout", DEFAULT_TIMEOUT)),
        max_full_zip_download_gb=float(dataset.get("max_full_zip_download_gb", 5.0)),
        allow_full_zip_download_when_size_unknown=bool(
            dataset.get("allow_full_zip_download_when_size_unknown", False)
        ),
    )
    catalog = make_numbered_lfq_catalog(raw_catalog)
    catalog.attrs.update({"accession": accession, "ID": dataset_id})
    print(f"{dataset_id} / {accession}: {len(tables)} table(s), {len(catalog)} LFQ columns.")

    # Display the analysis table ID together with its source URL/address.
    print("Table ID -> URL:")
    _display(
        tables[["table_id", "proteinGroups_address"]]
        .rename(columns={"proteinGroups_address": "URL"})
    )

    print("Select by row_number (1-based). Similar labels from all tables are shown together.")
    with pd.option_context("display.max_rows", None, "display.max_colwidth", None,
                           "display.max_columns", None):
        _display(catalog[["row_number", "sample_name", "table_id", "LFQ_column_name"]]
                 .set_index("row_number"))
    return tables, catalog



def _inclusive_ranges(selection: object, group_name: str, maximum: int) -> list[int]:
    """[1, 3] -> [1, 2, 3]; [[1, 3], [7, 8]] -> [1, 2, 3, 7, 8]."""
    def is_integer(value):
        return isinstance(value, Integral) and not isinstance(value, bool)

    if not isinstance(selection, (list, tuple)) or not selection:
        raise ValueError(f"Group {group_name!r}: use [start, end] or [[start, end], ...]")
    intervals = [selection] if len(selection) == 2 and all(map(is_integer, selection)) else selection
    rows: list[int] = []
    for interval in intervals:
        if (not isinstance(interval, (list, tuple)) or len(interval) != 2
                or not all(map(is_integer, interval))):
            raise ValueError(f"Group {group_name!r}: each range must be [start, end], both integers")
        start, end = map(int, interval)
        if start < 1 or end < start:
            raise ValueError(f"Group {group_name!r}: ranges must satisfy 1 <= start <= end")
        if end > maximum:
            raise ValueError(f"Group {group_name!r}: range [{start}, {end}] is outside 1..{maximum}")
        rows.extend(range(start, end + 1))
    if len(rows) != len(set(rows)):
        raise ValueError(f"Group {group_name!r}: overlapping ranges select a row more than once")
    return rows



def resolve_lfq_ranges(
    groups: Mapping[str, object],
    lfq_catalog: pd.DataFrame,
    *,
    allow_column_reuse: bool = False,
) -> dict[str, pd.DataFrame]:
    """Resolve global row numbers into exact source-table/column pairs."""
    required = {"row_number", "table_id", "LFQ_column_name"}
    if not required.issubset(lfq_catalog.columns):
        raise ValueError("Run Cell 1 to create the numbered LFQ catalogue")
    if not isinstance(groups, Mapping) or not groups:
        raise ValueError("Define at least one group in Cell 2, for example {'control': [1, 3]}")
    row_numbers = lfq_catalog["row_number"].tolist()
    if (any(not isinstance(value, Integral) or isinstance(value, bool) for value in row_numbers)
            or sorted(row_numbers) != list(range(1, len(lfq_catalog) + 1))):
        raise ValueError("Catalogue row numbers must be unique integers 1..N; rerun Cell 1")
    if lfq_catalog.duplicated(["table_id", "LFQ_column_name"]).any():
        raise ValueError("Catalogue contains a duplicate source-table/column pair")
    lookup = lfq_catalog.set_index("row_number", drop=False)
    resolved: dict[str, pd.DataFrame] = {}
    used: dict[int, str] = {}
    for group_name, selection in groups.items():
        if not isinstance(group_name, str) or not group_name.strip():
            raise ValueError("Group names must be non-empty strings")
        row_ids = _inclusive_ranges(selection, group_name, len(lookup))
        missing = [value for value in row_ids if value not in lookup.index]
        if missing:
            raise ValueError(f"Group {group_name!r}: row(s) {missing} are outside 1..{len(lookup)}")
        for row_id in row_ids:
            if not allow_column_reuse and row_id in used:
                raise ValueError(f"Row {row_id} is used in both {used[row_id]!r} and {group_name!r}")
            used[row_id] = group_name
        resolved[group_name] = lookup.loc[row_ids].reset_index(drop=True).copy()
    return resolved



def preview_lfq_range_groups(
    groups: Mapping[str, object],
    lfq_catalog: pd.DataFrame,
    *,
    allow_column_reuse: bool = False,
) -> pd.DataFrame:
    """Cell 2: show the exact columns selected by the user's inclusive ranges."""
    resolved = resolve_lfq_ranges(groups, lfq_catalog, allow_column_reuse=allow_column_reuse)
    preview = pd.DataFrame([
        {"group": group, "row_numbers": selected["row_number"].tolist(),
         "n_LFQ_columns": len(selected),
         "tables": "; ".join(dict.fromkeys(selected["table_id"])),
         "LFQ_columns": "; ".join(
             f"{row.row_number}: {row.LFQ_column_name} ({row.table_id})"
             for row in selected.itertuples(index=False))}
        for group, selected in resolved.items()
    ])
    with pd.option_context("display.max_colwidth", None, "display.max_rows", None):
        _display(preview)
    selected_count = len({number for selected in resolved.values() for number in selected["row_number"]})
    print(f"Selected {selected_count} of {len(lfq_catalog)} LFQ columns in {len(resolved)} group(s).")
    return preview


def _merge_range_report(new_report: pd.DataFrame, report_file: Path,
                        append_to_existing: bool) -> pd.DataFrame:
    """Add group columns by source file and row number, never by annotation."""
    if not append_to_existing or not report_file.exists():
        return new_report
    options = dict(keep_default_na=False,
                   dtype={column: str for column in ANNOTATION_COLUMNS + ["Source file"]})
    if report_file.suffix.casefold() == ".xlsx":
        existing = pd.read_excel(report_file, **options)
    else:
        existing = pd.read_csv(
            report_file, sep="\t" if report_file.suffix.casefold() == ".tsv" else ",",
            **options,
        )
    if (not set(REPORT_COLUMNS).issubset(existing.columns)
            or existing.columns.duplicated().any()):
        raise ValueError(
            "Existing report is not a protein-row report. Use a new report_file "
            "or set append_to_existing=False to replace it."
        )
    row_numbers = pd.to_numeric(existing["Source row"], errors="coerce")
    if (row_numbers.isna().any() or row_numbers.lt(1).any()
            or row_numbers.mod(1).ne(0).any()
            or existing["Source file"].str.strip().eq("").any()):
        raise ValueError("Existing report has invalid source file/row identifiers")
    existing["Source row"] = row_numbers.astype("int64")
    if existing.duplicated(SOURCE_COLUMNS).any():
        raise ValueError("Existing report repeats a source file/row identifier")
    for column in existing.columns:
        if column not in REPORT_COLUMNS:
            existing[column] = pd.to_numeric(
                existing[column].mask(existing[column].isin(["", "NA"])),
                errors="raise",
            )

    old = existing.set_index(SOURCE_COLUMNS)
    new = new_report.set_index(SOURCE_COLUMNS)
    old_sources = set(old.index.get_level_values("Source file"))
    for source in dict.fromkeys(new.index.get_level_values("Source file")):
        if source in old_sources:
            old_rows = old.loc[[source], ANNOTATION_COLUMNS]
            new_rows = new.loc[[source], ANNOTATION_COLUMNS]
            if not old_rows.equals(new_rows):
                raise ValueError(
                    "Source row order, count, or annotations changed. "
                    "Use a new report_file or set append_to_existing=False."
                )

    order = old.index.append(new.index[~new.index.isin(old.index)])
    merged = old.reindex(order).copy()
    for column in ANNOTATION_COLUMNS:
        merged.loc[new.index, column] = new[column]
    for column in new.columns.difference(ANNOTATION_COLUMNS, sort=False):
        # Also replace missing values, clearing stale results for a rerun group.
        merged[column] = new[column].reindex(order)
    merged = merged.reset_index()
    value_columns = [column for column in merged if column not in REPORT_COLUMNS]
    return merged[REPORT_COLUMNS + value_columns]


def _show_missing_as_na(frame: pd.DataFrame) -> pd.DataFrame:
    """Leave annotation cells empty; represent missing aggregate values as NA."""
    report = frame.copy()
    for column in ANNOTATION_COLUMNS:
        report[column] = report[column].fillna("")
    for column in report.columns:
        if column not in REPORT_COLUMNS:
            report[column] = report[column].astype(object).where(
                report[column].notna(), "NA"
            )
    return report


def build_range_representative_report(
    *,
    dataset: Mapping[str, object],
    groups: Mapping[str, object],
    lfq_catalog: pd.DataFrame,
    downloaded_tables: pd.DataFrame,
    report_file: str | Path = "processed_lfq_report.xlsx",
    group_aggregation: str = "median",
    presence_threshold: float | None = 0.10,
    treat_zero_as_missing: bool = False,
    append_to_existing: bool = True,
    allow_column_reuse: bool = False,
) -> dict[str, object]:
    """Aggregate LFQ columns within each original row, preserving every row.

    Presence uses every LFQ column in the source table and a strict > comparison.
    Rows failing the rule remain with NA values; None disables the presence rule.
    Rows from different tables are stacked, with each table's order preserved.
    No annotation is used for matching, filtering, expanding, or merging rows.
    """
    accession, dataset_id = _dataset_identity(dataset)
    if presence_threshold is not None and not 0 <= presence_threshold <= 1:
        raise ValueError("presence_threshold must be between 0 and 1, or None")
    _validate_aggregation(group_aggregation, "group_aggregation")
    report_file = Path(report_file)
    if report_file.suffix.casefold() not in {".xlsx", ".csv", ".tsv"}:
        raise ValueError("report_file must end with .xlsx, .csv, or .tsv")
    if downloaded_tables["table_id"].duplicated().any():
        raise ValueError("Downloaded tables have duplicate table IDs")
    if ("accession" in downloaded_tables
            and not downloaded_tables["accession"].astype(str).str.upper().eq(accession).all()):
        raise ValueError("Downloaded tables belong to a different accession; rerun Cell 1")
    if lfq_catalog.attrs.get("accession", accession) != accession:
        raise ValueError("The numbered catalogue belongs to another accession; rerun Cell 1")
    resolved = resolve_lfq_ranges(groups, lfq_catalog, allow_column_reuse=allow_column_reuse)
    output_names = {
        group: "_".join([_safe_label(dataset_id), _safe_label(accession), _safe_label(group)])
        for group in resolved
    }
    if (any(not _safe_label(group) for group in resolved)
            or len(set(output_names.values())) != len(output_names)):
        raise ValueError("Group names produce empty or duplicate output labels; choose distinct names")
    selected = pd.concat(list(resolved.values()), ignore_index=True).drop_duplicates("row_number")
    table_lookup = downloaded_tables.set_index("table_id").to_dict("index")
    for row in selected.itertuples(index=False):
        metadata = table_lookup.get(row.table_id)
        if metadata is None or row.LFQ_column_name not in metadata["lfq_columns"]:
            raise ValueError(f"Row {row.row_number} does not match the current tables; rerun Cell 1")
        for field in ("local_file", "proteinGroups_address"):
            if (field in selected and field in metadata
                    and str(getattr(row, field)) != str(metadata[field])):
                raise ValueError(f"Row {row.row_number} has changed source information; rerun Cell 1")

    table_reports, row_summary = [], []
    selected_table_ids = set(selected["table_id"])
    # Iterate in source-table order, never in gene, protein, or LFQ-sample order.
    for table_id, metadata in table_lookup.items():
        if table_id not in selected_table_ids:
            continue
        frame = create_presence_dataframe(metadata["local_file"])
        if list(lfq_columns(frame.columns)) != list(metadata["lfq_columns"]):
            raise ValueError(f"LFQ headers changed in {table_id}; rerun Cell 1 and review the ranges")
        eligible = (
            pd.Series(True, index=frame.index)
            if presence_threshold is None else frame["Presence"].gt(presence_threshold)
        )
        table_report = frame[ANNOTATION_COLUMNS].copy()
        address = metadata.get("proteinGroups_address", "")
        table_report["Source file"] = (
            str(address) if pd.notna(address) and str(address).strip()
            else str(Path(metadata["local_file"]).resolve())
        )
        table_report["Source row"] = frame.index + 1
        for group, group_selection in resolved.items():
            columns = group_selection.loc[
                group_selection["table_id"].eq(table_id), "LFQ_column_name"
            ].tolist()
            table_report[output_names[group]] = (
                _aggregate_rows(frame.loc[eligible], columns,
                                group_aggregation, treat_zero_as_missing).reindex(frame.index)
                if columns else float("nan")
            )
        table_reports.append(table_report)
        row_summary.append({
            "table_id": table_id,
            "source_rows": len(frame),
            "rows_passing_presence": int(eligible.sum()),
            "rows_failing_presence": int((~eligible).sum()),
            "selected_LFQ_columns": int(selected["table_id"].eq(table_id).sum()),
            "whole_table_LFQ_columns": len(metadata["lfq_columns"]),
            "presence_rule": "disabled" if presence_threshold is None else f"Presence > {presence_threshold:g}",
        })
    new_report = pd.concat(table_reports, ignore_index=True)
    if new_report.duplicated(SOURCE_COLUMNS).any():
        raise ValueError("Selected tables repeat a source file; each table must identify a distinct file")

    group_summary, selection_rows = [], []
    for group, group_selection in resolved.items():
        output_name = output_names[group]
        group_summary.append({
            "group": group,
            "row_numbers": group_selection["row_number"].tolist(),
            "tables": "; ".join(dict.fromkeys(group_selection["table_id"])),
            "LFQ_columns": len(group_selection),
            "aggregation": group_aggregation,
            "protein_rows_with_values": int(new_report[output_name].notna().sum()),
            "output_column": output_name,
        })
        selection_rows.extend(
            {"group": group, "output_column": output_name, **row}
            for row in group_selection.to_dict("records")
        )
    report = _show_missing_as_na(_merge_range_report(new_report, report_file, append_to_existing))
    report_file.parent.mkdir(parents=True, exist_ok=True)
    temporary = report_file.with_name(report_file.stem + ".partial" + report_file.suffix)
    try:
        _write_report(report, temporary)
        temporary.replace(report_file)
    finally:
        if temporary.exists():
            temporary.unlink()
    return {
        "report": report,
        "report_file": report_file.resolve(),
        "new_report": new_report,
        "group_summary": pd.DataFrame(group_summary),
        "row_summary": pd.DataFrame(row_summary),
        "selected_columns": pd.DataFrame(selection_rows),
        "lfq_catalog": lfq_catalog.copy(),
    }


def process_lfq_ranges_and_show(
    dataset: Mapping[str, object],
    groups: Mapping[str, object],
    settings: Mapping[str, object],
    downloaded_tables: pd.DataFrame,
    lfq_catalog: pd.DataFrame,
) -> dict[str, object]:
    """Cell 3: calculate, save, and display the protein-row report."""
    result = build_range_representative_report(
        dataset=dataset, groups=groups, lfq_catalog=lfq_catalog,
        downloaded_tables=downloaded_tables,
        report_file=settings.get("report_file", "processed_lfq_report.xlsx"),
        group_aggregation=str(settings.get("group_aggregation", "median")),
        presence_threshold=settings.get("presence_threshold", 0.10),
        treat_zero_as_missing=bool(settings.get("treat_zero_as_missing", False)),
        append_to_existing=bool(settings.get("append_to_existing", True)),
        allow_column_reuse=bool(settings.get("allow_column_reuse", False)),
    )
    print(f"Report saved to: {result['report_file']}")
    print("\nSource rows and whole-table presence:")
    _display(result["row_summary"])
    print("\nGroup representatives:")
    _display(result["group_summary"])
    print("\nReport preview:")
    _display(result["report"].head(20))
    return result


__all__ = [
    "annotation_columns",
    "create_presence_dataframe",
    "download_lfq_tables",
    "get_pride_files",
    "lfq_columns",
    "session_with_retries",
    "make_numbered_lfq_catalog",
    "prepare_lfq_selection",
    "resolve_lfq_ranges",
    "preview_lfq_range_groups",
    "build_range_representative_report",
    "process_lfq_ranges_and_show",
]

