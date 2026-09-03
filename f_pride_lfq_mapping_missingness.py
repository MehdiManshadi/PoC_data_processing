#!/usr/bin/env python3
"""Find PRIDE proteinGroups tables and report LFQ mappability + missingness.

Each output row represents one direct or ZIP-contained proteinGroups.txt file.
The gene matching and feature missingness logic comes from
``d_missingness_estimation.py``; the LFQ grouping logic comes from
``e_pride_lfq_mapping_checker.py``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import tempfile
import time
import urllib.parse
from collections import OrderedDict, defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from zipfile import ZipFile

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# Edit this small default list, pass accessions on the command line, or use
# --accessions-file for a long list.
PRIDE_ACCESSIONS = ['PXD000172',
'PXD000260',
'PXD000394',
'PXD000612',
'PXD000697',
'PXD000853',
'PXD000959',
'PXD000987',
'PXD001040',
'PXD001171',
'PXD001186',
'PXD001274',
'PXD001457',
'PXD001487',
'PXD001550',
'PXD001610',
'PXD001673',
'PXD001863',
'PXD002029',
'PXD002137',
'PXD002279',
'PXD002339',
'PXD002381',
'PXD002516',
'PXD002623',
'PXD002629',
'PXD002646',
'PXD002776',
'PXD002815',
'PXD002854',
'PXD003028',
'PXD003138',
'PXD003218',
'PXD003430',
'PXD003478',
'PXD003509',
'PXD003527',
'PXD003528',
'PXD003646',
'PXD003668',
'PXD003691',
'PXD003808',
]

PRIDE_FILES_API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
PRIDE_PROJECT_API = "https://www.ebi.ac.uk/pride/ws/archive/v2"
UNIPROT_API = "https://rest.uniprot.org"
REQUEST_TIMEOUT = 300
SMALL_TEXT_LIMIT = 20 * 1024**2

GENE_ALIASES = {
    "gene names", "gene name", "gene", "genes", "gene symbol",
    "gene symbols", "genesymbol", "genesymbols",
}
FASTA_ALIASES = {"fasta header", "fasta headers"}
MAJORITY_ALIASES = {"majority protein id", "majority protein ids"}
LFQ_RE = re.compile(
    r"^\s*LFQ(?:\s+|_)+intensity(?:\s+|_)*(.*?)\s*$",
    re.IGNORECASE,
)
UNIPROT_RE = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)"
)
MARKED_REPLICATE_RE = re.compile(
    r"(?i)^(.+?)(?:[\s_.-]*(?:bio(?:logical)?rep|tech(?:nical)?rep|"
    r"rep(?:licate)?|repl|experiment|exp|ex|run|r))[\s_.-]*0*(\d+)$"
)


@dataclass
class ProteinGroupsSource:
    proteinGroups_address: str
    local_file: str = ""
    archive_name: str = ""
    archive_member: str = ""
    download_error: str = ""


@dataclass
class MappingContext:
    documents: list[tuple[str, str]]
    raw_files: dict[str, list[str]]
    corpus: str
    warnings: list[str]


# ---------------------------------------------------------------------------
# HTTP and PRIDE discovery
# ---------------------------------------------------------------------------

def load_optional_dotenv() -> None:
    """Use a local .env when python-dotenv is installed; otherwise use os.environ."""
    try:
        from dotenv import load_dotenv
    except ImportError:
        return
    load_dotenv()


def http_session() -> requests.Session:
    retry = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods={"GET", "HEAD", "POST"},
    )
    session = requests.Session()
    session.headers.update({"User-Agent": "PRIDE-LFQ-mapping-missingness/1.0"})
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def api_files(payload) -> list[dict]:
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


def get_pride_files(accession: str, session: requests.Session) -> list[dict]:
    """Return every PRIDE file page, including datasets with more than 20 files."""
    files, seen = [], set()
    for page in range(1000):
        response = session.get(
            f"{PRIDE_FILES_API}/projects/{accession}/files",
            params={"page": page, "pageSize": 100},
            timeout=REQUEST_TIMEOUT,
        )
        response.raise_for_status()
        batch = api_files(response.json())
        if not batch:
            break

        added = 0
        for record in batch:
            name = str(record.get("fileName") or "").strip()
            url = safe_download_url(record)
            key = (name, url)
            if name and key not in seen:
                seen.add(key)
                files.append(record)
                added += 1
        if added == 0:
            break
    return files


def download_url(record: dict) -> str:
    locations = sorted(
        record.get("publicFileLocations") or [],
        key=lambda item: "FTP" in str(item.get("name", "")),
        reverse=True,
    )
    for location in locations:
        value = str(location.get("value") or "").strip()
        if value.startswith("ftp://ftp.pride.ebi.ac.uk/"):
            value = value.replace(
                "ftp://ftp.pride.ebi.ac.uk/",
                "https://ftp.pride.ebi.ac.uk/",
                1,
            )
        elif value.startswith("http://"):
            value = "https://" + value.removeprefix("http://")
        if value.startswith("https://"):
            return urllib.parse.quote(value, safe=":/@?&=%+;,$-_.!~*'()")
    raise ValueError(f"No HTTP/FTP URL for {record.get('fileName')}")


def safe_download_url(record: dict) -> str:
    try:
        return download_url(record)
    except ValueError:
        return ""


def file_size(record: dict) -> int:
    try:
        return int(record.get("fileSizeBytes") or 0)
    except (TypeError, ValueError):
        return 0


def project_metadata(accession: str, session: requests.Session) -> tuple[dict, str]:
    try:
        response = session.get(
            f"{PRIDE_PROJECT_API}/projects/{accession}",
            timeout=REQUEST_TIMEOUT,
        )
        response.raise_for_status()
        return response.json(), ""
    except Exception as exc:
        return {}, f"Project metadata unavailable: {exc}"


def normalized_name(path: str) -> str:
    return PurePosixPath(str(path).replace("\\", "/")).name


def normalized_text(value) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).casefold()).strip()


def compact_text(value: str) -> str:
    return "".join(re.findall(r"[a-z0-9]+", value.lower()))


def is_proteingroups_file(path: str) -> bool:
    return normalized_name(path).casefold() in {
        "proteingroups.txt",
        "proteingroups.txt.gz",
    }


def metadata_score(path: str) -> int:
    name = normalized_name(path).casefold()
    allowed = (".txt", ".tsv", ".csv", ".xml", ".json")
    if is_proteingroups_file(name) or not name.endswith(allowed):
        return 0
    words = [
        "sdrf", "experimentaldesign", "metadata", "manifest", "annotation",
        "design", "readme", "mqpar", "sample",
    ]
    return len(words) - next(
        (index for index, word in enumerate(words) if word in name),
        len(words),
    )


def decode_text(data: bytes) -> str:
    for encoding in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            pass
    return data.decode(errors="replace")


def read_url_prefix(session: requests.Session, url: str, limit: int) -> bytes:
    with session.get(
        url,
        headers={"Range": f"bytes=0-{limit - 1}", "Accept-Encoding": "identity"},
        stream=True,
        timeout=REQUEST_TIMEOUT,
    ) as response:
        response.raise_for_status()
        return response.raw.read(limit)


class RequestsRangeReader(io.RawIOBase):
    """Seekable HTTP range reader used by ZipFile without full ZIP download."""

    def __init__(
        self,
        session: requests.Session,
        url: str,
        size: int,
        chunk_size: int = 1024**2,
    ):
        self.session = session
        self.url = url
        self.size = size
        self.chunk_size = chunk_size
        self.position = 0
        self.cache: OrderedDict[int, bytes] = OrderedDict()

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.position

    def seek(self, offset: int, whence: int = os.SEEK_SET) -> int:
        if whence == os.SEEK_SET:
            position = offset
        elif whence == os.SEEK_CUR:
            position = self.position + offset
        elif whence == os.SEEK_END:
            position = self.size + offset
        else:
            raise ValueError(f"Unsupported whence: {whence}")
        self.position = max(0, position)
        return self.position

    def _block(self, index: int) -> bytes:
        if index in self.cache:
            self.cache.move_to_end(index)
            return self.cache[index]

        start = index * self.chunk_size
        end = min((index + 1) * self.chunk_size, self.size) - 1
        response = self.session.get(
            self.url,
            headers={"Range": f"bytes={start}-{end}", "Accept-Encoding": "identity"},
            timeout=REQUEST_TIMEOUT,
        )
        response.raise_for_status()
        if response.status_code != 206:
            raise OSError("Server does not support ZIP byte ranges")
        data = response.content
        if len(data) != end - start + 1:
            raise OSError("Incomplete ZIP byte range")

        self.cache[index] = data
        while len(self.cache) > 48:
            self.cache.popitem(last=False)
        return data

    def read(self, size: int = -1) -> bytes:
        if size is None or size < 0:
            size = self.size - self.position
        size = min(size, self.size - self.position)
        output = bytearray()
        while size > 0:
            block = self._block(self.position // self.chunk_size)
            offset = self.position % self.chunk_size
            take = min(size, len(block) - offset)
            if take <= 0:
                break
            output.extend(block[offset : offset + take])
            self.position += take
            size -= take
        return bytes(output)

    def readinto(self, buffer) -> int:
        data = self.read(len(buffer))
        buffer[: len(data)] = data
        return len(data)


def remote_file_size(
    session: requests.Session,
    url: str,
    size_hint: int = 0,
) -> int:
    if size_hint > 0:
        return size_hint
    with session.get(
        url,
        headers={"Range": "bytes=0-0", "Accept-Encoding": "identity"},
        stream=True,
        timeout=REQUEST_TIMEOUT,
    ) as response:
        response.raise_for_status()
        content_range = response.headers.get("Content-Range", "")
        if "/" in content_range:
            return int(content_range.rsplit("/", 1)[1])
        content_length = response.headers.get("Content-Length")
        if content_length:
            return int(content_length)
    raise OSError("Could not determine remote ZIP size")


@contextmanager
def remote_zip(
    session: requests.Session,
    url: str,
    size_hint: int = 0,
):
    size = remote_file_size(session, url, size_hint)
    archive = ZipFile(RequestsRangeReader(session, url, size))
    try:
        yield archive
    finally:
        archive.close()


def copy_stream(source, destination: Path) -> None:
    with destination.open("wb") as output:
        shutil.copyfileobj(source, output, length=4 * 1024**2)


def download_whole(
    session: requests.Session,
    url: str,
    destination: Path,
) -> None:
    with session.get(url, stream=True, timeout=REQUEST_TIMEOUT) as response:
        response.raise_for_status()
        with destination.open("wb") as output:
            for chunk in response.iter_content(4 * 1024**2):
                if chunk:
                    output.write(chunk)


def source_address(url: str, member: str = "") -> str:
    return f"{url}!/{member.lstrip('/')}" if member else url


def source_output_path(
    output_dir: Path,
    accession: str,
    address: str,
    archive_name: str,
    member: str = "",
) -> Path:
    if member:
        archive_label = re.sub(
            r"\.zip$", "", normalized_name(archive_name), flags=re.IGNORECASE
        )
        parent = str(PurePosixPath(member).parent)
        label = archive_label if parent == "." else f"{archive_label}_{parent}"
    else:
        label = "direct"
    label = re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("._-") or "source"
    digest = hashlib.sha256(address.encode("utf-8")).hexdigest()[:12]
    directory = output_dir / accession / "proteinGroups_files"
    directory.mkdir(parents=True, exist_ok=True)
    return directory / f"{accession}__{label[:80]}__{digest}.txt"


def extract_member(archive: ZipFile, member: str, destination: Path) -> None:
    with archive.open(member) as source:
        if member.casefold().endswith(".gz"):
            with gzip.GzipFile(fileobj=source) as decompressed:
                copy_stream(decompressed, destination)
        else:
            copy_stream(source, destination)


def inspect_archive(
    archive: ZipFile,
    accession: str,
    output_dir: Path,
    archive_name: str,
    url: str,
    already_seen: set[str],
    continue_on_error: bool = False,
) -> tuple[list[ProteinGroupsSource], list[tuple[str, str]]]:
    sources, documents, newly_seen = [], [], set()
    for info in archive.infolist():
        if info.is_dir():
            continue
        member = info.filename
        if is_proteingroups_file(member):
            address = source_address(url, member)
            if address in already_seen or address in newly_seen:
                continue
            destination = source_output_path(
                output_dir, accession, address, archive_name, member
            )
            try:
                extract_member(archive, member, destination)
                sources.append(
                    ProteinGroupsSource(
                        address,
                        str(destination),
                        archive_name,
                        member,
                    )
                )
            except Exception as exc:
                if not continue_on_error:
                    raise
                sources.append(
                    ProteinGroupsSource(
                        address,
                        archive_name=archive_name,
                        archive_member=member,
                        download_error=str(exc),
                    )
                )
            newly_seen.add(address)
        elif metadata_score(member) and info.file_size <= SMALL_TEXT_LIMIT:
            with archive.open(member) as source:
                content = source.read(SMALL_TEXT_LIMIT)
            documents.append((f"{archive_name}::{member}", decode_text(content)))
    already_seen.update(newly_seen)
    return sources, documents


def discover_accession(
    accession: str,
    output_dir: Path,
    session: requests.Session,
    max_full_zip_bytes: int,
) -> tuple[list[ProteinGroupsSource], MappingContext]:
    records = get_pride_files(accession, session)
    if not records:
        raise FileNotFoundError("PRIDE returned no files")

    sources: list[ProteinGroupsSource] = []
    documents: list[tuple[str, str]] = []
    warnings: list[str] = []
    seen_addresses: set[str] = set()

    for record in sorted(
        records,
        key=lambda item: -metadata_score(str(item.get("fileName") or "")),
    ):
        name = str(record.get("fileName") or "")
        url = safe_download_url(record)
        if not url:
            continue
        if is_proteingroups_file(name):
            address = source_address(url)
            if address in seen_addresses:
                continue
            destination = source_output_path(
                output_dir, accession, address, name
            )
            try:
                with session.get(
                    url, stream=True, timeout=REQUEST_TIMEOUT
                ) as response:
                    response.raise_for_status()
                    response.raw.decode_content = True
                    if name.casefold().endswith(".gz"):
                        with gzip.GzipFile(fileobj=response.raw) as decompressed:
                            copy_stream(decompressed, destination)
                    else:
                        copy_stream(response.raw, destination)
                sources.append(ProteinGroupsSource(address, str(destination)))
            except Exception as exc:
                sources.append(
                    ProteinGroupsSource(address, download_error=str(exc))
                )
            seen_addresses.add(address)
        elif metadata_score(name) and (
            file_size(record) == 0 or file_size(record) <= SMALL_TEXT_LIMIT
        ):
            try:
                documents.append(
                    (name, decode_text(read_url_prefix(session, url, SMALL_TEXT_LIMIT)))
                )
            except Exception as exc:
                warnings.append(f"Could not read metadata file {name}: {exc}")

    for record in records:
        archive_name = str(record.get("fileName") or "")
        if not normalized_name(archive_name).casefold().endswith(".zip"):
            continue
        url = safe_download_url(record)
        if not url:
            continue

        remote_error = None
        try:
            with remote_zip(session, url, file_size(record)) as archive:
                archive_sources, archive_documents = inspect_archive(
                    archive,
                    accession,
                    output_dir,
                    archive_name,
                    url,
                    seen_addresses,
                )
            sources.extend(archive_sources)
            documents.extend(archive_documents)
            continue
        except Exception as exc:
            remote_error = exc

        size = file_size(record)
        if size and size > max_full_zip_bytes:
            warnings.append(
                f"Skipped ZIP {archive_name}: byte-range inspection failed and "
                f"the archive exceeds the full-download limit ({remote_error})"
            )
            continue

        with tempfile.TemporaryDirectory() as temporary_directory:
            local_zip = Path(temporary_directory) / "archive.zip"
            try:
                download_whole(session, url, local_zip)
                with ZipFile(local_zip) as archive:
                    archive_sources, archive_documents = inspect_archive(
                        archive,
                        accession,
                        output_dir,
                        archive_name,
                        url,
                        seen_addresses,
                        continue_on_error=True,
                    )
                sources.extend(archive_sources)
                documents.extend(archive_documents)
            except Exception as exc:
                warnings.append(f"Could not inspect ZIP {archive_name}: {exc}")

    metadata, metadata_warning = project_metadata(accession, session)
    if metadata_warning:
        warnings.append(metadata_warning)
    raw_files = {
        compact_text(str(record.get("fileName") or "")): [
            str(record.get("fileName") or "")
        ]
        for record in records
        if normalized_name(str(record.get("fileName") or "")).casefold().endswith(
            (".raw", ".wiff", ".mzml", ".mzxml", ".d")
        )
    }
    corpus = " ".join(
        re.findall(
            r"[a-z0-9]+",
            (json.dumps(metadata) + " " + " ".join(x[1] for x in documents)).lower(),
        )
    )
    context = MappingContext(documents, raw_files, corpus, warnings)
    return sources, context


# ---------------------------------------------------------------------------
# LFQ mappability (same decision logic as e_pride_lfq_mapping_checker.py)
# ---------------------------------------------------------------------------

def lfq_columns(columns) -> list[tuple[str, str]]:
    return [
        (str(column), match.group(1).strip())
        for column in columns
        if (match := LFQ_RE.match(str(column)))
    ]


def useful_group(value: str) -> bool:
    key = compact_text(value)
    return len(key) >= 2 and not re.fullmatch(r"[a-z]\d*", key)


def remove_run_prefix(labels: list[str]) -> dict[str, str]:
    rows = [re.split(r"[\s_.-]+", label.strip(" _-")) for label in labels]
    if not rows or any(len(row) < 2 for row in rows):
        return dict(zip(labels, labels))
    first = [row[0] for row in rows]
    unique = len({x.lower() for x in first}) / len(first) >= 0.8
    numbered = all(re.search(r"\d", x) for x in first)
    skeletons = {re.sub(r"\d+", "", x.lower()) for x in first} - {""}
    remove = unique and numbered and len(skeletons) <= max(1, len(first) // 5)
    return {
        label: "_".join(row[1:] if remove else row)
        for label, row in zip(labels, rows)
    }


def split_rule(value: str, rule: str, width: int = 0):
    value = value.strip(" _-")
    if rule == "marked":
        match = MARKED_REPLICATE_RE.fullmatch(value)
        return (
            (match.group(1).strip(" _-"), str(int(match.group(2))))
            if match
            else None
        )
    if rule == "suffix_token":
        match = re.fullmatch(r"(.+?)[\s_.-]+([^\s_.-]+)", value)
        return (match.group(1), match.group(2)) if match else None
    if rule == "prefix_token":
        match = re.fullmatch(r"([^\s_.-]+)[\s_.-]+(.+)", value)
        return (match.group(2), match.group(1)) if match else None
    if len(value) <= width:
        return None
    if rule == "suffix_chars":
        return value[:-width].strip(" _-"), value[-width:].strip(" _-")
    return value[width:].strip(" _-"), value[:width].strip(" _-")


def grid_score(splits: dict, bonus: int):
    groups = defaultdict(list)
    for _, (group, member) in splits.items():
        group_key, member_key = compact_text(group), compact_text(member)
        if (
            not group_key
            or not member_key
            or group_key == member_key
            or len(member_key) > 8
        ):
            return None
        groups[group_key].append(member_key)
    if len(groups) < 2 or any(
        len(items) < 2 or len(items) != len(set(items))
        for items in groups.values()
    ):
        return None
    member_sets = [set(items) for items in groups.values()]
    overlaps = [
        len(left & right) / len(left | right)
        for index, left in enumerate(member_sets)
        for right in member_sets[index + 1 :]
    ]
    if not overlaps or min(overlaps) < 0.75:
        return None
    sizes = [len(items) for items in member_sets]
    if min(sizes) / max(sizes) < 0.5:
        return None
    return (
        50 * sum(overlaps) / len(overlaps)
        + 15 * min(sizes) / max(sizes)
        + 2 * min(len(groups), 6)
        + bonus
    )


def infer_grid(labels: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    values = remove_run_prefix(labels)
    rules = [("marked", 0, 40), ("suffix_token", 0, 30), ("prefix_token", 0, 15)]
    rules += [("suffix_chars", width, 10) for width in range(1, 7)]
    rules += [("prefix_chars", width, 5) for width in range(1, 7)]
    best = None
    for rule, width, bonus in rules:
        splits = {
            label: split_rule(value, rule, width)
            for label, value in values.items()
        }
        if any(value is None or not all(value) for value in splits.values()):
            continue
        score = grid_score(splits, bonus)
        if score is not None and (best is None or score > best[0]):
            best = score, splits
    if best is None:
        return {}, {}
    return (
        {label: pair[0] for label, pair in best[1].items()},
        {label: pair[1] for label, pair in best[1].items()},
    )


def structured_map(documents: list[tuple[str, str]]) -> dict[str, list[tuple[str, str]]]:
    output = defaultdict(list)
    for name, content in documents:
        lines = [line for line in content.splitlines() if line.strip()][:5001]
        if len(lines) < 2:
            continue
        delimiter = max(("\t", ",", ";"), key=lambda item: lines[0].count(item))
        rows = list(csv.reader(lines, delimiter=delimiter))
        heads = [re.sub(r"\W", "", item.lower()) for item in rows[0]]
        keys = [
            index
            for index, value in enumerate(heads)
            if any(
                word in value
                for word in (
                    "datafile", "rawfile", "filename", "assayname", "samplename"
                )
            )
            or value == "name"
        ]
        factors = [
            index
            for index, value in enumerate(heads)
            if any(
                word in value
                for word in (
                    "factorvalue", "condition", "group", "treatment",
                    "experiment", "phenotype", "genotype", "dose", "time",
                )
            )
        ]
        if not keys or not factors:
            continue
        width = len(rows[0])
        for row in rows[1:]:
            row += [""] * (width - len(row))
            values = [
                f"{rows[0][index]}={row[index]}"
                for index in factors
                if row[index].strip()
            ]
            for index in keys:
                if compact_text(row[index]) and values:
                    output[compact_text(row[index])].append(
                        ("; ".join(values), name)
                    )
    return output


def unique_match(label: str, mapping: dict) -> list:
    key = compact_text(label)
    if key in mapping:
        return mapping[key]
    matches = [
        value
        for candidate, value in mapping.items()
        if min(len(key), len(candidate)) >= 5
        and (key in candidate or candidate in key)
    ]
    return matches[0] if len(matches) == 1 else []


def informative(condition: str, label: str) -> bool:
    biological = (
        "condition", "group", "treatment", "factor", "phenotype",
        "genotype", "dose", "time",
    )
    for part in condition.split(";"):
        head, _, value = part.partition("=")
        if any(word in compact_text(head) for word in biological):
            return True
        if compact_text(value) != compact_text(label):
            return True
    return False


def assess_mappability(columns, context: MappingContext) -> dict:
    pairs = lfq_columns(columns)
    labels = [label for _, label in pairs]
    full_columns = {label: column for column, label in pairs}
    grid_groups, grid_members = infer_grid(labels)
    design = structured_map(context.documents)

    groups, members, raw_matches = {}, {}, {}
    mapped = 0
    for label in labels:
        hits = [
            hit
            for hit in unique_match(label, design)
            if informative(hit[0], label)
        ]
        raw_hit = unique_match(label, context.raw_files)
        if raw_hit:
            raw_matches[label] = raw_hit[0]
        if len({hit[0] for hit in hits}) == 1:
            groups[label] = hits[0][0]
            mapped += 1
        elif label in grid_groups:
            groups[label] = grid_groups[label]
            members[label] = grid_members[label]
            mapped += 1
        elif useful_group(label) and " ".join(
            re.findall(r"[a-z0-9]+", label.lower())
        ) in context.corpus:
            groups[label] = label
            mapped += 1

    total = len(labels)
    status = (
        "yes"
        if total and mapped == total
        else "partial"
        if mapped
        else "uncertain"
        if raw_matches
        else "no"
    )
    grouped_columns = defaultdict(list)
    for label, group in groups.items():
        grouped_columns[group].append(full_columns.get(label, label))

    return {
        "lfq_mapping_status": status,
        "lfq_mapping_reason": (
            f"{mapped}/{total} LFQ names mapped to groups"
            if total
            else "No LFQ intensity columns found"
        ),
        "total_lfq_columns": total,
        "mapped_lfq_columns": mapped,
        "mapped_lfq_percent": 100 * mapped / total if total else 0.0,
        "lfq_column_names": json.dumps([column for column, _ in pairs]),
        "lfq_sample_names": json.dumps(labels),
        "grouped_lfq_columns": json.dumps(dict(grouped_columns), sort_keys=True),
        "lfq_groups": json.dumps(groups, sort_keys=True),
        "lfq_replicates": json.dumps(members, sort_keys=True),
        "raw_file_matches": json.dumps(raw_matches, sort_keys=True),
    }


# ---------------------------------------------------------------------------
# Gene extraction and feature missingness (same logic as d_missingness_estimation.py)
# ---------------------------------------------------------------------------

def gene_sources(columns):
    normalized = {normalized_text(column): column for column in columns}
    gene_column = next(
        (normalized[alias] for alias in GENE_ALIASES if alias in normalized),
        None,
    )
    if gene_column is None:
        excluded = {
            "gene ontology", "gene ontology id", "gene ontology ids",
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
    lfq = [column for column, _ in lfq_columns(columns)]
    return gene_column, fasta_columns, majority_column, lfq


def fasta_genes(value) -> list[str]:
    if pd.isna(value):
        return []
    return re.findall(r"\bGN=([^\s;]+)", str(value), flags=re.IGNORECASE)


def uniprot_ids(value) -> list[str]:
    if pd.isna(value):
        return []
    result = []
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


def map_uniprot(
    accessions: list[str],
    session: requests.Session,
    cache: dict[str, list[str]],
) -> dict[str, list[str]]:
    pending = [item for item in dict.fromkeys(accessions) if item not in cache]
    for start in range(0, len(pending), 10000):
        batch = pending[start : start + 10000]
        response = session.post(
            f"{UNIPROT_API}/idmapping/run",
            data={
                "from": "UniProtKB_AC-ID",
                "to": "Gene_Name",
                "ids": ",".join(batch),
            },
            timeout=REQUEST_TIMEOUT,
        )
        response.raise_for_status()
        job = response.json()["jobId"]
        while True:
            status_response = session.get(
                f"{UNIPROT_API}/idmapping/status/{job}",
                timeout=REQUEST_TIMEOUT,
            )
            status_response.raise_for_status()
            status = status_response.json().get("jobStatus")
            if status in {"NEW", "RUNNING"}:
                time.sleep(3)
                continue
            if status == "FAILED":
                raise RuntimeError("UniProt mapping failed")
            break

        response = session.get(
            f"{UNIPROT_API}/idmapping/stream/{job}",
            params={"format": "json"},
            timeout=REQUEST_TIMEOUT,
        )
        response.raise_for_status()
        for accession in batch:
            cache.setdefault(accession, [])
        for item in response.json().get("results", []):
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
    file_path: Path,
    columns,
    session: requests.Session,
    uniprot_cache: dict[str, list[str]],
) -> tuple[pd.DataFrame, list[str]]:
    gene_column, fasta_columns, majority_column, lfq = gene_sources(columns)
    if not lfq:
        raise ValueError("No LFQ intensity columns found")
    if not (gene_column or fasta_columns or majority_column):
        raise ValueError("No usable gene source found")

    use_columns = list(
        dict.fromkeys(
            ([gene_column] if gene_column else [])
            + fasta_columns
            + ([majority_column] if majority_column else [])
            + lfq
        )
    )
    frame = pd.read_csv(file_path, sep="\t", usecols=use_columns, low_memory=False)
    genes_output, sources_output = [], []
    unresolved_rows, all_uniprot = {}, []

    for row_index, row in frame.iterrows():
        genes, seen, matching_sources = [], set(), []
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
            map_uniprot(all_uniprot, session, uniprot_cache)
            for row_index, ids in unresolved_rows.items():
                genes = []
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
            print(f"    UniProt fallback failed: {exc}")

    numeric_lfq = frame[lfq].apply(pd.to_numeric, errors="coerce")
    present = numeric_lfq.notna() & numeric_lfq.ne(0)
    frame["Presence_count"] = present.sum(axis=1)
    frame["Total_LFQ_columns"] = len(lfq)
    frame["Presence"] = frame["Presence_count"] / len(lfq)
    frame["Presence_percent"] = frame["Presence"] * 100

    output_columns = (
        ["Gene names", "Gene matching source"]
        + lfq
        + ["Presence_count", "Total_LFQ_columns", "Presence", "Presence_percent"]
    )
    return frame[output_columns], lfq


def read_feature_genes(path: Path) -> list[str]:
    genes = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        for gene in re.split(r"[;,\t]+", line):
            gene = gene.strip()
            if gene and gene.casefold() not in {
                "gene", "genes", "gene name", "gene names", "feature",
                "features", "gene symbol", "gene symbols",
            }:
                genes.append(gene.upper())
    unique = list(dict.fromkeys(genes))
    if not unique:
        raise ValueError(f"No feature genes found in {path}")
    return unique


def available_genes(frame: pd.DataFrame) -> set[str]:
    return set(
        frame["Gene names"]
        .dropna()
        .str.split(";")
        .explode()
        .str.strip()
        .str.upper()
    ) - {""}


def feature_missingness(frame: pd.DataFrame, feature_genes: list[str]) -> dict:
    present = available_genes(frame)
    target = {gene.upper() for gene in feature_genes}
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


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

REPORT_COLUMNS = [
    "PRIDE_accession",
    "proteinGroups_address",
    "analysis_status",
    "analysis_error",
    "lfq_mapping_status",
    "lfq_mapping_reason",
    "hyperglycemia_missingness_percent",
    "mitochondrial_myopathy_missingness_percent",
    "total_protein_groups",
    "protein_groups_with_gene_names",
    "number_of_available_gene_names",
    "total_lfq_columns",
    "mapped_lfq_columns",
    "mapped_lfq_percent",
    "lfq_column_names",
    "lfq_sample_names",
    "grouped_lfq_columns",
    "lfq_groups",
    "lfq_replicates",
    "raw_file_matches",
    "available_gene_names",
    "hyperglycemia_total_genes",
    "hyperglycemia_found_genes",
    "hyperglycemia_missing_genes",
    "hyperglycemia_found_gene_names",
    "hyperglycemia_missing_gene_names",
    "mitochondrial_myopathy_total_genes",
    "mitochondrial_myopathy_found_genes",
    "mitochondrial_myopathy_missing_genes",
    "mitochondrial_myopathy_found_gene_names",
    "mitochondrial_myopathy_missing_gene_names",
    "project_link",
    "downloaded_file",
    "source_archive",
    "source_member",
    "processed_file",
    "mapping_context_warnings",
]


def prefixed(prefix: str, values: dict) -> dict:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def analyze_source(
    accession: str,
    source: ProteinGroupsSource,
    context: MappingContext,
    output_dir: Path,
    hyperglycemia_genes: list[str],
    mitochondrial_genes: list[str],
    session: requests.Session,
    uniprot_cache: dict[str, list[str]],
    save_processed: bool,
) -> tuple[dict, pd.DataFrame | None]:
    row = {
        "PRIDE_accession": accession,
        "proteinGroups_address": source.proteinGroups_address,
        "analysis_status": "failed",
        "analysis_error": source.download_error,
        "project_link": f"https://www.ebi.ac.uk/pride/archive/projects/{accession}",
        "downloaded_file": source.local_file,
        "source_archive": source.archive_name,
        "source_member": source.archive_member,
        "processed_file": "",
        "mapping_context_warnings": " | ".join(context.warnings),
    }
    if source.download_error:
        return row, None

    try:
        columns = list(pd.read_csv(source.local_file, sep="\t", nrows=0).columns)
        row.update(assess_mappability(columns, context))
    except Exception as exc:
        row["analysis_error"] = f"Could not read table header: {exc}"
        return row, None

    try:
        frame, _ = create_presence_dataframe(
            Path(source.local_file), columns, session, uniprot_cache
        )
        if save_processed:
            processed_directory = output_dir / accession / "processed_tables"
            processed_directory.mkdir(parents=True, exist_ok=True)
            processed_file = (
                processed_directory / f"{Path(source.local_file).stem}_processed.txt"
            )
            frame.to_csv(processed_file, sep="\t", index=False)
            row["processed_file"] = str(processed_file)

        hyper = feature_missingness(frame, hyperglycemia_genes)
        mitochondrial = feature_missingness(frame, mitochondrial_genes)
        genes = sorted(available_genes(frame))
        row.update(
            {
                "analysis_status": "success",
                "analysis_error": "",
                "total_protein_groups": len(frame),
                "protein_groups_with_gene_names": int(
                    frame["Gene names"].notna().sum()
                ),
                "number_of_available_gene_names": len(genes),
                "available_gene_names": ";".join(genes),
                **prefixed("hyperglycemia", hyper),
                **prefixed("mitochondrial_myopathy", mitochondrial),
            }
        )
        return row, frame
    except Exception as exc:
        row["analysis_error"] = str(exc)
        return row, None


def parse_accessions_file(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8-sig")
    return re.findall(r"PXD\d+", text, flags=re.IGNORECASE)


def unique_accessions(values: list[str]) -> list[str]:
    return list(
        dict.fromkeys(
            value.strip().upper()
            for value in values
            if re.fullmatch(r"PXD\d+", value.strip(), flags=re.IGNORECASE)
        )
    )


def run_pipeline(
    accessions: list[str],
    output_dir: Path,
    hyperglycemia_file: Path,
    mitochondrial_file: Path,
    max_full_zip_gb: float = 20,
    save_processed: bool = True,
) -> tuple[dict, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    hyperglycemia_genes = read_feature_genes(hyperglycemia_file)
    mitochondrial_genes = read_feature_genes(mitochondrial_file)
    session = http_session()
    uniprot_cache: dict[str, list[str]] = {}
    dataframes: dict[str, dict[str, pd.DataFrame]] = {}
    rows = []

    try:
        for number, accession in enumerate(accessions, 1):
            print(f"[{number}/{len(accessions)}] {accession}")
            try:
                sources, context = discover_accession(
                    accession,
                    output_dir,
                    session,
                    int(max_full_zip_gb * 1024**3),
                )
            except Exception as exc:
                print(f"  Discovery failed: {exc}")
                continue
            if not sources:
                print("  No proteinGroups.txt files found")
                continue

            print(f"  Found {len(sources)} unique proteinGroups table(s)")
            dataframes.setdefault(accession, {})
            for source_number, source in enumerate(sources, 1):
                row, frame = analyze_source(
                    accession,
                    source,
                    context,
                    output_dir,
                    hyperglycemia_genes,
                    mitochondrial_genes,
                    session,
                    uniprot_cache,
                    save_processed,
                )
                rows.append(row)
                if frame is not None:
                    dataframes[accession][source.proteinGroups_address] = frame
                mapping = row.get("lfq_mapping_status", "not assessed")
                hyper = row.get("hyperglycemia_missingness_percent")
                mito = row.get("mitochondrial_myopathy_missingness_percent")
                scores = (
                    f"hyper={hyper:.2f}%, mito={mito:.2f}%"
                    if hyper is not None and mito is not None
                    else row.get("analysis_error", "analysis failed")
                )
                print(
                    f"  [{source_number}/{len(sources)}] "
                    f"mapping={mapping}; {scores}"
                )
    finally:
        session.close()

    report = pd.DataFrame(rows).reindex(columns=REPORT_COLUMNS)
    if not report.empty and report.duplicated(
        subset=["PRIDE_accession", "proteinGroups_address"], keep=False
    ).any():
        raise RuntimeError("Duplicate accession/proteinGroups address pairs in report")

    report_path = output_dir / "dataset_missingness_mapping_report.csv"
    report.to_csv(report_path, index=False)
    print(f"Saved report: {report_path.resolve()}")
    return dataframes, report


def self_test() -> None:
    coded = [
        f"{group}_{replicate}"
        for group in (
            "F1", "F2", "F3", "F4", "F5", "M1", "M2", "M3", "M4",
            "M5", "Top6_Top20",
        )
        for replicate in ("1", "2", "3")
    ]
    groups, members = infer_grid(coded)
    assert len(groups) == 33
    assert set(groups.values()) == {
        "F1", "F2", "F3", "F4", "F5", "M1", "M2", "M3", "M4", "M5",
        "Top6_Top20",
    }
    assert set(members.values()) == {"1", "2", "3"}
    assert infer_grid(["Alpha", "Beta", "Gamma", "Delta"]) == ({}, {})

    context = MappingContext([], {}, "", [])
    mapping = assess_mappability(
        [
            "Gene names",
            "LFQ intensity Control_1",
            "LFQ intensity Control_2",
            "LFQ intensity Drug_1",
            "LFQ intensity Drug_2",
        ],
        context,
    )
    assert mapping["lfq_mapping_status"] == "yes"
    assert mapping["mapped_lfq_columns"] == 4

    frame = pd.DataFrame({"Gene names": ["A;B", "D", pd.NA]})
    missingness = feature_missingness(frame, ["A", "C"])
    assert missingness["missingness_percent"] == 50
    assert missingness["found_gene_names"] == "A"
    assert missingness["missing_gene_names"] == "C"

    with tempfile.TemporaryDirectory() as temporary_directory:
        table_path = Path(temporary_directory) / "proteinGroups.txt"
        pd.DataFrame(
            {
                "Gene names": ["A;B", "D"],
                "LFQ intensity Control_1": [10, 0],
                "LFQ intensity Control_2": [11, 3],
                "LFQ intensity Drug_1": [5, 7],
                "LFQ intensity Drug_2": [6, 8],
            }
        ).to_csv(table_path, sep="\t", index=False)
        source = ProteinGroupsSource("test://proteinGroups.txt", str(table_path))
        row, processed = analyze_source(
            "PXD000000",
            source,
            context,
            Path(temporary_directory),
            ["A", "C"],
            ["B", "D"],
            None,
            {},
            False,
        )
        assert row["analysis_status"] == "success"
        assert row["lfq_mapping_status"] == "yes"
        assert row["hyperglycemia_missingness_percent"] == 50
        assert row["mitochondrial_myopathy_missingness_percent"] == 0
        assert row["total_protein_groups"] == 2
        assert processed is not None
    print("Self-test passed")


def build_parser() -> argparse.ArgumentParser:
    output_dir = Path(os.getenv("OUTPUT_DIR", "PRIDE_downloads"))
    hyperglycemia_file = os.getenv("HYPERGLYCEMIA_FEATURES_FILE")
    mitochondrial_file = os.getenv("MITOCHONDRIAL_MYOPATHY_FEATURES_FILE")
    parser = argparse.ArgumentParser(
        description=(
            "Create one report row per PRIDE proteinGroups.txt with LFQ "
            "mappability and two feature-set missingness scores."
        )
    )
    parser.add_argument("accessions", nargs="*", help="PRIDE accessions")
    parser.add_argument(
        "--accessions-file",
        type=Path,
        help="Text/CSV/Python file containing PXD accessions",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=output_dir
    )
    parser.add_argument(
        "--hyperglycemia-features",
        type=Path,
        default=Path(hyperglycemia_file) if hyperglycemia_file else None,
    )
    parser.add_argument(
        "--mitochondrial-features",
        type=Path,
        default=(
            Path(mitochondrial_file)
            if mitochondrial_file
            else None
        ),
    )
    parser.add_argument("--max-full-zip-gb", type=float, default=20)
    parser.add_argument("--no-save-processed", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser


def main():
    load_optional_dotenv()
    parser = build_parser()
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return None, None

    supplied = list(args.accessions)
    if args.accessions_file:
        supplied.extend(parse_accessions_file(args.accessions_file))
    accessions = unique_accessions(supplied or PRIDE_ACCESSIONS)
    if not accessions:
        parser.error("Supply PRIDE accessions or edit PRIDE_ACCESSIONS")
    if args.hyperglycemia_features is None:
        parser.error(
            "Provide --hyperglycemia-features or HYPERGLYCEMIA_FEATURES_FILE"
        )
    if args.mitochondrial_features is None:
        parser.error(
            "Provide --mitochondrial-features or "
            "MITOCHONDRIAL_MYOPATHY_FEATURES_FILE"
        )

    return run_pipeline(
        accessions=accessions,
        output_dir=args.output_dir,
        hyperglycemia_file=args.hyperglycemia_features,
        mitochondrial_file=args.mitochondrial_features,
        max_full_zip_gb=args.max_full_zip_gb,
        save_processed=not args.no_save_processed,
    )


if __name__ == "__main__":
    dataframes, dataset_missingness_mapping_report = main()
