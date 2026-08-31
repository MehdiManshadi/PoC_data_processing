#!/usr/bin/env python3
"""Find human or zebrafish PRIDE datasets with an exact proteinGroups.txt file.

Pipeline
--------
1. Discover human or zebrafish proteomics records through OmicsDI.
2. Keep only PRIDE datasets whose public PRIDE manifest contains the exact
   canonical MaxQuant basename ``proteinGroups.txt`` (case-sensitive), either
   as a direct file or inside a public ZIP archive. ZIP files are inspected
   remotely through their central directory; they are never downloaded in full.
3. Enrich each retained dataset with structured OmicsDI, PRIDE, and SDRF
   metadata.
4. Retrieve a related publication abstract from Europe PMC when a PMID or DOI
   is available.
5. Text-mine tissue/biofluid names and explicit sample/participant/replicate
   counts from repository text and the publication abstract.
6. Read only the beginning of proteinGroups.txt (also when it is stored in a
   ZIP) to check whether actual ``LFQ intensity ...`` columns exist and count
   the quantitative profiles.
7. Conservatively compare LFQ sample labels with PRIDE manifest file basenames
   and report whether all, some, or none are directly mappable.

The CSV keeps evidence and provenance separate. A filename match proves that
the public PRIDE manifest exposes proteinGroups.txt; it does not by itself
prove that the table contains LFQ columns or reveal the biological sample
count. Those questions are reported separately.

Requires Python 3.10+ and uses only the standard library.

Quick test (targeted discovery; stops after 5 verified datasets)::

    python omicsdi_pride_proteingroups.py --test --test-limit 5 \
        --output human_proteingroups_test.csv

Zebrafish LFQ discovery::

    python omicsdi_pride_proteingroups.py --organism zebrafish \
        --require-lfq-columns --output zebrafish_lfq_report.csv

Full discovery::

    python omicsdi_pride_proteingroups.py \
        --max-candidates 5000 \
        --output human_proteingroups_report.csv

Strict LFQ-only output::

    python omicsdi_pride_proteingroups.py \
        --require-lfq-columns \
        --output human_lfq_proteingroups_report.csv
"""

from __future__ import annotations

import argparse
import csv
import html
import io
import json
import re
import sys
import struct
import time
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path, PurePosixPath
from typing import Any, Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import quote, unquote, urlencode
from urllib.request import Request, urlopen


OMID_BASE = "https://www.omicsdi.org/ws"
OMID_SEARCH = f"{OMID_BASE}/dataset/search"
PRIDE_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v3"
EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

TARGET_FILE = "proteinGroups.txt"
ZIP_EOCD_SIGNATURE = b"PK\x05\x06"
ZIP_CENTRAL_SIGNATURE = b"PK\x01\x02"
ZIP_LOCAL_SIGNATURE = b"PK\x03\x04"
# Free-text organism queries are used for discovery because taxonomy filters
# are not consistently populated in OmicsDI. Every retained dataset is then
# reconfirmed from structured OmicsDI, PRIDE, or SDRF metadata.
ORGANISM_CONFIG: dict[str, dict[str, str]] = {
    "human": {
        "label": "Homo sapiens (human)",
        "query": 'human AND omics_type:"Proteomics" AND repository:"PRIDE"',
    },
    "zebrafish": {
        "label": "Danio rerio (zebrafish)",
        "query": '(zebrafish OR "Danio rerio" OR "NCBITaxon:7955") '
        'AND omics_type:"Proteomics" AND repository:"PRIDE"',
    },
}
USER_AGENT = "pride-proteingroups-discovery/1.1"

CSV_FIELDS = [
    "accession",
    "title",
    "publication_date",
    "target_organism",
    "organism_verified",
    "human_verified",
    "organisms",
    "protein_groups_verified",
    "protein_groups_location",
    "protein_groups_archive_name",
    "protein_groups_file_name",
    "protein_groups_file_category",
    "protein_groups_file_size_mb",
    "protein_groups_url",
    "lfq_columns_detected",
    "lfq_column_count",
    "lfq_sample_names",
    "lfq_names_mappable_to_pride_files",
    "lfq_names_mapped_count",
    "lfq_names_unmapped_count",
    "lfq_names_unmapped",
    "lfq_name_file_matches",
    "ordinary_intensity_column_count",
    "pride_structured_tissues",
    "omicsdi_structured_tissues",
    "sdrf_tissues",
    "text_mined_tissues",
    "combined_tissues",
    "tissue_count",
    "sdrf_row_count",
    "sdrf_unique_source_count",
    "sdrf_unique_sample_count",
    "sample_count_estimate",
    "sample_count_basis",
    "sample_count_confidence",
    "sample_count_hints",
    "replicate_count_hints",
    "sample_evidence_pride",
    "sample_evidence_abstract",
    "tissue_evidence_sources",
    "pmid",
    "doi",
    "publication_title",
    "publication_url",
    "pride_dataset_url",
    "omicsdi_dataset_url",
]


# Canonical tissue/biofluid name -> expressions commonly used in metadata.
# Longer, specific terms are intentionally retained alongside broad coverage.
TISSUE_TERMS: dict[str, tuple[str, ...]] = {
    "adipose tissue": ("adipose tissue", "adipose", "fat tissue"),
    "visceral adipose tissue": ("visceral adipose tissue", "visceral fat", "vat"),
    "subcutaneous adipose tissue": (
        "subcutaneous adipose tissue",
        "subcutaneous fat",
        "sat",
    ),
    "skeletal muscle": ("skeletal muscle", "musculus skeletal muscle"),
    "cardiac muscle": ("cardiac muscle", "myocardium", "myocardial tissue"),
    "muscle": ("muscle tissue", "muscle"),
    "serum": ("serum",),
    "plasma": ("blood plasma", "plasma"),
    "whole blood": ("whole blood",),
    "blood": ("peripheral blood", "blood"),
    "pbmc": ("peripheral blood mononuclear cells", "pbmcs", "pbmc"),
    "bone marrow": ("bone marrow",),
    "cerebrospinal fluid": ("cerebrospinal fluid", "csf"),
    "urine": ("urine", "urinary proteome"),
    "saliva": ("saliva", "salivary fluid"),
    "breast": ("breast tissue", "mammary gland", "breast"),
    "liver": ("liver tissue", "liver", "hepatic tissue"),
    "brain": ("brain tissue", "brain"),
    "cerebral cortex": ("cerebral cortex", "brain cortex"),
    "hippocampus": ("hippocampus", "hippocampal tissue"),
    "cerebellum": ("cerebellum", "cerebellar tissue"),
    "heart": ("heart tissue", "heart"),
    "kidney": ("kidney tissue", "kidney", "renal tissue"),
    "lung": ("lung tissue", "lung", "pulmonary tissue"),
    "pancreas": ("pancreatic tissue", "pancreas"),
    "pancreatic islet": ("pancreatic islets", "pancreatic islet", "islets of langerhans"),
    "colon": ("colon tissue", "colon", "colonic tissue"),
    "rectum": ("rectal tissue", "rectum"),
    "small intestine": ("small intestine", "small bowel"),
    "duodenum": ("duodenum", "duodenal tissue"),
    "jejunum": ("jejunum", "jejunal tissue"),
    "ileum": ("ileum", "ileal tissue"),
    "stomach": ("stomach tissue", "stomach", "gastric tissue"),
    "esophagus": ("esophagus", "oesophagus", "esophageal tissue", "oesophageal tissue"),
    "prostate": ("prostate tissue", "prostate", "prostatic tissue"),
    "ovary": ("ovarian tissue", "ovary"),
    "testis": ("testicular tissue", "testis", "testes"),
    "placenta": ("placental tissue", "placenta"),
    "spleen": ("splenic tissue", "spleen"),
    "thyroid": ("thyroid gland", "thyroid tissue", "thyroid"),
    "skin": ("skin tissue", "skin", "dermis", "epidermis"),
    "bone": ("bone tissue", "bone"),
    "cartilage": ("cartilage", "cartilaginous tissue"),
    "retina": ("retinal tissue", "retina"),
    "eye": ("ocular tissue", "eye tissue"),
    "artery": ("arterial tissue", "artery", "aorta"),
    "vein": ("venous tissue", "vein"),
    "lymph node": ("lymph node", "lymph nodes"),
    "tonsil": ("tonsillar tissue", "tonsil"),
    "endometrium": ("endometrial tissue", "endometrium"),
    "uterus": ("uterine tissue", "uterus"),
    "bladder": ("urinary bladder", "bladder tissue", "bladder"),
    "gallbladder": ("gallbladder", "gall bladder"),
    "oral mucosa": ("oral mucosa", "buccal mucosa"),
    "nasal tissue": ("nasal tissue", "nasal epithelium", "nasal mucosa"),
    "bronchoalveolar lavage": ("bronchoalveolar lavage", "bal fluid", "balf"),
    "synovial fluid": ("synovial fluid",),
    "seminal plasma": ("seminal plasma", "seminal fluid"),
    "amniotic fluid": ("amniotic fluid",),
    "cell line": ("cell line", "cell lines"),
    "primary cells": ("primary cells", "primary cell"),
    "organoid": ("organoids", "organoid"),
}

SUPPRESS_BROAD_TISSUES = {
    "skeletal muscle": {"muscle"},
    "cardiac muscle": {"muscle", "heart"},
    "visceral adipose tissue": {"adipose tissue"},
    "subcutaneous adipose tissue": {"adipose tissue"},
    "whole blood": {"blood"},
    "plasma": {"blood"},
    "pbmc": {"blood"},
    "pancreatic islet": {"pancreas"},
    "cerebral cortex": {"brain"},
    "hippocampus": {"brain"},
    "cerebellum": {"brain"},
    "cerebrospinal fluid": {"brain"},
    "seminal plasma": {"plasma"},
}

NUMBER_WORDS = {
    "one": 1,
    "two": 2,
    "three": 3,
    "four": 4,
    "five": 5,
    "six": 6,
    "seven": 7,
    "eight": 8,
    "nine": 9,
    "ten": 10,
    "eleven": 11,
    "twelve": 12,
    "thirteen": 13,
    "fourteen": 14,
    "fifteen": 15,
    "sixteen": 16,
    "seventeen": 17,
    "eighteen": 18,
    "nineteen": 19,
    "twenty": 20,
}

NUMBER_TOKEN = r"(?:\d{1,6}|" + "|".join(NUMBER_WORDS) + r")"
ENTITY_TOKEN = (
    r"patients?|subjects?|participants?|individuals?|donors?|volunteers?|"
    r"biological samples?|clinical samples?|samples?|specimens?|biops(?:y|ies)|"
    r"tissue samples?|tissues?|organs?|biological replicates?|technical replicates?|replicates?"
)

COUNT_PATTERNS = [
    re.compile(
        rf"\b(?P<count>{NUMBER_TOKEN})\s+(?P<entity>{ENTITY_TOKEN})\b",
        re.IGNORECASE,
    ),
    re.compile(
        rf"\b(?P<entity>{ENTITY_TOKEN})\s+(?:from|of|in|were|was|included|comprised)\s+"
        rf"(?P<count>{NUMBER_TOKEN})\b",
        re.IGNORECASE,
    ),
    re.compile(
        rf"\b(?:cohort|group|set|panel)\s+of\s+(?P<count>{NUMBER_TOKEN})\s+"
        rf"(?P<entity>{ENTITY_TOKEN})\b",
        re.IGNORECASE,
    ),
    re.compile(rf"\bn\s*[=:]\s*(?P<count>{NUMBER_TOKEN})\b", re.IGNORECASE),
]


class NetworkError(RuntimeError):
    """A retried HTTP request did not complete successfully."""


def request_bytes(
    url: str,
    *,
    params: dict[str, Any] | None = None,
    headers: dict[str, str] | None = None,
    timeout: float = 45.0,
    allow_404: bool = False,
    max_bytes: int | None = None,
) -> bytes | None:
    if params:
        separator = "&" if "?" in url else "?"
        url = f"{url}{separator}{urlencode(params)}"
    request_headers = {"Accept": "application/json", "User-Agent": USER_AGENT}
    request_headers.update(headers or {})
    request = Request(url, headers=request_headers)
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            with urlopen(request, timeout=timeout) as response:
                return response.read(max_bytes) if max_bytes else response.read()
        except HTTPError as error:
            if allow_404 and error.code == 404:
                return None
            last_error = error
            if error.code not in {408, 429, 500, 502, 503, 504}:
                break
            retry_after = error.headers.get("Retry-After") if error.headers else None
            delay = float(retry_after) if retry_after and retry_after.isdigit() else 0.7 * (2**attempt)
        except (URLError, TimeoutError, OSError) as error:
            last_error = error
            delay = 0.7 * (2**attempt)
        if attempt < 3:
            time.sleep(min(delay, 15.0))
    raise NetworkError(f"GET failed for {url}: {last_error}")


def get_json(
    url: str,
    *,
    params: dict[str, Any] | None = None,
    timeout: float = 45.0,
    allow_404: bool = False,
) -> Any:
    content = request_bytes(
        url,
        params=params,
        timeout=timeout,
        allow_404=allow_404,
    )
    if content is None:
        return None
    return json.loads(content.decode("utf-8-sig"))


def get_text(
    url: str,
    *,
    timeout: float = 45.0,
    allow_404: bool = False,
) -> str:
    content = request_bytes(
        url,
        timeout=timeout,
        headers={"Accept": "text/plain,text/tab-separated-values,*/*"},
        allow_404=allow_404,
    )
    if content is None:
        return ""
    return content.decode("utf-8-sig", errors="replace")


def clean_text(value: Any) -> str:
    text = html.unescape(str(value or ""))
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def unique(values: Iterable[Any]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        text = clean_text(value)
        key = text.casefold()
        if text and key not in seen:
            seen.add(key)
            result.append(text)
    return result


def flatten_values(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, (int, float)):
        return [str(value)]
    if isinstance(value, list):
        return unique(item for child in value for item in flatten_values(child))
    if isinstance(value, dict):
        preferred = []
        for key in ("name", "value", "label", "title"):
            if key in value and value[key] not in (None, ""):
                preferred.extend(flatten_values(value[key]))
        if preferred:
            return unique(preferred)
    return []


def values_at_keys(value: Any, wanted_keys: set[str]) -> list[str]:
    found: list[str] = []

    def walk(node: Any) -> None:
        if isinstance(node, dict):
            for key, child in node.items():
                normalized = re.sub(r"[^a-z0-9]", "", str(key).casefold())
                if normalized in wanted_keys:
                    found.extend(flatten_values(child))
                walk(child)
        elif isinstance(node, list):
            for child in node:
                walk(child)

    walk(value)
    return unique(found)


def project_files(payload: Any) -> list[dict[str, Any]]:
    """Accept current PRIDE v3 lists and a few older wrapper shapes."""
    if isinstance(payload, list):
        return [item for item in payload if isinstance(item, dict)]
    if not isinstance(payload, dict):
        return []
    for key in ("files", "list", "content"):
        value = payload.get(key)
        if isinstance(value, list):
            return [item for item in value if isinstance(item, dict)]
    embedded = payload.get("_embedded")
    if isinstance(embedded, dict) and isinstance(embedded.get("files"), list):
        return [item for item in embedded["files"] if isinstance(item, dict)]
    return []


def file_name(record: dict[str, Any]) -> str:
    return clean_text(record.get("fileName") or record.get("name"))


def file_category(record: dict[str, Any]) -> str:
    value = record.get("fileCategory")
    if isinstance(value, dict):
        return clean_text(value.get("value") or value.get("name"))
    return clean_text(value or record.get("category"))


def direct_file_url(accession: str, record: dict[str, Any]) -> str:
    locations = record.get("publicFileLocations") or []
    candidates: list[tuple[int, str]] = []
    if isinstance(locations, list):
        for location in locations:
            if not isinstance(location, dict):
                continue
            label = clean_text(location.get("name")).casefold()
            value = clean_text(location.get("value"))
            if not value:
                continue
            if value.startswith("https://"):
                priority = 0
            elif value.startswith("http://"):
                priority = 1
            elif value.startswith("ftp://"):
                value = "https://" + value[len("ftp://") :]
                priority = 2
            elif "aspera" in label:
                priority = 9
            else:
                priority = 8
            candidates.append((priority, value))
    if candidates:
        candidates.sort(key=lambda item: item[0])
        if candidates[0][0] < 8:
            return candidates[0][1]

    name = (
        clean_text(record.get("_archive_file_name"))
        if record.get("_zip_member")
        else file_name(record)
    )
    return (
        "https://www.ebi.ac.uk/pride/ws/archive-file-downloader/files/s3/"
        f"{quote(accession, safe='')}/{quote(name, safe='')}"
    )


def ranged_bytes(url: str, start: int, end: int, timeout: float) -> tuple[bytes, str]:
    """Read one byte range and return its Content-Range header.

    This deliberately refuses a server that ignores Range, so a very large ZIP
    can never be accidentally downloaded in full.
    """
    request = Request(
        url,
        headers={
            "Accept": "application/zip,application/octet-stream,*/*",
            "Range": f"bytes={start}-{end}",
            "User-Agent": USER_AGENT,
        },
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            content_range = response.headers.get("Content-Range", "")
            if response.status != 206 or not content_range:
                return b"", ""
            return response.read(end - start + 1), content_range
    except (HTTPError, URLError, TimeoutError, OSError):
        return b"", ""


def suffix_bytes(url: str, count: int, timeout: float) -> tuple[bytes, str]:
    """Read the final ``count`` bytes, requiring HTTP range support."""
    request = Request(
        url,
        headers={
            "Accept": "application/zip,application/octet-stream,*/*",
            "Range": f"bytes=-{count}",
            "User-Agent": USER_AGENT,
        },
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            content_range = response.headers.get("Content-Range", "")
            if response.status != 206 or not content_range:
                return b"", ""
            return response.read(count), content_range
    except (HTTPError, URLError, TimeoutError, OSError):
        return b"", ""


def zip_member_records(
    accession: str,
    files: list[dict[str, Any]],
    timeout: float,
    max_central_directory_bytes: int,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Find targets and filenames in public ZIPs using HTTP ranges only."""
    matches: list[dict[str, Any]] = []
    archived_names: list[str] = []
    for archive in files:
        archive_name = unquote(file_name(archive)).replace("\\", "/")
        if not archive_name.casefold().endswith(".zip"):
            continue
        url = direct_file_url(accession, archive)
        tail, content_range = suffix_bytes(url, 65_557, timeout)
        total_match = re.search(r"/(\d+)$", content_range)
        if not total_match:
            continue
        total = int(total_match.group(1))
        eocd_at = tail.rfind(ZIP_EOCD_SIGNATURE)
        if eocd_at < 0 or eocd_at + 22 > len(tail):
            continue
        eocd = tail[eocd_at : eocd_at + 22]
        central_size = struct.unpack_from("<I", eocd, 12)[0]
        central_offset = struct.unpack_from("<I", eocd, 16)[0]
        if central_size == 0xFFFFFFFF or central_offset == 0xFFFFFFFF:
            # ZIP64 locator is immediately before the ordinary EOCD record.
            locator_at = eocd_at - 20
            if locator_at < 0 or tail[locator_at : locator_at + 4] != b"PK\x06\x07":
                continue
            zip64_offset = struct.unpack_from("<Q", tail, locator_at + 8)[0]
            zip64, _ = ranged_bytes(url, zip64_offset, zip64_offset + 55, timeout)
            if zip64[:4] != b"PK\x06\x06" or len(zip64) < 56:
                continue
            central_size = struct.unpack_from("<Q", zip64, 40)[0]
            central_offset = struct.unpack_from("<Q", zip64, 48)[0]
        if central_size < 1 or central_size > max_central_directory_bytes:
            continue
        central, _ = ranged_bytes(url, central_offset, central_offset + central_size - 1, timeout)
        position = 0
        while position + 46 <= len(central) and central[position : position + 4] == ZIP_CENTRAL_SIGNATURE:
            flags = struct.unpack_from("<H", central, position + 8)[0]
            compression = struct.unpack_from("<H", central, position + 10)[0]
            compressed_size = struct.unpack_from("<I", central, position + 20)[0]
            name_length, extra_length, comment_length = struct.unpack_from("<HHH", central, position + 28)
            local_offset = struct.unpack_from("<I", central, position + 42)[0]
            end = position + 46 + name_length + extra_length + comment_length
            if end > len(central):
                break
            raw_name = central[position + 46 : position + 46 + name_length]
            encoding = "utf-8" if flags & 0x800 else "cp437"
            member_name = raw_name.decode(encoding, errors="replace").replace("\\", "/")
            if member_name and not member_name.endswith("/"):
                archived_names.append(f"{archive_name}!/{member_name}")
            if PurePosixPath(member_name).name == TARGET_FILE:
                record = dict(archive)
                record["_archive_file_name"] = archive_name
                record["fileName"] = f"{archive_name}!/{member_name}"
                record["fileCategory"] = "proteinGroups.txt inside ZIP archive"
                record["_zip_member"] = {
                    "name": member_name,
                    "local_offset": local_offset,
                    "compressed_size": compressed_size,
                    "compression": compression,
                    "encrypted": bool(flags & 1),
                }
                matches.append(record)
            position = end
    return matches, unique(archived_names)


def exact_target_files(files: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Require the canonical MaxQuant basename, including its capitalization."""
    matches = []
    for record in files:
        name = unquote(file_name(record)).replace("\\", "/")
        if PurePosixPath(name).name == TARGET_FILE:
            matches.append(record)
    return matches


def mapping_key(value: str) -> str:
    """Normalize an LFQ label or submitted filename for conservative matching."""
    name = unquote(value).replace("\\", "/").rsplit("/", 1)[-1]
    # Strip common mass-spectrometry and table/archive extensions repeatedly.
    extensions = {
        ".raw", ".mzml", ".mzxml", ".wiff", ".wiff2", ".d", ".mgf",
        ".txt", ".tsv", ".csv", ".zip", ".gz", ".tar",
    }
    while True:
        suffix = PurePosixPath(name).suffix.casefold()
        if suffix not in extensions:
            break
        name = name[: -len(suffix)]
    return re.sub(r"[^a-z0-9]+", "", name.casefold())


def map_lfq_names_to_pride_files(
    lfq_names: list[str], files: list[dict[str, Any]], archived_names: list[str]
) -> dict[str, Any]:
    """Map LFQ column labels to PRIDE manifest basenames by normalized equality."""
    if not lfq_names:
        return {
            "status": "unverified",
            "mapped_count": "",
            "unmapped_count": "",
            "unmapped": [],
            "matches": [],
        }
    index: dict[str, list[str]] = {}
    submitted_names = [file_name(record) for record in files] + archived_names
    for submitted_name in submitted_names:
        display_name = unquote(submitted_name).replace("\\", "/")
        key = mapping_key(display_name)
        if key:
            index.setdefault(key, []).append(display_name)

    matches: list[str] = []
    unmapped: list[str] = []
    for label in lfq_names:
        candidates = unique(index.get(mapping_key(label), []))
        if candidates:
            matches.append(f"{label} -> {', '.join(candidates)}")
        else:
            unmapped.append(label)
    mapped_count = len(lfq_names) - len(unmapped)
    if mapped_count == len(lfq_names):
        status = "yes"
    elif mapped_count:
        status = "partial"
    else:
        status = "no"
    return {
        "status": status,
        "mapped_count": mapped_count,
        "unmapped_count": len(unmapped),
        "unmapped": unmapped,
        "matches": matches,
    }


def omicsdi_search(query: str, max_candidates: int, timeout: float) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    start = 0
    while len(results) < max_candidates:
        size = min(100, max_candidates - len(results))
        payload = get_json(
            OMID_SEARCH,
            params={"query": query, "start": start, "size": size},
            timeout=timeout,
        )
        batch = payload.get("datasets") if isinstance(payload, dict) else None
        if not isinstance(batch, list) or not batch:
            break
        results.extend(item for item in batch if isinstance(item, dict))
        start += len(batch)
        total = int(payload.get("count") or 0)
        if len(batch) < size or (total and start >= total):
            break
    return results[:max_candidates]


def discover_candidates(args: argparse.Namespace) -> list[dict[str, Any]]:
    if args.accessions:
        return [
            {"id": accession.upper(), "source": "PRIDE", "_query": "manual accession"}
            for accession in unique(args.accessions)
        ]

    organism_query = ORGANISM_CONFIG[args.organism]["query"]
    queries = [organism_query]
    if args.test:
        queries = [
            f'"{TARGET_FILE}" AND {organism_query}',
            f"MaxQuant AND {organism_query}",
            organism_query,
        ]

    candidate_cap = args.max_candidates
    if args.test:
        # Keep test mode genuinely small while allowing enough candidates to
        # find exact-file matches after PRIDE verification.
        candidate_cap = min(candidate_cap, max(50, args.test_limit * 25))

    merged: dict[str, dict[str, Any]] = {}
    remaining = candidate_cap
    for query in queries:
        if remaining <= 0:
            break
        print(f"OmicsDI query: {query}", file=sys.stderr)
        for item in omicsdi_search(query, remaining, args.timeout):
            accession = clean_text(item.get("id") or item.get("accession")).upper()
            source = clean_text(item.get("source") or item.get("repository"))
            if not re.fullmatch(r"PXD\d+", accession):
                continue
            if source and source.casefold() != "pride":
                continue
            if accession not in merged:
                merged[accession] = dict(item)
                merged[accession]["_query"] = query
        remaining = candidate_cap - len(merged)
        if not args.test:
            break

    return list(merged.values())[:candidate_cap]


def omicsdi_detail(accession: str, timeout: float) -> dict[str, Any]:
    payload = get_json(
        f"{OMID_BASE}/dataset/pride/{quote(accession, safe='')}",
        timeout=timeout,
        allow_404=True,
    )
    return payload if isinstance(payload, dict) else {}


def pride_project(accession: str, timeout: float) -> dict[str, Any]:
    payload = get_json(
        f"{PRIDE_BASE}/projects/{quote(accession, safe='')}",
        timeout=timeout,
        allow_404=True,
    )
    return payload if isinstance(payload, dict) else {}


def pride_files(accession: str, timeout: float) -> list[dict[str, Any]]:
    base = f"{PRIDE_BASE}/projects/{quote(accession, safe='')}/files"
    try:
        payload = get_json(
            f"{base}/all",
            timeout=max(timeout, 120.0),
            allow_404=True,
        )
        files = project_files(payload)
        if files or payload == []:
            return files
    except NetworkError:
        pass

    # Compatibility fallback for deployments without /files/all.
    files: list[dict[str, Any]] = []
    page_size = 100
    for page in range(10_000):
        payload = get_json(
            base,
            params={"pageSize": page_size, "page": page},
            timeout=max(timeout, 90.0),
            allow_404=True,
        )
        batch = project_files(payload)
        if not batch:
            break
        files.extend(batch)
        if len(batch) < page_size:
            break
    return files


def additional_values(detail: dict[str, Any], key: str) -> list[str]:
    additional = detail.get("additional")
    if not isinstance(additional, dict):
        return []
    return unique(flatten_values(additional.get(key)))


def structured_organisms(
    search_item: dict[str, Any], detail: dict[str, Any], project: dict[str, Any]
) -> list[str]:
    values: list[str] = []
    values.extend(flatten_values(search_item.get("organisms")))
    values.extend(additional_values(detail, "species"))
    values.extend(additional_values(detail, "organism"))
    values.extend(values_at_keys(project, {"organisms", "species"}))
    return unique(values)


def is_target_organism(
    target: str,
    organisms: Iterable[str],
    detail: dict[str, Any],
    project: dict[str, Any],
) -> bool:
    text = " ".join(organisms)
    text += " " + json.dumps(detail, ensure_ascii=False)
    text += " " + json.dumps(project.get("organisms") or [], ensure_ascii=False)
    if target == "human":
        pattern = (
            r"\bHomo\s+sapiens\b|\bhuman\b|\bNCBITaxon\s*:\s*9606\b|"
            r"(?:NCBI\s*)?Taxon(?:omy)?\s*[:=]?\s*9606\b"
        )
    elif target == "zebrafish":
        pattern = (
            r"\bDanio\s+rerio\b|\bzebra[\s-]?fish\b|\bNCBITaxon\s*:\s*7955\b|"
            r"(?:NCBI\s*)?Taxon(?:omy)?\s*[:=]?\s*7955\b"
        )
    else:  # Protected by argparse choices; retained for direct function use.
        raise ValueError(f"Unsupported organism: {target}")
    return bool(re.search(pattern, text, re.I))


def structured_tissues_from_pride(project: dict[str, Any]) -> list[str]:
    return values_at_keys(
        project,
        {
            "organismpart",
            "organismparts",
            "tissue",
            "tissues",
            "bodypart",
            "bodyparts",
            "bodysite",
        },
    )


def structured_tissues_from_omicsdi(detail: dict[str, Any]) -> list[str]:
    return unique(
        additional_values(detail, "tissue")
        + additional_values(detail, "organism_part")
        + additional_values(detail, "body_site")
    )


def fetch_sdrf(accession: str, files: list[dict[str, Any]], timeout: float) -> str:
    try:
        text = get_text(
            f"{PRIDE_BASE}/files/sdrf/{quote(accession, safe='')}",
            timeout=timeout,
            allow_404=True,
        )
        if "\t" in text and "\n" in text:
            return text
    except NetworkError:
        pass

    for record in files:
        name = file_name(record).casefold()
        if name.endswith("sdrf.tsv") or name.endswith("sdrf.txt"):
            try:
                text = get_text(direct_file_url(accession, record), timeout=timeout)
                if "\t" in text and "\n" in text:
                    return text
            except NetworkError:
                continue
    return ""


def parse_sdrf(text: str) -> dict[str, Any]:
    result = {
        "row_count": 0,
        "unique_sources": [],
        "unique_samples": [],
        "tissues": [],
        "organisms": [],
    }
    if not text:
        return result

    reader = csv.DictReader(io.StringIO(text.lstrip("\ufeff")), delimiter="\t")
    if not reader.fieldnames:
        return result

    source_columns: list[str] = []
    sample_columns: list[str] = []
    tissue_columns: list[str] = []
    organism_columns: list[str] = []
    for column in reader.fieldnames:
        normalized = re.sub(r"\s+", " ", column.casefold()).strip()
        if normalized == "source name":
            source_columns.append(column)
        if normalized == "sample name":
            sample_columns.append(column)
        if any(term in normalized for term in ("organism part", "tissue", "body site")):
            tissue_columns.append(column)
        if "organism" in normalized and "part" not in normalized:
            organism_columns.append(column)

    sources: list[str] = []
    samples: list[str] = []
    tissues: list[str] = []
    organisms: list[str] = []
    missing = {"", "not available", "not applicable", "na", "n/a", "unknown"}
    for row in reader:
        result["row_count"] += 1
        for column, target in (
            *((column, sources) for column in source_columns),
            *((column, samples) for column in sample_columns),
            *((column, tissues) for column in tissue_columns),
            *((column, organisms) for column in organism_columns),
        ):
            value = clean_text(row.get(column))
            if value.casefold() not in missing:
                target.append(value)

    result["unique_sources"] = unique(sources)
    result["unique_samples"] = unique(samples)
    result["tissues"] = unique(tissues)
    result["organisms"] = unique(organisms)
    return result


def text_fields(project: dict[str, Any], detail: dict[str, Any]) -> dict[str, str]:
    pride_values: list[str] = []
    for key in (
        "title",
        "projectTitle",
        "projectDescription",
        "description",
        "sampleProcessingProtocol",
        "dataProcessingProtocol",
        "keywords",
    ):
        pride_values.extend(flatten_values(project.get(key)))

    omicsdi_values: list[str] = []
    for key in ("name", "title", "description", "keywords"):
        omicsdi_values.extend(flatten_values(detail.get(key)))
    return {
        "PRIDE metadata": clean_text(" ".join(unique(pride_values))),
        "OmicsDI metadata": clean_text(" ".join(unique(omicsdi_values))),
    }


def extract_publication_ids(project: dict[str, Any], detail: dict[str, Any]) -> tuple[str, str]:
    references = project.get("references") or []
    pmids = values_at_keys(references, {"pmid", "pubmedid", "pubmed"})
    dois = values_at_keys(references, {"doi"})
    if not pmids:
        pmids = values_at_keys(detail, {"pmid", "pubmedid", "pubmed"})
    if not dois:
        dois = values_at_keys(detail, {"doi"})
    if not dois:
        dois = flatten_values(project.get("doi"))

    # Some PRIDE records store identifiers only inside a formatted reference line.
    reference_text = json.dumps(references, ensure_ascii=False)
    if not pmids:
        pmids = re.findall(r"(?:PMID|PubMed)\D{0,12}(\d{5,10})", reference_text, flags=re.I)
    if not dois:
        dois = re.findall(r"10\.\d{4,9}/[^\s\"<>]+", reference_text, flags=re.I)

    def clean_pmid(value: str) -> str:
        match = re.search(r"\b\d{5,10}\b", value)
        return match.group(0) if match else ""

    def clean_doi(value: str) -> str:
        value = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", value, flags=re.I)
        value = re.sub(r"^doi:\s*", "", value, flags=re.I).strip()
        match = re.search(r"10\.\d{4,9}/\S+", value, flags=re.I)
        return match.group(0).rstrip(".,;)") if match else ""

    pmid = next((clean_pmid(value) for value in pmids if clean_pmid(value)), "")
    doi = next((clean_doi(value) for value in dois if clean_doi(value)), "")
    return pmid, doi


def fetch_publication(pmid: str, doi: str, timeout: float) -> dict[str, Any]:
    queries = []
    if pmid:
        queries.append(f"EXT_ID:{pmid} AND SRC:MED")
    if doi:
        escaped = doi.replace('"', "")
        queries.append(f'DOI:"{escaped}"')

    for query in queries:
        try:
            payload = get_json(
                EPMC_SEARCH,
                params={
                    "query": query,
                    "format": "json",
                    "resultType": "core",
                    "pageSize": 1,
                },
                timeout=timeout,
            )
        except NetworkError:
            continue
        results = (
            payload.get("resultList", {}).get("result", [])
            if isinstance(payload, dict)
            else []
        )
        if isinstance(results, list) and results and isinstance(results[0], dict):
            return results[0]
    return {}


def publication_link(publication: dict[str, Any], pmid: str, doi: str) -> str:
    pmid = clean_text(publication.get("pmid") or pmid)
    doi = clean_text(publication.get("doi") or doi)
    if pmid:
        return f"https://europepmc.org/article/MED/{quote(pmid, safe='')}"
    if doi:
        return f"https://doi.org/{quote(doi, safe='/()') }"
    return ""


def short_snippet(text: str, start: int, end: int, max_words: int = 18) -> str:
    left = max(text.rfind(".", 0, start), text.rfind(";", 0, start), text.rfind("\n", 0, start))
    right_candidates = [position for marker in (".", ";", "\n") if (position := text.find(marker, end)) >= 0]
    right = min(right_candidates) if right_candidates else min(len(text), end + 180)
    snippet = clean_text(text[left + 1 : right + 1])
    words = snippet.split()
    if len(words) > max_words:
        snippet = " ".join(words[:max_words]) + " …"
    return snippet


def tissue_mentions(text: str, source: str) -> list[dict[str, str]]:
    found: list[dict[str, str]] = []
    if not text:
        return found
    for canonical, aliases in TISSUE_TERMS.items():
        best: re.Match[str] | None = None
        for alias in sorted(aliases, key=len, reverse=True):
            pattern = re.compile(rf"(?<![A-Za-z0-9]){re.escape(alias)}(?![A-Za-z0-9])", re.I)
            match = pattern.search(text)
            if match and (best is None or len(match.group(0)) > len(best.group(0))):
                best = match
        if best:
            found.append(
                {
                    "tissue": canonical,
                    "source": source,
                    "snippet": short_snippet(text, best.start(), best.end()),
                }
            )

    names = {item["tissue"] for item in found}
    suppress: set[str] = set()
    for specific, broad_names in SUPPRESS_BROAD_TISSUES.items():
        if specific in names:
            suppress.update(broad_names)
    return [item for item in found if item["tissue"] not in suppress]


def number_value(token: str) -> int:
    token = token.casefold()
    return int(token) if token.isdigit() else NUMBER_WORDS[token]


def normalized_entity(value: str) -> str:
    text = value.casefold()
    if "replicate" in text:
        return "replicates"
    if text.startswith("patient"):
        return "patients"
    if text.startswith("subject"):
        return "subjects"
    if text.startswith("participant"):
        return "participants"
    if text.startswith("individual"):
        return "individuals"
    if text.startswith("donor"):
        return "donors"
    if text.startswith("volunteer"):
        return "volunteers"
    if "biops" in text:
        return "biopsies"
    if text.startswith("specimen"):
        return "specimens"
    if "tissue" in text:
        return "tissues"
    if text.startswith("organ"):
        return "organs"
    return "samples"


def count_mentions(text: str, source: str) -> list[dict[str, Any]]:
    found: list[dict[str, Any]] = []
    occupied: set[tuple[int, int]] = set()
    for pattern in COUNT_PATTERNS:
        for match in pattern.finditer(text or ""):
            span = match.span()
            if any(max(span[0], other[0]) < min(span[1], other[1]) for other in occupied):
                continue
            occupied.add(span)
            entity = normalized_entity(match.groupdict().get("entity") or "n")
            if match.groupdict().get("entity") is None:
                entity = "unspecified n"
            found.append(
                {
                    "count": number_value(match.group("count")),
                    "entity": entity,
                    "source": source,
                    "snippet": short_snippet(text, match.start(), match.end()),
                }
            )
    deduplicated: dict[tuple[int, str, str], dict[str, Any]] = {}
    for item in found:
        key = (item["count"], item["entity"], item["source"])
        deduplicated.setdefault(key, item)
    return list(deduplicated.values())


def canonicalize_structured_tissues(values: Iterable[str]) -> list[str]:
    result: list[str] = []
    for value in values:
        mentions = tissue_mentions(value, "structured metadata")
        if mentions:
            result.extend(item["tissue"] for item in mentions)
        else:
            result.append(value)
    return unique(result)


def parse_protein_groups_header(content: bytes) -> dict[str, Any]:
    result = {"status": "unverified", "lfq_columns": [], "ordinary_intensity_count": ""}
    if b"\n" not in content:
        return result
    first_line = content.splitlines()[0].decode("utf-8-sig", errors="replace")
    columns = [column.strip() for column in first_line.split("\t")]
    lfq_columns = [
        column for column in columns
        if re.fullmatch(r"LFQ\s+intensity(?:\s+.+)?", column, flags=re.I)
    ]
    ordinary = [
        column for column in columns
        if re.fullmatch(r"Intensity\s+.+", column, flags=re.I)
        and not re.fullmatch(r"LFQ\s+intensity\s+.+", column, flags=re.I)
    ]
    result["status"] = "yes" if lfq_columns else "no"
    result["lfq_columns"] = lfq_columns
    result["ordinary_intensity_count"] = len(ordinary)
    return result


def read_zip_member_header(url: str, member: dict[str, Any], timeout: float) -> dict[str, Any]:
    """Read a ZIP member's first table line without downloading its archive."""
    result = {"status": "unverified", "lfq_columns": [], "ordinary_intensity_count": ""}
    local_offset = member.get("local_offset")
    compressed_size = member.get("compressed_size")
    compression = member.get("compression")
    if member.get("encrypted"):
        return result
    if not isinstance(local_offset, int) or not isinstance(compressed_size, int):
        return result
    local, _ = ranged_bytes(url, local_offset, local_offset + 29, timeout)
    if len(local) < 30 or local[:4] != ZIP_LOCAL_SIGNATURE:
        return result
    name_length, extra_length = struct.unpack_from("<HH", local, 26)
    data_start = local_offset + 30 + name_length + extra_length
    # A 1 MiB compressed prefix is ample for a tabular header but keeps network
    # use bounded. Deflate can be decoded incrementally even when incomplete.
    chunk_size = min(compressed_size, 1_048_576)
    payload, _ = ranged_bytes(url, data_start, data_start + chunk_size - 1, timeout)
    try:
        if compression == 0:
            decoded = payload
        elif compression == 8:
            decoded = zlib.decompressobj(-15).decompress(payload, 2_000_000)
        else:
            return result
    except zlib.error:
        return result
    return parse_protein_groups_header(decoded)


def read_table_header(url: str, timeout: float, max_bytes: int = 512_000) -> dict[str, Any]:
    result = {
        "status": "unverified",
        "lfq_columns": [],
        "ordinary_intensity_count": "",
    }
    if not url.startswith(("http://", "https://")):
        return result
    try:
        content = request_bytes(
            url,
            headers={"Accept": "text/plain,*/*", "Range": f"bytes=0-{max_bytes - 1}"},
            timeout=timeout,
            max_bytes=max_bytes,
        )
        return parse_protein_groups_header(content or b"")
    except NetworkError:
        return result


def count_hint_text(items: Iterable[dict[str, Any]]) -> str:
    values = []
    for item in items:
        values.append(
            f'{item["source"]}: {item["count"]} {item["entity"]} [{item["snippet"]}]'
        )
    return " | ".join(unique(values))


def choose_sample_estimate(
    sdrf: dict[str, Any], count_evidence: list[dict[str, Any]]
) -> tuple[str, str, str]:
    samples = sdrf.get("unique_samples") or []
    sources = sdrf.get("unique_sources") or []
    if samples:
        return str(len(samples)), "SDRF unique 'sample name' values", "high"
    if sources:
        return str(len(sources)), "SDRF unique 'source name' values", "high"

    source_rank = {"PRIDE metadata": 0, "OmicsDI metadata": 1, "publication abstract": 2}
    entity_rank = {
        "samples": 0,
        "specimens": 1,
        "biopsies": 1,
        "donors": 2,
        "patients": 3,
        "subjects": 3,
        "participants": 3,
        "individuals": 3,
        "volunteers": 3,
        "tissues": 4,
        "organs": 4,
        "unspecified n": 9,
    }
    candidates = [item for item in count_evidence if item["entity"] != "replicates"]
    if not candidates:
        return "", "", ""
    candidates.sort(
        key=lambda item: (
            source_rank.get(item["source"], 9),
            entity_rank.get(item["entity"], 8),
        )
    )
    best = candidates[0]
    confidence = "medium" if best["source"] == "PRIDE metadata" and best["entity"] != "unspecified n" else "low"
    return (
        str(best["count"]),
        f'{best["source"]}: explicit {best["entity"]} statement',
        confidence,
    )


def process_candidate(item: dict[str, Any], args: argparse.Namespace) -> dict[str, Any] | None:
    accession = clean_text(item.get("id") or item.get("accession")).upper()
    try:
        files = pride_files(accession, args.timeout)
        targets = exact_target_files(files)
        archived_names: list[str] = []
        if args.inspect_zip_archives:
            zipped_targets, archived_names = zip_member_records(
                accession,
                files,
                args.timeout,
                int(args.zip_central_directory_max_mb * 1_000_000),
            )
            targets.extend(zipped_targets)
        if not targets:
            return None

        project = pride_project(accession, args.timeout)
        detail = omicsdi_detail(accession, args.timeout)
        organisms = structured_organisms(item, detail, project)

        sdrf_text = fetch_sdrf(accession, files, args.timeout)
        sdrf = parse_sdrf(sdrf_text)
        organisms = unique(organisms + list(sdrf.get("organisms") or []))
        organism_verified = is_target_organism(args.organism, organisms, detail, project)
        if not organism_verified and not args.allow_unconfirmed_organism:
            return None

        target_urls = [direct_file_url(accession, record) for record in targets]
        if args.skip_header_check:
            header = {"status": "unverified", "lfq_columns": [], "ordinary_intensity_count": ""}
        elif targets[0].get("_zip_member"):
            header = read_zip_member_header(target_urls[0], targets[0]["_zip_member"], args.timeout)
        else:
            header = read_table_header(target_urls[0], args.timeout)
        if args.require_lfq_columns and header["status"] != "yes":
            return None

        pmid, doi = extract_publication_ids(project, detail)
        publication = {} if args.skip_publications else fetch_publication(pmid, doi, args.timeout)
        pmid = clean_text(publication.get("pmid") or pmid)
        doi = clean_text(publication.get("doi") or doi)

        sources = text_fields(project, detail)
        abstract = clean_text(publication.get("abstractText"))
        if abstract:
            sources["publication abstract"] = abstract

        tissue_evidence: list[dict[str, str]] = []
        count_evidence: list[dict[str, Any]] = []
        for source, text in sources.items():
            tissue_evidence.extend(tissue_mentions(text, source))
            count_evidence.extend(count_mentions(text, source))

        pride_tissues = structured_tissues_from_pride(project)
        omicsdi_tissues = structured_tissues_from_omicsdi(detail)
        sdrf_tissues = list(sdrf.get("tissues") or [])
        mined_tissues = unique(item["tissue"] for item in tissue_evidence)
        combined_tissues = unique(
            canonicalize_structured_tissues(pride_tissues + omicsdi_tissues + sdrf_tissues)
            + mined_tissues
        )

        sample_counts = [item for item in count_evidence if item["entity"] != "replicates"]
        replicate_counts = [item for item in count_evidence if item["entity"] == "replicates"]
        estimate, basis, confidence = choose_sample_estimate(sdrf, count_evidence)

        lfq_columns: list[str] = list(header["lfq_columns"])
        lfq_names = [
            re.sub(r"^LFQ\s+intensity\s*", "", column, flags=re.I).strip()
            for column in lfq_columns
        ]
        lfq_names = [name for name in lfq_names if name]
        lfq_mapping = map_lfq_names_to_pride_files(lfq_names, files, archived_names)

        categories = unique(file_category(record) for record in targets)
        sizes = []
        for record in targets:
            try:
                size_bytes = (
                    record.get("_zip_member", {}).get("compressed_size")
                    if record.get("_zip_member")
                    else record.get("fileSizeBytes")
                )
                sizes.append(f'{float(size_bytes or 0) / 1_000_000:.3f}')
            except (TypeError, ValueError):
                sizes.append("")

        locations = ["inside ZIP" if record.get("_zip_member") else "direct PRIDE file" for record in targets]
        archive_names = [
            file_name(record).split("!/", 1)[0]
            for record in targets
            if record.get("_zip_member")
        ]

        publication_date = clean_text(
            project.get("publicationDate")
            or item.get("publicationDate")
            or detail.get("publicationDate")
        )
        title = clean_text(
            project.get("title")
            or project.get("projectTitle")
            or item.get("title")
            or detail.get("name")
        )

        pride_sample_evidence = [
            evidence["snippet"]
            for evidence in count_evidence
            if evidence["source"] == "PRIDE metadata"
        ]
        abstract_sample_evidence = [
            evidence["snippet"]
            for evidence in count_evidence
            if evidence["source"] == "publication abstract"
        ]

        return {
            "accession": accession,
            "title": title,
            "publication_date": publication_date,
            "target_organism": ORGANISM_CONFIG[args.organism]["label"],
            "organism_verified": "yes" if organism_verified else "OmicsDI query only",
            "human_verified": (
                ("yes" if organism_verified else "OmicsDI query only")
                if args.organism == "human"
                else "not applicable"
            ),
            "organisms": "; ".join(organisms),
            "protein_groups_verified": "yes — exact canonical basename",
            "protein_groups_location": "; ".join(unique(locations)),
            "protein_groups_archive_name": "; ".join(unique(archive_names)),
            "protein_groups_file_name": "; ".join(file_name(record) for record in targets),
            "protein_groups_file_category": "; ".join(categories),
            "protein_groups_file_size_mb": "; ".join(sizes),
            "protein_groups_url": "; ".join(target_urls),
            "lfq_columns_detected": header["status"],
            "lfq_column_count": len(lfq_columns) if header["status"] != "unverified" else "",
            "lfq_sample_names": "; ".join(lfq_names),
            "lfq_names_mappable_to_pride_files": lfq_mapping["status"],
            "lfq_names_mapped_count": lfq_mapping["mapped_count"],
            "lfq_names_unmapped_count": lfq_mapping["unmapped_count"],
            "lfq_names_unmapped": "; ".join(lfq_mapping["unmapped"]),
            "lfq_name_file_matches": " | ".join(lfq_mapping["matches"]),
            "ordinary_intensity_column_count": header["ordinary_intensity_count"],
            "pride_structured_tissues": "; ".join(pride_tissues),
            "omicsdi_structured_tissues": "; ".join(omicsdi_tissues),
            "sdrf_tissues": "; ".join(sdrf_tissues),
            "text_mined_tissues": "; ".join(mined_tissues),
            "combined_tissues": "; ".join(combined_tissues),
            "tissue_count": len(combined_tissues),
            "sdrf_row_count": sdrf["row_count"] or "",
            "sdrf_unique_source_count": len(sdrf["unique_sources"]) or "",
            "sdrf_unique_sample_count": len(sdrf["unique_samples"]) or "",
            "sample_count_estimate": estimate,
            "sample_count_basis": basis,
            "sample_count_confidence": confidence,
            "sample_count_hints": count_hint_text(sample_counts),
            "replicate_count_hints": count_hint_text(replicate_counts),
            "sample_evidence_pride": " | ".join(unique(pride_sample_evidence)),
            "sample_evidence_abstract": " | ".join(unique(abstract_sample_evidence)),
            "tissue_evidence_sources": "; ".join(
                unique(f'{evidence["tissue"]} [{evidence["source"]}]' for evidence in tissue_evidence)
            ),
            "pmid": pmid,
            "doi": doi,
            "publication_title": clean_text(publication.get("title")),
            "publication_url": publication_link(publication, pmid, doi),
            "pride_dataset_url": f"https://www.ebi.ac.uk/pride/archive/projects/{accession}",
            "omicsdi_dataset_url": f"https://www.omicsdi.org/dataset/pride/{accession}",
        }
    except (NetworkError, ValueError, json.JSONDecodeError) as error:
        print(f"Warning: {accession} could not be completed: {error}", file=sys.stderr)
        return None


def process_all(
    candidates: list[dict[str, Any]], args: argparse.Namespace
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    inspected = 0
    batch_size = max(args.workers * 2, 1)
    for start in range(0, len(candidates), batch_size):
        batch = candidates[start : start + batch_size]
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(process_candidate, item, args): item for item in batch}
            for future in as_completed(futures):
                inspected += 1
                accession = clean_text(
                    futures[future].get("id") or futures[future].get("accession")
                )
                try:
                    row = future.result()
                except Exception as error:  # Keep a long search running after one unexpected record.
                    print(f"Warning: {accession} failed: {error}", file=sys.stderr)
                    row = None
                if row:
                    rows.append(row)
                    print(
                        f"Verified {row['accession']} ({len(rows)} retained; {inspected} inspected)",
                        file=sys.stderr,
                    )

        if args.test and len(rows) >= args.test_limit:
            rows = rows[: args.test_limit]
            break
        if inspected % 50 == 0 or inspected == len(candidates):
            print(
                f"Progress: {inspected}/{len(candidates)} candidates; {len(rows)} retained",
                file=sys.stderr,
            )
        if args.request_delay:
            time.sleep(args.request_delay)

    return sorted(rows, key=lambda row: row["accession"])


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Discover human or zebrafish PRIDE proteomics datasets through OmicsDI, require an exact "
            "proteinGroups.txt file (directly or inside ZIP), and extract tissue/sample evidence."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--organism",
        choices=sorted(ORGANISM_CONFIG),
        default="human",
        help="Target organism. Human remains the default for backward compatibility.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("human_proteingroups_report.csv"),
        help="Final CSV report.",
    )
    parser.add_argument(
        "--max-candidates",
        type=int,
        default=5000,
        help="Maximum unique OmicsDI candidates whose PRIDE manifests may be checked.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help=(
            "Limited test mode: prioritize proteinGroups/MaxQuant searches and stop after "
            "--test-limit verified datasets."
        ),
    )
    parser.add_argument(
        "--test-limit",
        type=int,
        default=10,
        help="Number of verified datasets retained in --test mode.",
    )
    parser.add_argument(
        "--accessions",
        nargs="+",
        help="Optional PXD accessions for a focused test; bypasses OmicsDI discovery only.",
    )
    parser.add_argument(
        "--require-lfq-columns",
        action="store_true",
        help="Keep only tables whose header contains one or more 'LFQ intensity ...' columns.",
    )
    parser.add_argument(
        "--skip-header-check",
        action="store_true",
        help="Do not read the beginning of proteinGroups.txt to inspect quantitative columns.",
    )
    parser.add_argument(
        "--inspect-zip-archives",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Inspect public ZIP central directories for proteinGroups.txt. Use "
            "--no-inspect-zip-archives to restore direct-file-only behavior."
        ),
    )
    parser.add_argument(
        "--zip-central-directory-max-mb",
        type=float,
        default=25.0,
        help="Maximum ZIP central-directory size read per archive; full ZIPs are never downloaded.",
    )
    parser.add_argument(
        "--skip-publications",
        action="store_true",
        help="Skip Europe PMC publication lookup and abstract text mining.",
    )
    parser.add_argument(
        "--allow-unconfirmed-organism",
        "--allow-unconfirmed-human",
        dest="allow_unconfirmed_organism",
        action="store_true",
        help=(
            "Keep OmicsDI query hits even when structured PRIDE/OmicsDI/SDRF metadata "
            "cannot reconfirm the selected organism. The old --allow-unconfirmed-human "
            "name remains accepted as an alias."
        ),
    )
    parser.add_argument("--workers", type=int, default=4, help="Concurrent dataset workers (1-8).")
    parser.add_argument("--timeout", type=float, default=45.0, help="HTTP timeout in seconds.")
    parser.add_argument(
        "--request-delay",
        type=float,
        default=0.15,
        help="Polite delay in seconds between worker batches.",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.max_candidates < 1:
        raise ValueError("--max-candidates must be positive")
    if args.test_limit < 1:
        raise ValueError("--test-limit must be positive")
    if not 1 <= args.workers <= 8:
        raise ValueError("--workers must be between 1 and 8")
    if args.timeout <= 0:
        raise ValueError("--timeout must be positive")
    if args.zip_central_directory_max_mb <= 0:
        raise ValueError("--zip-central-directory-max-mb must be positive")
    if args.require_lfq_columns and args.skip_header_check:
        raise ValueError("--require-lfq-columns cannot be combined with --skip-header-check")
    if args.accessions:
        invalid = [value for value in args.accessions if not re.fullmatch(r"PXD\d+", value, re.I)]
        if invalid:
            raise ValueError(f"Invalid PRIDE accession(s): {', '.join(invalid)}")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2

    try:
        candidates = discover_candidates(args)
    except (NetworkError, ValueError, json.JSONDecodeError) as error:
        print(f"Discovery failed: {error}", file=sys.stderr)
        return 1

    print(f"Discovered {len(candidates)} unique PRIDE candidates", file=sys.stderr)
    rows = process_all(candidates, args)
    write_csv(args.output, rows)
    print(f"Saved {len(rows)} verified datasets to {args.output.resolve()}")
    if not rows:
        print(
            "No verified datasets were retained. In a small test this can happen if none of "
            "the inspected manifests or range-readable ZIP archives exposes the exact "
            "canonical proteinGroups.txt basename."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
