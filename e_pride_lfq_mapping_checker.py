#!/usr/bin/env python3
"""Minimal PRIDE LFQ mappability checker (standard library only)."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import re
import tempfile
import time
import urllib.parse
import urllib.request
import zipfile
from collections import OrderedDict, defaultdict
from contextlib import contextmanager
from pathlib import Path


PRIDE_ACCESSIONS = [
    "PXD002854",
]
OUTPUT_FILE = Path("lfq_mapping_report.csv")
MAX_FULL_ZIP_GB = 2.0
API = "https://www.ebi.ac.uk/pride/ws/archive/v2"
UA = "PRIDE-LFQ-minimal/1.0"
LFQ = re.compile(r"^\s*LFQ(?:\s+|_)+intensity(?:\s+|_)*(.*?)\s*$", re.I)
MARKED_REP = re.compile(
    r"(?i)^(.+?)(?:[\s_.-]*(?:bio(?:logical)?rep|tech(?:nical)?rep|"
    r"rep(?:licate)?|repl|experiment|exp|ex|run|r))[\s_.-]*0*(\d+)$"
)


def request(url, headers=None, retries=4):
    headers = {"User-Agent": UA, "Accept-Encoding": "identity", **(headers or {})}
    error = None
    for attempt in range(retries):
        try:
            return urllib.request.urlopen(
                urllib.request.Request(url, headers=headers), timeout=90
            )
        except Exception as exc:
            error = exc
            if attempt + 1 < retries:
                time.sleep(2**attempt)
    raise RuntimeError(f"Request failed: {url}: {error}")


def get_json(url):
    with request(url, {"Accept": "application/json"}) as response:
        return json.loads(response.read().decode("utf-8"), strict=False)


def file_url(item):
    locations = item.get("publicFileLocations") or []
    locations = sorted(
        locations,
        key=lambda x: "FTP" in str(x.get("name", "")),
        reverse=True,
    )
    if not locations:
        return ""
    url = str(locations[0].get("value", ""))
    url = url.replace("ftp://ftp.pride.ebi.ac.uk/", "https://ftp.pride.ebi.ac.uk/")
    url = url.replace("http://", "https://", 1)
    return urllib.parse.quote(url, safe=":/@?&=%+;,$-_.!~*'()")


def project_files(accession):
    output, seen = [], set()
    for page in range(10000):
        url = f"{API}/projects/{accession}/files?page={page}&pageSize=100"
        payload = get_json(url)
        items = payload if isinstance(payload, list) else payload.get("files", [])
        if not items:
            break
        for item in items:
            url, name = file_url(item), str(item.get("fileName", ""))
            key = (name, url)
            if name and url and key not in seen:
                seen.add(key)
                try:
                    size = int(item.get("fileSizeBytes"))
                except (TypeError, ValueError):
                    size = None
                output.append({"name": name, "url": url, "size": size})
        if len(items) < 100:
            break
    return output


def prefix(url, limit):
    with request(url, {"Range": f"bytes=0-{limit - 1}"}) as response:
        return response.read(limit)


class RangeReader(io.RawIOBase):
    def __init__(self, url, size, chunk=1024**2):
        self.url, self.size, self.chunk = url, size, chunk
        self.pos, self.cache = 0, OrderedDict()

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=os.SEEK_SET):
        self.pos = {
            os.SEEK_SET: offset,
            os.SEEK_CUR: self.pos + offset,
            os.SEEK_END: self.size + offset,
        }[whence]
        return self.pos

    def _get(self, index):
        if index in self.cache:
            self.cache.move_to_end(index)
            return self.cache[index]
        start, end = index * self.chunk, min((index + 1) * self.chunk, self.size) - 1
        with request(self.url, {"Range": f"bytes={start}-{end}"}) as response:
            if getattr(response, "status", response.getcode()) != 206:
                raise OSError("Server does not support ZIP byte ranges")
            data = response.read(end - start + 1)
        self.cache[index] = data
        while len(self.cache) > 48:
            self.cache.popitem(last=False)
        return data

    def read(self, size=-1):
        size = self.size - self.pos if size is None or size < 0 else size
        size = min(size, self.size - self.pos)
        out = bytearray()
        while size > 0:
            block = self._get(self.pos // self.chunk)
            offset = self.pos % self.chunk
            take = min(size, len(block) - offset)
            if take <= 0:
                break
            out += block[offset : offset + take]
            self.pos, size = self.pos + take, size - take
        return bytes(out)

    def readinto(self, buffer):
        data = self.read(len(buffer))
        buffer[: len(data)] = data
        return len(data)


@contextmanager
def remote_zip(item, full_limit):
    size = item["size"]
    if size is None:
        with request(item["url"], {"Range": "bytes=0-0"}) as response:
            size = int(response.headers["Content-Range"].rsplit("/", 1)[1])
    try:
        archive = zipfile.ZipFile(RangeReader(item["url"], size))
    except Exception:
        if size > full_limit:
            raise
        with tempfile.TemporaryFile() as handle:
            with request(item["url"]) as response:
                while block := response.read(1024**2):
                    handle.write(block)
            handle.seek(0)
            with zipfile.ZipFile(handle) as archive:
                yield archive
        return
    try:
        yield archive
    finally:
        archive.close()


def base(path):
    return path.replace("\\", "/").rsplit("/", 1)[-1].lower()


def is_pg(path):
    name = base(path)
    return "proteingroups" in name and name.endswith((".txt", ".tsv"))


def metadata_score(path):
    name = base(path)
    if is_pg(name) or not name.endswith((".txt", ".tsv", ".csv", ".xml", ".json")):
        return 0
    words = ["sdrf", "experimentaldesign", "metadata", "manifest", "annotation", "design", "readme", "mqpar", "sample"]
    return len(words) - next((i for i, word in enumerate(words) if word in name), len(words))


def text(data):
    for encoding in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            pass
    return data.decode(errors="replace")


def header(data):
    line = next((x for x in text(data).splitlines() if x.strip()), "")
    delimiter = "\t" if line.count("\t") >= line.count(",") else ","
    return next(csv.reader([line], delimiter=delimiter), [])


def norm(value):
    return "".join(re.findall(r"[a-z0-9]+", value.lower()))


def lfq_columns(columns):
    return [(column, match.group(1).strip()) for column in columns if (match := LFQ.match(column))]


def useful_group(value):
    key = norm(value)
    return len(key) >= 2 and not re.fullmatch(r"[a-z]\d*", key)


def remove_run_prefix(labels):
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


def split_rule(value, rule, width=0):
    value = value.strip(" _-")
    if rule == "marked":
        match = MARKED_REP.fullmatch(value)
        return (match.group(1).strip(" _-"), str(int(match.group(2)))) if match else None
    if rule == "suffix_token":
        match = re.fullmatch(r"(.+?)[\s_.-]+([^\s_.-]+)", value)
        return (match.group(1), match.group(2)) if match else None
    if rule == "prefix_token":
        match = re.fullmatch(r"([^\s_.-]+)[\s_.-]+(.+)", value)
        return (match.group(2), match.group(1)) if match else None
    if len(value) <= width:
        return None
    return (
        (value[:-width].strip(" _-"), value[-width:].strip(" _-"))
        if rule == "suffix_chars"
        else (value[width:].strip(" _-"), value[:width].strip(" _-"))
    )


def grid_score(splits, bonus):
    groups = defaultdict(list)
    for label, (group, member) in splits.items():
        group_key, member_key = norm(group), norm(member)
        # For structural mappability, coded names such as F1, M3, or even
        # numeric group IDs are valid. Biological interpretability is not
        # required; only a non-empty, consistent grouping structure is needed.
        if not group_key or not member_key or group_key == member_key or len(member_key) > 8:
            return None
        groups[group_key].append(member_key)
    if len(groups) < 2 or any(len(x) < 2 or len(x) != len(set(x)) for x in groups.values()):
        return None
    sets = [set(x) for x in groups.values()]
    overlaps = [
        len(a & b) / len(a | b)
        for i, a in enumerate(sets)
        for b in sets[i + 1 :]
    ]
    if not overlaps or min(overlaps) < 0.75:
        return None
    sizes = [len(x) for x in sets]
    if min(sizes) / max(sizes) < 0.5:
        return None
    return 50 * sum(overlaps) / len(overlaps) + 15 * min(sizes) / max(sizes) + 2 * min(len(groups), 6) + bonus


def infer_grid(labels):
    values = remove_run_prefix(labels)
    # Explicit replicate markers receive a strong bonus so ex1/rep1/R1 is not
    # inverted into a group merely because the design is mathematically square.
    rules = [("marked", 0, 40), ("suffix_token", 0, 30), ("prefix_token", 0, 15)]
    rules += [("suffix_chars", n, 10) for n in range(1, 7)]
    rules += [("prefix_chars", n, 5) for n in range(1, 7)]
    best = None
    for rule, width, bonus in rules:
        splits = {label: split_rule(value, rule, width) for label, value in values.items()}
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


def structured_map(documents):
    output = defaultdict(list)
    for name, content in documents:
        lines = [line for line in content.splitlines() if line.strip()][:5001]
        if len(lines) < 2:
            continue
        delimiter = max(("\t", ",", ";"), key=lambda x: lines[0].count(x))
        rows = list(csv.reader(lines, delimiter=delimiter))
        heads = [re.sub(r"\W", "", x.lower()) for x in rows[0]]
        keys = [i for i, x in enumerate(heads) if any(w in x for w in ("datafile", "rawfile", "filename", "assayname", "samplename")) or x == "name"]
        factors = [i for i, x in enumerate(heads) if any(w in x for w in ("factorvalue", "condition", "group", "treatment", "experiment", "phenotype", "genotype", "dose", "time"))]
        if not keys or not factors:
            continue
        width = len(rows[0])
        for row in rows[1:]:
            row += [""] * (width - len(row))
            values = [f"{rows[0][i]}={row[i]}" for i in factors if row[i].strip()]
            for i in keys:
                if norm(row[i]) and values:
                    output[norm(row[i])].append(("; ".join(values), name))
    return output


def unique_match(label, mapping):
    key = norm(label)
    if key in mapping:
        return mapping[key]
    matches = [value for candidate, value in mapping.items() if min(len(key), len(candidate)) >= 5 and (key in candidate or candidate in key)]
    return matches[0] if len(matches) == 1 else []


def informative(condition, label):
    biological = ("condition", "group", "treatment", "factor", "phenotype", "genotype", "dose", "time")
    for part in condition.split(";"):
        head, _, value = part.partition("=")
        if any(word in norm(head) for word in biological) or norm(value) != norm(label):
            return True
    return False


def analyze(accession, max_zip):
    metadata = get_json(f"{API}/projects/{accession}")
    files = project_files(accession)
    documents, sources = [], []

    for item in sorted(files, key=lambda x: -metadata_score(x["name"])):
        if is_pg(item["name"]):
            sources.append((item["name"], header(prefix(item["url"], 8 * 1024**2))))
        elif metadata_score(item["name"]) and (item["size"] or 0) <= 20 * 1024**2:
            documents.append((item["name"], text(prefix(item["url"], 20 * 1024**2))))

    for item in files:
        if not base(item["name"]).endswith(".zip"):
            continue
        try:
            with remote_zip(item, max_zip) as archive:
                for member in archive.infolist():
                    if member.is_dir():
                        continue
                    if is_pg(member.filename):
                        with archive.open(member) as handle:
                            sources.append((f"{item['name']}::{member.filename}", header(handle.read(8 * 1024**2))))
                    elif metadata_score(member.filename) and member.file_size <= 20 * 1024**2:
                        with archive.open(member) as handle:
                            documents.append((f"{item['name']}::{member.filename}", text(handle.read())))
        except Exception as exc:
            print(f"[{accession}] ZIP skipped: {item['name']}: {exc}")

    design = structured_map(documents)
    raw = {norm(item["name"]): [item["name"]] for item in files if base(item["name"]).endswith((".raw", ".wiff", ".mzml", ".mzxml", ".d"))}
    corpus = " ".join(re.findall(r"[a-z0-9]+", (json.dumps(metadata) + " ".join(x[1] for x in documents)).lower()))
    rows = []

    if not sources:
        return [{"accession": accession, "protein_groups_path": "", "lfq_count": 0, "status": "no", "mapped_count": 0, "lfq_names": "[]", "groups": "{}", "members": "{}", "raw_matches": "{}", "reason": "No proteinGroups text file found"}]

    for path, columns in sources:
        pairs = lfq_columns(columns)
        labels = [label for _, label in pairs]
        grid_groups, grid_members = infer_grid(labels)
        groups, members, raw_matches, mapped = {}, {}, {}, 0
        for label in labels:
            hits = [hit for hit in unique_match(label, design) if informative(hit[0], label)]
            raw_hit = unique_match(label, raw)
            if raw_hit:
                raw_matches[label] = raw_hit[0]
            if len({hit[0] for hit in hits}) == 1:
                groups[label] = hits[0][0]
                mapped += 1
            elif label in grid_groups:
                groups[label], members[label] = grid_groups[label], grid_members[label]
                mapped += 1
            elif useful_group(label) and " ".join(re.findall(r"[a-z0-9]+", label.lower())) in corpus:
                groups[label] = label
                mapped += 1

        total = len(labels)
        status = "yes" if total and mapped == total else "partial" if mapped else "uncertain" if raw_matches else "no"
        reason = (
            f"{mapped}/{total} LFQ names mapped to groups"
            if total
            else "No LFQ intensity columns found"
        )
        rows.append({
            "accession": accession,
            "protein_groups_path": path,
            "lfq_count": total,
            "status": status,
            "mapped_count": mapped,
            "lfq_names": json.dumps(labels),
            "groups": json.dumps(groups),
            "members": json.dumps(members),
            "raw_matches": json.dumps(raw_matches),
            "reason": reason,
        })
    return rows


def self_test():
    compact = [f"{group}{member}" for group in ("Deltacys", "Rev150", "Rev190", "Rev31", "Rev362", "Rev622") for member in "ABCD"]
    groups, members = infer_grid(compact)
    assert len(groups) == 24 and set(groups.values()) == {"Deltacys", "Rev150", "Rev190", "Rev31", "Rev362", "Rev622"}
    assert set(members.values()) == set("ABCD")
    earlier = [f"{65 + 2*i}rok_shA_ex{i+1}" for i in range(5)] + [f"{66 + 2*i}rok_shCtrl_ex{i+1}" for i in range(5)]
    groups, members = infer_grid(earlier)
    assert set(groups.values()) == {"shA", "shCtrl"} and set(members.values()) == {"1", "2", "3", "4", "5"}
    coded = [
        f"{group}_{member}"
        for group in ("F1", "F2", "F3", "F4", "F5", "M1", "M2", "M3", "M4", "M5", "Top6_Top20")
        for member in ("1", "2", "3")
    ]
    groups, members = infer_grid(coded)
    assert len(groups) == 33
    assert set(groups.values()) == {"F1", "F2", "F3", "F4", "F5", "M1", "M2", "M3", "M4", "M5", "Top6_Top20"}
    assert set(members.values()) == {"1", "2", "3"}
    # Prefer the conventional final-token member boundary even when the
    # transposed mathematical grid would also be possible.
    many_members = [f"{group}_{member}" for group in ("Control", "Drug") for member in range(1, 11)]
    groups, members = infer_grid(many_members)
    assert set(groups.values()) == {"Control", "Drug"}
    assert set(members.values()) == {str(i) for i in range(1, 11)}
    assert infer_grid(["Alpha", "Beta", "Gamma", "Delta"]) == ({}, {})
    print("Self-test passed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accessions", nargs="*")
    parser.add_argument("--output", type=Path, default=OUTPUT_FILE)
    parser.add_argument("--max-full-zip-gb", type=float, default=MAX_FULL_ZIP_GB)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    accessions = [x.upper() for x in (args.accessions or PRIDE_ACCESSIONS)]
    if not accessions:
        parser.error("Supply accessions or edit PRIDE_ACCESSIONS")
    rows = []
    for accession in accessions:
        print(f"Checking {accession}...")
        rows.extend(analyze(accession, int(args.max_full_zip_gb * 1024**3)))
    with args.output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {args.output.resolve()}")


if __name__ == "__main__":
    main()
