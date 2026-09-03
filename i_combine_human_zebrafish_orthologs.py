#!/usr/bin/env python3
"""Combine human-to-zebrafish ortholog calls from NCBI, Ensembl, and ZFIN."""
'''
python i_combine_human_zebrafish_orthologs.py   --genes /Users/mehman/Desktop/Test/hyperglycemia_features_of_in
terest.txt   --ncbi /Users/mehman/Projects/PoC_data_processing/Orthology/ncbi_hyperglycemia.tsv   --ensembl /User
s/mehman/Projects/PoC_data_processing/Orthology/ensembl_hyperglycemia.csv   --zfin /Users/mehman/Projects/PoC_dat
a_processing/Orthology/human_orthos_zfin.txt   --output combined_orthologs.csv
'''
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


DEFAULT_GENES = "hyperglycemia_features_of_interest.txt"
DEFAULT_NCBI = "zebrafish_ncbi.tsv"
DEFAULT_ENSEMBL = "ensembl.csv"
DEFAULT_ZFIN = "human_orthos_zfin.txt"
DEFAULT_OUTPUT = "human_to_zebrafish_orthologs_combined.csv"

# Positional layout of the standard ZFIN human orthology download.
ZFIN_COLUMNS = [
    "ZFIN_gene_ID",
    "Zebrafish_gene",
    "Zebrafish_gene_name",
    "Human_gene",
    "Human_gene_name",
    "Human_OMIM_ID",
    "Human_NCBI_ID",
    "Human_HGNC_ID",
    "Evidence_code",
    "Publication_ID",
    "Evidence_method",
    "ECO_ID",
    "Evidence_description",
]


def clean_header(value: object) -> str:
    """Normalize a column name while retaining its words."""
    text = str(value).replace("\ufeff", "").strip().casefold()
    return re.sub(r"[^a-z0-9]+", " ", text).strip()


def find_column(frame: pd.DataFrame, candidates: list[str], source: str) -> str:
    """Find a column using normalized exact names, then distinctive phrases."""
    normalized = {clean_header(column): column for column in frame.columns}

    for candidate in candidates:
        key = clean_header(candidate)
        if key in normalized:
            return normalized[key]

    for candidate in candidates:
        key = clean_header(candidate)
        partial = [original for norm, original in normalized.items() if key in norm]
        if len(partial) == 1:
            return partial[0]

    available = ", ".join(map(str, frame.columns))
    raise ValueError(
        f"Could not identify the required column in {source}. "
        f"Available columns: {available}"
    )


def read_table(path: Path, *, separator: str | None = None) -> pd.DataFrame:
    """Read a delimited table entirely as text."""
    kwargs = {
        "dtype": str,
        "keep_default_na": False,
        "encoding": "utf-8-sig",
    }
    if separator is None:
        return pd.read_csv(path, sep=None, engine="python", **kwargs)
    return pd.read_csv(path, sep=separator, **kwargs)


def normalize_gene(value: object) -> str:
    """Create a case-insensitive gene-symbol key."""
    return str(value).strip().casefold()


def clean_gene(value: object) -> str:
    """Clean a gene symbol for display."""
    value = str(value).strip()
    return "" if value.casefold() in {"", "nan", "none", "na", "n/a", "-"} else value


def unique_symbols(values: pd.Series) -> list[str]:
    """Return nonempty symbols, deduplicated case-insensitively and sorted."""
    symbols: dict[str, str] = {}
    for raw_value in values:
        # Usually there is one symbol per field, but this also tolerates lists.
        for value in re.split(r"[;|]", clean_gene(raw_value)):
            symbol = clean_gene(value)
            if symbol:
                symbols.setdefault(normalize_gene(symbol), symbol)
    return sorted(symbols.values(), key=lambda value: (value.casefold(), value))


def build_mapping(
    frame: pd.DataFrame,
    human_column: str,
    zebrafish_column: str,
) -> dict[str, list[str]]:
    """Map each normalized human symbol to its unique zebrafish symbols."""
    selected = frame[[human_column, zebrafish_column]].copy()
    selected["_human_key"] = selected[human_column].map(normalize_gene)
    selected = selected[selected["_human_key"] != ""]

    return {
        human_key: unique_symbols(group[zebrafish_column])
        for human_key, group in selected.groupby("_human_key", sort=False)
    }


def read_input_genes(path: Path) -> list[str]:
    """Read one human gene symbol per line, preserving input order."""
    genes: list[str] = []
    seen: set[str] = set()

    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            # Accept a plain list and tolerate a one-column CSV/TSV-like file.
            gene = clean_gene(re.split(r"[\t,]", line.strip(), maxsplit=1)[0])
            key = normalize_gene(gene)
            if not gene or key in {"gene", "gene symbol", "human gene", "human_gene"}:
                continue
            if key not in seen:
                genes.append(gene)
                seen.add(key)

    if not genes:
        raise ValueError(f"No gene symbols were found in {path}")
    return genes


def read_ncbi_mapping(path: Path) -> dict[str, list[str]]:
    frame = read_table(path, separator="\t")
    # Accept both the earlier NCBI table and the detailed table produced by
    # i_human_to_zebrafish_local.py. In the detailed table, human_input is the
    # original symbol from the requested list, so it is the preferred join key.
    human = find_column(
        frame,
        ["Human_gene", "Human gene", "human_input", "human_symbol"],
        path.name,
    )
    zebrafish = find_column(
        frame,
        ["Zebrafish_gene", "Zebrafish gene", "zebrafish_symbol"],
        path.name,
    )
    return build_mapping(frame, human, zebrafish)


def read_ensembl_mapping(path: Path) -> dict[str, list[str]]:
    frame = read_table(path)
    human = find_column(frame, ["Gene name"], path.name)
    zebrafish = find_column(frame, ["Zebrafish gene name"], path.name)
    confidence_column = find_column(
        frame,
        [
            "Zebrafish orthology confidence",
            "drerio homolog orthology confidence",
        ],
        path.name,
    )

    # Ensembl represents low confidence as 0 and high confidence as 1.
    # Missing or nonnumeric confidence values are excluded as well.
    confidence = pd.to_numeric(
        frame[confidence_column].astype(str).str.strip(),
        errors="coerce",
    )
    frame = frame.loc[confidence.notna() & confidence.ne(0)]

    return build_mapping(frame, human, zebrafish)


def read_zfin_mapping(path: Path) -> dict[str, list[str]]:
    # The standard ZFIN file has no header. Reading positionally also avoids
    # dependence on small header-name differences across ZFIN downloads.
    frame = pd.read_csv(
        path,
        sep=None,
        engine="python",
        header=None,
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig",
    )
    if frame.shape[1] < 4:
        raise ValueError(
            f"{path.name} has {frame.shape[1]} columns; at least 4 were expected. "
            "Check that it is a comma- or tab-separated ZFIN orthology file."
        )

    if frame.shape[1] > len(ZFIN_COLUMNS):
        frame = frame.iloc[:, : len(ZFIN_COLUMNS)]
    frame.columns = ZFIN_COLUMNS[: frame.shape[1]]

    # If the file happens to include a header row, it will not match real input
    # genes and is harmless; this explicit filter removes it for completeness.
    header_like = frame["Human_gene"].map(clean_header).isin(
        {"human gene", "human gene symbol", "human symbol"}
    )
    frame = frame.loc[~header_like]
    return build_mapping(frame, "Human_gene", "Zebrafish_gene")


def combine_orthologs(*source_lists: list[str]) -> list[str]:
    combined: dict[str, str] = {}
    for source_list in source_lists:
        for symbol in source_list:
            combined.setdefault(normalize_gene(symbol), symbol)
    return sorted(combined.values(), key=lambda value: (value.casefold(), value))


def join_symbols(symbols: list[str]) -> str:
    return "; ".join(symbols)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find zebrafish orthologs for a human gene list using NCBI, "
            "Ensembl, and ZFIN mapping files."
        )
    )
    parser.add_argument("--genes", type=Path, default=Path(DEFAULT_GENES))
    parser.add_argument("--ncbi", type=Path, default=Path(DEFAULT_NCBI))
    parser.add_argument("--ensembl", type=Path, default=Path(DEFAULT_ENSEMBL))
    parser.add_argument("--zfin", type=Path, default=Path(DEFAULT_ZFIN))
    parser.add_argument("--output", type=Path, default=Path(DEFAULT_OUTPUT))
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    for path in (args.genes, args.ncbi, args.ensembl, args.zfin):
        if not path.is_file():
            raise FileNotFoundError(f"Input file not found: {path}")

    input_genes = read_input_genes(args.genes)
    ncbi_map = read_ncbi_mapping(args.ncbi)
    ensembl_map = read_ensembl_mapping(args.ensembl)
    zfin_map = read_zfin_mapping(args.zfin)

    rows: list[dict[str, object]] = []
    for human_gene in input_genes:
        key = normalize_gene(human_gene)
        ncbi = ncbi_map.get(key, [])
        ensembl = ensembl_map.get(key, [])
        zfin = zfin_map.get(key, [])
        combined = combine_orthologs(ncbi, ensembl, zfin)

        rows.append(
            {
                "Human_gene": human_gene,
                "NCBI_orthologs": join_symbols(ncbi),
                "Ensembl_orthologs": join_symbols(ensembl),
                "ZFIN_orthologs": join_symbols(zfin),
                "All_unique_orthologs": join_symbols(combined),
            }
        )

    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)

    mapped = int(output["All_unique_orthologs"].ne("").sum())
    print(f"Saved {len(output)} human genes to: {args.output}")
    print(f"Genes with at least one zebrafish ortholog: {mapped}/{len(output)}")


if __name__ == "__main__":
    main()
