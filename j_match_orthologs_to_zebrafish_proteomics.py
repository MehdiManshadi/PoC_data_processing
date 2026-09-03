#!/usr/bin/env python3
"""Find human genes whose zebrafish orthologs occur in a proteomics table."""
'''
python j_match_orthologs_to_zebrafish_proteomics.py --orthologs "/Users/mehman/Projects/PoC_data_processing/Orthology/combined_orthologs_mitochondrial_myopathy.csv" --proteomics "/Users/mehman/3-IntellAif/zebrafish POC/26031
3_Vaparanta_Proteins (1).xlsx"
'''
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


DEFAULT_ORTHOLOGS = "combined_orthologs_hyperglycemia.csv"
DEFAULT_PROTEOMICS = "260313_Vaparanta_Proteins (1).xlsx"
DEFAULT_OUTPUT = "human_genes_with_detected_zebrafish_orthologs.csv"

HUMAN_COLUMN = "Human_gene"
ORTHOLOG_COLUMN = "All_unique_orthologs"
PROTEOMICS_GENE_COLUMN = "PG.Genes"
PROTEIN_GROUP_COLUMN = "PG.ProteinGroups"

MISSING_VALUES = {"", "-", "na", "n/a", "nan", "none", "null"}


def clean_text(value: object) -> str:
    """Convert a table value to clean text, treating common markers as empty."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in MISSING_VALUES else text


def normalize_gene(value: object) -> str:
    """Create a case-insensitive exact-match key for a gene symbol."""
    return clean_text(value).casefold()


def split_gene_symbols(value: object) -> list[str]:
    """Split a field containing one or more gene symbols."""
    symbols: list[str] = []
    seen: set[str] = set()

    for part in re.split(r"[;|,]", clean_text(value)):
        symbol = clean_text(part)
        key = normalize_gene(symbol)
        if symbol and key not in seen:
            symbols.append(symbol)
            seen.add(key)

    return symbols


def join_values(values: list[str] | set[str], separator: str = "; ") -> str:
    """Join unique nonempty values in a stable, case-insensitive order."""
    unique = {normalize_gene(value): clean_text(value) for value in values if clean_text(value)}
    return separator.join(sorted(unique.values(), key=lambda value: (value.casefold(), value)))


def read_table(path: Path, sheet: str | int = 0) -> pd.DataFrame:
    """Read an Excel, CSV, or TSV table as text where possible."""
    suffixes = [suffix.casefold() for suffix in path.suffixes]

    if path.suffix.casefold() in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=sheet, dtype=str, keep_default_na=False)

    if ".tsv" in suffixes or path.suffix.casefold() in {".txt", ".tab"}:
        return pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            encoding="utf-8-sig",
        )

    return pd.read_csv(
        path,
        sep=None,
        engine="python",
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig",
    )


def resolve_column(frame: pd.DataFrame, required_name: str, source: Path) -> str:
    """Resolve a required column after trimming whitespace and any BOM."""
    cleaned = {
        str(column).replace("\ufeff", "").strip().casefold(): column
        for column in frame.columns
    }
    key = required_name.casefold()

    if key not in cleaned:
        available = ", ".join(map(str, frame.columns))
        raise ValueError(
            f"Required column {required_name!r} was not found in {source.name}. "
            f"Available columns: {available}"
        )

    return cleaned[key]


def optional_column(frame: pd.DataFrame, name: str) -> str | None:
    """Resolve an optional column after normalizing its header."""
    target = name.casefold()
    for column in frame.columns:
        if str(column).replace("\ufeff", "").strip().casefold() == target:
            return column
    return None


def build_proteomics_index(
    proteomics: pd.DataFrame,
    gene_column: str,
    protein_group_column: str | None,
) -> dict[str, dict[str, set]]:
    """Index every zebrafish gene symbol and its supporting protein rows."""
    index: dict[str, dict[str, set]] = defaultdict(
        lambda: {
            "symbols": set(),
            "row_numbers": set(),
            "protein_groups": set(),
        }
    )

    # Row 1 is the spreadsheet header, so data rows begin at row 2.
    for spreadsheet_row, (_, row) in enumerate(proteomics.iterrows(), start=2):
        protein_group = (
            clean_text(row[protein_group_column]) if protein_group_column else ""
        )

        for symbol in split_gene_symbols(row[gene_column]):
            key = normalize_gene(symbol)
            index[key]["symbols"].add(symbol)
            index[key]["row_numbers"].add(spreadsheet_row)
            if protein_group:
                index[key]["protein_groups"].add(protein_group)

    return dict(index)


def build_report(
    orthologs: pd.DataFrame,
    ortholog_column: str,
    proteomics_index: dict[str, dict[str, set]],
) -> pd.DataFrame:
    """Keep human genes with at least one detected zebrafish ortholog."""
    rows: list[dict[str, object]] = []

    for _, source_row in orthologs.iterrows():
        detected_symbols: list[str] = []
        matched_rows: set[int] = set()
        matched_protein_groups: set[str] = set()

        for ortholog in split_gene_symbols(source_row[ortholog_column]):
            match = proteomics_index.get(normalize_gene(ortholog))
            if match is None:
                continue

            detected_symbols.append(ortholog)
            matched_rows.update(match["row_numbers"])
            matched_protein_groups.update(match["protein_groups"])

        if not detected_symbols:
            continue

        output_row = source_row.to_dict()
        output_row["Detected_zebrafish_orthologs"] = join_values(detected_symbols)
        output_row["N_detected_zebrafish_orthologs"] = len(
            {normalize_gene(symbol) for symbol in detected_symbols}
        )
        output_row["Matched_proteomics_row_count"] = len(matched_rows)
        output_row["Matched_proteomics_rows"] = "; ".join(
            str(row_number) for row_number in sorted(matched_rows)
        )
        output_row["Matched_PG.ProteinGroups"] = join_values(
            matched_protein_groups,
            separator=" | ",
        )
        rows.append(output_row)

    extra_columns = [
        "Detected_zebrafish_orthologs",
        "N_detected_zebrafish_orthologs",
        "Matched_proteomics_row_count",
        "Matched_proteomics_rows",
        "Matched_PG.ProteinGroups",
    ]
    return pd.DataFrame(rows, columns=[*orthologs.columns, *extra_columns])


def parse_sheet(value: str) -> str | int:
    """Allow --sheet to be either an Excel sheet name or a zero-based index."""
    return int(value) if value.isdigit() else value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Report human genes having at least one zebrafish ortholog in "
            "the PG.Genes column of a proteomics file."
        )
    )
    parser.add_argument(
        "--orthologs",
        type=Path,
        default=Path(DEFAULT_ORTHOLOGS),
        help=f"Combined human-zebrafish mapping (default: {DEFAULT_ORTHOLOGS})",
    )
    parser.add_argument(
        "--proteomics",
        type=Path,
        default=Path(DEFAULT_PROTEOMICS),
        help=f"Zebrafish proteomics Excel/CSV/TSV file (default: {DEFAULT_PROTEOMICS})",
    )
    parser.add_argument(
        "--sheet",
        type=parse_sheet,
        default=0,
        help="Excel sheet name or zero-based sheet index (default: first sheet)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(DEFAULT_OUTPUT),
        help=f"Output CSV file (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    for path in (args.orthologs, args.proteomics):
        if not path.is_file():
            raise FileNotFoundError(f"Input file not found: {path}")

    orthologs = read_table(args.orthologs)
    proteomics = read_table(args.proteomics, sheet=args.sheet)

    resolve_column(orthologs, HUMAN_COLUMN, args.orthologs)
    ortholog_column = resolve_column(
        orthologs,
        ORTHOLOG_COLUMN,
        args.orthologs,
    )
    proteomics_gene_column = resolve_column(
        proteomics,
        PROTEOMICS_GENE_COLUMN,
        args.proteomics,
    )
    protein_group_column = optional_column(proteomics, PROTEIN_GROUP_COLUMN)

    proteomics_index = build_proteomics_index(
        proteomics,
        proteomics_gene_column,
        protein_group_column,
    )
    report = build_report(
        orthologs,
        ortholog_column,
        proteomics_index,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.output, index=False)

    print(f"Human genes in ortholog table: {len(orthologs)}")
    print(f"Unique zebrafish genes in proteomics table: {len(proteomics_index)}")
    print(f"Human genes with at least one detected ortholog: {len(report)}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
