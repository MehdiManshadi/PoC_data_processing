#!/usr/bin/env python3
'''
python h_ncbi_human_to_zebrafish_local.py   mitochondrial_myopathy_features_of_interest.txt
'''
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


HUMAN_TAX_ID = "9606"
ZEBRAFISH_TAX_ID = "7955"

ORTHOLOG_FILE = Path("gene_orthologs.gz")
HUMAN_INFO_FILE = Path("Homo_sapiens.gene_info.gz")
ZEBRAFISH_INFO_FILE = Path("Danio_rerio.gene_info.gz")


def read_gene_list(path: Path) -> list[str]:
    """Read one human gene symbol per line."""
    genes = []

    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            gene = line.strip().strip("\"'")
            if gene and not gene.startswith("#"):
                genes.append(gene)

    # Remove duplicates while preserving the original order.
    return list(dict.fromkeys(genes))


def load_gene_info(path: Path) -> pd.DataFrame:
    """Read the required columns from an NCBI gene_info file."""
    return pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        dtype=str,
        usecols=["GeneID", "Symbol", "Synonyms"],
    )


def map_input_symbols(
    genes: list[str],
    human_info: pd.DataFrame,
) -> pd.DataFrame:
    """Map input symbols or synonyms to human NCBI GeneIDs."""

    official_lookup = defaultdict(list)
    synonym_lookup = defaultdict(list)
    id_to_symbol = {}

    for row in human_info.itertuples(index=False):
        gene_id = row.GeneID
        symbol = row.Symbol

        id_to_symbol[gene_id] = symbol
        official_lookup[symbol.casefold()].append(gene_id)

        if pd.notna(row.Synonyms) and row.Synonyms != "-":
            for synonym in row.Synonyms.split("|"):
                synonym = synonym.strip()
                if synonym:
                    synonym_lookup[synonym.casefold()].append(gene_id)

    records = []

    for gene in genes:
        key = gene.casefold()

        gene_ids = official_lookup.get(key, [])
        match_type = "official_symbol"

        if not gene_ids:
            gene_ids = synonym_lookup.get(key, [])
            match_type = "synonym"

        gene_ids = sorted(set(gene_ids), key=int)

        if not gene_ids:
            records.append(
                {
                    "human_input": gene,
                    "human_symbol": None,
                    "human_gene_id": None,
                    "input_match": "not_found",
                }
            )
            continue

        for gene_id in gene_ids:
            records.append(
                {
                    "human_input": gene,
                    "human_symbol": id_to_symbol[gene_id],
                    "human_gene_id": gene_id,
                    "input_match": match_type,
                }
            )

    return pd.DataFrame(records)


def find_zebrafish_orthologs(
    selected_human_ids: set[str],
) -> pd.DataFrame:
    """Scan gene_orthologs.gz once and retain human-zebrafish pairs."""

    columns = [
        "#tax_id",
        "GeneID",
        "relationship",
        "Other_tax_id",
        "Other_GeneID",
    ]

    results = []

    for chunk in pd.read_csv(
        ORTHOLOG_FILE,
        sep="\t",
        compression="gzip",
        dtype=str,
        usecols=columns,
        chunksize=500_000,
    ):
        # Human -> zebrafish orientation
        forward = chunk.loc[
            (chunk["#tax_id"] == HUMAN_TAX_ID)
            & (chunk["Other_tax_id"] == ZEBRAFISH_TAX_ID)
            & (chunk["GeneID"].isin(selected_human_ids)),
            ["GeneID", "Other_GeneID", "relationship"],
        ].rename(
            columns={
                "GeneID": "human_gene_id",
                "Other_GeneID": "zebrafish_gene_id",
            }
        )

        # Also check the reverse orientation, if present.
        reverse = chunk.loc[
            (chunk["#tax_id"] == ZEBRAFISH_TAX_ID)
            & (chunk["Other_tax_id"] == HUMAN_TAX_ID)
            & (chunk["Other_GeneID"].isin(selected_human_ids)),
            ["Other_GeneID", "GeneID", "relationship"],
        ].rename(
            columns={
                "Other_GeneID": "human_gene_id",
                "GeneID": "zebrafish_gene_id",
            }
        )

        if not forward.empty:
            results.append(forward)

        if not reverse.empty:
            results.append(reverse)

    if not results:
        return pd.DataFrame(
            columns=[
                "human_gene_id",
                "zebrafish_gene_id",
                "relationship",
            ]
        )

    return (
        pd.concat(results, ignore_index=True)
        .drop_duplicates()
        .reset_index(drop=True)
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Map human gene symbols to zebrafish NCBI orthologs."
    )
    parser.add_argument(
        "human_gene_list",
        type=Path,
        help="Text file containing one human gene symbol per line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("zebrafish_ncbi.tsv"),
        help="Output TSV file (default: zebrafish_ncbi.tsv).",
    )
    args = parser.parse_args()

    required_files = [
        args.human_gene_list,
        ORTHOLOG_FILE,
        HUMAN_INFO_FILE,
        ZEBRAFISH_INFO_FILE,
    ]

    for path in required_files:
        if not path.exists():
            raise FileNotFoundError(f"Required file not found: {path}")

    genes = read_gene_list(args.human_gene_list)

    print(f"Loaded {len(genes)} unique human gene symbols.")
    print("Reading human gene information...")

    human_info = load_gene_info(HUMAN_INFO_FILE)
    human_matches = map_input_symbols(genes, human_info)

    selected_human_ids = set(
        human_matches["human_gene_id"].dropna().astype(str)
    )

    print(f"Matched {len(selected_human_ids)} human NCBI GeneIDs.")
    print("Scanning local NCBI ortholog table...")

    orthologs = find_zebrafish_orthologs(selected_human_ids)

    print("Reading zebrafish gene information...")

    zebrafish_info = load_gene_info(ZEBRAFISH_INFO_FILE)

    zebrafish_symbol_lookup = dict(
        zip(zebrafish_info["GeneID"], zebrafish_info["Symbol"])
    )

    result = human_matches.merge(
        orthologs,
        how="left",
        on="human_gene_id",
    )

    result["zebrafish_symbol"] = result["zebrafish_gene_id"].map(
        zebrafish_symbol_lookup
    )

    result["status"] = "mapped"

    result.loc[
        result["human_gene_id"].isna(),
        "status",
    ] = "human_gene_not_found"

    result.loc[
        result["human_gene_id"].notna()
        & result["zebrafish_gene_id"].isna(),
        "status",
    ] = "no_zebrafish_ortholog"

    result.loc[
        result["zebrafish_gene_id"].notna()
        & result["zebrafish_symbol"].isna(),
        "status",
    ] = "zebrafish_symbol_not_found"

    ortholog_counts = (
        result.groupby("human_input")["zebrafish_gene_id"]
        .nunique()
        .to_dict()
    )

    result["n_zebrafish_orthologs"] = (
        result["human_input"].map(ortholog_counts).fillna(0).astype(int)
    )

    input_order = {gene: position for position, gene in enumerate(genes)}
    result["_input_order"] = result["human_input"].map(input_order)

    result = result.sort_values(
        ["_input_order", "human_gene_id", "zebrafish_gene_id"],
        na_position="last",
    ).drop(columns="_input_order")

    result = result[
        [
            "human_input",
            "human_symbol",
            "human_gene_id",
            "input_match",
            "zebrafish_symbol",
            "zebrafish_gene_id",
            "relationship",
            "n_zebrafish_orthologs",
            "status",
        ]
    ]

    result.to_csv(
        args.output,
        sep="\t",
        index=False,
        na_rep="",
    )

    mapped_human_genes = result.loc[
        result["zebrafish_gene_id"].notna(),
        "human_input",
    ].nunique()

    ortholog_pair_count = result.loc[
        result["zebrafish_gene_id"].notna(),
        ["human_gene_id", "zebrafish_gene_id"],
    ].drop_duplicates().shape[0]

    print(f"Human genes with zebrafish orthologs: {mapped_human_genes}")
    print(f"Human-zebrafish ortholog pairs: {ortholog_pair_count}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()