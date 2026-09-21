#!/usr/bin/env python3
"""Aggregate replicate intensities into a gene-by-tissue TSV matrix."""

from pathlib import Path
import argparse

import pandas as pd


REQUIRED_COLUMNS = {"Gene", "Gene name", "Tissue", "Intensity"}


def build_gene_tissue_matrix(
    input_file: Path,
    output_file: Path,
    presence_threshold: float = 0.10,
) -> None:
    df = pd.read_csv(input_file, sep="\t", low_memory=False)

    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    # An unreported intensity is defined as zero in this dataset. Thus the
    # source matrix contains numeric intensities or zeros, never missing data.
    df["Intensity"] = pd.to_numeric(df["Intensity"], errors="coerce").fillna(0.0)
    df["Tissue"] = df["Tissue"].astype("string").str.strip()

    gene_keys = ["Gene", "Gene name"]

    # Apply the presence rule once per gene across ALL tissue samples together.
    # A gene passes only when more than 10% of all its intensity rows are
    # measured and non-zero.
    detected = df["Intensity"].notna() & df["Intensity"].ne(0)
    detected_fraction = detected.groupby(
        [df[column] for column in gene_keys],
        dropna=False,
        observed=True,
    ).transform("mean")
    gene_passes = detected_fraction > presence_threshold

    # Calculate each tissue mean from its replicates. The 10% rule is not
    # applied again at tissue level. Genes that fail the global rule remain in
    # the output, but every tissue value for those genes is set to NA.
    stats = (
        df.assign(_gene_passes=gene_passes)
        .groupby(gene_keys + ["Tissue"], dropna=False, observed=True)
        .agg(
            gene_passes=("_gene_passes", "first"),
            aggregated_intensity=("Intensity", "mean"),
        )
        .reset_index()
    )
    stats.loc[
        ~stats["gene_passes"],
        "aggregated_intensity",
    ] = float("nan")

    # Each gene-tissue combination is already unique after groupby, so pivot
    # directly. This avoids the very large Cartesian product that
    # pivot_table(dropna=False) can create for Gene × Gene-name combinations.
    result = stats.pivot(
            index=gene_keys,
        columns="Tissue",
        values="aggregated_intensity",
    ).reset_index()
    result.columns.name = None
    result = result.rename(
        columns={
            column: f"H1040_HPA_{'_'.join(str(column).split())}"
            for column in result.columns
            if column not in gene_keys
        }
    )
    # Write missing or criterion-failed results explicitly as NA.
    result.to_csv(output_file, sep="\t", index=False, na_rep="NA")

    print(f"Saved {result.shape[0]} genes × {result.shape[1] - 2} tissues to {output_file}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a gene-by-tissue matrix from replicate-level intensities."
    )
    parser.add_argument("input_tsv", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("--threshold", type=float, default=0.10)
    args = parser.parse_args()

    build_gene_tissue_matrix(args.input_tsv, args.output_tsv, args.threshold)


if __name__ == "__main__":
    main()
