from pathlib import Path

import numpy as np
import pandas as pd


def collect_features(
    features_file, data_folder, mapping_file, ortholog_files, output_file=None
):
    """Process all CSVs in data_folder and save their columns in one CSV.

    Rows follow features_file order. Output includes feature and sample names.
    Returns the combined DataFrame. Files are processed in filename order.
    """
    with open(features_file, encoding="utf-8") as file:
        features = file.read().splitlines()
    gene_protein_map = pd.read_csv(mapping_file, sep="\t", usecols=["From", "Entry"])
    gene_protein_map = gene_protein_map[
        gene_protein_map["From"].isin(features)
    ].drop_duplicates()

    output_file = (
        Path(output_file) if output_file is not None
        else Path(__file__).resolve().parent / "e_processed_features" / "all_processed_features.csv"
    )
    data_files = sorted(
        path for path in Path(data_folder).iterdir()
        if path.is_file() and path.suffix.lower() == ".csv"
        and path.resolve() != output_file.resolve()
    )
    if not data_files:
        raise ValueError(f"No CSV files found in {data_folder}")

    ortholog_maps = {}
    for dataset_type, ortholog_file in ortholog_files.items():
        ortholog_map = pd.read_csv(
            ortholog_file, usecols=["Human_gene", "All_unique_orthologs"]
        )
        ortholog_map = ortholog_map[
            ortholog_map["Human_gene"].isin(features)
        ].assign(
            ortholog=lambda frame: frame["All_unique_orthologs"].str.split(";")
        ).explode("ortholog")
        ortholog_map["ortholog"] = ortholog_map["ortholog"].str.strip()
        ortholog_maps[dataset_type] = ortholog_map

    processed = [
        _process_data_file(features, file_path, gene_protein_map, ortholog_maps)
        for file_path in data_files
    ]
    combined = pd.concat(processed, axis=1)
    combined.insert(0, "feature", features)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, index=False)
    return combined


def _process_data_file(features, file_path, gene_protein_map, ortholog_maps):
    """Normalize one dataset and aggregate mapped proteins to genes by median."""
    prefix = file_path.stem

    if prefix == "H1040_HPA":
        # This one file reports gene-level HPA tissue intensities (tab-
        # separated, "Gene"/"Gene name" columns) rather than protein-group
        # rows, so it is matched straight to "Gene name" instead of going
        # through the UniProt protein-to-gene mapping used for every other
        # dataset. This branch is specific to this file only.
        df = pd.read_csv(file_path, sep="\t")
        matching_columns = [
            col for col in df.columns
            if str(col).startswith(prefix)
        ]
        measurements = df[matching_columns].apply(pd.to_numeric, errors="coerce")
        log2_measurements = np.log2(measurements.where(measurements.ne(0)))
        medians = log2_measurements.median(skipna=True)
        normalized = log2_measurements.subtract(medians, axis="columns").mask(
            measurements.eq(0), 0
        )
        return (
            normalized.assign(**{"Gene name": df["Gene name"]})
            .groupby("Gene name")[matching_columns]
            .median()
            .reindex(features)
            .reset_index(drop=True)
        )

    df = pd.read_csv(file_path)
    # Remove the entire row if any of these quality-control fields contains "+".
    exclusion_columns = ["Only identified by site", "Reverse", "Potential contaminant"]
    df = df[~df[exclusion_columns].eq("+").any(axis=1)].copy()

    # Select columns whose names start with the dataset identifier (filename stem).
    dataset_type = prefix[0]
    uses_orthologs = dataset_type != "H"
    if uses_orthologs:
        selected_ortholog_map = ortholog_maps[dataset_type]
    matching_columns = [
        col for col in df.columns
        if str(col).startswith(prefix)
    ]

    # Keep protein and gene identifiers with all selected measurement columns.
    identifier_columns = ["Protein IDs", "Gene names"]
    if uses_orthologs:
        identifier_columns.append("Fasta headers")
    LFQ = df[identifier_columns + matching_columns].copy()

    # Convert nonnumeric entries to NaN before computing log2 intensities.
    measurements = LFQ[matching_columns].apply(pd.to_numeric, errors="coerce")

    # Temporarily treat original zeros as missing so they do not enter the log2 median.
    log2_measurements = np.log2(measurements.where(measurements.ne(0)))
    # Calculate one median per sample column, ignoring NaN values.
    medians = log2_measurements.median(skipna=True)
    # Subtract each column's median, then restore only the original zero positions.
    # Missing values remain NaN; normalized nonzero measurements may also equal zero.
    LFQ[matching_columns] = log2_measurements.subtract(medians, axis="columns").mask(
        measurements.eq(0), 0
    )
    # Optional inspection: print a protein's normalized measurements or sample medians.
    # print(LFQ.loc[LFQ["Protein IDs"] == "A0AVI4", matching_columns])
    # print(medians)

    if uses_orthologs:
        # Use the union of listed gene names and GN= values in FASTA headers.
        listed_genes = LFQ["Gene names"].fillna("").str.split(";")
        fasta_genes = LFQ["Fasta headers"].fillna("").str.findall(
            r"\bGN=([^\s;]+)"
        )
        all_genes = [list(dict.fromkeys(left + right)) for left, right in zip(
            listed_genes, fasta_genes
        )]

        # Match animal gene names to the orthologs of each requested human gene.
        protein_measurements = LFQ.assign(
            _row=range(len(LFQ)),
            ortholog=all_genes,
        ).explode("ortholog")
        protein_measurements["ortholog"] = protein_measurements["ortholog"].str.strip()
        protein_measurements = protein_measurements.merge(
            selected_ortholog_map[["Human_gene", "ortholog"]],
            on="ortholog",
            how="inner",
        ).drop_duplicates(["_row", "Human_gene"])
        group_column = "Human_gene"
    else:
        # Expand semicolon-separated protein groups so each complete ID can match a
        # feature. Each group member receives that row's normalized measurements.
        protein_measurements = LFQ.assign(
            feature=LFQ["Protein IDs"].astype("string").str.split(";")
        ).explode("feature")
        protein_measurements["feature"] = protein_measurements["feature"].str.strip()
        # UniProt entries use sp|accession|name or tr|accession|name. Extract
        # the accession while preserving plain IDs and isoform suffixes.
        accessions = protein_measurements["feature"].str.split("|").str[1]
        protein_measurements["feature"] = accessions.fillna(
            protein_measurements["feature"]
        ).str.strip()
        protein_measurements = protein_measurements.merge(
            gene_protein_map, left_on="feature", right_on="Entry", how="inner"
        )
        group_column = "From"

    # Aggregate all mapped proteins per gene, preserving the requested gene order.
    # Unmatched genes retain NaN; original zeros are included in the median.
    feature_measurements = (
        protein_measurements.groupby(group_column)[matching_columns]
        .median()
        .reindex(features)
        # Remove feature labels only after establishing the feature-file row order.
        .reset_index(drop=True)
    )

    return feature_measurements


if __name__ == "__main__":
    collect_features(
        features_file=(
            "/Users/mehman/3-IntellAif/zebrafish POC/Integrated_features.txt"
        ),
        mapping_file=(
            "/Users/mehman/Projects/PoC_data_processing/"
            "isoforms/human/features_gene_protein_map.tsv"
        ),
        data_folder=(
            "/Users/mehman/Projects/PoC_data_processing/scripts/"
            "c&d_aggrigated_datasets"
        ),
        ortholog_files={
            "Z": (
                "/Users/mehman/Projects/PoC_data_processing/gene_orthology/"
                "combined_features_zebrafish_orthologs.csv"
            ),
            "M": (
                "/Users/mehman/Projects/PoC_data_processing/gene_orthology/"
                "combined_features_mouse_orthologs.csv"
            ),
        },
    )
