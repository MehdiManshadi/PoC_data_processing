from pathlib import Path

import numpy as np
import pandas as pd


EXCLUSION_COLUMNS = ["Only identified by site", "Reverse", "Potential contaminant"]

NON_SAMPLE_COLUMNS = [
    "Protein IDs", "Majority protein IDs", "Protein names", "Gene names",
    "Fasta headers", "Only identified by site", "Reverse",
    "Potential contaminant", "Source file", "Source row",
]

MISSING_VALUE_POLICIES = ("zero", "min_minus_one", "neg_log_median")


def _missing_value_floor(missing_value_policy, centered_real, medians):
    """Return the floor value used for zero/missing measurements, and for
    genes with no mapped protein at all, under missing_value_policy:

    - "zero": missing/zero measurements are simply 0.
    - "min_minus_one": one unit below the minimum centered value actually
      observed among real (nonzero, non-missing) measurements in that column.
    - "neg_log_median": -log2(median + 1) — the value a raw zero would get
      under ordinary log2(x+1) centering against the median.
    """
    if missing_value_policy == "zero":
        return 0
    if missing_value_policy == "min_minus_one":
        return centered_real.min(skipna=True) - 1
    if missing_value_policy == "neg_log_median":
        return -medians
    raise ValueError(
        f"Unknown missing_value_policy: {missing_value_policy!r}; "
        f"expected one of {MISSING_VALUE_POLICIES}"
    )


def assign_datasets(
    features_file,
    reference_column,
    test_columns=(),
    external_weight_columns=(),
    output_dir=None,
):
    """Assign every sample column in features_file to a model input file.

    reference_column is a single column; external_weight_columns must name
    exactly two columns (group1, group2); test_columns can name any number
    of columns, in the order given. Every sample column not explicitly
    assigned to one of these goes to train_dataset.txt, in features_file's
    original column order. Row (feature) order is preserved from
    features_file throughout. evaluation_dataset.txt is not produced here —
    see generate_pseudo_replicates.

    All outputs are written without row or column labels, semicolon-separated,
    with missing values written as 0 — the format expected by
    sample_script.py's load_semicolon_1d_to_NF. test_dataset.txt and
    external_weight_dataset.txt are only written when test_columns /
    external_weight_columns are supplied.
    """
    features_file = Path(features_file)
    df = pd.read_csv(features_file)

    if reference_column not in df.columns:
        raise ValueError(f"{reference_column!r} not found in {features_file}")

    test_columns = list(test_columns)
    external_weight_columns = list(external_weight_columns)

    if external_weight_columns and len(external_weight_columns) != 2:
        raise ValueError(
            "external_weight_columns must name exactly two columns, got "
            f"{len(external_weight_columns)}: {external_weight_columns}"
        )

    for label, columns in (
        ("test", test_columns),
        ("external_weight", external_weight_columns),
    ):
        missing = [column for column in columns if column not in df.columns]
        if missing:
            raise ValueError(f"{label} column(s) not found in {features_file}: {missing}")

    assigned = {"feature", reference_column, *test_columns, *external_weight_columns}
    train_columns = [column for column in df.columns if column not in assigned]

    outputs = {
        "train_dataset.txt": df[train_columns].fillna(0),
        "reference_dataset.txt": df[[reference_column]].fillna(0),
    }
    if test_columns:
        outputs["test_dataset.txt"] = df[test_columns].fillna(0)
    if external_weight_columns:
        outputs["external_weight_dataset.txt"] = df[external_weight_columns].fillna(0)

    output_dir = (
        Path(output_dir) if output_dir is not None
        else features_file.resolve().parent.parent / "f_input_data"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    for filename, data in outputs.items():
        data.to_csv(
            output_dir / filename, sep=";", header=False, index=False, encoding="utf-8"
        )

    return outputs


def generate_pseudo_replicates(
    features_file,
    dataset_file,
    ortholog_file,
    retention_levels=(1.0, 0.9, 0.8, 0.7, 0.6, 0.5),
    repeats_per_level=10,
    random_seed=42,
    output_file=None,
    evaluation_output=None,
    missing_value_policy="zero",
):
    """Simulate missingness in one dataset and re-run e_features_collector.py's
    normalization/gene-mapping on each simulated draw.

    dataset_file is a standard aggregated-dataset CSV (Protein IDs/Gene
    names/Fasta headers/QC flag columns + one intensity column per sample),
    the same shape e_features_collector.py consumes. For every
    retention_level, repeats_per_level independent draws are made (only one
    draw for retention_level == 1.0, since there is nothing to randomize):
    within each sample column separately, enough of that column's currently
    non-zero/non-missing measurements are randomly zeroed out so only
    retention_level's fraction of the ORIGINAL available measurements
    remain. QC-flagged rows are dropped first, exactly as in
    e_features_collector.py.

    The log2 + per-column-median centering is then computed on each draw
    AFTER the removal, so the median reflects the reduced set of available
    measurements. Proteins are mapped to genes via ortholog_file (same
    ortholog-matching logic as e_features_collector.py's non-human branch)
    and aggregated to one row per feature by median.

    All resulting columns (one per retention_level x repeat x original
    sample) are concatenated into a single report CSV: rows follow
    features_file's order (one "feature" column), the same shape as
    e_features_collector.py's own output. Column names encode the source
    sample, retention level, and repeat, e.g.
    "M1006_PXD000288_Mouse_Triceps_p90_r03".

    The same values (minus the "feature" column and any labels) are also
    written to f_input_data/evaluation_dataset.txt — semicolon-separated,
    no row or column labels — the format sample_script.py's
    load_semicolon_1d_to_NF expects.

    missing_value_policy selects how zero/missing measurements (and genes
    with no mapped protein at all) are represented — see
    MISSING_VALUE_POLICIES / _missing_value_floor.
    """
    dataset_file = Path(dataset_file)

    with open(features_file, encoding="utf-8") as file:
        features = file.read().splitlines()

    ortholog_map = pd.read_csv(
        ortholog_file, usecols=["Human_gene", "All_unique_orthologs"]
    )
    ortholog_map = ortholog_map[
        ortholog_map["Human_gene"].isin(features)
    ].assign(
        ortholog=lambda frame: frame["All_unique_orthologs"].str.split(";")
    ).explode("ortholog")
    ortholog_map["ortholog"] = ortholog_map["ortholog"].str.strip()

    df = pd.read_csv(dataset_file)
    df = df[~df[EXCLUSION_COLUMNS].eq("+").any(axis=1)].copy()

    prefix = dataset_file.stem
    matching_columns = [column for column in df.columns if column not in NON_SAMPLE_COLUMNS]
    if not matching_columns:
        raise ValueError(f"No sample columns found in {dataset_file}")
    identifier_columns = ["Gene names", "Fasta headers"]

    identifiers = df[identifier_columns].copy()
    raw_measurements = df[matching_columns].apply(pd.to_numeric, errors="coerce")

    listed_genes = identifiers["Gene names"].fillna("").str.split(";")
    fasta_genes = identifiers["Fasta headers"].fillna("").str.findall(r"\bGN=([^\s;]+)")
    all_genes = [
        list(dict.fromkeys(left + right)) for left, right in zip(listed_genes, fasta_genes)
    ]

    rng = np.random.default_rng(random_seed)
    pseudo_columns = []

    for level in retention_levels:
        level_pct = round(level * 100)
        repeats = 1 if level_pct == 100 else repeats_per_level

        for repeat in range(1, repeats + 1):
            masked = raw_measurements.copy()

            if level_pct < 100:
                for column in matching_columns:
                    available_index = masked.index[masked[column].notna() & masked[column].ne(0)]
                    n_remove = len(available_index) - round(len(available_index) * level)
                    if n_remove > 0:
                        remove_index = rng.choice(available_index, size=n_remove, replace=False)
                        masked.loc[remove_index, column] = 0

            # Real, actually measured values (nonzero, non-missing after the
            # simulated removal above) are centered with a log2(x+1)
            # pseudocount against their own column median. Centering is
            # computed from `masked`, so the removals above change the
            # per-column median used here.
            real_mask = masked.ne(0) & masked.notna()
            log2_real = np.log2(masked.where(real_mask) + 1)
            medians = log2_real.median(skipna=True)
            centered_real = log2_real.subtract(medians, axis="columns")

            # Missing-value / zero policy, isolated so it can be revisited on
            # its own later.
            floor_value = _missing_value_floor(missing_value_policy, centered_real, medians)
            centered = centered_real.fillna(floor_value)

            protein_measurements = pd.concat([identifiers, centered], axis=1).assign(
                _row=range(len(df)),
                ortholog=all_genes,
            ).explode("ortholog")
            protein_measurements["ortholog"] = protein_measurements["ortholog"].str.strip()
            protein_measurements = protein_measurements.merge(
                ortholog_map[["Human_gene", "ortholog"]], on="ortholog", how="inner"
            ).drop_duplicates(["_row", "Human_gene"])

            feature_measurements = (
                protein_measurements.groupby("Human_gene")[matching_columns]
                .median()
                .reindex(features)
                # A gene with no mapped protein at all is also assumed zero.
                .fillna(floor_value)
                .reset_index(drop=True)
            )

            suffix = f"p{level_pct}_r{repeat:02d}" if level_pct < 100 else "p100"
            pseudo_columns.append(
                feature_measurements.rename(
                    columns={column: f"{column}_{suffix}" for column in matching_columns}
                )
            )

    combined = pd.concat(pseudo_columns, axis=1)
    combined.insert(0, "feature", features)

    output_file = (
        Path(output_file) if output_file is not None
        else Path(__file__).resolve().parent
        / "f_pseudo_replicates"
        / f"{prefix}_pseudo_replicates.csv"
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, index=False)

    evaluation_output = (
        Path(evaluation_output) if evaluation_output is not None
        else Path(__file__).resolve().parent / "f_input_data" / "evaluation_dataset.txt"
    )
    evaluation_output.parent.mkdir(parents=True, exist_ok=True)
    combined.drop(columns="feature").fillna(0).to_csv(
        evaluation_output, sep=";", header=False, index=False, encoding="utf-8"
    )

    return combined


def build_pathway_indices(
    features_file, genes_file, indices_output=None, breaks_output=None
):
    """Map genes_file's gene symbols to their row index (0-based) in features_file.

    Index order follows genes_file's order. features_file's "feature" values
    are matched after stripping a stray leading BOM character and surrounding
    whitespace, since some upstream feature-list files carried one. Raises if
    a gene from genes_file is not found.

    All mapped genes form a single pathway/subset, so pathbreaks.txt is
    written as the CSR boundaries [0, len(indices)] for that one subset —
    the format SubsetCatalogCSR expects (it requires boundaries[-1] to equal
    len(indices) exactly, since it slices indices[start:end]).
    """
    features_file = Path(features_file)
    genes_file = Path(genes_file)

    features = pd.read_csv(features_file, usecols=["feature"])["feature"]
    normalized_features = features.str.lstrip("﻿").str.strip()
    feature_to_index = {}
    for index, name in enumerate(normalized_features):
        feature_to_index.setdefault(name, index)

    with open(genes_file, encoding="utf-8-sig") as file:
        genes = [line.strip() for line in file if line.strip()]

    missing = [gene for gene in genes if gene not in feature_to_index]
    if missing:
        raise ValueError(f"{len(missing)} gene(s) not found in {features_file}: {missing}")

    indices = [feature_to_index[gene] for gene in genes]
    breaks = [0, len(indices)]

    processed_data_dir = features_file.resolve().parent.parent / "f_input_data"
    indices_output = (
        Path(indices_output) if indices_output is not None
        else processed_data_dir / "pathinds.txt"
    )
    breaks_output = (
        Path(breaks_output) if breaks_output is not None
        else processed_data_dir / "pathbreaks.txt"
    )

    for output_file, values in ((indices_output, indices), (breaks_output, breaks)):
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(
            "\n".join(str(value) for value in values) + "\n", encoding="utf-8"
        )

    return indices, breaks


if __name__ == "__main__":
    features_file = (
        Path(__file__).resolve().parent
        / "e_processed_features"
        / "all_processed_features.csv"
    )

    assign_datasets(
        features_file=features_file,
        reference_column="H1041_PXD010489_100_COXneg",
        test_columns=["H1041_PXD010489_20_COXneg",	"H1041_PXD010489_100_COXpos",	"H1041_PXD010489_20_COXpos"],
        external_weight_columns=["H1041_PXD010489_100_COXneg", "H1040_HPA_skeletal_muscle"],
    )

    build_pathway_indices(
        features_file=features_file,
        genes_file=(
            "/Users/mehman/Projects/PoC_data_processing/selected_features/"
            "mitochondrial_myopathy_intersection_PXD010489.txt"
        ),
    )

    generate_pseudo_replicates(
        features_file=(
            "/Users/mehman/Projects/PoC_data_processing/selected_features/Integrated_features.txt"
        ),
        dataset_file=(
            "/Users/mehman/Projects/PoC_data_processing/scripts/"
            "c&d_aggrigated_datasets/M1006_PXD000288_missingness_effect.csv"
        ),
        ortholog_file=(
            "/Users/mehman/Projects/PoC_data_processing/gene_orthology/"
            "combined_features_mouse_orthologs.csv"
        ),
        missing_value_policy="min_minus_one",
    )

#mehman$ scp -r /Users/mehman/Projects/PoC_data_processing/scripts/f_input_data mehman@c-cb47.ad.cmm.se:~/PoC/