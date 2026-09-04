from pathlib import Path
import pandas as pd

FEATURES_FILE = Path("/Users/mehman/Projects/PoC_data_processing/selected_features/mitochondrial_myopathy_features_intersect.txt")
ORTHOLOGS_FILE = Path("/Users/mehman/Projects/PoC_data_processing/orthology/combined_orthologs_mitochondrial_myopathy.csv")
OUTPUT_FILE = Path("/Users/mehman/Projects/PoC_data_processing/mitochondrial_myopathy_zebrafish_features_orthologs.txt")

human_genes = {
    line.strip().upper()
    for line in FEATURES_FILE.read_text(encoding="utf-8-sig").splitlines()
    if line.strip()
}

table = pd.read_csv(ORTHOLOGS_FILE)

rows = table[
    table["Human_gene"].str.strip().str.upper().isin(human_genes)
]

zebrafish_orthologs = sorted({
    ortholog.strip()
    for cell in rows["All_unique_orthologs"].dropna()
    for ortholog in cell.split(";")
    if ortholog.strip()
})

OUTPUT_FILE.write_text(
    "\n".join(zebrafish_orthologs) + "\n",
    encoding="utf-8",
)

print(f"Saved {len(zebrafish_orthologs)} unique zebrafish orthologs")