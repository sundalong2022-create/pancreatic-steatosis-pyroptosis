#!/usr/bin/env python3
"""Marker-based single-cell localization of the human pancreatic pyroptosis axis."""
from __future__ import annotations
import argparse, gzip
from pathlib import Path
import numpy as np
import pandas as pd

CORE = ["TXNIP", "NLRP3", "CASP1", "GSDMD", "IL1B", "IL18"]
EXTENDED = CORE + ["CASP4", "CASP5", "PYCARD", "AIM2", "NLRC4", "NLRP1", "GSDME", "GBP1", "GBP5"]
MARKERS = {
    "macrophage": ["LST1", "TYROBP", "FCER1G", "CTSS", "LILRB1", "AIF1", "C1QC", "CD68", "MS4A7"],
    "acinar": ["PRSS1", "PRSS2", "REG1A", "REG1B", "AMY2A", "AMY2B", "CPA1", "CTRB1", "CTRB2"],
    "ductal": ["KRT19", "KRT8", "KRT18", "KRT7", "MUC1", "KRT17", "KRT23", "CFTR", "SPP1", "KRT16"],
    "stellate": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "SPARC", "DCN", "LUM", "DES"],
    "beta": ["INS", "IAPP", "PCSK1", "ISL1", "MAFA", "DLK1", "ADCYAP1"],
    "alpha": ["GCG", "TTR", "LOXL4"],
    "endothelial": ["KDR", "ESAM", "EMCN", "PECAM1", "VWF", "RAMP2"],
    "t_cell": ["CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "LTB"],
    "epithelial_tumor": ["EPCAM", "KRT8", "KRT18", "KRT19", "MSLN", "CEACAM6", "MUC1", "KRT17"],
}

def opener(path):
    return gzip.open(path, "rt", errors="replace") if str(path).endswith(".gz") else open(path, "r", errors="replace")

def canonical(x):
    return str(x).strip().split(".")[0].upper()

def read_selected_rows(path, wanted, chunksize=512):
    with opener(path) as fh:
        header = fh.readline().rstrip("\n\r").split(",")
        cells, selected = header[1:], {}
        library_size = np.zeros(len(cells), dtype=float)
        for chunk in pd.read_csv(fh, header=None, names=header, chunksize=chunksize, dtype={0: str}, low_memory=False):
            numeric = chunk.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0)
            library_size += numeric.sum(axis=0).to_numpy(dtype=float)
            genes = chunk.iloc[:, 0].map(canonical)
            for i in np.flatnonzero(genes.isin(wanted).to_numpy()):
                gene = genes.iloc[i]
                vals = numeric.iloc[i, :].to_numpy(float)
                selected[gene] = selected.get(gene, 0) + vals
    if not selected:
        raise ValueError("No requested genes found; check species and matrix identifiers")
    return cells, pd.DataFrame(selected, index=cells).T, library_size

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--outdir", default="results/sc_location")
    args = ap.parse_args(); out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    wanted = set(CORE + EXTENDED + [g for x in MARKERS.values() for g in x])
    cells, raw, library_size = read_selected_rows(args.matrix, wanted)
    lib = pd.Series(library_size, index=cells).replace(0, np.nan)
    logcpm = np.log1p(raw.div(lib, axis=1) * 1e4).fillna(0)
    scores = pd.DataFrame(index=cells)
    for name, genes in MARKERS.items():
        present = [g for g in genes if g in logcpm.index]
        scores[name + "_score"] = logcpm.loc[present].mean(axis=0) if present else 0.0
    for name, genes in [("PSPS_core", CORE), ("pyroptosis_extended", EXTENDED)]:
        present = [g for g in genes if g in logcpm.index]
        scores[name] = logcpm.loc[present].mean(axis=0) if present else 0.0
    marker_cols = [x + "_score" for x in MARKERS]
    scores["cell_type"] = scores[marker_cols].idxmax(axis=1).str.replace("_score", "", regex=False)
    scores["top_score"] = scores[marker_cols].max(axis=1)
    scores["second_score"] = scores[marker_cols].apply(lambda x: np.sort(x.to_numpy())[-2], axis=1)
    scores["assignment_margin"] = scores["top_score"] - scores["second_score"]
    scores["low_information"] = scores["top_score"] <= 0
    scores["sample_id"] = scores.index.to_series().str.split(":", n=1).str[0].to_numpy()
    scores.to_csv(out / "per_cell_scores.csv"); raw.to_csv(out / "selected_gene_counts.csv")
    group = scores.groupby("cell_type", observed=True)
    group[["PSPS_core", "pyroptosis_extended", "top_score", "assignment_margin"]].agg(["count", "mean", "median", "std"]).to_csv(out / "cell_type_score_summary.csv")
    (raw > 0).T.assign(cell_type=scores["cell_type"].to_numpy()).groupby("cell_type", observed=True).mean().to_csv(out / "cell_type_gene_detection.csv")
    logcpm.T.assign(cell_type=scores["cell_type"].to_numpy()).groupby("cell_type", observed=True).mean().to_csv(out / "cell_type_gene_mean_logcpm.csv")
    rank = group.agg(n_cells=("PSPS_core", "size"), mean_PSPS=("PSPS_core", "mean"), median_PSPS=("PSPS_core", "median"), mean_extended=("pyroptosis_extended", "mean"), frac_low_information=("low_information", "mean")).sort_values("mean_PSPS", ascending=False)
    rank.to_csv(out / "cell_type_pyroptosis_ranking.csv")
    sample_type = scores.groupby(["sample_id", "cell_type"], observed=True).agg(
        n_cells=("PSPS_core", "size"), mean_PSPS=("PSPS_core", "mean"), median_PSPS=("PSPS_core", "median"),
        mean_extended=("pyroptosis_extended", "mean"), frac_low_information=("low_information", "mean")
    ).reset_index()
    sample_type.to_csv(out / "sample_cell_type_scores.csv", index=False)
    props = scores.groupby(["sample_id", "cell_type"], observed=True).size().rename("n_cells").reset_index()
    props["fraction_within_sample"] = props["n_cells"] / props.groupby("sample_id")["n_cells"].transform("sum")
    props.to_csv(out / "sample_cell_type_proportions.csv", index=False)
    pd.DataFrame({"gene": EXTENDED, "present_in_matrix": [g in raw.index for g in EXTENDED]}).to_csv(out / "gene_presence.csv", index=False)
    (out / "run_manifest.txt").write_text(f"matrix={args.matrix}\nn_cells={len(cells)}\nselected_genes={len(raw)}\nassignment=highest marker-module log1p(CPM) score\ninterpretation=transcriptomic localization only; not proof of functional pyroptotic death\n")
    print(f"cells={len(cells)} selected_genes={len(raw)} output={out}")

if __name__ == "__main__": main()
