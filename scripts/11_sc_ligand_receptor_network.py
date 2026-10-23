#!/usr/bin/env python3
"""Sample-level ligand-receptor network for GSE154778 PDAC scRNA-seq.

This is a transparent, marker-free expression-support analysis. OmniPath
ligrecextra interactions are used as a prior; ligand and receptor support are
summarised within each sample and inferred cell type. Cells are not treated as
independent replicates.
"""
from __future__ import annotations
import argparse
import gzip
import io
import os
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from scipy.stats import spearmanr

TARGET_TYPES = ["macrophage", "epithelial_tumor", "stellate", "t_cell"]
OMNI_URL = "https://omnipathdb.org/interactions?datasets=ligrecextra&format=tsv&genesymbols=yes"
# A biologically focused subset keeps the network interpretable and centers
# inflammation, chemokines, myeloid recruitment, TGF-beta, growth-factor and
# Notch axes that can plausibly connect macrophages to the three receiver
# compartments in this PDAC reference dataset.
FOCUSED_LIGANDS = set("IL1A IL1B IL18 TNF CCL2 CCL3 CCL4 CCL5 CXCL1 CXCL2 CXCL3 CXCL5 CXCL8 CXCL9 CXCL10 CXCL11 CXCL12 CXCL13 TGFB1 TGFB2 SPP1 CSF1 CSF2 CSF3 IFNG IL6 IL10 IL15 IL2 IL7 IL12A IL12B IL23A VEGFA PDGFA PDGFB EGF FGF2 HGF GAS6 JAG1 JAG2 DLL1 DLL4 ICAM1 VCAM1 CD40LG LTB".split())
FOCUSED_RECEPTORS = set("IL1R1 IL1RAP IL18R1 TNFRSF1A TNFRSF1B CCR2 CCR5 CXCR1 CXCR2 CXCR3 CXCR4 CXCR5 TGFBR1 TGFBR2 CD44 CSF1R CSF2RA CSF3R IFNGR1 IFNGR2 IL6R IL10RA IL15RA IL2RA IL7R IL12RB1 IL12RB2 IL23R KDR FLT1 PDGFRA PDGFRB EGFR FGFR1 MET AXL JAG1 JAG2 NOTCH1 NOTCH2 NOTCH3 NOTCH4 ICAM1 ITGA4 ITGB1 CD40 LTBR".split())

def canonical(x: str) -> str:
    return str(x).strip().split(".")[0].upper()

def read_header(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        header = fh.readline().rstrip("\r\n").split(",")
    return header[1:]

def load_lr(path: Path) -> pd.DataFrame:
    if path.exists():
        x = pd.read_csv(path, sep="\t")
    else:
        r = requests.get(OMNI_URL, timeout=60)
        r.raise_for_status()
        path.write_bytes(r.content)
        x = pd.read_csv(io.BytesIO(r.content), sep="\t")
    x = x.rename(columns={"source_genesymbol": "ligand", "target_genesymbol": "receptor"})
    keep = ["ligand", "receptor"] + [c for c in ["is_directed", "is_stimulation", "is_inhibition"] if c in x.columns]
    x = x[keep].copy()
    x["ligand"] = x["ligand"].map(canonical); x["receptor"] = x["receptor"].map(canonical)
    x = x.dropna(subset=["ligand", "receptor"]).drop_duplicates(["ligand", "receptor"])
    return x

def read_selected_matrix(path: Path, wanted: set[str], chunksize: int = 256):
    # Parsing one sparse CSV line at a time avoids materialising the full
    # 30k x 15k matrix in pandas. np.fromstring performs the numeric parse in
    # compiled code and is substantially faster and more memory-stable here.
    cells = read_header(path); n = len(cells); lib = np.zeros(n, dtype=float); selected = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        _ = fh.readline()
        for line in fh:
            line = line.rstrip("\r\n")
            if not line: continue
            gene_raw, sep, rest = line.partition(",")
            if not sep: continue
            vals = np.fromstring(rest, sep=",", dtype=float)
            if vals.size != n:
                if vals.size < n: vals = np.pad(vals, (0, n - vals.size))
                else: vals = vals[:n]
            lib += vals
            g = canonical(gene_raw)
            if g in wanted: selected[g] = selected.get(g, np.zeros(n, dtype=float)) + vals
    raw = pd.DataFrame(selected, index=cells).T
    logcpm = np.log1p(raw.div(pd.Series(lib, index=cells).replace(0, np.nan), axis=1) * 1e4).fillna(0)
    return logcpm

def sample_type_means(logcpm: pd.DataFrame, cell_scores: pd.DataFrame):
    cs = cell_scores.copy(); cs.index = cs.index.astype(str)
    common = logcpm.columns.intersection(cs.index)
    logcpm = logcpm.loc[:, common]; cs = cs.loc[common]
    cs["sample_id"] = cs["sample_id"].astype(str)
    cs["cell_type"] = cs["cell_type"].astype(str)
    rows = []
    for (sample, ct), ix in cs.groupby(["sample_id", "cell_type"], observed=True).groups.items():
        if ct not in TARGET_TYPES: continue
        sub = logcpm.loc[:, ix]
        m = sub.mean(axis=1); d = (sub > 0).mean(axis=1)
        for g in logcpm.index:
            rows.append((sample, ct, g, float(m[g]), float(d[g]), int(len(ix))))
    return pd.DataFrame(rows, columns=["sample_id", "cell_type", "gene", "mean_logcpm", "detection_fraction", "n_cells"])

def safe_rho(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4 or np.nanstd(x[ok]) == 0 or np.nanstd(y[ok]) == 0: return np.nan, np.nan
    r, p = spearmanr(x[ok], y[ok]); return float(r), float(p)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default="data/raw/GSE154778/GSE154778_dgeMtx.csv.gz")
    ap.add_argument("--cell-scores", default="results/sc_location_GSE154778/per_cell_scores.csv")
    ap.add_argument("--outdir", default="results/network_ligand_receptor")
    args = ap.parse_args(); out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    lr_all = load_lr(out / "omnipath_ligrecextra.tsv")
    # Keep database provenance explicit in every prior-edge table.
    lr_all["evidence_type"] = "数据库先验"
    lr_all.to_csv(out / "omnipath_ligrecextra_all_edges.csv", index=False)
    lr = lr_all[lr_all.ligand.isin(FOCUSED_LIGANDS) & lr_all.receptor.isin(FOCUSED_RECEPTORS)].copy()
    cell_scores = pd.read_csv(args.cell_scores, index_col=0)
    wanted = set(lr["ligand"]) | set(lr["receptor"])
    logcpm = read_selected_matrix(Path(args.matrix), wanted)
    lr = lr[lr.ligand.isin(logcpm.index) & lr.receptor.isin(logcpm.index)].copy()
    lr["evidence_type"] = "数据库先验"
    lr.to_csv(out / "ligand_receptor_prior_edges.csv", index=False)
    st = sample_type_means(logcpm, cell_scores); st.to_csv(out / "sample_celltype_gene_support.csv", index=False)
    # Wide lookup for each sample and cell type.
    key = st.set_index(["sample_id", "cell_type", "gene"])
    samples = sorted(st.sample_id.unique())
    edge_rows = []
    for _, e in lr.iterrows():
        for target in TARGET_TYPES[1:]:
            vals = []
            for s in samples:
                try:
                    ls = key.loc[(s, TARGET_TYPES[0], e.ligand)]
                    rr = key.loc[(s, target, e.receptor)]
                except KeyError:
                    continue
                score = float(np.sqrt(max(ls.mean_logcpm, 0) * max(rr.mean_logcpm, 0)))
                supported = bool(ls.detection_fraction >= 0.10 and rr.detection_fraction >= 0.10)
                vals.append((s, score, supported, ls.mean_logcpm, rr.mean_logcpm, ls.detection_fraction, rr.detection_fraction, ls.n_cells, rr.n_cells))
            for v in vals:
                edge_rows.append(dict(ligand=e.ligand, receptor=e.receptor, sender=TARGET_TYPES[0], receiver=target, sample_id=v[0], interaction_score=v[1], supported=v[2], ligand_mean_logcpm=v[3], receptor_mean_logcpm=v[4], ligand_detection=v[5], receptor_detection=v[6], sender_n_cells=v[7], receiver_n_cells=v[8], evidence_type="候选配体–受体表达支持"))
    detail = pd.DataFrame(edge_rows); detail.to_csv(out / "sample_level_lr_scores.csv", index=False)
    # Sample-level macrophage PSPS is used only for exploratory association.
    psps = cell_scores.loc[cell_scores.cell_type.eq("macrophage")].groupby("sample_id", observed=True)["PSPS_core"].mean()
    summary = []
    for (lig, rec, sender, target), g in detail.groupby(["ligand", "receptor", "sender", "receiver"], observed=True):
        # Correlations are descriptive only and are suppressed when the edge is
        # unsupported in fewer than three samples; otherwise zero-inflated
        # detection indicators can create unstable ranks.
        if int(g.supported.sum()) < 3:
            r, p = np.nan, np.nan
        else:
            r, p = safe_rho(g.sample_id.map(psps).to_numpy(float), g.interaction_score.to_numpy(float))
        summary.append(dict(ligand=lig, receptor=rec, sender=sender, receiver=target, n_samples=int(g.shape[0]), n_supported_samples=int(g.supported.sum()), prevalence=float(g.supported.mean()), mean_score=float(g.interaction_score.mean()), median_score=float(g.interaction_score.median()), spearman_rho_psps=r, spearman_p_psps=p, evidence_type="候选配体–受体表达支持"))
    summ = pd.DataFrame(summary)
    if len(summ):
        from statsmodels.stats.multitest import multipletests
        valid = summ.spearman_p_psps.notna(); summ["spearman_fdr_psps"] = np.nan
        summ.loc[valid, "spearman_fdr_psps"] = multipletests(summ.loc[valid, "spearman_p_psps"], method="fdr_bh")[1]
        summ = summ.sort_values(["prevalence", "mean_score"], ascending=False)
    summ.to_csv(out / "lr_edge_summary.csv", index=False)
    # Keep a compact top-edge table for plotting.
    top = summ[(summ.prevalence > 0) & (summ.n_supported_samples >= 3)].head(40) if len(summ) else summ
    top.to_csv(out / "lr_top_edges.csv", index=False)
    with open(out / "run_manifest.txt", "w") as f:
        f.write(f"matrix={args.matrix}\nsource=OmniPath ligrecextra ({OMNI_URL})\n")
        f.write(f"n_cells={len(cell_scores)}\nn_lr_prior_all={len(lr_all)}\nn_lr_prior_focused={len(lr)}\nn_lr_receiver_pairs={len(summ)}\n")
        f.write("unit_of_replication=sample; expression_support=mean log1p(CPM) and detection fraction; score=geometric mean of sender ligand and receiver receptor expression\n")
        f.write("evidence_labels=prior tables: 数据库先验; sample-level and summary edges: 候选配体–受体表达支持\n")
        f.write("interpretation=exploratory candidate communication network; not proof of ligand secretion, receptor activation, or causality\n")
    print(f"cells={len(cell_scores)} prior_edges={len(lr)} receiver_pairs={len(summ)} top_edges={len(top)} out={out}")

if __name__ == "__main__": main()
