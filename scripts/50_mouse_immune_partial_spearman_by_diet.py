#!/usr/bin/env python3
"""Exploratory diet-adjusted partial-Spearman sensitivity for mouse PRII modules.

The method ranks PRII and each module, residualizes ranks against the five-level
 diet factor, correlates the residuals, and uses a t approximation with n-k-2 df.
BH adjustment is applied to the nine prespecified module tests.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata, t

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "results/12wk_primary/immune/12wk_immune_sample_scores.csv"
OUTPUT = ROOT / "submission_bmc_genomics/tables/Table_S11_mouse_diet_adjusted_immune_correlations.csv"
MODULES = [
    "leukocyte", "t_cell", "cytotoxic", "b_cell", "myeloid",
    "macrophage_TAM", "neutrophil", "NK", "inflammation_noncore",
]

def bh_adjust(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out

def residualize(values, design):
    return values - design @ np.linalg.lstsq(design, values, rcond=None)[0]

def main():
    dat = pd.read_csv(INPUT)
    levels = sorted(dat["group"].dropna().unique())
    if len(levels) != 5 or dat["PSPS_core"].isna().any():
        raise ValueError("Expected complete PRII data from five diet groups")
    group = pd.Categorical(dat["group"], categories=levels)
    dummies = pd.get_dummies(group, drop_first=True, dtype=float).to_numpy()
    design = np.column_stack([np.ones(len(dat)), dummies])
    covariate_df = dummies.shape[1]
    df = len(dat) - covariate_df - 2
    prii_resid = residualize(rankdata(dat["PSPS_core"].to_numpy()), design)
    rows = []
    for module in MODULES:
        y = dat[module].to_numpy(dtype=float)
        if not np.isfinite(y).all():
            raise ValueError(f"Missing/non-finite values for {module}")
        y_resid = residualize(rankdata(y), design)
        rho = float(np.corrcoef(prii_resid, y_resid)[0, 1])
        stat = rho * np.sqrt(df / (1 - rho * rho))
        pvalue = float(2 * t.sf(abs(stat), df))
        rows.append({
            "module": module,
            "n": len(dat),
            "diet_groups": ";".join(levels),
            "diet_covariate_df": covariate_df,
            "partial_spearman_rho": rho,
            "df": df,
            "p_value": pvalue,
        })
    out = pd.DataFrame(rows)
    out["BH_q"] = bh_adjust(out["p_value"])
    out["method"] = "Spearman partial correlation: rank each variable, residualize on diet-group indicators, correlate residuals; t approximation with n-k-2 df"
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, index=False, float_format="%.8g")
    print(OUTPUT)
    print(out[["module", "partial_spearman_rho", "p_value", "BH_q"]].to_string(index=False))

if __name__ == "__main__":
    main()
