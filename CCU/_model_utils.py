# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 15:47:12 2025

@author: IGB
"""


import pandas as pd


def get_significant_params(rho_df, p_df, 
                           cutoff_p=0.05,
                           cutoff_rho=0.0,
                           round_rho_to=3,
                           round_p_to=3,
                           indicator_filter=None):
                   
    """
    rho_df: DataFrame of Spearman rho (rows = input parameters, cols = outputs/indicators)
    p_df  : DataFrame of p-values (same shape as rho_df)

    indicator_filter: list of substrings to select outputs/indicators by their name
        Example:
            ["Adjusted minimum selling price", "Total gwp100a"]
        Matching is done on the *second element* if the column is a tuple.
    """

    rows = []

    for param in rho_df.index:              # input parameter (e.g. Ash disposal price)
        for out in rho_df.columns:          # output indicator (e.g. Adjusted MSP, Total GWP)

            # Resolve output/indicator name (if tuple, use 2nd element)
            if isinstance(out, tuple):
                out_name = out[1]
            else:
                out_name = out

            # If filtering: keep only selected outputs
            if indicator_filter is not None:
                if not any(key in out_name for key in indicator_filter):
                    continue

            rho = rho_df.loc[param, out]
            p   = p_df.loc[param, out]

            if pd.isna(rho) or pd.isna(p):
                continue

            # significant?
            if (p < cutoff_p) and (abs(rho) >= cutoff_rho):
                rows.append({
                    "indicator": param,   # input parameter (what’s on your Excel “indicator” col)
                    "parameter": out,    # output indicator (Adjusted MSP, Total gwp100a, etc.)
                    "rho": float(rho),
                    "p_value": float(p),
                })

    # Nothing significant
    if not rows:
        return pd.DataFrame(columns=["indicator", "parameter", "rank", "rho", "p_value"])

    sig = pd.DataFrame(rows)

    # rank |rho| within each output indicator
    sig["abs_rho"] = sig["rho"].abs()
    sig = sig.sort_values(["parameter", "abs_rho"], ascending=[True, False])

    # rounding
    sig["rho"] = sig["rho"].round(round_rho_to)
    sig["p_value"] = sig["p_value"].round(round_p_to)

    # rank: 1 = most influential input parameter for this output
    sig["rank"] = (
        sig.groupby("parameter")["abs_rho"]
           .rank(method="dense", ascending=False)
           .astype(int)
    )

    # final column order
    sig = sig[["indicator", "parameter", "rank", "rho", "p_value"]]

    return sig