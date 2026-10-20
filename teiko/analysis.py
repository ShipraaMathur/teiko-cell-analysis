"""Analytical functions for Parts 2 to 4.

Every function reads from the SQLite database created by load_data.py, so the
pipeline and the dashboard always report the same numbers.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from .db import query

ALPHA = 0.05

# Reusable SQL: one row per (sample, population) with the sample total and
# the relative frequency, joined to all sample and subject metadata.
_FREQUENCY_CTE = """
WITH totals AS (
    SELECT sample_id, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
),
frequencies AS (
    SELECT
        s.sample_id                       AS sample,
        t.total_count                     AS total_count,
        p.name                            AS population,
        p.population_id                   AS population_order,
        c.count                           AS count,
        100.0 * c.count / t.total_count   AS percentage,
        s.subject_id                      AS subject,
        s.sample_type                     AS sample_type,
        s.time_from_treatment_start       AS time_from_treatment_start,
        sub.project_id                    AS project,
        sub.condition                     AS condition,
        sub.age                           AS age,
        sub.sex                           AS sex,
        sub.treatment                     AS treatment,
        sub.response                      AS response
    FROM cell_counts c
    JOIN totals      t   ON t.sample_id = c.sample_id
    JOIN populations p   ON p.population_id = c.population_id
    JOIN samples     s   ON s.sample_id = c.sample_id
    JOIN subjects    sub ON sub.subject_id = s.subject_id
)
"""


# --------------------------------------------------------------------------
# Part 2: relative frequency summary table
# --------------------------------------------------------------------------
def frequency_summary() -> pd.DataFrame:
    """Return one row per sample and population with relative frequencies."""
    sql = _FREQUENCY_CTE + """
    SELECT sample, total_count, population, count, percentage
    FROM frequencies
    ORDER BY sample, population_order
    """
    return query(sql)


def frequencies_with_metadata(
    condition: str | None = None,
    treatment: str | None = None,
    sample_type: str | None = None,
) -> pd.DataFrame:
    """Return the frequency table joined to metadata, optionally filtered."""
    clauses, params = [], {}
    for column, value in (("condition", condition), ("treatment", treatment),
                          ("sample_type", sample_type)):
        if value is not None:
            clauses.append(f"{column} = :{column}")
            params[column] = value
    where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
    sql = _FREQUENCY_CTE + f"""
    SELECT * FROM frequencies
    {where}
    ORDER BY sample, population_order
    """
    return query(sql, params)


# --------------------------------------------------------------------------
# Part 3: responders versus non-responders
# --------------------------------------------------------------------------
def response_cohort() -> pd.DataFrame:
    """Melanoma PBMC samples from miraclib-treated subjects with known response."""
    df = frequencies_with_metadata(
        condition="melanoma", treatment="miraclib", sample_type="PBMC"
    )
    df = df[df["response"].isin(["yes", "no"])].copy()
    df["response_group"] = df["response"].map(
        {"yes": "Responder", "no": "Non-responder"}
    )
    return df


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (false discovery rate)."""
    p = np.asarray(p_values, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted


def _mann_whitney(responders: pd.Series, non_responders: pd.Series) -> dict:
    result = mannwhitneyu(responders, non_responders, alternative="two-sided")
    n1, n2 = len(responders), len(non_responders)
    # Rank-biserial correlation: positive values mean responders tend higher.
    rank_biserial = 2 * result.statistic / (n1 * n2) - 1
    return {"u_statistic": result.statistic, "p_value": result.pvalue,
            "rank_biserial": rank_biserial}


def compare_responders(cohort: pd.DataFrame, unit: str = "sample") -> pd.DataFrame:
    """Two-sided Mann-Whitney U test per population with BH correction.

    unit="sample": every sample is one observation (all time points).
    unit="subject": each subject's percentages are averaged across time
        points first, so repeated samples from one person count once. This
        guards against inflated significance from non-independent samples.
    """
    data = cohort
    if unit == "subject":
        data = (cohort.groupby(["subject", "response_group", "population",
                                "population_order"], as_index=False)
                ["percentage"].mean())
    elif unit != "sample":
        raise ValueError("unit must be 'sample' or 'subject'")

    rows = []
    for (population, _), group in data.groupby(["population", "population_order"],
                                               sort=True):
        resp = group.loc[group["response_group"] == "Responder", "percentage"]
        non = group.loc[group["response_group"] == "Non-responder", "percentage"]
        rows.append({
            "population": population,
            "population_order": _,
            "n_responders": len(resp),
            "n_non_responders": len(non),
            "median_responders": resp.median(),
            "median_non_responders": non.median(),
            "median_difference": resp.median() - non.median(),
            "mean_responders": resp.mean(),
            "mean_non_responders": non.mean(),
            **_mann_whitney(resp, non),
        })

    table = pd.DataFrame(rows).sort_values("population_order")
    table["q_value_bh"] = benjamini_hochberg(table["p_value"].to_numpy())
    table["significant_raw"] = table["p_value"] < ALPHA
    table["significant_fdr"] = table["q_value_bh"] < ALPHA
    return table.drop(columns="population_order").reset_index(drop=True)


def compare_by_timepoint(cohort: pd.DataFrame) -> pd.DataFrame:
    """Per time point Mann-Whitney U tests (BH correction within each time point)."""
    frames = []
    for time_point, group in cohort.groupby("time_from_treatment_start"):
        table = compare_responders(group, unit="sample")
        table.insert(0, "time_from_treatment_start", time_point)
        frames.append(table)
    return pd.concat(frames, ignore_index=True)


def interpret(sample_level: pd.DataFrame, subject_level: pd.DataFrame) -> list[str]:
    """Plain-language conclusions for each population."""
    merged = sample_level.merge(subject_level, on="population",
                                suffixes=("", "_subject"))
    lines = []
    for _, row in merged.iterrows():
        direction = "higher" if row["median_difference"] > 0 else "lower"
        base = (f"{row['population']}: median {row['median_responders']:.2f}% in "
                f"responders vs {row['median_non_responders']:.2f}% in "
                f"non-responders (p = {row['p_value']:.4f}, "
                f"BH q = {row['q_value_bh']:.4f}, "
                f"subject-level p = {row['p_value_subject']:.4f}).")
        if row["significant_fdr"]:
            verdict = (f" Significant after FDR correction; responders show "
                       f"{direction} relative frequency.")
        elif row["significant_raw"] and row["significant_raw_subject"]:
            verdict = (f" Nominally significant (p < {ALPHA}) at both the sample "
                       f"and subject level, with responders {direction}, but not "
                       f"significant after correcting for five comparisons.")
        elif row["significant_raw"]:
            verdict = (f" Nominally significant at the sample level only; it "
                       f"does not hold at the subject level or after correction.")
        else:
            verdict = " No significant difference."
        lines.append(base + verdict)
    return lines


# --------------------------------------------------------------------------
# Part 4: baseline subset analysis
# --------------------------------------------------------------------------
_BASELINE_FILTER = """
    sub.condition = 'melanoma'
    AND sub.treatment = 'miraclib'
    AND s.sample_type = 'PBMC'
    AND s.time_from_treatment_start = 0
"""


def baseline_samples() -> pd.DataFrame:
    """Melanoma PBMC baseline samples from miraclib-treated subjects."""
    sql = f"""
    SELECT s.sample_id AS sample, s.subject_id AS subject, sub.project_id AS project,
           sub.condition, sub.treatment, sub.response, sub.sex, sub.age,
           s.sample_type, s.time_from_treatment_start
    FROM samples s
    JOIN subjects sub ON sub.subject_id = s.subject_id
    WHERE {_BASELINE_FILTER}
    ORDER BY s.sample_id
    """
    return query(sql)


def baseline_samples_per_project() -> pd.DataFrame:
    sql = f"""
    SELECT sub.project_id AS project, COUNT(*) AS n_samples
    FROM samples s
    JOIN subjects sub ON sub.subject_id = s.subject_id
    WHERE {_BASELINE_FILTER}
    GROUP BY sub.project_id
    ORDER BY sub.project_id
    """
    return query(sql)


def baseline_subjects_by_response() -> pd.DataFrame:
    sql = f"""
    SELECT CASE sub.response WHEN 'yes' THEN 'responder'
                             WHEN 'no'  THEN 'non-responder'
                             ELSE 'unknown' END AS response,
           COUNT(DISTINCT sub.subject_id) AS n_subjects
    FROM samples s
    JOIN subjects sub ON sub.subject_id = s.subject_id
    WHERE {_BASELINE_FILTER}
    GROUP BY sub.response
    ORDER BY sub.response DESC
    """
    return query(sql)


def baseline_subjects_by_sex() -> pd.DataFrame:
    sql = f"""
    SELECT CASE sub.sex WHEN 'M' THEN 'male' WHEN 'F' THEN 'female' END AS sex,
           COUNT(DISTINCT sub.subject_id) AS n_subjects
    FROM samples s
    JOIN subjects sub ON sub.subject_id = s.subject_id
    WHERE {_BASELINE_FILTER}
    GROUP BY sub.sex
    ORDER BY sub.sex DESC
    """
    return query(sql)


@dataclass
class BCellAnswer:
    average_b_cells: float
    n_samples: int

    @property
    def formatted(self) -> str:
        return f"{self.average_b_cells:.2f}"


def average_b_cells_melanoma_male_responders_baseline() -> BCellAnswer:
    """Mean B cell count for melanoma male responders at time 0.

    Includes all sample types and all treatments, as specified.
    """
    sql = """
    SELECT AVG(c.count) AS average_b_cells, COUNT(*) AS n_samples
    FROM cell_counts c
    JOIN populations p   ON p.population_id = c.population_id
    JOIN samples     s   ON s.sample_id = c.sample_id
    JOIN subjects    sub ON sub.subject_id = s.subject_id
    WHERE p.name = 'b_cell'
      AND sub.condition = 'melanoma'
      AND sub.sex = 'M'
      AND sub.response = 'yes'
      AND s.time_from_treatment_start = 0
    """
    row = query(sql).iloc[0]
    return BCellAnswer(float(row["average_b_cells"]), int(row["n_samples"]))
