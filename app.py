"""Interactive dashboard for the Loblaw Bio immune cell analysis.

Start with:
    make dashboard
or:
    streamlit run app.py

Reads directly from cell_counts.db, which `make pipeline` creates.
"""

import pandas as pd
import streamlit as st

from teiko import analysis, plots
from teiko.db import DB_PATH

st.set_page_config(page_title="Loblaw Bio: Immune Cell Analysis", layout="wide")


@st.cache_data(show_spinner="Loading frequency table...")
def load_summary() -> pd.DataFrame:
    return analysis.frequencies_with_metadata()


@st.cache_data(show_spinner="Running statistics...")
def load_part3():
    cohort = analysis.response_cohort()
    sample_level = analysis.compare_responders(cohort, unit="sample")
    subject_level = analysis.compare_responders(cohort, unit="subject")
    by_time = analysis.compare_by_timepoint(cohort)
    return cohort, sample_level, subject_level, by_time


@st.cache_data(show_spinner="Querying baseline subset...")
def load_part4():
    return (
        analysis.baseline_samples(),
        analysis.baseline_samples_per_project(),
        analysis.baseline_subjects_by_response(),
        analysis.baseline_subjects_by_sex(),
        analysis.average_b_cells_melanoma_male_responders_baseline(),
    )


STAT_COLUMNS = {
    "population": "Population",
    "n_responders": "n (R)",
    "n_non_responders": "n (NR)",
    "median_responders": "Median % (R)",
    "median_non_responders": "Median % (NR)",
    "median_difference": "Median diff (R - NR)",
    "rank_biserial": "Rank-biserial r",
    "p_value": "p-value",
    "q_value_bh": "BH q-value",
    "significant_raw": "p < 0.05",
    "significant_fdr": "q < 0.05",
}


def format_stats(table: pd.DataFrame) -> pd.DataFrame:
    shown = table[[c for c in STAT_COLUMNS if c in table.columns]].rename(
        columns=STAT_COLUMNS)
    if "time_from_treatment_start" in table.columns:
        shown.insert(0, "Day", table["time_from_treatment_start"])
    return shown


def stats_table(table: pd.DataFrame) -> None:
    st.dataframe(
        format_stats(table),
        hide_index=True,
        width="stretch",
        column_config={
            "Median % (R)": st.column_config.NumberColumn(format="%.2f"),
            "Median % (NR)": st.column_config.NumberColumn(format="%.2f"),
            "Median diff (R - NR)": st.column_config.NumberColumn(format="%.2f"),
            "Rank-biserial r": st.column_config.NumberColumn(format="%.3f"),
            "p-value": st.column_config.NumberColumn(format="%.4f"),
            "BH q-value": st.column_config.NumberColumn(format="%.4f"),
        },
    )


# --------------------------------------------------------------------------
st.title("Loblaw Bio: Immune Cell Population Analysis")
st.caption(f"Data source: `{DB_PATH.name}` (SQLite), built from cell-count.csv")

if not DB_PATH.exists():
    # Normally created by `make pipeline`. Building it here as a fallback lets
    # the dashboard also run on hosts that start from a clean checkout
    # (for example Streamlit Community Cloud).
    import load_data

    with st.spinner("Database not found; building it from cell-count.csv..."):
        if load_data.main() != 0:
            st.error("Could not build the database. Run `make pipeline` and check the output.")
            st.stop()

tab2, tab3, tab4 = st.tabs([
    "Part 2: Data Overview",
    "Part 3: Responder Analysis",
    "Part 4: Baseline Subset",
])

# Part 2 -------------------------------------------------------------------
with tab2:
    st.subheader("Relative frequency of each cell population per sample")
    st.write(
        "For each sample, the total cell count is the sum of the five populations. "
        "Each population's percentage is its count divided by that total."
    )
    summary = load_summary()

    col1, col2, col3, col4 = st.columns(4)
    sample_query = col1.text_input("Sample ID contains", placeholder="e.g. sample0001")
    populations = col2.multiselect("Population", sorted(summary["population"].unique()))
    conditions = col3.multiselect("Condition", sorted(summary["condition"].unique()))
    sample_types = col4.multiselect("Sample type", sorted(summary["sample_type"].unique()))

    filtered = summary
    if sample_query:
        filtered = filtered[filtered["sample"].str.contains(sample_query.strip(), case=False)]
    if populations:
        filtered = filtered[filtered["population"].isin(populations)]
    if conditions:
        filtered = filtered[filtered["condition"].isin(conditions)]
    if sample_types:
        filtered = filtered[filtered["sample_type"].isin(sample_types)]

    table = filtered[["sample", "total_count", "population", "count", "percentage"]]
    m1, m2, m3 = st.columns(3)
    m1.metric("Rows", f"{len(table):,}")
    m2.metric("Samples", f"{filtered['sample'].nunique():,}")
    m3.metric("Subjects", f"{filtered['subject'].nunique():,}")

    st.dataframe(
        table,
        hide_index=True,
        width="stretch",
        height=460,
        column_config={"percentage": st.column_config.NumberColumn(format="%.2f")},
    )
    st.download_button(
        "Download table as CSV",
        table.to_csv(index=False).encode("utf-8"),
        file_name="frequency_summary.csv",
        mime="text/csv",
    )

# Part 3 -------------------------------------------------------------------
with tab3:
    cohort, sample_level, subject_level, by_time = load_part3()
    st.subheader("Responders vs non-responders: melanoma, miraclib, PBMC samples")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Samples", f"{cohort['sample'].nunique():,}")
    c2.metric("Subjects", f"{cohort['subject'].nunique():,}")
    c3.metric("Responder subjects",
              f"{cohort.loc[cohort['response'] == 'yes', 'subject'].nunique():,}")
    c4.metric("Non-responder subjects",
              f"{cohort.loc[cohort['response'] == 'no', 'subject'].nunique():,}")

    time_points = sorted(cohort["time_from_treatment_start"].unique())
    selected_times = st.multiselect(
        "Time points shown in the boxplot (days from treatment start)",
        time_points, default=time_points,
    )
    plot_data = cohort[cohort["time_from_treatment_start"].isin(selected_times)]
    if plot_data.empty:
        st.info("Select at least one time point.")
    else:
        plot_stats = analysis.compare_responders(plot_data, unit="sample")
        st.plotly_chart(plots.boxplot_figure(plot_data, plot_stats),
                        width="stretch")

    st.markdown("#### Findings")
    for line in analysis.interpret(sample_level, subject_level):
        st.markdown(f"- {line}")

    with st.expander("Statistical method", expanded=False):
        st.markdown(
            """
- **Test:** two-sided Mann-Whitney U test per population. It is non-parametric,
  so it does not assume that percentages are normally distributed.
- **Multiple comparisons:** five populations are tested, so p-values are adjusted
  with the Benjamini-Hochberg false discovery rate procedure (q-values).
- **Effect size:** rank-biserial correlation (from -1 to 1). Positive values mean
  responders tend to have a higher relative frequency.
- **Repeated measures:** each subject contributes samples at days 0, 7 and 14,
  which are not independent. The subject-level analysis averages each subject's
  time points before testing, so every person counts once.
- **Significance threshold:** 0.05.
            """
        )

    st.markdown("#### Sample-level results (all time points, primary analysis)")
    stats_table(sample_level)
    st.markdown("#### Subject-level results (time points averaged per subject)")
    stats_table(subject_level)
    st.markdown("#### Results by time point")
    stats_table(by_time)

# Part 4 -------------------------------------------------------------------
with tab4:
    baseline, per_project, by_response, by_sex, b_cells = load_part4()
    st.subheader("Baseline melanoma PBMC samples from miraclib-treated patients")
    st.write("Filter: condition = melanoma, treatment = miraclib, sample type = PBMC, "
             "time from treatment start = 0.")

    st.metric("Baseline samples", f"{len(baseline):,}")

    p1, p2, p3 = st.columns(3)
    with p1:
        st.markdown("**Samples per project**")
        st.dataframe(per_project, hide_index=True, width="stretch")
        st.bar_chart(per_project.set_index("project"))
    with p2:
        st.markdown("**Subjects by response**")
        st.dataframe(by_response, hide_index=True, width="stretch")
        st.bar_chart(by_response.set_index("response"))
    with p3:
        st.markdown("**Subjects by sex**")
        st.dataframe(by_sex, hide_index=True, width="stretch")
        st.bar_chart(by_sex.set_index("sex"))

    st.divider()
    st.markdown("#### Average B cell count: melanoma male responders at time 0")
    st.write("Includes all sample types and all treatments.")
    st.metric("Average B cells", b_cells.formatted,
              help=f"Mean over {b_cells.n_samples} samples")

    with st.expander("Baseline sample list"):
        st.dataframe(baseline, hide_index=True, width="stretch")
