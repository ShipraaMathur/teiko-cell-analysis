"""Boxplots comparing responders and non-responders (Part 3)."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import plotly.express as px  # noqa: E402
import plotly.graph_objects as go  # noqa: E402

GROUP_ORDER = ["Responder", "Non-responder"]
GROUP_COLORS = {"Responder": "#2a7f62", "Non-responder": "#b5543c"}
POPULATION_LABELS = {
    "b_cell": "B cells",
    "cd8_t_cell": "CD8 T cells",
    "cd4_t_cell": "CD4 T cells",
    "nk_cell": "NK cells",
    "monocyte": "Monocytes",
}


def _population_order(cohort: pd.DataFrame) -> list[str]:
    return (cohort[["population", "population_order"]].drop_duplicates()
            .sort_values("population_order")["population"].tolist())


def _significance_label(stats: pd.DataFrame | None, population: str) -> str:
    if stats is None:
        return ""
    row = stats.loc[stats["population"] == population].iloc[0]
    return f"p = {row['p_value']:.3g}, q = {row['q_value_bh']:.3g}"


def boxplot_figure(cohort: pd.DataFrame, stats: pd.DataFrame | None = None) -> go.Figure:
    """Interactive Plotly figure: one panel per population, independent y axes."""
    populations = _population_order(cohort)
    fig = px.box(
        cohort,
        x="response_group",
        y="percentage",
        color="response_group",
        facet_col="population",
        facet_col_spacing=0.05,
        category_orders={"response_group": GROUP_ORDER, "population": populations},
        color_discrete_map=GROUP_COLORS,
        points="outliers",
        hover_data={"sample": True, "subject": True,
                    "time_from_treatment_start": True, "response_group": False},
        labels={"percentage": "Relative frequency (%)", "response_group": ""},
    )
    fig.update_yaxes(matches=None, showticklabels=True)
    fig.update_xaxes(title_text="", showticklabels=False)

    def rename(annotation):
        population = annotation.text.split("=")[-1]
        label = POPULATION_LABELS.get(population, population)
        stat = _significance_label(stats, population)
        annotation.update(text=f"<b>{label}</b><br><sup>{stat}</sup>" if stat
                          else f"<b>{label}</b>")

    fig.for_each_annotation(rename)
    fig.update_layout(
        height=480,
        margin=dict(t=90, l=60, r=20, b=40),
        legend=dict(title="", orientation="h", y=-0.08, x=0.5, xanchor="center"),
        template="plotly_white",
    )
    return fig


def save_static_boxplot(cohort: pd.DataFrame, stats: pd.DataFrame, path) -> None:
    """Static PNG version of the boxplot for the report."""
    populations = _population_order(cohort)
    fig, axes = plt.subplots(1, len(populations), figsize=(4 * len(populations), 4.8))
    for ax, population in zip(axes, populations):
        subset = cohort[cohort["population"] == population]
        data = [subset.loc[subset["response_group"] == g, "percentage"]
                for g in GROUP_ORDER]
        box = ax.boxplot(data, patch_artist=True, widths=0.6,
                         medianprops=dict(color="black"),
                         flierprops=dict(marker="o", markersize=2, alpha=0.4))
        for patch, group in zip(box["boxes"], GROUP_ORDER):
            patch.set_facecolor(GROUP_COLORS[group])
            patch.set_alpha(0.75)
        ax.set_xticks([1, 2], GROUP_ORDER)
        ax.set_title(f"{POPULATION_LABELS.get(population, population)}\n"
                     f"{_significance_label(stats, population)}", fontsize=11)
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("Relative frequency (%)")
    fig.suptitle("Melanoma, miraclib, PBMC: responders vs non-responders "
                 "(two-sided Mann-Whitney U, BH-adjusted q)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
