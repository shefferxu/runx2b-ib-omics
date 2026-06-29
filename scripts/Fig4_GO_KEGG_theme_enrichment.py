from pathlib import Path
import shutil
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.lines import Line2D

# 设置你的自定义默认颜色循环
plt.rcParams['axes.prop_cycle'] = cycler(
    color=['#A1A9D0', '#F0988C', '#B883D4', '#9E9E9E',
           '#CFEAF1', '#C4A5DE', '#F6CAE5', '#96CCCB']
)

# 设置 viridis 作为默认连续色阶
plt.rcParams['image.cmap'] = 'viridis'

# 中文支持
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False


ROOT_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT_DIR / "output"
SOURCE_DIR = OUTPUT_DIR / "theme_enrichment_reanalysis"
SOURCE_XLSX = SOURCE_DIR / "theme_enrichment_reanalysis.xlsx"

FIG5_DIR = OUTPUT_DIR / "figures" / "Fig5_theme_enrichment"
PANEL_DIR = FIG5_DIR / "single_panels"
DATA_DIR = FIG5_DIR / "source_data"

THEME_ORDER = ["ECM", "Cytoskeleton", "Calcium", "Muscle", "Ossification"]
THEME_COLORS = {
    "ECM": "#A1A9D0",
    "Cytoskeleton": "#F0988C",
    "Calcium": "#B883D4",
    "Muscle": "#96CCCB",
    "Ossification": "#C4A5DE",
}

COMPARE_MARKERS = {
    "He vs Ctrl": "o",
    "Ho vs Ctrl": "^",
    "Ho vs He": "s",
}

COMPARE_Y_OFFSETS = {
    "He vs Ctrl": 0.00,
    "Ho vs Ctrl": 0.13,
    "Ho vs He": -0.13,
}


def read_theme_table() -> pd.DataFrame:
    """读取主题富集结果表。

    参数:
        无。

    返回:
        主题富集结果数据表。
    """
    table = pd.read_excel(SOURCE_XLSX, sheet_name="theme_terms")
    table["padjust"] = pd.to_numeric(table["padjust"], errors="coerce")
    table["pvalue"] = pd.to_numeric(table["pvalue"], errors="coerce")
    table["count"] = pd.to_numeric(table["count"], errors="coerce")
    table["rich_factor"] = pd.to_numeric(table["rich_factor"], errors="coerce")
    table["neg_log10_p"] = pd.to_numeric(table["neg_log10_p"], errors="coerce")
    table = table.replace([np.inf, -np.inf], np.nan)
    table["canonical_themes"] = table["theme"].map(assign_canonical_themes)
    table = table[table["canonical_themes"].map(len) > 0].copy()
    return table


def assign_canonical_themes(theme_value: object) -> list[str]:
    """将原始主题映射到论文主线主题。

    参数:
        theme_value: 原始主题名称。

    返回:
        规范化主题列表。
    """
    if pd.isna(theme_value):
        return []
    text = str(theme_value).lower()
    themes = []
    if "ecm" in text or "adhesion" in text:
        themes.append("ECM")
    if "cytoskeleton" in text:
        themes.append("Cytoskeleton")
    if "calcium" in text:
        themes.append("Calcium")
    if "muscle" in text:
        themes.append("Muscle")
    if "bone" in text or "ossification" in text:
        themes.append("Ossification")
    return themes


def expand_themes(table: pd.DataFrame) -> pd.DataFrame:
    """将复合主题拆分为多行。

    参数:
        table: 主题富集结果表。

    返回:
        拆分后的主题富集表。
    """
    rows = []
    for _, row in table.iterrows():
        for theme in row["canonical_themes"]:
            new_row = row.copy()
            new_row["canonical_theme"] = theme
            rows.append(new_row)
    return pd.DataFrame(rows)


def select_top_terms(table: pd.DataFrame, database: str, max_per_theme: int = 3) -> pd.DataFrame:
    """选择每个主题最代表性的富集条目。

    参数:
        table: 主题富集表。
        database: 富集数据库名称。
        max_per_theme: 每个主题最多保留条目数。

    返回:
        用于气泡图的富集条目。
    """
    expanded = expand_themes(table)
    subset = expanded[expanded["database"] == database].copy()
    if subset.empty:
        return subset
    subset["sort_p"] = subset["padjust"].fillna(subset["pvalue"])
    subset = subset.sort_values(["canonical_theme", "sort_p", "neg_log10_p"], ascending=[True, True, False])
    selected = subset.groupby("canonical_theme", group_keys=False).head(max_per_theme).copy()
    selected["theme_rank"] = selected["canonical_theme"].map({theme: idx for idx, theme in enumerate(THEME_ORDER)})
    selected = selected.sort_values(["theme_rank", "sort_p", "term"])
    return selected


def wrap_term(term: str, width: int = 30) -> str:
    """换行显示富集条目名称。

    参数:
        term: 富集条目名称。
        width: 每行最大字符数。

    返回:
        换行后的条目名称。
    """
    return "\n".join(textwrap.wrap(str(term), width=width))


def build_term_labels(plot_table: pd.DataFrame) -> list[str]:
    """生成富集条目标签。

    参数:
        plot_table: 用于绘图的富集条目表。

    返回:
        list[str]: 换行后的条目标签。
    """
    return [wrap_term(term) for term in plot_table["term"]]


def add_panel_label(ax, label: str) -> None:
    """添加合并图面板标签。

    参数:
        ax: matplotlib 坐标轴。
        label: 面板标签。

    返回:
        无返回值。
    """
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="normal",
    )


def draw_bubble_panel(ax, selected: pd.DataFrame, title: str) -> None:
    """绘制主题富集气泡图。

    参数:
        ax: matplotlib 坐标轴。
        selected: 已筛选的富集条目。
        title: 图标题。

    返回:
        无返回值。
    """
    if selected.empty:
        ax.text(0.5, 0.5, "No matched terms", ha="center", va="center", fontsize=9)
        ax.set_title(title, fontsize=10)
        ax.axis("off")
        return

    plot_table = selected.copy().reset_index(drop=True)
    plot_table["term_label"] = build_term_labels(plot_table)
    unique_terms = list(dict.fromkeys(plot_table["term_label"].tolist()))
    y_lookup = {term: len(unique_terms) - 1 - idx for idx, term in enumerate(unique_terms)}
    plot_table["y_position"] = plot_table["term_label"].map(y_lookup)
    sizes = 20 + plot_table["count"].fillna(1).to_numpy() * 12
    colors = plot_table["neg_log10_p"].fillna(0).to_numpy()

    norm = plt.Normalize(vmin=float(np.nanmin(colors)), vmax=float(np.nanmax(colors)))
    cmap = plt.get_cmap("viridis")
    for compare_label, marker in COMPARE_MARKERS.items():
        subset = plot_table[plot_table["compare_label"] == compare_label]
        if subset.empty:
            continue
        subset_sizes = 20 + subset["count"].fillna(1).to_numpy() * 12
        y_values = subset["y_position"] + COMPARE_Y_OFFSETS.get(compare_label, 0)
        ax.scatter(
            subset["rich_factor"],
            y_values,
            s=subset_sizes,
            c=subset["neg_log10_p"].fillna(0).to_numpy(),
            cmap=cmap,
            norm=norm,
            marker=marker,
            edgecolor="#333333",
            linewidth=0.35,
            alpha=0.9,
        )

    y_positions = [y_lookup[term] for term in unique_terms]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(unique_terms, fontsize=6.5)
    term_theme = plot_table.drop_duplicates("term_label").set_index("term_label")["canonical_theme"].to_dict()
    for tick_label in ax.get_yticklabels():
        theme = term_theme.get(tick_label.get_text(), "")
        tick_label.set_color(THEME_COLORS.get(theme, "black"))
    ax.set_xlabel("Rich factor", fontsize=8)
    ax.set_title(title, fontsize=10, pad=6)
    ax.tick_params(axis="x", labelsize=7)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    marker_handles = [
        Line2D([0], [0], marker=marker, color="none", markerfacecolor="#BDBDBD",
               markeredgecolor="#333333", markersize=5.5, label=compare_label)
        for compare_label, marker in COMPARE_MARKERS.items()
    ]
    legend_loc = "lower right" if title.startswith("KEGG") else "upper right"
    ax.legend(
        handles=marker_handles,
        title="Comparison",
        loc=legend_loc,
        ncol=1,
        frameon=False,
        fontsize=6.3,
        title_fontsize=6.5,
        handletextpad=0.4,
        borderaxespad=0.2,
    )

    scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])
    colorbar = plt.colorbar(scalar_mappable, ax=ax, fraction=0.045, pad=0.02)
    colorbar.set_label("-log10(adjusted P)", fontsize=7)
    colorbar.ax.tick_params(labelsize=6)


def build_theme_summary(table: pd.DataFrame) -> pd.DataFrame:
    """构建五类主题归类统计表。

    参数:
        table: 主题富集结果表。

    返回:
        按数据库、组学和主题统计的结果表。
    """
    expanded = expand_themes(table)
    summary = (
        expanded.groupby(["database", "omics", "canonical_theme"], dropna=False)
        .agg(term_count=("term", "nunique"), mean_neg_log10_p=("neg_log10_p", "mean"))
        .reset_index()
    )
    return summary


def draw_theme_bar_panel(ax, summary: pd.DataFrame) -> None:
    """绘制主题分组条形图。

    参数:
        ax: matplotlib 坐标轴。
        summary: 主题统计表。

    返回:
        无返回值。
    """
    plot_table = (
        summary.groupby("canonical_theme")
        .agg(term_count=("term_count", "sum"), mean_neg_log10_p=("mean_neg_log10_p", "mean"))
        .reindex(THEME_ORDER)
        .reset_index()
    )
    bars = ax.barh(
        plot_table["canonical_theme"],
        plot_table["term_count"],
        color=[THEME_COLORS[theme] for theme in plot_table["canonical_theme"]],
        edgecolor="#333333",
        linewidth=0.45,
    )
    for bar, value in zip(bars, plot_table["term_count"]):
        ax.text(value + max(plot_table["term_count"]) * 0.02, bar.get_y() + bar.get_height() / 2, str(int(value)),
                va="center", fontsize=8)
    ax.set_xlabel("Number of enriched terms", fontsize=8)
    ax.set_title("Theme-level enrichment summary", fontsize=10, pad=6)
    ax.tick_params(labelsize=8)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()


def draw_single_panels(go_terms: pd.DataFrame, kegg_terms: pd.DataFrame, theme_summary: pd.DataFrame) -> None:
    """导出 Fig. 5 单图。

    参数:
        go_terms: GO 富集条目。
        kegg_terms: KEGG 富集条目。
        theme_summary: 主题统计表。

    返回:
        无返回值。
    """
    fig_a, ax_a = plt.subplots(figsize=(5.6, 5.2))
    draw_bubble_panel(ax_a, go_terms, "GO enrichment")
    fig_a.tight_layout()
    fig_a.savefig(PANEL_DIR / "Fig5a_GO_theme_enrichment.png", dpi=300, bbox_inches="tight")
    fig_a.savefig(PANEL_DIR / "Fig5a_GO_theme_enrichment.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig_a)

    fig_b, ax_b = plt.subplots(figsize=(5.6, 5.2))
    draw_bubble_panel(ax_b, kegg_terms, "KEGG enrichment")
    fig_b.tight_layout()
    fig_b.savefig(PANEL_DIR / "Fig5b_KEGG_theme_enrichment.png", dpi=300, bbox_inches="tight")
    fig_b.savefig(PANEL_DIR / "Fig5b_KEGG_theme_enrichment.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig_b)

    fig_c, ax_c = plt.subplots(figsize=(4.6, 3.0))
    draw_theme_bar_panel(ax_c, theme_summary)
    fig_c.tight_layout()
    fig_c.savefig(PANEL_DIR / "Fig5c_theme_grouped_bar.png", dpi=300, bbox_inches="tight")
    fig_c.savefig(PANEL_DIR / "Fig5c_theme_grouped_bar.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig_c)


def draw_combined_figure(go_terms: pd.DataFrame, kegg_terms: pd.DataFrame, theme_summary: pd.DataFrame) -> None:
    """导出 Fig. 5 合并图。

    参数:
        go_terms: GO 富集条目。
        kegg_terms: KEGG 富集条目。
        theme_summary: 主题统计表。

    返回:
        无返回值。
    """
    fig = plt.figure(figsize=(10.5, 6.8))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.15, 1.15], height_ratios=[1.0, 0.65], hspace=0.42, wspace=0.58)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    draw_bubble_panel(ax_a, go_terms, "GO enrichment")
    draw_bubble_panel(ax_b, kegg_terms, "KEGG enrichment")
    draw_theme_bar_panel(ax_c, theme_summary)
    add_panel_label(ax_a, "(a)")
    add_panel_label(ax_b, "(b)")
    add_panel_label(ax_c, "(c)")

    fig.savefig(FIG5_DIR / "Fig5_theme_enrichment.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG5_DIR / "Fig5_theme_enrichment.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_outputs(table: pd.DataFrame, go_terms: pd.DataFrame, kegg_terms: pd.DataFrame, theme_summary: pd.DataFrame) -> None:
    """保存 Fig. 5 源数据和图例说明。

    参数:
        table: 全部主题富集结果。
        go_terms: GO 富集条目。
        kegg_terms: KEGG 富集条目。
        theme_summary: 主题统计表。

    返回:
        无返回值。
    """
    table.to_csv(DATA_DIR / "Fig5_all_theme_enrichment_terms.csv", index=False, encoding="utf-8-sig")
    go_terms.to_csv(DATA_DIR / "Fig5_GO_selected_terms.csv", index=False, encoding="utf-8-sig")
    kegg_terms.to_csv(DATA_DIR / "Fig5_KEGG_selected_terms.csv", index=False, encoding="utf-8-sig")
    theme_summary.to_csv(DATA_DIR / "Fig5_theme_summary.csv", index=False, encoding="utf-8-sig")
    if SOURCE_XLSX.exists():
        shutil.copy2(SOURCE_XLSX, DATA_DIR / SOURCE_XLSX.name)

    legend_text = """# Fig. 5 Theme enrichment analysis

**Figure legend:** (a) GO enrichment bubble plot showing representative terms related to ECM/cell adhesion, cytoskeleton, calcium signaling, muscle contraction/development, and ossification/bone mineralization. (b) KEGG enrichment bubble plot highlighting pathway-level enrichment related to these five biological themes, including ECM-receptor interaction, focal adhesion, calcium signaling, and muscle-related pathways when present. Bubble size represents the number of enriched molecules, and color indicates -log10(adjusted P value). (c) Theme-level summary showing the number of enriched terms assigned to each biological theme.

## Main files

- `Fig5_theme_enrichment.png`
- `Fig5_theme_enrichment.tiff`

## Single panels

- `single_panels/Fig5a_GO_theme_enrichment.png`
- `single_panels/Fig5a_GO_theme_enrichment.tiff`
- `single_panels/Fig5b_KEGG_theme_enrichment.png`
- `single_panels/Fig5b_KEGG_theme_enrichment.tiff`
- `single_panels/Fig5c_theme_grouped_bar.png`
- `single_panels/Fig5c_theme_grouped_bar.tiff`

## Source data

- `source_data/Fig5_all_theme_enrichment_terms.csv`
- `source_data/Fig5_GO_selected_terms.csv`
- `source_data/Fig5_KEGG_selected_terms.csv`
- `source_data/Fig5_theme_summary.csv`
"""
    (FIG5_DIR / "Fig5_legend.md").write_text(legend_text, encoding="utf-8")


def main() -> None:
    """生成并整理 Fig. 5 主题富集图。

    参数:
        无。

    返回:
        无返回值。
    """
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    table = read_theme_table()
    go_terms = select_top_terms(table, "GO", max_per_theme=2)
    kegg_terms = select_top_terms(table, "KEGG", max_per_theme=2)
    theme_summary = build_theme_summary(table)

    draw_single_panels(go_terms, kegg_terms, theme_summary)
    draw_combined_figure(go_terms, kegg_terms, theme_summary)
    write_outputs(table, go_terms, kegg_terms, theme_summary)


if __name__ == "__main__":
    main()
