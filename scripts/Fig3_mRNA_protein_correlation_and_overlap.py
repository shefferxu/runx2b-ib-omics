from pathlib import Path
import math
import shutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.patches import Circle

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
SOURCE_DIR = OUTPUT_DIR / "correlation_and_venn"
MATCHED_TABLE = OUTPUT_DIR / "nine_quadrant_plots" / "nine_quadrant_matched_gene_protein.csv"
SUMMARY_FILE = SOURCE_DIR / "correlation_and_venn_summary.xlsx"

FIG3_DIR = OUTPUT_DIR / "figures" / "Fig3_mRNA_protein_correlation_overlap"
PANEL_DIR = FIG3_DIR / "single_panels"
DATA_DIR = FIG3_DIR / "source_data"

COMPARE_ORDER = ["He|Ctrl", "Ho|Ctrl", "Ho|He"]
COMPARE_LABELS = {
    "He|Ctrl": "He vs Ctrl",
    "Ho|Ctrl": "Ho vs Ctrl",
    "Ho|He": "Ho vs He",
}


def read_matched_table() -> pd.DataFrame:
    """读取 mRNA-protein 匹配表。

    参数:
        无。

    返回:
        清洗后的匹配数据表。
    """
    table = pd.read_csv(MATCHED_TABLE)
    table["gene_log2fc"] = pd.to_numeric(table["gene_log2fc"], errors="coerce")
    table["protein_log2fc"] = pd.to_numeric(table["protein_log2fc"], errors="coerce")
    table["gene_significant"] = table["gene_significant"].astype(str).str.lower()
    table["protein_significant"] = table["protein_significant"].astype(str).str.lower()
    table = table.replace([np.inf, -np.inf], np.nan)
    return table.dropna(subset=["gene_log2fc", "protein_log2fc"])


def read_summary_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """读取相关性、Venn 和 overlap 汇总表。

    参数:
        无。

    返回:
        correlation、venn_counts 和 DEG_DEP_overlap 三个数据表。
    """
    correlation = pd.read_excel(SUMMARY_FILE, sheet_name="correlation")
    venn_counts = pd.read_excel(SUMMARY_FILE, sheet_name="venn_counts")
    overlap = pd.read_excel(SUMMARY_FILE, sheet_name="DEG_DEP_overlap")
    return correlation, venn_counts, overlap


def clip_limits(values: pd.Series, min_abs: float = 2.0) -> tuple[float, float]:
    """计算稳健坐标范围。

    参数:
        values: 数值序列。
        min_abs: 最小绝对范围。

    返回:
        坐标轴下限和上限。
    """
    clean_values = values.dropna()
    lower = float(clean_values.quantile(0.005))
    upper = float(clean_values.quantile(0.995))
    limit = max(abs(lower), abs(upper), min_abs)
    return -limit * 1.08, limit * 1.08


def format_p_value(value: float) -> str:
    """格式化 P 值。

    参数:
        value: P 值。

    返回:
        格式化后的字符串。
    """
    if pd.isna(value):
        return "NA"
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}"


def add_panel_label(ax, label: str) -> None:
    """添加子图面板标注。

    参数:
        ax: matplotlib 坐标轴。
        label: 面板标注文本。

    返回:
        无返回值。
    """
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="normal",
    )


def draw_correlation_panel(ax, subset: pd.DataFrame, stats_row: pd.Series, compare: str) -> None:
    """绘制 mRNA-protein log2FC 相关性散点图。

    参数:
        ax: matplotlib 坐标轴。
        subset: 当前比较的数据。
        stats_row: 当前比较的相关性统计。
        compare: 比较名称。

    返回:
        无返回值。
    """
    both = (subset["gene_significant"] == "yes") & (subset["protein_significant"] == "yes")
    one = ((subset["gene_significant"] == "yes") | (subset["protein_significant"] == "yes")) & ~both
    neither = ~(both | one)

    ax.scatter(
        subset.loc[neither, "gene_log2fc"],
        subset.loc[neither, "protein_log2fc"],
        s=5,
        color="#C7C7C7",
        alpha=0.28,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        subset.loc[one, "gene_log2fc"],
        subset.loc[one, "protein_log2fc"],
        s=9,
        color="#A1A9D0",
        alpha=0.62,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        subset.loc[both, "gene_log2fc"],
        subset.loc[both, "protein_log2fc"],
        s=28,
        color="#F0988C",
        edgecolor="#333333",
        linewidths=0.35,
        alpha=0.92,
        rasterized=True,
    )

    fit_table = subset[["gene_log2fc", "protein_log2fc"]].dropna()
    if len(fit_table) >= 2:
        slope, intercept = np.polyfit(fit_table["gene_log2fc"], fit_table["protein_log2fc"], 1)
        x_values = np.linspace(fit_table["gene_log2fc"].quantile(0.01), fit_table["gene_log2fc"].quantile(0.99), 100)
        ax.plot(x_values, slope * x_values + intercept, color="#333333", linewidth=0.9)

    ax.axvline(0, color="#777777", linestyle="--", linewidth=0.65)
    ax.axhline(0, color="#777777", linestyle="--", linewidth=0.65)
    ax.set_title(COMPARE_LABELS[compare], fontsize=9, pad=4)
    ax.set_xlabel("mRNA log2FC", fontsize=8)
    ax.set_ylabel("Protein log2FC", fontsize=8)
    ax.set_xlim(*clip_limits(subset["gene_log2fc"]))
    ax.set_ylim(*clip_limits(subset["protein_log2fc"]))
    ax.grid(True, color="#E6E6E6", linewidth=0.4)
    ax.tick_params(labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    text = (
        f"Pearson r={stats_row['pearson_r']:.3f}\n"
        f"Spearman r={stats_row['spearman_r']:.3f}\n"
        f"P={format_p_value(stats_row['pearson_p'])}\n"
        f"n={int(stats_row['n'])}"
    )
    ax.text(0.03, 0.97, text, transform=ax.transAxes, ha="left", va="top", fontsize=6.8)


def draw_venn_panel(ax, row: pd.Series) -> None:
    """绘制 DEG/DEP overlap 韦恩图。

    参数:
        ax: matplotlib 坐标轴。
        row: Venn 数量统计行。

    返回:
        无返回值。
    """
    left_circle = Circle((0.41, 0.54), 0.29, facecolor="#A1A9D0", edgecolor="#6E78A8", alpha=0.58, linewidth=1.0)
    right_circle = Circle((0.59, 0.54), 0.29, facecolor="#F0988C", edgecolor="#C96D65", alpha=0.58, linewidth=1.0)
    ax.add_patch(left_circle)
    ax.add_patch(right_circle)
    ax.text(0.27, 0.54, str(int(row["deg_only"])), ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0.50, 0.54, str(int(row["overlap"])), ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0.73, 0.54, str(int(row["dep_only"])), ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0.30, 0.18, "DEG only", ha="center", va="center", fontsize=7)
    ax.text(0.70, 0.18, "DEP only", ha="center", va="center", fontsize=7)
    ax.text(0.50, 0.87, "Overlap", ha="center", va="center", fontsize=7)
    ax.set_title(row["compare"], fontsize=9, pad=4)
    ax.set_xlim(0.08, 0.92)
    ax.set_ylim(0.08, 0.96)
    ax.set_aspect("equal")
    ax.axis("off")


def build_direction_counts(overlap: pd.DataFrame) -> pd.DataFrame:
    """统计 overlap 基因/蛋白上下调方向。

    参数:
        overlap: DEG-DEP overlap 明细表。

    返回:
        每个比较的方向一致性统计表。
    """
    rows = []
    for compare in ["He vs Ctrl", "Ho vs Ctrl", "Ho vs He"]:
        sub = overlap[overlap["compare_label"] == compare].copy()
        up_up = int(((sub["gene_log2fc"] > 0) & (sub["protein_log2fc"] > 0)).sum())
        down_down = int(((sub["gene_log2fc"] < 0) & (sub["protein_log2fc"] < 0)).sum())
        opposite = int(((sub["gene_log2fc"] * sub["protein_log2fc"]) < 0).sum())
        rows.append({"compare": compare, "Up-Up": up_up, "Down-Down": down_down, "Opposite": opposite})
    return pd.DataFrame(rows)


def draw_direction_panel(ax, direction_counts: pd.DataFrame) -> None:
    """绘制 DEG-DEP overlap 方向一致性统计图。

    参数:
        ax: matplotlib 坐标轴。
        direction_counts: 方向一致性统计表。

    返回:
        无返回值。
    """
    x_positions = np.arange(len(direction_counts))
    bottom = np.zeros(len(direction_counts))
    colors = {"Up-Up": "#F0988C", "Down-Down": "#A1A9D0", "Opposite": "#9E9E9E"}

    for category in ["Up-Up", "Down-Down", "Opposite"]:
        values = direction_counts[category].to_numpy()
        ax.bar(
            x_positions,
            values,
            bottom=bottom,
            color=colors[category],
            edgecolor="#333333",
            linewidth=0.45,
            width=0.62,
            label=category,
        )
        for x_value, y_value, base in zip(x_positions, values, bottom):
            if y_value > 0:
                ax.text(x_value, base + y_value / 2, str(int(y_value)), ha="center", va="center", fontsize=8)
        bottom += values

    ax.set_xticks(x_positions)
    ax.set_xticklabels(direction_counts["compare"], fontsize=8)
    ax.set_ylabel("Number of DEG-DEP overlaps", fontsize=8)
    ax.set_ylim(0, max(bottom.max() + 1, 4))
    ax.legend(frameon=False, fontsize=7, ncols=3, loc="upper right")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.4)
    ax.tick_params(axis="y", labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_single_panels(
    matched: pd.DataFrame,
    correlation: pd.DataFrame,
    venn_counts: pd.DataFrame,
    direction_counts: pd.DataFrame,
) -> None:
    """分别导出 Fig. 3 的单个面板图。

    参数:
        matched: mRNA-protein 匹配表。
        correlation: 相关性统计表。
        venn_counts: Venn 数量统计表。
        direction_counts: 方向一致性统计表。

    返回:
        无返回值。
    """
    fig_a, axes_a = plt.subplots(1, 3, figsize=(8.8, 2.8), sharex=False, sharey=False)
    for ax, compare in zip(axes_a, COMPARE_ORDER):
        subset = matched[matched["compare"] == compare]
        stats_row = correlation[correlation["compare"] == COMPARE_LABELS[compare]].iloc[0]
        draw_correlation_panel(ax, subset, stats_row, compare)
    fig_a.tight_layout()
    fig_a.savefig(PANEL_DIR / "Fig3a_mRNA_protein_log2FC_correlation.png", dpi=300)
    fig_a.savefig(PANEL_DIR / "Fig3a_mRNA_protein_log2FC_correlation.tiff", dpi=300)
    plt.close(fig_a)

    fig_b, axes_b = plt.subplots(1, 3, figsize=(7.2, 2.4))
    for ax, compare in zip(axes_b, ["He vs Ctrl", "Ho vs Ctrl", "Ho vs He"]):
        row = venn_counts[venn_counts["compare"] == compare].iloc[0]
        draw_venn_panel(ax, row)
    fig_b.tight_layout()
    fig_b.savefig(PANEL_DIR / "Fig3b_DEG_DEP_overlap_venn.png", dpi=300)
    fig_b.savefig(PANEL_DIR / "Fig3b_DEG_DEP_overlap_venn.tiff", dpi=300)
    plt.close(fig_b)

    fig_c, ax_c = plt.subplots(figsize=(4.8, 3.0))
    draw_direction_panel(ax_c, direction_counts)
    fig_c.tight_layout()
    fig_c.savefig(PANEL_DIR / "Fig3c_overlap_direction_consistency.png", dpi=300)
    fig_c.savefig(PANEL_DIR / "Fig3c_overlap_direction_consistency.tiff", dpi=300)
    plt.close(fig_c)


def draw_combined_figure(
    matched: pd.DataFrame,
    correlation: pd.DataFrame,
    venn_counts: pd.DataFrame,
    direction_counts: pd.DataFrame,
) -> None:
    """导出 Fig. 3 合并图。

    参数:
        matched: mRNA-protein 匹配表。
        correlation: 相关性统计表。
        venn_counts: Venn 数量统计表。
        direction_counts: 方向一致性统计表。

    返回:
        无返回值。
    """
    fig = plt.figure(figsize=(8.8, 8.2))
    grid = fig.add_gridspec(3, 3, height_ratios=[1.15, 0.95, 1.0], hspace=0.42, wspace=0.34)

    axes_a = [fig.add_subplot(grid[0, i]) for i in range(3)]
    for ax, compare in zip(axes_a, COMPARE_ORDER):
        subset = matched[matched["compare"] == compare]
        stats_row = correlation[correlation["compare"] == COMPARE_LABELS[compare]].iloc[0]
        draw_correlation_panel(ax, subset, stats_row, compare)
    add_panel_label(axes_a[0], "(a)")

    axes_b = [fig.add_subplot(grid[1, i]) for i in range(3)]
    for ax, compare in zip(axes_b, ["He vs Ctrl", "Ho vs Ctrl", "Ho vs He"]):
        row = venn_counts[venn_counts["compare"] == compare].iloc[0]
        draw_venn_panel(ax, row)
    add_panel_label(axes_b[0], "(b)")

    ax_c = fig.add_subplot(grid[2, :])
    draw_direction_panel(ax_c, direction_counts)
    add_panel_label(ax_c, "(c)")

    fig.savefig(FIG3_DIR / "Fig3_mRNA_protein_correlation_DEG_DEP_overlap.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG3_DIR / "Fig3_mRNA_protein_correlation_DEG_DEP_overlap.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig)


def copy_source_files() -> None:
    """复制 Fig. 3 相关源数据和既有图件。

    参数:
        无。

    返回:
        无返回值。
    """
    source_files = [
        SUMMARY_FILE,
        MATCHED_TABLE,
        SOURCE_DIR / "combined_mRNA_protein_correlation.png",
        SOURCE_DIR / "combined_mRNA_protein_correlation.tiff",
        SOURCE_DIR / "combined_DEG_DEP_venn.png",
        SOURCE_DIR / "combined_DEG_DEP_venn.tiff",
        SOURCE_DIR / "mRNA_protein_correlation_legend.md",
        SOURCE_DIR / "DEG_DEP_venn_legend.md",
    ]
    for src in source_files:
        if not src.exists():
            continue
        target_dir = DATA_DIR if src.suffix.lower() in [".xlsx", ".csv"] else PANEL_DIR
        shutil.copy2(src, target_dir / src.name)


def write_outputs(correlation: pd.DataFrame, venn_counts: pd.DataFrame, direction_counts: pd.DataFrame) -> None:
    """写出 Fig. 3 统计表和英文图例。

    参数:
        correlation: 相关性统计表。
        venn_counts: Venn 数量统计表。
        direction_counts: 方向一致性统计表。

    返回:
        无返回值。
    """
    stats_path = DATA_DIR / "Fig3_correlation_overlap_direction_summary.xlsx"
    with pd.ExcelWriter(stats_path) as writer:
        correlation.to_excel(writer, sheet_name="correlation", index=False)
        venn_counts.to_excel(writer, sheet_name="venn_counts", index=False)
        direction_counts.to_excel(writer, sheet_name="direction_counts", index=False)

    legend_text = """# Fig. 3 mRNA-protein correlation and DEG-DEP overlap

**Figure legend:** (a) Scatter plots showing the relationship between mRNA log2FC and protein log2FC for all matched mRNA-protein pairs in each pairwise comparison. Grey points indicate non-significant pairs, blue points indicate pairs significant in one omics layer, and red points indicate DEG-DEP overlaps. Pearson and Spearman correlation coefficients are shown in each panel. (b) Venn diagrams showing the overlap between DEGs and DEPs. Circle areas are schematic and are not scaled to set size. (c) Directional consistency of DEG-DEP overlaps, including concordant up-regulation, concordant down-regulation, and opposite regulation between mRNA and protein levels.

## Main files

- `Fig3_mRNA_protein_correlation_DEG_DEP_overlap.png`
- `Fig3_mRNA_protein_correlation_DEG_DEP_overlap.tiff`

## Single panels

- `single_panels/Fig3a_mRNA_protein_log2FC_correlation.png`
- `single_panels/Fig3a_mRNA_protein_log2FC_correlation.tiff`
- `single_panels/Fig3b_DEG_DEP_overlap_venn.png`
- `single_panels/Fig3b_DEG_DEP_overlap_venn.tiff`
- `single_panels/Fig3c_overlap_direction_consistency.png`
- `single_panels/Fig3c_overlap_direction_consistency.tiff`

## Source data

- `source_data/correlation_and_venn_summary.xlsx`
- `source_data/nine_quadrant_matched_gene_protein.csv`
- `source_data/Fig3_correlation_overlap_direction_summary.xlsx`
"""
    (FIG3_DIR / "Fig3_legend.md").write_text(legend_text, encoding="utf-8")


def main() -> None:
    """生成并整理 Fig. 3 图件。

    参数:
        无。

    返回:
        无返回值。
    """
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    matched = read_matched_table()
    correlation, venn_counts, overlap = read_summary_tables()
    direction_counts = build_direction_counts(overlap)

    draw_single_panels(matched, correlation, venn_counts, direction_counts)
    draw_combined_figure(matched, correlation, venn_counts, direction_counts)
    copy_source_files()
    write_outputs(correlation, venn_counts, direction_counts)


if __name__ == "__main__":
    main()
