from pathlib import Path
import shutil

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
SOURCE_DIR = OUTPUT_DIR / "nine_quadrant_plots"
MATCHED_TABLE = SOURCE_DIR / "nine_quadrant_matched_gene_protein.csv"
SOURCE_SUMMARY = SOURCE_DIR / "nine_quadrant_summary_and_candidates.xlsx"
PROTEIN_ANNOTATION_FILE = ROOT_DIR / (
    "MJ20251014053-ZX-P-251019025-郑建波-高通量蛋白组学-9个样本-翘嘴鲌-尾部肌肉-20251117-FX2025102300143"
) / "翘嘴鲌-尾部肌肉-9" / "workflow_results" / "03_Annotation" / "Stat" / "proteins_anno_detail.xls"

FIG4_DIR = OUTPUT_DIR / "figures" / "Fig4_nine_quadrant"
PANEL_DIR = FIG4_DIR / "single_panels"
DATA_DIR = FIG4_DIR / "source_data"

GENE_UP_CUTOFF = 1.0
GENE_DOWN_CUTOFF = -1.0
PROTEIN_UP_CUTOFF = np.log2(1.2)
PROTEIN_DOWN_CUTOFF = np.log2(1 / 1.2)

COMPARE_ORDER = ["He|Ctrl", "Ho|Ctrl", "Ho|He"]
COMPARE_LABELS = {
    "He|Ctrl": "He vs Ctrl",
    "Ho|Ctrl": "Ho vs Ctrl",
    "Ho|He": "Ho vs He",
}

PATTERN_ORDER = [
    "Genes(up)_Proteins(up)",
    "Genes(up)_Proteins(no)",
    "Genes(up)_Proteins(down)",
    "Genes(no)_Proteins(up)",
    "Genes(no)_Proteins(no)",
    "Genes(no)_Proteins(down)",
    "Genes(down)_Proteins(up)",
    "Genes(down)_Proteins(no)",
    "Genes(down)_Proteins(down)",
]

PATTERN_LABEL = {
    "Genes(up)_Proteins(up)": "mRNA up / protein up",
    "Genes(up)_Proteins(no)": "mRNA up / protein unchanged",
    "Genes(up)_Proteins(down)": "mRNA up / protein down",
    "Genes(no)_Proteins(up)": "mRNA unchanged / protein up",
    "Genes(no)_Proteins(no)": "mRNA unchanged / protein unchanged",
    "Genes(no)_Proteins(down)": "mRNA unchanged / protein down",
    "Genes(down)_Proteins(up)": "mRNA down / protein up",
    "Genes(down)_Proteins(no)": "mRNA down / protein unchanged",
    "Genes(down)_Proteins(down)": "mRNA down / protein down",
}

PATTERN_SHORT = {
    "Genes(up)_Proteins(up)": "Up / Up",
    "Genes(up)_Proteins(no)": "Up / Unchanged",
    "Genes(up)_Proteins(down)": "Up / Down",
    "Genes(no)_Proteins(up)": "Unchanged / Up",
    "Genes(no)_Proteins(no)": "Unchanged / Unchanged",
    "Genes(no)_Proteins(down)": "Unchanged / Down",
    "Genes(down)_Proteins(up)": "Down / Up",
    "Genes(down)_Proteins(no)": "Down / Unchanged",
    "Genes(down)_Proteins(down)": "Down / Down",
}

PATTERN_COLOR = {
    "Genes(up)_Proteins(up)": "#F0988C",
    "Genes(down)_Proteins(down)": "#A1A9D0",
    "Genes(up)_Proteins(down)": "#B883D4",
    "Genes(down)_Proteins(up)": "#96CCCB",
    "Genes(up)_Proteins(no)": "#F6CAE5",
    "Genes(down)_Proteins(no)": "#C4A5DE",
    "Genes(no)_Proteins(up)": "#CFEAF1",
    "Genes(no)_Proteins(down)": "#A8DADC",
    "Genes(no)_Proteins(no)": "#BDBDBD",
}

FOCUS_PATTERNS = [
    "Genes(up)_Proteins(up)",
    "Genes(down)_Proteins(down)",
    "Genes(no)_Proteins(up)",
    "Genes(no)_Proteins(down)",
    "Genes(up)_Proteins(no)",
    "Genes(down)_Proteins(no)",
]

LEGEND_PATTERNS = [
    "Genes(up)_Proteins(up)",
    "Genes(down)_Proteins(down)",
    "Genes(up)_Proteins(down)",
    "Genes(down)_Proteins(up)",
    "Genes(no)_Proteins(up)",
    "Genes(no)_Proteins(down)",
    "Genes(up)_Proteins(no)",
    "Genes(down)_Proteins(no)",
    "Genes(no)_Proteins(no)",
]

MECHANISM_KEYWORDS = {
    "muscle": ["muscle", "myosin", "skeletal", "actin", "troponin", "tropomyosin"],
    "ECM": ["extracellular", "matrix", "collagen", "microfibril", "fibronectin", "laminin", "integrin"],
    "cytoskeleton": ["cytoskeleton", "keratin", "tubulin", "actinin", "plakin", "envoplakin"],
    "calcium": ["calcium", "calmodulin", "calc", "inositol", "phosphatase", "phospholipase"],
    "ossification": ["ossification", "bone", "mineral", "osteoblast", "osteoclast", "osteogenesis"],
}


def read_matched_table() -> pd.DataFrame:
    """读取并清洗九象限 mRNA-protein 匹配表。

    参数:
        无。

    返回:
        清洗后的九象限数据表。
    """
    table = pd.read_csv(MATCHED_TABLE)
    table["gene_log2fc"] = pd.to_numeric(table["gene_log2fc"], errors="coerce")
    table["protein_log2fc"] = pd.to_numeric(table["protein_log2fc"], errors="coerce")
    table["gene_significant"] = table["gene_significant"].astype(str).str.lower()
    table["protein_significant"] = table["protein_significant"].astype(str).str.lower()
    table = table.replace([np.inf, -np.inf], np.nan)
    return table.dropna(subset=["gene_log2fc", "protein_log2fc"])


def extract_nr_name(nr_value: object) -> str:
    """从 NR 注释中提取简短蛋白名称。

    参数:
        nr_value: NR 注释字段。

    返回:
        简短蛋白名称。
    """
    if pd.isna(nr_value):
        return ""
    text = str(nr_value)
    if "(" in text and ")" in text:
        text = text.split("(", 1)[1].rsplit(")", 1)[0]
    text = text.replace(" isoform X1", "").replace(" isoform X2", "")
    text = text.replace(" [Megalobrama amblycephala]", "")
    text = text.replace(" [Ctenopharyngodon idella]", "")
    return text.strip()


def read_protein_annotations() -> pd.DataFrame:
    """读取蛋白注释表。

    参数:
        无。

    返回:
        蛋白编号及简短注释表。
    """
    annotation = pd.read_csv(PROTEIN_ANNOTATION_FILE, sep="\t")
    annotation["protein_annotation"] = annotation["nr"].map(extract_nr_name)
    annotation["protein_annotation"] = annotation["protein_annotation"].replace("", np.nan)
    fallback = annotation["KO_name"].where(annotation["KO_name"].notna(), annotation["EggNOG_description"])
    annotation["protein_annotation"] = annotation["protein_annotation"].fillna(fallback)
    annotation["short_annotation"] = annotation["KO_name"].where(annotation["KO_name"].notna(), annotation["protein_annotation"])
    annotation["short_annotation"] = annotation["short_annotation"].fillna("unannotated")
    return annotation[["Accession_id", "protein_annotation", "short_annotation", "KO_name", "nr"]]


def assign_mechanism_terms(row: pd.Series) -> str:
    """根据注释关键词判定候选分子的机制主题。

    参数:
        row: 候选分子注释行。

    返回:
        命中的机制主题，多个主题用分号连接。
    """
    search_text = " ".join(
        str(row.get(column, ""))
        for column in ["protein", "short_annotation", "protein_annotation", "KO_name", "nr"]
    ).lower()
    matched_terms = []
    for theme, keywords in MECHANISM_KEYWORDS.items():
        if any(keyword in search_text for keyword in keywords):
            matched_terms.append(theme)
    return "; ".join(matched_terms)


def clip_limits(values: pd.Series, min_abs: float = 2.0) -> tuple[float, float]:
    """按分位数计算坐标范围。

    参数:
        values: 数值序列。
        min_abs: 最小绝对坐标范围。

    返回:
        坐标轴下限和上限。
    """
    clean_values = values.dropna()
    lower = float(clean_values.quantile(0.005))
    upper = float(clean_values.quantile(0.995))
    limit = max(abs(lower), abs(upper), min_abs)
    return -limit * 1.08, limit * 1.08


def add_panel_label(ax, label: str) -> None:
    """为合并图添加面板标签。

    参数:
        ax: matplotlib 坐标轴。
        label: 面板标签。

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


def select_candidates(table: pd.DataFrame, max_per_compare: int = 4) -> pd.DataFrame:
    """选择需要在九象限图中标注的核心候选。

    参数:
        table: 九象限数据表。
        max_per_compare: 每个比较最多标注数量。

    返回:
        候选基因/蛋白数据表。
    """
    candidates = table[
        (table["gene_significant"] == "yes")
        & (table["protein_significant"] == "yes")
    ].copy()
    candidates["compare_label"] = candidates["compare"].map(COMPARE_LABELS)
    candidates["pattern"] = candidates["label_type_2"].map(PATTERN_LABEL)
    candidates["rank_score"] = candidates["gene_log2fc"].abs() + candidates["protein_log2fc"].abs()
    annotations = read_protein_annotations()
    candidates = candidates.merge(annotations, left_on="protein", right_on="Accession_id", how="left")
    candidates["protein_annotation"] = candidates["protein_annotation"].fillna("unannotated protein")
    candidates["short_annotation"] = candidates["short_annotation"].fillna(candidates["protein_annotation"])
    candidates["mechanism_terms"] = candidates.apply(assign_mechanism_terms, axis=1)
    candidates["mechanism_hit"] = candidates["mechanism_terms"].ne("")
    candidates = candidates.sort_values(
        ["compare", "mechanism_hit", "rank_score"],
        ascending=[True, False, False],
    )
    candidates = candidates.groupby("compare", group_keys=False).head(max_per_compare)
    return candidates


def build_panel_a_legend_handles() -> tuple[list[Line2D], list[str]]:
    """生成 Fig. 4a 的颜色图例。

    参数:
        无。

    返回:
        图例句柄和标签。
    """
    handles = []
    labels = []
    for pattern in LEGEND_PATTERNS:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markersize=4.5,
                markerfacecolor=PATTERN_COLOR[pattern],
                markeredgecolor="none",
                alpha=0.85,
            )
        )
        labels.append(PATTERN_SHORT[pattern])
    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=5.8,
            markerfacecolor="#FFD166",
            markeredgecolor="#333333",
            alpha=1.0,
        )
    )
    labels.append("Highlighted candidate")
    return handles, labels


def build_summary(table: pd.DataFrame) -> pd.DataFrame:
    """生成每个比较的九象限数量统计。

    参数:
        table: 九象限数据表。

    返回:
        九象限数量统计表。
    """
    records = []
    for compare in COMPARE_ORDER:
        subset = table[table["compare"] == compare]
        counts = subset["label_type_2"].value_counts()
        for pattern in PATTERN_ORDER:
            records.append({
                "compare": COMPARE_LABELS[compare],
                "raw_pattern": pattern,
                "pattern": PATTERN_LABEL[pattern],
                "short_pattern": PATTERN_SHORT[pattern],
                "count": int(counts.get(pattern, 0)),
            })
    return pd.DataFrame(records)


def draw_nine_quadrant_panel(
    ax,
    subset: pd.DataFrame,
    compare: str,
    candidates: pd.DataFrame | None = None,
) -> None:
    """绘制单个比较的九象限散点图。

    参数:
        ax: matplotlib 坐标轴。
        subset: 单个比较组的数据。
        compare: 比较名称。
        candidates: 需要标注的候选表。

    返回:
        无返回值。
    """
    for pattern in PATTERN_ORDER:
        pattern_subset = subset[subset["label_type_2"] == pattern]
        if pattern_subset.empty:
            continue
        alpha = 0.22 if pattern == "Genes(no)_Proteins(no)" else 0.68
        size = 5 if pattern == "Genes(no)_Proteins(no)" else 9
        ax.scatter(
            pattern_subset["gene_log2fc"],
            pattern_subset["protein_log2fc"],
            s=size,
            color=PATTERN_COLOR[pattern],
            alpha=alpha,
            linewidths=0,
            rasterized=True,
        )

    ax.axvline(GENE_DOWN_CUTOFF, color="#555555", linestyle="--", linewidth=0.8)
    ax.axvline(GENE_UP_CUTOFF, color="#555555", linestyle="--", linewidth=0.8)
    ax.axhline(PROTEIN_DOWN_CUTOFF, color="#555555", linestyle="--", linewidth=0.8)
    ax.axhline(PROTEIN_UP_CUTOFF, color="#555555", linestyle="--", linewidth=0.8)
    ax.set_title(COMPARE_LABELS[compare], fontsize=9, pad=4)
    ax.set_xlabel("mRNA log2FC", fontsize=8)
    ax.set_ylabel("Protein log2FC", fontsize=8)
    ax.grid(True, color="#E6E6E6", linewidth=0.4)
    ax.tick_params(labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(*clip_limits(subset["gene_log2fc"]))
    ax.set_ylim(*clip_limits(subset["protein_log2fc"]))

    counts = subset["label_type_2"].value_counts()
    text = (
        f"Up/Up: {int(counts.get('Genes(up)_Proteins(up)', 0))}\n"
        f"Down/Down: {int(counts.get('Genes(down)_Proteins(down)', 0))}\n"
        f"Unchanged/changed: "
        f"{int(counts.get('Genes(no)_Proteins(up)', 0) + counts.get('Genes(no)_Proteins(down)', 0))}\n"
        f"Changed/unchanged: "
        f"{int(counts.get('Genes(up)_Proteins(no)', 0) + counts.get('Genes(down)_Proteins(no)', 0))}"
    )
    ax.text(0.02, 0.98, text, transform=ax.transAxes, ha="left", va="top", fontsize=6.6)

    if candidates is not None and not candidates.empty:
        candidate_subset = candidates[candidates["compare"] == compare]
        for _, row in candidate_subset.iterrows():
            ax.scatter(
                row["gene_log2fc"],
                row["protein_log2fc"],
                s=34,
                color="#FFD166",
                edgecolor="#333333",
                linewidth=0.5,
                zorder=5,
            )


def draw_summary_heatmap(ax, summary: pd.DataFrame) -> None:
    """绘制九象限数量热图。

    参数:
        ax: matplotlib 坐标轴。
        summary: 九象限数量统计表。

    返回:
        无返回值。
    """
    matrix = summary.pivot(index="short_pattern", columns="compare", values="count")
    matrix = matrix.loc[[PATTERN_SHORT[p] for p in PATTERN_ORDER], [COMPARE_LABELS[c] for c in COMPARE_ORDER]]
    log_matrix = np.log10(matrix + 1)
    image = ax.imshow(log_matrix, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, fontsize=7)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=6.5)
    ax.tick_params(length=0)
    for row_idx in range(matrix.shape[0]):
        for col_idx in range(matrix.shape[1]):
            value = int(matrix.iloc[row_idx, col_idx])
            color = "white" if log_matrix.iloc[row_idx, col_idx] > log_matrix.to_numpy().max() * 0.55 else "#222222"
            ax.text(col_idx, row_idx, str(value), ha="center", va="center", fontsize=6.5, color=color)
    ax.set_title("Nine-quadrant counts", fontsize=9, pad=5)
    colorbar = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
    colorbar.ax.tick_params(labelsize=6)
    colorbar.set_label("log10(count + 1)", fontsize=7)


def draw_candidate_panel(ax, candidates: pd.DataFrame) -> None:
    """绘制核心候选基因/蛋白 mRNA-protein 双点图。

    参数:
        ax: matplotlib 坐标轴。
        candidates: 候选基因/蛋白表。

    返回:
        无返回值。
    """
    plot_table = candidates.copy()
    plot_table["candidate_label"] = plot_table["symbol"].fillna(plot_table["protein"]).astype(str)
    plot_table["display_label"] = (
        plot_table["protein"].astype(str)
        + " | "
        + plot_table["short_annotation"].astype(str)
    )
    plot_table["compare_label"] = plot_table["compare"].map(COMPARE_LABELS)
    plot_table["compare_short"] = plot_table["compare_label"].str.replace(" vs ", "/", regex=False)
    plot_table = plot_table.sort_values(["compare", "rank_score"], ascending=[True, False])
    plot_table["row_label"] = plot_table["display_label"] + " | " + plot_table["compare_short"]
    y_positions = np.arange(len(plot_table))

    for y_value, (_, row) in zip(y_positions, plot_table.iterrows()):
        ax.plot(
            [row["gene_log2fc"], row["protein_log2fc"]],
            [y_value, y_value],
            color="#BDBDBD",
            linewidth=0.75,
            zorder=1,
        )
    ax.scatter(
        plot_table["gene_log2fc"],
        y_positions,
        s=34,
        color="#A1A9D0",
        edgecolor="#333333",
        linewidth=0.45,
        label="mRNA log2FC",
        zorder=3,
    )
    ax.scatter(
        plot_table["protein_log2fc"],
        y_positions,
        s=34,
        color="#F0988C",
        marker="s",
        edgecolor="#333333",
        linewidth=0.45,
        label="Protein log2FC",
        zorder=3,
    )
    ax.axvline(0, color="#777777", linestyle="--", linewidth=0.7)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(plot_table["row_label"], fontsize=5.8)
    ax.set_xlabel("log2FC", fontsize=8)
    ax.set_title("Highlighted DEG-DEP candidates", fontsize=9, pad=5)
    ax.tick_params(axis="x", labelsize=7)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.4)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    min_value = min(plot_table["gene_log2fc"].min(), plot_table["protein_log2fc"].min())
    max_value = max(plot_table["gene_log2fc"].max(), plot_table["protein_log2fc"].max())
    ax.set_xlim(min_value - 0.7, max_value + 0.7)
    ax.invert_yaxis()


def draw_single_panels(table: pd.DataFrame, summary: pd.DataFrame, candidates: pd.DataFrame) -> None:
    """导出 Fig. 4 单图。

    参数:
        table: 九象限数据表。
        summary: 九象限数量统计表。
        candidates: 候选基因/蛋白表。

    返回:
        无返回值。
    """
    fig_a, axes_a = plt.subplots(1, 3, figsize=(8.6, 3.25))
    for ax, compare in zip(axes_a, COMPARE_ORDER):
        draw_nine_quadrant_panel(ax, table[table["compare"] == compare], compare, candidates)
    handles, labels = build_panel_a_legend_handles()
    fig_a.legend(
        handles,
        labels,
        frameon=False,
        fontsize=6.2,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=5,
    )
    fig_a.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.25, wspace=0.34)
    fig_a.savefig(PANEL_DIR / "Fig4a_nine_quadrant_scatter.png", dpi=300, bbox_inches="tight")
    plt.close(fig_a)

    fig_b, ax_b = plt.subplots(figsize=(5.8, 3.4))
    draw_summary_heatmap(ax_b, summary)
    fig_b.tight_layout()
    fig_b.savefig(PANEL_DIR / "Fig4b_nine_quadrant_counts.png", dpi=300)
    plt.close(fig_b)

    fig_c, ax_c = plt.subplots(figsize=(5.8, 3.4))
    draw_candidate_panel(ax_c, candidates)
    fig_c.tight_layout()
    fig_c.savefig(PANEL_DIR / "Fig4c_highlighted_candidates.png", dpi=300)
    plt.close(fig_c)


def draw_combined_figure(table: pd.DataFrame, summary: pd.DataFrame, candidates: pd.DataFrame) -> None:
    """导出 Fig. 4 合并图。

    参数:
        table: 九象限数据表。
        summary: 九象限数量统计表。
        candidates: 候选基因/蛋白表。

    返回:
        无返回值。
    """
    fig = plt.figure(figsize=(11.2, 6.4))
    grid = fig.add_gridspec(
        3,
        12,
        height_ratios=[1.05, 0.18, 1.15],
        hspace=0.48,
        wspace=0.60,
    )

    axes_a = [fig.add_subplot(grid[0, 0:4]), fig.add_subplot(grid[0, 4:8]), fig.add_subplot(grid[0, 8:12])]
    for ax, compare in zip(axes_a, COMPARE_ORDER):
        draw_nine_quadrant_panel(ax, table[table["compare"] == compare], compare, candidates)
    add_panel_label(axes_a[0], "(a)")
    handles, labels = build_panel_a_legend_handles()
    legend_ax = fig.add_subplot(grid[1, :])
    legend_ax.axis("off")
    legend_ax.legend(
        handles,
        labels,
        frameon=False,
        fontsize=6.2,
        loc="center",
        ncol=5,
        bbox_to_anchor=(0.5, -0.12),
    )

    ax_b = fig.add_subplot(grid[2, 0:5])
    draw_summary_heatmap(ax_b, summary)
    add_panel_label(ax_b, "(b)")

    ax_c = fig.add_subplot(grid[2, 7:12])
    draw_candidate_panel(ax_c, candidates)
    add_panel_label(ax_c, "(c)")

    fig.savefig(FIG4_DIR / "Fig4_nine_quadrant_analysis.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_outputs(summary: pd.DataFrame, candidates: pd.DataFrame) -> None:
    """保存 Fig. 4 源数据、统计表和英文图例。

    参数:
        summary: 九象限数量统计表。
        candidates: 候选基因/蛋白表。

    返回:
        无返回值。
    """
    summary.to_csv(DATA_DIR / "Fig4_nine_quadrant_counts.csv", index=False, encoding="utf-8-sig")
    candidate_columns = [
        "compare", "compare_label", "gene", "protein", "short_annotation", "protein_annotation",
        "mechanism_terms", "mechanism_hit", "KO_name", "nr", "gene_log2fc", "protein_log2fc", "label_type_2",
        "rank_score",
    ]
    candidates.to_csv(
        DATA_DIR / "Fig4_highlighted_candidates.csv",
        index=False,
        encoding="utf-8-sig",
        columns=[column for column in candidate_columns if column in candidates.columns],
    )

    for src in [MATCHED_TABLE, SOURCE_SUMMARY]:
        if src.exists():
            shutil.copy2(src, DATA_DIR / src.name)

    legend_text = """# Fig. 4 Nine-quadrant analysis

**Figure legend:** (a) Nine-quadrant scatter plots showing joint mRNA-protein regulation patterns. The x-axis represents mRNA log2FC and the y-axis represents protein log2FC. Vertical dashed lines indicate mRNA log2FC thresholds of -1 and 1, and horizontal dashed lines indicate protein log2FC thresholds of log2(1/1.2) and log2(1.2). Highlighted points indicate candidate DEG-DEP overlaps prioritized by mechanism-related annotations, including muscle, ECM, cytoskeleton, calcium, and ossification themes, followed by joint fold-change magnitude. (b) Counts of matched gene-protein pairs in each of the nine expression quadrants. (c) Highlighted candidate gene-protein pairs showing joint transcriptomic and proteomic changes.

## Main files

- `Fig4_nine_quadrant_analysis.png`

## Single panels

- `single_panels/Fig4a_nine_quadrant_scatter.png`
- `single_panels/Fig4b_nine_quadrant_counts.png`
- `single_panels/Fig4c_highlighted_candidates.png`

## Source data

- `source_data/nine_quadrant_matched_gene_protein.csv`
- `source_data/nine_quadrant_summary_and_candidates.xlsx`
- `source_data/Fig4_nine_quadrant_counts.csv`
- `source_data/Fig4_highlighted_candidates.csv`
"""
    (FIG4_DIR / "Fig4_legend.md").write_text(legend_text, encoding="utf-8")


def remove_stale_tiff_files() -> None:
    """删除 Fig. 4 目录中旧的 TIFF 文件。

    参数:
        无。

    返回:
        无返回值。
    """
    for path in FIG4_DIR.rglob("*.tiff"):
        path.unlink()
    for path in FIG4_DIR.rglob("*.tif"):
        path.unlink()


def main() -> None:
    """生成并整理 Fig. 4 九象限分析图。

    参数:
        无。

    返回:
        无返回值。
    """
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    remove_stale_tiff_files()

    table = read_matched_table()
    summary = build_summary(table)
    candidates = select_candidates(table)

    draw_single_panels(table, summary, candidates)
    draw_combined_figure(table, summary, candidates)
    write_outputs(summary, candidates)


if __name__ == "__main__":
    main()
