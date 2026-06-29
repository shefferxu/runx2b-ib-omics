from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.patches import Patch

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
MAIN_CANDIDATE_FILE = OUTPUT_DIR / "main_candidate_heatmap" / "main_core_candidate_heatmap_data.xlsx"
MATCHED_TABLE = OUTPUT_DIR / "nine_quadrant_plots" / "nine_quadrant_matched_gene_protein.csv"

FIG6_DIR = OUTPUT_DIR / "figures" / "Fig6_candidate_expression_heatmap"
PANEL_DIR = FIG6_DIR / "single_panels"
DATA_DIR = FIG6_DIR / "source_data"

GROUP_COLUMNS = ["Ctrl", "He", "Ho"]
THEME_ORDER = ["ECM", "Cytoskeleton", "Calcium", "Muscle", "Ossification"]
THEME_COLORS = {
    "ECM": "#A1A9D0",
    "Cytoskeleton": "#F0988C",
    "Calcium": "#B883D4",
    "Muscle": "#96CCCB",
    "Ossification": "#C4A5DE",
}


def assign_primary_theme(row: pd.Series) -> str:
    """为候选分子分配主功能分类。

    参数:
        row: 候选分子记录。

    返回:
        ECM、Cytoskeleton、Calcium、Muscle 或 Ossification。
    """
    text = " ".join(
        str(row.get(column, ""))
        for column in ["module", "phenotype_relevance", "best_function_description", "best_domain", "best_go_description", "best_kegg"]
    ).lower()
    module = str(row.get("module", "")).lower()
    if "ecm" in module or "adhesion" in module:
        return "ECM"
    if "calcium" in module:
        return "Calcium"
    if "muscle" in module:
        return "Muscle"
    if "ossification" in module:
        return "Ossification"
    if "cytoskeleton" in module:
        return "Cytoskeleton"
    if "ecm" in text or "adhesion" in text or "collagen" in text or "cadherin" in text:
        return "ECM"
    if "calcium" in text or "annexin" in text or "calpain" in text:
        return "Calcium"
    if "muscle" in text or "myosin" in text:
        return "Muscle"
    if "ossification" in text or "bone" in text:
        return "Ossification"
    return "Cytoskeleton"


def short_gene_label(display_label: str) -> str:
    """从 display_label 中提取基因/蛋白简称。

    参数:
        display_label: 候选显示名，如 MYH1s | KAK9975044.1。

    返回:
        str: 简化后的行名。
    """
    label = str(display_label).split("|")[0].strip()
    if label.lower().startswith("phha") or "PAH" in label:
        return "PAH"
    if label == "PLIN5_RAT":
        return "PLIN5"
    return label


def read_candidates() -> pd.DataFrame:
    """读取并整理最终候选表。

    参数:
        无。

    返回:
        pandas.DataFrame: 候选表。
    """
    candidates = pd.read_excel(MAIN_CANDIDATE_FILE, sheet_name="main_heatmap_candidates")
    candidates["primary_theme"] = candidates.apply(assign_primary_theme, axis=1)
    candidates["gene_label"] = candidates["display_label"].map(short_gene_label)
    return candidates


def read_matched_table() -> pd.DataFrame:
    """读取九象限匹配表中的表达均值。

    参数:
        无。

    返回:
        pandas.DataFrame: 匹配表。
    """
    table = pd.read_csv(MATCHED_TABLE)
    for column in [
        "gene_Ctrl", "gene_He", "gene_Ho",
        "protein_Ctrl", "protein_He", "protein_Ho",
        "gene_log2fc", "protein_log2fc",
    ]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    return table


def mean_available(values: pd.Series) -> float:
    """计算非缺失值均值。

    参数:
        values: 数值序列。

    返回:
        float: 非缺失值均值。
    """
    numeric_values = pd.to_numeric(values, errors="coerce").dropna()
    if numeric_values.empty:
        return np.nan
    return float(numeric_values.mean())


def format_row_label(gene_label: str, protein_ids: str) -> str:
    """生成带 KAK 编号的简短热图行名。

    参数:
        gene_label: 基因简称。
        protein_ids: 分号分隔的蛋白编号。

    返回:
        str: 两行显示的行名。
    """
    ids = [item.strip() for item in str(protein_ids).split(";") if item.strip() and item.strip() != "nan"]
    unique_ids = sorted(set(ids))
    if not unique_ids:
        return gene_label
    if len(unique_ids) > 1:
        return f"{gene_label}\n{unique_ids[0]} +{len(unique_ids) - 1}"
    return f"{gene_label}\n{unique_ids[0]}"


def build_gene_level_tables(candidates: pd.DataFrame, matched: pd.DataFrame):
    """按基因简称合并候选并构建表达矩阵。

    参数:
        candidates: 候选表。
        matched: 九象限匹配表。

    返回:
        tuple: mRNA log2 均值矩阵、protein log2 均值矩阵、注释表。
    """
    records = []
    for _, candidate in candidates.iterrows():
        subsets = matched[
            (matched["gene"] == candidate["gene_id"])
            & (matched["protein"] == candidate["protein_id"])
        ]
        if subsets.empty:
            continue
        record = {
            "gene_label": candidate["gene_label"],
            "display_label": candidate["display_label"],
            "gene_id": candidate["gene_id"],
            "protein_id": candidate["protein_id"],
            "primary_theme": candidate["primary_theme"],
            "module": candidate["module"],
            "phenotype_relevance": candidate["phenotype_relevance"],
            "mRNA_log2FC_values": candidate["mRNA_log2FC"],
            "protein_log2FC_values": candidate["protein_log2FC"],
            "mRNA_Ctrl": mean_available(subsets["gene_Ctrl"]),
            "mRNA_He": mean_available(subsets["gene_He"]),
            "mRNA_Ho": mean_available(subsets["gene_Ho"]),
            "protein_Ctrl": mean_available(subsets["protein_Ctrl"]),
            "protein_He": mean_available(subsets["protein_He"]),
            "protein_Ho": mean_available(subsets["protein_Ho"]),
        }
        records.append(record)

    raw = pd.DataFrame(records)
    theme_rank = {theme: idx for idx, theme in enumerate(THEME_ORDER)}
    grouped_records = []
    for gene_label, group in raw.groupby("gene_label", sort=False):
        theme = group["primary_theme"].mode().iloc[0]
        grouped_records.append({
            "gene_label": gene_label,
            "primary_theme": theme,
            "theme_rank": theme_rank.get(theme, 99),
            "display_labels": "; ".join(sorted(set(group["display_label"].astype(str)))),
            "gene_ids": "; ".join(sorted(set(group["gene_id"].astype(str)))),
            "protein_ids": "; ".join(sorted(set(group["protein_id"].astype(str)))),
            "modules": "; ".join(sorted(set(group["module"].astype(str)))),
            "phenotype_relevance": "; ".join(sorted(set(group["phenotype_relevance"].astype(str)))),
            "mRNA_Ctrl": group["mRNA_Ctrl"].mean(skipna=True),
            "mRNA_He": group["mRNA_He"].mean(skipna=True),
            "mRNA_Ho": group["mRNA_Ho"].mean(skipna=True),
            "protein_Ctrl": group["protein_Ctrl"].mean(skipna=True),
            "protein_He": group["protein_He"].mean(skipna=True),
            "protein_Ho": group["protein_Ho"].mean(skipna=True),
            "mRNA_log2FC_mean": group["mRNA_log2FC_values"].mean(skipna=True),
            "protein_log2FC_mean": group["protein_log2FC_values"].mean(skipna=True),
        })

    annotation = pd.DataFrame(grouped_records).sort_values(["theme_rank", "gene_label"]).reset_index(drop=True)
    annotation["row_label"] = annotation.apply(lambda row: format_row_label(row["gene_label"], row["protein_ids"]), axis=1)
    mrna_log = np.log2(annotation[["mRNA_Ctrl", "mRNA_He", "mRNA_Ho"]].astype(float) + 1)
    protein_log = np.log2(annotation[["protein_Ctrl", "protein_He", "protein_Ho"]].astype(float) + 1)
    mrna_log.index = annotation["row_label"]
    protein_log.index = annotation["row_label"]
    mrna_log.columns = GROUP_COLUMNS
    protein_log.columns = GROUP_COLUMNS
    return mrna_log, protein_log, annotation


def row_zscore(matrix: pd.DataFrame) -> pd.DataFrame:
    """计算按三组行标准化的 z-score。

    参数:
        matrix: log2(mean+1) 表达矩阵。

    返回:
        pandas.DataFrame: 行 z-score 矩阵。
    """
    values = matrix.astype(float)
    mean = values.mean(axis=1, skipna=True)
    std = values.std(axis=1, skipna=True).replace(0, np.nan)
    zscore = values.sub(mean, axis=0).div(std, axis=0)
    return zscore.clip(-2, 2)


def add_panel_label(ax, label: str, x_offset: float = -0.12, fontsize: int = 16) -> None:
    """添加面板标签。

    参数:
        ax: 坐标轴。
        label: 标签。

    返回:
        None。
    """
    ax.text(
        x_offset,
        1.07,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fontsize,
        fontweight="bold",
        fontfamily="Times New Roman",
    )


def draw_heatmap(ax, matrix: pd.DataFrame, annotation: pd.DataFrame, title: str, show_y=True):
    """绘制候选表达热图。

    参数:
        ax: 坐标轴。
        matrix: row z-score 矩阵。
        annotation: 候选注释。
        title: 标题。
        show_y: 是否显示 y 轴标签。

    返回:
        matplotlib.image.AxesImage: 热图对象。
    """
    masked = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#D0D0D0")
    image = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-2, vmax=2)
    ax.set_xticks(np.arange(len(GROUP_COLUMNS)))
    ax.set_xticklabels(GROUP_COLUMNS, fontsize=8)
    ax.set_yticks(np.arange(matrix.shape[0]))
    if show_y:
        ax.set_yticklabels(matrix.index.tolist(), fontsize=5.6)
        for tick_label, theme in zip(ax.get_yticklabels(), annotation["primary_theme"]):
            tick_label.set_color(THEME_COLORS.get(theme, "black"))
            tick_label.set_linespacing(0.9)
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
    ax.set_title(title, fontsize=9, pad=5)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for x in np.arange(-0.5, len(GROUP_COLUMNS), 1):
        ax.axvline(x, color="white", linewidth=0.8)
    for y in np.arange(-0.5, matrix.shape[0], 1):
        ax.axhline(y, color="white", linewidth=0.8)
    return image


def draw_logfc_panel(ax, annotation: pd.DataFrame) -> None:
    """绘制候选 mRNA/protein 平均 log2FC 对照图。

    参数:
        ax: 坐标轴。
        annotation: 候选注释。

    返回:
        None。
    """
    plot_table = annotation.iloc[::-1].reset_index(drop=True)
    y_positions = np.arange(len(plot_table))
    for y_value, row in zip(y_positions, plot_table.itertuples()):
        ax.plot([row.mRNA_log2FC_mean, row.protein_log2FC_mean], [y_value, y_value], color="#CFCFCF", linewidth=0.7)
    ax.scatter(plot_table["mRNA_log2FC_mean"], y_positions, s=22, color="#A1A9D0", edgecolor="#333333", linewidth=0.35, label="mRNA log$_2$FC")
    ax.scatter(plot_table["protein_log2FC_mean"], y_positions, s=22, color="#F0988C", marker="s", edgecolor="#333333", linewidth=0.35, label="Protein log$_2$FC")
    ax.axvline(0, color="#333333", linewidth=0.7)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(plot_table["gene_label"], fontsize=7)
    for tick_label, theme in zip(ax.get_yticklabels(), plot_table["primary_theme"]):
        tick_label.set_color(THEME_COLORS.get(theme, "black"))
    ax.set_xlabel("log$_2$FC", fontsize=8)
    ax.set_title("Mean matched mRNA-protein log$_2$FC", fontsize=9, pad=5)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.4)


def write_outputs(mrna_log, protein_log, mrna_z, protein_z, annotation):
    """写出 Fig.6 重绘后的源数据和图例。

    参数:
        mrna_log: mRNA log2(mean+1) 矩阵。
        protein_log: protein log2(mean+1) 矩阵。
        mrna_z: mRNA row z-score 矩阵。
        protein_z: protein row z-score 矩阵。
        annotation: 注释表。

    返回:
        None。
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    annotation.to_csv(DATA_DIR / "Fig6_gene_level_candidate_annotation.csv", index=False, encoding="utf-8-sig")
    mrna_log.to_csv(DATA_DIR / "Fig6_gene_level_mRNA_log2_mean_expression.csv", encoding="utf-8-sig")
    protein_log.to_csv(DATA_DIR / "Fig6_gene_level_protein_log2_mean_expression.csv", encoding="utf-8-sig")
    mrna_z.to_csv(DATA_DIR / "Fig6_gene_level_mRNA_row_zscore.csv", encoding="utf-8-sig")
    protein_z.to_csv(DATA_DIR / "Fig6_gene_level_protein_row_zscore.csv", encoding="utf-8-sig")

    legend = """# Fig. 6. Candidate mRNA-protein expression patterns

**Figure legend:** (a) Heatmap showing row-scaled mRNA expression patterns of gene-level candidate molecules across Ctrl, He, and Ho groups. (b) Heatmap showing row-scaled protein abundance patterns of the same candidate set. For both heatmaps, group mean values were first transformed as log2(mean + 1), followed by row-wise z-score normalization across Ctrl, He, and Ho. Missing values are shown in light gray. (c) Mean matched mRNA and protein log2FC values for each candidate. Candidate rows were collapsed by gene/protein abbreviation, and multiple protein accessions assigned to the same abbreviation were averaged for visualization.
"""
    (FIG6_DIR / "Fig6_legend.md").write_text(legend, encoding="utf-8")


def draw_single_panels(mrna_z, protein_z, annotation):
    """导出 Fig.6 单图面板。

    参数:
        mrna_z: mRNA row z-score。
        protein_z: protein row z-score。
        annotation: 注释表。

    返回:
        None。
    """
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    fig_a, ax_a = plt.subplots(figsize=(3.2, 4.6), dpi=300)
    img = draw_heatmap(ax_a, mrna_z, annotation, "mRNA expression", show_y=True)
    fig_a.colorbar(img, ax=ax_a, fraction=0.045, pad=0.02).set_label("Row z-score", fontsize=7)
    fig_a.savefig(PANEL_DIR / "Fig6a_candidate_mRNA_heatmap.png", dpi=300, bbox_inches="tight")
    fig_a.savefig(PANEL_DIR / "Fig6a_candidate_mRNA_heatmap.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig_a)

    fig_b, ax_b = plt.subplots(figsize=(3.2, 4.6), dpi=300)
    img = draw_heatmap(ax_b, protein_z, annotation, "Protein abundance", show_y=True)
    fig_b.colorbar(img, ax=ax_b, fraction=0.045, pad=0.02).set_label("Row z-score", fontsize=7)
    fig_b.savefig(PANEL_DIR / "Fig6b_candidate_protein_heatmap.png", dpi=300, bbox_inches="tight")
    fig_b.savefig(PANEL_DIR / "Fig6b_candidate_protein_heatmap.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig_b)

    fig_c, ax_c = plt.subplots(figsize=(4.6, 4.6), dpi=300)
    draw_logfc_panel(ax_c, annotation)
    fig_c.savefig(PANEL_DIR / "Fig6c_candidate_mRNA_protein_log2FC.png", dpi=300, bbox_inches="tight")
    fig_c.savefig(PANEL_DIR / "Fig6c_candidate_mRNA_protein_log2FC.tiff", dpi=300, bbox_inches="tight")
    plt.close(fig_c)


def draw_combined(mrna_z, protein_z, annotation):
    """绘制 Fig.6 合并图。

    参数:
        mrna_z: mRNA row z-score。
        protein_z: protein row z-score。
        annotation: 注释表。

    返回:
        None。
    """
    FIG6_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(13.2, 5.6), dpi=300)
    gs = fig.add_gridspec(1, 5, width_ratios=[1.10, 1.10, 0.07, 0.24, 1.45], wspace=0.18)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])
    ax_c = fig.add_subplot(gs[0, 4])

    image = draw_heatmap(ax_a, mrna_z, annotation, "mRNA expression", show_y=True)
    draw_heatmap(ax_b, protein_z, annotation, "Protein abundance", show_y=False)
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label("Row z-score", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)
    draw_logfc_panel(ax_c, annotation)

    add_panel_label(ax_a, "a", x_offset=-0.25, fontsize=32)
    add_panel_label(ax_b, "b", x_offset=-0.12, fontsize=32)
    add_panel_label(ax_c, "c", x_offset=-0.18, fontsize=32)

    theme_handles = [Patch(facecolor=color, edgecolor="none", label=theme) for theme, color in THEME_COLORS.items()]
    fig.legend(handles=theme_handles, loc="lower center", bbox_to_anchor=(0.56, 0.005), ncol=5, frameon=False, fontsize=7)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.88, bottom=0.19)
    fig.savefig(FIG6_DIR / "Fig6_candidate_expression_heatmap.png", dpi=300)
    fig.savefig(FIG6_DIR / "Fig6_candidate_expression_heatmap.tiff", dpi=300)
    plt.close(fig)


def main():
    """按基因简称重绘 Fig.6 候选表达热图。

    参数:
        无。

    返回:
        None。
    """
    candidates = read_candidates()
    matched = read_matched_table()
    mrna_log, protein_log, annotation = build_gene_level_tables(candidates, matched)
    mrna_z = row_zscore(mrna_log)
    protein_z = row_zscore(protein_log)
    write_outputs(mrna_log, protein_log, mrna_z, protein_z, annotation)
    draw_single_panels(mrna_z, protein_z, annotation)
    draw_combined(mrna_z, protein_z, annotation)


if __name__ == "__main__":
    main()
