from pathlib import Path
import math
import re
import textwrap

import matplotlib.pyplot as plt
from cycler import cycler

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

import numpy as np
import pandas as pd


OUTPUT_DIR = Path("output") / "theme_enrichment_reanalysis"
COMPARES = ["He_vs_Ctrl", "Ho_vs_Ctrl", "Ho_vs_He"]
COMPARE_BAR = {"He_vs_Ctrl": "He|Ctrl", "Ho_vs_Ctrl": "Ho|Ctrl", "Ho_vs_He": "Ho|He"}
COMPARE_TITLE = {"He_vs_Ctrl": "He vs Ctrl", "Ho_vs_Ctrl": "Ho vs Ctrl", "Ho_vs_He": "Ho vs He"}

THEME_PATTERNS = {
    "ECM/cell adhesion": [
        "ecm-receptor", "ecm receptor", "extracellular matrix", "collagen",
        "cell adhesion", "adhesion molecule", "focal adhesion", "tight junction",
        "anchoring junction", "cell-substrate adhesion",
    ],
    "Cytoskeleton": [
        "actin", "cytoskeleton", "intermediate filament", "microtubule",
        "myofibril", "sarcomere",
    ],
    "MAPK/PI3K/Wnt/mTOR": [
        "mapk", "pi3k-akt", "pi3k akt", "wnt", "mtor", "tgf-beta", "tgf beta", "bmp",
    ],
    "Calcium signaling": [
        "calcium", "ca2", "ca(2", "calmodulin", "camk", "calcineurin",
    ],
    "Muscle contraction": [
        "muscle contraction", "muscle", "myosin", "troponin", "ryanodine",
        "contractile", "actomyosin",
    ],
    "Bone/ossification": [
        "bone", "ossification", "osteoblast", "osteoclast", "chondrocyte",
        "cartilage", "skeletal", "mineralization", "runx",
    ],
}


def read_tsv_or_excel(path):
    """读取制表符文本或 Excel 结果表。

    参数:
        path: 文件路径。

    返回:
        pandas DataFrame。
    """
    if path.suffix.lower() == ".xlsx":
        return pd.read_excel(path)
    try:
        return pd.read_csv(path, sep="\t")
    except UnicodeDecodeError:
        return pd.read_csv(path, sep="\t", encoding="gbk")


def find_one(predicate):
    """按条件查找单个文件。

    参数:
        predicate: 接收 Path 并返回布尔值的函数。

    返回:
        第一个匹配文件 Path。
    """
    for path in Path(".").rglob("*"):
        if path.is_file() and predicate(path):
            return path
    return None


def parse_ratio_count(value):
    """从 5/137 形式的比例中提取分子数量。

    参数:
        value: 比例字符串或数值。

    返回:
        分子数量。
    """
    if pd.isna(value):
        return np.nan
    text = str(value)
    if "/" in text:
        try:
            return float(text.split("/")[0])
        except ValueError:
            return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def identify_theme(term):
    """根据通路或 GO 名称识别主题类别。

    参数:
        term: GO/KEGG 名称。

    返回:
        匹配到的主题字符串。
    """
    text = str(term).lower()
    hits = []
    for theme, patterns in THEME_PATTERNS.items():
        if any(pattern in text for pattern in patterns):
            hits.append(theme)
    return "; ".join(hits)


def normalize_enrichment(df, source, omics, db, compare, term_col, id_col=None, p_col=None, padj_col=None,
                         count_col=None, ratio_col=None, rich_col=None, gene_list_col=None):
    """统一不同来源的富集表字段。

    参数:
        df: 原始富集表。
        source: 结果来源。
        omics: DEG/DEP/Joint 等层级。
        db: GO 或 KEGG。
        compare: 比较组。
        term_col: term 名称列。
        id_col: ID 列。
        p_col: P 值列。
        padj_col: 校正 P 值列。
        count_col: 分子数量列。
        ratio_col: 比例列。
        rich_col: 富集因子列。
        gene_list_col: 分子列表列。

    返回:
        标准化 DataFrame。
    """
    out = pd.DataFrame(index=df.index)
    out["source"] = source
    out["omics"] = omics
    out["database"] = db
    out["compare"] = compare
    out["compare_label"] = COMPARE_TITLE.get(compare, compare)
    out["term"] = df[term_col].astype(str)
    out["term_id"] = df[id_col].astype(str) if id_col and id_col in df.columns else ""
    out["pvalue"] = pd.to_numeric(df[p_col], errors="coerce") if p_col and p_col in df.columns else np.nan
    out["padjust"] = pd.to_numeric(df[padj_col], errors="coerce") if padj_col and padj_col in df.columns else np.nan
    if count_col and count_col in df.columns:
        out["count"] = pd.to_numeric(df[count_col], errors="coerce")
    elif ratio_col and ratio_col in df.columns:
        out["count"] = df[ratio_col].apply(parse_ratio_count)
    else:
        out["count"] = np.nan
    out["rich_factor"] = pd.to_numeric(df[rich_col], errors="coerce") if rich_col and rich_col in df.columns else np.nan
    out["molecules"] = df[gene_list_col].astype(str) if gene_list_col and gene_list_col in df.columns else ""
    out["theme"] = out["term"].apply(identify_theme)
    p_for_plot = out["padjust"].where(out["padjust"].notna(), out["pvalue"]).clip(lower=1e-300)
    out["neg_log10_p"] = -np.log10(p_for_plot)
    return out


def load_deg_enrichment():
    """读取转录组 DEG GO/KEGG 富集表。

    参数:
        无。

    返回:
        标准化富集表列表。
    """
    records = []
    for compare in COMPARES:
        go_path = find_one(lambda p: p.name == "go_enrich_stat.xlsx" and f"{compare}_G" in str(p))
        kegg_path = find_one(lambda p: p.name == "kegg_enrich_stat.xls" and f"{compare}_G" in str(p))
        if go_path:
            df = read_tsv_or_excel(go_path)
            records.append(normalize_enrichment(
                df, "DEG GO enrichment", "DEG", "GO", compare,
                term_col="Description", id_col="GO ID", p_col="Pvalue", padj_col="Padjust",
                count_col="Number", ratio_col="Ratio_in_study", rich_col="Rich Factor", gene_list_col="Gene IDs",
            ))
        if kegg_path:
            df = read_tsv_or_excel(kegg_path)
            records.append(normalize_enrichment(
                df, "DEG KEGG enrichment", "DEG", "KEGG", compare,
                term_col="Description", id_col="Pathway ID", p_col="Pvalue", padj_col="Padjust",
                count_col="Number", ratio_col="Ratio_in_study", rich_col="Rich Factor", gene_list_col="Gene IDs",
            ))
    return records


def load_dep_enrichment():
    """读取蛋白组 DEP GO/KEGG 富集表。

    参数:
        无。

    返回:
        标准化富集表列表。
    """
    records = []
    for compare in COMPARES:
        go_path = find_one(lambda p: p.name == f"{compare}_all_go_enrich_protein.xlsx")
        kegg_path = find_one(lambda p: p.name == f"{compare}_all_kegg_enrichment_stat.xls")
        if go_path:
            df = read_tsv_or_excel(go_path)
            records.append(normalize_enrichment(
                df, "DEP GO enrichment", "DEP", "GO", compare,
                term_col="discription", id_col="go_id", p_col="Pvalue", padj_col="Adjusted Pvalue",
                count_col="Protein Number", ratio_col="ratio_in_study", rich_col="Rich factor",
                gene_list_col="accession_list",
            ))
        if kegg_path:
            df = read_tsv_or_excel(kegg_path)
            records.append(normalize_enrichment(
                df, "DEP KEGG enrichment", "DEP", "KEGG", compare,
                term_col="KEGG Description", id_col="Pathway_ID", p_col="P-Value", padj_col="Adjusted Pvalue",
                count_col="Protein Number", ratio_col="Ratio_in_study", rich_col="Rich factor",
                gene_list_col="Proteins",
            ))
    return records


def load_joint_enrichment():
    """读取联合分析 GO/KEGG 富集表。

    参数:
        无。

    返回:
        标准化富集表列表。
    """
    records = []
    for compare in COMPARES:
        go_path = find_one(lambda p: p.name == f"{compare}__total.go_enrich.xls" and "3_Function_Analysis" in str(p)
                           and "2_Enrich" in str(p) and "1_GO_Enrich" in str(p))
        kegg_path = find_one(lambda p: p.name == f"{compare}__total.kegg_enrich.xls" and "3_Function_Analysis" in str(p)
                             and "2_Enrich" in str(p) and "4_KEGG_Enrich" in str(p))
        if go_path:
            df = read_tsv_or_excel(go_path)
            records.append(normalize_enrichment(
                df, "Joint GO enrichment", "Joint DEG-DEP", "GO", compare,
                term_col="go_term", id_col="go_id", p_col="pvalue", padj_col="padjust",
                count_col="pr_num", ratio_col=None, rich_col="enrich_factor_gene", gene_list_col="pr_list",
            ))
        if kegg_path:
            df = read_tsv_or_excel(kegg_path)
            records.append(normalize_enrichment(
                df, "Joint KEGG enrichment", "Joint DEG-DEP", "KEGG", compare,
                term_col="pathway_term", id_col="pathway_id", p_col="pvalue", padj_col="padjust",
                count_col="pr_num", ratio_col=None, rich_col="enrich_factor_gene", gene_list_col="pr_list",
            ))
    return records


def load_nine_quadrant_concordant_hits(theme_table):
    """从九象限同向变化分子中统计主题富集表命中。

    参数:
        theme_table: 已标准化富集/注释表。

    返回:
        命中统计 DataFrame。
    """
    matched_path = Path("output") / "nine_quadrant_plots" / "nine_quadrant_matched_gene_protein.csv"
    if not matched_path.exists():
        return pd.DataFrame()
    matched = pd.read_csv(matched_path)
    concordant = matched[matched["label_type_2"].isin([
        "Genes(up)_Proteins(up)", "Genes(down)_Proteins(down)"
    ])].copy()
    rows = []
    for _, term_row in theme_table[theme_table["theme"].astype(bool)].iterrows():
        molecules = set(re.split(r"[;|,\s]+", str(term_row["molecules"])))
        molecules.discard("")
        subset = concordant[concordant["compare"].map(COMPARE_TITLE) == term_row["compare_label"]]
        if subset.empty or not molecules:
            continue
        hit_mask = subset["gene"].astype(str).isin(molecules) | subset["protein"].astype(str).isin(molecules) | subset["seq_id"].astype(str).isin(molecules)
        hit_count = int(hit_mask.sum())
        if hit_count > 0:
            rows.append({
                "compare_label": term_row["compare_label"],
                "database": term_row["database"],
                "term": term_row["term"],
                "term_id": term_row["term_id"],
                "theme": term_row["theme"],
                "concordant_hit_count": hit_count,
                "hit_pairs": ";".join(subset.loc[hit_mask, "seq_id"].astype(str).tolist()),
                "note": "Annotation overlap in concordant nine-quadrant pairs; not an independent enrichment test.",
            })
    return pd.DataFrame(rows)


def filter_theme_terms(enrichment):
    """筛选主题相关富集条目。

    参数:
        enrichment: 标准化富集表。

    返回:
        主题条目 DataFrame。
    """
    theme = enrichment[enrichment["theme"].astype(str).str.len() > 0].copy()
    valid_p = theme["neg_log10_p"].notna()
    valid_count = theme["count"].fillna(0) > 0
    theme = theme[valid_p & valid_count].copy()
    theme = theme.sort_values(["database", "omics", "compare_label", "padjust", "pvalue", "term"])
    return theme


def format_panel_label(panel):
    """将气泡图横轴标签改为两行短标签。
    参数:
        panel: 原始 panel 标签。
    返回:
        适合作图显示的短标签。
    """
    parts = [part.strip() for part in str(panel).split("|")]
    if len(parts) >= 2:
        omics_label = "Joint" if parts[0] == "Joint DEG-DEP" else parts[0]
        compare_label = parts[1].replace(" vs ", "/")
        return f"{omics_label}\n{compare_label}"
    return str(panel)


def wrap_term_label(term, width=42):
    """将较长的 GO/KEGG 术语换行，避免左侧标签被压缩。
    参数:
        term: 术语名称。
        width: 每行最大字符数。
    返回:
        换行后的术语标签。
    """
    return "\n".join(textwrap.wrap(str(term), width=width, break_long_words=False, max_lines=2))


def plot_theme_bubble(theme, database, output_stem):
    """绘制主题富集气泡图。

    参数:
        theme: 主题富集表。
        database: GO 或 KEGG。
        output_stem: 输出文件名前缀。

    返回:
        无。
    """
    plot_df = theme[theme["database"] == database].copy()
    if plot_df.empty:
        return
    plot_df["panel"] = plot_df["omics"] + " | " + plot_df["compare_label"]
    # 每个 panel 最多保留显著性最强的 8 个主题条目，避免图过密。
    plot_df = plot_df.sort_values(["panel", "neg_log10_p"], ascending=[True, False]).groupby("panel").head(8)
    plot_df["term_short"] = plot_df["term"].str.replace(" signaling pathway", "", regex=False)
    plot_df["term_short"] = plot_df["term_short"].str.replace("Regulation of ", "Reg. ", regex=False)
    terms = list(dict.fromkeys(plot_df.sort_values("term_short")["term_short"].tolist()))
    panels = list(dict.fromkeys(plot_df["panel"].tolist()))
    x_map = {panel: i for i, panel in enumerate(panels)}
    y_map = {term: i for i, term in enumerate(terms)}
    fig_width = max(11.5, 0.86 * len(panels) + 5.0)
    term_height = 0.62 if database == "GO" else 0.44
    fig_height = max(5.8, term_height * len(terms) + 2.2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    x = plot_df["panel"].map(x_map)
    y = plot_df["term_short"].map(y_map)
    size = plot_df["count"].fillna(1).clip(lower=1) * 18
    scatter = ax.scatter(
        x, y, s=size, c=plot_df["neg_log10_p"], cmap="viridis",
        edgecolor="#555555", linewidth=0.35, alpha=0.88
    )
    ax.set_xticks(range(len(panels)))
    ax.set_xticklabels([format_panel_label(panel) for panel in panels], rotation=0, ha="center", fontsize=8)
    ax.set_yticks(range(len(terms)))
    y_fontsize = 7 if database == "GO" else 8
    ax.set_yticklabels([wrap_term_label(term, width=60) for term in terms], fontsize=y_fontsize)
    ax.set_title(f"Theme-focused {database} enrichment", fontsize=11, pad=8)
    ax.grid(True, color="#E6E6E6", linewidth=0.4)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.subplots_adjust(left=0.34, right=0.86, bottom=0.18, top=0.92)
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.035)
    cbar.set_label("-log10(P adjusted/P)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.savefig(OUTPUT_DIR / f"{output_stem}.png", dpi=300)
    fig.savefig(OUTPUT_DIR / f"{output_stem}.tiff", dpi=300)
    plt.close(fig)


def plot_theme_heatmap(theme):
    """绘制主题类别命中数量热图。

    参数:
        theme: 主题富集表。

    返回:
        无。
    """
    rows = []
    for _, row in theme.iterrows():
        for theme_name in str(row["theme"]).split("; "):
            if theme_name:
                rows.append({
                    "theme": theme_name,
                    "panel": f"{row['omics']} | {row['compare_label']} | {row['database']}",
                    "count": 1,
                })
    if not rows:
        return
    count_df = pd.DataFrame(rows).groupby(["theme", "panel"], as_index=False)["count"].sum()
    matrix = count_df.pivot(index="theme", columns="panel", values="count").fillna(0)
    fig_width = max(8, 0.55 * matrix.shape[1] + 2.5)
    fig_height = max(3.8, 0.45 * matrix.shape[0] + 1.2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(matrix.values, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=55, ha="right", fontsize=8)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=9)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = int(matrix.iloc[i, j])
            if value:
                ax.text(j, i, str(value), ha="center", va="center", fontsize=8, color="white" if value > matrix.values.max() / 2 else "black")
    ax.set_title("Theme term counts across omics layers", fontsize=11, pad=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Number of theme terms", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "theme_term_count_heatmap.png", dpi=300)
    fig.savefig(OUTPUT_DIR / "theme_term_count_heatmap.tiff", dpi=300)
    plt.close(fig)


def write_legend():
    """写出图例说明。

    参数:
        无。

    返回:
        无。
    """
    text = """# Figure legend

Figure X. Theme-focused GO and KEGG enrichment reanalysis of transcriptomic, proteomic, and joint transcriptome-proteome results. Terms related to ECM/cell adhesion, cytoskeleton, MAPK/PI3K/Wnt/mTOR signaling, calcium signaling, muscle contraction, and bone/ossification were selected from the original enrichment tables. Bubble size represents the number of genes, proteins, or matched gene-protein pairs assigned to each term, and color indicates -log10 adjusted P value when available or -log10 P value otherwise. DEG and DEP enrichment panels are based on the original transcriptome and proteome enrichment outputs, whereas joint panels summarize terms from the original transcriptome-proteome association enrichment results.
"""
    (OUTPUT_DIR / "theme_enrichment_figure_legend.md").write_text(text, encoding="utf-8")


def main():
    """执行主题功能富集重分析。

    参数:
        无。

    返回:
        无。
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_tables = []
    for loader in [load_deg_enrichment, load_dep_enrichment, load_joint_enrichment]:
        all_tables.extend(loader())
    enrichment = pd.concat(all_tables, ignore_index=True)
    theme = filter_theme_terms(enrichment)
    concordant_hits = load_nine_quadrant_concordant_hits(theme)
    plot_theme_bubble(theme, "KEGG", "theme_KEGG_enrichment_bubble")
    plot_theme_bubble(theme, "GO", "theme_GO_enrichment_bubble")
    plot_theme_heatmap(theme)
    write_legend()
    with pd.ExcelWriter(OUTPUT_DIR / "theme_enrichment_reanalysis.xlsx") as writer:
        enrichment.to_excel(writer, sheet_name="all_enrichment", index=False)
        theme.to_excel(writer, sheet_name="theme_terms", index=False)
        if not concordant_hits.empty:
            concordant_hits.to_excel(writer, sheet_name="concordant_hits", index=False)
    print("Theme terms by omics/database")
    print(theme.groupby(["omics", "database"]).size().reset_index(name="n").to_string(index=False))
    if not concordant_hits.empty:
        print("\nConcordant nine-quadrant annotation hits")
        print(concordant_hits.head(20).to_string(index=False))
    print(f"\n输出目录: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
