from pathlib import Path

import pandas as pd


OUTPUT_DIR = Path("output") / "figure_source_manifest"
OUTPUT_XLSX = OUTPUT_DIR / "figure_source_manifest.xlsx"
OUTPUT_MD = OUTPUT_DIR / "figure_source_manifest.md"


FIGURE_ROWS = [
    {
        "figure_group": "Volcano plots",
        "main_output_png": "output/volcano_plots/combined_DEG_DEP_volcano.png",
        "main_output_tiff": "output/volcano_plots/combined_DEG_DEP_volcano.tiff",
        "generated_by_script": "output/draw_volcano_plots.py",
        "processed_data_output": "output/volcano_plots/volcano_summary.xlsx",
        "primary_data_source": "Transcriptome DEG result tables and proteome DEP result tables in service-provider workflow_results directories.",
        "source_identification_rule": "Script searches original DEG/DEP files for He_vs_Ctrl, Ho_vs_Ctrl, and Ho_vs_He; proteome significance follows original significant/regulate fields.",
        "plotted_values": "log2 fold change and -log10 adjusted P value/P value; colors indicate up-, down-, or non-significant molecules.",
        "reliability_note": "Counts were cross-checked against the service-provider DEG/DEP classification logic.",
    },
    {
        "figure_group": "mRNA-protein correlation",
        "main_output_png": "output/correlation_and_venn/combined_mRNA_protein_correlation.png",
        "main_output_tiff": "output/correlation_and_venn/combined_mRNA_protein_correlation.tiff",
        "generated_by_script": "output/draw_correlation_and_venn.py",
        "processed_data_output": "output/correlation_and_venn/correlation_and_venn_summary.xlsx",
        "primary_data_source": "output/nine_quadrant_plots/nine_quadrant_matched_gene_protein.csv generated from matched transcriptome-proteome nine-quadrant table.",
        "source_identification_rule": "Uses matched gene-protein pairs with gene_log2fc and protein_log2fc for each comparison.",
        "plotted_values": "mRNA log2FC versus protein log2FC; Pearson and Spearman correlations reported.",
        "reliability_note": "Interpret as global matched-pair correlation; weak correlation does not imply absence of pathway-level convergence.",
    },
    {
        "figure_group": "DEG-DEP Venn",
        "main_output_png": "output/correlation_and_venn/combined_DEG_DEP_venn.png",
        "main_output_tiff": "output/correlation_and_venn/combined_DEG_DEP_venn.tiff",
        "generated_by_script": "output/draw_correlation_and_venn.py",
        "processed_data_output": "output/correlation_and_venn/correlation_and_venn_summary.xlsx",
        "primary_data_source": "Matched gene-protein table derived from original DEG/DEP and nine-quadrant association results.",
        "source_identification_rule": "DEG set: gene_significant == yes; DEP set: protein_significant == yes in each comparison.",
        "plotted_values": "DEG-only, DEP-only, and DEG-DEP overlap counts.",
        "reliability_note": "Overlap is small; use as descriptive integration evidence rather than pathway enrichment evidence.",
    },
    {
        "figure_group": "Nine-quadrant plots",
        "main_output_png": "output/nine_quadrant_plots/combined_nine_quadrant.png",
        "main_output_tiff": "output/nine_quadrant_plots/combined_nine_quadrant.tiff",
        "generated_by_script": "output/draw_nine_quadrant_plots.py",
        "processed_data_output": "output/nine_quadrant_plots/nine_quadrant_summary_and_candidates.xlsx; output/nine_quadrant_plots/nine_quadrant_matched_gene_protein.csv",
        "primary_data_source": "Original transcriptome-proteome association Nine_quadrant result.csv tables from service-provider interaction/association workflow.",
        "source_identification_rule": "Uses matched gene/protein log2FC, significance labels, and label_type_2 quadrant classification.",
        "plotted_values": "mRNA log2FC and protein log2FC; quadrants classify concordant, discordant, and single-layer changes.",
        "reliability_note": "Concordant quadrants are used for candidate screening; discordant quadrants indicate possible post-transcriptional regulation.",
    },
    {
        "figure_group": "Theme-focused KEGG/GO enrichment bubbles",
        "main_output_png": "output/theme_enrichment_reanalysis/theme_KEGG_enrichment_bubble.png; output/theme_enrichment_reanalysis/theme_GO_enrichment_bubble.png",
        "main_output_tiff": "output/theme_enrichment_reanalysis/theme_KEGG_enrichment_bubble.tiff; output/theme_enrichment_reanalysis/theme_GO_enrichment_bubble.tiff",
        "generated_by_script": "output/theme_enrichment_reanalysis.py",
        "processed_data_output": "output/theme_enrichment_reanalysis/theme_enrichment_reanalysis.xlsx",
        "primary_data_source": "Original DEG GO/KEGG enrichment tables, DEP GO/KEGG enrichment tables, and joint transcriptome-proteome enrichment tables.",
        "source_identification_rule": "Terms are selected by predefined theme keywords: ECM/cell adhesion, cytoskeleton, MAPK/PI3K/Wnt/mTOR, calcium signaling, muscle contraction, and bone/ossification.",
        "plotted_values": "Bubble size: number of DEGs or DEPs assigned to the term; color: -log10 adjusted P value or -log10 P value when adjusted P is unavailable.",
        "reliability_note": "Joint DEG-DEP entries with pr_num = 0 or missing P values were excluded; do not claim significant DEG-DEP intersection enrichment for those terms.",
    },
    {
        "figure_group": "Theme term count heatmap",
        "main_output_png": "output/theme_enrichment_reanalysis/theme_term_count_heatmap.png",
        "main_output_tiff": "output/theme_enrichment_reanalysis/theme_term_count_heatmap.tiff",
        "generated_by_script": "output/theme_enrichment_reanalysis.py",
        "processed_data_output": "output/theme_enrichment_reanalysis/theme_enrichment_reanalysis.xlsx",
        "primary_data_source": "Theme-filtered DEG/DEP GO and KEGG enrichment terms.",
        "source_identification_rule": "Counts terms per theme, omics layer, database, and comparison after validity filtering.",
        "plotted_values": "Number of theme-related enrichment terms.",
        "reliability_note": "This is a term-count summary, not a statistical enrichment test by itself.",
    },
    {
        "figure_group": "Candidate DEG-DEP heatmaps",
        "main_output_png": "output/candidate_heatmaps/combined_candidate_mRNA_protein_heatmap.png",
        "main_output_tiff": "output/candidate_heatmaps/combined_candidate_mRNA_protein_heatmap.tiff",
        "generated_by_script": "output/draw_candidate_heatmaps.py",
        "processed_data_output": "output/candidate_heatmaps/candidate_gene_protein_heatmap_data.xlsx",
        "primary_data_source": "output/nine_quadrant_plots/nine_quadrant_matched_gene_protein.csv",
        "source_identification_rule": "Selects gene-protein pairs where both gene_significant and protein_significant are yes.",
        "plotted_values": "Group mean mRNA expression and protein abundance after row-wise z-score scaling.",
        "reliability_note": "Shows expression pattern of strict DEG-DEP candidate pairs only.",
    },
    {
        "figure_group": "Core pathway candidate heatmap",
        "main_output_png": "output/core_pathway_candidates/core_pathway_mRNA_protein_heatmap.png",
        "main_output_tiff": "output/core_pathway_candidates/core_pathway_mRNA_protein_heatmap.tiff",
        "generated_by_script": "output/build_core_pathway_candidates.py",
        "processed_data_output": "output/core_pathway_candidates/core_pathway_candidate_table.xlsx",
        "primary_data_source": "Nine-quadrant matched table plus theme-focused enrichment annotation table.",
        "source_identification_rule": "Candidates are matched gene-protein pairs with at least one significant layer and annotation to core pathway terms.",
        "plotted_values": "Group mean mRNA expression and protein abundance after row-wise z-score scaling.",
        "reliability_note": "Complete candidate table is preserved; plotted subset prioritizes High/Medium candidates.",
    },
    {
        "figure_group": "Core pathway-candidate network",
        "main_output_png": "output/core_pathway_candidates/core_pathway_network.png",
        "main_output_tiff": "output/core_pathway_candidates/core_pathway_network.tiff",
        "generated_by_script": "output/build_core_pathway_candidates.py",
        "processed_data_output": "output/core_pathway_candidates/core_pathway_candidate_table.xlsx; output/core_pathway_candidates/annotated_priority_candidates.xlsx",
        "primary_data_source": "Core pathway candidate table and annotation-derived pathway membership.",
        "source_identification_rule": "Edges connect selected candidate gene/protein pairs to KEGG/GO theme pathways from annotation hits.",
        "plotted_values": "Network topology of candidate-pathway associations; node colors indicate concordant up/down, discordant, or other patterns.",
        "reliability_note": "Network is an explanatory visualization, not a network inference or causality test.",
    },
    {
        "figure_group": "Main representative candidate heatmap",
        "main_output_png": "output/main_candidate_heatmap/main_core_candidate_heatmap.png",
        "main_output_tiff": "output/main_candidate_heatmap/main_core_candidate_heatmap.tiff",
        "generated_by_script": "output/draw_main_candidate_heatmap.py",
        "processed_data_output": "output/main_candidate_heatmap/main_core_candidate_heatmap_data.xlsx",
        "primary_data_source": "output/core_pathway_candidates/candidate_annotation_from_existing_genome.xlsx and output/core_pathway_candidates/core_pathway_candidate_table.xlsx.",
        "source_identification_rule": "Selects representative High/Medium candidates with shared He/Ho regulation, concordant transcript-protein changes, and relevance to ECM, cytoskeleton, calcium signaling, muscle contraction, or ossification.",
        "plotted_values": "Group mean mRNA expression and protein abundance after row-wise z-score scaling; labels use resolved KO/ortholog/domain annotation when available.",
        "reliability_note": "This is the recommended main-text heatmap; full High/Medium candidates remain available as supplementary table/figure.",
    },
]


def path_exists_text(value):
    """检查分号分隔路径是否存在。
    参数:
        value: 路径字符串。
    返回:
        存在状态说明。
    """
    paths = [item.strip() for item in str(value).split(";") if item.strip()]
    statuses = []
    for path in paths:
        statuses.append(f"{path}: {'exists' if Path(path).exists() else 'missing'}")
    return " | ".join(statuses)


def main():
    """生成图件数据来源溯源表。
    参数:
        无。
    返回:
        无。
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.DataFrame(FIGURE_ROWS)
    manifest["png_status"] = manifest["main_output_png"].apply(path_exists_text)
    manifest["tiff_status"] = manifest["main_output_tiff"].apply(path_exists_text)
    manifest["script_status"] = manifest["generated_by_script"].apply(path_exists_text)
    manifest["processed_data_status"] = manifest["processed_data_output"].apply(path_exists_text)

    with pd.ExcelWriter(OUTPUT_XLSX) as writer:
        manifest.to_excel(writer, sheet_name="figure_source_manifest", index=False)

    lines = ["# 图件数据来源与可靠性溯源表", ""]
    for _, row in manifest.iterrows():
        lines.append(f"## {row['figure_group']}")
        lines.append("")
        lines.append(f"- 输出图：`{row['main_output_png']}`")
        lines.append(f"- 生成脚本：`{row['generated_by_script']}`")
        lines.append(f"- 处理后数据：`{row['processed_data_output']}`")
        lines.append(f"- 主要数据来源：{row['primary_data_source']}")
        lines.append(f"- 数据选择规则：{row['source_identification_rule']}")
        lines.append(f"- 图中数值含义：{row['plotted_values']}")
        lines.append(f"- 可靠性说明：{row['reliability_note']}")
        lines.append("")
    OUTPUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(OUTPUT_XLSX.resolve())
    print(OUTPUT_MD.resolve())


if __name__ == "__main__":
    main()
