from pathlib import Path
import csv
import math
from collections import defaultdict

import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

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


WORK_DIR = Path(__file__).resolve().parents[1]
FIG_DIR = WORK_DIR / "output" / "figures" / "Fig7_candidate_mechanism_network"
SINGLE_DIR = FIG_DIR / "single_panels"
SOURCE_DIR = FIG_DIR / "source_data"
INPUT_CSV = (
    WORK_DIR
    / "output"
    / "figures"
    / "Fig6_candidate_expression_heatmap"
    / "source_data"
    / "Fig6_candidate_annotation.csv"
)

THEME_LABELS = {
    "ECM": "ECM remodeling",
    "Cytoskeleton": "Cytoskeleton\norganization",
    "Calcium": "Calcium\nregulation",
    "Muscle": "Muscle\ndevelopment",
    "Ossification": "Ossification",
}

THEME_COLORS = {
    "ECM": "#6BAED6",
    "Cytoskeleton": "#9E9AC8",
    "Calcium": "#74C476",
    "Muscle": "#FDAE6B",
    "Ossification": "#E7969C",
}

REGULATION_COLORS = {
    "Up": "#D55E00",
    "Down": "#0072B2",
    "Mixed": "#7F7F7F",
}

REGULATION_LABELS = {
    "Up": "concordant up-regulation",
    "Down": "concordant down-regulation",
    "Mixed": "mixed or comparison-dependent",
}


def read_candidate_rows(input_csv):
    """读取 Fig.6 候选注释表。

    参数:
        input_csv: 候选注释 CSV 路径。

    返回:
        list[dict]: 候选记录列表。
    """
    with input_csv.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def theme_set_from_row(row):
    """根据主分类和功能关键词提取候选所属机制主题。

    参数:
        row: 单条候选记录。

    返回:
        set[str]: 机制主题集合。
    """
    text = " ".join([
        row.get("primary_theme", ""),
        row.get("module", ""),
        row.get("phenotype_relevance", ""),
    ]).lower()
    themes = set()
    if "ecm" in text or "adhesion" in text or "matrix" in text:
        themes.add("ECM")
    if "cytoskeleton" in text or "keratin" in text or "actin" in text:
        themes.add("Cytoskeleton")
    if "calcium" in text or "annexin" in text:
        themes.add("Calcium")
    if "muscle" in text or "myosin" in text or "contraction" in text:
        themes.add("Muscle")
    if "ossification" in text or "bone" in text or "mineral" in text:
        themes.add("Ossification")

    primary = row.get("primary_theme", "").strip()
    if primary in THEME_LABELS:
        themes.add(primary)
    return themes or {"Cytoskeleton"}


def get_regulation(rows):
    """汇总同一候选在不同比较中的表达方向。

    参数:
        rows: 同一候选的所有比较记录。

    返回:
        str: Up、Down 或 Mixed。
    """
    directions = []
    for row in rows:
        mrna = float(row.get("mRNA_log2FC") or 0)
        protein = float(row.get("protein_log2FC") or 0)
        if mrna > 0 and protein > 0:
            directions.append("Up")
        elif mrna < 0 and protein < 0:
            directions.append("Down")
        else:
            directions.append("Mixed")

    unique_directions = set(directions)
    if unique_directions == {"Up"}:
        return "Up"
    if unique_directions == {"Down"}:
        return "Down"
    return "Mixed"


def aggregate_candidates(rows):
    """把比较层面的候选合并为分子层面的网络节点。

    参数:
        rows: Fig.6 候选记录列表。

    返回:
        list[dict]: 聚合后的候选节点。
    """
    grouped = defaultdict(list)
    for row in rows:
        display_label = row["display_label"].strip()
        grouped[display_label].append(row)

    candidates = []
    for display_label, group_rows in grouped.items():
        gene_symbol, protein_id = [x.strip() for x in display_label.split("|", 1)]
        themes = set()
        comparisons = []
        for row in group_rows:
            themes.update(theme_set_from_row(row))
            comparisons.append(row["compare_label"])

        candidates.append({
            "display_label": display_label,
            "short_label": f"{gene_symbol}\n{protein_id}",
            "gene_symbol": gene_symbol,
            "protein_id": protein_id,
            "themes": sorted(themes),
            "primary_theme": group_rows[0].get("primary_theme", "Cytoskeleton"),
            "regulation": get_regulation(group_rows),
            "comparisons": "; ".join(sorted(set(comparisons))),
            "n_records": len(group_rows),
        })

    theme_order = {theme: i for i, theme in enumerate(THEME_LABELS)}
    return sorted(candidates, key=lambda x: (theme_order.get(x["primary_theme"], 99), x["gene_symbol"], x["protein_id"]))


def write_source_tables(candidates):
    """导出 Fig.7 网络节点和边信息。

    参数:
        candidates: 聚合后的候选节点。

    返回:
        None。
    """
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    node_path = SOURCE_DIR / "Fig7_candidate_network_nodes.csv"
    edge_path = SOURCE_DIR / "Fig7_candidate_network_edges.csv"

    with node_path.open("w", encoding="utf-8-sig", newline="") as f:
        fieldnames = ["node_id", "label", "node_type", "regulation", "themes", "comparisons", "n_records"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({
            "node_id": "mechanism",
            "label": "Candidate mechanism",
            "node_type": "central_theme",
            "regulation": "",
            "themes": "",
            "comparisons": "",
            "n_records": "",
        })
        for theme, label in THEME_LABELS.items():
            writer.writerow({
                "node_id": theme,
                "label": label.replace("\n", " "),
                "node_type": "mechanism_theme",
                "regulation": "",
                "themes": theme,
                "comparisons": "",
                "n_records": "",
            })
        for candidate in candidates:
            writer.writerow({
                "node_id": candidate["display_label"],
                "label": candidate["short_label"].replace("\n", " | "),
                "node_type": "overlap",
                "regulation": candidate["regulation"],
                "themes": "; ".join(candidate["themes"]),
                "comparisons": candidate["comparisons"],
                "n_records": candidate["n_records"],
            })

    with edge_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["source", "target", "edge_type"])
        writer.writeheader()
        for theme in THEME_LABELS:
            writer.writerow({"source": "mechanism", "target": theme, "edge_type": "theme"})
        for candidate in candidates:
            for theme in candidate["themes"]:
                writer.writerow({
                    "source": theme,
                    "target": candidate["display_label"],
                    "edge_type": "candidate_theme",
                })


def draw_round_box(ax, xy, text, width, height, color, fontsize=9):
    """绘制机制主题圆角节点。

    参数:
        ax: Matplotlib 坐标轴。
        xy: 节点中心坐标。
        text: 节点文本。
        width: 节点宽度。
        height: 节点高度。
        color: 节点颜色。
        fontsize: 文本字号。

    返回:
        None。
    """
    x, y = xy
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        facecolor=color,
        edgecolor="#4D4D4D",
        linewidth=0.9,
        alpha=0.85,
        zorder=4,
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, color="black", zorder=5)


def candidate_positions(candidates, theme_positions):
    """计算候选节点在主题外侧的布局。

    参数:
        candidates: 聚合后的候选节点。
        theme_positions: 主题节点坐标。

    返回:
        dict: 候选显示标签到坐标的映射。
    """
    by_theme = defaultdict(list)
    for candidate in candidates:
        primary = candidate["primary_theme"] if candidate["primary_theme"] in THEME_LABELS else candidate["themes"][0]
        by_theme[primary].append(candidate)

    positions = {}
    for theme, items in by_theme.items():
        tx, ty = theme_positions[theme]
        angle = math.atan2(ty, tx)
        outward_x = math.cos(angle)
        outward_y = math.sin(angle)
        tangent_x = -outward_y
        tangent_y = outward_x
        spacing = 0.95 if len(items) > 3 else 1.05
        start = -(len(items) - 1) / 2
        for index, candidate in enumerate(items):
            offset = (start + index) * spacing
            x = tx + outward_x * 1.55 + tangent_x * offset
            y = ty + outward_y * 1.15 + tangent_y * offset
            positions[candidate["display_label"]] = (x, y)
    return positions


def draw_network(output_png, output_tiff=None):
    """绘制候选机制网络图。

    参数:
        output_png: PNG 输出路径。
        output_tiff: TIFF 输出路径，可为空。

    返回:
        None。
    """
    rows = read_candidate_rows(INPUT_CSV)
    candidates = aggregate_candidates(rows)
    write_source_tables(candidates)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    SINGLE_DIR.mkdir(parents=True, exist_ok=True)

    theme_positions = {
        "ECM": (-3.6, 1.7),
        "Cytoskeleton": (0.0, 2.25),
        "Calcium": (3.6, 1.7),
        "Muscle": (-2.35, -2.0),
        "Ossification": (2.35, -2.0),
    }
    center = (0, 0)
    candidate_pos = candidate_positions(candidates, theme_positions)

    fig, ax = plt.subplots(figsize=(12.5, 8.2), dpi=300)
    ax.set_xlim(-6.6, 6.6)
    ax.set_ylim(-4.8, 4.15)
    ax.axis("off")

    draw_round_box(
        ax,
        center,
        "Reduced intermuscular\nbone phenotype",
        width=2.25,
        height=0.82,
        color="#F7F7F7",
        fontsize=9.5,
    )

    for theme, pos in theme_positions.items():
        ax.plot([center[0], pos[0]], [center[1], pos[1]], color="#BDBDBD", lw=1.1, zorder=1)
        draw_round_box(
            ax,
            pos,
            THEME_LABELS[theme],
            width=1.75,
            height=0.62,
            color=THEME_COLORS[theme],
            fontsize=8.6,
        )

    for candidate in candidates:
        x, y = candidate_pos[candidate["display_label"]]
        for theme in candidate["themes"]:
            tx, ty = theme_positions[theme]
            ax.plot([tx, x], [ty, y], color=THEME_COLORS[theme], lw=0.85, alpha=0.58, zorder=1)

        color = REGULATION_COLORS[candidate["regulation"]]
        ax.scatter(
            x,
            y,
            s=230,
            marker="D",
            facecolor=color,
            edgecolor="white",
            linewidth=1.0,
            alpha=0.92,
            zorder=3,
        )
        if x < 0 and y < 0:
            text_y = y - 0.31
            va = "top"
        else:
            text_y = y + 0.31
            va = "bottom"

        ax.text(
            x,
            text_y,
            candidate["short_label"],
            ha="center",
            va=va,
            fontsize=7.1,
            color="black",
            linespacing=0.92,
            zorder=6,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=0.7),
        )

    legend_items = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#FFFFFF",
               markeredgecolor="#4D4D4D", markersize=7, label="mRNA node"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#FFFFFF",
               markeredgecolor="#4D4D4D", markersize=7, label="protein node"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="#FFFFFF",
               markeredgecolor="#4D4D4D", markersize=7, label="mRNA-protein overlap"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=REGULATION_COLORS["Up"],
               markeredgecolor="white", markersize=8, label=REGULATION_LABELS["Up"]),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=REGULATION_COLORS["Down"],
               markeredgecolor="white", markersize=8, label=REGULATION_LABELS["Down"]),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=REGULATION_COLORS["Mixed"],
               markeredgecolor="white", markersize=8, label=REGULATION_LABELS["Mixed"]),
    ]
    ax.legend(
        handles=legend_items,
        loc="lower left",
        bbox_to_anchor=(0.015, 0.015),
        ncol=2,
        frameon=False,
        fontsize=8.2,
        handletextpad=0.6,
        columnspacing=1.2,
    )

    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.10)
    fig.savefig(output_png, dpi=300)
    if output_tiff is not None:
        fig.savefig(output_tiff, dpi=300)
    plt.close(fig)


def write_legend_md():
    """生成 Fig.7 英文图例说明。

    参数:
        无。

    返回:
        None。
    """
    legend_text = """# Fig. 7. Candidate mechanism network

Candidate mechanism network linking the final prioritized mRNA-protein candidates to five phenotype-relevant biological themes: ECM remodeling, cytoskeleton organization, calcium regulation, muscle development, and ossification. Diamond nodes indicate mRNA-protein overlap candidates. Node color represents the summarized direction of matched mRNA and protein changes across the available pairwise comparisons: orange, concordant up-regulation; blue, concordant down-regulation; gray, mixed or comparison-dependent regulation. Edges connect each candidate to the biological theme(s) supported by functional annotation and phenotype-relevance screening.
"""
    (FIG_DIR / "Fig7_legend.md").write_text(legend_text, encoding="utf-8")


def main():
    """生成 Fig.7 候选机制网络图及说明文件。

    参数:
        无。

    返回:
        None。
    """
    output_png = FIG_DIR / "Fig7_candidate_mechanism_network.png"
    output_tiff = FIG_DIR / "Fig7_candidate_mechanism_network.tiff"
    single_png = SINGLE_DIR / "Fig7_candidate_mechanism_network_single.png"
    single_tiff = SINGLE_DIR / "Fig7_candidate_mechanism_network_single.tiff"

    draw_network(output_png, output_tiff)
    draw_network(single_png, single_tiff)
    write_legend_md()


if __name__ == "__main__":
    main()
