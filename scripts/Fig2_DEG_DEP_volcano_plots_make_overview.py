from pathlib import Path
import shutil

import pandas as pd
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


ROOT_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT_DIR / "output"
FIG2_DIR = OUTPUT_DIR / "figures" / "Fig2_DEG_DEP_overview"
PANEL_DIR = FIG2_DIR / "single_panels"
DATA_DIR = FIG2_DIR / "source_data"

MRNA_PCA_FILE = ROOT_DIR / (
    "MJ20251014051-ZX-R-251014142-郑建波-真核有参转录组测序-9个样本-翘嘴鲌-20251027-FX2025102000108"
) / "翘嘴鲌(MJ20251014051) - 9" / "workflow_results" / "06_Express" / "ExpPCA" / "PCA.xls"

MRNA_VARIANCE_FILE = MRNA_PCA_FILE.with_name("Explained_variance_ratio.xls")

PROTEIN_PCA_FILE = ROOT_DIR / (
    "MJ20251014053-ZX-P-251019025-郑建波-高通量蛋白组学-9个样本-翘嘴鲌-尾部肌肉-20251117-FX2025102300143"
) / "翘嘴鲌-尾部肌肉-9" / "workflow_results" / "02_SampleComp" / "SamPca" / "pca_sites.xls"

PROTEIN_VARIANCE_FILE = PROTEIN_PCA_FILE.with_name("pca_importance.xls")


def infer_group(sample_name: str) -> str:
    """根据样本名推断分组。

    参数:
        sample_name: 样本名称。

    返回:
        样本分组名称。
    """
    if sample_name.startswith("WT"):
        return "Ctrl"
    if sample_name.startswith("He"):
        return "He"
    if sample_name.startswith("Ho"):
        return "Ho"
    return "Other"


def read_pca_table(path: Path, sample_col: str) -> pd.DataFrame:
    """读取 PCA 坐标表并补充分组列。

    参数:
        path: PCA 坐标表路径。
        sample_col: 样本列名称。

    返回:
        包含 sample、PC1、PC2 和 group 的数据表。
    """
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={sample_col: "sample"})
    df["group"] = df["sample"].astype(str).map(infer_group)
    return df


def read_variance_table(path: Path) -> dict[str, float]:
    """读取 PCA 方差解释率。

    参数:
        path: 方差解释率表路径。

    返回:
        PC 名称到方差解释率的映射。
    """
    df = pd.read_csv(path, sep="\t")
    pc_col = df.columns[0]
    value_col = df.columns[1]
    return dict(zip(df[pc_col].astype(str), df[value_col].astype(float)))


def draw_pca(
    df: pd.DataFrame,
    variance: dict[str, float],
    title: str,
    out_stem: Path,
    label_offsets: dict[str, tuple[int, int]] | None = None,
) -> None:
    """绘制统一风格 PCA 散点图。

    参数:
        df: PCA 坐标与分组数据。
        variance: 方差解释率映射。
        title: 子图标题。
        out_stem: 输出文件路径，不含扩展名。
        label_offsets: 样本标签偏移设置。

    返回:
        无返回值，保存 PNG 和 TIFF。
    """
    label_offsets = label_offsets or {}
    colors = {"Ctrl": "#A1A9D0", "He": "#F0988C", "Ho": "#B883D4", "Other": "#9E9E9E"}
    markers = {"Ctrl": "o", "He": "s", "Ho": "^", "Other": "D"}
    fig, ax = plt.subplots(figsize=(4.6, 3.8))

    for group in ["Ctrl", "He", "Ho", "Other"]:
        sub_df = df[df["group"] == group]
        if sub_df.empty:
            continue
        ax.scatter(
            sub_df["PC1"],
            sub_df["PC2"],
            s=52,
            color=colors[group],
            marker=markers[group],
            edgecolor="#333333",
            linewidth=0.6,
            label=group,
            zorder=3,
        )
        for _, row in sub_df.iterrows():
            offset = label_offsets.get(row["sample"], (4, 4))
            ax.annotate(
                row["sample"],
                (row["PC1"], row["PC2"]),
                xytext=offset,
                textcoords="offset points",
                fontsize=7.5,
            )

    pc1_label = f"PC1 ({variance.get('PC1', 0) * 100:.1f}%)" if "PC1" in variance else "PC1"
    pc2_label = f"PC2 ({variance.get('PC2', 0) * 100:.1f}%)" if "PC2" in variance else "PC2"
    ax.set_xlabel(pc1_label, fontsize=9)
    ax.set_ylabel(pc2_label, fontsize=9)
    ax.set_title(title, fontsize=9.5, pad=6)
    ax.axhline(0, color="#CCCCCC", linewidth=0.8, zorder=1)
    ax.axvline(0, color="#CCCCCC", linewidth=0.8, zorder=1)
    ax.grid(True, color="#E6E6E6", linewidth=0.6, zorder=0)
    ax.legend(frameon=False, fontsize=7.5, loc="best")
    ax.tick_params(labelsize=7.5)
    x_range = df["PC1"].max() - df["PC1"].min()
    y_range = df["PC2"].max() - df["PC2"].min()
    ax.set_xlim(df["PC1"].min() - x_range * 0.12, df["PC1"].max() + x_range * 0.18)
    ax.set_ylim(df["PC2"].min() - y_range * 0.15, df["PC2"].max() + y_range * 0.15)
    fig.tight_layout()
    fig.savefig(out_stem.with_suffix(".png"), dpi=300)
    fig.savefig(out_stem.with_suffix(".tiff"), dpi=300)
    plt.close(fig)


def copy_if_exists(src: Path, dst: Path) -> None:
    """复制存在的文件。

    参数:
        src: 源文件路径。
        dst: 目标文件路径。

    返回:
        无返回值。
    """
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def collect_fig2_panels() -> list[Path]:
    """整理 Fig. 2 已有单图。

    参数:
        无。

    返回:
        已复制文件的目标路径列表。
    """
    copied_files = []
    volcano_dir = OUTPUT_DIR / "volcano_plots"
    source_files = [
        volcano_dir / "DEG_He_vs_Ctrl_volcano.png",
        volcano_dir / "DEG_He_vs_Ctrl_volcano.tiff",
        volcano_dir / "DEG_Ho_vs_Ctrl_volcano.png",
        volcano_dir / "DEG_Ho_vs_Ctrl_volcano.tiff",
        volcano_dir / "DEG_Ho_vs_He_volcano.png",
        volcano_dir / "DEG_Ho_vs_He_volcano.tiff",
        volcano_dir / "DEP_He_vs_Ctrl_volcano.png",
        volcano_dir / "DEP_He_vs_Ctrl_volcano.tiff",
        volcano_dir / "DEP_Ho_vs_Ctrl_volcano.png",
        volcano_dir / "DEP_Ho_vs_Ctrl_volcano.tiff",
        volcano_dir / "DEP_Ho_vs_He_volcano.png",
        volcano_dir / "DEP_Ho_vs_He_volcano.tiff",
        volcano_dir / "combined_DEG_DEP_volcano.png",
        volcano_dir / "combined_DEG_DEP_volcano.tiff",
        volcano_dir / "combined_DEG_DEP_volcano_legend.md",
        volcano_dir / "volcano_summary.xlsx",
    ]

    for src in source_files:
        if not src.exists():
            continue
        dst = PANEL_DIR / src.name
        copy_if_exists(src, dst)
        copied_files.append(dst)

    return copied_files


def write_legend(copied_files: list[Path]) -> None:
    """生成 Fig. 2 图例说明文件。

    参数:
        copied_files: 已整理到 Fig. 2 目录的文件列表。

    返回:
        无返回值。
    """
    legend_path = FIG2_DIR / "Fig2_legend.md"
    unique_files = sorted({path.name for path in copied_files})
    file_lines = "\n".join(f"- `single_panels/{name}`" for name in unique_files)
    legend_text = f"""# Fig. 2 DEG/DEP expression overview

**Figure legend:** Principal component analysis (PCA) and volcano plots summarize the transcriptomic and proteomic expression profiles among Ctrl, He, and Ho groups. The mRNA PCA was redrawn from `PCA.xls`, with PC1 and PC2 explaining 37.7% and 19.4% of the variance, respectively. The protein PCA was redrawn from `pca_sites.xls` and `pca_importance.xls`, with PC1 and PC2 explaining 33.3% and 16.3% of the variance, respectively. Volcano plots show differentially expressed genes (DEGs) and differentially expressed proteins (DEPs) across pairwise comparisons.

## Single panels

{file_lines}

## Source data

- `source_data/mRNA_PCA.xls`
- `source_data/mRNA_Explained_variance_ratio.xls`
- `source_data/protein_pca_sites.xls`
- `source_data/protein_pca_importance.xls`
"""
    legend_path.write_text(legend_text, encoding="utf-8")


def main() -> None:
    """生成并整理 Fig. 2 所需图件。

    参数:
        无。

    返回:
        无返回值。
    """
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    mrna_pca = read_pca_table(MRNA_PCA_FILE, "sample")
    mrna_variance = read_variance_table(MRNA_VARIANCE_FILE)
    protein_pca = read_pca_table(PROTEIN_PCA_FILE, "Sample_ID")
    protein_variance = read_variance_table(PROTEIN_VARIANCE_FILE)

    draw_pca(
        mrna_pca,
        mrna_variance,
        "mRNA PCA",
        PANEL_DIR / "Fig2A_mRNA_PCA",
        {"Ho2": (4, -12), "WT2": (-18, 4), "WT3": (4, 8), "Ho3": (4, 10)},
    )
    draw_pca(
        protein_pca,
        protein_variance,
        "Protein PCA",
        PANEL_DIR / "Fig2B_protein_PCA",
        {"WT2": (4, -12), "WT1": (-18, -18), "Ho1": (20, 8), "He3": (-34, -10)},
    )

    copied_files = collect_fig2_panels()

    copy_if_exists(MRNA_PCA_FILE, DATA_DIR / "mRNA_PCA.xls")
    copy_if_exists(MRNA_VARIANCE_FILE, DATA_DIR / "mRNA_Explained_variance_ratio.xls")
    copy_if_exists(PROTEIN_PCA_FILE, DATA_DIR / "protein_pca_sites.xls")
    copy_if_exists(PROTEIN_VARIANCE_FILE, DATA_DIR / "protein_pca_importance.xls")

    copied_files.extend([
        PANEL_DIR / "Fig2A_mRNA_PCA.png",
        PANEL_DIR / "Fig2A_mRNA_PCA.tiff",
        PANEL_DIR / "Fig2B_protein_PCA.png",
        PANEL_DIR / "Fig2B_protein_PCA.tiff",
    ])
    write_legend(copied_files)


if __name__ == "__main__":
    main()
