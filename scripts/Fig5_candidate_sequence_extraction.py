from pathlib import Path
import csv
import re
from collections import defaultdict


WORK_DIR = Path(__file__).resolve().parents[1]
FIG7_DIR = WORK_DIR / "output" / "figures" / "Fig7_candidate_mechanism_network"
OUT_DIR = WORK_DIR / "output" / "figures" / "Fig7_candidate_mechanism_network" / "source_data" / "candidate_sequences"

FIG7_NODES = FIG7_DIR / "source_data" / "Fig7_candidate_network_nodes.csv"
FIG6_ANNOTATION = (
    WORK_DIR
    / "output"
    / "figures"
    / "Fig6_candidate_expression_heatmap"
    / "source_data"
    / "Fig6_candidate_annotation.csv"
)

TRANSCRIPT_FASTA = (
    WORK_DIR
    / "MJ20251014051-ZX-R-251014142-郑建波-真核有参转录组测序-9个样本-翘嘴鲌-20251027-FX2025102000108"
    / "翘嘴鲌(MJ20251014051) - 9"
    / "workflow_results"
    / "04_Assemble"
    / "Sequence"
    / "all_transcripts.fa"
)
CDS_FASTA = TRANSCRIPT_FASTA.with_name("all_cds.fa")
PEP_FASTA = TRANSCRIPT_FASTA.with_name("all_pep.fa")
PROTEOME_FASTA = (
    WORK_DIR
    / "MJ20251014053-ZX-P-251019025-郑建波-高通量蛋白组学-9个样本-翘嘴鲌-尾部肌肉-20251117-FX2025102300143"
    / "翘嘴鲌-尾部肌肉-9"
    / "workflow_results"
    / "ProteinFasta"
    / "protein.fa"
)


def read_csv_rows(path):
    """读取 CSV 文件。

    参数:
        path: CSV 文件路径。

    返回:
        list[dict]: CSV 记录列表。
    """
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def normalize_gene_id(gene_id):
    """统一 gene_id 格式。

    参数:
        gene_id: 原始 gene_id。

    返回:
        str: 去掉 gene- 前缀后的 gene_id。
    """
    return gene_id.strip().replace("gene-", "")


def get_fig7_candidate_ids():
    """获得 Fig.7 中实际展示的候选节点。

    参数:
        无。

    返回:
        set[str]: Fig.7 overlap 节点 display_label 集合。
    """
    nodes = read_csv_rows(FIG7_NODES)
    return {
        row["node_id"].strip()
        for row in nodes
        if row.get("node_type", "").strip() == "overlap"
    }


def build_candidate_table():
    """把 Fig.7 节点回连到 Fig.6 注释表，生成候选序列提取清单。

    参数:
        无。

    返回:
        list[dict]: 去重后的候选清单。
    """
    fig7_ids = get_fig7_candidate_ids()
    rows = read_csv_rows(FIG6_ANNOTATION)
    grouped = {}
    for row in rows:
        display_label = row["display_label"].strip()
        if display_label not in fig7_ids:
            continue
        key = (display_label, row["gene_id"].strip(), row["protein_id"].strip())
        if key not in grouped:
            gene_symbol = display_label.split("|", 1)[0].strip()
            grouped[key] = {
                "display_label": display_label,
                "gene_symbol": gene_symbol,
                "gene_id": row["gene_id"].strip(),
                "gene_id_normalized": normalize_gene_id(row["gene_id"]),
                "protein_id": row["protein_id"].strip(),
                "primary_theme": row.get("primary_theme", ""),
                "comparisons": set(),
                "regulation_types": set(),
            }
        grouped[key]["comparisons"].add(row.get("compare_label", ""))
        grouped[key]["regulation_types"].add(row.get("nine_quadrant_type", ""))

    candidates = []
    for item in grouped.values():
        item["comparisons"] = "; ".join(sorted(x for x in item["comparisons"] if x))
        item["regulation_types"] = "; ".join(sorted(x for x in item["regulation_types"] if x))
        candidates.append(item)
    return sorted(candidates, key=lambda x: (x["primary_theme"], x["gene_symbol"], x["protein_id"]))


def parse_fasta(path):
    """解析 FASTA 文件。

    参数:
        path: FASTA 文件路径。

    返回:
        list[tuple[str, str]]: header 和 sequence 组成的列表。
    """
    records = []
    header = None
    seq_parts = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def index_by_gene(records):
    """按 gene= 字段建立 FASTA 索引。

    参数:
        records: FASTA 记录。

    返回:
        dict[str, list[tuple[str, str]]]: gene_id 到 FASTA 记录的映射。
    """
    index = defaultdict(list)
    for header, seq in records:
        match = re.search(r"\bgene=([^\s]+)", header)
        if not match:
            continue
        gene_id = normalize_gene_id(match.group(1))
        index[gene_id].append((header, seq))
    return index


def index_by_first_token(records):
    """按 FASTA header 第一个字段建立索引。

    参数:
        records: FASTA 记录。

    返回:
        dict[str, tuple[str, str]]: ID 到 FASTA 记录的映射。
    """
    index = {}
    for header, seq in records:
        token = header.split()[0]
        index[token] = (header, seq)
    return index


def wrap_sequence(seq, width=70):
    """按固定宽度折行序列。

    参数:
        seq: 序列字符串。
        width: 每行字符数。

    返回:
        str: 折行后的序列。
    """
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_fasta(records, path):
    """写出 FASTA 文件。

    参数:
        records: FASTA 记录列表。
        path: 输出路径。

    返回:
        None。
    """
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for header, seq in records:
            f.write(f">{header}\n{wrap_sequence(seq)}\n")


def make_output_header(candidate, source_header, sequence_type):
    """生成便于 PCR 追踪的 FASTA header。

    参数:
        candidate: 候选清单记录。
        source_header: 原始 FASTA header。
        sequence_type: 序列类型。

    返回:
        str: 新 FASTA header。
    """
    safe_symbol = candidate["gene_symbol"].replace(" ", "_")
    return (
        f"{safe_symbol}|gene_id={candidate['gene_id']}|protein_id={candidate['protein_id']}"
        f"|theme={candidate['primary_theme']}|type={sequence_type}|source={source_header}"
    )


def main():
    """提取 Fig.7 候选基因/蛋白序列。

    参数:
        无。

    返回:
        None。
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    candidates = build_candidate_table()

    transcript_index = index_by_gene(parse_fasta(TRANSCRIPT_FASTA))
    cds_index = index_by_gene(parse_fasta(CDS_FASTA))
    pep_index = index_by_first_token(parse_fasta(PEP_FASTA))
    proteome_index = index_by_first_token(parse_fasta(PROTEOME_FASTA))

    transcript_records = []
    cds_records = []
    pep_records = []
    proteome_records = []
    summary_rows = []

    for candidate in candidates:
        gene_key = candidate["gene_id_normalized"]
        protein_key = candidate["protein_id"]
        transcripts = transcript_index.get(gene_key, [])
        cds_list = cds_index.get(gene_key, [])
        pep = pep_index.get(protein_key)
        proteome = proteome_index.get(protein_key)

        for idx, (header, seq) in enumerate(transcripts, start=1):
            out_header = make_output_header(candidate, header, f"transcript_isoform_{idx}")
            transcript_records.append((out_header, seq))
        for idx, (header, seq) in enumerate(cds_list, start=1):
            out_header = make_output_header(candidate, header, f"cds_isoform_{idx}")
            cds_records.append((out_header, seq))
        if pep:
            header, seq = pep
            pep_records.append((make_output_header(candidate, header, "transcriptome_pep"), seq))
        if proteome:
            header, seq = proteome
            proteome_records.append((make_output_header(candidate, header, "proteome_protein"), seq))

        row = dict(candidate)
        row.update({
            "transcript_records": len(transcripts),
            "cds_records": len(cds_list),
            "transcriptome_pep_found": "yes" if pep else "no",
            "proteome_protein_found": "yes" if proteome else "no",
        })
        summary_rows.append(row)

    write_fasta(transcript_records, OUT_DIR / "Fig7_candidate_transcript_sequences.fa")
    write_fasta(cds_records, OUT_DIR / "Fig7_candidate_cds_sequences.fa")
    write_fasta(pep_records, OUT_DIR / "Fig7_candidate_pep_from_transcriptome.fa")
    write_fasta(proteome_records, OUT_DIR / "Fig7_candidate_protein_from_proteome.fa")

    summary_path = OUT_DIR / "Fig7_candidate_sequence_summary.csv"
    with summary_path.open("w", encoding="utf-8-sig", newline="") as f:
        fieldnames = [
            "display_label", "gene_symbol", "gene_id", "gene_id_normalized",
            "protein_id", "primary_theme", "comparisons", "regulation_types",
            "transcript_records", "cds_records", "transcriptome_pep_found",
            "proteome_protein_found",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"candidates={len(candidates)}")
    print(f"transcript_records={len(transcript_records)}")
    print(f"cds_records={len(cds_records)}")
    print(f"transcriptome_pep_records={len(pep_records)}")
    print(f"proteome_protein_records={len(proteome_records)}")
    print(OUT_DIR)


if __name__ == "__main__":
    main()
