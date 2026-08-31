# PHap v2 开发测试记录

更新：2026-08-31
代码基线：`testing/v2-multispecies-20260823` 分支，提交
[`3cded6f`](https://github.com/pxxiao-hz/PHap.v2/commit/3cded6f)

本文单独记录 PHap v2 开发期间实际使用的输入、代码、运行目录、参数和结果。
它是开发工作站的可追溯索引，不是原始数据发布清单：原始 FASTQ、BAM、p_utg、
Hi-C links 和运行结果均未复制进 Git 仓库，也不应被修改。

开发工作目录记为：

```text
WORKSPACE=/home/pxxiao/test/PHap.cluster
```

路径以 `$WORKSPACE` 表示。每次 phase_reads 运行的精确输入签名、工具路径、
线程与 checkpoint 状态保存在相应 `run_manifest.json`；详细数值保存在各阶段
`*.summary.json/tsv`。下文的汇总不替代这些原始审计文件。

## 1. 数据集和共享输入

本次开发数据为四倍体马铃薯（12 条染色体、4 个 haplotype groups）。

| 用途 | 开发工作站文件 | 大小 | 说明 |
| --- | --- | ---: | --- |
| p_utg | `$WORKSPACE/hs.100k.asm.bp.p_utg.gfa.remove.plastid.contamination.fa` | 2,623,806,054 bp | allelic table、cluster 与 rescue 的组装输入 |
| mT2T | `$WORKSPACE/hs.mT2T.v2.fa` | 758,144,122 bp | 构建染色体局部 allelic table 的参考 |
| contig dosage/type | `$WORKSPACE/contig_depth.txt` | 170,502 B | dosage 与 collapsed unitig 类型 |
| Hi-C unitig links | `$WORKSPACE/full_links.pkl` | 22,277,344 B | cluster、recluster、rescue 的 Hi-C 证据 |
| HiFi→p_utg BAM | `$WORKSPACE/01.hifi.data/HiFi.p_utg.sort.bam` | 41,699,502,565 B | HiFi reads 分配 |
| ONT→p_utg BAM | `$WORKSPACE/01.hifi.data/ONT.p_utg.sort.bam` | 44,532,834,625 B | ONT reads 分配 |
| Hi-C→p_utg BAM | `$WORKSPACE/01.hifi.data/HiC.p_utg.sort.bam` | 242,111,604,113 B | fast-v1 Hi-C pairs 分配 |

原始 HiFi、ONT、Hi-C FASTQ 均保持在开发工作站的只读数据目录。其绝对路径、
文件大小和修改时间记录在以下 manifest，而不在本文重复公开：

- HiFi/ONT：`03.phase_reads.v2.nomapq_all_hifi_ont_assemble_20260828/run_manifest.json`
- Hi-C：`03.phase_reads.v2.hic_fast_v1_all_20260828/run_manifest.json`

用于 cluster 的原始 PAF、no-sequence GFA、每条染色体 FASTA 和对应 allelic
table 均由运行脚本显式记录，见第 3 节的两个 shell 脚本。

## 2. 使用的 v2 代码与测试

本次整理后的代码入口为：

| 功能 | 文件 |
| --- | --- |
| 命令行入口 | `PHap.py`、`scripts/phap_allelic_table.py`、`scripts/phap_cluster.py`、`scripts/phap_phase_reads.py` |
| allelic table | `utils/allelic_table_generate_v2.py`、`utils/find_collinear_chains.py` |
| initial cluster | `utils/cluster_allelic_unitigs_v2.py` |
| recluster / rescue | `utils/chr_uncluster_recluster.py`、`utils/unchr_recluster.py` |
| phase_reads | `utils/phase_reads_assignment.py`、`utils/phase_reads_assemble_anchor.py`、`utils/phase_reads_demux.awk` |
| 自动化测试 | `tests/test_allelic_table_v2.py`、`test_chromosome_extraction.py`、`test_cluster_v2.py`、`test_recluster_v2.py`、`test_rescue_v2.py`、`test_phase_reads_v2.py` |

提交 `3cded6f` 上执行：

```shell
python -m unittest discover -s tests -v
```

结果为 **58 tests passed**。测试覆盖 allelic table、染色体提取、cluster、
recluster、rescue、phase_reads 的参数校验、checkpoint/input-signature、
fast extraction、corrected fast-v1 Hi-C 和 scaffold-only 小规模接口。

## 3. Cluster 开发验证

正式测试目录：

```text
$WORKSPACE/02.cluster.v2.directevidence.validation_20260824
```

可复核的实际运行脚本：

- `run_all_cluster.sh`：12 条染色体的 initial cluster；四条染色体并行。
- `run_remaining_04_05.sh`：12 条染色体 recluster 后执行一次全局 rescue。

核心配置为 ploidy 4、`--hic-link-normalization dosage`、最小 adjusted links 5、
group/chromosome margin 0.10、低置信度 `defer`，并保留 allelic block 联合分配和
长 unitig（>=5 Mb）保护。算法不使用亲本信息；yak 仅在流程后评估。

结果：

- `03.cluster`、`04.recluster` 均完成 12 条染色体；`05.rescue` 输出 48 个
  `chr*_group*.reassignment.fa`。
- 最终分配 2,895 个 unitig、4,053 个 dosage memberships、2,493,576,229 bp；
  1,669 个 unitig（87,160,532 bp）保持显式未分配。
- `rescue_summary.json` 的 dosage、长度、partition、违反约束及已固定 assignment
  改动均为 0。

主要结果与审计文件：

```text
$WORKSPACE/02.cluster.v2.directevidence.validation_20260824/03.cluster/
$WORKSPACE/02.cluster.v2.directevidence.validation_20260824/04.recluster/
$WORKSPACE/02.cluster.v2.directevidence.validation_20260824/05.rescue/
```

其中最终 `group.reassignment.cluster.txt` 是后续 phase_reads 的 group 输入。
更早的 allelic-table/cluster 验证细节见
[VALIDATION_CURRENT_DATA.md](VALIDATION_CURRENT_DATA.md)。

## 4. HiFi + ONT reads 分配、提取和组装

运行目录：

```text
$WORKSPACE/03.phase_reads.v2.nomapq_all_hifi_ont_assemble_20260828
```

输入为第 1 节的 HiFi/ONT BAM、原始 HiFi/ONT FASTQ、`contig_depth.txt` 和第 3 节
最终 group 文件。该运行不按 MAPQ 过滤或加权 reads；HiFi/ONT 仍使用长度和 identity
阈值，collapsed reads 以完整 unitig 为单位，只在其已有目标 groups 中互斥均衡。

主要运行参数：`--data-types hifi ont`、`--collapsed-policy balanced`、
`--extract-backend seqkit`、`--extract-threads 10`、`--jobs 4`、
`--threads-per-job 10`。HiFi/ONT FASTQ 各扫描一次，经 `seqkit -> gawk -> FIFO -> pigz`
直接分流；原始 FASTQ 不被写入或修改。

| 项目 | HiFi | ONT |
| --- | ---: | ---: |
| 原始 FASTQ records | 5,199,674 | 285,720 |
| 已请求/写入 reads | 4,869,083 / 4,869,083 | 260,590 / 260,590 |
| 缺失已分配 reads | 0 | 0 |
| 输出 groups | 48 | 48 |

48 个 hifiasm primary-contig GFA 均已完成，并具有 per-group checkpoint；组装工具为
hifiasm 0.19.8，运行资源为 4 jobs × 10 threads。48 个组装结果随后完成 yak 后验
评估；详细表在 `04.yak/yak_summary.tsv`。本次记录不将 yak 用于任何方法选择。

关键输出：

```text
$WORKSPACE/03.phase_reads.v2.nomapq_all_hifi_ont_assemble_20260828/01.assignments/
$WORKSPACE/03.phase_reads.v2.nomapq_all_hifi_ont_assemble_20260828/02.reads/
$WORKSPACE/03.phase_reads.v2.nomapq_all_hifi_ont_assemble_20260828/03.assembly/
$WORKSPACE/03.phase_reads.v2.nomapq_all_hifi_ont_assemble_20260828/04.yak/yak_summary.tsv
```

## 5. corrected fast-v1 Hi-C 分配与提取

运行目录：

```text
$WORKSPACE/03.phase_reads.v2.hic_fast_v1_all_20260828
```

配置为 `--data-types hic --hic-assignment-backend fast-v1`。该模式只扫描 primary
read1，并从同一 BAM 记录的 RNEXT 获取 mate unitig；MAPQ 不参与过滤或权重。两个
complete unitig 的已有 group 集合取交集；collapsed pair 只在该 unitig 已定义的
groups 中互斥分配。详细 SQLite Hi-C 模式保留为对照接口，但不是该大数据运行的后端。

| 指标 | 结果 |
| --- | ---: |
| 扫描 BAM alignment records | 2,288,859,376 |
| 有 group 证据的 pairs | 780,828,588 |
| 已分配 pairs | 439,688,014 |
| mate conflict 而 defer 的 pairs | 341,140,574 |
| collapsed 互斥均衡 pairs | 122,286,947 |
| FASTQ 写入 pairs | 439,688,014 |
| 缺失已分配 pairs | 0 |
| 输出 Hi-C group pairs | 48 |

提取使用 `seqkit grep -> gawk count -> pigz`，资源为 24 个并发 groups、每个 seqkit
1 thread。它为每个 group 扫描一次每端原始 FASTQ（48 次/端）；这是修正后 v1 快速
模式保留的 I/O 特征，运行成功不代表其 I/O 成本最优。

关键输出：

```text
$WORKSPACE/03.phase_reads.v2.hic_fast_v1_all_20260828/01.assignments/hic.fast_v1/summary.json
$WORKSPACE/03.phase_reads.v2.hic_fast_v1_all_20260828/02.reads/hic.summary.json
$WORKSPACE/03.phase_reads.v2.hic_fast_v1_all_20260828/02.reads/chr*_group*.Hi-C.[12].fq.gz
```

## 6. scaffold-only 小规模真实数据测试

运行目录：

```text
$WORKSPACE/03.phase_reads.v2.scaffold_test_20260831
```

该测试不重跑 reads 分配、提取或组装，而是连接第 4 节的 `03.assembly` 与第 5 节的
`02.reads`。测试 groups 为 `chr03_group1` 和 `chr09_group2`，参数为 2 jobs、
5 threads/job、1 HapHiC process；BWA、samtools、samblaster、filter_bam 与 HapHiC
均成功完成。

| Group | 输入 contig | 输入 bp | 输出 scaffold | 主 scaffold bp | 未接入主 scaffold 的 contig |
| --- | ---: | ---: | ---: | ---: | ---: |
| chr03_group1 | 20 | 63,891,697 | 2 | 63,851,385 | 1（42,112 bp） |
| chr09_group2 | 36 | 67,781,544 | 8 | 67,264,589 | 7（共519,755 bp） |

所有输入 contig 均保留；输出额外的 100-bp `N` 为 scaffold 连接间隙。两组均有
完整 checkpoint，且运行日志没有 error/warning/traceback。chr03 的排序一致性
指标低于 chr09，若后续扩大 scaffold 测试，应优先检查其 contact map 和连接方向。

关键输出：

```text
$WORKSPACE/03.phase_reads.v2.scaffold_test_20260831/04.scaffold/chr03_group1/02.haphic/04.build/
$WORKSPACE/03.phase_reads.v2.scaffold_test_20260831/04.scaffold/chr09_group2/02.haphic/04.build/
$WORKSPACE/03.phase_reads.v2.scaffold_test_20260831/run_manifest.json
```

## 7. 当前验证边界

已完成：12 条染色体的 cluster/recluster/rescue、48 groups 的 HiFi/ONT
assignment/extraction/assembly/yak、48对 Hi-C group FASTQ，以及2个 group 的
scaffold-only 真实数据测试。

本轮明确未做：

- 全部 48 groups 的 HapHiC scaffolding；
- Juicebox 输出及人工 contact-map 校正；
- 其他物种或非四倍体的端到端验证；
- variant-based read phasing。

因此该分支是可供测试的 v2 开发版本，而不是已经完成全部生产验证的通用版本。

## 8. 全流程复测脚本

当前马铃薯数据的完整复测脚本已经随代码保存：

1. [`examples/run_cluster_full_potato_v2.sh`](examples/run_cluster_full_potato_v2.sh)
   从 allelic table 运行至 `05.rescue`，并检查最终48个 group FASTA。
2. [`examples/run_phase_reads_full_potato_v2.sh`](examples/run_phase_reads_full_potato_v2.sh)
   使用前一步的最终 group 文件，运行 HiFi/ONT/Hi-C assignment、FASTQ extraction、
   hifiasm 和48个 HapHiC scaffold，并检查48个结果和 checkpoint。

两个脚本均拒绝写入非空输出目录。第二个脚本需要在包含 hifiasm、HapHiC 及其依赖
工具的环境中运行；其默认资源是4个并发任务、每任务10线程、Hi-C seqkit每任务4线程。
若更改资源或原始数据路径，应创建新的输出目录并将实际命令与 manifest 一并保存。
