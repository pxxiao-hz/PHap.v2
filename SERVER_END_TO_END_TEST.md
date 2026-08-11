# PHap.v2 服务器从头测试流程

本文用于测试 GitHub 分支 `codex/bootstrap-minimal-pr`。所有形如
`__NAME__` 的字符串都是占位符，运行前必须替换；不要原样提交计算任务。

测试覆盖 PHap 当前可执行主流程：

```text
原始 HiFi/ONT/Hi-C
  -> hifiasm 初始组装
  -> p_ctg/p_utg FASTA
  -> dosage
  -> mT2T
  -> Hi-C 比对与 links
  -> cluster/recluster/rescue
  -> reads phasing、各 haplotype 重组装与 HapHiC scaffolding
```

Juicebox 人工校正、TGS-GapCloser 和 T2T-Polish 属于下游人工/可选步骤，
不纳入本次自动验收。

## 1. 一次性填写占位符

先创建服务器测试目录，并把下面内容保存为 `config.sh` 后逐项替换：

```bash
export PHAP_REPO="__PHAP_REPOSITORY_DIRECTORY__"
export PHAP_RUN="__EMPTY_ABSOLUTE_TEST_DIRECTORY__"

# HapHiC Conda 环境根目录。PHap 仍使用自己的 .venv；这里只借用其中的外部命令。
export HAPHIC_ENV="__HAPHIC_CONDA_ENV_DIRECTORY__"

export HIFI_READS="__ABSOLUTE_HIFI_FASTQ_GZ__"
export ONT_READS="__ABSOLUTE_ONT_FASTQ_GZ__"
export HIC_R1="__ABSOLUTE_HIC_R1_FASTQ_GZ__"
export HIC_R2="__ABSOLUTE_HIC_R2_FASTQ_GZ__"

export PLOIDY="__PLOIDY_INTEGER__"
export CHR_NUM="__HAPLOID_CHROMOSOME_COUNT__"
export THREADS="__THREADS_PER_PROCESS__"
export PROCESSES="__MAX_PARALLEL_PROCESSES__"
export CPU_BUDGET="__TOTAL_ALLOCATED_CPUS__"
export SEED="__FIXED_RANDOM_SEED__"

# Hi-C 文库实际限制性内切酶识别序列，例如 GATC；不可照抄未知值。
export RE_SITE="__RECOGNITION_SEQUENCE__"

# 若 auto 无法唯一识别单倍体深度，将它替换为有证据支持的数字，例如 32.5。
export HAPLOID_DEPTH="auto"
```

每次登录或提交调度任务后执行：

```bash
set -euo pipefail
source __ABSOLUTE_PATH_TO_CONFIG_SH__
test "$PHAP_RUN" != "/"
test "$PHAP_RUN" != "$HOME"
test "$PHAP_REPO" != "$PHAP_RUN"
test -d "$HAPHIC_ENV/bin"
test $((THREADS * PROCESSES)) -le "$CPU_BUDGET"
```

`PHAP_RUN` 必须是独立的新目录。PHap 的 `02.cluster/` 和 `03.phase_reads/`
由当前工作目录决定，不要在源码仓库或已有生产结果目录中运行。

## 2. 获取指定分支并建立环境

服务器能访问 GitHub 时：

```bash
git clone --branch codex/bootstrap-minimal-pr --single-branch \
  https://github.com/pxxiao-hz/PHap.v2.git "$PHAP_REPO"
cd "$PHAP_REPO"
git rev-parse --abbrev-ref HEAD
git rev-parse HEAD
```

预期分支为 `codex/bootstrap-minimal-pr`。本流程编写时对应提交为
`c465b20f895b42e73e39702920c3c707a0c15c76`；若提交不同，应记录实际哈希，
并确认是否有意测试更新后的分支。

离线服务器可上传仓库根目录已有的 `PHap.v2-c465b20.bundle`，然后运行：

```bash
git clone __ABSOLUTE_BUNDLE_PATH__ "$PHAP_REPO"
cd "$PHAP_REPO"
git switch codex/bootstrap-minimal-pr
```

创建环境：

```bash
python3 --version
python3 -m venv "$PHAP_REPO/.venv"
source "$PHAP_REPO/.venv/bin/activate"
python -m pip install --upgrade pip
python -m pip install -e "$PHAP_REPO[dev]"
```

PHap 使用独立 `.venv`，但 `filter_bam`、`haphic` 和 `samblaster` 可以来自
HapHiC 的 Conda 环境。不要在 `.venv` 激活后再运行 `conda activate haphic`，否则
`python` 和 `phap` 也可能被切换到 Conda 环境。按以下方式只加入外部命令：

```bash
source "$PHAP_REPO/.venv/bin/activate"
export PATH="$VIRTUAL_ENV/bin:$HAPHIC_ENV/bin:$PATH"
hash -r

test "$(command -v python)" = "$VIRTUAL_ENV/bin/python"
test "$(command -v phap)" = "$VIRTUAL_ENV/bin/phap"
command -v filter_bam
command -v haphic
command -v samblaster
```

对于本文对应的服务器，可在 `config.sh` 中填写：

```bash
export HAPHIC_ENV="/home/pxxiao/tools/Anaconda3/envs/haphic"
```

后续每个新的 SSH 会话或调度脚本都应依次执行：

```bash
source __ABSOLUTE_PATH_TO_CONFIG_SH__
source "$PHAP_REPO/.venv/bin/activate"
export PATH="$VIRTUAL_ENV/bin:$HAPHIC_ENV/bin:$PATH"
hash -r
```

Python 必须为 3.9 或更高。Linux 离线安装时应准备 Linux 对应的 wheel，
不要上传 macOS 的二进制 wheel。

## 3. 软件与代码基线验收

```bash
source "$PHAP_REPO/.venv/bin/activate"

phap --version
phap --help
phap dosage --help
phap mt2t --help
phap cluster --help
phap phase_reads --help

cd "$PHAP_REPO"
ruff check phap_core PHap.py tests scripts/phap_cluster.py
mypy
pytest -q
python -m compileall -q PHap.py phap_core scripts utils shell
git diff --check
```

所有命令必须以 0 退出。测试数量以后可能变化，因此以“0 failed”为验收标准，
不要硬编码测试总数。

检查完整流程所需外部命令。下面的检查只汇总缺失项，不调用 `exit`，因此直接在
SSH 登录 shell 中粘贴也不会断开连接：

```bash
missing_tools=""
for tool in awk bash bwa filter_bam haphic hifiasm mash minimap2 \
  pandepth pigz Rscript samblaster samtools seqkit sort; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "MISSING required tool: $tool" >&2
    missing_tools="$missing_tools $tool"
  else
    command -v "$tool"
  fi
done

if [ -n "$missing_tools" ]; then
  echo "Required tools not ready:$missing_tools" >&2
else
  echo "All required external tools found"
fi

# 仅当后续 phase_reads 设置了非默认 --ont_length 或 --ont_quality 时才需要。
if command -v chopper >/dev/null 2>&1; then
  command -v chopper
else
  echo "OPTIONAL missing: chopper (not needed with --ont_length 1 --ont_quality 0)"
fi

python - <<'PY'
import numpy, pandas, pysam, scipy
print("Python scientific dependencies: OK")
PY
```

`chopper` 只在 phase_reads 启用 ONT 长度/质量过滤时使用。当前文档使用默认
`--ont_length 1 --ont_quality 0`，因此缺少它不会阻止本次流程。
记录所有实际版本；无法用 `--version` 的工具按其帮助文档查询：

```bash
{
  git -C "$PHAP_REPO" rev-parse HEAD
  python --version
  phap --version
  hifiasm --version
  minimap2 --version
  mash --version
  samtools --version | head -n 1
  bwa 2>&1 | head -n 3
  seqkit version
} > "$PHAP_RUN/tool_versions.txt" 2>&1
```

## 4. 输入预检

```bash
mkdir -p "$PHAP_RUN"
cd "$PHAP_RUN"

for file in "$HIFI_READS" "$ONT_READS" "$HIC_R1" "$HIC_R2"; do
  test -s "$file" || { echo "missing/empty input: $file" >&2; exit 1; }
done

gzip -t "$HIFI_READS"
gzip -t "$ONT_READS"
gzip -t "$HIC_R1"
gzip -t "$HIC_R2"

# 只做少量抽查，不解压整个数据集。head 会让上游 gzip 正常收到 SIGPIPE，
# 因此只在这四条预览命令周围暂时关闭 pipefail。
set +o pipefail
zcat "$HIFI_READS" | head -n 8
zcat "$ONT_READS" | head -n 8
zcat "$HIC_R1" | head -n 8
zcat "$HIC_R2" | head -n 8
set -o pipefail

sha256sum "$HIFI_READS" "$ONT_READS" "$HIC_R1" "$HIC_R2" \
  > input_reads.sha256
```

确认 Hi-C R1/R2 来自同一文库、read name 可配对，且四个文件不是同一路径。

## 5. 初始 hifiasm 组装

```bash
mkdir -p "$PHAP_RUN/00.initial_assembly"
cd "$PHAP_RUN/00.initial_assembly"

hifiasm \
  -o initial.asm \
  -t "$CPU_BUDGET" \
  --ul "$ONT_READS" \
  "$HIFI_READS" \
  > hifiasm.stdout.log \
  2> hifiasm.stderr.log
```

从 GFA 导出 FASTA，ID 保持为 GFA `S` 行第二列：

```bash
awk 'BEGIN{OFS=""} $1=="S"{print ">",$2,"\n",$3}' \
  initial.asm.bp.p_ctg.gfa > p_ctg.fa.tmp
mv p_ctg.fa.tmp p_ctg.fa

awk 'BEGIN{OFS=""} $1=="S"{print ">",$2,"\n",$3}' \
  initial.asm.bp.p_utg.gfa > p_utg.fa.tmp
mv p_utg.fa.tmp p_utg.fa

test -s p_ctg.fa
test -s p_utg.fa
samtools faidx p_ctg.fa
samtools faidx p_utg.fa
test -z "$(cut -f1 p_ctg.fa.fai | LC_ALL=C sort | uniq -d)"
test -z "$(cut -f1 p_utg.fa.fai | LC_ALL=C sort | uniq -d)"
```

后续统一使用：

```bash
export P_CTG="$PHAP_RUN/00.initial_assembly/p_ctg.fa"
export P_UTG="$PHAP_RUN/00.initial_assembly/p_utg.fa"
```

## 6. dosage：HiFi 深度与剂量判定

```bash
mkdir -p "$PHAP_RUN/01.dosage"

bash "$PHAP_REPO/shell/px_shell_dosage.sh" \
  -g "$P_UTG" \
  -i "$HIFI_READS" \
  -t "$CPU_BUDGET" \
  -o "$PHAP_RUN/01.dosage"

# PanDepth 默认留下压缩表；PHap dosage 读取明确的纯文本文件。
gzip -cd "$PHAP_RUN/01.dosage/aln.sort.clean.pandepth.win.stat.gz" \
  | awk '$0 !~ /RegionLength/' \
  > "$PHAP_RUN/01.dosage/aln.sort.clean.pandepth.win.stat.txt.tmp"
mv "$PHAP_RUN/01.dosage/aln.sort.clean.pandepth.win.stat.txt.tmp" \
  "$PHAP_RUN/01.dosage/aln.sort.clean.pandepth.win.stat.txt"

phap dosage \
  --input-file "$PHAP_RUN/01.dosage/aln.sort.clean.pandepth.win.stat.txt" \
  --pandepth \
  --ploidy "$PLOIDY" \
  --haploid-depth "$HAPLOID_DEPTH" \
  --output-dir "$PHAP_RUN/01.dosage"
```

验收：

```bash
test -s "$PHAP_RUN/01.dosage/contig_depth.txt"
test -s "$PHAP_RUN/01.dosage/dosage_windows.tsv"
python -m json.tool "$PHAP_RUN/01.dosage/dosage_model.json" >/dev/null
head "$PHAP_RUN/01.dosage/contig_depth.txt"
```

如果 `auto` 报 `weakly identified`，这是安全停止。检查深度分布后，把
`HAPLOID_DEPTH` 改成证据支持的正数并重跑，不要猜测阈值。

```bash
export CONTIG_TYPE="$PHAP_RUN/01.dosage/contig_depth.txt"
```

## 7. 构建 mT2T

```bash
cd "$PHAP_RUN"

phap mt2t \
  --p_ctg "$P_CTG" \
  --output-directory "$PHAP_RUN/01.mT2T" \
  --threads "$THREADS" \
  --process "$PROCESSES" \
  --cpu-budget "$CPU_BUDGET" \
  --min-identity 0.8 \
  --min-shorter-coverage 0.2 \
  --min-longer-coverage 0.05 \
  > "$PHAP_RUN/01.mT2T.stdout.log" \
  2> "$PHAP_RUN/01.mT2T.stderr.log"

export MT2T="$PHAP_RUN/01.mT2T/04.remove.redundancy/mT2T.fa"
test -s "$MT2T"
samtools faidx "$MT2T"
```

审计每条输入 contig 恰有一条 routing：

```bash
awk '/^>/{sub(/^>/,""); split($0,a,/[^[:graph:]]/); print a[1]}' "$P_CTG" \
  | LC_ALL=C sort > "$PHAP_RUN/mt2t.input.ids"
tail -n +2 "$PHAP_RUN/01.mT2T/04.remove.redundancy/mt2t_sequence_routing.tsv" \
  | cut -f1 | LC_ALL=C sort > "$PHAP_RUN/mt2t.routed.ids"
diff -u "$PHAP_RUN/mt2t.input.ids" "$PHAP_RUN/mt2t.routed.ids"
test -z "$(uniq -d "$PHAP_RUN/mt2t.routed.ids")"
```

同时人工查看：

```bash
head -n 30 "$PHAP_RUN/01.mT2T/04.remove.redundancy/mt2t_join_audit.tsv"
head -n 30 "$PHAP_RUN/01.mT2T/04.remove.redundancy/mt2t_containment_audit.tsv"
```

branch、cycle、orientation conflict 不应被程序猜测性连接；应保留并进入审计。

## 8. Hi-C 比对并生成 full_links.pkl / paired_links.clm

该脚本在当前目录产生 `HiC.bam` 和 `HiC.filtered.bam`，因此在独立目录运行：

```bash
mkdir -p "$PHAP_RUN/01.hic_links"
cd "$PHAP_RUN/01.hic_links"

bash "$PHAP_REPO/shell/haphic_shell_data-prepare.sh" \
  "$P_UTG" "$HIC_R1" "$HIC_R2" --threads "$CPU_BUDGET" \
  > hic_mapping.stdout.log \
  2> hic_mapping.stderr.log

samtools quickcheck -v HiC.bam HiC.filtered.bam

python "$PHAP_REPO/shell/parse.hic.bam.py" \
  --bam HiC.filtered.bam \
  --fasta "$P_UTG" \
  --flank 500000 \
  --threads "$THREADS" \
  > parse_hic.stdout.log \
  2> parse_hic.stderr.log

test -s full_links.pkl
test -s paired_links.clm
export FULL_LINKS="$PHAP_RUN/01.hic_links/full_links.pkl"
export CLM="$PHAP_RUN/01.hic_links/paired_links.clm"
```

注意：`parse.hic.bam.py` 当前通过只读取 read1 来保证每个 Hi-C pair 计数一次。
保留原 BAM、过滤后 BAM 和日志，以便审计过滤策略。

## 9. cluster / recluster / rescue

建议先按 mT2T FASTA 的真实目标 ID 创建明确的 locus 列表，而不是仅依赖“最长
`CHR_NUM` 条序列”。文件必须恰好有 `CHR_NUM` 个不重复 ID，每行一个：

```bash
export LOCUS_TARGETS="__ABSOLUTE_LOCUS_TARGET_ID_FILE__"
test "$(awk 'NF && $1 !~ /^#/ {n++} END {print n+0}' "$LOCUS_TARGETS")" \
  -eq "$CHR_NUM"
```

若还没有可靠列表，可暂时不传 `--locus-targets`，但必须在报告中记录 PHap 使用了
最长 `CHR_NUM` 条 mT2T 序列作为目标。

```bash
mkdir -p "$PHAP_RUN/02.clustering_run"
cd "$PHAP_RUN/02.clustering_run"

phap cluster \
  --p_utg "$P_UTG" \
  --mT2T "$MT2T" \
  --contig_type "$CONTIG_TYPE" \
  --full_links "$FULL_LINKS" \
  --clm "$CLM" \
  --ploidy "$PLOIDY" \
  --chr_num "$CHR_NUM" \
  --threads "$THREADS" \
  --process "$PROCESSES" \
  --RE "$RE_SITE" \
  --hic-score-mode re_density \
  --min-locus-identity 0.8 \
  --min-locus-query-coverage 0.05 \
  --min-locus-score-margin 0.05 \
  --min-bin-support-bases 1 \
  --min-bin-coverage 0 \
  --min-hic-score 0 \
  --min-hic-margin 0 \
  --locus-targets "$LOCUS_TARGETS" \
  > cluster.run1.stdout.log \
  2> cluster.run1.stderr.log
```

第一次探索参数时，`--min-bin-*` 和 `--min-hic-*` 的宽松值仅用于验证流程可运行，
不是生产阈值。生产分析必须依据数据分布做敏感性测试。

如果不使用 locus 文件，删除命令中的整行
`--locus-targets "$LOCUS_TARGETS"`。

基本验收：

```bash
grep 'alignment cache' cluster.run1.stdout.log
test -s 02.cluster/01.putg_vs_mT2T/p_utg_vs_mT2T.paf
python -m json.tool \
  02.cluster/01.putg_vs_mT2T/p_utg_vs_mT2T.stage_manifest.json >/dev/null
test -s 02.cluster/01.putg_vs_mT2T/paf_alignment_audit.tsv
test -s 02.cluster/01.putg_vs_mT2T/unitig_locus_candidates.tsv
test -s 02.cluster/01.putg_vs_mT2T/locus_rescue_decisions.tsv
test -s 02.cluster/01.putg_vs_mT2T/allelic_bin_candidates.tsv
test -s 02.cluster/01.putg_vs_mT2T/allelic_bin_decisions.tsv
test -s 02.cluster/02.chr_seq/locus_sequence_manifest.tsv
test -s 02.cluster/02.chr_seq/locus_sequence_routing.tsv
test -s 02.cluster/05.rescue/recluster_merge_sources.tsv
test -s 02.cluster/05.rescue/group.reassignment.cluster.txt
```

确认每个 p_utg 恰有一条 locus routing，未定位序列保留到 `un_chr.fa`：

```bash
cut -f1 "$P_UTG.fai" | LC_ALL=C sort > input.unitig.ids
tail -n +2 02.cluster/02.chr_seq/locus_sequence_routing.tsv \
  | cut -f1 | LC_ALL=C sort > routed.unitig.ids
diff -u input.unitig.ids routed.unitig.ids
test -z "$(uniq -d routed.unitig.ids)"
test -f 02.cluster/02.chr_seq/un_chr.fa
```

确认 rescue merge 只读取当前 manifest 声明的 locus：

```bash
awk -F '\t' 'NR>1 && $1!="un_chr.fa" {sub(/\.putg\.fa$/, "", $1); print $1}' \
  02.cluster/02.chr_seq/locus_sequence_manifest.tsv \
  | LC_ALL=C sort > current.loci
tail -n +2 02.cluster/05.rescue/recluster_merge_sources.tsv \
  | cut -f1 | LC_ALL=C sort > merged.loci
diff -u current.loci merged.loci
```

检查 group 内和 group 间是否存在意外重复；collapsed unitig 按 dosage 出现在多个 group
是预期行为，所以不要仅用全局 `uniq -d` 判错，应与 `contig_depth.txt` 的 dosage 对照。

## 10. cluster 缓存和确定性复测

先记录机器可读输出的内容哈希。哈希清单本身放到 `02.cluster` 外，避免把自己纳入：

```bash
find 02.cluster -type f \
  \( -name '*.tsv' -o -name '*.txt' -o -name '*.json' -o -name '*.paf' \) \
  -print0 | LC_ALL=C sort -z | xargs -0 sha256sum \
  > checksums.run1
```

用第 9 节完全相同的 `phap cluster` 命令重跑，仅把日志改为 `run2`。然后：

```bash
grep 'alignment cache: hit' cluster.run2.stdout.log

find 02.cluster -type f \
  \( -name '*.tsv' -o -name '*.txt' -o -name '*.json' -o -name '*.paf' \) \
  -print0 | LC_ALL=C sort -z | xargs -0 sha256sum \
  > checksums.run2

diff -u checksums.run1 checksums.run2
```

相同输入、参数、PHap/外部工具版本和 seed 下，机器可读输出应字节稳定。

## 11. 为 phase_reads 准备三类 BAM

HiFi 和 ONT 都比对到同一个 `p_utg`；使用 primary alignment 并排序、索引：

```bash
mkdir -p "$PHAP_RUN/03.phase_input"
cd "$PHAP_RUN/03.phase_input"

minimap2 -ax map-hifi -t "$CPU_BUDGET" --secondary=no "$P_UTG" "$HIFI_READS" \
  | samtools sort -@ "$THREADS" -o hifi.p_utg.sorted.bam -
samtools index hifi.p_utg.sorted.bam

minimap2 -ax map-ont -t "$CPU_BUDGET" --secondary=no "$P_UTG" "$ONT_READS" \
  | samtools sort -@ "$THREADS" -o ont.p_utg.sorted.bam -
samtools index ont.p_utg.sorted.bam

samtools quickcheck -v hifi.p_utg.sorted.bam ont.p_utg.sorted.bam \
  "$PHAP_RUN/01.hic_links/HiC.filtered.bam"
```

Hi-C BAM 已在第 8 节相对于同一 `p_utg` 生成。PHap 会过滤 unmapped、secondary、
supplementary、duplicate、QC-fail 和低 MAPQ alignment，并按 pair 生成一次决策。

## 12. phase_reads、haplotype 重组装与 HapHiC

该阶段会同时运行多个 `seqkit`、`hifiasm` 和 HapHiC 作业，实际 CPU 峰值约受
`PROCESSES × THREADS` 影响。必须满足调度器分配，不要用 `nohup` 隐藏退出码。

```bash
mkdir -p "$PHAP_RUN/03.phasing_run"
cd "$PHAP_RUN/03.phasing_run"

phap phase_reads \
  --bam_hifi "$PHAP_RUN/03.phase_input/hifi.p_utg.sorted.bam" \
  --bam_hic "$PHAP_RUN/01.hic_links/HiC.filtered.bam" \
  --bam_ont "$PHAP_RUN/03.phase_input/ont.p_utg.sorted.bam" \
  --contig_type "$CONTIG_TYPE" \
  --group "$PHAP_RUN/02.clustering_run/02.cluster/05.rescue/group.reassignment.cluster.txt" \
  --ploidy "$PLOIDY" \
  --hifi "$HIFI_READS" \
  --ont "$ONT_READS" \
  --hic1 "$HIC_R1" \
  --hic2 "$HIC_R2" \
  --threads "$THREADS" \
  --process "$PROCESSES" \
  --seed "$SEED" \
  --min_mapq 1 \
  > phase_reads.stdout.log \
  2> phase_reads.stderr.log
```

若需要启用 ONT 过滤，再增加有证据支持的 `--ont_length` 和 `--ont_quality`；默认
`1/0` 不触发 chopper。

验收：

```bash
test -d 03.phase_reads
test -s 03.phase_reads/read_assignments.tsv
test -s 03.phase_reads/read_assignment_summary.tsv
test -s 03.phase_reads/unitig_candidates.tsv

for group in $(seq 1 "$PLOIDY"); do
  test -s "03.phase_reads/group${group}.HiFi.fq.gz"
  test -s "03.phase_reads/group${group}.ONT.fq.gz"
  test -s "03.phase_reads/group${group}.asm/group${group}.asm.bp.p_ctg.gfa"
  test -d "03.phase_reads/group${group}.asm/scaffolding/02_haphic/04.build"
done
```

如果实际 group 名不是 `group1...groupN`，按
`group.reassignment.cluster.txt` 第一列替换上述循环。重点检查：

- 每个 read entity 最多分配到一个最终 haplotype；
- `ambiguous` 和 `unassigned` 有明确数量且未被复制进多个组；
- collapsed-unitig reads 的分区在固定 seed 下可复现且各组互斥；
- 每组 hifiasm、HapHiC 和 Juicebox 脚本日志均无非零退出或截断输出。

## 13. 最终报告与科学指标

保存完整命令、输入哈希、版本和日志，并至少汇总：

```text
代码提交：
输入数据 SHA-256：
PLOIDY / CHR_NUM / RE_SITE：
线程、进程、CPU budget、seed：
dosage 模型与各状态 unitig 数：
mT2T 输入/保留/contained/joined/conflict contig 数：
locus assigned / ambiguous / unassigned 数：
各 group 的 unitig copy 数守恒检查：
各 group HiFi/ONT/Hi-C read 数：
跨 group read overlap：
cluster run2 cache 是否 hit：
run1/run2 哈希是否一致：
失败、跳过和人工调整项：
```

组装效果不能只报告 N50。至少比较 completeness、continuity、switch/misjoin
证据、haplotype duplication、原始 reads 支持和 Hi-C 支持。若调整科学阈值，保留
调整前后机器可读指标及完整命令。

## 14. 失败时的最小回传材料

不要上传原始 reads。回传以下小文件即可定位多数问题：

```text
tool_versions.txt
实际 git commit hash
失败命令及退出码
对应 stdout/stderr log 的末尾 100 行
dosage_model.json
p_utg_vs_mT2T.stage_manifest.json
locus_evidence_manifest.tsv
locus_rescue_decisions.tsv
allelic_bin_decisions.tsv
locus_sequence_manifest.tsv
recluster_merge_sources.tsv
read_assignment_summary.tsv（若已进入 phase_reads）
```

任何失败阶段都应停止，不要把部分输出当作有效缓存继续下游。
