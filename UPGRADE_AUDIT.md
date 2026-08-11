# PHap v2 Upgrade Audit

审查基线：commit `a43272d`。本轮进行了静态代码审查、Python 语法编译和
CLI smoke test；未使用生产数据，也未运行 hifiasm、HapHiC、minimap2 等外部工具。

## 实施状态

第一批最小 PR 已完成以下工程基线工作：

- 新增 `pyproject.toml`、可安装的 `phap` console entry point 和依赖声明；
- 新增 CLI/preflight/worker failure smoke tests、tiny FASTA/PAF/group/SAM fixtures
  和 GitHub Actions CI；
- Python 子脚本改为 `sys.executable -m ...` 调用，并统一传播外部命令失败；
- 移除 read assembly/scaffolding 路径中的个人绝对工具路径，修复 `pigp` 拼写；
- 修正文档中的仓库地址、重复 subcommand、Hi-C mate 文件和不可移植示例路径。

第二批本地 correctness 修复已完成：

- dosage 从固定 `base_value=28` 改为 window-level、ploidy-aware 的自动间距拟合，
  同时保留显式 `--haploid-depth` 覆盖；弱识别模型默认停止，不静默猜测；
- collapsed read 改为每种数据类型全局一次、确定性且互斥的分配，冲突和无支持记录
  进入显式审计输出；
- 新增 `phap dosage`、window/model 审计、unitig candidate/read assignment/group
  summary 审计及相应回归测试。

当前 dosage 分类策略按后续测试反馈作了显式过渡调整：低于 `0.5 × haploid_depth`
的 window 归入 `dosage_1`；unitig 的最终 dosage/type 仅使用全部 window 的
`average_depth`，并由同一个 dosage model 和 `--min-confidence` 阈值分类。
`dominant_class`、`dominant_fraction`、`class_counts` 和 `mixed_dosage` 只作为
window 组成的审计字段，不影响最终判定。原 `--min-unitig-support` 已移除，避免
暴露不再生效的比例参数。机器可读的合成 before/after 证据保存在
`tests/fixtures/dosage_policy_before_after.tsv`，验证命令为：

```shell
python -m pytest -q tests/test_dosage.py
```

第三批本地 clustering correctness 修复已完成：

- clustering/re-clustering 改用显式 `--ploidy` 和动态 group，dosage `d` 必须
  完整进入 `d` 个 group，否则不产生部分分配；
- 显式 legacy/external `low_coverage`、`ambiguous`、`high_copy` 和缺失 dosage
  不会在下游静默回退为 haplotig；
- raw Hi-C count 与 group restriction-site density 同时计算和审计，
  `--hic-score-mode` 决定实际排序；
- 无 Hi-C 或证据并列的 unitig 不再删除或猜测；有 mT2T chromosome/bin 的记录
  保留为 `locus_assigned_haplotype_unresolved`，并写出 ID、FASTA、reason；
- group 文件解析显式区分带 count 和不带 count 两个阶段，输出稳定排序。

第四批本地 PAF/locus correctness 修复已完成：

- 新增 lossless、tag-aware PAF parser；`tp` 按标签名读取，缺失、secondary、
  malformed tag 和越界半开区间均显式审计或报错；
- locus 候选从完整 PAF 构建，不再先丢弃 next-best target；正向/反向链同时检查
  query 与 oriented target 的共线性；
- query coverage 使用 interval union，identity 使用
  `sum(nmatch) / sum(alignment_block_length)`，并实际应用 identity、coverage 和
  best-minus-next score margin；
- mT2T 仅决定 locus routing，不直接分配 haplotype group；low-coverage 只有在
  可选独立 read support 达标时才标记为候选，仍不强制 dosage 1；
- allelic-bin target support 改用 clipped target interval union；首个 group seed
  也不再允许完全无 Hi-C 的非全倍性 unitig 绕过证据检查。

后续本地 correctness 修复还已完成：

- 旧 target-gap-only LIS 已被严格、可审计的 alignment chain 替换；
- mT2T reverse join 和 pairwise-first-win 已替换为冲突感知 oriented overlap graph；
- allelic-bin 候选从 contig-count `top_n` 截断改为 dosage-capacity 优化，
  目标为 `union_support_bases * dosage`；并列最优不按 ID 猜测，每个
  candidate/bin 都有独立审计记录。

当前仍不能宣称任意倍性生产数据端到端已验证：尚未运行真实全流程 benchmark；
严格 locus chain、mT2T overlap graph 和 dosage-capacity allelic selection 仍需以
相同生产输入比较 completeness、switch/misjoin、duplication、read/Hi-C support
和资源使用。

下面的问题清单仍保留为基线审计记录，实施状态以各条目内的
`状态` 为准。

## 1. 当前基线

- 仓库约 5,600 行，核心逻辑分散在 `scripts/`、`utils/`、`shell/` 的脚本中。
- 没有 `pyproject.toml`/requirements、自动化测试、CI、示例小数据或可验证的
  benchmark。
- `python3 -m compileall` 通过，`python3 PHap.py --version` 可运行。
- `mt2t --help` 和 `cluster --help` 在当前干净 Python 环境中因缺少 `pandas`
  导入失败。帮助信息不应依赖完整科学计算环境才能显示。
- 当前环境中仅检测到 Rscript；其余外部组装工具未安装，因此本轮不能进行
  端到端运行验证。

## 2. P0：先修复，否则结果可能错误或流程无法运行

### P0-1 collapsed unitig 的 reads 没有被互斥拆分

**状态：已在第二批本地修复。**

`utils/phase_reads_assemble_anchor.py:249-265` 对每个 group 独立调用
`split_reads()`，每次使用相同 seed。一个 diplotig/triplotig/tetraplotig 同时出现
在多个 group 时，各 group 会得到同一个 `1/dosage` 子集，而不是不同的互斥子集；
其余 reads 被丢弃。这会直接造成 haplotype 间污染和覆盖不足。

当前实现先验证 unitig dosage 与候选 group 数守恒，再以 read/read-pair 为单位全局
处理一次。兼容候选组使用交集解析，collapsed reads 使用带 seed 的稳定 BLAKE2
分配；各 group 互斥，ambiguous/unassigned 保留在审计表中。

### P0-2 read phasing 接收所有长读长比对

**状态：已在第二批本地修复。**

`utils/phase_reads_assemble_anchor.py:114-129` 未过滤 unmapped、secondary、
supplementary、duplicate、MAPQ 或最佳比对。一个多重比对 read 可进入多个 unitig，
随后进入多个 haplotype。Hi-C 路径也没有完整、可配置的质量策略。

当前策略仅使用 mapped、primary、non-supplementary、non-duplicate、non-QC-fail
且达到 `--min_mapq` 的比对；Hi-C mates 作为一个 pair entity 处理。跨 group
冲突单独报告，不允许通过 `set.update()` 静默复制。non-proper Hi-C alignment
有意保留，因为有效 Hi-C contact 不要求传统 paired-end 的 proper-pair 几何关系。

### P0-3 dosage 归一化计算后未使用

**状态：活动 clustering/re-clustering 路径已在第三批本地修复。**

`scripts/phap_cluster.py:138-150` 调用 `adjust_hic_signals()` 后丢弃返回值，后续仍
使用原始 links。当前 dosage 权重实际上没有生效。

修复方向：将 Hi-C scoring 提取为纯函数，用合成 links 验证 dosage、RE sites、
group size 归一化，并明确使用 raw 还是 normalized score。

### P0-4 无 Hi-C 信号的 unitig 被当作错误分型并删除

**状态：活动 clustering/re-clustering 路径已在第三批本地修复。**

`scripts/phap_cluster.py:645-653,681-685` 从四个 group 中直接过滤不在
`dic_contig_hic` 的 unitig。缺少 Hi-C links 不等于序列错误，会降低完整性并导致
无审计的数据丢失。

当前实现不再调用删除无 Hi-C unitig 的旧过滤路径。无信号、弱信号和分数边界并列
分别进入结构化 decision audit；有 mT2T locus 的序列保留为
`locus_assigned_haplotype_unresolved`，无 locus 的序列进入 `unassigned`，并输出
原始 ID 的 FASTA。旧函数将在后续独立的 behavior-preserving cleanup commit 删除。

### P0-5 mT2T reverse join 坐标转换错误，单 alignment 还可能崩溃

**状态：已在第六批本地修复。** reverse overlap 先保留原始区间，再按
`[length-old_end, length-old_start)` 一次性转换；单 alignment 直接作为合法链处理，
不再调用 `linregress()`。join 只接受两侧精确到端点的 overlap，避免把未比对 prefix
混入 trim；正向、反向、反向 canonical traversal、单 alignment 和非零 overhang
均有精确序列或保留断言。

`scripts/util.py:798-801` 先覆盖 `tstarts`，再用已覆盖值计算 `tends`，反向坐标不是
正确的 `[tlen-old_end, tlen-old_start)`。`calc_stats()` 在 LIS 只有一个点时调用
`linregress()`，也没有边界处理。

修复方向：先保留 old starts/ends 后成对转换；为正向、反向、单 alignment、
contained 和冲突 overlap 建立精确序列断言。

### P0-6 mT2T 目前不是冲突感知的多 contig 组装

**状态：已在第六批本地修复。** 新的 `phap_core/mt2t_overlap.py` 以 contig 物理端
构建 bidirected oriented overlap graph。只有无歧义 simple path 会被物化；branch、
cycle 和 orientation conflict 的整组 contig 原样保留，并写入 edge、join、routing、
ID mapping 和 manifest 审计。每个 source contig 恰好一条 routing 记录。
Containment 只有在 child interval-union coverage 为 100%、两端完整且 host 唯一时
才删除；partial/competing child 均原样保留。identity、shorter/longer coverage 都是
显式、方向无关的 CLI 参数。

`scripts/util.py:871-884` 按 PAF 字典顺序尝试独立 pairwise merge。每条 contig
只能碰巧加入第一个成功 pair，不能处理 A-B-C 链、分支、环、方向冲突或多个候选
overlap；结果还依赖遍历顺序。

修复方向：构建带方向和置信度的 overlap graph，先去 contained contig，再进行
冲突检测与 path extraction；有歧义时输出 graph/报告，不猜测连接。merge 后必须
以原始 reads、覆盖连续性和 Hi-C 支持复核。

### P0-7 PAF 覆盖率和 primary 标记解析不可靠

**状态：活动 cluster locus 路径已在第四批本地修复；mT2T overlap backend 已在
第六批迁移。**

多个模块把 alignment match length 直接累加，重叠 alignment 会重复计数，比例可
被放大；多个位置使用 `parts[16][-1]` 判断 primary/supplementary，但 PAF optional
tag 没有固定列号。

涉及：

- `scripts/util.py`
- `utils/find_longest_subsequence.py`
- `utils/allelic_table_generate.py`
- `utils/extract_chr_from_putg.py`

当前 cluster 流程从完整 PAF 生成 alignment、candidate、decision 和 routing
审计；optional tag 顺序、缺失 `tp`、secondary、overlap union、正反向 chain、
竞争 locus 和确定性输出均有回归测试。mT2T overlap backend 现在复用同一 PAF
parser 和 chain 规则；旧的 fixed-column、summed-coverage 和 pairwise-first-win
实现已从 `scripts/util.py` 删除。

mT2T 上游不再按文件存在性复用 split/Mash/minimap2 结果：每次在隔离 staging
目录重建，使用安全数字文件名、参数列表 subprocess 和原子 `merge.paf`；0/1 个
达到长度阈值的 contig 直接生成空 PAF 后进入完整 routing。upstream manifest 记录
输入 hash、参数、工具版本、精确命令和 PAF hash。

### P0-8 流程存在明确的可执行性故障

- Python 依赖未声明：至少涉及 numpy、pandas、scipy、pysam、intervaltree，
  Excel 输出还需要对应 writer。
- `phap_cluster.py:783-821,831-838,876-897` 直接执行 mode `100644` 的 Python
  文件，干净 checkout 会 `Permission denied`；应统一使用 `sys.executable`。
- `utils/phase_reads_assemble_anchor.py:274` 使用个人 chopper 路径且把 `pigz`
  写成 `pigp`。
- 同文件 `:328-340` 硬编码 `/home/pxxiao/...`，其他机器无法运行。
- `PHap.py` 注册了不存在的 `scripts/phap_fill_gaps.py`。
- 并行 worker 失败只记录或被忽略，父流程可能以成功状态继续使用残缺输出。

修复方向：先建立可安装包、依赖清单、外部工具 preflight、统一 subprocess
runner 和 CLI smoke tests。

## 3. P1：显著影响组装质量、泛化能力或可复现性

### P1-1 “LIS” 实现没有计算 query/target 共线性

**状态：已在第五批本地修复。** 活动 cluster 路径现在直接使用严格 locus 阶段
输出的完整 alignment chain。该阶段按 strand 检查 query 与 oriented-target
坐标单调性，以 interval union 计算 query coverage，并在逐 alignment、逐 candidate
和逐 unitig 审计中记录选择依据。旧的 target-gap-only “LIS” 脚本及不再生效的
`--min_lis_*` 参数已删除。

`utils/find_longest_subsequence.py:147-173` 只是按 target gap 切段，未检查 query
坐标递增/递减、链方向或 overlap；这不是最长递增子序列。错误 alignment chain
会进入 allelic table。

建议：按 strand 分组，排序 query 坐标，在 target 坐标上计算带权 LIS/collinear
chain，并以 union coverage、identity、gap penalty 评分。

### P1-2 `top_n` 没有贯穿全流程

CLI 声称支持 ploidy/top_n，但 clustering、reclustering、rescue、输出命名及 dosage
上限多处硬编码 4；`utils/allelic_table_refresh.py:126,136` 也固定按四倍体过滤。
`utils/unchr_recluster.py:345-346` 更硬编码了 chr01-chr12 × group1-group4。

建议：引入一个 `AssemblyConfig(ploidy, chromosome_count, ...)`，所有 group 和 dosage
逻辑由它生成；先明确当前版本“只支持 autotetraploid”，或完成真正泛化后再宣称
支持其他倍性。

### P1-3 allelic table 的 dosage 逻辑前后不一致

**状态：活动 generate 路径已在第七批本地修复。** 完整 bin 候选集现在
通过 0/1 dosage-capacity 动态规划选择，约束为 `sum(dosage) <= ploidy`，
目标为 copy-weighted target interval-union support。缺失/非法 dosage 不参与优化；
多个最优 subset 同分时整个 bin 进入 `ambiguous`，ID 仅用于稳定输出。
每个 candidate 和 bin 都输出 status/reason、copy count、剩余容量和 objective。
bin 内最小 union support bases/覆盖率也已暴露为 CLI 参数并写入审计；
默认值是明示的宽松值，不在代码中暗藏数据集阈值。

旧 refresh bridge 启发式仍保留；它没有 bin alignment support，因此如果补入后
仍超容量，当前会保守地删除该 bin 并记录 `ambiguous_bin_exceeds_ploidy`，
不按 FASTA 长度猜测科学选择。

`utils/allelic_table_generate.py` 的 active path 只是选 top N contig，传入的 contig
type 未参与 active ranking；refresh 先按 contig 数量裁剪，再按 dosage 和固定 4
二次过滤。跨 chromosome 的 `seen_unitigs` 和前后 search window 也没有显式隔离。

建议：以 copy-count capacity 为核心一次完成约束优化，按 chromosome 独立处理，
输出每个 bin 的候选、score、选择及拒绝原因。

### P1-4 中间结果缓存可能静默复用错误输入

**状态：活动 cluster 的两个 existence-only gate 已在第八批本地修复。**
`p_utg_vs_mT2T.paf` 现在只在 stage manifest 完整验证通过时复用；
cache key 包含两个 FASTA 的路径/字节大小/SHA-256、完整命令参数、PHap
版本、minimap2 解析路径和版本。PAF 哈希、manifest schema 或任一签名
不匹配都会重建。输出使用同目录临时文件原子替换，manifest 最后提交；
运行中输入变化时拒绝写 completion marker，并发窗口通过 stage lock 串行。

每染色体 `corrected_allelic_table` 不再缓存：它每次都从本轮 global
table 精确匹配 chromosome 字段并原子覆盖，包括覆盖为空文件。
step 5 的 recluster merge 也只接受本轮 locus 列表显式对应的输出；历史
染色体目录不再通过 glob 混入 rescue，合并来源、大小和 SHA-256 会写入
审计表。其他旧 workflow 的 stage-level manifest 仍需在后续按活动调用
路径逐项迁移。

多个 stage 只判断固定文件是否存在；read phasing 还使用固定 pickle 名。更换输入、
参数或代码后，旧结果仍会被加载。

建议：每个 stage 写 manifest（输入 fingerprint、参数、PHap/tool 版本）并原子提交；
不匹配则明确重算或要求 `--force`。

### P1-5 结果不完全确定

- 多处 `list(set(...))` 后直接输出，顺序受 hash seed 影响。
- mT2T 使用 set 构造输入顺序。
- 多个进程同时 `>> distances.txt`，存在顺序不稳定和并发写风险。
- tie-breaking 多数只比较一个 score，没有稳定的 ID/长度次级键。

建议：所有输出排序，显式稳定 tie-break；worker 各写独立文件，由父进程按 key 合并。

### P1-6 资源配置会严重超订阅

`process_seqkit = process * threads`，随后每个 seqkit 又使用 `threads`，例如
`process=12, threads=10` 最坏会创建约 120 个进程、每个再用 10 threads。cluster
reassignment 还固定 `run_in_parallel(..., 12)`。

建议：用总 CPU budget 统一推导 stage concurrency，记录 peak memory/threads，
不允许多层乘法并行。

## 4. P2：冗余、维护性和性能

- `chr_uncluster_recluster.py` 与 `unchr_recluster.py` 约 500 行高度重复，应合并为一个
  library + 两个薄入口；目前两份实现已出现行为漂移。
- `scripts/util.py` 的两套旧 mT2T redundancy/join 实现已删除；文件中其他历史 helper
  仍需按调用关系继续拆分和清理。
- FASTA、contig type、pickle、subprocess、并行 helper 被重复实现多次。
- `unitig_overlap_ratio()` 为全部 unitig 两两建表，时间和内存均为 `O(U²)`；只需在
  同 bin 的 unitig 间建立 sparse pair。
- rescue 对每个待分配 contig 扫描全部 Hi-C edge，接近 `O(C×E)`；应先建 adjacency。
- 大量逐 record `print("Debugs...")` 会产生巨型日志并拖慢运行。
- `PHap.py:parse_arguments()` 是未使用的第二套 CLI 定义，与实际透传入口已经漂移。
- tracked `scripts/.DS_Store`、无 pycache ignore、示例脚本不可执行等属于仓库卫生问题。

## 5. 文档问题

- README clone URL 仍指向旧仓库 `JiaoLab2021/PHap`。
- README 未说明 Python 包依赖、版本范围、安装方式和外部工具发现规则。
- `Pipeline.md:80` 重复写了 `phap phase_reads phase_reads`。
- `Pipeline.md:88` 把 `--hic2` 也写成了 hic_1。
- 文档给出的 mt2t 参数默认值与当前 CLI 已有差异。
- shell 示例依赖个人绝对路径和手动 conda 环境，不能从干净环境复现。

## 6. 建议的升级顺序

### Phase 0：冻结基线

1. 准备不含敏感数据的小型 synthetic fixture。
2. 保存当前一套代表性 potato 数据的完整命令、工具版本和输出指标。
3. 建立 sequence/read 去向审计脚本，先量化当前丢失、重复和跨 group overlap。

### Phase 1：工程可运行

1. 添加 `pyproject.toml`、console entry point、依赖与可选 integration 依赖。
2. 合并 CLI，修复帮助、缺失 subcommand、路径、权限、退出码和 preflight。
3. 建统一 command runner、日志、output directory、manifest/resume。
4. 建 CI：lint、unit、CLI smoke、synthetic end-to-end。

### Phase 2：先修科学正确性

按 P0-1 至 P0-7 逐项修复，每项一个 regression test；这一阶段不同时大改默认阈值。

### Phase 3：重构核心数据层

建立共享 FASTA/PAF/BAM parser、interval coverage、Hi-C adjacency、dosage model 和
动态 group model；删除重复脚本与旧实现。

### Phase 4：提升组装算法

1. mT2T：oriented overlap graph、contained filtering、conflict-aware path、join QA。
2. allelic table：strand-aware weighted chaining、dosage-capacity selection。
3. clustering：可解释的 multi-evidence score 与 ambiguity margin。
4. read phasing：best-hit/likelihood 分配、collapsed reads 互斥分区、污染审计。

### Phase 5：真实数据 benchmark

新旧版本必须在相同输入上比较：

- completeness（例如 BUSCO/Compleasm 类指标）；
- k-mer completeness、QV、haplotype duplication/switch；
- contig/scaffold continuity、gap 数、端粒完整性；
- Hi-C contact map misjoin 证据；
- reads 回贴覆盖、跨 haplotype contamination；
- peak memory、CPU time、磁盘占用。

N50 只能作为连续性指标之一，不能单独作为“组装效果提升”的结论。

## 7. 第一批建议实施的最小 PR

第一批仅做“可运行 + 可测”，不改科学默认行为：

1. `pyproject.toml`、依赖、CLI entry point；
2. CLI/help/version smoke tests；
3. 统一 Python 子脚本调用和外部工具 preflight；
4. 统一 worker 失败传播；
5. 修复硬编码路径、`pigp`、缺失脚本入口和文档命令；
6. 加入 tiny FASTA/PAF/group fixtures，为下一批 P0 correctness fixes 建安全网。

## 8. mT2T 与 MGA 后续设计依据

PHap 论文中的 mT2T 和 MGA 都生成 haplotype-mixed/mosaic reference。两者在
PHap 流水线中的**功能角色可以互换**：都可为后续 p_utg 定位、allelic table
构建和分型提供一套混合单倍型伪 T2T 坐标。它们不是同一算法，输入契约也不同：

- PHap 当前 mT2T 以 `p_ctg` 为输入，先用 Mash 筛选候选 contig pair，再依据
  minimap2 overlap 构图并连接序列。
- [MGA](https://github.com/ZhangZhenmiao/MGA) 直接以 HiFi reads 为输入，基于 LJA
  的 multiplex de Bruijn graph，执行 read/graph cleaning、detouring、dewhirling、
  decoupling、broken-tip repair、short-edge contraction、scaffolding 和 cognate
  contig deduplication。
- [MGA 论文](https://link.springer.com/article/10.1186/s13059-026-04128-5)
  讨论的是 diploid haplotype-mixed consensus assembly；尚不能据此假定其 bubble
  简化和去冗余规则适用于 autotetraploid dosage。
- MGA 使用 BSD 3-Clause License，但 PHap 不应复制其实现作为本轮修复的一部分。
  后续应将 MGA 作为显式、可选的外部 backend，记录版本与完整命令，并与 PHap
  overlap-graph backend 在相同数据上比较完整性、phase switch/misjoin、
  duplication、read support、连续性和资源使用。

建议将 `phap mt2t` 改为统一 provider 接口：

- `--backend overlap`：输入 `p_ctg`，执行 PHap 原生 Mash/minimap2/oriented
  overlap graph 流程；
- `--backend mga`：输入 HiFi reads，调用固定版本的 MGA；
- 两个 backend 都必须产出规范化的 `mT2T.fa`、运行 manifest、序列 ID mapping
  和质量统计，之后进入完全相同的 allelic-table/cluster 流程。

原生 overlap backend 仍需完成共享 PAF parser、conflict-aware path extraction
和序列物化测试。MGA backend 则以外部工具适配、输出规范化和四倍体数据 benchmark
为主，不要求复制 MGA 内部图算法。
