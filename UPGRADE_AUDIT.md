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

下面的问题清单仍保留为基线审计记录。P0-1 至 P0-7 的科学正确性修复不属于第一批
工程 PR，尚未实施。

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

`utils/phase_reads_assemble_anchor.py:249-265` 对每个 group 独立调用
`split_reads()`，每次使用相同 seed。一个 diplotig/triplotig/tetraplotig 同时出现
在多个 group 时，各 group 会得到同一个 `1/dosage` 子集，而不是不同的互斥子集；
其余 reads 被丢弃。这会直接造成 haplotype 间污染和覆盖不足。

修复方向：以 `(unitig, data_type)` 为单位只洗牌一次，再按目标 group 数量做完整、
互斥、确定性的分区；输出 union、intersection、unassigned 统计并加回归测试。

### P0-2 read phasing 接收所有长读长比对

`utils/phase_reads_assemble_anchor.py:114-129` 未过滤 unmapped、secondary、
supplementary、duplicate、MAPQ 或最佳比对。一个多重比对 read 可进入多个 unitig，
随后进入多个 haplotype。Hi-C 路径也没有完整、可配置的质量策略。

修复方向：先定义并测试 primary/best-hit/ambiguous 策略；对跨 group 的 read 冲突
单独报告，不允许通过 `set.update()` 静默复制。

### P0-3 dosage 归一化计算后未使用

`scripts/phap_cluster.py:138-150` 调用 `adjust_hic_signals()` 后丢弃返回值，后续仍
使用原始 links。当前 dosage 权重实际上没有生效。

修复方向：将 Hi-C scoring 提取为纯函数，用合成 links 验证 dosage、RE sites、
group size 归一化，并明确使用 raw 还是 normalized score。

### P0-4 无 Hi-C 信号的 unitig 被当作错误分型并删除

`scripts/phap_cluster.py:645-653,681-685` 从四个 group 中直接过滤不在
`dic_contig_hic` 的 unitig。缺少 Hi-C links 不等于序列错误，会降低完整性并导致
无审计的数据丢失。

修复方向：保留为 `unassigned/low_support`，在最终审计表中记录原因；只有明确的
生物学或比对证据才能删除序列。

### P0-5 mT2T reverse join 坐标转换错误，单 alignment 还可能崩溃

`scripts/util.py:798-801` 先覆盖 `tstarts`，再用已覆盖值计算 `tends`，反向坐标不是
正确的 `[tlen-old_end, tlen-old_start)`。`calc_stats()` 在 LIS 只有一个点时调用
`linregress()`，也没有边界处理。

修复方向：先保留 old starts/ends 后成对转换；为正向、反向、单 alignment、
contained 和冲突 overlap 建立精确序列断言。

### P0-6 mT2T 目前不是冲突感知的多 contig 组装

`scripts/util.py:871-884` 按 PAF 字典顺序尝试独立 pairwise merge。每条 contig
只能碰巧加入第一个成功 pair，不能处理 A-B-C 链、分支、环、方向冲突或多个候选
overlap；结果还依赖遍历顺序。

修复方向：构建带方向和置信度的 overlap graph，先去 contained contig，再进行
冲突检测与 path extraction；有歧义时输出 graph/报告，不猜测连接。merge 后必须
以原始 reads、覆盖连续性和 Hi-C 支持复核。

### P0-7 PAF 覆盖率和 primary 标记解析不可靠

多个模块把 alignment match length 直接累加，重叠 alignment 会重复计数，比例可
被放大；多个位置使用 `parts[16][-1]` 判断 primary/supplementary，但 PAF optional
tag 没有固定列号。

涉及：

- `scripts/util.py`
- `utils/find_longest_subsequence.py`
- `utils/allelic_table_generate.py`
- `utils/extract_chr_from_putg.py`

修复方向：实现一个共享 PAF parser，按 tag 名解析 `tp`，覆盖率使用 interval
union，并对缺失 tag、乱序 tag、overlap alignment 建 fixture。

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

`utils/allelic_table_generate.py` 的 active path 只是选 top N contig，传入的 contig
type 未参与 active ranking；refresh 先按 contig 数量裁剪，再按 dosage 和固定 4
二次过滤。跨 chromosome 的 `seen_unitigs` 和前后 search window 也没有显式隔离。

建议：以 copy-count capacity 为核心一次完成约束优化，按 chromosome 独立处理，
输出每个 bin 的候选、score、选择及拒绝原因。

### P1-4 中间结果缓存可能静默复用错误输入

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
- `scripts/util.py` 同时保留 `remove_redundancy()` 和
  `remove_redundancy_v2()`，并含大量失效注释、重复 import 和未使用变量。
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
