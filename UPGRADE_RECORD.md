# PHap v2 代码与算法升级记录

更新时间：2026-08-22

本文汇总当前工作区中 PHap v2 的代码结构、已经保留的算法改进、验证后放弃的尝试、现有结果及后续运行建议。各阶段的详细算法和输出定义仍以对应的专题文档为准。

## 1. 当前结论

当前代码的推荐方案是：

1. 用 mT2T 投影、真实区间重叠、dosage 和可选 GFA 构建 allelic unitig table。
2. 在 `03.cluster` 中把 allelic table 作为互斥约束，求解 exact-dosage 分组，再用 Hi-C 做保持约束的局部、Kempe、phase-block 和 phase-interval 优化。
3. 在 `04.recluster` 中用 Hi-C 复核弱 seed、分配其余染色体 unitig，并用 **allelic block 联合分配**恢复单 unitig margin 难以处理的局部互斥关系。
4. 在 `05.rescue` 中只处理 `un_chr.fa`，不改变 `04.recluster` 的结果。
5. Hi-C 默认按两端 dosage 的乘积矫正；保留 `raw` 参数用于严格对照实验。
6. `phase_reads` 采用 read-centric 比对证据；只在 collapsed reads 无法区分时按估算深度互斥均衡分配，再分阶段完成 FASTQ 提取、组装和挂载。

“统计置信度 + allelic block”中的 Poisson/p-value 统计置信度已经从当前源码、参数和新输出中删除。当前 `04.recluster` 只保留 allelic block 的 `constraint_forced` 或 `margin` 两种接受依据。

## 2. 不可改变的方法边界

### 2.1 亲本信息只用于评估

聚类和重分配只允许使用：

- mT2T 投影及 allelic table；
- unitig dosage；
- Hi-C links 和 restriction-site 数量；
- unitig 长度、table 坐标及组装图证据。

亲本 yak/k-mer 信息只能在一次运行全部完成后评估结果，不得参与选组、阈值设定、move 接受、染色体特异参数调整或 allelic table 构建。本文列出的 yak 结果都是后验评估，不是算法输入。

### 2.2 独立证据不循环使用

- Hi-C 不用于创建 allelic table，避免随后用同一证据验证或满足自己创建的约束。
- GFA 直接连接可以否定一对候选 allele，但没有图连接不能证明两者一定等位。
- group bp 只作 QC 和轻量平衡项，不能为了长度相等覆盖强 Hi-C 或 exact-dosage 约束。

### 2.3 保守分配和完整审计

- 低置信度序列默认进入显式 unassigned 输出，不强制猜测。
- 每一阶段都检查 dosage、互斥约束、ID 唯一性和输入/输出 partition。
- 新输出先写临时目录，全部验证通过后再替换正式结果。

## 3. 当前代码结构

| 阶段 | 主要源码 | 职责 |
| --- | --- | --- |
| 总入口 | `PHap.py`、`scripts/phap_cluster.py`、`scripts/phap_allelic_table.py` | 参数解析、阶段编排、并行染色体任务及 `--stop_after` |
| 01 allelic table | `utils/find_collinear_chains.py`、`utils/allelic_table_generate_v2.py`、`utils/validate_allelic_table_gfa.py` | 聚合 PAF、局部投影、atomic interval、dosage/GFA 过滤和审计 |
| 02 chromosome extraction | `utils/extract_chr_from_putg.py` | 按最佳 mT2T 染色体证据分配，保留 `un_chr.fa` |
| 03 cluster | `utils/cluster_allelic_unitigs_v2.py` | exact-dosage 约束求解及多层 Hi-C 优化 |
| 04 recluster | `utils/chr_uncluster_recluster.py` | seed 复核、批量传播、同步 refinement、allelic block 联合分配 |
| 05 rescue | `utils/unchr_recluster.py` | 仅救援染色体未定位 unitig，保持 04 不变 |
| phase reads | `utils/phase_reads_assignment.py`、`utils/phase_reads_assemble_anchor.py` | read-centric 分配、collapsed 深度均衡、单次 FASTQ 分发、分阶段组装和挂载 |
| 测试 | `tests/test_*.py` | 表构建、约束求解、raw/dosage、block、recluster/rescue、原子写出等回归测试 |

旧目录 `PHap` 未被本轮升级修改；新实现集中在 `PHap.v2`。

## 4. 分阶段算法改进

### 4.1 Allelic unitig table

旧方法容易把碎片化、方向改变或相邻而不重叠的投影误当成等位关系。v2 的主要改进为：

1. PAF 先按 primary、alignment length、identity 和 MAPQ 过滤，再按 unitig/染色体聚合局部证据；重叠 PAF 记录不会重复计分。
2. 允许同一 unitig 出现 mixed strand、倒位和非单调局部块，不再要求全长单链共线。
3. 稀疏投影按最大无支持间隙递归切分，保留可靠的局部 `segmented` block。
4. 通过投影端点 sweep 生成 half-open atomic intervals；只有参考坐标真实重叠的 unitig 才能进入同一 row，相邻区间不构成 allele pair。
5. 每个 interval 的 dosage 总和不得超过 ploidy；超过时只有最佳可行子集相对次优解达到 margin 才保留，否则省略。
6. 可选 GFA 检查会删除存在直接 `L` edge 的候选等位 interval，并记录 bubble-like、shared-read 和 neutral 证据。
7. 输出 projection、pair、rejection、summary 和 run manifest，保证每个决定可追踪。

当前默认最小 query 长度仍为 20 kb，没有全局提升到 50 kb。详见 [ALLELIC_TABLE_V2.md](ALLELIC_TABLE_V2.md)。

### 4.2 染色体提取

`02.chr_seq` 将“能否定位到染色体”和“能否成为高置信等位约束”分开处理：

- allelic projection 覆盖不足不会自动把 unitig 从染色体中删除；
- 最佳染色体相对第二名 margin 至少为 0.05 时可救援定位；
- 每条输入序列只能进入一个 chromosome FASTA 或 `un_chr.fa`；
- FASTA 流式读取，并核对 PAF、QC、mT2T 与 p_utg 的长度一致性。

详见 [CHROMOSOME_EXTRACTION_V2.md](CHROMOSOME_EXTRACTION_V2.md)。

### 4.3 初始聚类 `03.cluster`

当前 `03.cluster` 不再是简单按 Hi-C 贪心分组，而是：

1. 将 haplotig、diplotig、triplotig、tetraplotig 转成 1–4 的 exact dosage。
2. 将 allelic row 转成互斥 conflict graph。
3. 用确定性的 DSATUR 风格顺序、forward checking 和 backtracking 求 exact-dosage graph coloring。
4. 对数学上不可行的超 ploidy clique，默认仅放松方向重叠比例和 overlap bp 最弱的 edge；也可用 `fail` 直接停止。
5. 以无量纲目标 `Hi-C cohesion - balance_weight × group-bp relative MSE` 做单 unitig 和 Kempe-component refinement。
6. 对与 Hi-C 强烈矛盾的分配，只允许在 link、group margin、全局目标及弱 edge overlap 条件都通过时释放有限约束。
7. 增加 suffix phase-block 和双边界 phase-interval permutation，修复长距离 phase switch；unitig 始终作为完整记录移动，不切断序列。
8. 每个放松 edge、unitig move 和 block/interval move 都写入独立 TSV。

当前关键默认值包括：`balance_weight=1.0`、refinement 10 rounds、block 最小 gain 0.005、约束重分配最小 links 5、margin 0.10。详见 [CLUSTER_V2.md](CLUSTER_V2.md)。

### 4.4 染色体内重分配 `04.recluster`

这一阶段解决两个问题：弱的 `03.cluster` seed 可能分错，以及大量非 table unitig 尚未进入 haplotype。

1. 只有满足阈值的 `hic_supported` seed 成为不可变 anchor。
2. 其余 seed 先移出传播 anchor，用染色体内 Hi-C 重评；证据不足时默认保留原组，但不让它传播不确定性。
3. 新候选按批次同步接受，避免按 unitig 长度或输入顺序逐个更新造成级联偏差。
4. 传播后进行同步 refinement，检测稳定状态和 oscillation。
5. 常规接受条件为：exact dosage、每个选中组有正支持、selected links 至少 5，并且 weakest-selected 与 strongest-unselected 的 density margin 至少 0.10。
6. 最后读取染色体局部 allelic rows，对弱成员联合枚举 dosage-exact、互斥的配置：
   - 只有一个合法补全时，以 `constraint_forced` 接受；
   - 有多个配置时，最佳配置必须同时达到 block margin 0.10 和最小 links 5，才以 `margin` 接受。
7. 含有 `03.cluster` 已放松 pair 的 allelic row 直接跳过，避免在 04 中重新引入已经判弱的约束。
8. 每个 accepted block 写入 `allelic_block_reassignments.tsv`；summary 记录接受依据和 changed unitig 数。

当前 block exact-enumeration 默认至少需要 1 个 anchor，上限 256 个配置。详见 [RECLUSTER_V2.md](RECLUSTER_V2.md)。

### 4.5 染色体外救援 `05.rescue`

- 候选严格限定为 `02.chr_seq/un_chr.fa`；04 内已分配和显式 deferred unitig 均不可变。
- 先比较最佳与次佳 chromosome density，再在最佳染色体内比较 dosage 边界的 group density。
- 默认 chromosome margin 和 group margin 都为 0.10，最小 links 为 5。
- `other`、`replotig` 或缺少 dosage 的候选默认 defer。
- 最终 assigned 与 unassigned 必须互斥且并集等于原始 p_utg。

05 不使用 chromosome-local allelic block，因为这些候选本来就没有可靠的本地染色体归属。详见 [RESCUE_V2.md](RESCUE_V2.md)。

### 4.6 工程可靠性

- v2 不再接收或生成 CLM。03–05 的 Hi-C 分配只使用 `full_links.pkl`；后续
  `phase_reads` 会把分组后的 Hi-C reads 重新比对到新组装，不消费历史拆分 CLM。
- 03–05 的 summary 均记录输入、参数、计数和 validation。
- 失败时保留旧正式输出，减少半成品覆盖。
- `--stop_after table|chr_seq|cluster|recluster|rescue` 支持分阶段运行和复核。

## 5. Hi-C 阈值和 dosage/raw 参数

统一参数为：

```text
--hic_link_normalization dosage   # 默认
--hic_link_normalization raw
```

standalone utility 使用连字符形式 `--hic-link-normalization`。03、04、05 使用同一个模式，并将模式写入各自 summary JSON。

dosage 模式的 link 为：

```text
adjusted_links(u, v) = raw_read_pairs(u, v) / (dosage(u) × dosage(v))
```

因此阈值 5 表示所选归一化尺度上的 5，而不总是 5 个原始 read pairs：

| 两端 dosage | 达到 adjusted=5 所需 raw links |
| --- | ---: |
| 1 × 1 | 5 |
| 1 × 2 | 10 |
| 2 × 2 | 20 |

`raw` 模式直接令 adjusted value 等于原始 pair 数。为兼容旧下游，部分 TSV 字段仍包含 `adjusted` 名称，但数值遵循 summary 中记录的模式。

group margin 定义为：

```text
(weakest selected group density - strongest unselected group density)
/ weakest selected group density
```

短、近乎相同或 Hi-C 稀疏的 unitig 很容易达不到固定 margin。当前解决办法不是降低全局阈值，而是：强证据正常分配；有 anchor 的 allelic block 联合解；其余显式 defer。这样不会因为放宽一个全局阈值而扩大系统性误分配。

## 6. 输入顺序和 50 kb 长度阈值的设计结论

### 6.1 Allelic table 中 unitig/row 顺序

- row 内 unitig 在主要搜索中会按 ID、长度或约束度排序，成员书写顺序通常不改变约束含义。
- `03.cluster` 的 anchor row 在证据和区间长度完全相同时以原 row 顺序作最终 tie-break；启发式搜索在完全等价解中也可能选择不同但同分的 group label。
- `04.recluster` 的 block 先按已有 anchor 数、坐标和 block ID 排序；一个已接受 block 会锁定其 movable 成员。因此重叠 block 在完全打平的特殊情况下可能有顺序效应。

结论是：相同输入下当前实现是确定性的，普通 row/member 顺序不会系统性改变结果；但不应把“任意打乱 allelic table 后 bitwise 完全一致”当成算法保证。对顺序敏感性的正式检查应随机置换 table，多次运行后对齐任意 haplotype label 再比较。

### 6.2 是否只保留大于 50 kb 的 unitig

当前 12 条染色体的 `03.cluster` seed 共 2,601 个 unitig；其中 541 个小于 50 kb，占数量 20.80%，但仅 23,134,920 bp，占 seed 序列 0.942%。当前 03 初始 anchor 中最短 unitig 也超过 2 Mb，所以 50 kb 过滤不会改变本次数据的 anchor 本身，但可能删除短 unitig 提供的局部互斥边、减少 table recall，并把问题推迟到 Hi-C 更弱的 04 阶段。

因此目前未加入全局 50 kb 硬过滤。如果以后比较，应新增“table seed 专用最小长度”参数，而不是直接提高共享 PAF 的 `min_query_length`；后者会同时影响 chromosome extraction，混淆实验因素。

## 7. 已完成验证和算法取舍

### 7.1 当前全染色体基线

目录：`../02.cluster.v2.adaptive.nogfa`

- 12 条染色体的 03、04、05 均已运行，最终有 48 个 chromosome-haplotype groups。
- 03 yak 汇总：345,594 / 24,197,556，错误率 1.428219%。
- 04 yak 汇总：213,608 / 24,304,274，错误率 0.878891%。
- 这些是加入 allelic block 之前的完整基线；不能冒充当前 block-only 源码的全染色体结果。

### 7.2 “统计置信度 + allelic block”试验：放弃

目录：`../02.cluster.v2.statblock.validation/04.recluster`，测试染色体 chr02、chr05、chr08。

| 方案 | wrong / informative | 错误率 |
| --- | ---: | ---: |
| 原 04 基线 | 155,529 / 5,471,327 | 2.842619% |
| 统计置信度 + block | 157,144 / 5,479,835 | 2.867678% |

chr05 和 chr08 略有改善，但 chr02 明显退化。`utg000927l` 和 `utg001320l` 即使得到很小的 Poisson p-value，仍分别增加约 6,463 和 3,520 个错误 marker。这说明独立 Poisson 假设低估了系统性、相关性和过度离散的 Hi-C 噪声；单纯收紧 p-value 也不能解决模型失配。

最终处理：删除 `hic_confidence.py` 源码、p-value/z-score 计算、相关 CLI 参数和新输出字段。历史试验目录可能仍含旧格式统计列，仅用于追溯，不代表当前代码接口。

### 7.3 仅 allelic block 试验：保留

目录：`../02.cluster.v2.blockonly.validation/04.recluster`，阈值仍为 links 5、margin 0.10，没有用统计显著性作为接受条件。

| 染色体 | 原 04 | block-only | wrong 变化 |
| --- | ---: | ---: | ---: |
| chr02 | 12,851 / 322,593 (3.983657%) | 11,544 / 322,593 (3.578503%) | -1,307 |
| chr05 | 46,602 / 2,217,088 (2.101946%) | 42,827 / 2,217,088 (1.931678%) | -3,775 |
| chr08 | 96,076 / 2,931,646 (3.277203%) | 95,158 / 2,931,646 (3.245890%) | -918 |
| 合计 | 155,529 / 5,471,327 (2.842619%) | 149,529 / 5,471,327 (2.732957%) | -6,000 |

三个染色体各接受 6 个 block、改变 6 个 unitig，validation violation 均为 0：

- chr02：4 个 `constraint_forced`，2 个 `margin`；
- chr05：6 个 `margin`；
- chr08：1 个 `constraint_forced`，5 个 `margin`。

相对错误 marker 降低 3.86%，绝对错误率下降 0.109663 个百分点。因此 allelic block 已保留到当前代码。清理统计代码后重新运行 chr02，assignment 与 block-only 验证结果完全一致。

### 7.4 dosage/raw 对照：尚未完成

目录：`../02.cluster.v2.norm_compare_20260822`

该实验在代码快照上分别用 `dosage` 和 `raw` 从 03 运行到 05，并计划对 03、04 运行 yak。它用于隔离 Hi-C normalization 这一变量；快照早于后续 block-only 最终清理，不能用来验证 block 改进。

截至本文核对时：dosage/raw 的 03 和 04 summary 均已生成，dosage 的 05 已生成；raw 的 05 summary 不存在，yak 也未完整结束，后台进程已经不在运行，但 status 文件仍停留在 `RUNNING`。所以当前不报告 dosage/raw 优劣，需完整续跑或重跑并确认两个模式的 03、04、05 和 yak 都结束后再比较。

## 8. 当前数据中的已知个案

- chr04 g4 偏长主要由 `utg000213l` 造成。该 unitig 与 mT2T 只能比对到一小部分，现阶段按用户决定视为组装异常并暂时忽略，不针对它修改聚类规则。
- chr10 旧的 g3 yak 较差主要与 `utg000982l` 的组别有关；04 的 seed review 将其从 g3 调整到 g1后明显改善，说明弱 seed 复核是必要的。
- chr08 错误呈区间成簇，且 chr05 也有相似的高相似重复/近等位区域；很多相关 unitig 较短、Hi-C 弱或接近对半。原因更符合上游 assembly/dosage 模糊与局部 Hi-C 不可分，而不是单纯 allelic-table 坐标整体错误。
- chr02 的错误区间更长、更分散，既有 projection/dosage 问题，也有 04 丢失局部 allelic context 的影响；block-only 在 chr02 的改善支持后一部分判断。

这些判断使用亲本 marker 定位错误来源，但没有把亲本信息反馈到算法分配中。

## 9. 当前推荐运行和验收顺序

1. 先运行单元测试。
2. 用当前代码完整重跑 03–05，输出到新的结果根目录，不覆盖历史验证目录。
3. 逐染色体检查 03 的三个 validation 计数及 relaxed constraints。
4. 检查 04 的 dosage/partition/fixed-seed/refinement validation，并汇总 `allelic_block_reassignments.tsv`。
5. 检查 05 的 complete partition 和 immutable-04 validation。
6. 最后运行 yak；先对齐任意 haplotype group label，再比较错误率和错误区间。
7. dosage/raw 对照需使用同一代码快照、同一输入、同一阈值和独立输出目录，不能把不完整的当前目录作为结论。

详细的当前工作区基础验证见 [VALIDATION_CURRENT_DATA.md](VALIDATION_CURRENT_DATA.md)。

## 10. 尚未完成的工作

- 用当前“无统计置信度、保留 allelic block”的源码完整重跑全部 12 条染色体 03–05，并生成新的正式 yak 汇总。
- 完成或从干净目录重跑 dosage/raw 对照。
- 如需验证 table 顺序鲁棒性，增加随机 row/member permutation 的重复试验和 group-label 对齐比较。
- 如需验证 50 kb 方案，增加 table-seed-only 参数并做控制变量实验，不改共享 PAF/染色体定位阈值。
- 使用后续提供的真实 BAM/FASTQ 校准 `phase_reads` 的比对过滤与 margin 参数，并评估内存、I/O、分组深度、组装连续性和 Hi-C 挂载结果。
- variant-based read phasing 暂缓；待无亲本信息的现有方案完成真实数据验证后再单独升级。

## 11. phase_reads 分步升级（已实现，待真实数据验证）

详细方案见 [PHASE_READS_V2.md](PHASE_READS_V2.md)。当前实现遵循以下边界：

- 亲本信息不参与 reads 分配，只能在整套方法固定后用 yak 评估；
- 一个 HiFi/ONT read 最多进入一个 group，Hi-C 两端始终共同保留或共同 defer；
- collapsed-only reads 不再按 unitig 重复抽样，而是在候选 haplotype 间互斥分配，并以“已分配 read bp / group unitig bp”估算深度做确定性均衡；
- FASTQ 每类只扫描一次，按 assignment table 同时写出所有 group；
- `assign`、`extract`、`assemble`、`scaffold` 可分阶段运行和安全续跑，输入或已完成阶段参数变化时拒绝复用旧结果；
- 外部工具失败会终止阶段，临时结果不会替换上一次完整输出。

phase_reads 的合成测试已经扩展到磁盘 evidence、BAM offset 续跑、group
选择时保留 sibling haplotypes，以及单 group 组装 checkpoint。chr01 HiFi
assignment 已完成首次真实 BAM 工程验证；ONT、Hi-C、FASTQ extraction、组装
和挂载仍待真实数据验收，因此尚不能作为最终生物学效果结论。

2026-08-23 增加以下工程升级：

- coordinate-sorted BAM 按记录写入 SQLite evidence，支持从已提交的 BAM
  virtual offset 继续；按 read ID 在磁盘聚合，不再同时保存三类数据的全部
  Python alignment objects；
- HiFi、ONT、Hi-C 分别产生 assignment SQLite index、TSV 和完成标记，完成
  一类后立即释放其 evidence 数据库；extract 同样按数据类型 checkpoint；
- hifiasm 和 HapHiC 改为每个 group 独立原子输出和 checkpoint，一个 group
  失败不再清除其他已经完成的 group；
- 新增 `--groups`、`--chromosomes`。单独请求某个 group 时仍使用同染色体
  全部 sibling groups 进行公平分配，仅限制最终输出；
- 新增 `phase_reads.log`、四阶段进度、BAM/FASTQ 记录速率和 group 完成比例。
- 新增 `--resume --rerun-from assign|extract|assemble|scaffold`，可保留指定
  阶段之前的 checkpoint 并主动重跑后续流程；只调整 jobs/threads 不会使已有
  结果 checkpoint 失效。
- 新增 `--data-types hifi|ont|hic`，允许三类 BAM 按完成顺序独立测试；只有
  被选择的数据类型才要求对应 BAM/FASTQ。
- 新增 `--temp-dir`：大型 evidence 数据库、FASTQ 提取暂存目录以及组装/挂载
  暂存目录均可放到指定位置；同时为 SQLite 和外部子进程设置临时目录。该参数
  属于运行时参数，不会使已经完成的生物学 assignment checkpoint 失效。输入
  BAM/FASTQ 始终只读，跨文件系统发布结果时先复制到输出目录的隐藏文件再原子
  替换，避免暴露不完整结果。

### 11.1 chr01 HiFi-only 首次真实数据测试（2026-08-23）

输入为 `01.hifi.data/HiFi.p_utg.sort.bam`（39 GB，coordinate-sorted），使用
`02.cluster.v2/05.rescue/group.reassignment.cluster.txt` 的 chr01 四个 groups，
输出在 `03.phase_reads.v2.hifi_chr01_test_20260823`。本次目的是验证工程流程和
reads 深度均衡，不代表当前最新版聚类算法的最终生物学结果。

- BAM 共读取 9,346,596 条 alignment records；保留 512,714 条 chr01 group
  alignment evidence，聚合为 512,575 个 reads；
- 357,296 个 reads（69.7061%）由 alignment evidence 明确分配；155,272 个
  reads（30.2925%）按 collapsed depth balance 分配；仅7个 reads defer；
- 运行中人为中断后从第150,000条 BAM checkpoint 成功续跑，证明 virtual
  offset checkpoint 可用；最终临时 evidence SQLite 已自动清理；
- 最终输出约192 MB，其中 assignment SQLite 约188 MB、压缩TSV约13 MB；
- 四组估算HiFi深度分别为24.3960、24.3960、25.6827、25.6827×；总体深度
  CV为2.569%，最大/最小为1.0527。group1/2及group3/4内部几乎完全一致，
  剩余组间差异由不同 collapsed candidate sets 的可分配范围约束；
- 7个 defer 均为 `ambiguous_unique_group_evidence`，没有发生强制猜测。
