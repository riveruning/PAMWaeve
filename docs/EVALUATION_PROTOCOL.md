# 性能评价方案（冻结版）

状态日期：**2026-10-03**
冻结提交：`54384e6`（本轮改动前的仓库 HEAD；本轮分析修复后的实际源码版本
记录在 `results/results.csv` 的 `code_version` 与 `code_source_sha256`）
冻结清单：`benchmarks/multi_evidence_v2/manifest_final22.json`（**22 系统**，
本轮完整参考运行所用）
清单 SHA-256：见 `data/parsed/bench_final22/benchmark_report.json` 的
`manifest_sha256` 字段（由运行器在运行时重算并写入）。

> **在本提交包内复跑时用哪份清单？**
> 用 `manifest_redistributable21.json`（**21 系统**）——它的全部工件都随包分发。
> `manifest_final22.json` 多出的 SpCas9 `paired_evidence` 输入为
> **CC-BY-NC-ND-4.0**，不随包分发，因此对它跑 `--check-files` 会**非零退出**。
> 两者差异**恰为这 1 个系统**；其余 21 个系统的指标与参考结果**逐一相同**
> （实测 0 项差异）。详见 [SUBMISSION.md](SUBMISSION.md) 第 2、3 节。

本文件在**查看本轮任何模型分数之前**写成，用于固定样本纳入／排除条件、
候选空间、参数、参考答案、指标、并列规则与弃权处理。
下文所有设置在本轮运行后**未做任何修改**。

## 1. 本轮回答什么、不回答什么

回答：

- 在已冻结的 22 个 Cas9 系统上，`family_prior` / `protein_only` /
  `spacer_only` / `fusion` 四路方法各自的候选覆盖与排名表现；
- 哪些系统上两路证据一致、冲突或不足；
- 现有材料是否支持「融合带来改善」。

不回答（**尚未证明**）：

- 模型泛化能力；
- 融合带来的准确率提升；
- 切割活性、编辑效率或成功概率。

最后一条尤其重要：本仓库的分数是**证据兼容度与排序值**，
不是活性或效率的估计，本轮不为其重新命名。

## 2. 样本单位与关联结构

**评价单位 = 一个 Cas9 系统（system）**，不是「一行 gold」、
不是「一个 PAM」、不是「一个命中」、不是「一个 contig」。

- 同一蛋白的多个可接受 PAM 属于同一个系统，按 PAM 长度组成 spectrum，
  **每个系统只计一次**；
- `published_flanks` 的每条记录是「一条 spacer 一条聚合侧翼」，
  `support_target_count` 因此**刻意不填**（不是 0）：见第 5 节；
- 同一病毒的多个 contig、同一 contig 上的多个命中都**不**算独立样本。

已知关联（必须在读表时一起看）：

| 关联 | 系统 | 说明 |
|---|---|---|
| 同一论文 | `asp-acidiphilium-21-60-14-protein-only`、`cba-caulobacterales-protein-only` | DOI `10.1089/crispr.2024.0013`，**不是两篇独立论文** |
| 同一论文来源 | 18 个 `published_flanks` 系统 | 聚合侧翼来自 CRISPRCasDB（DOI `10.1186/s13059-021-02495-9`） |
| 同源蛋白 | `nme-*` 与 `nme2-*` | 同属 *N. meningitidis* Cas9，非独立系统对 |

## 3. 赛道划分与纳入／排除条件

按输入条件分四个赛道，**严格独立资格按赛道分别判定**，
蛋白独立不代表 spacer/target 也独立。

| 赛道 | 本轮系统数 | 纳入条件 | 排除条件 |
|---|---:|---|---|
| `protein_only` | 22 | 有 Cas9 全序列 + 实验 PAM 谱或功能共识 | 无蛋白序列 |
| `paired_evidence` | 1 | 同一系统同时有定向 spacer **和**原始 target contig | target 未聚类到病毒层级 |
| `published_flanks` | 18 | 有论文发布的 per-spacer 聚合侧翼 | 无原始 target 标识 |
| `end_to_end_genome` | 0 | 需完整宿主基因组 + 冻结注释流程 | **本轮无符合条件系统**，不评价 |

`end_to_end_genome` 本轮为 0 个合格系统，这是**容量缺口**，不是失败：
基因组工作流（`scripts/run_genome_agent.py`）已在真实基因组上跑通并留证，
但没有冻结的实验 gold 谱与之配对，因此不进入排名评价。
它只作为**工程演示**记录，见 `CASE_SELECTION.md`。

### 3.1 全部保留，不做选择性剔除

- 预期表现差的系统**不剔除**：`cj-*`、`wvi-*`、`tde-*`、`sth1a-*` 等
  全部保留在正式表中；
- 弃权（abstain）与失败的样本**保留并单独报告**，不删除以提高指标；
- 用于探索的 `spacer_fallback` / `legacy` 融合策略**不是**正式结果，
  正式结果一律 `--fusion-policy abstain`。

### 3.2 严格独立资格定义（沿用已冻结的审计规则）

一个 system 在某赛道上严格独立，当且仅当：

1. `status=ready`；
2. gold 为实验 activity 或实验 functional spectrum；
3. 有 DOI；
4. `protein_training_exposure=none`；
5. target 已聚类到独立病毒层级（`paired_evidence` 赛道要求）；
6. protein / spacer / target 与方向信息齐全。

**本轮实际结果**：严格独立系统 3 个，全部是 `protein_only`
（Cj4Cas9、AspCas9、CbaCas9）；严格独立 `paired_evidence` 为 **0**，
低于预注册门槛 5。因此本轮**不产生泛化结论**。

## 4. 候选空间与参考答案

- **正式排名评价**使用完整具体 A/C/G/T 候选空间，按 PAM 长度分别生成：
  长度 3→64、4→256、5→1024、6→4096、7→16384、8→65536，上限 100000；
- **不同长度之间不比较名次**，绝不把三碱基分数与六碱基分数放同一榜；
- IUPAC gold 使用**集合语义**：gold `NGG` 时 `AGG/CGG/GGG/TGG` 均算命中；
- **展示用小候选清单（如 `NGG,NAG,NGA,NRG`）不得代替正式排名评价**，
  只用于 `results.csv` 的展示与人工阅读；
- 用参考答案长度限定候选空间（如已知 gold 是 7 nt 就只生成 7 nt 候选）
  意味着评价是在**已知 PAM 长度**条件下进行的，**不能**称为完整 PAM
  发现能力。本报告在每处涉及排名的地方都注明长度条件。

参考答案来源分两类，**在结果中必须区分标签**：

| 来源标签 | 含义 | 系统 |
|---|---|---|
| `engineering_known_spectrum` | 工程阳性对照的已知 PAM | `spcas9-*` |
| `experimental_functional_spectrum` | 论文实验功能谱 | 严格独立 3 个 |
| `published_aggregate_flank` | **文献侧翼汇总，非原始匹配** | 18 个 `*-published-flanks` |

## 5. 参数与缺失值规则

固定参数（本轮未改）：

```text
--model            cas9_full
--device           cpu
--fusion-policy    abstain
--bootstrap        2000
--max-concrete-candidates 100000
pam_side           downstream（全部系统）
max_mismatches     4（paired_evidence）；published_flanks 用 --min-effective-spacers
```

缺失值规则：

- **缺失证据一律表示为缺失／不适用**，**不填成 0**；
- `protein_only` 系统的 `spacer_only` 与 `fusion` 标记 `not_applicable`，
  其指标为 `null`（TSV 中为空单元格），并**从该方法的覆盖率分母中剔除**；
- `published_flanks` 的 `support_target_count` 缺失，因为一条聚合侧翼是
  **一条记录**而不是独立观测到的靶标；此处置为缺失而非 0；
- 只有真的算了且得 0 分的才写 0。

## 6. 指标与并列规则

主指标：

- **MRR**；
- **Recall@1 / Recall@3 / Recall@5**；
- **覆盖率（coverage）**：该方法在适用系统上真正给出排序的比例；
- **弃权率与弃权原因**。

并列处理（保守，沿用已冻结实现）：与最佳 gold **同分**的非 gold 候选
一律视为排在 gold 之前；同时保留最乐观名次用于审计字段
`best_gold_rank_optimistic`。这样「全部候选同分」不能免费拿到 Recall@1。

区间估计：仅在样本量与独立单位支持时给出。本轮严格组 n=3 且其中 2 个
共享同一论文，**不给 bootstrap 区间**，只报逐系统数值。

## 7. 方法比较的合格样本一致性

- 方法对比**只在适用系统上**进行；
- `protein_only` 在 22 个系统上都适用；
- `spacer_only` / `fusion` 只在 19 个有 spacer 通道的系统上适用
  （1 个 `paired_evidence` + 18 个 `published_flanks`）；
- 融合相对基线的配对比较**只在两路都适用且都给出排序**的系统上进行；
- **不删除弃权或失败样本来提高指标**。

## 8. 运行后不得修改的项目

以下内容在本轮结果产生后**不得**再调整：

候选空间、gold 定义、`abstain` 门控阈值、并列规则、样本纳排条件、
指标定义、方法适用性判定。

若将来要改，必须新建 manifest 版本与新的冻结日期，并同时保留旧结果。

## 9. 复跑命令

```bash
cd <repo>
P=python3                     # 见 SUBMISSION 环境说明
export PYTHONPATH=.pylibs:.reference/Protein2PAM:src:.
export HF_HOME=$PWD/data/checkpoints/hf

# 1) 不加载模型即可校验清单与全部 60 个工件哈希
$P scripts/validate_benchmark_manifest.py \
   --manifest benchmarks/multi_evidence_v2/manifest_final22.json --check-files

# 2) 正式评价（真实推理，CPU）
$P -u scripts/run_multi_evidence_benchmark.py \
   --manifest benchmarks/multi_evidence_v2/manifest_final22.json \
   --device cpu --outdir data/parsed/bench_final22

# 已有推理产物时只重算指标
$P -u scripts/run_multi_evidence_benchmark.py \
   --manifest benchmarks/multi_evidence_v2/manifest_final22.json \
   --device cpu --outdir data/parsed/bench_final22 --reuse-existing
```
