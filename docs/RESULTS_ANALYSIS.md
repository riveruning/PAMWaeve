# 性能评价结果与限制

状态日期：**2026-10-03**
方案：[EVALUATION_PROTOCOL.md](EVALUATION_PROTOCOL.md)（先于本轮分数冻结）
结果目录：`data/parsed/bench_final22/`
运行器：`scripts/run_multi_evidence_benchmark.py`（`--fusion-policy abstain`）

## 1. 结论摘要

```text
FRAMEWORK_VALIDATED_NO_STRICT_COHORT
```

- 22 个系统全部完成真实模型推理，无一失败；
- 严格独立系统 **3 个，全部是 `protein_only`**；严格独立
  `paired_evidence` 为 **0**（门槛 5）；
- **没有证据支持「融合带来改善」**（理由见第 4 节，这是本轮最重要的负结论）；
- 严格独立蛋白单路 MRR **0.008 ~ 0.067**，是诚实的负结果。

## 2. 样本数量与独立资格

| 项 | 数量 | 说明 |
|---|---:|---|
| 清单系统 | 22 | `manifest_final22.json` |
| 实际完成评价 | 22 | 全部真实推理 |
| `protein_only` 赛道 | 22 | 全部适用 |
| `paired_evidence` 赛道 | 1 | 仅 SpCas9 工程对照 |
| `published_flanks` 赛道 | 18 | 文献聚合侧翼 |
| `end_to_end_genome` 赛道 | **0** | 无配对冻结 gold，不评价 |
| **严格独立（全部赛道）** | **3** | 均为 `protein_only` |
| **严格独立 `paired_evidence`** | **0** | **阻塞原因见下** |

### 严格独立样本不足的阻塞原因

1. 22 个系统中 **19 个蛋白**与 Protein2PAM 官方训练序列完全相同或近邻
   ≥90%，基础模型的既有暴露**不能**通过后续排除消除；
2. 18 个 `published_flanks` 系统的证据来自 CRISPRCasDB 的**聚合侧翼**
   （一条 spacer 一条记录），**没有原始 target 标识**，
   无法做病毒层级聚类，因此 target 独立性不可评估；
3. 唯一有原始 target 的 SpCas9 系统，蛋白是**精确训练暴露**；
4. 3 个严格独立系统**没有任何 spacer/target 材料**，不构成配对证据。

**未为凑数量开展无期限数据采集**。已有样本继续用于工程回归，
但不据此声称泛化。

## 3. 方法汇总（22 系统，真实推理）

`n/a` = 该方法对该系统不适用（该系统无 spacer 通道），
**已从覆盖率分母中剔除**，不是 0 分。

| 方法 | 覆盖率 | MRR (all) | Recall@1 | Recall@3 | Recall@5 | 不适用 | 适用系统 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `family_prior` | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0 | **1** |
| `protein_only` | 1.000 | 0.683 | 0.636 | 0.682 | 0.727 | 0 | 22 |
| `spacer_only` | 0.947 | 0.700 | 0.579 | 0.842 | 0.842 | 3 | 19 |
| `fusion` | 0.474 | 0.447 | 0.421 | 0.474 | 0.474 | 3 | 19 |

> **`family_prior` 的 1.000 不具信息量**：只有 1 个系统预注册了
> `family_prior_pam`（SpCas9 的 NGG），n=1 且该 gold 是经典已知 PAM。
> 它**不能**与 protein/spacer 的 22/19 系统结果并列比较。

## 4. 融合是否改善？——没有证据支持

这是本轮必须讲清楚的一条。

### 4.1 配对比较实际只有一个非零样本

融合与蛋白的配对比较 **n=9**（两路都适用且都给出排序的系统）：

| 统计量 | 值 |
|---|---:|
| 平均 MRR 差 | **+0.1083** |
| **中位数 MRR 差** | **0.0000** |
| 融合更好 / 持平 / 更差 | **1 / 8 / 0** |

平均值的 **+0.1083 完全来自单个系统** `ssa-jim8777-published-flanks`
（+0.974，其余 8 个恰好都是 0.000）：

```text
0.974 / 9 = 0.1082
```

`ssa-jim8777` 的细节：gold `NNAGAAA`（长度 7，16384 个候选），
蛋白单路把 gold 排到 **第 39 位**（top-1 是 `AAAAAAA`），
spacer 单路排 **第 1**，融合因此也排第 1。

**因此「融合提升 MRR」实际上是一个系统的单点效应**，
bootstrap 95% 区间 `[0.000, 0.325]` **包含 0**。
按本项目的样本纪律，这**不构成**融合改善的证据。

### 4.2 门控更多是在「不合并」，而不是在提升

在 19 个有 spacer 通道的系统中，正式 `abstain` 策略的处置：

| 门控结果 | 系统数 |
|---|---:|
| `joint_ranking_top_agreement`（联合排序） | 9 |
| `abstain_severe_top_disagreement`（严重冲突，弃权） | 9 |
| `abstain_missing_scored_source`（源不足，弃权） | 1 |

其中 **7 个被弃权的系统，蛋白单路本来就把 gold 排在第 1 位**
（`sth1a`、`ain`、`sdy`、`tsp`、`nme2`、`psp`、`wvi`）。
也就是说，如果把这些系统强行纳入融合统计，融合的覆盖率与指标会**更低**；
门控的作用是**避免错误合并**，而不是制造性能提升。

`spacer_only` 的 MRR 0.700 高于 `protein_only` 的 0.683，
但这是**同源 published aggregate flanks** 的结果，
不能替代独立噬菌体 target 或切割实验证据。

## 5. 严格独立组结果（逐系统，不汇总）

| system_id | PAM 长度 | gold | 候选数 | gold 名次 | MRR |
|---|---:|---|---:|---:|---:|
| `cj4-campylobacter-jejuni-414-protein-only` | 6 | `NNNGRY` | 4096 | 28 | 0.036 |
| `asp-acidiphilium-21-60-14-protein-only` | 7 | `NNNNGCA` | 16384 | 122 | 0.008 |
| `cba-caulobacterales-protein-only` | 7 | `NGWNCCA` | 16384 | 15 | 0.067 |

（候选数、名次与 MRR 均取自 `data/parsed/bench_final22/<system_id>/evaluation.json`
的 `methods.protein_only.by_length`。）

**限制**：n=3；Asp 与 Cba 共享同一论文（独立来源实际为 2）；
Cj4Cas9 是 0.8944 一致性的**门槛边界**样本。
**不给区间估计，不做泛化结论。**

三者在 6–7 nt 全候选空间中均未把 gold 排进前列，
这是本项目在**目前唯一可用的独立样本**上的真实表现，如实保留。

## 6. 失败与冲突案例分析

### 6.1 严重冲突（两路 top 候选差异过大）— 9 个

典型是 `wvi-dsm16922`：蛋白 MRR 1.000、spacer MRR 0.125，两路 top 不一致，
门控弃权。这类系统的正确读法是**保留两路各自结论**，而不是合成一个分数。

### 6.2 方向需复核 — 1 个

`nme2-de10444-published-flanks` 的 `orientation_pass=false`
（`orientation_warning` 触发），其 spacer MRR 为 0.000。
在方向未确认前，该系统的 spacer 侧结论**不可用**，已在表中标注 `review`。

### 6.3 源不足弃权 — 1 个

`cj-nctc11168-published-flanks` 只有 **2 条**有效 spacer，
低于预注册阈值 5，`spacer_only` 与 `fusion` 均弃权（属契约行为，不是 bug）。

### 6.4 蛋白侧真实排序失败 — 1 个

`ssa-jim8777`：蛋白把 gold 排到第 39/16384。这是真实失败，不是方向或标签问题。

### 6.5 训练暴露与 gold 不一致

19 个系统蛋白为精确训练暴露却仍有失败排序（如 `nme`、`tde`），
说明这些失败**不能用「模型没见过」解释**，更可能是训练标签与实验 gold 定义不一致。
本轮不做进一步归因（超出范围），如实记录。

## 7. 结果复用与校验

本轮结果**没有**直接沿用旧文档数字，而是重新运行。复用前已核对：

- 清单与全部 **60 个工件 SHA-256 全部通过**（`--check-files`）；
- 模型 `cas9_full`，权重快照 `407f7fc32146a4c0db13c05284f9f3a7cf0ff612`；
- 参数与 [EVALUATION_PROTOCOL.md](EVALUATION_PROTOCOL.md) 第 5 节一致；
- 代码版本见 [SUBMISSION.md](SUBMISSION.md) 的 commit 记录。

旧文档 `docs/MULTI_EVIDENCE_BENCHMARK.md` 的 19 系统表与本轮 22 系统表
**数字不同**，原因是本轮新增 3 个严格蛋白系统并修正了
「无 spacer 通道被当成 0 分计入分母」的口径。旧表保留为历史记录，不删除。

## 8. 复跑命令

见 [EVALUATION_PROTOCOL.md](EVALUATION_PROTOCOL.md) 第 9 节。

## 9. 尚未证明

- 模型泛化能力（严格配对系统 0 个）；
- 融合带来的准确率提升（唯一非零增益来自 1 个系统，CI 含 0）；
- 切割活性、编辑效率、成功概率（本项目**从未**输出这类量）；
- 端到端基因组赛道的性能（合格系统 0 个）。
