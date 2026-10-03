# 参赛案例选择与候选清单说明

状态日期：**2026-10-03**
对应评测方案：[EVALUATION_PROTOCOL.md](EVALUATION_PROTOCOL.md)
对应结果目录：`data/parsed/bench_final22/`

## 1. 主展示案例：`spcas9-pampredict-example`

**层级标注：工程展示（engineering showcase），存在训练暴露，不是独立评价。**

| 项 | 值 |
|---|---|
| system_id | `spcas9-pampredict-example` |
| Cas9 | SpCas9（*Streptococcus pyogenes*，Type II-A） |
| 蛋白序列 SHA-256 | 清单 `provenance.protein.sha256`，作用域 `selected_fasta_record_sequence` |
| 蛋白来源 | `data/raw/cas9_full.fasta` 中的 SpCas9 记录 |
| spacer | PAMpredict 示例 `Example/spacers.fna`，37 条，sha256 `7f4c3de1…85a1` |
| target | `Example/Phages/phages.fna`，67 contig，sha256 `1bcb4311…b37c` |
| 参考 PAM | `NGG`（工程已知 PAM，`gold_evidence_type=engineering_known_pam`） |
| 赛道 | `paired_evidence`（仓库中唯一有原始 spacer＋target 的系统） |
| 方向 | `spacer_orientation=reverse`（已独立重算验证） |
| 训练暴露 | **`exact`** — 该蛋白与 Protein2PAM 官方训练序列完全相同 |
| 严格独立资格 | **否**（暴露 + target 仅 contig 级） |

### 选择理由

1. **它是唯一同时具备「蛋白序列 + 定向 spacer + 原始 target contig」的系统**，
   另外 18 个序列侧系统只有文献聚合侧翼，**不能**写成重新完成了原始 spacer 匹配，
   证据层级明显更低；
2. 它能同时展示本项目的两条证据通路与后融合链路；
3. 它的方向曾被上游隐式决定，本项目**独立重算**并用负对照证明方向敏感
   （方向错误时 xGG 信号归零、破坏同源时命中为 0），是可信的工程校验。

### 必须一起说明的局限（不做美化）

- 蛋白为 Protein2PAM **精确训练暴露**，因此**不能**据此说模型泛化；
- 37 条 spacer 中 **20 条命中同一 contig**（`uvig_394886`），
  所以「55 个靶标」**不代表** 55 个独立观测；按 contig 折叠后是 **7/15 = 46.7%**；
- 67 个 contig 中 64 个是无宿主标注的 `uvig_`/`ivig_` 记录，来自人肠道病毒组，
  **不是** *S. pyogenes* 的生态位，属技术性重建；
- gold 是物种级经典 NGG；
- **第二路证据没有改变 top-1**：蛋白单路已把 NGG 排第 1，第二路是印证而非提升。

## 2. 补充案例：`sp7f7-published-flanks`

**层级标注：回顾性生物学（retrospective biological），无原始 target。**

| 项 | 值 |
|---|---|
| Cas9 | WP_038431314.1（*S. pyogenes* 7F7，Type II-A） |
| spacer | 8 条（sha256 `32e2c799…5329`） |
| 证据 | `published_flanks.tsv`，8 行→7 行可用（sha256 `a4c8f117…7076`） |
| 来源 | Vink et al., CRISPRCasDB, DOI `10.1186/s13059-021-02495-9`，Additional file 2 |
| 许可 | **CC-BY-4.0**（可再分发，已入 Git） |
| 方向 | 源表逐行声明 `orientation_PAMbased`（7/7），无需推断 |
| 训练暴露 | `near`（最近官方训练蛋白一致性 0.99415205） |

**为什么选它作补充**：它**已入 Git**，因此在没有 `data/` 的机器上仍可复跑，
是唯一「不依赖受限第三方数据」的完整双证据案例；同时它清楚展示了
聚合侧翼的证据层级限制（`support_target_count` 缺失而非 0）。

## 3. 严格独立评价系统（3 个，均为 `protein_only`）

这 3 个是**仅有的**满足冻结独立性定义的样本，全部如实保留，包括差结果：

| system_id | Cas9 | PAM 长度 | 参考 PAM | 蛋白 MRR | 来源 |
|---|---|---:|---|---:|---|
| `cj4-campylobacter-jejuni-414-protein-only` | EFC33367.1 | 6 | `NNNGRY` | **0.036** | DOI `10.1038/s42003-025-09430-9` |
| `asp-acidiphilium-21-60-14-protein-only` | OYV68964.1 | 7 | `NNNNGCA` | **0.008** | DOI `10.1089/crispr.2024.0013` |
| `cba-caulobacterales-protein-only` | MBN8552000.1 | 7 | `NGWNCCA` | **0.067** | 同上 DOI |

**如实标注的限制**：

- Cj4Cas9 最近官方训练近邻一致性 **0.89441624**，只比 0.90 门槛低 0.0056，
  属**门槛边界**样本而非远离分布的 OOD；
- Asp 与 Cba **共享同一论文**，**不是两篇独立论文**，因此
  「3 个严格系统」在论文层面只有 **2 个独立来源**；
- n=3 且其中 2 个同源，**不给 bootstrap 区间、不做泛化结论**；
- 这 3 个**没有** spacer/target 材料，因此不参与 `spacer_only` / `fusion` 评价。

## 4. 工程演示案例（不进入性能评价）

`scripts/run_genome_agent.py` 的基因组→注释→Cas9→PAM 流程已在真实基因组上
跑通（`data/parsed/genome_agent_real_20260928/`），但**没有冻结的实验 gold 谱
与之配对**，因此 `end_to_end_genome` 赛道本轮**合格系统数为 0**。
它只作工程演示，不作排名，也不计入任何指标。

## 5. 候选清单 `results.csv`

### 字段说明

| 字段 | 含义 |
|---|---|
| `candidate_id` | 候选编号，`<system_id>-<run_index>` |
| `track` | 赛道：`paired_evidence` / `published_flanks` / `protein_only` |
| `cas9_id` | Cas9 标识（清单 `protein_record_id`） |
| `protein_sequence` | 蛋白序列（FASTA 单行） |
| `protein_sequence_sha256` | 蛋白序列校验值 |
| `candidate_pam` | 候选 PAM（IUPAC） |
| `pam_length` | PAM 长度 |
| `protein_compatibility_score` | 蛋白侧**兼容度**（非活性） |
| `protein_rank_in_length` | 同长度内蛋白侧名次 |
| `spacer_evidence_score` | spacer 侧证据分（缺失则留空） |
| `spacer_rank_in_length` | 同长度内 spacer 侧名次（缺失则留空） |
| `support_spacer_count` | 支持该候选的独立 spacer 数 |
| `support_target_count` | 支持该候选的 target contig 数；**聚合侧翼证据此列留空** |
| `evidence_status` | 一致／冲突／不足 |
| `ranking_basis` | 排序依据或弃权原因 |
| `model_version` | 使用的模型 |
| `code_version` | 代码 commit |
| `evidence_source` | `raw_spacer_target_match` 或 `published_aggregate_flank` |
| `notes` | 备注 |

### 关键规则（已在实现中强制）

- **缺失一律留空**，不填 0；
- 聚合侧翼证据的 `evidence_source` 明确标为 `published_aggregate_flank`，
  与 `raw_spacer_target_match` **区分开**；
- 兼容度与融合排序值**不**改名成切割活性／编辑效率／成功概率；
- **不跨蛋白、不跨 PAM 长度比较名次**（长度内排名）；
- 本轮两个证据通道**没有校准总分**，因此**不拼接综合分**；
  `evidence_status` 是状态标签，不是综合分。

### 生成命令

```bash
cd <repo>
python3 scripts/build_competition_results.py \
  --benchmark-dir data/parsed/bench_final22 \
  --manifest benchmarks/multi_evidence_v2/manifest_final22.json \
  --out results.csv \
  --display-candidates results/display_candidates.json
```

`--display-candidates` 用小候选清单（`NGG,NAG,NGA,NRG` 等）用于**展示**；
它**不代表**正式排名评价，正式排名始终用完整具体 A/C/G/T 候选空间。

## 6. 许可与再分发

| 材料 | 许可 | 本包是否收录 |
|---|---|---|
| `benchmarks/multi_evidence_v2/**`（18 个侧翼系统 + 3 个 strict 蛋白） | 见各 `source_summary.json`；多数 CC-BY-4.0 | 是 |
| Protein2PAM cas9_full 权重 | CC-BY-NC-4.0 | **否**，给出获取步骤 |
| PAMpredict 示例（spacers/phages） | CC-BY-NC-ND-4.0 | **否**，给出 clone 步骤 |
| Protein2PAM 官方训练 TSV | CC BY-NC 4.0 | **否**，给出 URL |
| 全部 `data/` 产物 | 混合，部分**待核实** | 否 |

许可证未知项在 [THIRD_PARTY.md](THIRD_PARTY.md) 中明确列为**待核实**，
不默认可再分发。
