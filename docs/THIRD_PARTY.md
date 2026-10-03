# 数据与第三方软件来源、版本、用途与许可

状态日期：**2026-10-03**

**可再分发性是逐项判断的**：材料进入本地 `data/` 不代表可以再次分发。
本包只收录许可清楚、允许再分发的部分；其余只给获取说明。
**许可未知的项目明确列为「待核实」，不默认可再分发。**

## 1. 随本包分发的数据

| 材料 | 位置 | 来源 | 版本／校验 | 许可 | 可再分发 |
|---|---|---|---|---|---|
| 18 个 `published_flanks` 系统输入 | `benchmarks/multi_evidence_v2/systems/*-published-flanks/` | Vink et al., CRISPRCasDB, DOI `10.1186/s13059-021-02495-9`, Additional file 2 | 每个文件 SHA-256 记录在 manifest 与 `source_summary.json` | **CC-BY-4.0** | **是** |
| 3 个严格独立蛋白系统 | `benchmarks/.../cj4-*`, `asp-*`, `cba-*` | NCBI Protein / 出版商补充材料，DOI 见各系统 `provenance` | `protein.faa` SHA-256 见 manifest | Asp/Cba：**CC-BY-4.0**（文章补充材料）；Cj4：NCBI 分子数据 | **是** |
| SpCas9 工程对照输入 | `benchmarks/dual_evidence_cases/`（案例定义） | PAMpredict 示例 | 定义文件 SHA-256 见案例 JSON | 案例定义为本项目原创 | **是** |
| 清单与审计表 | `benchmarks/multi_evidence_v2/*.tsv/json` | 本项目生成 | 含 SHA-256 溯源列 | 本项目 | **是** |
| `results/` 下的结果 | 本项目运行产物 | 本项目 | 见 `benchmark_report.json` | 本项目 | **是** |

## 2. 仅给获取说明、**不随包分发**的材料

| 材料 | 来源 | 固定版本 | 许可 | 不再分发原因 |
|---|---|---|---|---|
| Protein2PAM `cas9_full` 权重 | HF `Profluent-Bio/protein2pam-cas9_full` | rev `407f7fc32146a4c0db13c05284f9f3a7cf0ff612`，权重 SHA-256 `aa1dc5c0…41c2` | **CC-BY-NC-4.0** | 体积 2.5 GiB + 非商业限制 |
| Protein2PAM 上游代码 | `github.com/Profluent-Bio/Protein2PAM` | commit `887026ce058d6a53b04b3c9205990a65d5a92306` | **PolyForm Noncommercial 1.0.0** | 许可限制，需自行 clone |
| Protein2PAM 官方训练 TSV | `https://storage.googleapis.com/protein2pam-x83y9z7q4k/protein2pam_train_seqs.tsv` | 58,356,844 B，md5 `a76d28712a0796a85eb6a53825f9acf3` | **CC BY-NC 4.0** | 非商业限制；仅用于训练暴露审计 |
| PAMpredict 示例（spacers/phages） | `github.com/Matteo-Ciciani/PAMpredict` | git `4ad096ce…8593` | **CC-BY-NC-ND-4.0** | 不随此包分发；获取时遵守原许可 |
| 基因组注释库 | NCBI RefSeq/GenBank 组装 | 见各运行目录 | NCBI 条款 | 体积 |
| 用户提供的 PDF 补充材料 | 见 `data/raw/supp_*` | — | **待核实** | 第三方版权，未确认可再分发 |

## 3. 第三方软件

| 软件 | 版本 | 用途 | 许可 | 是否随包 |
|---|---|---|---|---|
| PyTorch | `2.12.1+cu126` | 模型推理 | BSD-3-Clause | 否（pip 安装） |
| transformers | `4.44.2` | 模型加载 | Apache-2.0 | 否 |
| tokenizers | `0.19.1` | 分词 | Apache-2.0 | 否 |
| safetensors | `0.8.0` | 权重读取 | Apache-2.0 | 否 |
| huggingface_hub | `0.36.2` | 权重获取 | Apache-2.0 | 否 |
| NumPy | `1.23.5` | 数值计算 | BSD-3-Clause | 否 |
| pandas | `2.0.3` | 审计辅助 | BSD-3-Clause | 否 |
| Matplotlib | `3.7.2` | 历史绘图 | PSF-based | 否 |
| Protein2PAM（上游代码） | commit `887026ce…` | 模型定义与推理适配 | **PolyForm Noncommercial 1.0.0** | 否 |
| CRISPRCasTyper | `1.9.0` | **可选**基因组注释（`end_to_end_genome` 赛道） | 见上游 | 否 |

## 4. 上游数据许可与再分发说明（重要）

- 18 个 `published_flanks` 系统派生自 **CC-BY-4.0** 的 CRISPRCasDB 补充表，
  本包**保留来源与许可标注**，符合署名要求；
- **Protein2PAM 相关的一切（权重、代码、训练集）均为非商业许可**，
  本包只提供获取步骤，不分发，并已在
  [MODEL_CARD.md](../models/MODEL_CARD.md) 标注 CC-BY-NC-4.0；
- **PAMpredict 示例为 CC-BY-NC-ND-4.0（限制改编材料分发）**，
  因此 SpCas9 工程案例的 **spacer/target 原始文件不随包分发**，
  只提供案例定义与 clone 步骤；离线复跑请改用 `sp7f7` 案例
  （CC-BY-4.0，已随包）。

## 5. 待核实清单（**不默认可再分发**）

| 材料 | 位置 | 状态 |
|---|---|---|
| 用户提供的补充 PDF（Asp/Cba 蛋白来源） | `data/raw/` | **待核实**：文章为 CC-BY-4.0，但用户提供副本的再分发范围未确认 |
| `data/raw/epmc_cache*`（Europe PMC 缓存） | `data/raw/` | **待核实**：逐篇许可不同 |
| `data/parsed/batch_v*`（文献抽取产物） | `data/parsed/` | **待核实**：混合来源，未逐条确认 |
| `data/corpus/gold_*` | `data/corpus/` | **待核实**：派生自混合来源 |
| 项目整体许可证 | 仓库根 | 原创代码 MIT；第三方材料保留原许可 |

**处置**：以上材料**均未进入本提交包**。

## 6. 复现所需最小步骤

```bash
# 1) 本包内已含的清单与可再分发输入可直接校验
python3 scripts/validate_benchmark_manifest.py \
  --manifest benchmarks/multi_evidence_v2/manifest_redistributable21.json

# 2) 需要真实推理时，另行获取模型与上游代码（见 MODEL_CARD.md 第 1 节）
# 3) 需要 SpCas9 paired_evidence 赛道时，另行 clone PAMpredict 示例
```
