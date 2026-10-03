# Detailed reproduction guide / 详细复现说明

本包是第一届全球大学生生命科学挑战赛「AI 基因编辑与核酸工具设计」赛道的
代码提交包，对应**附件 5《代码提交要求》**。

它做三件事：

1. 给定一条 **Cas9 蛋白序列**，在完整具体候选空间中给出**候选 PAM 兼容度排序**；
2. 叠加**第二条独立证据**（spacer → protospacer 定向侧翼），
   并**只做透明后融合**，输出一致／冲突状态，**不伪造联合概率**；
3. 在**按 Cas9 系统聚合**的评价框架下报告 MRR / Recall@k / 覆盖率与弃权，
   并对每个赛道**分别审计严格独立资格**。

> **一句话边界**：本项目输出的是**证据兼容度与排序值**，
> **不是**切割活性、编辑效率或实验成功概率。
> 本轮**没有**完成湿实验，严格独立样本不足，**不宣称泛化或融合提升**。

---

## 0. 先读这三份文件

| 文件 | 内容 |
|---|---|
| [docs/EVALUATION_PROTOCOL.md](EVALUATION_PROTOCOL.md) | **先于跑分冻结**的评价方案：样本纳排、候选空间、参数、指标、并列与弃权规则 |
| [docs/RESULTS_ANALYSIS.md](RESULTS_ANALYSIS.md) | 实际结果、失败与冲突分析、**为什么说「融合改善」尚无证据** |
| [docs/CASE_SELECTION.md](CASE_SELECTION.md) | 最终参赛案例与选择理由、`results.csv` 字段说明 |

---

## 1. 环境配置

已在以下环境验证：

| 项 | 值 |
|---|---|
| 操作系统 | Ubuntu 24.04.4 LTS（WSL2，内核 `6.6.87.2-microsoft-standard-WSL2`） |
| Python | CPython **3.10.20**（conda-forge） |
| CPU | 22 逻辑核 |
| 内存 | 23 GiB |
| GPU | **未使用**（本轮全部 CPU 推理；CUDA 路径未测） |
| 依赖 | 见 [requirements.txt](../requirements.txt)（固定版本） |

```bash
# 1) 建立环境（示例用 venv；conda 亦可）
python3.10 -m venv .venv
source .venv/bin/activate

# 2) 安装固定版本依赖
python -m pip install -r requirements.txt

# 3) 获取模型权重与上游代码（非商业许可，不随包分发）
#    完整命令与校验值见 models/MODEL_CARD.md 第 1 节
export HF_HOME=$PWD/data/checkpoints/hf
python -m pip install "huggingface_hub==0.36.2"
python - <<'PY'
from huggingface_hub import snapshot_download
snapshot_download(
    repo_id="Profluent-Bio/protein2pam-cas9_full",
    revision="407f7fc32146a4c0db13c05284f9f3a7cf0ff612",
)
PY

# 4) 获取上游推理适配层（非商业许可，不随包分发）
git clone https://github.com/Profluent-Bio/Protein2PAM .reference/Protein2PAM
git -C .reference/Protein2PAM checkout 887026ce058d6a53b04b3c9205990a65d5a92306

# 5) 让本包代码可导入
export PYTHONPATH=src:$PYTHONPATH
```

**上游适配层放在哪里？** 按以下顺序自动查找，**不写死开发机路径**：

1. 环境变量 `PAMPRIDICT_PROTEIN2PAM_HOME`（指向 clone 出来的仓库根目录，
   或直接指向其中的 `protein2pam` 包目录）；
2. `PYTHONPATH` 上任何含 `protein2pam/` 的目录；
3. 包内约定位置 `<package-root>/.reference/Protein2PAM`。

若三者都没有，推理会**明确报错并列出这三种提供方式**，而不是抛出难懂的
`FileNotFoundError`。

如果上游树不在包内（例如放在 `/opt/Protein2PAM`），这样即可：

```bash
export PAMPRIDICT_PROTEIN2PAM_HOME=/opt/Protein2PAM
export PYTHONPATH=src:$PYTHONPATH
```


**预计耗时与资源**：

| 步骤 | 时间 | 资源 |
|---|---|---|
| 不加载模型的清单校验 | < 5 秒 | CPU |
| 22 系统完整评价（CPU，含真实推理） | **约 25–40 分钟** | 单核为主，< 4 GiB 内存 |
| 仅重算指标（`--reuse-existing`） | < 30 秒 | CPU |
| 生成 `results.csv` | < 10 秒 | CPU |

CPU 推理是主要瓶颈：每个系统都要加载一次 2.5 GiB 权重（**每个系统一个子进程**）。

---

## 2. 主运行入口：一条命令生成 `results.csv`

### 2.0 先分清两个清单（**重要**）

本包带**两份**清单，它们用途不同，**默认一律使用可再分发的那份**：

| 清单 | 系统数 | 工件 | 用途 |
|---|---:|---|---|
| **`manifest_redistributable21.json`** | **21** | 57/57 齐全 | **默认**。所有输入都随包分发，开箱即可校验通过并复跑 |
| `manifest_final22.json` | 22 | 57/60（缺 3） | 仅供追溯：它是本轮**完整参考结果**所用的清单 |

`manifest_final22.json` 比默认清单多出 **`spcas9-pampredict-example`**
（唯一 `paired_evidence` 系统）。它的 3 个输入
（`cas9_full.fasta`、`spacers.fna`、`phages.fna`）为上游
**CC-BY-NC-ND-4.0（禁止演绎）**，**按许可不随包分发**，
因此对该清单跑 `--check-files` 会以**非零状态退出**，这是**预期行为**
而不是缺陷。获取方式见 [docs/THIRD_PARTY.md](THIRD_PARTY.md) 第 2 节。

> **复跑得到的是 21 系统子集，不是完整的 22 系统参考结果。**
> 二者的区别是精确的：少的恰好是 SpCas9 这 1 个系统，
> 其余 **21 个系统的指标与参考结果逐一相同**
> （第 4 节验证记录中已实测：0 项指标差异）。
> 提交的 `results/results.csv` 与 `results/benchmark_report.*` 来自
> **完整 22 系统参考运行**，其中 SpCas9 那 5 行在复跑时会被跳过——
> 这是**数据许可导致的可复现性边界**，已如实标注，不掩饰。

### 2.1 第一步：校验输入（**不需要模型**）

```bash
cd <package-root>
python3 scripts/validate_benchmark_manifest.py \
  --manifest benchmarks/multi_evidence_v2/manifest_redistributable21.json \
  --check-files
```

预期：`"valid": true`，`checked_artifacts: 57`，`failed_artifacts: 0`，
**退出码 0**。`--check-files` 会逐一重算蛋白 / spacer / flank 文件的 SHA-256。

（如需确认完整清单的缺失项，可对 `manifest_final22.json` 跑同一命令：
它会报 60 个工件中 3 个缺失，且**全部**是上述 SpCas9 输入，退出码非 0。）

### 2.2 第二步：运行评价（**真实模型推理**）

```bash
export PYTHONPATH=src:$PYTHONPATH   # 上游适配层见第 1 节
export HF_HOME=$PWD/data/checkpoints/hf

python -u scripts/run_multi_evidence_benchmark.py \
  --manifest benchmarks/multi_evidence_v2/manifest_redistributable21.json \
  --device cpu \
  --fusion-policy abstain \
  --bootstrap 2000 \
  --outdir data/parsed/bench_repro21
```

产物在 `data/parsed/bench_repro21/`：
`benchmark_report.{json,md}`、`benchmark_systems.tsv`、`manifest_audit.json`、
`<system_id>/evaluation.json`、`<system_id>/{protein,spacer,fusion}_scores.tsv`。

### 2.3 第三步：生成候选清单

```bash
python3 scripts/build_competition_results.py \
  --benchmark-dir data/parsed/bench_repro21 \
  --manifest benchmarks/multi_evidence_v2/manifest_redistributable21.json \
  --display-candidates results/display_candidates.json \
  --top-n 5 \
  --systems sp7f7-published-flanks,cj4-campylobacter-jejuni-414-protein-only,asp-acidiphilium-21-60-14-protein-only,cba-caulobacterales-protein-only \
  --out results/results_repro21.csv
```

这样得到的是 **20 行**（4 个系统 × 5 个候选），对应 21 系统子集中
本包可完整复跑的那 4 个冻结案例。
包内随附的 `results/results.csv` 是 **25 行**，含 SpCas9 的 5 行。

`--top-n 5` 是一个**对所有系统统一施加的固定展示规则**
（列出蛋白侧前 5 名，并并入显示清单），不是逐系统挑好看的结果。

### 2.4 一键复跑

```bash
bash run_all.sh
```

`run_all.sh` 默认使用 `manifest_redistributable21.json`；
仅当检测到 PAMpredict SpCas9 输入确实存在时，才切换到完整 22 系统清单。

### 2.5 测试

```bash
PYTHONPATH=src:$PYTHONPATH python3 scripts/run_tests.py
```

预期 **168 passed, 0 failed，退出码 0**，不需要权重、网络或 CRISPRCasTyper。
这是**提交包专用**入口：主仓库的 `scripts/run_tests.py`（504 项）还包含
本地网页界面、智能助手与复现底座的用例，那些模块**不在本包内**，
因此本入口只列出**确实随包分发**的模块，并显式打印被排除的模块清单。

---

## 3. `results.csv` 字段

UTF-8、逗号分隔，1 行表头 + 25 行候选（覆盖 5 个冻结案例），24 列。

> **`competition_track` 与 `track` 是两个不同的东西，不要混用。**
> `competition_track` 是**比赛赛道**；`track` 是本项目内部的**证据类型**标签。

| 字段 | 含义 |
|---|---|
| `candidate_id` | 候选编号 `<system_id>-<序号>` |
| **`competition_track`** | **比赛赛道**：`赛道二：AI 基因编辑与核酸工具设计` |
| `competition_track_id` | 该赛道的机器可读标识 `track-2-ai-gene-editing-nucleic-acid-tools` |
| `track` | **证据类型**（不是比赛赛道）：`paired_evidence` / `published_flanks` / `protein_only` |
| `track_definition` | 上一列取值的文字解释 |
| `cas9_id` | Cas9 标识（如 `SpCas9`、`EFC33367.1`） |
| `protein_sequence` | 蛋白氨基酸序列 |
| `protein_sequence_sha256` | 蛋白序列校验值 |
| `candidate_pam` | 候选 PAM（IUPAC） |
| `pam_length` | PAM 长度 |
| `protein_compatibility_score` | 蛋白侧**兼容度**（**非**活性／效率／成功概率） |
| `protein_rank_in_length` | **同长度内**蛋白侧名次 |
| `spacer_evidence_score` | spacer 侧证据分（不适用则**留空**） |
| `spacer_rank_in_length` | **同长度内** spacer 侧名次（不适用则留空） |
| `support_spacer_count` | 支持该候选的独立 spacer 数 |
| `support_target_count` | 支持该候选的 target contig 数；**聚合侧翼证据此列留空** |
| `evidence_status` | 证据状态标签（一致／冲突／不足）；**不是综合分** |
| `ranking_basis` | 排序依据或弃权原因（含 `fusion_gate`） |
| `model_version` | 使用的模型（`cas9_full`） |
| `code_version` | **实际产生该分析的源码版本**，形如 `<git-HEAD>+src.<源码哈希前12位>` |
| `code_version_source` | 上一列是怎么得到的（Git HEAD？源码哈希？） |
| `code_source_sha256` | 对分析源码集合计算的 SHA-256（**解压后没有 Git 也保留**） |
| `evidence_source` | `raw_spacer_target_match` 或 `published_aggregate_flank` |
| `notes` | 该行证据层级的限制说明 |

### 代码版本是怎么记录的

`code_version` 记录的是**真正产生该 CSV 的源码**，而不是某个固定字符串：

- 本机在 Git 检出内运行时：`<commit>+src.<12 位源码哈希>`，
  这样**即使 commit 之后又改了代码也能看出来**；
- 包被解压到**没有 Git** 的位置时：退化为 `src.<12 位源码哈希>`，
  并明确写「no Git checkout enclosing this package, so no commit is claimed」——
  **不冒用一个无关仓库的 commit**；
- 哈希覆盖 11 个分析源文件（运行器、评分器、融合、benchmark 模块、
  推理适配层），清单见生成脚本的 `CODE_IDENTITY_FILES`，
  构建时会打印在 `code_identity_files` 里。

一致性校验示例：

```bash
# CSV 里记录的源码哈希，应与对同一组文件重新计算的结果一致
python3 - <<'PY'
import csv, hashlib, pathlib
files = ["scripts/run_multi_evidence_benchmark.py", "scripts/build_competition_results.py",
         "scripts/score_candidate_pams.py", "scripts/score_spacer_pams.py",
         "scripts/score_published_flanks.py", "scripts/fuse_pam_evidence.py",
         "src/pamdict/benchmark/multi_evidence.py", "src/pamdict/score/candidate.py",
         "src/pamdict/score/spacer.py", "src/pamdict/score/fusion.py",
         "src/pamdict/infer/p2pam.py"]
d = hashlib.sha256()
for rel in sorted(files):
    d.update(rel.encode()); d.update(b"\0")
    d.update(pathlib.Path(rel).read_bytes()); d.update(b"\0")
recorded = next(csv.DictReader(open("results/results.csv")))["code_source_sha256"]
print("recomputed:", d.hexdigest())
print("recorded  :", recorded)
print("MATCH" if d.hexdigest() == recorded else "MISMATCH")
PY
```

### 已强制的规则

- **缺失一律留空，绝不写 0**（`support_target_count` 对聚合侧翼就是留空，
  因为原始 target 标识**从未发表**，写 0 会变成一句假话）；
- 两类来源标签**明确区分**：文献聚合侧翼**不会**被写成「重新完成了原始
  spacer 匹配」；
- 兼容度与融合排序值**不改名**成切割活性／编辑效率／成功概率；
- **不跨蛋白、不跨 PAM 长度**比较名次；
- 两个证据通道**没有校准总分**，因此**不拼接综合分**；
- 比赛赛道与内部证据类型**分列两个字段**，不互相顶替。

### 赛事模板映射

当前 `request/` 附件**未提供**完整候选清单模板，
因此本包按附件 5「未另行规定时默认 UTF-8 `results.csv`」输出，
并加上 `competition_track` 字段显式声明所属赛道。
**未声称满足任何尚未发布的字段规范。**

---

## 4. 最终案例

| 案例 | 层级 | 说明 |
|---|---|---|
| `spcas9-pampredict-example` | **工程展示** | 唯一有原始 spacer＋target 的系统；**蛋白为精确训练暴露** |
| `sp7f7-published-flanks` | 回顾性生物学 | 已入包、CC-BY-4.0，无 `data/` 也能跑；聚合侧翼证据 |
| `cj4-*`、`asp-*`、`cba-*` | **严格独立**（`protein_only`） | 仅有的 3 个合格样本，**差结果也全部保留** |

理由与局限见 [docs/CASE_SELECTION.md](CASE_SELECTION.md)。

---

## 5. 实际结果摘要（22 系统，真实 CPU 推理）

| 方法 | 覆盖率 | MRR (all) | Recall@1 | 适用系统 |
|---|---:|---:|---:|---:|
| `protein_only` | 1.000 | 0.683 | 0.636 | 22 |
| `spacer_only` | 0.947 | 0.700 | 0.579 | 19 |
| `fusion`（`abstain`） | 0.474 | 0.447 | 0.421 | 19 |

- 严格独立系统 **3 个，全部 `protein_only`**；严格独立 `paired_evidence` **0 个**；
- 严格独立蛋白单路 MRR **0.008 / 0.036 / 0.067**（n=3）——**诚实的负结果**；
- **没有证据支持「融合改善」**：配对比较 n=9，平均 MRR 差 +0.108，
  但**中位数为 0.000**，且 +0.108 **完全来自单个系统**
  （其余 8 个恰为 0），bootstrap 区间含 0；
- 门控在 10 个系统上弃权，其中 **7 个**蛋白单路本来就把 gold 排第 1 ——
  门控的作用是**避免错误合并**，不是制造提升。

完整分析见 [docs/RESULTS_ANALYSIS.md](RESULTS_ANALYSIS.md)。

---

## 6. 目录结构

```text
.
├── README.md                     # 本文件
├── requirements.txt              # 固定版本依赖
├── run_all.sh                    # 一键复跑
├── src/pamdict/                  # 核心源代码
│   ├── score/                    #   候选兼容度、spacer 证据、透明后融合
│   ├── benchmark/                #   系统级评价、指标、契约审计
│   ├── infer/                    #   Protein2PAM 推理适配层
│   └── collect/ finetune/ ...    #   数据收集与近邻诊断
├── scripts/                      # 运行入口
│   ├── run_multi_evidence_benchmark.py   # 主评价入口
│   ├── build_competition_results.py      # 生成 results.csv
│   ├── validate_benchmark_manifest.py    # 输入校验（不需模型）
│   └── score_candidate_pams.py 等        # 单步工具
├── benchmarks/
│   ├── multi_evidence_v2/        # 冻结清单 + 系统输入 + 审计表
│   └── dual_evidence_cases/      # 案例定义
├── models/MODEL_CARD.md          # 模型说明（权重版本、许可、局限）
├── results/
│   ├── results.csv               # ★ 最终候选清单
│   ├── benchmark_report.{md,json}
│   └── benchmark_systems.tsv
├── docs/
│   ├── EVALUATION_PROTOCOL.md    # 冻结的评价方案
│   ├── RESULTS_ANALYSIS.md       # 结果与限制
│   ├── CASE_SELECTION.md         # 案例选择
│   ├── THIRD_PARTY.md            # 数据与第三方来源、许可
│   └── SUBMISSION.md             # 提交包说明与验证记录
└── logs/                         # 运行日志
```

---

## 7. 可选：基因组注释（**本提交包未使用**）

`src/pamdict/genome.py` 与 `scripts/run_genome_agent.py` 支持从组装基因组 FASTA
出发做 CRISPR-Cas 注释 → Cas9 提取 → PAM 报告，需要额外安装
**CRISPRCasTyper 1.9.0**。

该赛道（`end_to_end_genome`）本轮**合格评价系统为 0**：
流程已在真实基因组上跑通，但**没有冻结的实验 gold 谱与之配对**，
因此**不进入排名评价**，只作工程演示。

---

## 8. 测试

见第 2.5 节。摘要：

```bash
# 离线回归（不需要权重、网络或 CRISPRCasTyper）
export PYTHONPATH=src:$PYTHONPATH   # 上游适配层见第 1 节
python3 scripts/run_tests.py
```

本包实测 **168 passed, 0 failed，退出码 0**（隔离解压目录中同样通过）。
主仓库的合并回归为 **504 项**，额外包含本地网页界面、智能助手与复现底座的
用例；那些模块**未进入本包**，本入口会显式打印被排除的模块清单，
不用「跳过」来掩盖差异，也**不把两处计数相互替代**。

---

## 9. 引用与许可

- 模型权重：`Profluent-Bio/protein2pam-cas9_full`，**CC-BY-NC-4.0**；
- 上游代码：`Profluent-Bio/Protein2PAM`，**PolyForm Noncommercial 1.0.0**；
- 文献侧翼数据：CRISPRCasDB / Vink et al.，**CC-BY-4.0**；
- 详见 [docs/THIRD_PARTY.md](THIRD_PARTY.md) 与
  [models/MODEL_CARD.md](../models/MODEL_CARD.md)。

本仓库原创代码采用 MIT 许可；第三方数据、模型及上游软件不受本仓库 MIT 许可覆盖。
第三方数据保留来源与原许可证。
