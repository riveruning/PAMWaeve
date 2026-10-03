# 提交包说明与验证记录

状态日期：**2026-10-03**

## 1. 这个包是什么

对应**附件 5《代码提交要求》**的代码提交包，赛道为
**二、AI 基因编辑与核酸工具设计**。

包根目录：`deliverables/PAMPRIDICT_submission/`

### 收录了什么

| 类别 | 内容 |
|---|---|
| 核心源代码 | `src/pamdict/{score,benchmark,infer,collect,finetune,schema,genome_tests}` |
| 主运行入口 | `scripts/run_multi_evidence_benchmark.py`（主评价）、`scripts/build_competition_results.py`（出候选清单）、`scripts/validate_benchmark_manifest.py`（输入校验） |
| 单步工具 | 蛋白/spacer/侧翼评分与融合：`score_candidate_pams.py`、`score_spacer_pams.py`、`score_published_flanks.py`、`fuse_pam_evidence.py` |
| 案例与审计脚本 | `run_dual_evidence_case.py`、`verify_dual_evidence_value.py`、`audit_*`、`build_pammla_benchmark.py` 等——**包内测试会直接调用它们**，故必须随包 |
| 测试入口 | `scripts/run_tests.py`（**提交包专用**，168 项，见第 3.2 节第 8 项） |
| 依赖清单 | `requirements.txt`（固定版本） |
| 冻结清单与输入 | `benchmarks/multi_evidence_v2/`（21 系统可再分发清单 + 22 系统全量清单 + 各系统输入） |
| 案例定义 | `benchmarks/dual_evidence_cases/` |
| 示例输入 | `examples/example_protein.faa` |
| 结果样例 | `results/results.csv`、`benchmark_report.{md,json}`、`benchmark_systems.tsv` |
| 模型说明 | `models/MODEL_CARD.md` |
| 说明文档 | `docs/{EVALUATION_PROTOCOL,RESULTS_ANALYSIS,CASE_SELECTION,THIRD_PARTY,SUBMISSION}.md` |
| 一键复跑 | `run_all.sh`（默认 21 系统可再分发清单） |

### 明确排除了什么

| 排除项 | 原因 |
|---|---|
| 本地网页界面（`pamdict/app`）、智能助手（`pamdict/agent`） | **本轮不新增网页功能**；且核心结果**无需大模型 API** |
| 复现底座（`pamdict/repro`）及应用层测试 | 与核心分析无关，属调试／交付脚手架 |
| 模型权重（2.5 GiB） | 非商业许可，只给固定版本获取步骤 |
| Protein2PAM 上游代码 | PolyForm Noncommercial 1.0.0，只给固定 commit |
| PAMpredict 示例（SpCas9 paired 输入） | **CC-BY-NC-ND-4.0（禁止演绎）**，不可再分发 |
| `data/` 全部本地产物、缓存、用户样本 | 混合来源，可再分发性**待核实** |
| 密钥、私有配置、环境缓存 | 见第 4 节检查 |
| 既有 PPT、`request/` 比赛附件、历史展示材料 | 与本次代码提交无关 |
| 中间试跑清单（`manifest_stage_*.json`） | 系统数与命名易误导；已删除，只保留两份正式清单 |

## 2. 最终模型与历史探索的区分

| 项 | 结论 |
|---|---|
| **最终采用权重** | `Profluent-Bio/protein2pam-cas9_full`，rev `407f7fc32146a4c0db13c05284f9f3a7cf0ff612` |
| 是否使用本项目微调权重 | **否**。全部结果来自干净预训练权重 |
| 历史微调（F0-step20 等） | **未采用**，项目自查判定 `PROMISING_BUT_NOT_ACCEPTED`，**未包装成模型改进** |

详见 [MODEL_CARD.md](../models/MODEL_CARD.md) 第 2 节。

## 3. 验证记录（**本机验证，不是独立新机验收**）

### 3.1 验证环境

| 项 | 值 |
|---|---|
| 机器 | 本开发机（WSL2），**非独立新机** |
| 系统 | Ubuntu 24.04.4 LTS，内核 `6.6.87.2-microsoft-standard-WSL2` |
| Python | CPython 3.10.20（conda-forge） |
| 设备 | CPU（`--device cpu`），无 GPU |

### 3.2 已执行的验证

| # | 验证项 | 命令 | 结果 |
|---|---|---|---|
| 1 | **默认清单**（README 第 2 节） | `validate_benchmark_manifest.py --manifest manifest_redistributable21.json --check-files` | `valid=true`，21 系统，**57 工件，0 失败，退出码 0** |
| 2 | 全量清单（追溯用） | 同上换 `manifest_final22.json` | 22 系统；60 工件中 **57 通过，3 未通过**，退出码非 0 |
| 3 | 未通过项定位 | 同上 | 3 项**全部**是 SpCas9 的 `cas9_full.fasta` / `spacers.fna` / `phages.fna`——**按 CC-BY-NC-ND-4.0 刻意不随包分发**，属预期 |
| 4 | **`run_all.sh` 默认流程** | 隔离副本中 `bash run_all.sh` | 三步全过；产 **20 行** `results_repro21.csv`；退出码 0 |
| 5 | 隔离目录完整评价 | 隔离副本跑 21 系统真实推理 | **21/21 完成**；与参考运行 **0 项指标差异** |
| 6 | 隔离目录候选清单 | 隔离副本重新生成 | 20 个共享候选 **0 项字段差异** |
| 7 | 压缩包解压后复跑 | 解压 zip → 跑 21 系统 | **21/21 完成**；**0 项指标差异** |
| 8 | **包内测试入口** | `PYTHONPATH=src:$PYTHONPATH python3 scripts/run_tests.py` | **168 passed, 0 failed，退出码 0**（隔离解压目录中同样通过） |
| 9 | 主仓库合并回归 | 主仓库 `scripts/run_tests.py` | **504 passed, 0 failed** |
| 10 | `code_version` 无 Git 行为 | 在无 `.git` 的副本中生成清单 | 记为 `src.6f3d90da486c` 并声明「no commit is claimed」；**源码哈希与含 Git 时一致** |
| 11 | 结果不变量 | 脚本校验 `results.csv` | **11/11 通过**，见 3.4 |

> 第 9 项的 504 = 基线 499 + 本轮新增 5 项回归
> （蛋白单路不计入分母、配对比较剔除不适用、`n/a` 渲染、`rel()` 路径、
> 聚合侧翼计数不得读成 0 target）。第 8 项的 168 是**提交包内**可运行的
> 子集：网页界面、智能助手与复现底座的模块未随包分发，
> 包内入口会**显式打印被排除的模块名**，不用「跳过」掩盖差异，
> 也不把 168 与 504 相互替代。
>
> 一次中间运行曾出现 `test_assistant` 1 项失败；**经基线对照确认为既有偶发**
> （该用例不导入任何本轮改动模块），随后完整重跑为 504/504 全通过。
> 此处如实记录，不隐藏。


### 3.3 隔离目录端到端运行

在**只含本包文件**的目录（无 `.pylibs`、无 `.reference`、无 `data/`）中执行
主入口，模型权重与上游适配层按 [MODEL_CARD.md](../models/MODEL_CARD.md) 的
说明通过环境变量指向包外位置——这正是文档所声明的部署方式。

**这一步发现并修复了两个真实缺陷**（此前在开发机内始终被掩盖）：

1. 运行器**硬编码** `PYTHONPATH=.pylibs:.reference/Protein2PAM:...`，
   会覆盖操作者自己的 `PYTHONPATH`；现改为「总是前置包内 `src`，
   仅在本机确实存在时才追加约定目录，否则保留继承的 `PYTHONPATH`」；
2. `pamdict/infer/p2pam.py` **硬编码**上游适配层相对路径；现按
   `PAMPRIDICT_PROTEIN2PAM_HOME` → `sys.path` 上的 `protein2pam` → 包内
   `.reference/Protein2PAM` 顺序查找，找不到时给出**可操作的报错**。

**结论**：入口在包内可运行，**不依赖未说明的开发机绝对路径**
（`run_all.sh` 用 `$(dirname "${BASH_SOURCE[0]}")` 定位包根，
用 `PYTHON`/`HF_HOME`/`PAMPRIDICT_PROTEIN2PAM_HOME` 定位解释器、权重与适配层）。

### 3.4 `results.csv` 不变量校验

| 检查 | 结果 |
|---|---|
| `published_aggregate_flank` 行的 `support_target_count` 为空而非 0 | **通过** |
| `raw_spacer_target_match` 行的 `support_target_count` 有值 | **通过** |
| `protein_only` 行的 spacer 字段为空而非 0 | **通过** |
| 无 `activity` / `efficiency` / `success_prob` / `cleavage` 字段名 | **通过** |
| 5 个冻结案例的 Cas9 标识齐全 | **通过** |
| 每行都有蛋白序列与 SHA-256 | **通过** |
| 每行都有 `model_version` 与 `code_version` | **通过** |
| 名次从 1 开始 | **通过** |

### 3.5 结果与清单一致性

`results/results.csv` 由 `scripts/build_competition_results.py` 直接投影
`data/parsed/bench_final22/` 的冻结产物生成，**不重新排序、不重新融合**。
其 `benchmark_manifest_sha256` 与运行器写入
`benchmark_report.json` 的值一致
（`ab3ea376df2bad64bca44937113b7213c51d0d23e601afbbd862f4538edaf0ec`）。

### 3.6 **未做**的验证（不得当成已完成）

- **独立新机安装：由用户负责，本轮未做**；第 3 节全部结论仅限本机；
- **真实用户试用：未做**；
- GPU / CUDA 推理：**未测**；
- macOS / Windows 原生：**未测**；
- `end_to_end_genome` 赛道的性能评价：**无合格系统**，未做。

## 4. 打包前安全检查（实际执行记录）

扫描命令（在包根目录执行）：

```bash
# 明文密钥／令牌
grep -rInE "sk-[A-Za-z0-9]{20,}|api[_-]?key.{0,4}[:=].{0,4}['\"][A-Za-z0-9]{12,}|AKIA[0-9A-Z]{16}|-----BEGIN [A-Z ]*PRIVATE KEY-----" .
# 开发机绝对路径
grep -rIln "/home/<dev-user>" .
```

| 检查项 | 结果 |
|---|---|
| 明文密钥／令牌 | **0 个真实命中**。唯一的模式匹配是**蛋白氨基酸序列**中偶然出现的 `AKIA` 字样（如 `...AAKIATVSK...`），属假阳性，非 AWS 凭据 |
| 开发机绝对路径 | **0 个**。初次扫描发现 11 个文件带有开发机路径（7 个溯源 summary、1 个测试脚本文档串、1 个报告 JSON、1 个源码注释），**已全部修正**：溯源记录改为包内相对路径，运行器新增 `rel()` 只写仓库相对路径 |
| 用户样本／上传文件 | 未收录（`data/app_jobs*`、`assistant_uploads` 均不在包内） |
| 环境缓存 | 未收录（无 `__pycache__`、`.venv`、`hf` 缓存） |
| 无关材料 | 未收录（无 PPT、无 `request/` 比赛附件、无历史展示 zip） |
| 依赖未声明路径 | 入口只用包内相对路径 + `PYTHON`/`HF_HOME`/`PAMPRIDICT_PROTEIN2PAM_HOME` 环境变量；隔离目录实测通过（3.3 节） |

## 5. 压缩包

| 项 | 值 |
|---|---|
| 路径 | `deliverables/PAMPRIDICT_code_submission.zip` |
| 顶层目录 | `PAMPRIDICT_submission/` |
| 文件数 | 247 |
| 大小 | 约 0.5 MB |

打包命令：

```bash
cd deliverables
zip -rq PAMPRIDICT_code_submission.zip PAMPRIDICT_submission \
  -x '*/__pycache__/*' -x '*.pyc' -x '*/.DS_Store'
```

**解压后复跑已验证**（3.2 节第 8 项）：解压出的目录可直接跑通 21 系统评价，
与参考运行 **0 项指标差异**。

## 6. 复跑

```bash
cd deliverables/PAMPRIDICT_submission
PYTHON=python3 bash run_all.sh
```

或按 [README.md](../README.md) 第 2 节逐步执行。

## 7. 未完成事项（不得当成已完成）

| 事项 | 状态 | 责任 |
|---|---|---|
| **独立新机安装** | **未做** | 由用户执行 |
| **真实用户试用** | **未做** | 由用户执行 |
| GPU / CUDA 推理 | **未测** | — |
| macOS / Windows 原生 | **未测** | — |
| `end_to_end_genome` 赛道性能 | **无合格系统**，未评价 | — |
| 严格独立配对证据 | **0 个**系统（门槛 5），未达成 | — |
| 湿实验验证 | **未做**（本轮范围外） | — |
