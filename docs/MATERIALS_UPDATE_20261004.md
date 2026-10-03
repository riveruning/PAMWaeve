# 2026-10-04 交付材料补充

本次保留核心评分、模型权重和基准评价逻辑，补充可复现运行的交付材料。

## 文件与入口

- `data/README.md`：说明 `benchmarks/` 中实际分发的数据、来源、许可及受限输入获取方式。
- `scripts/run_submission.py`：执行输入校验、真实推理/评价与候选导出，保存实际执行记录。
- `run_all.sh`：保留原主入口命令，增加 `--fresh`、`--outdir`、`--logdir` 参数。
- `logs/README.md`：解释日志与种子字段；普通用户运行写入不自动提交的 `logs/runs/`。
- `logs/verification_20261004/`：分发本次实际本地验证的输出与元数据。
- `notebooks/01_reproduce_candidates.ipynb`：逐步执行同一入口并检查候选，附实际执行输出与 Notebook 环境依赖。

## 目录与必备材料对应

| 附件5板块 | 本项目提供方式 |
|---|---|
| 项目、环境、输入输出说明 | `README.md`、`docs/REPRODUCTION.md`、固定版本 `requirements.txt` |
| 数据与许可 | `data/README.md`、`benchmarks/`、`docs/THIRD_PARTY.md` |
| 核心源码 | `src/pamdict/` |
| 最终模型 | `models/MODEL_CARD.md`：固定版本、获取步骤、校验值和局限 |
| 主预测/筛选入口 | `bash run_all.sh`，等效入口名称按附件5正文允许 |
| Notebook | `notebooks/01_reproduce_candidates.ipynb`，实际执行同一入口并检查候选 |
| 训练入口 | 最终结果采用预训练模型；复现无需训练。未采用的历史微调探索在 Model Card 中披露，相关探索代码保留 |
| 候选清单 | 完整22系统参考 `results/results.csv`；本次21系统复跑 `results/results_repro21.csv` |
| 日志、参数、种子 | `logs/verification_20261004/` 及运行时生成的 `run.json` |

## 验证边界

本次验证在开发机执行，不能解释为独立新机、GPU 或湿实验验证。
可再分发队列为21系统；完整22系统参考包含额外的受限 SpCas9 案例。
验证完成后，以日志中的退出码、实际耗时及交付检查记录为准。

## 已完成验证

- 21系统全新CPU模型推理、输入校验和候选导出全部成功，主流程耗时约263秒。
- 57个可分发输入工件全部通过哈希校验，生成20行候选。
- 与完整22系统参考结果中的共享系统比较：21个系统的评价指标、20行候选的科学字段均无差异。Git运行版本字段随实际检出不同，单独标注。
- 168项离线回归全部通过。
- 故意缺少上游适配层的隔离测试确认：错误退出码被保留，运行记录标为失败，候选导出不会继续执行。
- 模型权重、配置的实际SHA-256与Model Card的固定版本一致。
- Notebook使用真实Jupyter内核执行四个代码单元，复用前述21系统的真实模型输出；不宣称第二次全新推理。

上述原始输出与检查结果在 `logs/verification_20261004/`，其中 `run.json` 是全新推理记录，`result_comparison.json` 是科学结果逐项比较记录。
