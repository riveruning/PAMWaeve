# 可执行流程展示

`01_reproduce_candidates.ipynb` 使用同一套主运行入口，逐步完成输入校验、真实推理/评价、候选导出和清单检查。

先按主 README 配好分析环境、固定版本权重与上游代码。Notebook 本身的执行依赖见 `requirements.txt`；这些是本次验证使用的版本，不是核心 CLI 的必需依赖。

```bash
python -m pip install -r notebooks/requirements.txt
jupyter execute notebooks/01_reproduce_candidates.ipynb \
  --output 01_reproduce_candidates.executed.ipynb
```

也可在 Jupyter 中打开并从头运行。首次运行没有缓存时会执行模型推理；`FRESH=True` 可强制重跑。通过环境变量 `PAMWEAVE_NOTEBOOK_OUTDIR` 可指定已有分析目录，默认 `data/parsed/bench_repro21`。

随包 Notebook 的输出来自本次实际内核运行；其中模型缓存复用明确显示为21个系统。此前同一目录的21系统全新推理记录见 `logs/verification_20261004/run.json`，不将 Notebook 的缓存运行说成第二次全新推理。
