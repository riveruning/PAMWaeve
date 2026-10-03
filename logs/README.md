# 运行日志与复现记录

`bash run_all.sh` 自动将实际执行记录写入 `logs/runs/<UTC时间>/`：

- `01_validate.log`：输入校验 stdout/stderr。
- `02_benchmark.log`：模型推理与系统评价 stdout/stderr。
- `03_candidates.log`：候选清单生成 stdout/stderr。
- `run.json`：命令、Python/依赖/系统、设备、参数、种子、源码与产物哈希。

每个阶段记录实际退出码和耗时。普通运行日志不自动提交 Git；经检查的 `verification_20261004/` 验证记录随仓库与压缩包分发。
日志仅将本机绝对路径替换为 `<PACKAGE_ROOT>`、`<MODEL_CACHE>` 等，不改动数值和执行结论。

Bootstrap 使用分析函数的实际默认种子 **17**，重复 **2000** 次。入口对子进程设置 `PYTHONHASHSEED=17`。预训练模型在 `eval()` 模式推理；没有声称设置新的 PyTorch 推理种子，也没有新增训练或微调。

本次验证是开发机本地复跑，不等于独立新机安装、GPU 验证或湿实验。历史完整 22 系统结果与本次 21 系统验证分别标注。
