# 数据入口与获取说明

可再分发输入位于 `benchmarks/multi_evidence_v2/systems/`，不重复复制到本目录。默认清单 `manifest_redistributable21.json` 包含 21 个系统、57 个工件；路径、用途与 SHA-256 均写在清单中。

| 数据 | 包内位置或获取方式 | 用途与来源 |
|---|---|---|
| Cas9 蛋白 | `systems/*/protein.faa`；`examples/example_protein.faa` | 模型输入；NCBI 或文献补充材料，来源见清单 |
| 聚合侧翼 | `systems/*-published-flanks/published_flanks.tsv` | 18 个系统，CRISPRCasDB / Vink et al.，CC-BY-4.0 |
| spacer | `systems/*-published-flanks/spacers.fna` | 来源关联；聚合侧翼没有原始 target 标识 |
| SpCas9 配对输入 | 按 `docs/THIRD_PARTY.md` 获取 PAMpredict 示例，并提供完整清单登记的蛋白输入 | 原始配对工程案例；受限输入不随包分发 |
| 权重与适配层 | `models/MODEL_CARD.md` 的固定版本获取命令 | Protein2PAM 预训练模型；分别保留模型与上游软件许可 |

完整 22 系统清单 `manifest_final22.json` 用于追溯参考结果。缺少的三个 SpCas9 输入不视为已经分发的数据；默认 21 系统运行生成 20 行候选，完整参考 CSV 有 25 行。

```bash
python scripts/validate_benchmark_manifest.py \
  --manifest benchmarks/multi_evidence_v2/manifest_redistributable21.json \
  --check-files
```

`data/checkpoints/` 是默认权重缓存；`data/parsed/` 存放运行生成的模型输出与评价报告。这些目录由运行过程创建，不进入源码包。许可、数据来源和训练暴露审计见 `docs/THIRD_PARTY.md`、`docs/EVALUATION_PROTOCOL.md` 与 `results/manifest_audit.json`。
