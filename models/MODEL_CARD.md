# Model Card — Protein2PAM `cas9_full`（本项目最终采用的权重）

状态日期：**2026-10-03**

## 1. 最终采用哪个权重

**本项目最终结果（`results.csv`、`benchmark_report`）全部由
`Profluent-Bio/protein2pam-cas9_full` 的干净预训练权重产生，
没有使用任何本项目微调权重。**

| 项 | 值 |
|---|---|
| 模型名称 | `Profluent-Bio/protein2pam-cas9_full` |
| 类型 | 蛋白语言模型 + PAM 分类头（`EsmForSequenceClassification`） |
| Hub revision（快照） | `407f7fc32146a4c0db13c05284f9f3a7cf0ff612` |
| `model.safetensors` SHA-256 | `aa1dc5c017ddd885bf26c8126675075af0007d1af3ab1a2893e62aaffe0341c2` |
| `model.safetensors` 大小 | 2.5 GiB |
| `config.json` SHA-256 | `7fe91bd35342572b6dcd23f937013ec2603b2e136970f4e3a21a9af5b5c41cb5` |
| 来源 | HuggingFace Hub `Profluent-Bio/protein2pam-cas9_full` |
| 上游代码 pin | `Profluent-Bio/Protein2PAM` commit `887026ce058d6a53b04b3c9205990a65d5a92306` |
| 许可（模型权重） | **CC-BY-NC-4.0**（非商业） |
| 许可（上游代码） | **PolyForm Noncommercial 1.0.0**（非商业） |

### 获取方式（不随包分发）

```bash
export HF_HOME=$PWD/data/checkpoints/hf
python -m pip install "huggingface_hub==0.36.2"
python - <<'PY'
from huggingface_hub import snapshot_download
snapshot_download(
    repo_id="Profluent-Bio/protein2pam-cas9_full",
    revision="407f7fc32146a4c0db13c05284f9f3a7cf0ff612",
    local_dir_use_symlinks=True,
)
PY
```

上游模型适配层（`protein2pam` 包）需单独获取：

```bash
git clone https://github.com/Profluent-Bio/Protein2PAM
cd Protein2PAM && git checkout 887026ce058d6a53b04b3c9205990a65d5a92306
```

下载后请核对 `model.safetensors` 的 SHA-256 与上表一致；
**不一致则本轮结果不可复现**。

## 2. 历史微调探索（**不是**最终模型，未用于任何提交结果）

仓库中确实做过微调探索，如实列出，**但它们不是本项目的模型改进**：

| 产物 | 位置 | 状态 |
|---|---|---|
| F0-step20（Cas9 均衡重训） | `data/checkpoints/finetune_v3_cas9_F0_equal/` | **未采用** |
| `finetune_pilot` | `data/checkpoints/finetune_pilot/` | **未采用** |
| `finetune_v2_mixedfamily_invalid` | `data/checkpoints/finetune_v2_mixedfamily_invalid/` | **已作废**（`invalid`） |
| `finetune_v3_cas9_smoke*` | `data/checkpoints/` | 冒烟测试，**未采用** |

项目自身文档（`docs/SYSTEM_SPECTRUM_EVAL.md`）对 F0 的判定为
**`PROMISING_BUT_NOT_ACCEPTED`**，且指出 48 个可用于概率 logo 比较的
system **全部来自同一 DOI**，属同来源校准信号，不能解释为跨论文泛化。

**因此**：

- 本包**不包含**任何微调权重；
- 本轮**没有**做「微调提高」的主张；
- 若未来要采用微调权重，必须重新冻结 manifest 并重跑全部评价。

## 3. 输入与输出

- **输入**：单条 Cas9 蛋白氨基酸序列（FASTA 或字面序列）；
- **输出**：10×4 的 A/C/G/T 概率矩阵（每条蛋白一行），
  按 PAM 长度逐位给出碱基概率；
- 本项目在其上派生：
  - `specificity_adjusted_score`（蛋白侧兼容度）：
    IUPAC 允许概率质量相对均匀四碱基背景的富集；
  - `predicted_pam`：按信息量阈值 `CONSENSUS_THRESHOLD_BITS` 取共识。

## 4. 本项目的实际贡献（相对上游）

1. **不改模型权重**，而是在其输出之上建立
   「候选 PAM 兼容度排序」的可审计口径（`pamdict/score/candidate.py`）；
2. **第二条独立证据通路**：spacer → protospacer 定向侧翼搜索
   （`pamdict/score/spacer.py`），并支持文献聚合侧翼输入；
3. **透明后融合**：只输出一致／冲突状态与审查优先级，
   **不训练联合概率**（`pamdict/score/fusion.py`）；
4. **系统级评价框架**：按 Cas9 system 聚合、完整具体候选空间、
   并列保守处理、弃权门控与逐赛道独立性审计
   （`pamdict/benchmark/multi_evidence.py`、`scripts/run_multi_evidence_benchmark.py`）；
5. **独立性审计**：显式区分训练暴露、近邻与严格独立，并拒绝在样本不足时
   输出泛化结论。

## 5. 已知局限（必须与任何结果一起读）

- **分数不是活性**：输出是**证据兼容度与排序值**，
  **不是**切割活性、编辑效率或实验成功概率；
- **训练暴露普遍**：22 个评价系统中 19 个蛋白与官方训练序列完全相同或
  ≥90% 近邻，基础模型的既有暴露**不能**通过后续排除消除；
- **严格独立样本不足**：仅 3 个，且其中 2 个共享同一论文；
- **未见泛化证据**：严格独立蛋白单路 MRR 0.008–0.067（n=3，
  其中一个还是 0.8944 一致性的门槛边界样本）；
- **PAM 长度分离**：不同长度的分数不可直接比较；
- **无 GPU 验证**：本轮全部为 **CPU** 推理，CUDA 路径未测；
- **无湿实验**：本项目未做任何湿实验验证；
- **上游许可为非商业**：CC-BY-NC-4.0，商业使用需另行取得授权。

## 6. 复现校验值

```text
model  : Profluent-Bio/protein2pam-cas9_full
rev    : 407f7fc32146a4c0db13c05284f9f3a7cf0ff612
weights: aa1dc5c017ddd885bf26c8126675075af0007d1af3ab1a2893e62aaffe0341c2
config : 7fe91bd35342572b6dcd23f937013ec2603b2e136970f4e3a21a9af5b5c41cb5
```
