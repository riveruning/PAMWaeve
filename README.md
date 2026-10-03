# PAMWeave

**Evidence-aware candidate PAM ranking for Cas9.**

PAMWeave combines Protein2PAM-based candidate compatibility with spacer/protospacer flank evidence. It keeps the evidence channels visible, checks disagreements, and can abstain from joint ranking. It supports both raw spacer–target alignments and published aggregate flanks, with distinct provenance labels.

The release contains the analysis code and reproducible examples from the competition submission. The original internal package name `pamdict` and configuration prefix `PAMPRIDICT_` are retained for compatibility. The optional web interface and LLM assistant are not included.

## What the scores mean

Scores express model compatibility or sequence-evidence support. They are **not cleavage activity, editing efficiency, or experimental success probabilities**. The current evaluation does not establish independent generalization or an improvement from fusion. No wet-lab validation has been completed.

## Quick start

Python 3.10 is the tested interpreter. CPU inference is supported; GPU inference has not been validated for this release.

```bash
git clone https://github.com/riveruning/PAMWaeve.git
cd PAMWaeve
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
export PYTHONPATH="$PWD/src:${PYTHONPATH:-}"

# Check the 21 systems distributed with this repository; no weights needed.
python scripts/validate_benchmark_manifest.py \
  --manifest benchmarks/multi_evidence_v2/manifest_redistributable21.json \
  --check-files

# Run the package's offline regression suite.
python scripts/run_tests.py
```

For real predictions, download the pinned model weights and upstream adapter following [the reproduction guide](docs/REPRODUCTION.md#1-环境配置), then run:

```bash
bash run_all.sh
```

To rerun inference in a separate output directory:

```bash
bash run_all.sh --fresh --outdir data/parsed/fresh_repro21
```

Every run saves actual stdout/stderr, exit codes, timing, environment versions,
commands, source/output hashes and the bootstrap seed in `logs/runs/<UTC time>/`.
See [data instructions](data/README.md) and [execution records](logs/README.md).
Checked-in local verification logs are in `logs/verification_20261004/`.

The default run evaluates the 21 redistributable systems and generates `results/results_repro21.csv` (20 candidate rows from four showcase systems). The supplied `results/results.csv` contains 25 rows from the full 22-system reference evaluation. Reproducing its additional SpCas9 case requires separately obtaining the omitted inputs; the guide explains the distinction. Existing artifacts are reused by `run_all.sh`; use the explicit evaluation command in the guide for a fresh run in a new output directory.

## Documentation

- [Detailed setup and reproduction](docs/REPRODUCTION.md)
- [Model versions, downloads and limitations](models/MODEL_CARD.md)
- [Frozen evaluation protocol](docs/EVALUATION_PROTOCOL.md)
- [Results and failure analysis](docs/RESULTS_ANALYSIS.md)
- [Case selection](docs/CASE_SELECTION.md)
- [Data and software attribution](docs/THIRD_PARTY.md)
- [Original submission verification record](docs/SUBMISSION.md)
- [Delivery materials update](docs/MATERIALS_UPDATE_20261004.md)
- [Executable reproduction Notebook](notebooks/01_reproduce_candidates.ipynb)

The directory layout follows the competition's functional requirements. Data is
distributed under `benchmarks/`, with `data/README.md` explaining its location.
`run_all.sh` is the equivalent main prediction/screening entry; filenames such as
`predict.py` are not mandatory. The final results use pretrained weights, with
historical, unused fine-tuning experiments disclosed in the Model Card. No new
training is needed to reproduce them. An executable Notebook is also provided
to walk through the same entry and check the candidate list.

The original evaluation used 22 systems; three qualify for strict protein-only evaluation, while none qualify for strict paired-evidence evaluation. Training-exposed and retrospective examples are reported separately from independent samples. See the analysis for cohort sizes, denominators and abstention handling.

## License

Original project code and original software documentation are released under the [MIT License](LICENSE). Redistributed third-party datasets keep their original licenses and attribution; see [NOTICE](NOTICE.md). Model weights and the upstream Protein2PAM implementation are downloaded separately and have noncommercial terms. The MIT license does not change those terms.

## Contributing

Please include a minimal reproducer, input provenance where shareable, Python/dependency versions, and expected versus observed behavior when reporting an issue. See [CONTRIBUTING.md](CONTRIBUTING.md).

## 中文说明

PAMWeave 面向 Cas9 候选 PAM 筛选，分别保留蛋白预测和 spacer 侧翼证据，并展示一致、冲突及弃权。它用于安排和核查候选，不输出切割活性结论。安装、真实推理和结果字段的中文说明见[详细复现指南](docs/REPRODUCTION.md)。
