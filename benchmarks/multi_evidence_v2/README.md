# Multi-evidence benchmark v2

v2 separates biological evidence tier from available product input track.

Tracks:

- `protein_only`: a Cas protein and experimental PAM spectrum;
- `paired_evidence`: the same biological system also has oriented spacers and a
  frozen target collection;
- `published_flanks`: the same system has published per-spacer aggregate flank
  evidence, without raw target identifiers;
- `end_to_end_genome`: a host genome is supplied and Cas/array discovery is
  part of the evaluated pipeline.

A system may support more than one track. Strict independence is audited per
track: a protein can be independent while its target collection is not
virus-clustered, but that system must not be counted as strict paired evidence.

Validate without loading a model:

```bash
python scripts/validate_benchmark_manifest.py
```

需要同时确认本地输入没有被替换时：

```bash
python scripts/validate_benchmark_manifest.py --check-files
```

它会分别校验选中的蛋白 FASTA record、spacer、target 和 published-flank 文件，
不会加载 Protein2PAM 或使用 GPU。

To require the planned minimum of five strict paired systems:

```bash
python scripts/validate_benchmark_manifest.py --require-strict-paired
```

The current manifest contains 20 systems: one SpCas9 engineering control,
18 retrospective biological published-flank systems, and one strict-independent
protein-only system (Cj4Cas9, EFC33367.1). It validates all 58 registered
artifacts. The strict protein-only count is one; the strict paired-evidence count
remains zero because no independent virus-clustered target collection is yet
available. Copy `system_template.json` during curation; never insert placeholder
systems into the frozen manifest.

Expansion audit artifacts:

- `expansion_candidates.tsv`: all 644 qualified accession/subtype cohorts;
- `expansion_priority_A.tsv`: the 106 unregistered A-tier curation targets;
- `gasiunas_experimental_catalog.tsv`: 79 experimental rows with protein sequence;
- `expansion_all_experimental_crossmatch.tsv`: all-cohort conservative species/strain matching;
- `exact_strain_protein_verification_all.tsv`: exact-strain GenBank translation checks;
- `species_only_exact_protein_screen.tsv`: 170 species-only pairs screened by exact
  experimental-protein SHA256, yielding six deduplicated systems now in the manifest.
- `recent_experimental_cas9_audit.tsv`: four Cas9s with empirical randomized-7N
  PAM profiles from Becker et al. (2025), kept outside the manifest because two
  are exact Protein2PAM training exposures and two have >=90% training neighbours.
- `recent_strict_candidate_audit.tsv`: CoCas9, AalCas9, and Cj4Cas9 checked in
  the latest literature pass. Cj4Cas9 is registered as the first strict
  protein-only system; CoCas9 has a >=90% training neighbour, while AalCas9 is
  held out pending an exact paper-to-accession link. Cj4Cas9's nearest official
  neighbor is 0.89441624 identical, so it is explicitly treated as a cutoff-
  boundary case rather than a distant-OOD protein.

Reproduce the recent-paper audit with:

```bash
PYTHONPATH=src python scripts/audit_recent_experimental_cas9.py
PYTHONPATH=src python scripts/audit_recent_strict_candidates.py
```
