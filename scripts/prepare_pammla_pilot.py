"""Verify source-linked Cas9 cores and prepare the frozen quantitative pilot.

The small GenBank parser is intentionally specific to the verified BPK848 file;
it fails closed on any other file. No arbitrary joined feature parser is used.
"""
import argparse
import csv
import hashlib
import itertools
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from audit_pp_df_training_neighbors import audit, rows, sha
from audit_recent_strict_candidates import is_cas9

OUT = ROOT / 'data/parsed/pammla_pilot_v1'
RAW = ROOT / 'data/raw/pammla_backbone'
POSITIONS = (1135, 1136, 1218, 1219, 1335, 1337)
GENBANK_SHA = '499cc828be9e4a0325ea3284b38a42ead62a64481f3a04fc36f1ad3a02ab5e17'
REF_URL = 'https://raw.githubusercontent.com/RachelSilverstein/multiplex_seq_analysis/main/reference_MSP2582_extraction.fasta'
TABLE = dict(zip(map(''.join, itertools.product('TCAG', repeat=3)),
                 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


def save(path, content):
    """Never replace different existing evidence."""
    if path.exists() and path.read_bytes() != content:
        raise ValueError(f'Existing artifact differs: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def save_json(path, value):
    save(path, (json.dumps(value, indent=2) + '\n').encode())


def translate(dna):
    if len(dna) % 3 or set(dna) - set('ACGT'):
        raise ValueError('Invalid coding DNA')
    return ''.join(TABLE[dna[i:i+3]] for i in range(0, len(dna), 3))


def verify_core(gbk, reference):
    if sha(gbk) != GENBANK_SHA:
        raise ValueError('Unexpected GenBank; requires a new explicit coordinate audit')
    text = gbk.decode()
    dna = ''.join(re.findall('[acgt]+', text.split('ORIGIN')[1].split('//')[0])).upper()
    assert len(dna) == 8021
    feature = re.search(r'     CDS +761\.\.4864\n(.*?)(?=\n     \S)', text, re.S)
    assert feature and '/label=Cas9\n' in feature.group(1)
    annotation = ''.join(re.search(r'/translation="([^"]+)"', feature.group(1)).group(1).split())
    core = translate(dna[760:4864])
    assert core == annotation and len(core) == 1368 and '*' not in core
    assert ''.join(core[p-1] for p in POSITIONS) == 'DSGERT'
    lines = reference.decode().splitlines()
    assert sum(line.startswith('>') for line in lines) == 1 and 'MSP2582' in lines[0]
    downstream = ''.join(lines[1:]).upper()
    hits = []
    for frame in range(3):
        coding = downstream[frame:frame + ((len(downstream)-frame)//3)*3]
        protein = translate(coding)
        for match in re.finditer(re.escape(core), protein):
            hits.append(frame + match.start()*3)
    assert len(hits) == 1, 'Downstream reference must contain one exact full Cas9 core'
    return core, {'genbank_sha256': sha(gbk), 'genbank_core_1based_inclusive': [761,4864],
                  'core_length': len(core), 'core_sha256': sha(core.encode()),
                  'six_residues': 'DSGERT', 'reference_sha256': sha(reference),
                  'reference_url': REF_URL, 'reference_header': lines[0],
                  'reference_core_1based_inclusive': [hits[0]+1, hits[0]+len(core)*3],
                  'downstream_core_exact_match': True,
                  'annotation_warning': 'Overlapping VQR/VRER labels do not describe the actual translated six residues; actual DNA and plain Cas9 feature agree.',
                  'model_input_policy': '1368-aa Cas9 core only; vector NLS/FLAG/P2A/EGFP excluded; no claim of cell-activity calibration',
                  'mutation_evidence': 'Paper methods: whole-ORF sequencing excluded indels, mixed colonies and unintended point mutations; six-site variants used in HT-PAMDA.',
                  'paper_url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12449813/',
                  'library_spelling': 'Methods use both MNW94 (library construction) and MMW94 (random sampling); CSV MMW94 retained literally.'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--genbank', required=True)
    ap.add_argument('--assay-reference', required=True)
    args = ap.parse_args()
    gbk = Path(args.genbank).read_bytes()
    reference = Path(args.assay_reference).read_bytes()
    core, provenance = verify_core(gbk, reference)
    save(RAW/'BPK848_addgene_181745.gbk', gbk)
    save(RAW/'reference_MSP2582_extraction.fasta', reference)
    frozen_path = ROOT/'data/parsed/pammla_pilot_v1.json'
    frozen = json.loads(frozen_path.read_text())
    source = ROOT/'data/raw/pam_prediction_papers/PAMmla/220924_6pos_selected_random_only.csv'
    assert sha(source.read_bytes()) == frozen['source_sha256']
    with source.open(encoding='utf-8-sig') as f:
        dataset = {r['name']: r for r in csv.DictReader(f)}
    records = []
    for sample in frozen['samples']:
        name = sample['variant']; row = dataset[name]
        assert row['source library'] == 'MMW94 - 6 position library'
        assert row['status'] == 'good'
        assert row['extra mutations (other than 14 positions)?'] == 'no'
        assert ''.join(row[f'{wt}{p}'] for wt,p in zip('DSGERT', POSITIONS)) == name
        sequence = list(core)
        for p, residue in zip(POSITIONS, name):
            assert residue in 'ACDEFGHIKLMNPQRSTVWY'
            sequence[p-1] = residue
        sequence = ''.join(sequence)
        assert all(a==b or i+1 in POSITIONS for i,(a,b) in enumerate(zip(core,sequence)))
        records.append({'protein_id':name, 'protein_sequence':sequence,
                        'pam_consensus':'NNNN', 'group':sample['group'],
                        'sequence_sha256':sha(sequence.encode())})
    assert len(records) == 13 and len({r['protein_sequence'] for r in records}) == 13
    save_json(OUT/'backbone_evidence.json', provenance)
    import io
    buf = io.StringIO(newline='')
    writer = csv.DictWriter(buf, fieldnames=list(records[0]), delimiter='\t')
    writer.writeheader(); writer.writerows(records)
    save(OUT/'bench.tsv', buf.getvalue().encode())
    # Triangle inequality: variants differ from WT by <=6 edits. With max
    # length >=1368, any variant match >=.85 must match WT at >=.84.
    # Exhaustively screen WT at .84, then verify every retained sequence
    # for each variant at .85; no heuristic k-mer cap.
    exposure = {'method':'Exhaustive WT global-edit identity >=.84 then per-variant >=.85; safe triangle bound for <=6 substitutions and length1368',
                'references':{}, 'variants':{r['protein_id']:{} for r in records}}
    for relative in ['data/raw/protein2pam_train_seqs.tsv',
                     'data/corpus/augmented_train_v2_cas9_train.tsv']:
        path = ROOT/relative; refs = rows(path)
        if relative.startswith('data/raw'):
            refs = [r for r in refs if is_cas9(r)]
        print('Exhaustive WT screening:',relative,len(refs),flush=True)
        screen = audit(core, refs, threshold=.84)
        hashes = {m['sequence_sha256'] for m in screen['matches_at_or_above_85']}
        nearby = [r for r in refs if sha(r['protein_sequence'].strip().upper().encode()) in hashes]
        exposure['references'][relative] = {'sha256':sha(path.read_bytes()),'rows':len(refs),
            'wt_screen_threshold':.84, 'screened_unique':screen['length_eligible_sequences_checked'],
            'retained_rows':len(nearby)}
        for row in records:
            result = audit(row['protein_sequence'], nearby)
            exposure['variants'][row['protein_id']][relative] = result
            print(row['protein_id'],'exact',result['exact'],'near85',result['near85_count'],flush=True)
    save_json(OUT/'exposure.json', exposure)
    checkpoint = ROOT/'data/checkpoints/finetune_v3_cas9_F0_equal/checkpoint_step000020.pt'
    save_json(OUT/'inputs.json', {'status':'READY_FOR_FIXED_MODEL_PILOT',
        'frozen_selection_sha256':sha(frozen_path.read_bytes()),'csv_sha256':frozen['source_sha256'],
        'backbone_evidence_sha256':sha((OUT/'backbone_evidence.json').read_bytes()),
        'exposure_sha256':sha((OUT/'exposure.json').read_bytes()),
        'bench_sha256':sha((OUT/'bench.tsv').read_bytes()),'checkpoint_sha256':sha(checkpoint.read_bytes()),
        'scope':'13 closely related engineered SpCas9 cores; not independent natural systems',
        'score':'Existing specificity_adjusted_score, full precision, downstream NXYZ at positions1..4; N ignored, no expansion to 256 independent observations',
        'metrics':frozen['metrics'],
        'top5_ties':'Expected intersection under independent uniform tie-breaking: dot(fractional_memberships)/5',
        'constant_predictions_or_targets':'Undefined Spearman (JSON null), excluded only from defined-pair aggregate; report counts',
        'placeholder_label':'bench pam_consensus NNNN is a schema placeholder, not a gold label; do not run motif-accuracy scorer'} )
    print('READY:',OUT,flush=True)


if __name__ == '__main__':
    main()
