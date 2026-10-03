"""Validate packaged source linkage, measurements, splits and exposure flags."""
import csv
import json
import math
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from build_pammla_benchmark import DEST, SOURCE, splits
from prepare_pammla_pilot import sha


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    manifest=json.loads((DEST/'manifest.json').read_text())
    for name,digest in manifest['files'].items():assert sha((DEST/name).read_bytes())==digest,name
    assert sha(SOURCE.read_bytes())==manifest['source_sha256']
    with SOURCE.open(encoding='utf-8-sig') as f:source={r['name']:r for r in csv.DictReader(f)}
    variants=read(DEST/'variants.tsv');measurements=read(DEST/'measurements.tsv')
    quarantined=read(DEST/'quarantine.tsv');names={r['variant_id'] for r in variants}
    assert len(variants)==len(names)==739 and len(measurements)==739*64
    assert names.isdisjoint({r['variant_id'] for r in quarantined})
    assert names|{r['variant_id'] for r in quarantined}==set(source)
    assert len(quarantined)==31
    assert len({r['sequence_sha256'] for r in variants})==739
    exposure=json.loads((DEST/'exposure.json').read_text())
    for ref in exposure['references'].values():assert sha((ROOT/ref['path']).read_bytes())==ref['sha256']
    lookup={r['variant_id']:r for r in variants};groups={n:set() for n in names}
    for r in variants:
        assert len(r['protein_sequence'])==1368
        assert sha(r['protein_sequence'].encode())==r['sequence_sha256']
        assert r['source_library']==source[r['variant_id']]['source library']
        for label in ['official','F0']:
            found=exposure['variants'][r['variant_id']][label]
            assert r[label+'_exact_match']==str(found['exact'])
            assert int(r[label+'_near85_count'])==found['near85_count']
            assert int(r[label+'_near90_count'])==found['near90_count']
    for r in measurements:
        name,pam=r['variant_id'],r['pam_group'];assert name in names and pam not in groups[name]
        assert len(pam)==4 and pam[0]=='N' and set(pam[1:])<=set('ACGT')
        groups[name].add(pam);raw=float(r['raw_log10_k']);target=float(r['ranking_target_log10_k'])
        assert math.isfinite(raw) and raw==float(source[name][pam]) and target==max(raw,-5)
        assert r['at_or_below_floor']==str(raw<=-5)
    assert all(len(v)==64 for v in groups.values())
    pilot=json.loads((ROOT/'data/parsed/pammla_pilot_v1.json').read_text())
    seen={s['variant'] for s in pilot['samples']}
    assignments=splits(sorted(names),seen)
    for name,(group,role) in assignments.items():
        assert lookup[name]['hamming2_component']==group and lookup[name]['future_task_split']==role
        if name in seen:assert role=='development'
    print('PASS: 739 unique proteins, 47296 source-linked measurements, 31 quarantined; hashes, exposure and component splits verified')


if __name__=='__main__':main()
