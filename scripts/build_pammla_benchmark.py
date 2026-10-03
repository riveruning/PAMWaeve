"""Package source-verified six-site variants; no training or model inference."""
import collections
import csv
import io
import itertools
import json
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from prepare_pammla_pilot import RAW, POSITIONS, verify_core, save, save_json, sha
from audit_pp_df_training_neighbors import audit, rows
from audit_recent_strict_candidates import is_cas9

DEST=ROOT/'benchmarks/pammla_quantitative_v1'
SOURCE=ROOT/'data/raw/pam_prediction_papers/PAMmla/220924_6pos_selected_random_only.csv'
SOURCE_SHA='701aa1f34a75c107e23722ca2da91b97e501f34ba4041a93a41adeb8ce62f92c'


def components(names, distance=2):
    if len(set(names))!=len(names) or any(len(n)!=6 for n in names):
        raise ValueError('Unique six-residue codes required')
    parent=list(range(len(names)))
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,name in enumerate(names):
        for j in range(i):
            if sum(a!=b for a,b in zip(name,names[j]))<=distance:
                parent[find(i)]=find(j)
    grouped=collections.defaultdict(list)
    for i,name in enumerate(names):grouped[find(i)].append(name)
    return sorted((sorted(g) for g in grouped.values()),key=lambda g:g[0])


def splits(names, previously_used):
    assignments={}
    for group in components(names):
        fingerprint=sha('|'.join(group).encode())
        # Never promote an already evaluated pilot or its close-neighbor
        # component to an untouched reserve. Select only by identifiers.
        role='development'
        if not set(group)&set(previously_used) and int(fingerprint[:8],16)%5==0:
            role='reserved_evaluation'
        for name in group:assignments[name]=(fingerprint,role)
    return assignments


def table(path, data, fields=None):
    fields=fields or list(data[0])
    buf=io.StringIO(newline='');w=csv.DictWriter(buf,fieldnames=fields,delimiter='\t')
    w.writeheader();w.writerows(data);save(path,buf.getvalue().encode())


def main():
    assert sha(SOURCE.read_bytes())==SOURCE_SHA
    with SOURCE.open(encoding='utf-8-sig') as f:all_rows=list(csv.DictReader(f))
    assert len(all_rows)==len({r['name'] for r in all_rows})==770
    eligible=[r for r in all_rows if r['source library']=='MMW94 - 6 position library']
    assert len(eligible)==739
    core,evidence=verify_core((RAW/'BPK848_addgene_181745.gbk').read_bytes(),
                             (RAW/'reference_MSP2582_extraction.fasta').read_bytes())
    pilot=json.loads((ROOT/'data/parsed/pammla_pilot_v1.json').read_text())
    seen={s['variant'] for s in pilot['samples']}
    assignment=splits([r['name'] for r in eligible],seen)
    pams=['N'+''.join(p) for p in itertools.product('ACGT',repeat=3)]
    variants=[];measurements=[];seqs={};censored=0;constant=0
    for r in eligible:
        name=r['name'];assert r['status']=='good' and r['extra mutations (other than 14 positions)?']=='no'
        assert ''.join(r[f'{aa}{p}'] for aa,p in zip('DSGERT',POSITIONS))==name
        assert set(name)<=set('ACDEFGHIKLMNPQRSTVWY')
        seq=list(core)
        for p,aa in zip(POSITIONS,name):seq[p-1]=aa
        seq=''.join(seq);seqs[name]=seq
        assert len(seq)==1368 and all(a==b or i+1 in POSITIONS for i,(a,b) in enumerate(zip(core,seq)))
        values=[float(r[p]) for p in pams];assert all(math.isfinite(v) for v in values)
        clipped=[max(v,-5) for v in values];is_constant=len(set(clipped))==1
        constant+=is_constant
        component,role=assignment[name]
        variants.append({'variant_id':name,'protein_sequence':seq,'sequence_sha256':sha(seq.encode()),
            'protein_length':len(seq),'sample':r['sample'],'prep':r['prep'],
            'selection_group':r['selected/rational'],'source_library':r['source library'],
            'is_wt':name=='DSGERT','previously_evaluated_pilot':name in seen,
            'distance_from_wt_6sites':sum(a!=b for a,b in zip(name,'DSGERT')),
            'constant_target_after_floor':is_constant,'family_group':'SpCas9_BPK848_MSP2582',
            'hamming2_component':component,'future_task_split':role,
            'strict_natural_generalization':False,'published_pammla_training_exposure':True})
        for pam,value,target in zip(pams,values,clipped):
            censored+=value<=-5
            measurements.append({'variant_id':name,'pam_group':pam,'raw_log10_k':value,
                                 'ranking_target_log10_k':target,'at_or_below_floor':value<=-5})
    assert len(set(seqs.values()))==739
    exposure={'method':'Exhaustive WT global-edit >=.84 then each variant >=.85; safe triangle bound <=6 substitutions /1368 aa',
              'identity_definition':'1-Levenshtein_distance/max_length; near counts are unique sequences',
              'references':{},'variants':{name:{} for name in seqs}}
    gap={'benchmark':{'rows':739,'independent_backbone_families':1,'assay':'HT-PAMDA',
        'target':'64 PAM-group log10 rate constants; floor -5 for rank evaluation',
        'condition':'Published processing averages two distinct spacer target assays; not guide-free intrinsic activity or in-cell editing efficiency',
        'input':'1368-aa engineered SpCas9 core only; at most six substitutions',
        'selected_random_counts':dict(collections.Counter(r['selected/rational'] for r in eligible)),
        'pammla_exposure':'Source CSV is published PAMmla training data; released PAMmla is not an independent baseline here'},
        'training_references':{}}
    for label,relative in [('official','data/raw/protein2pam_train_seqs.tsv'),
                           ('F0','data/corpus/augmented_train_v2_cas9_train.tsv')]:
        path=ROOT/relative;refs=rows(path)
        if label=='official':refs=[r for r in refs if is_cas9(r)]
        print('screen',label,len(refs),flush=True)
        screen=audit(core,refs,.84)
        hashes={m['sequence_sha256'] for m in screen['matches_at_or_above_85']}
        nearby=[r for r in refs if sha(r['protein_sequence'].strip().upper().encode()) in hashes]
        exposure['references'][label]={'path':relative,'sha256':sha(path.read_bytes()),'rows':len(refs),
            'wt_screen_threshold':.84,'wt_screen_unique_checked':screen['length_eligible_sequences_checked'],
            'retained_rows':len(nearby)}
        gap['training_references'][label]={**exposure['references'][label],
            'source_counts':dict(collections.Counter(r.get('source','') for r in refs)),
            'label_kind_counts':dict(collections.Counter(r.get('label_kind','not_recorded') for r in refs)),
            'label_semantics':'PAM consensus/per-position nucleotide logo, not PAM-specific kinetic rate constants',
            'scope_caveat':'Public Cas9 training catalog; checkpoint-specific original row consumption not independently proven' if label=='official' else 'Actual F0 input file, fingerprint checked against frozen pilot; F0 inherits base pretraining exposure'}
        if label=='F0':
            expected=json.loads((ROOT/'data/parsed/pammla_pilot_v1/exposure.json').read_text())['references'][relative]['sha256']
            assert expected==sha(path.read_bytes())
        for i,v in enumerate(variants):
            result=audit(seqs[v['variant_id']],nearby)
            exposure['variants'][v['variant_id']][label]=result
            v[label+'_exact_match']=result['exact'];v[label+'_near90_count']=result['near90_count']
            v[label+'_near85_count']=result['near85_count']
            v[label+'_best_identity_ge85']=result['matches_at_or_above_85'][0]['identity'] if result['near85_count'] else ''
            if (i+1)%100==0:print(label,i+1,'/739',flush=True)
    gap['sequence_exposure_summary']={label:{'exact_variants':sum(v[label+'_exact_match'] for v in variants),
        'variants_with_near90':sum(v[label+'_near90_count']>0 for v in variants),
        'variants_with_near85':sum(v[label+'_near85_count']>0 for v in variants)} for label in ['official','F0']}
    dev=[n for n,(_,role) in assignment.items() if role=='development']
    reserve=[n for n,(_,role) in assignment.items() if role=='reserved_evaluation']
    minimum=min(sum(a!=b for a,b in zip(x,y)) for x in dev for y in reserve)
    assert minimum>=3 and seen<=set(dev)
    manifest={'format_version':1,'status':'DATASET_READY_NOT_MODEL_VALIDATED','source_sha256':SOURCE_SHA,
        'source_url':'https://github.com/RachelSilverstein/PAMmla','paper_doi':'10.1038/s41586-025-09021-y',
        'variants':739,'pam_groups_per_variant':64,'measurements':len(measurements),
        'quarantined_rows':31,'censored_measurements':censored,'constant_target_variants':constant,
        'natural_independent_systems_added':0,'backbone_families':1,
        'model_evaluation_policy':'Score frozen Protein2PAM on all eligible non-WT rows by group; development vs reserve is for future task-model work, not an independence claim',
        'split_policy':'Hamming distance <=2 connected components over six sites; identifier hash modulo5 reserve; all previously evaluated pilot components forced development',
        'split_counts':dict(collections.Counter(role for _,role in assignment.values())),
        'split_component_count':len(set(group for group,_ in assignment.values())),
        'minimum_cross_split_6site_hamming':minimum,
        'reserve_caveat':'Not evaluated in the local 13-sample pilot, but public labels and PAMmla training exposure mean this is not a secret or universally unseen test set',
        'scoring_contract':{'unit':'one variant; 64 groups not independent proteins','primary':'per-variant Spearman with average ties and -5 floor',
            'secondary':'fractional Top5 overlap; separate selected/random','constant_target':'Spearman null, report counts; never silently discard',
            'control':'same-model WT prediction transferred to variants; WT excluded from gain average',
            'interpretation':'model compatibility ranking, not calibrated cleavage or cell-editing efficiency'},
        'exposure_summary':gap['sequence_exposure_summary'],'files':{}}
    table(DEST/'variants.tsv',variants);table(DEST/'measurements.tsv',measurements)
    table(DEST/'quarantine.tsv',[{'variant_id':r['name'],'source_library':r['source library'],
          'reason':'non-six-site-library; full mutation provenance not verified'} for r in all_rows if r not in eligible])
    save(DEST/'proteins.faa',''.join(f'>{n}\n{s}\n' for n,s in seqs.items()).encode())
    save_json(DEST/'backbone_evidence.json',evidence);save_json(DEST/'exposure.json',exposure)
    save_json(DEST/'training_gap.json',gap)
    for p in sorted(DEST.iterdir()):
        if p.name not in {'manifest.json','README.md'}:manifest['files'][p.name]=sha(p.read_bytes())
    save_json(DEST/'manifest.json',manifest)
    print(json.dumps(manifest,indent=2),flush=True)


if __name__=='__main__':main()
