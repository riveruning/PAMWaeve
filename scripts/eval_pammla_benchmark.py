"""Score external 64-PAM predictions without training or inference."""
import argparse
import csv
import json
import math
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from build_pammla_benchmark import DEST
from prepare_pammla_pilot import sha,save_json
from score_pammla_pilot import spearman,overlap


def score_vectors(predicted, target, wt_prediction):
    rho=spearman(predicted,target);control=spearman(wt_prediction,target)
    top=overlap(predicted,target);wt_top=overlap(wt_prediction,target)
    return {'spearman':rho,'wt_transfer_spearman':control,
            'delta_spearman':rho-control if rho is not None and control is not None else None,
            'top5_overlap':top,'wt_transfer_top5_overlap':wt_top,'delta_top5_overlap':top-wt_top}


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def select_scope(variants, scope):
    if scope not in {'development','all'}:raise ValueError('Unknown scope')
    selected=[r for r in variants if scope=='all' or r['future_task_split']=='development']
    if not any(r['variant_id']=='DSGERT' for r in selected):raise ValueError('WT control missing')
    return selected


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pred',type=Path,required=True)
    ap.add_argument('--model-name',required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--scope',choices=['development','all'],default='development',
                    help='Default protects reserved evaluation; all must be explicitly requested')
    args=ap.parse_args()
    manifest=json.loads((DEST/'manifest.json').read_text())
    for name,digest in manifest['files'].items():assert sha((DEST/name).read_bytes())==digest
    variants=select_scope(read(DEST/'variants.tsv'),args.scope)
    names={r['variant_id'] for r in variants}
    measurements=[r for r in read(DEST/'measurements.tsv') if r['variant_id'] in names]
    targets={(r['variant_id'],r['pam_group']):float(r['ranking_target_log10_k']) for r in measurements}
    predictions={}
    for r in read(args.pred):
        key=r['variant_id'],r['pam_group']
        if key in predictions:raise ValueError(f'Duplicate prediction: {key}')
        value=float(r['score'])
        if not math.isfinite(value):raise ValueError(f'Non-finite prediction: {key}')
        predictions[key]=value
    if predictions.keys()!=targets.keys():
        raise ValueError(f'Prediction coverage mismatch: missing={len(targets.keys()-predictions.keys())}, extra={len(predictions.keys()-targets.keys())}')
    pams=sorted({p for _,p in targets});assert len(pams)==64
    wt=[predictions['DSGERT',p] for p in pams]
    per_variant=[]
    for row in variants:
        name=row['variant_id'];target=[targets[name,p] for p in pams]
        per_variant.append({'variant_id':name,'is_wt':row['is_wt']=='True',
            'selection_group':row['selection_group'],'future_task_split':row['future_task_split'],
            'hamming2_component':row['hamming2_component'],
            'official_best_identity_ge85':float(row['official_best_identity_ge85']),
            'F0_best_identity_ge85':float(row['F0_best_identity_ge85']),
            'constant_target':len(set(target))==1,
            **score_vectors([predictions[name,p] for p in pams],target,wt)})
    groups={'all_non_wt':[r for r in per_variant if not r['is_wt']]}
    for column in ['selection_group','future_task_split']:
        for value in sorted({r[column] for r in per_variant}):
            groups[f'{column}={value}']=[r for r in per_variant if not r['is_wt'] and r[column]==value]
    summaries={}
    for name,items in groups.items():
        paired=[r for r in items if r['delta_spearman'] is not None]
        summary={'n':len(items),'constant_targets':sum(r['constant_target'] for r in items),
                 'defined_spearman_pairs':len(paired),'undefined_spearman_pairs':len(items)-len(paired)}
        summary['positive_spearman_gains']=sum(r['delta_spearman']>0 for r in paired)
        summary['negative_spearman_gains']=sum(r['delta_spearman']<0 for r in paired)
        for key in ['spearman','wt_transfer_spearman','delta_spearman','top5_overlap','wt_transfer_top5_overlap','delta_top5_overlap']:
            values=[r[key] for r in (paired if 'spearman' in key else items)]
            summary['mean_'+key]=float(np.mean(values)) if values else None
        if paired:
            clusters={r['hamming2_component'] for r in paired}
            sums=np.array([sum(r['delta_spearman'] for r in paired if r['hamming2_component']==c) for c in sorted(clusters)])
            counts=np.array([sum(r['hamming2_component']==c for r in paired) for c in sorted(clusters)])
            rng=np.random.default_rng(20260927)
            indices=rng.integers(len(clusters),size=(2000,len(clusters)))
            boot=sums[indices].sum(axis=1)/counts[indices].sum(axis=1)
            summary['exploratory_component_bootstrap_delta_spearman_95pct']=np.quantile(boot,[.025,.975]).tolist()
            summary['components']=len(clusters)
        summaries[name]=summary
    save_json(args.out,{'model_name':args.model_name,'scope':args.scope,
        'uncertainty':'Exploratory 2000 component-bootstrap resamples, seed20260927; development data, not confirmatory significance',
        'benchmark_manifest_sha256':sha((DEST/'manifest.json').read_bytes()),
        'predictions_sha256':sha(args.pred.read_bytes()),'prediction_rows':len(predictions),
        'limitations':'One engineered SpCas9 family, official-catalog near-neighbor exposure for all; published PAMmla training data; not cell activity calibration or natural-family independence',
        'aggregate':summaries,'per_variant':per_variant})
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
