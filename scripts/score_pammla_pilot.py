"""Evaluate frozen 64-group quantitative PAM ranking, never fit on labels."""
import csv
import hashlib
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
from pamdict.score.candidate import score_candidate_pam
from pamdict.score.spectrum import stable_benchmark_row_ids
from prepare_pammla_pilot import OUT, save_json, sha


def ranks(values):
    a = np.asarray(values, dtype=float)
    if a.ndim != 1 or not len(a) or not np.isfinite(a).all():
        raise ValueError('Expected nonempty finite vector')
    order = np.argsort(a, kind='stable')
    result = np.empty(len(a), dtype=float)
    i = 0
    while i < len(a):
        j = i + 1
        while j < len(a) and a[order[j]] == a[order[i]]:
            j += 1
        result[order[i:j]] = (i+j-1)/2 + 1
        i = j
    return result


def spearman(left, right):
    if len(left) != len(right):
        raise ValueError('Length mismatch')
    a, b = ranks(left), ranks(right)
    a -= a.mean(); b -= b.mean()
    denominator = np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.dot(a,b)/denominator) if denominator else None


def top_membership(values, k=5):
    a = np.asarray(values, dtype=float)
    if a.ndim != 1 or not np.isfinite(a).all() or not 1 <= k <= len(a):
        raise ValueError('Invalid top-k input')
    cutoff = np.sort(a)[-k]
    result = (a > cutoff).astype(float)
    tied = a == cutoff
    result[tied] = (k-result.sum())/tied.sum()
    assert abs(result.sum()-k) < 1e-10
    return result


def overlap(left, right, k=5):
    if len(left) != len(right):
        raise ValueError('Length mismatch')
    return float(np.dot(top_membership(left,k), top_membership(right,k))/k)


def main():
    inputs = json.loads((OUT/'inputs.json').read_text())
    checks = {'bench.tsv':'bench_sha256','backbone_evidence.json':'backbone_evidence_sha256',
              'exposure.json':'exposure_sha256'}
    for file, key in checks.items():
        assert sha((OUT/file).read_bytes()) == inputs[key]
    frozen_path = ROOT/'data/parsed/pammla_pilot_v1.json'
    assert sha(frozen_path.read_bytes()) == inputs['frozen_selection_sha256']
    frozen = json.loads(frozen_path.read_text())
    source = ROOT/'data/raw/pam_prediction_papers/PAMmla/220924_6pos_selected_random_only.csv'
    assert sha(source.read_bytes()) == inputs['csv_sha256']
    checkpoint = ROOT/'data/checkpoints/finetune_v3_cas9_F0_equal/checkpoint_step000020.pt'
    assert sha(checkpoint.read_bytes()) == inputs['checkpoint_sha256']
    with source.open(encoding='utf-8-sig') as f:
        data = {r['name']:r for r in csv.DictReader(f)}
    with (OUT/'bench.tsv').open() as f:
        bench = list(csv.DictReader(f, delimiter='\t'))
    ids = stable_benchmark_row_ids(bench)
    pams = ['N'+''.join(x) for x in itertools.product('ACGT',repeat=3)]
    score_vectors = {}; pred_hashes = {}
    for model in ['base','F0_step20']:
        file = OUT/(model+'.tsv'); pred_hashes[model] = sha(file.read_bytes())
        with file.open() as f:
            predictions = list(csv.DictReader(f,delimiter='\t'))
        by_id = {p['row_id']:p for p in predictions}
        assert len(predictions)==len(by_id)==len(bench) and set(by_id)==set(ids)
        score_vectors[model] = {}
        for row_id, row in zip(ids,bench):
            p = by_id[row_id]
            assert p['sequence_sha256'] == sha(row['protein_sequence'].encode())
            assert p['protein_id'] == row['protein_id']
            matrix = json.loads(p['pred_probability_json'])
            # Unrounded internal scores: display rounding must not introduce ties.
            vector = [score_candidate_pam(matrix,pam,side='downstream').specificity_adjusted_score for pam in pams]
            score_vectors[model][row['protein_id']] = vector
    result = {'inputs':inputs, 'prediction_sha256':pred_hashes,'pams':pams,
              'score_vectors':score_vectors, 'variants':[], 'aggregate':{}}
    snapshot = ROOT/'data/checkpoints/hf/hub/models--Profluent-Bio--protein2pam-cas9_full/snapshots/407f7fc32146a4c0db13c05284f9f3a7cf0ff612'
    model_hash = hashlib.sha256()
    with (snapshot/'model.safetensors').open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):
            model_hash.update(chunk)
    result['implementation'] = {
        'base_snapshot':snapshot.name, 'base_weights_sha256':model_hash.hexdigest(),
        'base_config_sha256':sha((snapshot/'config.json').read_bytes()),
        'inference':'CUDA BF16, batch_size=1, eval mode, local weights only',
        'code_sha256':{str(p):sha((ROOT/p).read_bytes()) for p in [
            'scripts/prepare_pammla_pilot.py','scripts/score_pammla_pilot.py',
            'scripts/predict_benchmark_base.py','scripts/eval_finetuned.py',
            'src/pamdict/infer/p2pam.py','src/pamdict/score/candidate.py']}}
    for sample in frozen['samples']:
        name = sample['variant']
        raw = np.array([float(data[name][pam]) for pam in pams])
        assert np.isfinite(raw).all()
        target = np.maximum(raw,-5)
        item = {**sample,'clipped_target':target.tolist(),
                'censored_groups':int((raw <= -5).sum()),'unique_target_values':len(set(target)),
                'constant_target':len(set(target))==1,'models':{}}
        for model, vectors in score_vectors.items():
            predicted = vectors[name]; control = vectors['DSGERT']
            rho = spearman(predicted,target); wt_rho = spearman(control,target)
            top = overlap(predicted,target); wt_top = overlap(control,target)
            item['models'][model] = {'spearman':rho,'wt_transfer_spearman':wt_rho,
                'delta_spearman':rho-wt_rho if rho is not None and wt_rho is not None else None,
                'top5_overlap':top,'wt_transfer_top5_overlap':wt_top,'delta_top5_overlap':top-wt_top,
                'score_vector_correlation_to_wt':spearman(predicted,control),
                'max_absolute_score_change_from_wt':float(np.max(np.abs(np.asarray(predicted)-control)))}
        result['variants'].append(item)
    for group in ['selected','random','all_non_wt']:
        cohort = [r for r in result['variants'] if not r['is_wt'] and (group=='all_non_wt' or r['group']==group)]
        result['aggregate'][group] = {}
        for model in score_vectors:
            ms = [r['models'][model] for r in cohort]
            values = [r['delta_spearman'] for r in ms if r['delta_spearman'] is not None]
            summary = {'n':len(ms),'defined_spearman_pairs':len(values),
                       'undefined_spearman_pairs':len(ms)-len(values),
                       'positive_spearman_deltas':sum(v>0 for v in values)}
            for key in ['spearman','wt_transfer_spearman','delta_spearman','top5_overlap','wt_transfer_top5_overlap','delta_top5_overlap']:
                defined = [r[key] for r in ms if r[key] is not None]
                summary['mean_'+key] = float(np.mean(defined)) if defined else None
                summary['median_'+key] = float(np.median(defined)) if defined else None
            result['aggregate'][group][model] = summary
    save_json(OUT/'results.json',result)
    print(json.dumps(result['aggregate'],indent=2))


if __name__ == '__main__': main()
