"""Fixed-weight precision and representation diagnostics for the frozen pilot."""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'src'))
from scripts.prepare_pammla_pilot import OUT, POSITIONS, save_json, sha
from scripts.score_pammla_pilot import spearman, overlap
from pamdict.infer.p2pam import P2PAMPredictor
from pamdict.score.candidate import score_candidate_pam
from scripts.finetune_p2pam import apply_lora


def vector_change(value, reference):
    a=np.asarray(value,dtype=float).ravel(); b=np.asarray(reference,dtype=float).ravel()
    delta=np.linalg.norm(a-b); denom=np.linalg.norm(b)
    cosine_denom=np.linalg.norm(a)*denom
    return {'l2':float(delta),'relative_l2':float(delta/denom) if denom else None,
            'cosine':float(np.dot(a,b)/cosine_denom) if cosine_denom else None,
            'max_abs':float(np.max(np.abs(a-b)))}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--model',choices=['base','F0_step20'],required=True)
    ap.add_argument('--precision',choices=['fp32','bf16'],required=True)
    args=ap.parse_args()
    inputs=json.loads((OUT/'inputs.json').read_text())
    original=json.loads((OUT/'results.json').read_text())
    assert sha((OUT/'bench.tsv').read_bytes())==inputs['bench_sha256']
    with (OUT/'bench.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    destination=OUT/'diagnostics'/f'{args.model}_{args.precision}.json'
    if destination.exists():
        raise ValueError(f'Diagnostic already exists: {destination}; do not overwrite')
    assert torch.cuda.is_available()
    torch.backends.cuda.matmul.allow_tf32=False
    torch.backends.cudnn.allow_tf32=False
    torch.manual_seed(0)
    predictor=P2PAMPredictor('cas9_full',device='cpu',local_files_only=True)
    model=predictor.model
    checkpoint=ROOT/'data/checkpoints/finetune_v3_cas9_F0_equal/checkpoint_step000020.pt'
    if args.model=='F0_step20':
        assert sha(checkpoint.read_bytes())==inputs['checkpoint_sha256']
        ck=torch.load(checkpoint,map_location='cpu')
        assert ck['format_version']==3 and ck['classifier_mode']=='freeze' and ck['model_name']=='cas9_full'
        apply_lora(model,rank=int(ck['lora_rank']),alpha=int(ck['lora_alpha']))
        state=ck['trainable_state']; current=model.state_dict()
        expected={k for k in current if 'lora_' in k}
        assert set(state)==expected
        assert all(state[k].shape==current[k].shape for k in state)
        incompatible=model.load_state_dict(state,strict=False)
        assert not incompatible.unexpected_keys
        del ck,state,current
    dtype=torch.float32 if args.precision=='fp32' else torch.bfloat16
    model.to(dtype=dtype).to('cuda').eval()
    torch.cuda.reset_peak_memory_stats()
    capture={}; token_positions=[]
    def capture_input(module,args_):
        features=args_[0].detach()
        capture['cls']=features[0,0].float().cpu().tolist()
        capture['six_site_hidden']=features[0,token_positions].float().cpu().tolist()
    def capture_head(module,args_,output):
        capture['head_hidden']=output[0].detach().float().cpu().tolist()
    hooks=[model.classifier.register_forward_pre_hook(capture_input),
           model.classifier.layers.register_forward_hook(capture_head)]
    samples=[]
    try:
        with torch.inference_mode():
            for row in rows:
                encoded=predictor.tokenizer.encode(row['protein_sequence'])
                assert sum(encoded.attention_mask)==1370
                token_positions[:]=[next(i for i,(start,end) in enumerate(encoded.offsets)
                                         if start==p-1 and end==p and not encoded.special_tokens_mask[i]) for p in POSITIONS]
                ids=torch.tensor([encoded.ids],device='cuda')
                mask=torch.tensor([encoded.attention_mask],device='cuda')
                output=model(input_ids=ids,attention_mask=mask)
                logits=output.logits[0].float()
                probability=torch.softmax(logits,dim=-1).cpu().numpy()
                centered=(logits-logits.mean(dim=-1,keepdim=True)).cpu().numpy()
                score=[score_candidate_pam(probability,pam).specificity_adjusted_score for pam in original['pams']]
                samples.append({'variant':row['protein_id'],'sequence_sha256':sha(row['protein_sequence'].encode()),
                    'probability':probability.tolist(),'centered_logits':centered.tolist(),
                    'scores':score,**capture})
                print(args.model,args.precision,len(samples),len(rows),flush=True)
                del output,logits,probability,centered,ids,mask
    finally:
        for hook in hooks: hook.remove()
    wt=samples[0]; assert wt['variant']=='DSGERT'
    targets={v['variant']:v['clipped_target'] for v in original['variants']}
    for sample in samples:
        target=targets[sample['variant']]
        old=original['score_vectors'][args.model][sample['variant']]
        sample['changes_from_wt']={key:vector_change(sample[key],wt[key]) for key in
                                   ['six_site_hidden','cls','head_hidden','centered_logits','probability']}
        rho=spearman(sample['scores'],target); control=spearman(wt['scores'],target)
        sample['metrics']={'spearman':rho,'wt_spearman':control,
            'delta_spearman':rho-control if rho is not None and control is not None else None,
            'top5_overlap':overlap(sample['scores'],target),
            'delta_top5_overlap':overlap(sample['scores'],target)-overlap(wt['scores'],target),
            'rank_correlation_to_wt':spearman(sample['scores'],wt['scores']),
            'rank_correlation_to_original_bf16':spearman(sample['scores'],old),
            'max_score_difference_from_original_bf16':float(np.max(np.abs(np.asarray(sample['scores'])-old)))}
    gain=[s['metrics']['delta_spearman'] for s in samples[1:] if s['metrics']['delta_spearman'] is not None]
    summary={'mean_delta_spearman':float(np.mean(gain)) if gain else None,
        'defined_non_wt_pairs':len(gain),
        'mean_delta_top5_overlap':float(np.mean([s['metrics']['delta_top5_overlap'] for s in samples[1:]])),
        'min_rank_correlation_to_wt':min(s['metrics']['rank_correlation_to_wt'] for s in samples[1:]),
        'representation_mean_relative_l2':{k:float(np.mean([s['changes_from_wt'][k]['relative_l2'] for s in samples[1:]]))
            for k in ['six_site_hidden','cls','head_hidden','centered_logits','probability']}}
    save_json(destination,{'model':args.model,'precision':args.precision,'tf32':False,
        'torch_version':torch.__version__,'gpu':torch.cuda.get_device_name(),
        'peak_allocated_bytes':torch.cuda.max_memory_allocated(),'bench_sha256':inputs['bench_sha256'],
        'original_results_sha256':sha((OUT/'results.json').read_bytes()),
        'script_sha256':sha(Path(__file__).read_bytes()),
        'checkpoint_sha256':inputs['checkpoint_sha256'] if args.model=='F0_step20' else None,
        'samples':samples,'summary':summary})
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
