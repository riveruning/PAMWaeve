"""Resumable FP32 inference on DEVELOPMENT ONLY; never load reserved proteins."""
import argparse
import csv
import hashlib
import io
import itertools
import json
import sys
from pathlib import Path
import torch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'src'))
from scripts.build_pammla_benchmark import DEST
from scripts.prepare_pammla_pilot import sha,save_json,save
from pamdict.infer.p2pam import P2PAMPredictor,_local_snapshot_dir
from pamdict.score.candidate import score_candidate_pam
from scripts.finetune_p2pam import apply_lora

OUT=ROOT/'data/parsed/pammla_development_fp32_v1'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--model',choices=['base','F0_step20'],required=True)
    args=ap.parse_args()
    manifest=json.loads((DEST/'manifest.json').read_text())
    assert sha((DEST/'variants.tsv').read_bytes())==manifest['files']['variants.tsv']
    with (DEST/'variants.tsv').open() as f:
        rows=[r for r in csv.DictReader(f,delimiter='\t') if r['future_task_split']=='development']
    assert len(rows)==670 and sum(r['variant_id']=='DSGERT' for r in rows)==1
    snapshot=Path(_local_snapshot_dir('Profluent-Bio/protein2pam-cas9_full'))
    hasher=hashlib.sha256()
    with (snapshot/'model.safetensors').open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):hasher.update(block)
    checkpoint=ROOT/'data/checkpoints/finetune_v3_cas9_F0_equal/checkpoint_step000020.pt'
    context={'model':args.model,'precision':'fp32','tf32':False,'scope':'development_only',
        'variants':len(rows),'reserved_variants_predicted':0,
        'manifest_sha256':sha((DEST/'manifest.json').read_bytes()),
        'base_snapshot':snapshot.name,'base_weights_sha256':hasher.hexdigest(),
        'checkpoint_sha256':sha(checkpoint.read_bytes()) if args.model=='F0_step20' else None,
        'tokenizer_sha256':sha((ROOT/'.reference/Protein2PAM/protein2pam/huggingface/tokenizer.json').read_bytes()),
        'code_sha256':{p:sha((ROOT/p).read_bytes()) for p in ['scripts/predict_pammla_development.py',
            'src/pamdict/infer/p2pam.py','src/pamdict/score/candidate.py',
            '.reference/Protein2PAM/protein2pam/huggingface/modeling_esm.py']},
        'torch_version':torch.__version__}
    assert context['base_weights_sha256']=='aa1dc5c017ddd885bf26c8126675075af0007d1af3ab1a2893e62aaffe0341c2'
    cache=OUT/args.model;save_json(cache/'run_context.json',context)
    run_hash=sha((cache/'run_context.json').read_bytes())
    done={}
    for row in rows:
        path=cache/(row['variant_id']+'.json')
        if path.exists():
            item=json.loads(path.read_text())
            assert item['run_context_sha256']==run_hash and item['sequence_sha256']==row['sequence_sha256']
            assert item['variant_id']==row['variant_id']
            done[row['variant_id']]=item
    pams=['N'+''.join(p) for p in itertools.product('ACGT',repeat=3)]
    if len(done)<len(rows):
        assert torch.cuda.is_available()
        torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
        torch.manual_seed(0)
        predictor=P2PAMPredictor('cas9_full',device='cpu',local_files_only=True)
        model=predictor.model
        if args.model=='F0_step20':
            expected=json.loads((ROOT/'data/parsed/pammla_pilot_v1/inputs.json').read_text())['checkpoint_sha256']
            assert context['checkpoint_sha256']==expected
            ck=torch.load(checkpoint,map_location='cpu')
            assert ck['format_version']==3 and ck['model_name']=='cas9_full' and ck['classifier_mode']=='freeze'
            apply_lora(model,rank=int(ck['lora_rank']),alpha=int(ck['lora_alpha']))
            state=ck['trainable_state'];current=model.state_dict()
            assert set(state)=={k for k in current if 'lora_' in k}
            assert all(state[k].shape==current[k].shape for k in state)
            assert not model.load_state_dict(state,strict=False).unexpected_keys
            del ck,state,current
        model.to(dtype=torch.float32).to('cuda').eval();predictor._device='cuda'
        assert all(b.dtype==torch.float32 for name,b in model.named_buffers() if name.endswith('inv_freq'))
        print(args.model,'FP32 development only; resume',len(done),'/670',flush=True)
        for row in rows:
            name=row['variant_id']
            if name in done:continue
            seq=row['protein_sequence'];assert sha(seq.encode())==row['sequence_sha256']
            encoded=predictor.tokenizer.encode(seq);assert sum(encoded.attention_mask)==1370
            probability=predictor.predict_probability_matrix([seq])[0]
            scores=[score_candidate_pam(probability,p).specificity_adjusted_score for p in pams]
            item={'variant_id':name,'sequence_sha256':row['sequence_sha256'],
                'run_context_sha256':run_hash,'probability':probability.tolist(),'scores':scores}
            save_json(cache/(name+'.json'),item);done[name]=item
            if len(done)%50==0 or len(done)==670:print(args.model,len(done),'/670',flush=True)
    buf=io.StringIO(newline='');w=csv.DictWriter(buf,fieldnames=['variant_id','pam_group','score'],delimiter='\t')
    w.writeheader()
    for row in rows:
        for pam,score in zip(pams,done[row['variant_id']]['scores']):
            w.writerow({'variant_id':row['variant_id'],'pam_group':pam,'score':score})
    save(OUT/(args.model+'_scores.tsv'),buf.getvalue().encode())
    print('Complete:',args.model,'670 development; 0 reserved',flush=True)


if __name__=='__main__':main()
