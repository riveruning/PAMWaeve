"""Bounded Cas9 report workflow. Run --demo without model weights or network."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pamdict.delivery import make_report, write_report
from pamdict.score.spectrum import valid_pam


def validate_protein(path):
    records = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if not line[1:].strip():
                raise ValueError("蛋白 FASTA 标题不能为空")
            records.append('')
        elif not records:
            raise ValueError("蛋白序列必须位于 FASTA 标题之后")
        else:
            records[-1] += ''.join(line.split()).upper()
    if len(records) != 1 or not records[0]:
        raise ValueError("本入口一次只接受一条非空完整 Cas9 蛋白序列")
    if set(records[0]) - set('ACDEFGHIKLMNPQRSTVWY'):
        raise ValueError("仅接受 20 种标准氨基酸；请先处理终止符或未知残基")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument('--demo', action='store_true', help='合成数据工程演示，不加载模型')
    source.add_argument('--protein', type=Path, help='单条完整 Cas9 FASTA')
    ap.add_argument('--candidates', default='NGG,NAG,NGA,NRG', help='同长度 IUPAC 候选；默认仅为示例，不代表所有 Cas9')
    ap.add_argument('--auto-candidates', action='store_true', help='只报告模型共识；不用默认三碱基候选')
    ap.add_argument('--spacers', type=Path)
    ap.add_argument('--targets', type=Path, action='append')
    ap.add_argument('--reverse-spacers', action='store_true')
    ap.add_argument('--max-mismatches', type=int, default=2)
    ap.add_argument('--device', choices=['cpu', 'cuda'], default='cpu')
    ap.add_argument('--timeout', type=int, default=600, help='每个计算阶段最多秒数')
    ap.add_argument('--outdir', type=Path, required=True, help='必须是尚不存在的新目录')
    args = ap.parse_args()
    motifs = list(dict.fromkeys(x.strip().upper() for x in args.candidates.split(',') if x.strip()))
    try:
        if args.demo and args.auto_candidates:
            raise ValueError('合成演示不支持自动模型共识')
        if not motifs or any(not valid_pam(m) or not 1 <= len(m) <= 10 for m in motifs):
            raise ValueError('候选必须是长度 1–10 的 IUPAC PAM')
        if len({len(m) for m in motifs}) != 1:
            raise ValueError('本入口要求候选 PAM 长度一致；不同长度请分别运行')
        if bool(args.spacers) != bool(args.targets):
            raise ValueError('--spacers 和 --targets 必须一起提供')
        if args.demo and (args.spacers or args.reverse_spacers):
            raise ValueError('演示模式不能混用真实 spacer 输入或反向参数')
        if args.reverse_spacers and not args.spacers:
            raise ValueError('--reverse-spacers 需要 --spacers')
        if args.timeout <= 0 or args.max_mismatches < 0:
            raise ValueError('timeout 必须为正，max-mismatches 不能为负')
        paths = [p.resolve() for p in [args.protein, args.spacers, *(args.targets or [])] if p]
        for p in paths:
            if not p.is_file():
                raise ValueError(f'输入文件不存在：{p}')
        if args.protein:
            validate_protein(args.protein)
        out = args.outdir.resolve()
        out.mkdir(parents=True, exist_ok=False)
    except (ValueError, OSError) as exc:
        ap.error(str(exc))

    env = os.environ.copy()
    env['PYTHONPATH'] = os.pathsep.join([str(ROOT / '.pylibs'), str(ROOT / '.reference/Protein2PAM'), str(ROOT / 'src'), str(ROOT), env.get('PYTHONPATH', '')])
    env.setdefault('HF_HOME', str(ROOT / 'data/checkpoints/hf'))
    env['HF_HUB_OFFLINE'] = '1'
    env['TRANSFORMERS_OFFLINE'] = '1'
    provenance = {'started_utc': datetime.now(timezone.utc).isoformat(), 'python': sys.executable,
                  'arguments': {k: str(v) for k, v in vars(args).items()},
                  'input_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  'commands': [], 'status': 'running'}
    provenance['source_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in [Path(__file__), ROOT/'src/pamdict/delivery.py',
                                            ROOT/'src/pamdict/score/candidate.py',
                                            ROOT/'src/pamdict/score/spacer.py', ROOT/'src/pamdict/infer/p2pam.py']}
    def run(script, options):
        cmd = [sys.executable, '-u', str(ROOT / 'scripts' / script), *options]
        provenance['commands'].append(cmd)
        print('运行阶段：'+script, flush=True)
        with (out / (script+'.log')).open('w', encoding='utf-8') as log:
            process = subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = process.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                # Kill this stage's whole process group so no grandchildren survive.
                try:
                    pgid = os.getpgid(process.pid)
                except ProcessLookupError:
                    pgid = None
                if pgid is not None:
                    for sig in (signal.SIGTERM, signal.SIGKILL):
                        try:
                            os.killpg(pgid, sig)
                        except ProcessLookupError:
                            break
                        try:
                            process.wait(timeout=10)
                            break
                        except subprocess.TimeoutExpired:
                            continue
                raise TimeoutError(f'{script} 超过 {args.timeout}s，其进程组已终止')
            if code != 0:
                raise RuntimeError(f'{script} 退出码 {code}；见 {script}.log')
    try:
        spacer = None
        if args.demo:
            from pamdict.score.candidate import rank_candidate_pams
            matrix = [[.25]*4, [.02,.02,.94,.02], [.02,.02,.94,.02]] + [[.25]*4 for _ in range(7)]
            protein = {'model': 'synthetic_fixture_NOT_Protein2PAM',
                       'proteins': [{'protein_id':'synthetic_demo', 'predicted_pam':'NGG', 'probability_matrix':matrix}],
                       'scores':[s.to_dict() for s in rank_candidate_pams(matrix, motifs, side='downstream')]}
            (out/'protein.json').write_text(json.dumps(protein, indent=2), encoding='utf-8')
        else:
            candidate_args = ['--auto-candidates'] if args.auto_candidates else ['--candidates', ','.join(motifs)]
            run('score_candidate_pams.py', ['--protein', str(args.protein.resolve()), '--model', 'cas9_full',
                *candidate_args, '--device', args.device, '--train-fasta', '',
                '--out', str(out/'protein.tsv'), '--json-out', str(out/'protein.json')])
            protein = json.loads((out/'protein.json').read_text())
            if args.auto_candidates:
                motifs = [r['candidate_pam'] for r in protein['scores']]
                if not motifs:
                    print('模型未给出明确 PAM（无共识）。不生成候选表，结果标记为“未得到明确 PAM”。',
                          file=sys.stderr)
        if args.spacers and not motifs:
            raise ValueError(
                'spacer 证据需要明确的候选 PAM 长度，但本次未得到明确 PAM；'
                '请用 --candidates 指定候选集合，或省略 --spacers 只做蛋白单路输出'
            )
        if args.spacers:
            opts = ['--spacers', str(args.spacers.resolve()), '--pam-side', 'downstream',
                    '--pam-length', str(len(motifs[0])), '--max-mismatches', str(args.max_mismatches),
                    '--candidates', ','.join(motifs), '--outdir', str(out/'spacer')]
            for p in args.targets:
                opts.extend(['--targets', str(p.resolve())])
            if args.reverse_spacers:
                opts.append('--reverse-spacers')
            run('score_spacer_pams.py', opts)
            spacer = json.loads((out/'spacer/summary.json').read_text())
        provenance['status'] = 'complete'
        write_report(make_report(protein, spacer, mode='synthetic_demo' if args.demo else 'real_inference', provenance=provenance), out)
    except Exception as exc:
        provenance.update(status='failed', error=str(exc))
        (out/'run.json').write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding='utf-8')
        print(f'运行失败，未完成报告。请查看 {out} 中的 run.json 和日志：{exc}', file=sys.stderr)
        return 1
    (out/'run.json').write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'报告：{out / "report.html"}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
