"""Offline delivery checks."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from pamdict.delivery import evidence_status, make_report, write_report

ROOT = Path(__file__).resolve().parents[4]


def test_delivery_evidence_states():
    p = {'specificity_adjusted_score': 90}
    s = {'effective_spacer_hits': 3, 'spacer_evidence_score': 90}
    assert '未提供' in evidence_status(p, None)
    assert '均支持' in evidence_status(p, s)
    assert '方向警告' in evidence_status(p, s, orientation_warning=True)
    assert '不足' in evidence_status(p, dict(s, effective_spacer_hits=0))
    assert '冲突' in evidence_status(p, dict(s, spacer_evidence_score=20))
    assert '中性' in evidence_status(p, dict(s, spacer_evidence_score=50))


def test_delivery_report_escapes_user_text():
    p = {'model': '<script>bad()</script>', 'proteins': [],
         'scores': [{'candidate_pam': 'NGG', 'specificity_adjusted_score': 90}]}
    with tempfile.TemporaryDirectory() as d:
        report = make_report(p, mode='real_inference', provenance={'user':'<script>'})
        write_report(report, d)
        page = (Path(d)/'report.html').read_text()
        assert '<script>' not in page
        assert '&lt;script&gt;' in page
        assert report['training_exposure'] == 'unknown'
        assert report['candidates'][0]['spacer_score'] is None


def test_delivery_cli_demo_and_refuse_overwrite():
    with tempfile.TemporaryDirectory() as d:
        out = Path(d)/'report'
        cmd = [sys.executable, str(ROOT/'scripts/run_pam_report.py'), '--demo', '--outdir', str(out)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stderr
        report = json.loads((out/'report.json').read_text())
        assert report['mode'] == 'synthetic_demo'
        assert report['model'] == 'synthetic_fixture_NOT_Protein2PAM'
        assert '没有运行蛋白模型' in (out/'report.html').read_text()
        before = (out/'report.json').read_bytes()
        assert subprocess.run(cmd, capture_output=True, timeout=30).returncode != 0
        assert before == (out/'report.json').read_bytes()


def test_delivery_cli_invalid_inputs():
    with tempfile.TemporaryDirectory() as d:
        for extras in [['--candidates','NGG,N'], ['--candidates','XYZ'],
                       ['--spacers','missing.fna'], ['--timeout','0'], ['--reverse-spacers']]:
            out = Path(d)/'invalid'
            r = subprocess.run([sys.executable,str(ROOT/'scripts/run_pam_report.py'),
                                '--demo','--outdir',str(out),*extras],capture_output=True,timeout=30)
            assert r.returncode != 0
            assert not out.exists()


def test_delivery_cli_multiple_proteins_rejected():
    with tempfile.TemporaryDirectory() as d:
        fasta = Path(d)/'proteins.faa'
        fasta.write_text('>a\nMALW\n>b\nMALW\n')
        r = subprocess.run([sys.executable,str(ROOT/'scripts/run_pam_report.py'),
                            '--protein',str(fasta),'--outdir',str(Path(d)/'out')],
                           capture_output=True,text=True,timeout=30)
        assert r.returncode != 0
        assert '一次只接受一条' in r.stderr
