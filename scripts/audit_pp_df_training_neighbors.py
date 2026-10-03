"""Audit Pp/Df leads; exhaustive global edit checks after safe length filter.

Outputs are candidate evidence only, never automatic benchmark registration.
"""
import csv
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from pamdict.finetune.similarity import global_edit_identity

URL = ('https://rest.uniprot.org/uniprotkb/search?'
       'query=xref%3Arefseq-WP_081825100.1&format=fasta')


def rows(path):
    with path.open(encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def audit(query, reference_rows, threshold=0.85):
    grouped = {}
    for row in reference_rows:
        seq = row.get('protein_sequence', '').strip().upper()
        if seq:
            grouped.setdefault(seq, []).append(row)
    matches = []
    checked = 0
    for seq, records in grouped.items():
        if min(len(query), len(seq)) / max(len(query), len(seq)) < threshold:
            continue  # Length difference alone exceeds permitted edit count.
        checked += 1
        identity = global_edit_identity(query, seq, threshold)
        if identity is not None:
            matches.append({'identity': identity, 'length': len(seq),
                            'sequence_sha256': sha(seq.encode()),
                            'records': [{'protein_id': r.get('protein_id', ''),
                                         'pam': r.get('pam_consensus', ''),
                                         'source': r.get('source', '')} for r in records]})
    matches.sort(key=lambda m: (-m['identity'], m['sequence_sha256']))
    return {'reference_rows': len(reference_rows), 'unique_sequences': len(grouped),
            'length_eligible_sequences_checked': checked,
            'exact': query in grouped,
            'near90_count': sum(m['identity'] >= .90 for m in matches),
            'near85_count': len(matches), 'matches_at_or_above_85': matches}


def main():
    raw = ROOT / 'data/raw/dfcas9_uniprot_WP_081825100.faa'
    if raw.exists():
        data = raw.read_bytes()
    else:
        data = urllib.request.urlopen(URL, timeout=20).read()
        if data.count(b'>') != 1:
            raise ValueError('Expected one UniProt record; re-audit accession mapping')
        raw.write_bytes(data)
    text = data.decode()
    if text.count('>') != 1 or not text.startswith('>'):
        raise ValueError('Invalid single-record FASTA')
    dfseq = ''.join(text.splitlines()[1:]).strip().upper()
    if len(dfseq) != 1079 or set(dfseq) - set('ACDEFGHIKLMNPQRSTVWY'):
        raise ValueError('Unexpected Df lead length/alphabet')
    gold = rows(ROOT / 'data/corpus/gold_master.tsv')
    queries = {'DfCas9_UniProt_lead': dfseq}
    for r in gold:
        if r['protein_id'] == 'A0ABY6H3A4':
            seq = r['protein_sequence'].strip().upper()
            queries[f'legacy_A0ABY6H3A4_length{len(seq)}'] = seq
    paths = ['data/raw/protein2pam_train_seqs.tsv',
             'data/corpus/augmented_train_v2_cas9_train.tsv']
    report = {'timestamp_utc': datetime.now(timezone.utc).isoformat(),
              'identity_definition': '1-global_levenshtein_distance/max_length',
              'search': 'exhaustive unique sequences after safe length filter; no kmer candidate cap',
              'df_source_url': URL, 'df_source_sha256': sha(data),
              'df_header': text.splitlines()[0],
              'df_provenance_status': 'database lead; exact experimental construct link pending',
              'reference_files': {}, 'candidates': {}}
    for path in paths:
        p = ROOT / path
        refs = rows(p)
        if path.startswith('data/raw'):
            refs = [r for r in refs if r.get('cas_family', '').upper().startswith('CAS9')
                    or r.get('crispr_type', '').replace(' ', '').upper() in {'II', 'TYPEII'}]
        report['reference_files'][path] = {'sha256': sha(p.read_bytes()), 'rows': len(refs)}
        for name, seq in queries.items():
            print(f'audit {name} against {path}', flush=True)
            candidate = report['candidates'].setdefault(name, {
                'length': len(seq), 'sequence_sha256': sha(seq.encode()), 'references': {}})
            candidate['references'][path] = audit(seq, refs)
    out = ROOT / 'data/parsed/pp_df_training_neighbors.json'
    out.write_text(json.dumps(report, indent=2) + '\n')
    for name, result in report['candidates'].items():
        for path, found in result['references'].items():
            print(name, path, 'exact', found['exact'], 'near90', found['near90_count'],
                  'best', found['matches_at_or_above_85'][:1])
    print(out)


if __name__ == '__main__':
    main()
