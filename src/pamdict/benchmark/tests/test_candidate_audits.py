"""Regression checks for candidate evidence, without network or training data."""
import io
import zipfile

from scripts.audit_pp_df_training_neighbors import audit
from scripts.extract_pp_df_pam_evidence import extract


def test_exhaustive_audit_includes_identity_boundary_and_deduplicates():
    query = 'A' * 100
    rows = [{'protein_sequence': query, 'protein_id': 'first'},
            {'protein_sequence': query, 'protein_id': 'duplicate'},
            {'protein_sequence': 'A' * 90 + 'C' * 10},
            {'protein_sequence': 'A' * 85 + 'C' * 15},
            {'protein_sequence': 'A' * 84}]
    result = audit(query, rows)
    assert result['exact']
    assert result['unique_sequences'] == 4
    assert result['length_eligible_sequences_checked'] == 3
    assert result['near90_count'] == 2
    assert result['near85_count'] == 3
    assert len(result['matches_at_or_above_85'][0]['records']) == 2


def test_exhaustive_audit_no_hit_does_not_invent_nearest_neighbor():
    result = audit('A' * 100, [{'protein_sequence': 'C' * 100}])
    assert not result['exact']
    assert result['matches_at_or_above_85'] == []
    assert result['length_eligible_sequences_checked'] == 1


def workbook(pams):
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    strings = ['DfCas9 depleted PAM'] + pams
    shared = f'<sst xmlns="{ns}">' + ''.join(
        f'<si><t>{s}</t></si>' for s in strings) + '</sst>'
    sheet = f'<worksheet xmlns="{ns}"><sheetData>' + ''.join(
        f'<row><c r="A{i+1}" t="s"><v>{i}</v></c></row>'
        for i in range(len(strings))) + '</sheetData></worksheet>'
    data = io.BytesIO()
    with zipfile.ZipFile(data, 'w') as archive:
        archive.writestr('xl/sharedStrings.xml', shared)
        archive.writestr('xl/worksheets/sheet1.xml', sheet)
    return data.getvalue()


def test_pam_attachment_extraction_preserves_cell_provenance():
    header, cells = extract(workbook(['AAAAAAA', 'ACGTACG']), 2)
    assert header == 'DfCas9 depleted PAM'
    assert cells == [('A2', 'AAAAAAA'), ('A3', 'ACGTACG')]


def test_pam_attachment_rejects_duplicates_bad_alphabet_and_count():
    for pams, count in [(['AAAAAAA', 'AAAAAAA'], 2),
                        (['AAAAAAN'], 1), (['AAAAAAA'], 2)]:
        try:
            extract(workbook(pams), count)
        except AssertionError:
            continue
        raise AssertionError('Invalid attachment was accepted')
