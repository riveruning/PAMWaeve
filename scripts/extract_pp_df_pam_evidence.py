"""Extract published positive PAM lists; NOT benchmark admission or activity data."""
import csv
import hashlib
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def extract(data, expected):
    with zipfile.ZipFile(io.BytesIO(data)) as book:
        strings = [''.join(n.itertext()) for n in
                   ET.fromstring(book.read('xl/sharedStrings.xml')).findall('m:si', NS)]
        sheet = ET.fromstring(book.read('xl/worksheets/sheet1.xml'))
        cells = []
        header = None
        for cell in sheet.findall('.//m:sheetData/m:row/m:c', NS):
            value = cell.find('m:v', NS)
            if value is None:
                continue
            assert cell.get('t') == 's', 'Unexpected cell encoding'
            text = strings[int(value.text)].strip()
            if cell.get('r') == 'A1':
                header = text
            else:
                assert len(text) == 7 and set(text) <= set('ACGT'), text
                cells.append((cell.get('r'), text))
        assert len(cells) == expected
        assert len({pam for _, pam in cells}) == expected
        return header, cells


def main():
    source = ROOT / 'data/raw/supp_char/PMC7708072.zip'
    raw = source.read_bytes()
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
        nested = outer.read('gkaa998_supplemental_files.zip')
    records, evidence = [], []
    with zipfile.ZipFile(io.BytesIO(nested)) as supplements:
        for system, number, count, assay in [
            ('DfCas9', 2, 122, 'E_coli_plasmid_transformation_interference'),
            ('PpCas9', 3, 79, 'in_vitro_plasmid_library_cleavage')]:
            name = f'Supplementary File S{number}.xlsx'
            data = supplements.read(name)
            header, cells = extract(data, count)
            assert system in header, header
            evidence.append({'member': name, 'sha256': digest(data),
                             'header': header, 'unique_pams': len(cells)})
            for cell, pam in cells:
                records.append({'system': system, 'pam': pam, 'assay': assay,
                                'doi': '10.1093/nar/gkaa998', 'member': name,
                                'sheet': 'sheet1', 'cell': cell,
                                'evidence': 'published_significantly_depleted_PAM',
                                'activity_value': '',
                                'admission': 'pending_exact_experimental_protein_link'})
    target = ROOT / 'data/parsed/pp_df_positive_pam_candidates.tsv'
    with target.open('w', newline='', encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)
    report = {'source': str(source.relative_to(ROOT)), 'sha256': digest(raw),
              'nested_zip_sha256': digest(nested), 'members': evidence,
              'rows': len(records), 'systems': 2, 'strict_systems_added': 0,
              'limitations': ['Positive lists only; missing PAM is not a negative.',
                              'No read counts or quantitative activity in these sheets.',
                              'Do not combine the two assay contexts.',
                              'Sequence provenance pending; Df lead fails near90 gate.'],
              'output_sha256': digest(target.read_bytes())}
    (ROOT / 'data/parsed/pp_df_positive_pam_evidence.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
