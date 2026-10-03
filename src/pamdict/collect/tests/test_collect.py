"""Offline tests for the collect package (no network).

Run:  PYTHONPATH=src pytest src/pamdict/collect/tests -q
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from pamdict.collect.discover import (  # noqa: E402
    Candidate,
    classify,
    deduplicate,
    _norm_title,
)
from pamdict.collect.parse import (  # noqa: E402
    extract_pam_motifs,
    extract_uniprot,
    infer_cas_family,
    infer_crispr_type,
    parse_fulltext_xml,
    parse_supp_table,
    sanitize_iupac,
)
from pamdict.schema.record import pam_matrix_from_logo, pam_logo_from_matrix  # noqa: E402
from pamdict.collect.parse_readcount import (  # noqa: E402
    parse_readcount,
    parse_readcount_weighted,
    readcount_freq_to_matrix,
)
from pamdict.collect.parse_pam_table import parse_pam_table  # noqa: E402
from pamdict.collect.parse_position_matrix import parse_position_matrix  # noqa: E402


def make_record(**kw):
    base = {
        "id": "1", "title": "Title", "authorString": "A B", "pubYear": "2025",
        "journalInfo": {"journal": {"title": "J"}}, "doi": "10.1/x",
        "pmid": "", "pmcid": "", "abstractText": "PAM depletion assay",
        "fullTextIdList": {"fullTextId": ["x"]},
    }
    base.update(kw)
    return base


def test_norm_title():
    assert _norm_title("Hello, World!") == "helloworld"


def test_classify_experimental():
    c = classify(make_record(abstractText="PAM depletion assay for Cas9"), ["type_ii"])
    assert c.evidence_status == "candidate_experimental"


def test_classify_review_takes_priority():
    c = classify(make_record(abstractText="review of CRISPR PAM depletion"), ["type_i"])
    assert c.evidence_status == "review_or_secondary"


def test_classify_computational_only():
    c = classify(make_record(abstractText="deep learning prediction of PAM"), ["type_v"])
    assert c.evidence_status == "computational_only"


def test_dedup_by_doi():
    a = classify(make_record(id="1", doi="10.0/a", title="First"), ["type_i"])
    b = classify(make_record(id="2", doi="10.0/A", title="First Two"), ["type_v"])
    out = deduplicate([a, b])
    assert len(out) == 1
    assert "type_v" in out[0].inferred_crispr_types  # merged


def test_dedup_by_title_key():
    a = classify(make_record(id="1", doi="", title="Same Paper"), ["type_i"])
    b = classify(make_record(id="2", doi="", title="Same Paper!!"), ["type_ii"])
    assert len(deduplicate([a, b])) == 1


def test_sanitize_iupac():
    assert sanitize_iupac("Ngg") == "NGG"
    assert sanitize_iupac("5'TTCN3'") == "TTCN"


def test_extract_pam_motifs():
    text = "The PAM is NGG for this Cas9. Another system has TTTN PAM."
    motifs = extract_pam_motifs(text)
    assert "NGG" in motifs and "TTTN" in motifs


def test_extract_uniprot():
    assert "P12345" in extract_uniprot("UniProt P12345 and Q9XYZ1")


def test_infer_family_type():
    assert infer_cas_family("Cas12a cleaves") == "Cas12a"
    assert infer_crispr_type("Type I-F system") == "Type I"


def test_parse_supp_table():
    rows = [
        {"pam": "NGG", "organism": "Streptococcus pyogenes", "accession": "P0D2D1"},
        {"pam": "TTN", "organism": "Pyrococcus furiosus"},
    ]
    samples = parse_supp_table(rows, doi="10.0/x", citation="Paper 2024")
    assert len(samples) == 2
    assert samples[0].pam_consensus == "NGG"
    assert samples[0].organism == "Streptococcus pyogenes"
    assert samples[0].method == "supp_table"


def test_schema_roundtrip():
    m = [[0.0] * 4 for _ in range(10)]
    m[0][2] = 1.2
    logo = pam_logo_from_matrix(m)
    back = pam_matrix_from_logo(logo)
    assert back[0][2] == 1.2 and len(back) == 10 and all(len(r) == 4 for r in back)


def _make_readcount_rows():
    """Minimal HT-PAMDA table: 3-nt library, Bound (selected) channel."""
    header = [
        ["PAM", "Read Counts", "", "", "Normalized Read Counts"],
        ["", "#1", "", "#2", ""],
        ["", "Unbound", "Bound", "Unbound", "Bound"],
    ]
    # positions 1..3 all "G" only -> pure G logo at each position.
    data = [[ "GGG", "100", "900", "100", "900"]]
    return header + data


def test_parse_readcount_positive_channel():
    rows = _make_readcount_rows()
    freq = parse_readcount(rows)
    assert set(freq) == {"A", "C", "G", "T"}
    # Only G contributes -> freq G == 1.0 at every position, others 0.
    for b in "ACT":
        assert all(abs(v) < 1e-9 for v in freq[b]), f"{b} should be ~0: {freq[b]}"
    assert all(abs(v - 1.0) < 1e-9 for v in freq["G"])


def test_parse_readcount_mixed_bases_normalize():
    # Two sequences sharing a position with equal weights -> 0.5 / 0.5.
    header = [
        ["PAM", "Read Counts", "", ""],
        ["", "#1", "", ""],
        ["", "Unbound", "Bound", "Unbound", "Bound"],
    ]
    # AGT and CGT: pos1 A=100, C=100; pos2 G=200; pos3 T=200.
    data = [
        ["AGT", "50", "100", "50", "100"],
        ["CGT", "50", "100", "50", "100"],
    ]
    freq = parse_readcount(header + data)
    assert abs(freq["A"][0] - 0.5) < 1e-9
    assert abs(freq["C"][0] - 0.5) < 1e-9
    assert abs(freq["G"][1] - 1.0) < 1e-9
    assert abs(freq["T"][2] - 1.0) < 1e-9


def test_readcount_rejects_non_table():
    assert parse_readcount([["foo", "bar"], ["x", "y"]]) == {}


def test_readcount_freq_to_matrix_shape():
    freq = {"A": [1.0], "C": [0.0], "G": [0.0], "T": [0.0]}
    m = readcount_freq_to_matrix(freq)
    assert len(m) == 10 and all(len(r) == 4 for r in m)
    assert m[0][0] == 1.0  # pos1 A == 1.0
    assert all(m[i][0] == 0.0 for i in range(1, 10))  # padded


def test_readcount_mean_enrichment_preferred():
    # A table with an Enrichment 'Mean' column: the mean ratio should win over
    # raw counts as the weight.  Here AAT has huge Mean enrichment, so pos2 A
    # should dominate; raw counts would have made a different call.
    header = [
        ["PAM", "Read Counts", "Enrichment", "Mean"],
        ["", "Unbound", "Bound", ""],
    ]
    data = [
        ["AAT", "10", "5", "100.0"],   # strong enrichment -> pos2 A
        ["ACT", "100", "90", "1.0"],   # weak enrichment
        ["AGT", "10", "5", "1.0"],
        ["ATT", "10", "5", "1.0"],
    ]
    freq, mode = parse_readcount_weighted(header + data)
    assert mode == "mean_enrichment"
    # pos2: A carries 100 vs others ~1 -> A ~ ~0.97
    assert freq["A"][1] > 0.9
    # pos1: AAT(100) vs ACT(1)+AGT(1)+ATT(1) -> still A-heavy but we only assert
    # pos2 cleanly.
    assert freq["A"][0] > freq["C"][0]


def test_readcount_selected_counts_fallback():
    # No Mean column -> fall back to summed Bound counts, mode selected_counts.
    rows = _make_readcount_rows()
    freq, mode = parse_readcount_weighted(rows)
    assert mode == "selected_counts"
    assert all(abs(v - 1.0) < 1e-9 for v in freq["G"])


def test_parse_pam_table_basic():
    rows = [
        ["Target ID", "Species", "Description", "PAM"],
        ["8N", "N/A", "depletion", "NNNNNNNN"],
        ["EMX1", "H. sapiens", "targeting", "ATG"],
        ["Amp", "N/A", "interference", "ATG"],
    ]
    recs = parse_pam_table(rows)
    # The all-N control is filtered out.
    assert len(recs) == 2
    assert all(r["pam_consensus"] == "ATG" for r in recs)
    assert recs[0]["raw_ref"].startswith("org=")


def test_parse_pam_table_with_protein_sequence():
    seq = "M" + "A" * 60  # 61 aa
    rows = [
        ["Name", "Protein sequence", "PAM"],
        ["CoCas9", seq, "NGAG"],
    ]
    recs = parse_pam_table(rows)
    assert len(recs) == 1
    assert recs[0]["pam_consensus"] == "NGAG"
    assert recs[0]["protein_sequence"] == seq
    assert recs[0]["protein_id"] == "CoCas9"


def test_parse_position_matrix_nrg():
    # 3 positions; pos2 & pos3 are pure G (like NGG), pos1 uniform.
    rows = [
        ["b\\p", "1", "2", "3"],
        ["A", "100", "0", "0"],
        ["C", "100", "0", "0"],
        ["G", "100", "100", "100"],
        ["T", "100", "0", "0"],
    ]
    m = parse_position_matrix(rows)
    assert m and len(m["A"]) == 3
    # pos1 uniform -> ~0.25 each
    assert abs(m["A"][0] - 0.25) < 1e-9
    # pos2 pure G
    assert abs(m["G"][1] - 1.0) < 1e-9
    assert abs(m["G"][2] - 1.0) < 1e-9
    assert m["C"][2] == 0.0


def test_parse_position_matrix_rejects():
    assert parse_position_matrix([["foo", "bar"], ["1", "2"]]) == {}


def test_clean_pam_motif_mutation_annotation():
    from pamdict.collect.parse_pam_table import _clean_pam_motif
    # (X>Y) mutation annotations are stripped; IUPAC ambiguity codes kept.
    assert _clean_pam_motif("NR(A>G)TTTT") == "NRTTTT"
    assert _clean_pam_motif("NAR(G>A)H(W>C)H(A>T>C)GN(C>T>R)") == "NARHHGN"
    assert _clean_pam_motif("N(C>D)M(A>C)RN(A>B)AY(C>T)") == "NMRNAY"
    # quotes stripped, U->T
    assert _clean_pam_motif('"NNNNGAAA"') == "NNNNGAAA"
    # full IUPAC ambiguity codes preserved
    assert _clean_pam_motif("BRTTTTT") == "BRTTTTT"
    assert _clean_pam_motif("NRRWC") == "NRRWC"


def test_clean_pam_motif_rejects_empty():
    from pamdict.collect.parse_pam_table import _clean_pam_motif, _is_iupac_pam
    assert _clean_pam_motif("-") == ""
    assert not _is_iupac_pam("-")
    assert not _is_iupac_pam("NNNNNNNN")  # all-N control


def test_parse_pam_table_complex_quadruple():
    # CoCas9-style single-table quadruple: ortholog / organism / PAM / AA seq.
    seq = "M" + "A" * 80
    rows = [
        ["Ortholog", "Subtype", "Organism", "PAM", "AA sequence"],
        ["CoCas9", "II-C", "Capnocytophaga ochracea DSM 7271", "NR(A>G)TTTT", seq],
    ]
    recs = parse_pam_table(rows)
    assert len(recs) == 1
    r = recs[0]
    assert r["pam_consensus"] == "NRTTTT"
    assert r["protein_sequence"] == seq
    assert r["protein_id"] == "CoCas9"
    assert "ochracea" in r["raw_ref"]


def test_parse_pam_table_truncated_header_and_degenerate_consensus():
    # IscB-style: "Amino acid sequenc" (truncated), "Degenerate_consensus" (PAM),
    # "Key name" (protein id).  Also a leading title row + empty cells must not
    # cause the empty-header v. empty-synonym prefix-match bug.
    seq = "M" + "A" * 90
    rows = [
        ["Supplementary table 1", "", "", ""],
        ["", "", "", ""],
        ["Key name", "Amino acid sequenc", "Degenerate_consensus", "scaffold"],
        ["IscB.m1", seq, "AGARGA", "GGCT..."],
    ]
    recs = parse_pam_table(rows)
    assert len(recs) == 1
    r = recs[0]
    assert r["pam_consensus"] == "AGARGA"
    assert r["protein_sequence"] == seq
    assert r["protein_id"] == "IscB.m1"


def test_find_col_skips_empty_headers():
    from pamdict.collect.parse_pam_table import _find_col, _PAM_HEADERS
    # Empty header cells must not prefix-match any synonym; "PAM" exact matches.
    assert _find_col(["", "", "", "PAM"], _PAM_HEADERS) == 3
    assert _find_col(["supplementary table", "", "", ""], _PAM_HEADERS) is None
    # "PAM variant median count" must NOT match (only "degenerate consensus" is
    # prefix-eligible; "pam" is exact-only).
    assert _find_col(["PAM variant median count", "789"], _PAM_HEADERS) is None
    # "Degenerate consensus TAM (if functional)" IS prefix-eligible.
    assert _find_col(["Key", "Amino acid sequenc", "Degenerate_consensus TAM"],
                     _PAM_HEADERS, prefix_synonyms=("degenerate consensus",)) == 2
    # A short data value like "DE" (Direction column) must NOT match
    # "degenerate consensus" by the reverse prefix direction.
    assert _find_col(["SYMBOL", "Status", "Direction"],
                     _PAM_HEADERS, prefix_synonyms=("degenerate consensus",)) is None


def test_parse_pam_table_rejects_gene_named_pam():
    # A GSEA-style table where a *data* cell is the gene name "Pam" next to a
    # column of numbers must NOT be treated as a PAM table.
    rows = [
        ["Fig.5a", "", ""],
        ["", "NAME", "SYMBOL", "RANK"],
        ["", "row_0", "St8sia5", "42"],
        ["", "row_1", "Far2", "152"],
        ["", "row_2", "Pam", "37"],
    ]
    assert parse_pam_table(rows) == []
