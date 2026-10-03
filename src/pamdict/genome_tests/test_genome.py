"""Offline tests for the genome-to-Cas9 adapter.

Fixtures deliberately reproduce the *real* CRISPRCasTyper 1.9.0 / pyrodigal-gv
output shape (headers without ``partial=``, ``genes.tab`` with ``Pos`` last,
``cas_operons_putative.tab`` when no good operon exists).  No network, no model
weights and no CCTyper install are required.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile

from pamdict.genome import (
    check_genome,
    collect_systems,
    fasta,
    literal_list,
    reverse_complement,
    safe_id,
    select_system,
    sha256,
    translate,
    verify_protein,
)

ROOT = Path(__file__).resolve().parents[3]

_STANDARD = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _codons_for(residue):
    return sorted(c for c, r in _STANDARD.items() if r == residue)


def _random_protein(rng, length):
    """A random protein beginning with M so it can be encoded with a start codon."""
    return "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(length - 1))


def _dna_for(protein, rng, start_codon="ATG", stop="TAA"):
    """Encode ``protein`` as DNA.

    ``start_codon`` replaces the first codon, so callers must pass a protein
    whose first residue is consistent with it: 'M' for ATG/GTG/TTG (pyrodigal
    reports M for an alternative start), or the naive residue otherwise.
    """
    if start_codon in {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}:
        assert protein[0] == "M", "alternative start codons are reported as M"
    else:
        assert _STANDARD[start_codon] == protein[0], "start codon encodes a different residue"
    parts = [start_codon]
    for residue in protein[1:]:
        parts.append(rng.choice(_codons_for(residue)))
    parts.append(stop)
    return "".join(parts)


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _make_genome_file(directory, contig="CT1", dna=None):
    """Write a genome FASTA and return its path."""
    dna = dna or ("ACGT" * 50)
    path = Path(directory) / "genome.fna"
    _write(path, f">{contig} test contig\n{dna}\n")
    return path


def _annotation_dir(root, *, cas9_dna, contig, pool_dna, table="cas_operons.tab",
                    prediction="II-A", trusted="True", spacer_count=3,
                    array_start=1, array_end=60, repeat="GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC"):
    """Build a CCTyper-like output directory around one Cas9 gene.

    Returns the annotation directory path.
    """
    directory = Path(root) / "annotation"
    directory.mkdir(parents=True, exist_ok=True)

    protein = translate(cas9_dna).rstrip("*")
    if protein.startswith("V") and cas9_dna[:3] == "GTG":
        protein = "M" + protein[1:]
    gene_start = pool_dna.index(cas9_dna) + 1
    gene_end = gene_start + len(cas9_dna) - 1
    contig_length = len(pool_dna)

    # pyrodigal-gv header: contig_pos # start # end # strand -- NO partial field.
    _write(directory / "proteins.faa",
           f">{contig}_1 # {gene_start} # {gene_end} # 1\n{protein}\n")
    _write(directory / "genes.tab",
           "Contig\tStart\tEnd\tStrand\tPos\n"
           f"{contig}\t{gene_start}\t{gene_end}\t1\t1\n")

    header = ("Contig\tOperon\tStart\tEnd\tPrediction\tComplete_Interference\t"
              "Complete_Adaptation\tBest_type\tBest_score\tGenes\tPositions\t"
              "E-values\tCoverageSeq\tCoverageHMM\tStrand_Interference\tStrand_Adaptation\n")
    row = (f"{contig}\t{contig}@1\t1\t{contig_length}\t{prediction}\t100%\t100%\t"
           f"II-A\t9.0\t['Cas9_1_II', 'Cas1_4_CAS-I-II-III-IV-V-VI']\t[1, 2]\t"
           "['0.00e+00', '1.00e-50']\t[1.0, 0.9]\t[1.0, 0.9]\t1\t1\n")
    _write(directory / table, header + row)

    _write(directory / "CRISPR_Cas.tab",
           "Contig\tOperon\tOperon_Pos\tPrediction\tCRISPRs\tDistances\t"
           "Prediction_Cas\tPrediction_CRISPRs\n"
           f"{contig}\t{contig}@1\t[1, {contig_length}]\t{prediction}\t"
           f"['{contig}_1']\t[100]\t{prediction}\t['{prediction}']\n")

    _write(directory / "crisprs_all.tab",
           "Contig\tCRISPR\tStart\tEnd\tConsensus_repeat\tN_repeats\tRepeat_len\t"
           "Spacer_len_avg\tRepeat_identity\tSpacer_identity\tSpacer_len_sem\t"
           "Trusted\tPrediction\tSubtype\tSubtype_probability\n"
           f"{contig}\t{contig}_1\t{array_start}\t{array_end}\t{repeat}\t3\t"
           f"{len(repeat)}\t32\t98.1\t57.5\t0.0\t{trusted}\tII-A\tII-A\t0.99\n")

    spacers = "".join(f">{contig}_1:{i}\n{'ACGT' * 8}\n" for i in range(1, spacer_count + 1))
    _write(directory / "spacers" / f"{contig}_1.fa", spacers)

    return directory


# --- FASTA validation ----------------------------------------------------

def test_fasta_rejects_empty_file():
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "empty.fna"
        path.write_text("")
        try:
            check_genome(path)
        except ValueError as exc:
            assert "empty" in str(exc).lower()
            return
        raise AssertionError("empty FASTA was accepted")


def test_fasta_rejects_duplicate_identifiers():
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "dup.fna"
        path.write_text(">c1 a\nACGT\n>c1 b\nACGT\n")
        try:
            check_genome(path)
        except ValueError as exc:
            assert "duplicate" in str(exc).lower()
            return
        raise AssertionError("duplicate FASTA identifiers were accepted")


def test_fasta_rejects_bad_symbols_and_sequence_before_header():
    with tempfile.TemporaryDirectory() as d:
        bad = Path(d) / "bad.fna"
        bad.write_text(">c1\nACGTZZZ\n")
        try:
            check_genome(bad)
        except ValueError as exc:
            assert "symbol" in str(exc).lower()
        else:
            raise AssertionError("invalid symbols were accepted")
        early = Path(d) / "early.fna"
        early.write_text("ACGT\n>c1\nACGT\n")
        try:
            check_genome(early)
        except ValueError as exc:
            assert "before" in str(exc).lower()
            return
        raise AssertionError("sequence before header was accepted")


def test_fasta_rejects_empty_record_and_oversize():
    with tempfile.TemporaryDirectory() as d:
        empty = Path(d) / "e.fna"
        empty.write_text(">c1\n>c2\nACGT\n")
        try:
            check_genome(empty)
        except ValueError as exc:
            assert "empty" in str(exc).lower()
        else:
            raise AssertionError("empty record was accepted")
        big = Path(d) / "big.fna"
        big.write_text(">c1\n" + "ACGT" * 20 + "\n")
        try:
            check_genome(big, max_bytes=10)
        except ValueError as exc:
            assert "30 MB" in str(exc) or "up to" in str(exc)
            return
        raise AssertionError("oversize genome was accepted")


def test_check_genome_reports_hash_and_contigs():
    with tempfile.TemporaryDirectory() as d:
        path = _make_genome_file(d, contig="CT1", dna="ACGT" * 30)
        dna = "ACGT" * 30
        summary = check_genome(path)
        assert summary["contigs"] == 1
        assert summary["bases"] == len(dna)
        assert summary["sha256"] == sha256(path)
        assert summary["contig_lengths"] == {"CT1": len(dna)}


# --- translation / integrity --------------------------------------------

def test_verify_protein_accepts_exact_forward_translation():
    rng = random.Random(7)
    protein = _random_protein(rng, 40)
    dna = _dna_for(protein, rng)
    verified, notes = verify_protein(protein, contig_sequence=dna,
                                     start=1, end=len(dna), strand=1)
    assert verified, notes
    assert "translation_matches_exactly" in notes


def test_verify_protein_accepts_reverse_strand_translation():
    rng = random.Random(11)
    protein = _random_protein(rng, 40)
    coding = _dna_for(protein, rng)
    contig = "TTTT" + reverse_complement(coding) + "GGGG"
    start = 5
    end = start + len(coding) - 1
    verified, notes = verify_protein(protein, contig_sequence=contig,
                                     start=start, end=end, strand=-1)
    assert verified, notes
    assert "translation_matches_exactly" in notes


def test_verify_protein_accepts_alternative_start_codon():
    """Real case: pyrodigal reports M for a GTG start; naive table gives V."""
    rng = random.Random(13)
    protein = "M" + _random_protein(rng, 30)[1:]
    coding = _dna_for(protein, rng, start_codon="GTG")
    naive = translate(coding).rstrip("*")
    assert naive[0] == "V" and naive[1:] == protein[1:]
    verified, notes = verify_protein(protein, contig_sequence=coding,
                                     start=1, end=len(coding), strand=1)
    assert verified, notes
    assert "translation_matches_with_alternative_start_codon" in notes


def test_verify_protein_rejects_truncated_and_mismatched():
    rng = random.Random(17)
    protein = _random_protein(rng, 30)
    dna = _dna_for(protein, rng)
    # Truncated call: the reported protein is a prefix only.
    verified, notes = verify_protein(protein[:20], contig_sequence=dna,
                                     start=1, end=len(dna), strand=1)
    assert not verified
    assert "translation_mismatch" in notes
    # Wrong coordinates entirely.
    verified, notes = verify_protein(protein, contig_sequence=dna,
                                     start=1, end=len(dna) - 3, strand=1)
    assert not verified


def test_verify_protein_flags_bad_stop_and_internal_stop_and_range():
    rng = random.Random(19)
    protein = _random_protein(rng, 20)
    coding = _dna_for(protein, rng, stop="TAA")
    no_stop = coding[:-3] + "GGG"
    _, notes = verify_protein(protein, contig_sequence=no_stop,
                              start=1, end=len(no_stop), strand=1)
    assert "no_valid_stop_codon" in notes
    _, notes = verify_protein("M" + protein[:5] + "*" + protein[6:],
                              contig_sequence=coding, start=1,
                              end=len(coding), strand=1)
    assert "internal_stop_codon" in notes
    verified, notes = verify_protein(protein, contig_sequence=coding,
                                     start=1, end=len(coding) + 500, strand=1)
    assert not verified and "gene_coordinates_out_of_range" in notes
    verified, notes = verify_protein(protein, contig_sequence=coding,
                                     start=1, end=len(coding), strand="?")
    assert not verified and "unknown_strand" in notes


def test_verify_protein_without_genome_coordinates_is_unverified():
    verified, notes = verify_protein("M" + "A" * 20)
    assert not verified
    assert "no_genome_coordinates" in notes


# --- system collection ---------------------------------------------------

def test_collect_systems_reads_pyrodigal_headers_without_partial_field():
    """The regression that motivated this module: no partial=00 in the header."""
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(23)
        protein = _random_protein(rng, 50)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        summary = check_genome(genome_path)
        systems = collect_systems(Path(d) / "annotation", summary,
                                  genome_path=genome_path)
        assert len(systems) == 1
        system = systems[0]
        assert "partial" not in system["protein_header"]
        assert system["integrity_verified"], system["integrity_notes"]
        assert "partial_or_unverified_gene_boundary" not in system["issues"]
        assert system["issues"] == [], system["issues"]
        assert select_system(systems) is system


def test_collect_systems_reads_putative_table_when_good_table_absent():
    """A genome whose only Cas9 sits in an Ambiguous operon has no good table."""
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(29)
        protein = _random_protein(rng, 50)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool,
                        table="cas_operons_putative.tab", prediction="Ambiguous")
        summary = check_genome(genome_path)
        systems = collect_systems(Path(d) / "annotation", summary,
                                  genome_path=genome_path)
        assert len(systems) == 1, "putative-only Cas9 was silently dropped"
        assert systems[0]["operon_table"] == "cas_operons_putative.tab"
        assert "ambiguous_or_putative_subtype" in systems[0]["issues"]
        assert select_system(systems) is None, "ambiguous system was auto-selected"
        # It is still reachable through an explicit, informed user choice.
        chosen = select_system(systems, systems[0]["system_id"], allow_flagged=True)
        assert chosen["system_id"] == systems[0]["system_id"]


def test_collect_systems_rejects_unsupported_subtype_and_missing_protein():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(31)
        protein = _random_protein(rng, 30)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        directory = _annotation_dir(d, cas9_dna=cas9, contig="CT1",
                                    pool_dna=pool, prediction="I-E")
        summary = check_genome(genome_path)
        systems = collect_systems(directory, summary, genome_path=genome_path)
        assert "unsupported_subtype" in systems[0]["issues"]
        # Case-insensitive Cas9 detection still fires for cas9/csn1 HMM names.
        assert systems[0]["gene_name"] == "Cas9_1_II"


def test_collect_systems_no_cas9_returns_empty():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(37)
        cas9 = _dna_for(_random_protein(rng, 20), rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        directory = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        # Replace the gene list with a non-Cas9 HMM name.
        table = (directory / "cas_operons.tab").read_text().replace("'Cas9_1_II'", "'Cas3_1_CAS-I'")
        _write(directory / "cas_operons.tab", table)
        summary = check_genome(genome_path)
        assert collect_systems(directory, summary, genome_path=genome_path) == []


def test_collect_systems_rejects_corrupt_tables():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(41)
        cas9 = _dna_for(_random_protein(rng, 20), rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        directory = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        summary = check_genome(genome_path)

        # Missing required column.
        _write(directory / "cas_operons.tab", "Contig\tOperon\nCT1\tCT1@1\n")
        try:
            collect_systems(directory, summary, genome_path=genome_path)
        except ValueError as exc:
            assert "column" in str(exc).lower()
        else:
            raise AssertionError("wrong columns were accepted")

        _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        _write(directory / "cas_operons.tab",
               "Contig\tOperon\tPrediction\tGenes\tPositions\n"
               "CT1\tCT1@1\tII-A\t['Cas9_1_II']\t[1, 2]\n")
        try:
            collect_systems(directory, summary, genome_path=genome_path)
        except ValueError as exc:
            assert "mismatch" in str(exc).lower()
        else:
            raise AssertionError("gene/position length mismatch was accepted")

        _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        _write(directory / "cas_operons.tab",
               "Contig\tOperon\tPrediction\tGenes\tPositions\n"
               "CT1\tCT1@1\tII-A\tnot_a_list\t[1]\n")
        try:
            collect_systems(directory, summary, genome_path=genome_path)
        except ValueError as exc:
            assert "malformed" in str(exc).lower()
            return
        raise AssertionError("malformed list field was accepted")


def test_collect_systems_requires_proteins_faa():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(43)
        cas9 = _dna_for(_random_protein(rng, 20), rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        directory = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        (directory / "proteins.faa").unlink()
        summary = check_genome(genome_path)
        try:
            collect_systems(directory, summary, genome_path=genome_path)
        except ValueError as exc:
            assert "proteins.faa" in str(exc)
            return
        raise AssertionError("missing proteins.faa was accepted")


# --- arrays --------------------------------------------------------------

def test_array_issues_are_reported_not_hidden():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(47)
        protein = _random_protein(rng, 30)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        directory = _annotation_dir(d, cas9_dna=cas9, contig="CT1",
                                    pool_dna=pool, array_end=10 ** 9)
        (directory / "spacers" / "CT1_1.fa").unlink()
        summary = check_genome(genome_path)
        arrays = collect_systems(directory, summary, genome_path=genome_path)[0]["arrays"]
        assert len(arrays) == 1
        assert "array_coordinates_out_of_range" in arrays[0]["issues"]
        assert "spacer_file_missing" in arrays[0]["issues"]
        assert arrays[0]["spacer_file_present"] is False


def test_multiple_arrays_are_all_listed_and_orientation_unconfirmed():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(53)
        protein = _random_protein(rng, 30)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        directory = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        # Add a second array linked to the same operon.
        crisprs = (directory / "crisprs_all.tab").read_text().rstrip("\n")
        extra = crisprs.splitlines()[-1].replace("CT1_1", "CT1_2").replace("\tTrue\t", "\tFalse\t")
        _write(directory / "crisprs_all.tab", crisprs + "\n" + extra + "\n")
        _write(directory / "spacers" / "CT1_2.fa", ">CT1_2:1\nACGTACGT\n")
        crispr_cas = (directory / "CRISPR_Cas.tab").read_text().replace(
            "['CT1_1']", "['CT1_1', 'CT1_2']")
        _write(directory / "CRISPR_Cas.tab", crispr_cas)

        summary = check_genome(genome_path)
        arrays = collect_systems(directory, summary, genome_path=genome_path)[0]["arrays"]
        assert len(arrays) == 2
        assert {a["array_id"] for a in arrays} == {"CT1_1", "CT1_2"}
        assert all(a["orientation"] == "unconfirmed" for a in arrays)
        assert [a["trusted"] for a in arrays if a["array_id"] == "CT1_2"] == [False]


def test_unsafe_array_identifier_is_rejected():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(59)
        cas9 = _dna_for(_random_protein(rng, 20), rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        directory = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        table = (directory / "crisprs_all.tab").read_text().replace("CT1_1", "../../etc/passwd")
        _write(directory / "crisprs_all.tab", table)
        crispr_cas = (directory / "CRISPR_Cas.tab").read_text().replace("CT1_1", "../../etc/passwd")
        _write(directory / "CRISPR_Cas.tab", crispr_cas)
        summary = check_genome(genome_path)
        try:
            collect_systems(directory, summary, genome_path=genome_path)
        except ValueError as exc:
            assert "unsafe" in str(exc).lower()
            return
        raise AssertionError("path-escaping array identifier was accepted")


def test_safe_id_and_literal_list_helpers():
    for bad in ["../x", "a/b", "", ".", "..", "a\\b"]:
        try:
            safe_id(bad)
        except ValueError:
            continue
        raise AssertionError(f"unsafe id accepted: {bad!r}")
    assert safe_id("NC_002163.1_1") == "NC_002163.1_1"
    assert literal_list("['a', 'b']") == ["a", "b"]
    assert literal_list("[1, 2]") == [1, 2]
    try:
        literal_list("not_a_list")
    except ValueError:
        pass
    else:
        raise AssertionError("malformed literal accepted")


def test_select_system_auto_only_when_single_eligible():
    one = {"system_id": "A", "issues": []}
    two = {"system_id": "B", "issues": []}
    bad = {"system_id": "C", "issues": ["ambiguous_or_putative_subtype"]}
    assert select_system([one])["system_id"] == "A"
    assert select_system([one, bad])["system_id"] == "A"
    assert select_system([one, two]) is None
    assert select_system([bad]) is None
    assert select_system([one, two], "B")["system_id"] == "B"
    for args in (([one], "NOPE"), ([bad], "C"), ([one, two], "C")):
        try:
            select_system(*args)
        except ValueError:
            continue
        raise AssertionError(f"select_system accepted {args}")


def test_explicit_selection_overrides_annotation_flags_but_never_integrity():
    """A user may accept an ambiguous subtype call, but never a broken protein."""
    flagged = {"system_id": "C", "issues": ["ambiguous_or_putative_subtype"]}
    broken = {"system_id": "D", "issues": ["gene_integrity_unverified"]}
    # Overridable annotation-confidence flag.
    chosen = select_system([flagged], "C", allow_flagged=True)
    assert chosen["system_id"] == "C"
    # Protein that failed genome-level verification stays refused.
    try:
        select_system([broken], "D", allow_flagged=True)
    except ValueError as exc:
        assert "integrity" in str(exc)
    else:
        raise AssertionError("integrity failure was overridable")
    # A system with both a flag and an integrity failure is also refused.
    both = {"system_id": "E", "issues": ["ambiguous_or_putative_subtype",
                                         "gene_integrity_unverified"]}
    try:
        select_system([both], "E", allow_flagged=True)
    except ValueError:
        pass
    else:
        raise AssertionError("mixed integrity failure was overridable")


def test_single_contig_with_two_systems_keeps_arrays_separate():
    """Spacers from different systems on one contig must not be merged."""
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(61)
        protein = _random_protein(rng, 30)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        directory = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        # Two operons, each linked to its own array.
        header = "Contig\tOperon\tPrediction\tGenes\tPositions\n"
        _write(directory / "cas_operons_putative.tab",
               header
               + "CT1\tCT1@9\tII-C\t['Cas9_1_II']\t[1]\n")
        crisprs = (directory / "crisprs_all.tab").read_text().rstrip("\n")
        extra = crisprs.splitlines()[-1].replace("CT1_1", "CT1_9")
        _write(directory / "crisprs_all.tab", crisprs + "\n" + extra + "\n")
        _write(directory / "spacers" / "CT1_9.fa", ">CT1_9:1\nACGTACGT\n")
        _write(directory / "CRISPR_Cas.tab",
               "Contig\tOperon\tPrediction\tCRISPRs\n"
               "CT1\tCT1@1\tII-A\t['CT1_1']\n"
               "CT1\tCT1@9\tII-C\t['CT1_9']\n")
        # Now the good operon must not pick up the other operon's array.
        _write(directory / "cas_operons.tab",
               "Contig\tOperon\tPrediction\tGenes\tPositions\n"
               "CT1\tCT1@1\tII-A\t['Cas9_1_II']\t[1]\n")
        summary = check_genome(genome_path)
        systems = {s["operon"]: s for s in
                   collect_systems(directory, summary, genome_path=genome_path)}
        assert set(systems) == {"CT1@1", "CT1@9"}
        assert [a["array_id"] for a in systems["CT1@1"]["arrays"]] == ["CT1_1"]
        assert [a["array_id"] for a in systems["CT1@9"]["arrays"]] == ["CT1_9"]


# --- CLI / state machine -------------------------------------------------

def _agent(args, timeout=120, stub_predictor=None):
    """Run the agent as a child process.

    ``stub_predictor`` swaps in a fake ``run_pam_report.py`` so tests can drive
    the whole state machine (completion, resume, integrity) without the real
    model weights. Normal use leaves it ``None`` and runs the real inference.
    """
    env = dict(os.environ)
    if stub_predictor is not None:
        env["PAMDICT_PREDICTOR_SCRIPT"] = str(stub_predictor)
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts/run_genome_agent.py"), *args],
        capture_output=True, text=True, timeout=timeout, cwd=str(ROOT), env=env,
    )


def _stub_predictor(directory, *, consensus="NGG", positions=(2, 3),
                    fail=False):
    """Write a stand-in predictor that emits a realistic protein.json.

    Keeping the offline suite free of model weights matters: without this the
    tests silently required a multi-GB checkpoint and failed on a clean machine.
    """
    directory = Path(directory)
    matrix = [[0.25] * 4 for _ in range(10)]
    matrix[1] = [0.02, 0.01, 0.95, 0.02]
    matrix[2] = [0.03, 0.02, 0.93, 0.02]
    payload = {
        "model": "stub_not_protein2pam",
        "pam_side": "downstream",
        "proteins": [{
            "protein_id": "p1",
            "predicted_pam": consensus,
            "no_clear_pam": not consensus.replace("N", ""),
            "partial_consensus": len(positions) < len(consensus),
            "determined_positions": list(positions),
            "probability_matrix": matrix,
            "information": {"signal_positions": list(positions)},
        }],
        "scores": [{
            "candidate_pam": consensus, "pam_length": len(consensus),
            "rank_within_length": 1, "specificity_adjusted_score": 1.0,
            "allowed_probability_geomean": 0.9,
        }] if consensus.replace("N", "") else [],
    }
    script = directory / "stub_predictor.py"
    script.write_text(
        "import json, sys, pathlib\n"
        "args = sys.argv[1:]\n"
        "out = pathlib.Path(args[args.index('--outdir') + 1])\n"
        "out.mkdir(parents=True, exist_ok=True)\n"
        f"payload = {payload!r}\n"
        + ("raise SystemExit(3)\n" if fail else "")
        + "(out / 'report.json').write_text(json.dumps(payload))\n"
        "(out / 'report.html').write_text('<html>stub</html>')\n"
        "(out / 'candidates.tsv').write_text('candidate_pam\\n')\n"
    )
    script.chmod(0o755)
    return script


def test_agent_rejects_bad_genome_before_writing_output():
    with tempfile.TemporaryDirectory() as d:
        bad = Path(d) / "bad.fna"
        bad.write_text(">c1\nACGTZZ\n")
        out = Path(d) / "out"
        result = _agent(["--genome", str(bad), "--outdir", str(out)])
        assert result.returncode != 0
        assert not out.exists(), "output directory was created for an invalid genome"


def test_agent_reports_missing_annotator_as_blocked_not_success():
    with tempfile.TemporaryDirectory() as d:
        genome = _make_genome_file(d)
        out = Path(d) / "out"
        result = _agent(["--genome", str(genome), "--outdir", str(out),
                         "--annotator", "/nonexistent/cctyper"])
        assert result.returncode == 3, result.stderr
        job = json.loads((out / "job.json").read_text())
        assert job["status"] == "blocked_dependency"
        assert "No prediction was made" in job["message"]
        # No success artefacts.
        assert not (out / "pam" / "report.html").exists()
        assert not (out / "report.html").exists()
        page = (out / "index.html").read_text()
        assert "依赖缺失" in page


def test_agent_refuses_to_overwrite_existing_outdir():
    with tempfile.TemporaryDirectory() as d:
        genome = _make_genome_file(d)
        out = Path(d) / "out"
        out.mkdir()
        (out / "keep.txt").write_text("existing")
        result = _agent(["--genome", str(genome), "--outdir", str(out),
                         "--annotator", "/nonexistent/cctyper"])
        assert result.returncode != 0
        assert (out / "keep.txt").read_text() == "existing"


def test_agent_rejects_candidates_with_auto_candidates():
    with tempfile.TemporaryDirectory() as d:
        genome = _make_genome_file(d)
        out = Path(d) / "out"
        result = _agent(["--genome", str(genome), "--outdir", str(out),
                         "--candidates", "NGG", "--auto-candidates"])
        assert result.returncode != 0
        assert "mutually exclusive" in (result.stderr + result.stdout)


def test_agent_rejects_nonpositive_settings():
    with tempfile.TemporaryDirectory() as d:
        genome = _make_genome_file(d)
        for flag, value in (("--threads", "0"), ("--annotation-timeout", "0"),
                            ("--inference-timeout", "-1")):
            out = Path(d) / f"out{flag}"
            result = _agent(["--genome", str(genome), "--outdir", str(out),
                             flag, value])
            assert result.returncode != 0, flag


def test_agent_timeout_kills_process_group_and_leaves_no_success():
    """A stage that spawns a grandchild must not leak it on timeout."""
    with tempfile.TemporaryDirectory() as d:
        genome = _make_genome_file(d)
        out = Path(d) / "out"
        # A fake annotator that starts a long-lived grandchild, then hangs.
        marker = Path(d) / "grandchild.pid"
        fake = Path(d) / "fake_cctyper"
        fake.write_text(
            "#!/usr/bin/env python3\n"
            "import subprocess, sys, time\n"
            "p = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(600)'])\n"
            f"open({str(marker)!r}, 'w').write(str(p.pid))\n"
            "time.sleep(600)\n"
        )
        fake.chmod(0o755)
        result = _agent(["--genome", str(genome), "--outdir", str(out),
                         "--annotator", str(fake), "--annotation-timeout", "3"],
                        timeout=90)
        assert result.returncode != 0
        job = json.loads((out / "job.json").read_text())
        assert job["status"] == "failed"
        assert "Timeout" in job["message"] or "terminated" in job["error"]
        assert not (out / "pam" / "report.html").exists()
        # The grandchild must be gone: kill(pid, 0) must fail.
        import os as _os
        import time as _time
        assert marker.exists(), "fake annotator never started"
        pid = int(marker.read_text())
        for _ in range(40):
            try:
                _os.kill(pid, 0)
            except ProcessLookupError:
                return
            _time.sleep(0.25)
        raise AssertionError(f"grandchild process {pid} survived the timeout")


def test_agent_flags_annotation_that_exits_zero_without_output():
    """Exit 0 with no usable tables is not a reliable negative result."""
    with tempfile.TemporaryDirectory() as d:
        genome = _make_genome_file(d)
        out = Path(d) / "out"
        fake = Path(d) / "fake_cctyper"
        fake.write_text("#!/usr/bin/env python3\nimport sys, pathlib\n"
                        "pathlib.Path(sys.argv[2]).mkdir(parents=True, exist_ok=True)\n"
                        "sys.exit(0)\n")
        fake.chmod(0o755)
        result = _agent(["--genome", str(genome), "--outdir", str(out),
                         "--annotator", str(fake)])
        assert result.returncode == 1, result.stderr
        job = json.loads((out / "job.json").read_text())
        assert job["status"] == "failed"
        assert "not a reliable negative" in job["error"]
        assert not (out / "pam" / "report.html").exists()


def test_agent_reports_no_cas9_without_claiming_no_crispr():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(67)
        cas9 = _dna_for(_random_protein(rng, 20), rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        annotation = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        table = (annotation / "cas_operons.tab").read_text().replace(
            "'Cas9_1_II'", "'Cas3_1_CAS-I'")
        _write(annotation / "cas_operons.tab", table)

        out = Path(d) / "out"
        fake = Path(d) / "fake_cctyper"
        fake.write_text(
            "#!/usr/bin/env python3\nimport shutil, sys\n"
            f"shutil.copytree({str(annotation)!r}, sys.argv[2], dirs_exist_ok=True)\n"
        )
        fake.chmod(0o755)
        stub = _stub_predictor(Path(d))
        result = _agent(["--genome", str(genome_path), "--outdir", str(out),
                         "--annotator", str(fake)],
                        stub_predictor=stub)
        assert result.returncode == 0, result.stderr
        job = json.loads((out / "job.json").read_text())
        assert job["status"] == "no_supported_cas9"
        assert "does not prove" in job["message"]
        assert not (out / "pam" / "report.html").exists()


def test_agent_requires_array_confirmation_for_spacer_evidence():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(71)
        protein = _random_protein(rng, 30)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        annotation = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        target = Path(d) / "phage.fna"
        target.write_text(">phage1\n" + "ACGT" * 40 + "\n")

        out = Path(d) / "out"
        fake = Path(d) / "fake_cctyper"
        fake.write_text(
            "#!/usr/bin/env python3\nimport shutil, sys\n"
            f"shutil.copytree({str(annotation)!r}, sys.argv[2], dirs_exist_ok=True)\n"
        )
        fake.chmod(0o755)
        result = _agent(["--genome", str(genome_path), "--outdir", str(out),
                         "--annotator", str(fake), "--targets", str(target)])
        assert result.returncode == 2, result.stderr
        job = json.loads((out / "job.json").read_text())
        assert job["status"] == "needs_array_confirmation"
        assert "NOT inferred" in job["message"]
        assert not (out / "pam" / "report.html").exists()


def test_agent_resume_rejects_changed_genome_and_changed_annotation():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(73)
        protein = _random_protein(rng, 30)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        annotation = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)

        out = Path(d) / "out"
        fake = Path(d) / "fake_cctyper"
        fake.write_text(
            "#!/usr/bin/env python3\nimport shutil, sys\n"
            f"shutil.copytree({str(annotation)!r}, sys.argv[2], dirs_exist_ok=True)\n"
        )
        fake.chmod(0o755)
        stub = _stub_predictor(Path(d))
        first = _agent(["--genome", str(genome_path), "--outdir", str(out),
                        "--annotator", str(fake)],
                       stub_predictor=stub)
        assert first.returncode == 0, first.stderr
        job = json.loads((out / "job.json").read_text())
        assert job["status"] == "complete"
        assert (out / "pam" / "report.html").exists()

        # Changing the genome must be rejected.
        genome_path.write_text(">CT1\n" + "ACGT" * 20 + "\n")
        changed = _agent(["--genome", str(genome_path), "--outdir", str(out), "--resume"])
        assert changed.returncode != 0
        assert "genome file changed" in (changed.stderr + changed.stdout)

        # Restoring the genome but mutating the annotation must also be rejected.
        genome_path.write_text(f">CT1 test contig\n{pool}\n")
        (out / "annotation" / "genes.tab").write_text("Contig\tStart\tEnd\tStrand\tPos\n")
        mutated = _agent(["--genome", str(genome_path), "--outdir", str(out), "--resume"])
        assert mutated.returncode != 0
        assert "annotation output changed" in (mutated.stderr + mutated.stdout)

        # A complete job cannot be resumed into a new input.
        complete = _agent(["--genome", str(genome_path), "--outdir", str(out), "--resume"])
        assert complete.returncode != 0


def test_rejected_resume_does_not_destroy_completed_job_record():
    """A refused --resume must leave the existing successful job untouched."""
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(83)
        protein = _random_protein(rng, 30)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        annotation = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        out = Path(d) / "out"
        fake = Path(d) / "fake_cctyper"
        fake.write_text(
            "#!/usr/bin/env python3\nimport shutil, sys\n"
            f"shutil.copytree({str(annotation)!r}, sys.argv[2], dirs_exist_ok=True)\n"
        )
        fake.chmod(0o755)
        stub = _stub_predictor(Path(d))
        assert _agent(["--genome", str(genome_path), "--outdir", str(out),
                       "--annotator", str(fake)],
                      stub_predictor=stub).returncode == 0
        job_path = out / "job.json"
        before = job_path.read_text()
        assert json.loads(before)["status"] == "complete"

        # Reject on a changed genome and confirm the on-disk record is intact.
        genome_path.write_text(">CT1\n" + "ACGT" * 20 + "\n")
        result = _agent(["--genome", str(genome_path), "--outdir", str(out), "--resume"])
        assert result.returncode != 0
        after = json.loads(job_path.read_text())
        assert after["status"] == "complete", "rejected resume clobbered the job status"
        assert after["report"] == "pam/report.html"
        # The report itself must still be present and readable.
        assert (out / "pam" / "report.html").is_file()


def test_agent_resume_refuses_when_partial_prediction_exists():
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(79)
        protein = _random_protein(rng, 30)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        annotation = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool)
        out = Path(d) / "out"
        fake = Path(d) / "fake_cctyper"
        fake.write_text(
            "#!/usr/bin/env python3\nimport shutil, sys\n"
            f"shutil.copytree({str(annotation)!r}, sys.argv[2], dirs_exist_ok=True)\n"
        )
        fake.chmod(0o755)
        stub = _stub_predictor(Path(d))
        assert _agent(["--genome", str(genome_path), "--outdir", str(out),
                       "--annotator", str(fake)],
                      stub_predictor=stub).returncode == 0
        # Simulate a run that crashed part-way through prediction: the status
        # never advanced past 'predicting' and pam/ holds only partial output.
        job_path = out / "job.json"
        job = json.loads(job_path.read_text())
        job["status"] = "predicting"
        job_path.write_text(json.dumps(job))
        shutil.rmtree(out / "pam")
        (out / "pam").mkdir()
        (out / "pam" / "half.tsv").write_text("partial")
        result = _agent(["--genome", str(genome_path), "--outdir", str(out), "--resume"])
        assert result.returncode != 0
        assert "partial prediction" in (result.stderr + result.stdout)
        # The partial output must be preserved, not deleted.
        assert (out / "pam" / "half.tsv").read_text() == "partial"


def test_report_marks_no_clear_pam_without_fake_candidates():
    """An all-N / empty consensus must not be presented as a PAM result."""
    from pamdict.delivery import make_report, write_report
    protein = {
        "model": "cas9_full",
        "candidate_scope_note": "consensus only",
        "no_clear_pam": True,
        "proteins": [{"protein_id": "p1", "predicted_pam": None}],
        "scores": [],
    }
    with tempfile.TemporaryDirectory() as d:
        report = make_report(protein, mode="real_inference", provenance={})
        assert report["no_clear_pam"] is True
        assert report["candidates"] == []
        write_report(report, d)
        page = (Path(d) / "report.html").read_text()
        assert "未得到明确 PAM" in page
        assert "无明确共识" in page
        tsv = (Path(d) / "candidates.tsv").read_text()
        assert tsv.startswith("candidate_pam\t")
        assert "\n" in tsv and len(tsv.strip().splitlines()) == 1


def test_report_renders_rows_without_relying_on_dict_order():
    from pamdict.delivery import make_report, write_report
    protein = {
        "model": "cas9_full",
        "proteins": [{"protein_id": "p1", "predicted_pam": "NGG"}],
        "scores": [{"candidate_pam": "NGG", "specificity_adjusted_score": 91.0,
                    "rank_within_length": 1}],
    }
    with tempfile.TemporaryDirectory() as d:
        report = make_report(protein, mode="real_inference", provenance={})
        write_report(report, d)
        lines = (Path(d) / "candidates.tsv").read_text().strip().splitlines()
        assert len(lines) == 2
        assert "NGG" in lines[1]


def test_partial_consensus_is_not_reported_as_a_strict_pam_call():
    """A consensus with undetermined N positions represents a set, not one PAM."""
    from pamdict.delivery import make_report, write_report
    protein = {
        "model": "cas9_full",
        "proteins": [{"protein_id": "p1", "predicted_pam": "NNNNANA",
                      "partial_consensus": True, "determined_positions": [2, 4, 6]}],
        "scores": [],
    }
    with tempfile.TemporaryDirectory() as d:
        report = make_report(protein, mode="real_inference", provenance={})
        assert report["partial_consensus"] == ["NNNNANA"]
        write_report(report, d)
        page = (Path(d) / "report.html").read_text()
        assert "一组 PAM" in page
        assert "未判定" in page
        # The determined positions are surfaced so N is never read as a base.
        assert "[2, 4, 6]" in page


def test_agent_grades_partial_consensus_distinctly():
    """_summarise_pam_result must not label a partly-unknown consensus as clean."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "run_genome_agent", str(ROOT / "scripts/run_genome_agent.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    with tempfile.TemporaryDirectory() as d:
        pam = Path(d)
        # determined_positions is published by score_candidate_pams using the
        # consensus probability rule, and consumed here -- never re-derived from
        # the separate information diagnostic.
        cases = [
            ({"predicted_pam": "NNNNNNNNNN", "determined_positions": []},
             "no_clear_pam"),
            ({"predicted_pam": "", "determined_positions": []},
             "no_clear_pam"),
            # Real observed output for a divergent Cas9: only positions 5 and 7.
            ({"predicted_pam": "NNNNANA", "determined_positions": [5, 7]},
             "partial_consensus"),
            # Real SpCas9-like output: position 1 is not called, so this denotes
            # a set of PAMs even though it is the expected consensus.
            ({"predicted_pam": "NGG", "determined_positions": [2, 3]},
             "partial_consensus"),
            ({"predicted_pam": "AAA", "determined_positions": [1, 2, 3]},
             "model_consensus_only"),
        ]
        for entry, expected in cases:
            (pam / "report.json").write_text(json.dumps({"proteins": [entry]}))
            job = {}
            mod._summarise_pam_result(pam, job)
            assert job["pam_conclusion"] == expected, (
                entry["predicted_pam"], job["pam_conclusion"])
            assert job["pam_conclusion_note"]
        # The published positions are passed through unchanged.
        (pam / "report.json").write_text(json.dumps(
            {"proteins": [{"predicted_pam": "NGG", "determined_positions": [2, 3]}]}))
        job = {}
        mod._summarise_pam_result(pam, job)
        assert job["determined_positions"] == [[2, 3]], job["determined_positions"]
        # A fully called consensus carries no warning.
        (pam / "report.json").write_text(json.dumps(
            {"proteins": [{"predicted_pam": "AAA", "determined_positions": [1, 2, 3]}]}))
        job = {}
        mod._summarise_pam_result(pam, job)
        assert not job.get("warnings")
        # An artifact from the older scoring code is refused rather than graded
        # on a value the agent would have to recompute incorrectly.
        (pam / "report.json").write_text(json.dumps(
            {"proteins": [{"predicted_pam": "NGG"}]}))
        try:
            mod._summarise_pam_result(pam, {})
        except RuntimeError as exc:
            assert "determined_positions" in str(exc)
        else:
            raise AssertionError("stale protein.json was graded silently")


def test_summarise_pam_result_requires_report_json():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "run_genome_agent2", str(ROOT / "scripts/run_genome_agent.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    with tempfile.TemporaryDirectory() as d:
        try:
            mod._summarise_pam_result(Path(d), {})
        except RuntimeError as exc:
            assert "report.json" in str(exc)
            return
        raise AssertionError("missing report.json was accepted")


# --- regressions for the three reported acceptance failures ----------------

def _stub_protein_json(path, *, consensus, matrix, scores=None, determined=None,
                       partial=False, no_clear=False):
    import numpy as np
    _write(Path(path), json.dumps({
        "model": "cas9_full", "pam_side": "downstream",
        "proteins": [{
            "protein_id": "p1", "predicted_pam": consensus,
            "no_clear_pam": no_clear, "partial_consensus": partial,
            "determined_positions": determined or [],
            "probability_matrix": matrix,
            "information": {"signal_positions": []},
        }],
        "scores": scores or [],
    }))


def test_regression_ngg_auto_mode_still_yields_a_candidate():
    """A canonical NGG consensus must not be dropped just because it contains N.

    Dropping it emptied the candidate list, so the spacer stage never ran and the
    dual-evidence path was unreachable. The gate is exercised directly on the
    real observed NGG matrix; no model weights are needed.
    """
    from pamdict.infer.p2pam import CONSENSUS_THRESHOLD_BITS, consensus_from_info
    from pamdict.score.candidate import rank_candidate_pams
    import numpy as np

    # Real observed cas9_full matrix for NZ_CP007240.1_685: only positions 2 and
    # 3 reach the 0.70 consensus rule, so the consensus is NGG.
    matrix = np.array([
        [0.246, 0.238, 0.165, 0.351],
        [0.028, 0.007, 0.953, 0.012],
        [0.036, 0.025, 0.925, 0.014],
        [0.211, 0.330, 0.122, 0.337],
    ] + [[0.25] * 4 for _ in range(6)])
    # consensus_from_info expects an information matrix; build it the same way
    # score_candidate_pams does and confirm the consensus really is NGG.
    from pamdict.infer.p2pam import prob_to_info_numpy
    consensus = consensus_from_info(prob_to_info_numpy(matrix), side="downstream")
    assert consensus == "NGG", consensus

    determined = [i + 1 for i, row in enumerate(matrix[:len(consensus)])
                  if float(max(row)) >= CONSENSUS_THRESHOLD_BITS]
    assert determined == [2, 3], determined
    # It is 'partial' because position 1 was not called...
    assert len(determined) < len(consensus)
    # ...but it is NOT unresolved, so it must still be scored.
    unresolved = not consensus.replace("N", "")
    assert not unresolved
    scored = rank_candidate_pams(matrix, [consensus], side="downstream")
    assert scored, "NGG consensus produced no candidate (spacer stage would break)"
    assert scored[0].candidate_pam == "NGG"
    assert scored[0].pam_length == 3


def test_regression_determined_positions_use_consensus_threshold():
    """determined_positions must follow the 0.70 consensus rule, not 0.15 bits.

    Real case: consensus NNNNANA has concrete base calls only at positions 5 and
    7, but the information diagnostic (0.15-bit floor) marked 2-7 as signal and
    the report wrongly claimed all of them were called.
    """
    from pamdict.infer.p2pam import CONSENSUS_THRESHOLD_BITS
    assert CONSENSUS_THRESHOLD_BITS == 0.70
    matrix = [
        [0.400, 0.137, 0.197, 0.266],
        [0.478, 0.130, 0.204, 0.188],
        [0.684, 0.043, 0.150, 0.123],
        [0.656, 0.063, 0.150, 0.130],
        [0.915, 0.041, 0.029, 0.014],
        [0.071, 0.423, 0.006, 0.500],
        [0.832, 0.024, 0.074, 0.070],
    ] + [[0.25] * 4 for _ in range(3)]
    called = [
        i + 1 for i, row in enumerate(matrix[:7])
        if max(row) >= CONSENSUS_THRESHOLD_BITS
    ]
    assert called == [5, 7], called
    # Positions 2,3,4 are genuinely below 0.70 and must not be reported as calls.
    assert 2 not in called and 3 not in called and 4 not in called


def test_regression_orphan_cas9_without_arrays_is_not_a_failure():
    """A Cas operon with no CRISPR array is a normal orphan result.

    CCTyper logs "No CRISPRs found." and deletes crisprs_all.tab; the run must
    still reach protein-only inference instead of being marked failed.
    """
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(103)
        protein = _random_protein(rng, 40)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        annotation = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool,
                                     table="cas_operons_putative.tab",
                                     prediction="Ambiguous")
        # Reproduce the real orphan layout: no array tables at all.
        (annotation / "crisprs_all.tab").unlink()
        (annotation / "CRISPR_Cas.tab").unlink()
        for f in (annotation / "spacers").glob("*"):
            f.unlink()

        # The emptiness check must accept this as a valid annotation.
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "rga_orphan", str(ROOT / "scripts/run_genome_agent.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        assert mod.annotation_is_scientifically_empty(annotation) is None, \
            "orphan Cas9 annotation was treated as empty/failed"

        systems = collect_systems(annotation, check_genome(genome_path),
                                  genome_path=genome_path)
        assert len(systems) == 1
        assert systems[0]["arrays"] == []
        assert "gene_integrity_unverified" not in systems[0]["issues"]


def test_regression_orphan_table_is_read_when_present():
    """cas_operons_orphan.tab must be parsed, not ignored."""
    with tempfile.TemporaryDirectory() as d:
        rng = random.Random(107)
        protein = _random_protein(rng, 40)
        cas9 = _dna_for(protein, rng)
        pool = "ACGT" * 5 + cas9 + "ACGT" * 5
        genome_path = _make_genome_file(d, contig="CT1", dna=pool)
        annotation = _annotation_dir(d, cas9_dna=cas9, contig="CT1", pool_dna=pool,
                                     table="cas_operons_orphan.tab",
                                     prediction="Ambiguous")
        summary = check_genome(genome_path)
        systems = collect_systems(annotation, summary, genome_path=genome_path)
        assert len(systems) == 1, "cas_operons_orphan.tab was ignored"
        assert systems[0]["operon_table"] == "cas_operons_orphan.tab"


def test_regression_annotation_still_requires_gene_and_protein_outputs():
    """Relaxing the array requirement must not weaken the real emptiness check."""
    with tempfile.TemporaryDirectory() as d:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "rga_empty", str(ROOT / "scripts/run_genome_agent.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        empty_dir = Path(d)
        assert mod.annotation_is_scientifically_empty(empty_dir) is not None
        _write(empty_dir / "genes.tab", "Contig\tStart\tEnd\tStrand\tPos\n")
        assert mod.annotation_is_scientifically_empty(empty_dir) is not None
        _write(empty_dir / "proteins.faa", ">a\nMAA\n")
        # genes + proteins but no operon table at all is still not a Cas call.
        assert mod.annotation_is_scientifically_empty(empty_dir) is not None
        _write(empty_dir / "cas_operons_orphan.tab",
               "Contig\tOperon\tPrediction\tGenes\tPositions\n")
        assert mod.annotation_is_scientifically_empty(empty_dir) is None


# --- end-of-round fixes: portable tests + the no-Cas negative --------------

def test_regression_upstream_crash_is_unknown_not_negative():
    """An upstream exception cannot establish a biological negative."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "rga_nocas", str(ROOT / "scripts/run_genome_agent.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    real_trace = (
        "Traceback (most recent call last):\n"
        "KeyError: \"['Hmm', 'ORF'] not in index\"\n"
    )
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        for f in ("genes.tab", "proteins.faa"):
            _write(d / f, "x")
        (d / "hmmer").mkdir()          # CCTyper creates this before crashing
        _write(d / "annotating.log", real_trace)
        reason = mod.annotation_upstream_error(d, d / "annotating.log")
        assert reason is not None
        assert "无法判断" in reason and "不是阴性结果" in reason


def test_regression_upstream_error_detection_is_narrow():
    """Only the matching crash receives this diagnostic; it remains failure."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "rga_narrow", str(ROOT / "scripts/run_genome_agent.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    real_trace = "KeyError: \"['Hmm', 'ORF'] not in index\"\n"
    cases = [
        # (files, log, should receive the known-crash diagnostic)
        ((), real_trace, True),
        (("hmmer.tab",), real_trace, False),
        (("cas_operons.tab",), real_trace, False),
        (("cas_operons_putative.tab",), real_trace, False),
        (("cas_operons_orphan.tab",), real_trace, False),
        (("crisprs_all.tab",), real_trace, False),
        (("CRISPR_Cas.tab",), real_trace, False),
        (("hmmer/hit.txt",), real_trace, False),
        ((), "KeyError: 'some other column'\n", False),
        ((), "no error here at all\n", False),
        ((), "", False),
    ]
    with tempfile.TemporaryDirectory() as d:
        for i, (files, log, expected) in enumerate(cases):
            case = Path(d) / str(i)
            case.mkdir()
            for f in files:
                path = case / f
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("x")
            _write(case / "annotating.log", log)
            got = mod.annotation_upstream_error(
                case, case / "annotating.log") is not None
            assert got is expected, (files, log[:40], got, expected)


def test_upstream_crash_workflow_preserves_unknown_and_nonzero_exit():
    from unittest.mock import patch
    from scripts import run_genome_agent as mod
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        genome, out = d / 'genome.fna', d / 'job'
        genome.write_text('>contig\nACGTACGT\n')
        def crash(command, output, job, stage, timeout, env=None):
            (output / 'annotating.log').write_text('KeyError: "[\'Hmm\', \'ORF\'] not in index"\n')
            raise RuntimeError('Stage annotating exited with code 1')
        with patch.object(sys, 'argv', ['run_genome_agent.py', '--genome', str(genome),
                '--outdir', str(out), '--annotator', sys.executable]), patch.object(mod, 'run_stage', side_effect=crash):
            assert mod.main() == 1
        job = json.loads((out / 'job.json').read_text())
        assert job['status'] == 'failed'
        assert job['annotation_outcome'] == 'unknown'
        assert job['annotated_system_count'] is None
        assert '不是阴性结果' in job['message']
        assert (out / 'annotating.log').is_file()
        assert not (out / 'pam').exists()


def test_offline_suite_does_not_need_the_model():
    """Guard against environment-dependent tests creeping back into this suite.

    The offline suite must run on a clean machine. Any test here that needs the
    protein2pam weights or a historical run artifact belongs in
    ``test_genome_integration`` instead, so this asserts that this module does
    not reference either.
    """
    source = Path(__file__).read_text(encoding="utf-8")
    # Split so this check does not match its own literals.
    banned_markers = [
        "data/parsed/" + "genome_agent_real",
        "HF" + "_HOME",
        "checkpoints/" + "hf",
    ]
    for banned in banned_markers:
        for line in source.splitlines():
            if banned in line:
                assert line.lstrip().startswith("#"), (
                    f"offline test suite references environment-dependent "
                    f"resource {banned!r} outside a comment: {line.strip()}"
                )
