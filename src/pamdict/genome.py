"""CRISPRCasTyper adapter for a bounded, local genome-to-Cas9 workflow.

The parsers here were written against the *actually installed* CRISPRCasTyper
1.9.0 output, not against a floating master README.  The details that matter:

* ``proteins.faa`` headers are produced by ``pyrodigal_gv`` and look like
  ``>NC_002163.1_1464 # 1456880 # 1459834 # -1``.  They contain **no**
  ``partial=`` field, unlike the older Prodigal-based releases.  Gene
  integrity is therefore established from genome coordinates, strand, the
  translation table, the stop codon and translation agreement with the
  reported protein -- never by searching for a header string.
* ``cas_operons.tab`` only holds operons CCTyper considers good.  Genomes whose
  only Cas9 sits in an ``Ambiguous``/``False`` operon produce
  ``cas_operons_putative.tab`` **and no** ``cas_operons.tab``.  Reading only the
  good table silently discards real Cas9 genes, so both are read.
* ``genes.tab`` columns are ``Contig, Start, End, Strand, Pos``.

Nothing here invents a repeat-finding algorithm: all annotation is CCTyper's.
"""
from __future__ import annotations

import ast
import csv
import hashlib
from pathlib import Path

# CCTyper 1.9.0 writes NCBI-style contig names and gene IDs derived from them.
# These are used to build file names, so they are restricted to a safe charset.
_SAFE_ID = set(
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-"
)

# Codons that pyrodigal-gv may use to start a gene.  When one of these is used,
# the translated first residue is methionine regardless of the naive table.
START_CODONS = frozenset({"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"})

STANDARD_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

# Defects that make a called gene boundary untrustworthy.  A translation that
# merely "lines up" is not enough while one of these is present.
BLOCKING_DEFECTS = (
    "internal_stop_codon",
    "nonstandard_residues",
    "gene_length_not_multiple_of_three",
    "no_valid_start_codon",
    "no_valid_stop_codon",
)

# System-level issues that a user's explicit --system-id may override, because
# they describe CCTyper's confidence rather than a defect in the protein.
OVERRIDABLE_ISSUES = frozenset({
    "ambiguous_or_putative_subtype",
    "unsupported_subtype",
})

# System-level issues that must never be overridden: these mean the protein
# could not be reproduced from the genome, so it must not be scored.
BLOCKING_ISSUES = frozenset({
    "gene_integrity_unverified",
    "missing_protein_sequence",
    "gene_coordinates_missing",
    "contig_not_in_input_genome",
    "gene_coordinates_invalid",
})

_COMPLEMENT = str.maketrans("ACGTNRYKMSWBDHV", "TGCANYRMKSWVHDB")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta(path, alphabet, *, allow_duplicates=False):
    """Read a FASTA file, validating identifiers and sequence symbols.

    Returns ``(records, headers)`` where both are dicts keyed by the first
    whitespace-delimited token of each header.
    """
    records: dict[str, str] = {}
    headers: dict[str, str] = {}
    name = None
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            fields = line[1:].split()
            if not fields:
                raise ValueError("FASTA contains an empty identifier")
            if fields[0] in records and not allow_duplicates:
                raise ValueError(f"FASTA contains duplicate identifier: {fields[0]}")
            name = fields[0]
            records.setdefault(name, "")
            headers.setdefault(name, line[1:])
        else:
            if name is None:
                raise ValueError("FASTA sequence before the first header")
            seq = "".join(line.split()).upper()
            invalid = set(seq) - set(alphabet)
            if invalid:
                raise ValueError(
                    f"Invalid sequence symbol(s) {''.join(sorted(invalid))} in {name}"
                )
            records[name] += seq
    if not records:
        raise ValueError("Empty FASTA: no records found")
    empty = [k for k, v in records.items() if not v]
    if empty:
        raise ValueError(f"Empty FASTA record: {empty[0]}")
    return records, headers


def check_genome(path, max_bytes=30_000_000, min_bases=1):
    """Validate an assembled nucleotide FASTA and summarise it."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"Genome FASTA not found: {path}")
    size = path.stat().st_size
    if size == 0:
        raise ValueError("Genome FASTA is empty")
    if size > max_bytes:
        raise ValueError(
            f"First release accepts assembled genomes up to {max_bytes // 1_000_000} MB only"
        )
    records, _ = fasta(path, "ACGTRYSWKMBDHVN")
    bases = sum(len(v) for v in records.values())
    if bases < min_bases:
        raise ValueError("Genome FASTA contains no usable sequence")
    return {
        "sha256": sha256(path),
        "bytes": size,
        "contigs": len(records),
        "bases": bases,
        "contig_lengths": {k: len(v) for k, v in records.items()},
    }


def load_sequences(path, alphabet="ACGTRYSWKMBDHVN"):
    """Load sequences for coordinate-level validation (uppercased)."""
    records, headers = fasta(path, alphabet)
    return records, headers


def table(path, required):
    """Read a CCTyper TSV, tolerating a missing file (returns [])."""
    path = Path(path)
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not set(required) <= fields:
            missing = sorted(set(required) - fields)
            raise ValueError(f"Unsupported CCTyper columns in {path.name}: missing {missing}")
        return list(reader)


def literal_list(value, *, field="CCTyper list"):
    """Parse a Python-literal list field written by pandas/CCTyper."""
    try:
        result = ast.literal_eval(value)
    except (ValueError, SyntaxError) as exc:
        raise ValueError(f"Malformed {field}: {value!r}") from exc
    if not isinstance(result, (list, tuple)):
        raise ValueError(f"Expected a list in {field}, got {type(result).__name__}")
    return list(result)


def safe_id(value, *, what="identifier"):
    """Reject identifiers that could escape the annotation directory."""
    text = str(value)
    if not text or set(text) - _SAFE_ID or text in {".", ".."}:
        raise ValueError(f"Unsafe {what}: {text!r}")
    return text


def _as_int(value, *, what):
    try:
        return int(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Non-integer {what}: {value!r}") from exc


def translate(dna, table=STANDARD_TABLE):
    """Translate a nucleotide string, ignoring a trailing partial codon."""
    return "".join(
        table.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3)
    )


def reverse_complement(dna):
    return dna.translate(_COMPLEMENT)[::-1]


def verify_protein(
    sequence,
    *,
    contig_sequence=None,
    start=None,
    end=None,
    strand=None,
):
    """Check a reported Cas protein against the genome it was called from.

    Returns ``(verified: bool, notes: list[str])``.  ``verified`` is only True
    when the protein can be reproduced from the genome coordinates.  Anything
    that cannot be confirmed is reported as unverified rather than assumed
    complete -- the caller must not silently promote it to an eligible system.
    """
    notes: list[str] = []
    if not sequence:
        return False, ["protein_sequence_missing"]

    # CCTyper strips the stop, so an internal stop is always a real defect.
    if "*" in sequence:
        notes.append("internal_stop_codon")
    if set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        notes.append("nonstandard_residues")

    if contig_sequence is None or start is None or end is None or strand is None:
        notes.append("no_genome_coordinates")
        return False, notes

    contig_length = len(contig_sequence)
    if not 1 <= start <= end <= contig_length:
        notes.append("gene_coordinates_out_of_range")
        return False, notes

    if str(strand) not in {"1", "-1", "+", "-"}:
        notes.append("unknown_strand")
        return False, notes
    forward = str(strand) in {"1", "+"}

    nucleotide = contig_sequence[start - 1:end]
    coding = nucleotide if forward else reverse_complement(nucleotide)

    if len(coding) % 3 != 0:
        notes.append("gene_length_not_multiple_of_three")
    if coding[:3] not in START_CODONS:
        notes.append("no_valid_start_codon")
    if coding[-3:] not in STOP_CODONS:
        notes.append("no_valid_stop_codon")

    translated = translate(coding)
    without_stop = translated.rstrip("*")

    def blocking():
        """A defect that makes the called gene boundary untrustworthy.

        A missing start or stop codon means the reported coordinates do not
        describe a complete gene, even if the residues happen to line up, so it
        blocks verification just like an internal stop does.
        """
        return any(n in notes for n in BLOCKING_DEFECTS)

    if without_stop == sequence:
        notes.append("translation_matches_exactly")
        return not blocking(), notes

    # pyrodigal rewrites the first residue to M when an alternative start codon
    # (GTG/TTG/...) is used, mapping to methionine in the reported protein.
    if (
        len(without_stop) == len(sequence)
        and without_stop[1:] == sequence[1:]
        and sequence[0] == "M"
        and coding[:3] in START_CODONS
    ):
        notes.append("translation_matches_with_alternative_start_codon")
        return not blocking(), notes

    notes.append("translation_mismatch")
    return False, notes


def collect_systems(directory, genome, *, genome_path=None, sequences=None):
    """Extract candidate Cas9 systems from a CCTyper output directory.

    ``genome`` is the summary from :func:`check_genome`.  ``sequences`` (or
    ``genome_path``) supplies contig sequences for coordinate-level checks.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError(f"CCTyper output directory not found: {directory}")

    if sequences is None and genome_path is not None:
        sequences, _ = load_sequences(genome_path)
    if sequences is None:
        sequences = {}

    # Read every operon table CCTyper may write, because which one exists depends
    # on the genome: good operons go to cas_operons.tab, an ambiguous subtype to
    # cas_operons_putative.tab (with no good table at all), and a Cas operon with
    # no nearby array to cas_operons_orphan.tab.  Reading only the good table
    # silently drops real Cas9 genes.
    operon_columns = ["Contig", "Operon", "Prediction", "Genes", "Positions"]
    operons = table(directory / "cas_operons.tab", operon_columns)
    for row in operons:
        row["_table"] = "cas_operons.tab"
    seen_operons = {row["Operon"] for row in operons}
    for name in ("cas_operons_putative.tab", "cas_operons_orphan.tab"):
        extra = table(directory / name, operon_columns)
        for row in extra:
            row["_table"] = name
            if row["Operon"] not in seen_operons:
                operons.append(row)
                seen_operons.add(row["Operon"])

    loci = table(directory / "CRISPR_Cas.tab", ["Operon", "CRISPRs"])
    arrays = table(
        directory / "crisprs_all.tab",
        ["CRISPR", "Contig", "Start", "End", "Consensus_repeat", "Trusted"],
    )
    genes = table(
        directory / "genes.tab", ["Contig", "Pos", "Start", "End", "Strand"]
    )
    gene_map = {(r["Contig"], str(r["Pos"])): r for r in genes}

    if not operons:
        return []

    proteins_path = directory / "proteins.faa"
    if not proteins_path.is_file():
        raise ValueError("CCTyper did not write proteins.faa")
    proteins, protein_headers = fasta(
        proteins_path,
        "ACDEFGHIKLMNPQRSTVWYBXZJUO*",
        allow_duplicates=True,
    )

    systems = []
    for op in operons:
        contig = op["Contig"]
        operation = op["Operon"]
        names = literal_list(op["Genes"], field="Genes")
        positions = literal_list(op["Positions"], field="Positions")
        if len(names) != len(positions):
            raise ValueError(
                f"Cas gene / position length mismatch in operon {operation}"
            )
        prediction = str(op["Prediction"])
        for name, pos in zip(names, positions):
            if not str(name).lower().startswith(("cas9", "csn1")):
                continue
            position = str(pos)
            system_id = f"{contig}_{position}"
            sequence = proteins.get(system_id, "").rstrip("*")
            gene = gene_map.get((contig, position))

            issues: list[str] = []
            if prediction.startswith("II-"):
                pass
            elif prediction in {"False", "Ambiguous", "Unknown"} or "Putative" in prediction:
                issues.append("ambiguous_or_putative_subtype")
            else:
                issues.append("unsupported_subtype")

            if not sequence:
                issues.append("missing_protein_sequence")

            if gene is None:
                issues.append("gene_coordinates_missing")
            if contig not in genome["contig_lengths"]:
                issues.append("contig_not_in_input_genome")

            verified = False
            notes: list[str] = []
            if gene is not None and contig in sequences:
                verified, notes = verify_protein(
                    sequence,
                    contig_sequence=sequences[contig],
                    start=_as_int(gene["Start"], what="gene Start"),
                    end=_as_int(gene["End"], what="gene End"),
                    strand=gene["Strand"],
                )
            if not verified:
                issues.append("gene_integrity_unverified")

            gene_info = None
            if gene is not None:
                gene_info = {
                    "contig": contig,
                    "pos": position,
                    "start": _as_int(gene["Start"], what="gene Start"),
                    "end": _as_int(gene["End"], what="gene End"),
                    "strand": str(gene["Strand"]),
                }

            linked = set()
            for row in loci:
                if row["Operon"] == operation:
                    linked.update(str(x) for x in literal_list(row["CRISPRs"], field="CRISPRs"))

            linked_arrays = []
            for a in arrays:
                if a["CRISPR"] not in linked or a["Contig"] != contig:
                    continue
                array_id = safe_id(a["CRISPR"], what="array identifier")
                start = _as_int(a["Start"], what="array Start")
                end = _as_int(a["End"], what="array End")
                array_issues = []
                contig_length = genome["contig_lengths"].get(contig)
                if contig_length is not None and not 1 <= start <= end <= contig_length:
                    array_issues.append("array_coordinates_out_of_range")
                spacer_file = directory / "spacers" / f"{array_id}.fa"
                if not spacer_file.is_file():
                    array_issues.append("spacer_file_missing")
                repeat = str(a["Consensus_repeat"]).upper()
                if not repeat or set(repeat) - set("ACGT"):
                    array_issues.append("invalid_consensus_repeat")
                linked_arrays.append({
                    "array_id": array_id,
                    "trusted": str(a["Trusted"]).strip().lower() == "true",
                    "start": start,
                    "end": end,
                    "repeat": repeat,
                    "repeat_length": len(repeat),
                    "spacer_file": str(spacer_file),
                    "spacer_file_present": spacer_file.is_file(),
                    "orientation": "unconfirmed",
                    "issues": array_issues,
                })

            systems.append({
                "system_id": system_id,
                "gene_name": str(name),
                "operon": operation,
                "operon_table": op["_table"],
                "prediction": prediction,
                "subtype": prediction,
                "operon_start": _as_int(op["Start"], what="operon Start") if op.get("Start") else None,
                "operon_end": _as_int(op["End"], what="operon End") if op.get("End") else None,
                "contig": contig,
                "gene": gene_info,
                "protein_length": len(sequence),
                "protein_sequence": sequence,
                "protein_header": protein_headers.get(system_id, ""),
                "integrity_verified": verified,
                "integrity_notes": notes,
                "issues": issues,
                "arrays": linked_arrays,
            })
    return systems


def select_system(systems, requested=None, *, allow_flagged=False):
    """Choose a Cas9 system, or report why an explicit choice is needed.

    With no ``requested`` ID only a single *fully eligible* system is returned
    automatically; anything ambiguous falls through to the caller for a user
    decision.

    With ``requested`` set, the user has inspected the reported issues and is
    making an explicit, informed choice.  ``allow_flagged`` therefore accepts a
    system whose only problems are annotation-ambiguity flags, while always
    refusing one that failed protein-integrity verification: an unverifiable
    protein must never be sent to the model merely because it was named.
    """
    if requested:
        matches = [s for s in systems if s["system_id"] == requested]
        if not matches:
            raise ValueError(f"Unknown --system-id: {requested}")
        selected = matches[0]
        if not selected["issues"]:
            return selected
        if not allow_flagged:
            raise ValueError(
                f"--system-id {requested} failed quality checks: "
                f"{', '.join(selected['issues'])}"
            )
        blocking = [i for i in selected["issues"] if i in BLOCKING_ISSUES]
        if blocking:
            raise ValueError(
                f"--system-id {requested} cannot be used even with an explicit "
                f"override; its protein did not pass integrity verification: "
                f"{', '.join(blocking)}"
            )
        return selected
    eligible = [s for s in systems if not s["issues"]]
    if len(eligible) == 1:
        return eligible[0]
    return None
