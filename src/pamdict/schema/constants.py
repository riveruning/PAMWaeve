"""Fixed constants for the PAMdict schema.

These encode the PAM representation contract that the reference audit
established from Protein2PAM's own consumers:

- ``pams.py``      : channel order ``["A","C","G","T"]``, positions 1..10.
- ``oracle.py``    : ``predictions = prob_to_info(...)`` -> units are
                     information content (bits), NOT 0..1 probabilities.
- ``common.py``    : ``prob_to_info`` definition (per-position weighted info).
- ``config.json``  : ``output_shape: [10, 4]``, ``id2label`` A/C/G/T.

Any code that compares a matrix value against the 0.70 threshold must use
the *information content* unit, or first convert explicitly.
"""

from __future__ import annotations

PAM_NUCLEOTIDES: tuple[str, ...] = ("A", "C", "G", "T")
PAM_POSITIONS: list[int] = list(range(1, 11))  # 1..10

# Official training TSV column order (verified via range request, 2026-09-05).
CORE_COLUMNS: tuple[str, ...] = (
    "crispr_type",
    "cas_family",
    "source",
    "protein_id",
    "citation",
    "doi",
    "protein_sequence",
    "pid_sequence",
    "pam_consensus",
    "pam_logo_acgt",
)

# Units of ``pam_logo_acgt`` as consumed by the reference implementation.
LOGO_UNITS: str = "information content (bits)"

# Extension columns appended by this production line (do NOT reorder/adjust the
# core columns; extension columns are additive and record provenance).
EXTENSION_COLUMNS: tuple[str, ...] = (
    "created_at",
    "raw_ref",
    "sha256",
)

ALL_COLUMNS: tuple[str, ...] = CORE_COLUMNS + EXTENSION_COLUMNS

# Consensus threshold used by the reference implementation (bits).
CONSENSUS_THRESHOLD_BITS: float = 0.70
