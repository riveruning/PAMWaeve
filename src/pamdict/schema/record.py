"""PAM record model and matrix<->logo conversion helpers.

The schema is a plain dict/dataclass mapping to the core 10 columns (plus
extension columns) from the official Protein2PAM training TSV. We intentionally
avoid a heavy schema library so that the parser/collector stays dependency-light
and the corpus can be round-tripped as TSV/Parquet without surprises.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any

from .constants import (
    ALL_COLUMNS,
    CORE_COLUMNS,
    LOGO_UNITS,
    PAM_NUCLEOTIDES,
    PAM_POSITIONS,
)


@dataclass
class PAMRecord:
    """One (protein -> PAM) training sample.

    Core fields mirror the official 10-column TSV. ``pam_logo_acgt`` is stored
    as a JSON array of shape [10 x 4] (row-major, positions then A/C/G/T) whose
    values are information content in bits.
    """

    crispr_type: str
    cas_family: str
    source: str
    protein_id: str
    citation: str
    doi: str
    protein_sequence: str
    pid_sequence: str
    pam_consensus: str
    pam_logo_acgt: Any  # list[list[float]] shape [10, 4] in bits

    # Extension columns (provenance), optional.
    created_at: str = ""
    raw_ref: str = ""
    sha256: str = ""

    def to_dict(self, include_extensions: bool = True) -> dict[str, Any]:
        cols = ALL_COLUMNS if include_extensions else CORE_COLUMNS
        return {c: getattr(self, c) for c in cols}

    @property
    def matrix(self) -> list[list[float]]:
        """Return the 10x4 info-content matrix (deep copy)."""
        m = json.loads(self.pam_logo_acgt) if isinstance(self.pam_logo_acgt, str) else self.pam_logo_acgt
        return [list(row) for row in m]


def pam_matrix_from_logo(logo_json: str | list) -> list[list[float]]:
    """Parse a ``pam_logo_acgt`` value into a 10x4 matrix."""
    raw = json.loads(logo_json) if isinstance(logo_json, str) else logo_json
    rows = [list(map(float, row)) for row in raw]
    assert len(rows) == len(PAM_POSITIONS), f"expected {len(PAM_POSITIONS)} positions, got {len(rows)}"
    assert all(len(r) == len(PAM_NUCLEOTIDES) for r in rows), "each position must have A/C/G/T"
    return rows


def pam_logo_from_matrix(matrix: list[list[float]]) -> str:
    """Serialize a 10x4 matrix back to the official logo JSON array string."""
    return json.dumps([[float(v) for v in row] for row in matrix])
