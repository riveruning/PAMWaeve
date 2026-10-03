"""PAMdict unified schema.

This is the single source of truth for the PAM training-data contract. Every
parser, scorer, and finetuner in this production line reads/writes this schema.

The ten core columns mirror the official Protein2PAM training TSV
(``protein2pam_train_seqs.tsv``) so that new data can be merged with the
official baseline. See the prior provenance audit in the reference project
(PAM_PREDICT) for the origin of these field facts.
"""

from .record import PAMRecord, pam_logo_from_matrix, pam_matrix_from_logo
from .constants import (
    PAM_NUCLEOTIDES,   # ["A", "C", "G", "T"]
    PAM_POSITIONS,     # list(range(1, 11))
    CORE_COLUMNS,      # official 10-column order
    LOGO_UNITS,        # "information content (bits)"
)

__all__ = [
    "PAMRecord",
    "pam_logo_from_matrix",
    "pam_matrix_from_logo",
    "PAM_NUCLEOTIDES",
    "PAM_POSITIONS",
    "CORE_COLUMNS",
    "LOGO_UNITS",
]
