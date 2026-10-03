"""Curated CAS alias -> source organism mapping (offline, conservative).

Why this exists
---------------
The `species` field in `cas_pam_pairs.tsv` is a co-occurrence string from
full-text; it mixes the *true* source organism with *host / target* organisms
("Mus musculus", "Escherichia coli", "Oryza sativa", "Arabidopsis thaliana") and
generic descriptors ("In vitro", "Lachnospiraceae bacterium"). Using that field
directly as the UniProt anchor produces wrong protein sequences.

The reliable anchor is the Cas *alias* itself, whose source organism is fixed
domain knowledge (e.g. "SpCas9" == Streptococcus pyogenes Cas9). This module
records exactly that: alias -> (source organism, family), **only for entries
that are uncontroversial**. Anything we are not confident about is intentionally
left out (resolved via live UniProt query under explicit review instead).

``VERIFIED_ACCESSION`` records alias -> UniProt accession that we have actually
resolved against UniProt and sanity-checked (length, family, organism). These
are safe to use offline.
"""

from __future__ import annotations

# alias -> source organism.
#
# ONLY the uncontroversial, universally-agreed CRISPR names are listed here.
# Anything we are not 100% sure about is intentionally ABSENT — it must be
# resolved via a live UniProt query and (ideally) human review, never guessed.
# Short prefixes are ambiguous by nature: "CdCas9" vs "SauriCas9" vs "CbCas9"
# are exactly the kind of alias that is easy to confuse with a lookalike, so
# they are not hard-coded here.
ALIAS_ORGANISM: dict[str, str] = {
    # Type II-A Cas9 (unambiguous, consensus names)
    "SpCas9": "Streptococcus pyogenes",
    "SpyCas9": "Streptococcus pyogenes",
    "dSpCas9": "Streptococcus pyogenes",  # catalytically-dead SpCas9
    "SaCas9": "Staphylococcus aureus",
    "dSaCas9": "Staphylococcus aureus",
    "FnCas9": "Francisella novicida",
    "CjCas9": "Campylobacter jejuni",
    "NmCas9": "Neisseria meningitidis",
    "NmeCas9": "Neisseria meningitidis",
    "StCas9": "Streptococcus thermophilus",
    # Type V-A / Cas12a (Cpf1)
    "AsCas12a": "Acidaminococcus sp. BV3L6",
    "AsCpf1": "Acidaminococcus sp. BV3L6",
    "FnCas12a": "Francisella novicida",
    "FnCpf1": "Francisella novicida",
    "LbCas12a": "Lachnospiraceae bacterium ND2006",
    "LbCpf1": "Lachnospiraceae bacterium ND2006",
    # Cold-start aliases whose source organism was confirmed in-session
    # (2026-09), fixing earlier paper-species mis-attribution. These are the
    # uncontroversial source organisms discovered during the out-of-domain
    # generalization eval — see docs/P2PAM_NEW_DATA_EVAL.md §4.
    "CoCas9": "Capnocytophaga ochracea",
    "SmacCas9": "Capnocytophaga canimorsus",
    "SmutCas9": "Streptococcus mutans",
    "SaHyCas9": "Staphylococcus hyicus",
    "CdCas9": "Corynebacterium diphtheriae",
    "DfCas9": "Rodentibacter pneumotropicus",
    "ScCas9": "Streptococcus macacae",
    # 2026-09-17: newly added source organisms (live-verified via UniProt)
    "MbCas12a": "Moraxella bovis",
    "EbCas12a": "Eubacterium album",
}

# alias -> UniProt accession, *verified* against UniProt (length+family+organism
# sanity-checked in-session). These are safe to use without further queries.
VERIFIED_ACCESSION: dict[str, str] = {
    "SpCas9": "Q99ZW2",
    # SaCas9: J7RUA5 (1053 aa, Swiss-Prot, S. aureus) — re-verified 2026-09-17
    # via live UniProt query. The previous value A0A386IRL6 was the contaminated
    # entry flagged in docs/DATA_QUALITY_ALERT_A0A386IRL6.md (sequence identical
    # to SpCas9 Q99ZW2 despite the S. aureus label); see also the fix applied by
    # scripts/fix_a0a386irl6.py to the data files.
    "SaCas9": "J7RUA5",
    "FnCas9": "A0Q5Y3",
    "NmeCas9": "C9X1G5",
    "CjCas9": "Q0P897",
    "StCas9": "G3ECR1",
    "AsCas12a": "U2UMQ6",
    "FnCas12a": "A0Q7Q2",
    # 2026-09-17: alias-equivalence verifications (live UniProt queries, same
    # protein under a different name — zero-risk additions).
    "AsCpf1": "U2UMQ6",     # same protein as AsCas12a (Acidaminococcus sp. BV3L6, 1307 aa)
    "FnCpf1": "A0Q7Q2",     # same protein as FnCas12a (F. novicida U112, 1300 aa)
    "LbCpf1": "C4Z1Q1",     # Lachnospira eligens Cpf1, 1282 aa (LbCas12a family ortholog)
    "Nme1Cas9": "C9X1G5",   # same protein as NmeCas9 (N. meningitidis 8013, 1082 aa)
    # 2026-09-17: newly resolved natural orthologs (Moraxella/Eubacterium Cas12a)
    "MbCas12a": "A0AAQ2Q2F9",  # Moraxella bovis Cas12a, 1261 aa (TrEMBL)
    "EbCas12a": "A0ABT2M290",  # Eubacterium album Cas12a, 1260 aa (TrEMBL; candidate pair with A0ABT2M0G8 1175 aa — PMC11139299 mapping needs confirm)
    # Cold-start aliases verified against UniProt in-session (length+family+
    # organism sanity-checked); fixes earlier paper-species mis-resolution.
    "CoCas9": "C7M7G9",     # Capnocytophaga ochracea Cas9, 1426 aa
    "SmacCas9": "F9YQX1",   # Capnocytophaga canimorsus Cas9, 1430 aa (was wrongly C7M7G9)
    "SmutCas9": "Q8DTE3",   # Streptococcus mutans Cas9, 1345 aa
    "SaHyCas9": "A0ACD5FQ62",  # Staphylococcus hyicus Cas9, 1055 aa
    "CdCas9": "Q6NKI3",     # Corynebacterium diphtheriae Cas9, 1084 aa (was wrongly S. auricularis)
    "DfCas9": "A0AAW5LCH4", # Rodentibacter pneumotropicus Cas9, 1055 aa
    "ScCas9": "G5JVJ9",     # Streptococcus macacae Cas9, 1338 aa (PAM itself unreliable)
}


def family_of(alias: str) -> str | None:
    """Infer the Cas family keyword from an alias (same logic as resolve_uniprot)."""
    n = alias.lower().strip()
    if "cas12f" in n:
        return "cas12f"
    if "cas12" in n or "cpf1" in n:
        return "cas12"
    if "cas13" in n or "c2c2" in n:
        return "cas13"
    if "cas9" in n:
        return "cas9"
    if "iscb" in n:
        return "iscb"
    return None


def organism_of(alias: str) -> str | None:
    """Return the source organism for a Cas alias, or None if unknown."""
    return ALIAS_ORGANISM.get(alias)
