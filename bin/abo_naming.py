#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ISBT allele naming: internal allele tokens (A1, A2.01, O1, ...) to ISBT-style names, narrowed by the panel markers that fired or were ruled out.

Part of aggregate_abo_reports.py (nf-core/abotyper); split out of that script so each part can be read and
tested on its own. This file is part of nf-core/abotyper and is released under the MIT License (see LICENSE).
"""
from __future__ import annotations

import re





ISBT_ALLELE_NAMES = {
    "A1": "ABO*A1.01",
    "A1.02": "ABO*A1.02",
    "A2.01": "ABO*A2.01",
    "A2.06": "ABO*A2.06",
    "A2": "ABO*A2",
    "A": "ABO*A",
    "B": "ABO*B.01",
    "O1": "ABO*O.01.01",
    "O2": "ABO*O.02.01",
    "O3": "ABO*O.03",
}


# Family prefix an allele name must start with to be a candidate for a token.
_ISBT_FAMILY = {
    "A1": "ABO*A1.", "A1.02": "ABO*A1.", "A2.01": "ABO*A2.", "A2.06": "ABO*A2.", "A2": "ABO*A2.",
    "B": "ABO*B.", "O1": "ABO*O.01.", "O2": "ABO*O.02.", "O3": "ABO*O.03",
}


_ISBT_NAME_RE = re.compile(r"(?:ABO\*)?([A-Za-z]+[0-9]*\.\d+(?:\.\d+)?)(?![\d-])")


_ALT_NOTE_RE = re.compile(r"\(c\.\d+(?:[ACGT]>([ACGT])|(dup|del)[ACGT]*)\)")


# Ael alleles sit outside the A1/A2 families but are refinements of them: c.804dupG alone is ABO*AEL.01,
# with c.467C>T it is ABO*AEL.08 (ISBT v9); c.804delG alone is ABO*AEL.03
_EXTRA_FAMILY = {"A1": ("ABO*AEL.01", "ABO*AEL.03"), "A1.02": ("ABO*AEL.08",)}


def _marker_candidate_alleles(marker, alt=None) -> set:
    """Exact (non-wildcard) ISBT allele names listed in a panel marker's associated_alleles.

    A marker with several alt bases can annotate each entry with the variant it belongs to, e.g.
    'ABO*A3.04 (c.539G>A)'. When the alt base that actually fired is known, entries annotated with a
    different alt base are skipped; unannotated entries always count.
    """
    out = set()
    for entry in (getattr(marker, "associated_alleles", "") or "").split(";"):
        names = ["ABO*" + m.group(1).upper() for m in _ISBT_NAME_RE.finditer(entry)]
        if not names:
            continue
        note = _ALT_NOTE_RE.search(entry)
        if alt and note and (note.group(1) or note.group(2)).upper() != str(alt).upper():
            continue
        out.update(names)
    return out


def _unpack_fired(mk):
    """Fired-marker entries are (marker, alt_base) pairs; plain markers are accepted too."""
    return mk if isinstance(mk, tuple) else (mk, None)


# What the primary A markers say about the A allele itself: (c.467C>T, c.1061delC) carried or not. A candidate allele
# whose ISBT definition disagrees (A2.20 is c.829G>A alone, so it cannot sit on an A2.01 background) is dropped.
_SAME_FAMILY_SPLIT = ("O1", "O2")


def _cdna_tokens(cdna: str):
    """The c. changes of a panel row ('c.539G>C / c.539G>A' -> two tokens); a descriptive row such as
    'c.804 region (homopolymer, legacy O3/O4 marker)' gives its leading 'c.804'."""
    out = []
    for tk in re.split(r"\s*/\s*|,\s*(?=c\.)", cdna or ""):
        m = re.match(r"c\.[^\s,]+", tk.strip())
        if m:
            out.append(m.group(0))
    return out


_A_PRIMARY_CONTEXT = {"A2.01": (True, True), "A2.06": (False, True), "A1.02": (True, False), "A1": (False, False)}


def primary_carriers(panel):
    """(alleles carrying c.467C>T, alleles carrying c.1061delC) from the two primary markers' associated alleles."""
    out = []
    for role in ("a1a2_467", "a2_1061del"):
        m = panel.role(role) if panel is not None else None
        out.append(_marker_candidate_alleles(m) if m is not None else None)
    return tuple(out)


def encode_allele(token: str, fired_markers=(), absent_markers=(), carriers=None, pooled: bool = True, other_token: str = "") -> str:
    """Internal allele token (A1, O2, ...) -> ISBT-style name.

    Starts from the reference allele of the token's family, then narrows using the
    non-primary panel markers that actually fired in this sample: the allele names
    those markers are associated with, restricted to the token's family. Several
    candidates are narrowed further by dropping every allele that is defined by a marker
    which was scanned with enough coverage and did not fire (absent_markers), e.g.
    ABO*O.01.75 needs c.106G>T, so it is ruled out when c.106 reads as reference.
    One surviving candidate -> that allele; several -> the family name only (not resolved
    further); none -> the reference allele. Unknown tokens pass through.
    """
    t = str(token).strip()
    default = ISBT_ALLELE_NAMES.get(t, t)
    fam = _ISBT_FAMILY.get(t)
    if not fam:
        return default
    sets = []
    for mk in fired_markers:
        marker, alt = _unpack_fired(mk)
        cand = {a for a in _marker_candidate_alleles(marker, alt)
                if a.startswith(fam) or a == fam.rstrip(".") or a in _EXTRA_FAMILY.get(t, ())}
        if cand:
            sets.append(cand)
    if not sets:
        return default
    common = set.intersection(*sets)
    cands = common if common else set.union(*sets)
    ctx = _A_PRIMARY_CONTEXT.get(t)
    if ctx and carriers and all(c is not None for c in carriers):
        consistent = {a for a in cands if (a in carriers[0]) == ctx[0] and (a in carriers[1]) == ctx[1]}
        cands = consistent or cands
        if pooled:
            # the markers of both alleles were pooled, so a marker may belong to the other allele. When every marker
            # behind the candidates also lists alleles of the other allele's family (c.829G>A is on A2.20 and on the
            # O.01 cassette alleles), the background allele itself (A2.01 for an A2.01 token) stays a candidate
            other_fam = _ISBT_FAMILY.get(str(other_token).strip(), "")
            if other_fam:
                behind = [(mk, alt) for mk, alt in (_unpack_fired(x) for x in fired_markers)
                          if {a for a in _marker_candidate_alleles(mk, alt) if a.startswith(fam)}]
                if behind and all(any(a.startswith(other_fam) for a in _marker_candidate_alleles(mk, alt))
                                  for mk, alt in behind):
                    cands = cands | {default}
    if len(cands) > 1 and absent_markers:
        ruled_out = set()
        for mk in absent_markers:
            ruled_out |= _marker_candidate_alleles(mk)
        cands = (cands - ruled_out) or cands
    if len(cands) == 1:
        return next(iter(cands))
    return fam.rstrip(".")   # ambiguous: family only, e.g. ABO*O.01


def encode_extended_genotype(ext: str, fired_markers=(), absent_markers=(), carriers=None, fired_pct=None) -> str:
    """'A2.01/O1' -> 'ABO*A2.01/ABO*O.01.01'. 'Unknown', 'Error' and empty values are left untouched."""
    if not isinstance(ext, str) or ext in ("", "Unknown", "Error"):
        return ext
    tokens = ext.split("/")
    markers = list(fired_markers)
    fired_ids = {getattr(_unpack_fired(m)[0], "id", None) for m in markers}
    absent = [m for m in absent_markers if getattr(m, "id", None) not in fired_ids]
    if len(tokens) == 2 and tokens[0] == tokens[1] and tokens[0] in _SAME_FAMILY_SPLIT and fired_pct is not None:
        # Two alleles of the same family (O1/O1) and no phasing. A marker at a heterozygous fraction sits on one of
        # the two alleles, so it must not rename both; a homozygous marker is on both. The heterozygous markers are
        # given to one allele (assumed in cis) and the other allele is named from the homozygous markers alone.
        hom = [m for m in markers if (fired_pct.get(getattr(_unpack_fired(m)[0], "id", None)) or 0) >= 80]
        het = [m for m in markers if m not in hom]
        if het:
            tok = tokens[0]
            with_het = encode_allele(tok, hom + het, absent, carriers, pooled=False)
            without = encode_allele(tok, hom, absent + [_unpack_fired(m)[0] for m in het], carriers, pooled=False)
            return "/".join(sorted([with_het, without]))
    if "B" in tokens:
        # B-lineage backbone SNPs (associated alleles include ABO*B.*) come from the B
        # allele of an AB/BO sample; they must not rename the non-B allele (e.g. to A2.04)
        backbone = [m for m in markers if "ABO*B.*" in (getattr(_unpack_fired(m)[0], "associated_alleles", "") or "")]
        non_b = [m for m in markers if m not in backbone]
        return "/".join(encode_allele(t, markers if t == "B" else non_b, absent, carriers,
                                      other_token=tokens[1 - i] if len(tokens) == 2 else "")
                        for i, t in enumerate(tokens))
    return "/".join(encode_allele(t, markers, absent, carriers, other_token=tokens[1 - i] if len(tokens) == 2 else "")
                    for i, t in enumerate(tokens))
