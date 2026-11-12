"""The genotype table gives exactly what the original 15-branch if/elif chain gave, for every combination of the five
primary-marker states and every outcome of the A-subtype step.

legacy_chain() is the chain as it was in aggregate_abo_reports.py before the table (nf-core/abotyper 9adba82), kept
here verbatim as the reference.
"""
import itertools

import pytest

import aggregate_abo_reports as agg

STATES = ("ref", "het", "alt", "none")
# how the original code turned each marker's state into the strings the chain matched (primary_type)
LEGACY_STRINGS = [
    {"ref": "A or B or O", "alt": "O1", "het": "O1 and (A or B or O)", "none": ""},          # c.261delG
    {"ref": "A or O", "alt": "B", "het": "(A or O) and B", "none": ""},                      # c.796C>A
    {"ref": "O and (A or B)", "alt": "O2", "het": "O2 and (O or A or B)", "none": ""},       # c.802G>A
    {"ref": "A or O", "alt": "B", "het": "(A or O) and B", "none": ""},                      # c.803G>C
    {"ref": "O and (A or B)", "alt": "O3", "het": "O3 and (O or A or B)", "none": ""},       # c.804
]
SUBTYPE_OUTCOMES = [("", None, "hom"), ("A2.01", "warn", "het"), ("A2.01", None, "hom"), ("A1.02", "w2", "het"), ("A1", None, "hom")]


def legacy_chain(type_exon6, type_exon7_422, type_exon7_428, type_exon7_429, type_exon7_431, determine_a_subtype):
    a_subtype_warning = None

    if (type_exon6 == "O1 and (A or B or O)" and type_exon7_422 == "A or O"
            and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
            and type_exon7_431 == "O and (A or B)"):
        a_subtype, a_subtype_warning, _ = determine_a_subtype()
        Phenotype, Genotype = "A", "AO"
        ExtendedGenotype = f"{a_subtype}/O1" if a_subtype else "A/O1"

    elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
            and type_exon7_428 == "O2 and (O or A or B)" and type_exon7_429 == "A or O"
            and type_exon7_431 == "O and (A or B)"):
        a_subtype, a_subtype_warning, _ = determine_a_subtype()
        Phenotype, Genotype = "A", "AO"
        ExtendedGenotype = f"{a_subtype}/O2" if a_subtype else "A/O2"

    elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
            and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
            and type_exon7_431 == "O3 and (O or A or B)"):
        a_subtype, a_subtype_warning, _ = determine_a_subtype()
        Phenotype, Genotype = "A", "AO"
        ExtendedGenotype = f"{a_subtype}/O3" if a_subtype else "A/O3"

    elif (type_exon6 == "O1 and (A or B or O)" and type_exon7_422 == "(A or O) and B"
            and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "(A or O) and B"
            and type_exon7_431 == "O and (A or B)"):
        Phenotype, Genotype, ExtendedGenotype = "B", "BO", "B/O1"

    elif (type_exon6 == "A or B or O" and type_exon7_422 == "(A or O) and B"
            and type_exon7_428 == "O2 and (O or A or B)" and type_exon7_429 == "(A or O) and B"
            and type_exon7_431 == "O and (A or B)"):
        Phenotype, Genotype, ExtendedGenotype = "B", "BO", "B/O2"

    elif (type_exon6 == "A or B or O" and type_exon7_422 == "(A or O) and B"
            and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "(A or O) and B"
            and type_exon7_431 == "O3 and (O or A or B)"):
        Phenotype, Genotype, ExtendedGenotype = "B", "BO", "B/O3"

    elif (type_exon6 == "O1 and (A or B or O)" and type_exon7_422 == "A or O"
            and type_exon7_428 == "O2 and (O or A or B)" and type_exon7_429 == "A or O"
            and type_exon7_431 == "O and (A or B)"):
        Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O1/O2"

    elif (type_exon6 == "O1 and (A or B or O)" and type_exon7_422 == "A or O"
            and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
            and type_exon7_431 == "O3 and (O or A or B)"):
        Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O1/O3"

    elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
            and type_exon7_428 == "O2 and (O or A or B)" and type_exon7_429 == "A or O"
            and type_exon7_431 == "O3 and (O or A or B)"):
        Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O2/O3"

    elif (type_exon6 == "O1" and type_exon7_422 == "A or O"
            and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
            and type_exon7_431 == "O and (A or B)"):
        Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O1/O1"

    elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
            and type_exon7_428 == "O2" and type_exon7_429 == "A or O"
            and type_exon7_431 == "O and (A or B)"):
        Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O2/O2"

    elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
            and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
            and type_exon7_431 == "O3"):
        Phenotype, Genotype, ExtendedGenotype = "O", "OO", "O3/O3"

    elif (type_exon6 == "A or B or O" and type_exon7_422 == "A or O"
            and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "A or O"
            and type_exon7_431 == "O and (A or B)"):
        a_subtype, a_subtype_warning, a_subtype_zygosity = determine_a_subtype()
        if a_subtype and a_subtype_zygosity == "het":
            Phenotype, Genotype, ExtendedGenotype = "A", "AA", f"A1/{a_subtype}"
        elif a_subtype:
            Phenotype, Genotype, ExtendedGenotype = a_subtype, "AA", f"{a_subtype}/{a_subtype}"
        else:
            Phenotype, Genotype, ExtendedGenotype = "A", "AA", "A/A"

    elif (type_exon6 == "A or B or O" and type_exon7_422 == "B"
            and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "B"
            and type_exon7_431 == "O and (A or B)"):
        Phenotype, Genotype, ExtendedGenotype = "B", "BB", "B/B"

    elif (type_exon6 == "A or B or O" and type_exon7_422 == "(A or O) and B"
            and type_exon7_428 == "O and (A or B)" and type_exon7_429 == "(A or O) and B"
            and type_exon7_431 == "O and (A or B)"):
        a_subtype, a_subtype_warning, _ = determine_a_subtype()
        Phenotype, Genotype = "AB", "AB"
        ExtendedGenotype = f"{a_subtype}/B" if a_subtype else "A/B"

    else:
        Phenotype, Genotype, ExtendedGenotype = "Unknown", "Unknown", "Unknown"
    return Phenotype, Genotype, ExtendedGenotype, a_subtype_warning


@pytest.mark.parametrize("outcome", SUBTYPE_OUTCOMES)
def test_table_matches_the_original_chain(outcome):
    for states in itertools.product(STATES, repeat=5):
        calls = {"legacy": 0, "table": 0}

        def subtype(which):
            def f():
                calls[which] += 1
                return outcome
            return f

        strings = [m[s] for m, s in zip(LEGACY_STRINGS, states)]
        expected = legacy_chain(*strings, subtype("legacy"))
        got = agg.resolve_primary_genotype(states, subtype("table"))
        assert got == expected, (states, outcome)
        # the A-subtype scan has side effects (it records fired markers), so it must run in exactly the same cases
        assert calls["table"] == calls["legacy"], states


def test_table_rows_are_reachable():
    assert len(agg.GENOTYPE_TABLE) == 15
    assert all(len(k) == 5 and set(k) <= {"ref", "het", "alt"} for k in agg.GENOTYPE_TABLE)
