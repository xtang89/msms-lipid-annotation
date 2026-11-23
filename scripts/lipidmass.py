"""Shared mass functions for the MS/MS annotation skill.

All masses are monoisotopic. Formulas are written without charge (for example C5H15NO4P).
The electron mass is applied when a charge is given.
"""
import re

ELEMENTS = {
    "C": 12.0,
    "H": 1.00782503207,
    "N": 14.0030740048,
    "O": 15.99491461956,
    "P": 30.97376163,
    "S": 31.97207100,
    "Na": 22.9897692809,
    "K": 38.96370668,
    "Li": 7.01600455,
}
ELECTRON = 0.00054858

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula):
    """Return a dict element -> count. Raise ValueError for unknown symbols or extra text."""
    counts = {}
    pos = 0
    for m in _TOKEN.finditer(formula):
        if m.start() != pos:
            raise ValueError(f"cannot read formula: {formula!r}")
        el, n = m.group(1), m.group(2)
        if el not in ELEMENTS:
            raise ValueError(f"unknown element {el!r} in {formula!r}")
        counts[el] = counts.get(el, 0) + (int(n) if n else 1)
        pos = m.end()
    if pos != len(formula):
        raise ValueError(f"cannot read formula: {formula!r}")
    return counts


def mono_mass(formula):
    """Monoisotopic mass of a neutral formula."""
    return sum(ELEMENTS[el] * n for el, n in parse_formula(formula).items())


def ion_mz(formula, charge):
    """m/z of an ion with this elemental formula and this charge (electron mass applied)."""
    return (mono_mass(formula) - charge * ELECTRON) / abs(charge)


def precursor_mz(neutral_formula, adduct_def):
    """m/z of a precursor from the neutral formula and an adduct definition from the rule file.

    adduct_def example: {"add": "Na", "charge": 1} or {"remove": "H", "charge": -1}
    """
    mass = mono_mass(neutral_formula)
    if adduct_def.get("add"):
        mass += mono_mass(adduct_def["add"])
    if adduct_def.get("remove"):
        mass -= mono_mass(adduct_def["remove"])
    z = adduct_def["charge"]
    return (mass - z * ELECTRON) / abs(z)


def neutral_loss_mz(prec_mz, charge, loss_formula):
    """m/z of a fragment after a neutral loss (the fragment keeps the precursor charge)."""
    return prec_mz - mono_mass(loss_formula) / abs(charge)


def expand_template(template_formula, n, d):
    """Fill a chain template such as 'C{n}H{2n-2d}O2' for a chain with n carbons and d double bonds."""
    def repl(m):
        expr = m.group(1).replace(" ", "")
        expr = re.sub(r"(\d)([nd])", r"\1*\2", expr)
        value = eval(expr, {"__builtins__": {}}, {"n": n, "d": d})
        if value < 0:
            raise ValueError(f"negative atom count in {template_formula} for {n}:{d}")
        return str(value)
    formula = re.sub(r"\{([^}]+)\}", repl, template_formula)
    # remove zero counts such as H0
    return re.sub(r"([A-Z][a-z]?)0(?!\d)", "", formula)
