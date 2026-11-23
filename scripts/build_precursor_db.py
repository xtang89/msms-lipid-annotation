"""Build the precursor database from LIPID MAPS COMP_DB and the rule files.

For each class + adduct in the rule files, the script takes every COMP_DB species of that
head group, calculates the precursor m/z from the formula, and writes one row per species
and adduct. Prostaglandins (fixed_precursors in negative.yaml) are added from their formulas.

Usage (from the skill folder):
    python scripts/build_precursor_db.py
    python scripts/build_precursor_db.py --comp-db COMP_DB.tsv --rules rules --out data/precursor_db.csv
"""
import argparse
import csv
import re
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lipidmass import mono_mass, precursor_mz  # noqa: E402

SKILL_ROOT = Path(__file__).resolve().parent.parent

# "PC 34:1", "PC O-34:1", "PE P-36:4", "SM 34:1;O2", "TG 52:2;O"
# Note: COMP_DB puts lyso species (LPC, LPE, ...) under the same head group as the diacyl class,
# so the class is taken from the name, not from the head-group column.
ABBREV = re.compile(r"^(?P<cls>\S+)\s+(?P<ether>[OP]-)?(?P<c>\d+):(?P<d>\d+)(?P<mod>;\S+)?$")

# Oxygen suffix that is normal for a class (not an oxidation). SM ;O2 = dihydroxy sphingoid base.
NORMAL_MOD = {"SM": ";O2"}

COLUMNS = ["class", "abbrev", "formula", "neutral_mass", "adduct", "polarity", "charge",
           "precursor_mz", "sum_carbons", "sum_db", "ether", "rules", "source"]


def load_rules(rules_dir):
    rules = []
    for fn in ("positive.yaml", "negative.yaml"):
        with open(rules_dir / fn, encoding="utf-8") as f:
            rules.append(yaml.safe_load(f))
    return rules


def load_comp_db(path):
    with open(path, encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def build(comp_db_path, rules_dir, out_path, config_path, include_oxidized=False, include_ether=True, mass_check=0.001):
    rules = load_rules(rules_dir)
    comp = load_comp_db(comp_db_path)
    by_hg = {}
    for row in comp:
        by_hg.setdefault(row["headgroup"], []).append(row)

    out_rows, warnings = [], []
    n_impossible = 0
    for rule in rules:
        pol = rule["meta"]["polarity"]
        adducts = rule["meta"]["adduct_formulas"]
        cmin, cmax = rule["chain_range"]["carbons"]
        dmax = rule["chain_range"]["double_bonds"][1]
        for entry in rule["entries"]:
            hg, adduct = entry["comp_db_headgroup"], entry["adduct"]
            species = by_hg.get(hg, [])
            if not species:
                warnings.append(f"no COMP_DB species for head group {hg}")
            for row in species:
                m = ABBREV.match(row["abbrev"].strip())
                if not m:
                    warnings.append(f"skipped (cannot read name): {row['abbrev']}")
                    continue
                if m.group("cls") != entry["class"]:
                    continue  # for example LPC rows under head group PC
                mod = m.group("mod")
                if mod and mod != NORMAL_MOD.get(entry["class"]) and not include_oxidized:
                    continue
                if m.group("ether") and not include_ether:
                    continue
                n_ch = entry.get("n_chains", 0)
                if n_ch > 0:
                    c_sum, d_sum = int(m.group("c")), int(m.group("d"))
                    if not (n_ch * cmin <= c_sum <= n_ch * cmax and d_sum <= n_ch * dmax):
                        n_impossible += 1  # cannot be made from n_ch chains in the chain range
                        continue
                formula = row["formula"].strip()
                neutral = mono_mass(formula)
                if abs(neutral - float(row["mass"])) > mass_check:
                    warnings.append(f"mass differs from COMP_DB for {row['abbrev']}: "
                                    f"{neutral:.4f} vs {row['mass']}")
                out_rows.append({
                    "class": entry["class"], "abbrev": row["abbrev"].strip(), "formula": formula,
                    "neutral_mass": f"{neutral:.6f}", "adduct": adduct, "polarity": pol,
                    "charge": adducts[adduct]["charge"],
                    "precursor_mz": f"{precursor_mz(formula, adducts[adduct]):.6f}",
                    "sum_carbons": m.group("c"), "sum_db": m.group("d"),
                    "ether": "yes" if m.group("ether") else "no",
                    "rules": "yes", "source": "COMP_DB",
                })
        for fx in rule.get("fixed_precursors", []):
            adduct = fx["adduct"]
            out_rows.append({
                "class": fx["name"], "abbrev": fx["name"], "formula": fx["formula"],
                "neutral_mass": f"{mono_mass(fx['formula']):.6f}", "adduct": adduct, "polarity": pol,
                "charge": adducts[adduct]["charge"],
                "precursor_mz": f"{precursor_mz(fx['formula'], adducts[adduct]):.6f}",
                "sum_carbons": "", "sum_db": "", "ether": "no", "rules": "yes", "source": "fixed",
            })

    # species of classes (or class + adduct) without rules: accurate-mass annotation only
    with open(config_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    ruled = {(e["class"], e["adduct"]) for r in rules for e in r["entries"]}
    adduct_defs = {r["meta"]["polarity"]: r["meta"]["adduct_formulas"] for r in rules}
    class_chains = {e["class"]: e.get("n_chains", 0) for r in rules for e in r["entries"]}
    crange = rules[0]["chain_range"]
    n_ms1 = 0
    for row in comp:
        m = ABBREV.match(row["abbrev"].strip())
        if not m:
            continue
        cls, mod = m.group("cls"), m.group("mod")
        if mod and mod != NORMAL_MOD.get(cls) and not include_oxidized:
            continue
        if m.group("ether") and not include_ether:
            continue
        n_ch = class_chains.get(cls, 0)
        if n_ch > 0:  # class known from the rules: use the same chain-range filter
            c_sum, d_sum = int(m.group("c")), int(m.group("d"))
            if not (n_ch * crange["carbons"][0] <= c_sum <= n_ch * crange["carbons"][1]
                    and d_sum <= n_ch * crange["double_bonds"][1]):
                continue
        formula = row["formula"].strip()
        for pol, adducts in cfg["adducts"].items():
            for adduct in adducts:
                if (cls, adduct) in ruled:
                    continue
                ad = adduct_defs[pol][adduct]
                out_rows.append({
                    "class": cls, "abbrev": row["abbrev"].strip(), "formula": formula,
                    "neutral_mass": f"{mono_mass(formula):.6f}", "adduct": adduct, "polarity": pol,
                    "charge": ad["charge"], "precursor_mz": f"{precursor_mz(formula, ad):.6f}",
                    "sum_carbons": m.group("c"), "sum_db": m.group("d"),
                    "ether": "yes" if m.group("ether") else "no", "rules": "no", "source": "COMP_DB",
                })
                n_ms1 += 1

    out_rows.sort(key=lambda r: (r["polarity"], float(r["precursor_mz"])))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(out_rows)

    counts = {}
    for r in [r for r in out_rows if r["rules"] == "yes"]:
        counts[(r["class"], r["adduct"])] = counts.get((r["class"], r["adduct"]), 0) + 1
    print(f"wrote {len(out_rows)} rows to {out_path}")
    print(f"excluded {n_impossible} species outside the chain range of the rule files")
    print(f"added {n_ms1} rows without rules (accurate-mass annotation only)")
    for (c, a), n in sorted(counts.items()):
        print(f"  {c:20s} {a:14s} {n}")
    if warnings:
        print(f"{len(warnings)} warnings (first 20):")
        for wmsg in warnings[:20]:
            print("  " + wmsg)
    return out_rows, warnings


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--comp-db", type=Path, default=SKILL_ROOT / "COMP_DB.tsv")
    ap.add_argument("--rules", type=Path, default=SKILL_ROOT / "rules")
    ap.add_argument("--out", type=Path, default=SKILL_ROOT / "data" / "precursor_db.csv")
    ap.add_argument("--config", type=Path, default=SKILL_ROOT / "config.yaml")
    ap.add_argument("--include-oxidized", action="store_true",
                    help="also include oxidized species (names with ;O, ;O2 ...); SM ;O2 is always included")
    ap.add_argument("--no-ether", action="store_true", help="exclude ether species (O-, P-)")
    args = ap.parse_args()
    build(args.comp_db, args.rules, args.out, args.config, include_oxidized=args.include_oxidized,
          include_ether=not args.no_ether)


if __name__ == "__main__":
    main()
