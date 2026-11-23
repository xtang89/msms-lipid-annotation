"""Match one averaged MS/MS spectrum against the precursor database and the rule files.

Decision tree (per window)
1. Polarity from the scan filter -> adduct list from config.yaml.
2. Candidates = precursor_db rows with that polarity and adduct within tolerance of the precursor m/z.
   No candidate -> "no match".
3. For each candidate, find the class fragments (diagnostic / supporting) and the chain fragments.
4. Level: high (>=2 diagnostic, or 1 diagnostic + >=1 supporting), medium (1 diagnostic),
   low (supporting only), ms1_only (no fragment). Classes with class_rule "chains" (TG, DG, MG):
   high if >=2 chain fragments fit one chain combination (1 is enough if all its chains are the same).
   high and medium = MS/MS confirmed.
5. Chain combination reported only if fragments for all its chains are present.
"""
import csv
import itertools
from pathlib import Path

import numpy as np
import yaml

from lipidmass import ion_mz, mono_mass, neutral_loss_mz, expand_template

SKILL_ROOT = Path(__file__).resolve().parent.parent
LEVEL_ORDER = {"high": 0, "medium": 1, "low": 2, "ms1_only": 3}
CONFIRMED = {"high", "medium"}


# ---------------------------------------------------------------- loading
def load_rules(rules_dir=SKILL_ROOT / "rules"):
    """Return {polarity: rule dict}, and an index {(class, adduct): entry}."""
    rules, index = {}, {}
    for fn in ("positive.yaml", "negative.yaml"):
        with open(Path(rules_dir) / fn, encoding="utf-8") as f:
            r = yaml.safe_load(f)
        pol = r["meta"]["polarity"]
        rules[pol] = r
        for e in r["entries"]:
            index[(e["class"], e["adduct"])] = e
        for fx in r.get("fixed_precursors", []):
            index[(fx["name"], fx["adduct"])] = dict(fx, **{"class": fx["name"], "n_chains": 0,
                                                            "chain_fragments": []})
    return rules, index


def load_precursor_db(path=SKILL_ROOT / "data" / "precursor_db.csv"):
    rows = []
    with open(path, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            r["precursor_mz"] = float(r["precursor_mz"])
            r["charge"] = int(r["charge"])
            rows.append(r)
    rows.sort(key=lambda r: r["precursor_mz"])
    return rows


# ---------------------------------------------------------------- helpers
def tol_da(mz, cfg):
    return max(mz * cfg["tolerance_ppm"] * 1e-6, cfg["min_tolerance_da"])


def ppm_err(obs, theo):
    return (obs - theo) / theo * 1e6


class Peaks:
    """Centroided peaks above the intensity threshold (relative to the base peak, full range)."""

    def __init__(self, mz, intensity, threshold):
        mz, it = np.asarray(mz, float), np.asarray(intensity, float)
        base = float(it.max()) if it.size else 0.0
        rel = it / base if base > 0 else it
        keep = rel >= threshold
        order = np.argsort(mz[keep])
        self.mz, self.rel = mz[keep][order], rel[keep][order]

    def find(self, target, tol):
        """Most intense peak within target +- tol. Return (index, mz, rel) or None."""
        if self.mz.size == 0:
            return None
        lo, hi = np.searchsorted(self.mz, target - tol), np.searchsorted(self.mz, target + tol, side="right")
        if hi <= lo:
            return None
        k = int(lo) + int(np.argmax(self.rel[lo:hi]))
        return k, float(self.mz[k]), float(self.rel[k])


def fragment_mz(frag, prec_theo, charge):
    """Theoretical m/z of a class fragment."""
    if frag["type"] == "ion":
        z = frag.get("charge", 1 if charge > 0 else -1)
        return ion_mz(frag["formula"], z)
    return neutral_loss_mz(prec_theo, charge, frag["formula"])


def chain_list(rule):
    c0, c1 = rule["chain_range"]["carbons"]
    d0, d1 = rule["chain_range"]["double_bonds"]
    return [(c, d) for c in range(c0, c1 + 1) for d in range(d0, d1 + 1) if d <= c // 2]


def chain_combinations(n_chains, sum_c, sum_d, chains):
    """All multisets of n_chains chains whose carbons and double bonds add to the sum."""
    out = []
    for combo in itertools.combinations_with_replacement(sorted(chains), n_chains):
        if sum(c for c, _ in combo) == sum_c and sum(d for _, d in combo) == sum_d:
            out.append(combo)
    return out


def chain_name(ch):
    return f"{ch[0]}:{ch[1]}"


# ---------------------------------------------------------------- candidate scoring
def score_candidate(cand, entry, rule, peaks, cfg):
    charge = cand["charge"]
    prec_theo = cand["precursor_mz"]
    matched = []

    for frag in entry.get("fragments", []):
        theo = fragment_mz(frag, prec_theo, charge)
        hit = peaks.find(theo, tol_da(theo, cfg))
        if hit:
            k, obs, rel = hit
            matched.append({"name": frag["name"], "role": frag["role"], "type": frag["type"],
                            "formula": frag["formula"],
                            "theo_mz": round(theo, 5), "obs_mz": round(obs, 5),
                            "ppm": round(ppm_err(obs, theo), 2), "rel_int": round(rel, 4), "peak": k})

    # chains
    chain_hits = {}          # chain -> list of matched chain fragments
    complete, chain_rule_met = [], False
    n_chains = entry.get("n_chains", 0)
    can_chain = (n_chains > 0 and entry.get("chain_fragments") and cand.get("ether") != "yes"
                 and cand.get("sum_carbons"))
    if can_chain:
        combos = chain_combinations(n_chains, int(cand["sum_carbons"]), int(cand["sum_db"]), chain_list(rule))
        needed = sorted({ch for combo in combos for ch in combo})
        for ch in needed:
            for cf in entry["chain_fragments"]:
                t = rule["chain_templates"][cf["template"]]
                formula = expand_template(t["formula"], *ch)
                name = t["name"].replace("{n}", str(ch[0])).replace("{d}", str(ch[1]))
                if t["type"] == "ion":
                    theo = ion_mz(formula, t.get("charge", 1 if charge > 0 else -1))
                else:
                    loss = mono_mass(formula) + (mono_mass(cf["plus"]) if cf.get("plus") else 0.0)
                    theo = prec_theo - loss / abs(charge)
                    if cf.get("plus"):
                        name += f" + {cf['plus']}"
                hit = peaks.find(theo, tol_da(theo, cfg))
                if hit:
                    k, obs, rel = hit
                    chain_hits.setdefault(ch, []).append(
                        {"name": name, "role": "chain", "type": t["type"], "chain": chain_name(ch),
                         "formula": formula + (f"+{cf['plus']}" if cf.get("plus") and t["type"] == "nl" else ""),
                         "theo_mz": round(theo, 5), "obs_mz": round(obs, 5),
                         "ppm": round(ppm_err(obs, theo), 2), "rel_int": round(rel, 4), "peak": k})
        for combo in combos:
            distinct = set(combo)
            if all(ch in chain_hits for ch in distinct):
                complete.append("_".join(chain_name(ch) for ch in combo))
            # class rule "chains": fragments (distinct peaks) from chains of this one combination
            n_frag = len({f["peak"] for ch in distinct for f in chain_hits.get(ch, [])})
            if n_frag >= 2 or (len(distinct) == 1 and n_frag >= 1):
                chain_rule_met = True

    chain_frags = [f for lst in chain_hits.values() for f in lst]
    n_diag = sum(1 for f in matched if f["role"] == "diagnostic")
    n_supp = sum(1 for f in matched if f["role"] == "supporting")

    # level
    if entry.get("class_rule") == "chains":
        if chain_rule_met:
            level = "high"
        elif chain_frags:
            level = "low"
        else:
            level = "ms1_only"
    else:
        if n_diag >= 2 or (n_diag >= 1 and n_supp >= 1):
            level = "high"
        elif n_diag == 1:
            level = "medium"
        elif n_supp >= 1:
            level = "low"
        else:
            level = "ms1_only"

    species = f"{cand['abbrev']} {cand['adduct']}"
    if level in CONFIRMED and complete:
        annotation = f"{species} ({', '.join(complete)})"
    else:
        annotation = species
    explained_peaks = {f["peak"] for f in matched + chain_frags}
    return {
        "class": cand["class"], "abbrev": cand["abbrev"], "adduct": cand["adduct"],
        "formula": cand["formula"], "theo_mz": round(prec_theo, 5),
        "ether": cand.get("ether") == "yes",
        "level": level, "msms_confirmed": level in CONFIRMED,
        "n_diagnostic": n_diag, "n_supporting": n_supp, "n_chain_fragments": len(chain_frags),
        "fragments": matched, "chain_fragments": chain_frags,
        "chains_observed": sorted({f["chain"] for f in chain_frags}),
        "complete_chain_combinations": complete,
        "annotation": annotation,
        "explained_intensity": round(float(sum(peaks.rel[k] for k in explained_peaks)), 4),
        "_peaks": explained_peaks,
    }


def match_window(window, rules, index, db, cfg):
    pol = window["polarity"]
    peaks = Peaks(window["mz"], window["intensity"], cfg["intensity_threshold"])
    prec = window["precursor_mz"]
    ptol = (cfg["filter_precursor_tolerance_da"] if window["precursor_source"] == "filter"
            else tol_da(prec, cfg))
    allowed = set(cfg["adducts"][pol])

    cands = [r for r in db if r["polarity"] == pol and r["adduct"] in allowed
             and abs(r["precursor_mz"] - prec) <= ptol]
    results = []
    for c in cands:
        entry = index.get((c["class"], c["adduct"]))
        if entry is None:
            res = {
                "class": c["class"], "abbrev": c["abbrev"], "adduct": c["adduct"], "formula": c["formula"],
                "theo_mz": round(c["precursor_mz"], 5), "ether": c.get("ether") == "yes",
                "level": "ms1_only", "msms_confirmed": False, "rules_available": False,
                "n_diagnostic": 0, "n_supporting": 0, "n_chain_fragments": 0,
                "fragments": [], "chain_fragments": [], "chains_observed": [],
                "complete_chain_combinations": [],
                "annotation": f"{c['abbrev']} {c['adduct']} (accurate mass only; no rules)",
                "explained_intensity": 0.0, "_peaks": set(),
            }
        else:
            res = score_candidate(c, entry, rules[pol], peaks, cfg)
            res["rules_available"] = True
        res["precursor_ppm"] = round(ppm_err(prec, c["precursor_mz"]), 2)
        results.append(res)

    results.sort(key=lambda r: (LEVEL_ORDER[r["level"]], not r["rules_available"], -r["n_diagnostic"], -r["n_supporting"],
                                -r["n_chain_fragments"], -r["explained_intensity"], abs(r["precursor_ppm"])))

    # unexplained peaks (not matched by any candidate, not the precursor)
    explained = set().union(*[r["_peaks"] for r in results]) if results else set()
    unexplained = []
    for k in np.argsort(-peaks.rel):
        k = int(k)
        if k in explained or abs(peaks.mz[k] - prec) <= tol_da(prec, cfg):
            continue
        unexplained.append({"mz": round(float(peaks.mz[k]), 5), "rel_int": round(float(peaks.rel[k]), 4)})
        if len(unexplained) >= cfg["n_unexplained"]:
            break
    for r in results:
        r.pop("_peaks")

    best_level = results[0]["level"] if results else "no_match"
    top = [r for r in results if r["level"] == best_level] if results else []
    return {
        "candidates": results,
        "best_level": best_level,
        "top_annotations": [r["annotation"] for r in top],
        "unexplained_peaks": unexplained,
        "n_peaks_above_threshold": int(peaks.mz.size),
    }
