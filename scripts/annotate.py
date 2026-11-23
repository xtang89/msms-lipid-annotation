"""Annotate the MS/MS windows of one or more Thermo .raw files.

For each .raw file, the results go into <rawfile>_annotation/ next to the .raw file:
    results.json   (for the LLM)    results.xlsx   (for people)    plots/ (made later by plot_origin.py)

Usage:
    python scripts/annotate.py file1.raw [file2.raw ...]
    python scripts/annotate.py folder_with_raw_files
"""
import argparse
import datetime
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from average import read_windows  # noqa: E402
from match import load_rules, load_precursor_db, match_window  # noqa: E402
from report import write_all  # noqa: E402

SKILL_ROOT = Path(__file__).resolve().parent.parent


def annotate_file(raw_path, cfg, rules, index, db):
    windows = read_windows(raw_path, cfg)
    out = []
    for w in windows:
        res = match_window(w, rules, index, db, cfg)
        inten = np.asarray(w["intensity"], float)
        base = float(inten.max()) if inten.size else 1.0
        keep = inten / base >= cfg["save_floor"]
        w_out = {k: v for k, v in w.items() if k not in ("mz", "intensity")}
        w_out["spectrum"] = {
            "mz": [round(x, 5) for x in np.asarray(w["mz"])[keep].tolist()],
            "rel_int": [round(x, 5) for x in (inten[keep] / base).tolist()],
        }
        w_out["match"] = res
        w_out["status"] = "pending"
        w_out["final_annotation"] = None
        w_out["comment"] = ""
        out.append(w_out)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", type=Path)
    ap.add_argument("--config", type=Path, default=SKILL_ROOT / "config.yaml")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config, encoding="utf-8"))
    rules, index = load_rules()
    db = load_precursor_db()

    raws = []
    for p in args.inputs:
        raws += sorted(q for q in p.iterdir() if q.suffix.lower() == ".raw") if p.is_dir() else [p]
    for raw in raws:
        windows = annotate_file(raw, cfg, rules, index, db)
        results = {
            "run_info": {
                "raw_file": str(raw),
                "date": datetime.datetime.now().isoformat(timespec="seconds"),
                "rules_version": {pol: r["meta"]["version"] for pol, r in rules.items()},
                "precursor_db": "data/precursor_db.csv (LIPID MAPS COMP_DB)",
                "settings": cfg,
            },
            "windows": windows,
        }
        out_dir = raw.parent / f"{raw.stem}_annotation"
        jpath, xpath = write_all(results, out_dir)
        n_conf = sum(1 for w in windows if w["match"]["best_level"] in ("high", "medium"))
        print(f"{raw.name}: {len(windows)} windows, {n_conf} MS/MS confirmed -> {jpath.parent}")


if __name__ == "__main__":
    main()
