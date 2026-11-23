"""Plot annotated MS/MS spectra in Origin (Windows only). Export OPJU and TIF per window.

Design (rules/../SKILL.md, "Plot design")
- Sticks: x = m/z, y = relative intensity (%), 0 to y_max_percent.
- Colors (Okabe-Ito): unannotated peaks gray; class fragments (diagnostic + supporting) blue;
  chain fragments vermillion; precursor marked black.
- Labels only on matched fragments: line 1 = observed m/z, line 2 = assignment.
- Header: annotation, precursor m/z, activation and collision energy.

Usage:
    python scripts/plot_origin.py <annotation_folder>            (windows with status confirmed or changed)
    python scripts/plot_origin.py <annotation_folder> --all      (all windows)
"""
import argparse
import json
import math
import sys
from pathlib import Path

import originpro as op

SKILL_ROOT = Path(__file__).resolve().parent.parent
COLORS = {"other": "#999999", "class": "#0072B2", "chain": "#D55E00", "precursor": "#000000"}
PLOT_DEFAULTS = {"width_cm": 8.5, "height_cm": 6.0, "font": "Arial", "font_size_pt": 9,
                 "label_decimals": 4, "y_max_percent": 110, "tif_dpi": 600}


def load_plot_settings():
    """Plot settings from config.yaml; Origin's own Python may not have PyYAML, then the defaults are used."""
    try:
        import yaml
    except ImportError:
        return dict(PLOT_DEFAULTS)
    with open(SKILL_ROOT / "config.yaml", encoding="utf-8") as f:
        return {**PLOT_DEFAULTS, **yaml.safe_load(f).get("plot", {})}


def sticks(mz, rel):
    x, y = [], []
    for m, r in zip(mz, rel):
        x += [m, m, math.nan]
        y += [0.0, 100.0 * r, math.nan]
    return x, y


def short_label(f, prec_theo, polarity):
    """Assignment text: '-183.066' for a neutral loss, 'C5H15NO4P+' for an ion, 'FA 18:1' for a chain ion."""
    sign = "+" if polarity == "positive" else "-"
    if f["role"] == "chain":
        if f["type"] == "nl":
            return "-FA " + f["chain"] + (" (ketene)" if "ketene" in f["name"] else "")
        return "FA " + f["chain"]
    if f["type"] == "nl":
        return f"-{prec_theo - f['theo_mz']:.3f}"
    return f["formula"] + sign


def place_labels(items, lo, hi, ymax, line_h):
    """items: (x, y_peak, line1, line2). Return (x_center, y_low_line, line1, line2).
    Labels are centred over the peak; a label that overlaps a placed label moves up."""
    span = hi - lo
    placed = []
    for x, y, l1, l2 in sorted(items, key=lambda t: -t[1]):
        half_w = 0.0105 * span * max(len(l1), len(l2)) / 2
        xc = min(max(x, lo + half_w), hi - half_w)       # keep the text inside the frame
        yb = y + 0.6 * line_h
        for _ in range(15):
            clash = [p for p in placed
                     if abs(p[0] - xc) < (p[4] + half_w) and abs(p[1] - yb) < 2.1 * line_h]
            if not clash:
                break
            yb = max(p[1] for p in clash) + 2.1 * line_h
        placed.append((xc, yb, l1, l2, half_w))
    return [(xc, yb, l1, l2) for xc, yb, l1, l2, _ in placed]


def pick_candidate(w):
    cands = w["match"]["candidates"]
    if not cands:
        return None
    final = w.get("final_annotation") or ""
    for c in cands:
        if final and (c["annotation"] == final or f"{c['abbrev']} {c['adduct']}" in final):
            return c
    return cands[0]


def plot_window(w, out_dir, pcfg):
    cand = pick_candidate(w)
    mz, rel = w["spectrum"]["mz"], w["spectrum"]["rel_int"]
    frags = (cand["fragments"] + cand["chain_fragments"]) if cand else []
    class_mz = {f["obs_mz"] for f in frags if f["role"] != "chain"}
    chain_mz = {f["obs_mz"] for f in frags if f["role"] == "chain"}

    groups = {"other": ([], []), "class": ([], []), "chain": ([], [])}
    for m, r in zip(mz, rel):
        key = "class" if round(m, 5) in class_mz else "chain" if round(m, 5) in chain_mz else "other"
        groups[key][0].append(m)
        groups[key][1].append(r)

    op.new(asksave=False)
    wks = op.new_sheet("w", lname=w["window_id"][:30])
    col = 0
    datasets = []
    for key in ("other", "class", "chain"):
        x, y = sticks(*groups[key])
        if not x:
            continue
        wks.from_list(col, x, lname="m/z", axis="X")
        wks.from_list(col + 1, y, lname=f"{key} (%)", axis="Y")
        datasets.append((key, col, col + 1))
        col += 2

    gp = op.new_graph(template="line", lname=w["window_id"][:30])
    gl = gp[0]
    for key, cx, cy in datasets:
        p = gl.add_plot(wks, coly=cy, colx=cx, type="l")
        p.color = COLORS[key]
    lo, hi = (w.get("scan_range") or [min(mz) - 10, max(mz) + 10])
    gl.set_xlim(lo, hi)
    gl.set_ylim(0, pcfg["y_max_percent"])
    gl.axis("x").title = "m/z"
    gl.axis("y").title = "Relative intensity (%)"

    op.lt_exec("label -r legend;")  # the design has no legend

    gp.activate()
    fs = pcfg["font_size_pt"]
    op.lt_exec(f"layer.x.label.pt={fs}; layer.y.label.pt={fs}; xb.fsize={fs}; yl.fsize={fs};")
    d = pcfg["label_decimals"]
    ymax = pcfg["y_max_percent"]
    line_h = 0.045 * ymax
    rel_of = {round(m, 5): rr for m, rr in zip(mz, rel)}
    prec_theo = cand["theo_mz"] if cand else w["precursor_mz"]
    items = [(f["obs_mz"], 100.0 * rel_of.get(f["obs_mz"], f["rel_int"]),
              f"{f['obs_mz']:.{d}f}", short_label(f, prec_theo, w["polarity"])) for f in frags]
    n = 0

    def text_at(xc, yc, text, size):
        nonlocal n
        n += 1
        name = f"pk{n}"
        op.lt_exec(f'label -n {name} "{text}"; {name}.fsize={size}; '
                   f'{name}.x={xc:.5f}; {name}.y={yc:.3f};')

    lab_fs = max(fs - 2, 6)
    for xc, yb, l1, l2 in place_labels(items, lo, hi, ymax, line_h):
        text_at(xc, yb, l2, lab_fs)               # assignment (lower line)
        text_at(xc, yb + line_h, l1, lab_fs)      # m/z (upper line)
    # precursor marker above the precursor peak
    prec = w["precursor_mz"]
    p_rel = max([100.0 * rr for m, rr in zip(mz, rel) if abs(m - prec) < 0.01] or [0.0])
    text_at(prec, min(p_rel + 0.6 * line_h, ymax - line_h / 2), "\\x(25BC)", lab_fs)
    # header above the frame
    title = w.get("final_annotation") or (cand["annotation"] if cand else "no match")
    header = f"{title}, m/z {prec:.{d}f}, {w.get('activation') or ''} {w.get('collision_energy') or ''}".strip()
    op.lt_exec(f'label -n hdr "{header}"; hdr.fsize={fs}; hdr.y={ymax * 1.06:.3f}; '
               f'hdr.x={lo:.5f} + hdr.dx / 2;')

    plots = Path(out_dir).resolve() / "plots"
    plots.mkdir(exist_ok=True)
    width_px = int(pcfg["width_cm"] / 2.54 * pcfg["tif_dpi"])
    gp.save_fig(str(plots / f"{w['window_id']}.tif"), type="tif", width=width_px, replace=True)
    op.save(str(plots / f"{w['window_id']}.opju"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", type=Path)
    ap.add_argument("--all", action="store_true", help="all windows, not only confirmed/changed")
    ap.add_argument("--windows", default=None, help="comma-separated window_ids to plot")
    ap.add_argument("--confirmed-level", action="store_true",
                    help="plot windows with level high or medium (before the user decision)")
    args = ap.parse_args()
    args.folder = args.folder.resolve()  # Origin is a separate program: give it absolute paths
    pcfg = load_plot_settings()
    results = json.load(open(args.folder / "results.json", encoding="utf-8"))
    if args.windows:
        ids = set(args.windows.split(","))
        windows = [w for w in results["windows"] if w["window_id"] in ids]
    elif args.confirmed_level:
        windows = [w for w in results["windows"] if w["match"]["best_level"] in ("high", "medium")]
    else:
        windows = [w for w in results["windows"] if args.all or w["status"] in ("confirmed", "changed")]
    if op.oext:
        op.set_show(False)
    try:
        for w in windows:
            plot_window(w, args.folder, pcfg)
            print(f"plotted {w['window_id']}")
    finally:
        if op.oext:
            op.exit()


if __name__ == "__main__":
    main()
