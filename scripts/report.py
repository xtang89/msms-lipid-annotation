"""Write results.json and results.xlsx from the same data (they always agree)."""
import json
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter


def _frag_text(frags):
    return "; ".join(f"{f['name']} ({f['obs_mz']:.4f}, {f['ppm']:+.1f} ppm, {100 * f['rel_int']:.1f}%)" for f in frags)


def write_json(results, out_dir):
    path = Path(out_dir) / "results.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1)
    return path


def write_xlsx(results, out_dir):
    wb = Workbook()
    ws = wb.active
    ws.title = "results"
    head = ["window_id", "raw_file", "polarity", "precursor_mz", "precursor_source", "activation", "collision_energy",
            "n_scans", "level", "msms_confirmed", "top_annotations", "diagnostic_fragments", "supporting_fragments",
            "chain_fragments", "chain_combinations", "other_candidates", "unexplained_peaks",
            "status", "final_annotation", "comment", "filter"]
    ws.append(head)
    fr = wb.create_sheet("fragments")
    fr.append(["window_id", "candidate", "level", "fragment", "role", "type", "theo_mz", "obs_mz", "ppm", "rel_int_%"])
    for w in results["windows"]:
        cands = w["match"]["candidates"]
        top = cands[0] if cands else None
        level = w["match"]["best_level"]
        same = [c for c in cands if c["level"] == level] if cands else []
        others = [c for c in cands if c not in same]
        ws.append([
            w["window_id"], w["raw_file"], w["polarity"], round(w["precursor_mz"], 5), w["precursor_source"],
            w["activation"], w["collision_energy"], w["n_scans"], level,
            "yes" if top and top["msms_confirmed"] else "no",
            " | ".join(w["match"]["top_annotations"]),
            " | ".join(_frag_text([f for f in c["fragments"] if f["role"] == "diagnostic"]) for c in same),
            " | ".join(_frag_text([f for f in c["fragments"] if f["role"] == "supporting"]) for c in same),
            " | ".join(_frag_text(c["chain_fragments"]) for c in same),
            " | ".join(", ".join(c["complete_chain_combinations"]) for c in same),
            "; ".join(f"{c['annotation']} [{c['level']}]" for c in others),
            "; ".join(f"{p['mz']:.4f} ({100 * p['rel_int']:.1f}%)" for p in w["match"]["unexplained_peaks"]),
            w.get("status", "pending"), w.get("final_annotation") or "", w.get("comment") or "", w["filter"],
        ])
        for c in cands:
            for f in c["fragments"] + c["chain_fragments"]:
                fr.append([w["window_id"], c["annotation"], c["level"], f["name"], f["role"], f["type"],
                           f["theo_mz"], f["obs_mz"], f["ppm"], round(100 * f["rel_int"], 2)])
    st = wb.create_sheet("settings")
    for k, v in results["run_info"].items():
        st.append([k, json.dumps(v) if isinstance(v, (dict, list)) else v])
    for sheet in wb.worksheets:
        for cell in sheet[1]:
            cell.font = Font(bold=True)
        for i, col in enumerate(sheet.columns, 1):
            width = max(len(str(c.value)) if c.value is not None else 0 for c in col)
            sheet.column_dimensions[get_column_letter(i)].width = min(max(10, width + 2), 60)
        sheet.freeze_panes = "B2"
    path = Path(out_dir) / "results.xlsx"
    wb.save(path)
    return path


def write_all(results, out_dir):
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    return write_json(results, out_dir), write_xlsx(results, out_dir)
