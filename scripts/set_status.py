"""Record the user's decision for a window and rewrite results.json and results.xlsx.

Usage:
    python scripts/set_status.py <annotation_folder> <window_id> confirmed
    python scripts/set_status.py <annotation_folder> <window_id> changed --annotation "PC 34:1 [M+K]+" --comment "..."
    python scripts/set_status.py <annotation_folder> <window_id> rejected --comment "..."
    python scripts/set_status.py <annotation_folder> all confirmed      (all windows that are still pending)
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from report import write_all  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", type=Path)
    ap.add_argument("window_id")
    ap.add_argument("status", choices=["confirmed", "changed", "rejected", "pending"])
    ap.add_argument("--annotation", default=None)
    ap.add_argument("--comment", default=None)
    args = ap.parse_args()

    path = args.folder / "results.json"
    results = json.load(open(path, encoding="utf-8"))
    targets = [w for w in results["windows"]
               if (args.window_id == "all" and w["status"] == "pending") or w["window_id"] == args.window_id]
    if not targets:
        sys.exit(f"window not found: {args.window_id}")
    if args.status == "changed" and not args.annotation:
        sys.exit("status 'changed' needs --annotation")
    for w in targets:
        w["status"] = args.status
        if args.status == "confirmed":
            w["final_annotation"] = args.annotation or " | ".join(w["match"]["top_annotations"]) or None
        elif args.status == "changed":
            w["final_annotation"] = args.annotation
        elif args.status == "rejected":
            w["final_annotation"] = None
        if args.comment is not None:
            w["comment"] = args.comment
    write_all(results, args.folder)
    print(f"updated {len(targets)} window(s)")


if __name__ == "__main__":
    main()
