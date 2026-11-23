"""Check the environment for the MS/MS annotation skill and print a JSON report.

The LLM runs this first, shows the report to the user, and asks before it installs anything.
Usage: python scripts/check_env.py
"""
import importlib
import json
import platform
import shutil
import sys
from pathlib import Path

SKILL_ROOT = Path(__file__).resolve().parent.parent

PACKAGES = {
    "numpy": "numpy",
    "yaml": "pyyaml",
    "openpyxl": "openpyxl",
    "fisher_py": "fisher_py",
    "native_fisher_py": "native-fisher-py",
    "originpro": "originpro",
}


def check_package(module):
    try:
        mod = importlib.import_module(module)
        return {"installed": True, "version": getattr(mod, "__version__", "unknown")}
    except Exception as exc:  # report any import problem, do not stop
        return {"installed": False, "error": f"{type(exc).__name__}: {exc}"}


def find_origin():
    """Return the path of an Origin executable if one is found in the usual folders."""
    names = ["Origin64.exe", "Origin.exe"]
    for name in names:
        hit = shutil.which(name)
        if hit:
            return hit
    for base in [Path("C:/Program Files/OriginLab"), Path("C:/Program Files (x86)/OriginLab")]:
        if base.exists():
            for exe in base.rglob("Origin64.exe"):
                return str(exe)
    return None


def main():
    report = {
        "os": platform.system(),
        "os_version": platform.version(),
        "python": sys.version.split()[0],
        "python_executable": sys.executable,
        "python_ok": sys.version_info >= (3, 9),
        "windows": platform.system() == "Windows",
        "packages": {},
        "missing_pip_names": [],
        "origin_executable": find_origin(),
        "data_files": {
            "rules/positive.yaml": (SKILL_ROOT / "rules" / "positive.yaml").exists(),
            "rules/negative.yaml": (SKILL_ROOT / "rules" / "negative.yaml").exists(),
            "data/precursor_db.csv": (SKILL_ROOT / "data" / "precursor_db.csv").exists(),
            "config.yaml": (SKILL_ROOT / "config.yaml").exists(),
        },
    }
    for module, pip_name in PACKAGES.items():
        res = check_package(module)
        report["packages"][pip_name] = res
        if not res["installed"]:
            report["missing_pip_names"].append(pip_name)
    report["ready_for_annotation"] = (
        report["python_ok"]
        and all(report["packages"][p]["installed"] for p in ["numpy", "pyyaml", "openpyxl"])
        and (report["packages"]["fisher_py"]["installed"] or report["packages"]["native-fisher-py"]["installed"])
        and all(report["data_files"].values())
    )
    report["ready_for_plot"] = (report["windows"] and report["packages"]["originpro"]["installed"]
                                and report["origin_executable"] is not None)
    # raw reader: fisher_py gives the Thermo averager (needs .NET, present on Windows)
    report["raw_reader"] = ("fisher_py" if report["packages"]["fisher_py"]["installed"] else
                            "native-fisher-py (own averaging)" if report["packages"]["native-fisher-py"]["installed"]
                            else None)
    if report["missing_pip_names"]:
        report["install_command"] = f'"{sys.executable}" -m pip install ' + " ".join(report["missing_pip_names"])
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
