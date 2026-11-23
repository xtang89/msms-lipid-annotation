"""Read a Thermo .raw file and make one averaged, centroided MS/MS spectrum per precursor window.

Steps
1. Group the MS2 scans by their scan filter string (one filter = one precursor window).
2. Average all scans of each group (Thermo average, or own average of centroided scans).
3. Centroid the averaged spectrum if it is profile data.
4. Find the precursor m/z: MS1 peak -> precursor peak in the MS2 spectrum -> filter value (2 decimals).

This module is used by annotate.py. It can also run alone to list the windows of a file:
    python scripts/average.py path/to/file.raw
"""
import re
import sys
from pathlib import Path

import numpy as np

# Reader: fisher_py (pythonnet + Thermo .NET libraries; needs .NET, present on Windows) gives the Thermo
# scan averager. native-fisher-py (no .NET needed) reads the scans, but its averager does not work.
try:
    from fisher_py import RawFile as _FisherRaw
    READER = "fisher_py"
except Exception:  # no fisher_py, or no .NET runtime
    _FisherRaw = None
    from native_fisher_py import RawFile as _NativeRaw
    READER = "native-fisher-py"


class Raw:
    """Small wrapper with the same methods for both readers."""

    def __init__(self, path):
        self._r = _FisherRaw(path) if _FisherRaw else _NativeRaw(path)
        self.first_scan, self.last_scan = self._r.first_scan, self._r.last_scan

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        acc = getattr(self._r, "_raw_file_access", None)
        if acc is not None:
            acc.dispose()
        elif hasattr(self._r, "__exit__"):
            self._r.__exit__(*exc)

    def filter_string(self, s):
        if _FisherRaw:
            return self._r._raw_file_access.get_scan_event_string_for_scan_number(s)
        return str(self._r.get_filter_for_scan_number(s))

    def scan(self, s):
        return self._r.get_scan_from_scan_number(s)

    def rt(self, s):
        return self._r.get_retention_time_from_scan_number(s)

    def thermo_average(self, scans):
        """Thermo scan averager (fisher_py only). Returns (mz, intensity)."""
        if not _FisherRaw:
            raise RuntimeError("the Thermo averager needs fisher_py")
        a = self._r._raw_file_access.average_scans(list(scans), None, None, True)
        return np.array(a.preferred_masses, float), np.array(a.preferred_intensities, float)

# Example filters:
#   "FTMS + p NSI Full ms2 760.58@hcd30.00 [50.0000-800.0000]"
#   "FTMS - c ESI d Full ms2 834.53@cid35.00 [225.0000-850.0000]"
FILTER_MS2 = re.compile(r"ms2\s+(?P<prec>\d+\.\d+)@(?P<act>[a-z]+)(?P<ce>\d+(\.\d+)?)", re.I)
FILTER_POL = re.compile(r"\s([+-])\s")
FILTER_RANGE = re.compile(r"\[(?P<lo>\d+(\.\d+)?)-(?P<hi>\d+(\.\d+)?)\]")
FILTER_MS1 = re.compile(r"Full ms(?!\d)", re.I)


def parse_filter(filt):
    """Return a dict with ms_level, polarity, precursor (filter value), activation, collision energy."""
    pol = FILTER_POL.search(" " + filt + " ")
    info = {"filter": filt, "polarity": None, "ms_level": None}
    if pol:
        info["polarity"] = "positive" if pol.group(1) == "+" else "negative"
    m = FILTER_MS2.search(filt)
    if m:
        info.update(ms_level=2, precursor_filter=float(m.group("prec")),
                    activation=m.group("act").upper(), collision_energy=float(m.group("ce")))
    elif FILTER_MS1.search(filt):
        info["ms_level"] = 1
    r = FILTER_RANGE.search(filt)
    if r:
        info["scan_range"] = [float(r.group("lo")), float(r.group("hi"))]
    return info


def tol_da(mz, ppm, min_da):
    return max(mz * ppm * 1e-6, min_da)


def centroid_profile(mz, inten):
    """Simple centroiding of profile data: local maxima, intensity-weighted m/z over each peak."""
    mz, inten = np.asarray(mz, float), np.asarray(inten, float)
    if mz.size < 3:
        return mz, inten
    peaks_mz, peaks_i = [], []
    i, n = 0, mz.size
    while i < n:
        if inten[i] <= 0:
            i += 1
            continue
        j = i
        while j + 1 < n and inten[j + 1] > 0:
            j += 1
        seg_mz, seg_i = mz[i:j + 1], inten[i:j + 1]
        # split the segment at local minima
        cuts = [0] + [k for k in range(1, len(seg_i) - 1)
                      if seg_i[k] < seg_i[k - 1] and seg_i[k] <= seg_i[k + 1]] + [len(seg_i) - 1]
        for a, b in zip(cuts[:-1], cuts[1:]):
            sm, si = seg_mz[a:b + 1], seg_i[a:b + 1]
            top = int(np.argmax(si))
            lo, hi = max(top - 1, 0), min(top + 2, len(si))
            w = si[lo:hi]
            peaks_mz.append(float(np.sum(sm[lo:hi] * w) / np.sum(w)))
            peaks_i.append(float(si[top]))
        i = j + 1
    return np.array(peaks_mz), np.array(peaks_i)


def looks_like_profile(mz):
    """Profile data has many closely spaced points."""
    if len(mz) < 50:
        return False
    return float(np.median(np.diff(np.sort(mz)))) < 0.01


def own_average(raw, scans, ppm, min_da):
    """Average centroided scans: cluster peaks within tolerance, mean m/z (intensity weighted), sum/n intensity."""
    all_mz, all_i = [], []
    for s in scans:
        m, i, _, _ = raw.scan(s)
        all_mz.append(np.asarray(m, float))
        all_i.append(np.asarray(i, float))
    mz = np.concatenate(all_mz) if all_mz else np.array([])
    it = np.concatenate(all_i) if all_i else np.array([])
    if mz.size == 0:
        return mz, it
    order = np.argsort(mz)
    mz, it = mz[order], it[order]
    out_mz, out_i = [], []
    start = 0
    for k in range(1, mz.size + 1):
        if k == mz.size or mz[k] - mz[k - 1] > tol_da(mz[k - 1], ppm, min_da):
            w = it[start:k]
            out_mz.append(float(np.sum(mz[start:k] * w) / np.sum(w)) if np.sum(w) > 0 else float(mz[start]))
            out_i.append(float(np.sum(w) / len(scans)))
            start = k
    return np.array(out_mz), np.array(out_i)


def find_peak(mz, inten, target, window):
    """Most intense peak within target +- window. Return m/z or None."""
    if len(mz) == 0:
        return None
    mz, inten = np.asarray(mz), np.asarray(inten)
    sel = np.abs(mz - target) <= window
    if not np.any(sel):
        return None
    idx = np.where(sel)[0][np.argmax(inten[sel])]
    return float(mz[idx])


def read_windows(raw_path, cfg):
    """Return a list of averaged MS/MS windows for one .raw file."""
    ppm, min_da = cfg["tolerance_ppm"], cfg["min_tolerance_da"]
    search = cfg["precursor_search_da"]
    windows = []
    with Raw(str(raw_path)) as raw:
        ms2, ms1_scans, rts = [], [], {}
        for s in range(raw.first_scan, raw.last_scan + 1):
            filt = raw.filter_string(s)
            info = parse_filter(filt)
            if info["ms_level"] == 2:
                ms2.append((info["polarity"], info.get("activation"), info.get("collision_energy"),
                            info["precursor_filter"], s, filt))
                rts[s] = raw.rt(s)
            elif info["ms_level"] == 1:
                ms1_scans.append(s)

        # one window = same polarity, activation and collision energy, precursor within tolerance
        groups = []
        ms2.sort(key=lambda t: (t[0] or "", t[1] or "", t[2] or 0, t[3]))
        for item in ms2:
            last = groups[-1] if groups else None
            if (last and last["key"] == item[:3]
                    and item[3] - last["precs"][-1] <= tol_da(item[3], ppm, min_da)):
                last["precs"].append(item[3]); last["scans"].append(item[4]); last["filters"].add(item[5])
            else:
                groups.append({"key": item[:3], "precs": [item[3]], "scans": [item[4]], "filters": {item[5]}})

        ms1_mz, ms1_i = (own_average(raw, ms1_scans, ppm, min_da) if ms1_scans
                         else (np.array([]), np.array([])))

        thermo_ok = True
        groups.sort(key=lambda g: (g["key"][0] or "", float(np.median(g["precs"]))))
        for n, g in enumerate(groups, 1):
            scans = sorted(g["scans"])
            filt = sorted(g["filters"])[0]
            info = parse_filter(filt)
            info["precursor_filter"] = float(np.median(g["precs"]))
            method = cfg.get("average_method", "thermo")
            if method == "thermo" and thermo_ok:
                try:
                    mz, inten = raw.thermo_average(scans)
                    if looks_like_profile(mz):
                        mz, inten = centroid_profile(mz, inten)
                except RuntimeError:
                    thermo_ok = False  # native-fisher-py: use own average
            if method != "thermo" or not thermo_ok:
                method = "own"
                mz, inten = own_average(raw, scans, ppm, min_da)

            target = info["precursor_filter"]
            prec = find_peak(ms1_mz, ms1_i, target, search)
            source = "MS1"
            if prec is None:
                prec = find_peak(mz, inten, target, search)
                source = "MS2"
            if prec is None:
                prec, source = target, "filter"

            keep = inten > 0
            windows.append({
                "window_id": f"{Path(raw_path).stem}_W{n:04d}_{target:.4f}",
                "raw_file": Path(raw_path).name,
                "filter": filt if len(g["filters"]) == 1 else f"{filt} (+{len(g['filters']) - 1} similar filters)",
                "polarity": info["polarity"],
                "precursor_filter": target,
                "precursor_mz": prec,
                "precursor_source": source,
                "activation": info.get("activation"),
                "collision_energy": info.get("collision_energy"),
                "scan_range": info.get("scan_range"),
                "n_scans": len(scans),
                "scans": [scans[0], scans[-1]],
                "rt_min": [round(min(rts[s] for s in scans), 4), round(max(rts[s] for s in scans), 4)],
                "average_method": method,
                "mz": mz[keep].tolist(),
                "intensity": inten[keep].tolist(),
            })
    return windows


if __name__ == "__main__":
    import yaml
    root = Path(__file__).resolve().parent.parent
    cfg = yaml.safe_load(open(root / "config.yaml", encoding="utf-8"))
    for w in read_windows(sys.argv[1], cfg):
        print(f"{w['window_id']}: {w['filter']} | scans {w['n_scans']} | precursor {w['precursor_mz']:.4f} "
              f"({w['precursor_source']}) | {len(w['mz'])} peaks")
