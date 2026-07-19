"""Per-machine runtime calibration for heavy jobs. Records (job, input size, elapsed) samples
on this PC and estimates the next run's duration, so the UI can show an accurate ETA.
Calibration lives in .teagle/timings.json and is refined as more runs complete."""
from __future__ import annotations
import os, json
from . import appdirs

_FILE = os.path.join(appdirs.user_data_dir(), "timings.json")
_MAXKEEP = 40
# fallback until this PC has samples: (base_seconds, seconds_per_kb)
_HEURISTIC = {"annotate": (25.0, 3.0), "splice": (3.0, 0.05), "install": (300.0, 0.0)}


def _load() -> dict:
    try:
        with open(_FILE, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def _save(d: dict):
    try:
        os.makedirs(os.path.dirname(_FILE), exist_ok=True)
        with open(_FILE, "w", encoding="utf-8") as f:
            json.dump(d, f)
    except Exception:
        pass


def record(job: str, size: int, secs: float):
    d = _load()
    arr = d.get(job, [])
    arr.append({"size": int(size), "secs": round(float(secs), 2)})
    d[job] = arr[-_MAXKEEP:]
    _save(d)


def estimate(job: str, size: int) -> dict:
    """Predict this PC's runtime for `job` on an input of `size` bp.
    Least-squares line from samples when available, else a heuristic."""
    arr = _load().get(job, [])
    size = int(size)
    base, per_kb = _HEURISTIC.get(job, (10.0, 1.0))
    heuristic = base + per_kb * (size / 1000.0)
    if len(arr) >= 2:
        n = len(arr)
        sx = sum(p["size"] for p in arr); sy = sum(p["secs"] for p in arr)
        sxx = sum(p["size"] ** 2 for p in arr); sxy = sum(p["size"] * p["secs"] for p in arr)
        denom = n * sxx - sx * sx
        sizes = [p["size"] for p in arr]
        if denom > 0:
            b = (n * sxy - sx * sy) / denom
            a = (sy - b * sx) / n
            est = a + b * size
            # only trust the line within/near the sampled size range with a positive slope;
            # a flat/negative slope or an extrapolation far past the samples reverts to the heuristic
            if b > 0 and size <= 1.5 * max(sizes) and est >= 1.0:
                return {"eta_s": round(est), "basis": "calibrated", "n": n}
            mean_rate = (sy / n) / max(sum(sizes) / n, 1)      # calibrated per-bp rate for extrapolation
            return {"eta_s": round(max(1.0, mean_rate * size)), "basis": "calibrated", "n": n}
        return {"eta_s": round(max(1.0, sy / n)), "basis": "calibrated", "n": n}   # all samples same size
    if arr:
        p = arr[-1]
        rate = p["secs"] / max(p["size"], 1)
        return {"eta_s": round(max(1.0, rate * size)), "basis": "calibrated", "n": 1}
    return {"eta_s": round(heuristic), "basis": "heuristic", "n": 0}
