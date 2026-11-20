"""Full confusion matrices of model runs against several Sentinel-1 flood maps.

Reads the per-cell state at the radar pass saved by structural.py or
time_series.py (radar_state.npz) and scores it against each radar map given, so
the scores can be checked for the radar classification threshold and for the
cell-fraction threshold without rerunning the model.

Usage:
    python -m research.periyar.radar_check OUT_CSV --states NAME=NPZ ... --sar NAME=TIF ...
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .validate import OPEN_CLASSES, _sample


def confusion(model_wet, observed_wet):
    hits = int((model_wet & observed_wet).sum())
    misses = int((~model_wet & observed_wet).sum())
    false_alarms = int((model_wet & ~observed_wet).sum())
    negatives = int((~model_wet & ~observed_wet).sum())
    n = hits + misses + false_alarms + negatives
    expected = ((hits + misses) * (hits + false_alarms)) / n if n else 0.0
    return {
        "hits": hits, "misses": misses, "false_alarms": false_alarms, "correct_negatives": negatives,
        "hit_rate": hits / max(hits + misses, 1),
        "false_alarm_ratio": false_alarms / max(hits + false_alarms, 1),
        "csi": hits / max(hits + misses + false_alarms, 1),
        "bias": (hits + false_alarms) / max(hits + misses, 1),
        # Equitable threat score: CSI corrected for hits expected by chance.
        "ets": (hits - expected) / max(hits + misses + false_alarms - expected, 1e-9),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out_csv", type=Path)
    parser.add_argument("--states", nargs="+", required=True, help="NAME=radar_state.npz")
    parser.add_argument("--sar", nargs="+", required=True, help="NAME=flood_fraction.tif")
    parser.add_argument("--fractions", type=float, nargs="+", default=[0.5, 0.3])
    args = parser.parse_args()

    rows = []
    for state_arg in args.states:
        run, path = state_arg.split("=", 1)
        s = np.load(path)
        for sar_arg in args.sar:
            sar_name, tif = sar_arg.split("=", 1)
            fraction = _sample(Path(tif), s["x"], s["y"], -1.0)
            land = s["land"] & (fraction >= 0)
            open_land = land & np.isin(s["landcover"], OPEN_CLASSES)
            for cut in args.fractions:
                observed = fraction > cut
                for domain, mask in (("all_land", land), ("open_land", open_land)):
                    rows.append({"run": run, "sar": sar_name, "fraction_cut": cut, "domain": domain,
                                 **confusion(s["model_wet"][mask], observed[mask])})
    table = pd.DataFrame(rows)
    table.to_csv(args.out_csv, index=False)
    print(table.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
