"""Audit of the inflow perturbation: is any water counted twice?

Writes the perturbed Neeleswaram inflow of the base model for factor combinations
at the ends of their ranges, under both release models of sensitivity.py, and
reports total inflow volume, peak discharge and the share of the lagged dam
releases in the gauged hydrograph. The Chalakudy inflow is not perturbed.

Usage:
    python -m research.periyar.forcing_audit BASE_DIR OUT_DIR
"""

import argparse
import itertools
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from hydrolib.core.dflowfm.bc.models import ForcingModel

from .sensitivity import IDAMALAYAR_LAG_HOURS, IDUKKI_LAG_HOURS, NOMINAL, _release_hours, \
    _rewrite_forcing


def inflows(bc_file: Path) -> dict:
    model = ForcingModel(bc_file)
    out = {}
    for forcing in model.forcing:
        if forcing.name.startswith(("upstream", "chalakudy")):
            block = np.asarray(forcing.datablock, float)
            out[forcing.name] = (block[:, 0], block[:, 1])
    return out


def volume_mcm(t_hours, q):
    return float(np.trapezoid(q, t_hours * 3600.0) / 1e6)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base_dir", type=Path)
    parser.add_argument("out_dir", type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    base = inflows(args.base_dir / "boundaries.bc")
    t, q = base["upstream_0001"]
    (t_idk, q_idk), (t_ida, q_ida) = _release_hours((IDUKKI_LAG_HOURS, IDAMALAYAR_LAG_HOURS))
    r_idk, r_ida = np.interp(t, t_idk, q_idk), np.interp(t, t_ida, q_ida)
    share = pd.DataFrame({"hours": t, "neeleswaram_m3s": q, "idukki_lagged_m3s": r_idk,
                          "idamalayar_lagged_m3s": r_ida})
    share.to_csv(args.out_dir / "release_share.csv", index=False)
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for model, (qr, mi, ma) in itertools.product(
                ["additive", "multiplicative"], itertools.product([0.85, 1.0, 1.15], [0.8, 1.0, 1.2], [0.8, 1.0, 1.2])):
            run_dir = Path(tmp) / "run"
            shutil.rmtree(run_dir, ignore_errors=True)
            run_dir.mkdir()
            shutil.copy(args.base_dir / "boundaries.bc", run_dir / "boundaries.bc")
            row = pd.Series({**NOMINAL, "q_rating": qr, "q_idukki": mi, "q_idamalayar": ma})
            _rewrite_forcing(run_dir, row, model)
            new = inflows(run_dir / "boundaries.bc")
            tq, qq = new["upstream_0001"]
            rows.append({"release_model": model, "q_rating": qr, "q_idukki": mi, "q_idamalayar": ma,
                         "upstream_volume_mcm": volume_mcm(tq, qq), "upstream_peak_m3s": qq.max(),
                         "chalakudy_volume_mcm": volume_mcm(*new["chalakudy_0001"]),
                         "chalakudy_unchanged": bool(np.allclose(new["chalakudy_0001"][1], base["chalakudy_0001"][1]))})
    table = pd.DataFrame(rows)
    nominal = table.query("release_model == 'additive' and q_rating == 1 and q_idukki == 1 and q_idamalayar == 1").iloc[0]
    table["volume_vs_nominal"] = table["upstream_volume_mcm"] / nominal["upstream_volume_mcm"]
    table.to_csv(args.out_dir / "forcing_audit.csv", index=False)
    summary = {
        "nominal_upstream_volume_mcm": nominal["upstream_volume_mcm"],
        "nominal_peak_m3s": float(q.max()),
        "idukki_lagged_volume_mcm": volume_mcm(t, r_idk),
        "idamalayar_lagged_volume_mcm": volume_mcm(t, r_ida),
        "release_share_of_volume": (volume_mcm(t, r_idk) + volume_mcm(t, r_ida)) / volume_mcm(t, q),
        "release_share_at_peak": float((r_idk + r_ida)[q.argmax()] / q.max()),
        "max_difference_between_models_m3s": float(np.max(np.abs(
            0.15 * 0.2 * (r_idk + r_ida)))),
    }
    pd.Series(summary).to_csv(args.out_dir / "summary.csv", header=["value"])
    print(pd.Series(summary).round(3).to_string())
    print(table.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
