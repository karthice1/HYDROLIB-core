from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hydrolib.core.dflowfm.bc.models import ForcingModel
from hydrolib.core.dflowfm.friction.models import FrictionModel
from hydrolib.core.dflowfm.structure.models import StructureModel
from research.gsa import analysis, ensemble, sampling
from research.gsa.factors import PERIYAR_2018, Factor, to_salib_problem

BASE = (Path(__file__).parents[3] / "tests/data/input/e02/f152_1d2d_projectmodels_rhu"
        / "c04_DHydamo-MGB-initialisation/fm")

FACTORS = [
    Factor("n_main", 0.5, 2.0, "friction_multiplier", ("Strickler_70.0",)),
    Factor("q_upstream", 0.5, 1.5, "discharge_multiplier", ("140703_391887",)),
    Factor("h_downstream", -0.5, 0.5, "waterlevel_offset", ("141134_395446",)),
]
FILES = {"friction": ["roughness/r001.ini"], "forcing": ["boundaries.bc"]}


def test_periyar_factors_are_valid():
    problem = to_salib_problem(PERIYAR_2018)
    assert problem["num_vars"] == len(PERIYAR_2018)


def test_invalid_factor_rejected():
    with pytest.raises(ValueError):
        to_salib_problem([Factor("a", 1.0, 0.5, "friction_multiplier")])


def test_sample_sizes():
    assert len(sampling.morris(FACTORS, trajectories=5)) == 5 * (len(FACTORS) + 1)
    assert len(sampling.sobol(FACTORS, base_samples=8)) == 8 * (len(FACTORS) + 2)


def test_ensemble_applies_perturbations(tmp_path):
    samples = pd.DataFrame([[2.0, 1.5, 0.25]], columns=[f.name for f in FACTORS])
    manifest = ensemble.build(BASE, tmp_path, samples, FACTORS, FILES)
    run = Path(manifest.loc[0, "path"])

    # Strickler is a conveyance: twice as rough means half the coefficient.
    assert FrictionModel(run / "roughness/r001.ini").global_[0].frictionvalue == pytest.approx(35.0)

    base = {f.name: f for f in ForcingModel(BASE / "boundaries.bc").forcing}
    new = {f.name: f for f in ForcingModel(run / "boundaries.bc").forcing}
    q0 = np.array(base["140703_391887"].datablock)[:, 1]
    q1 = np.array(new["140703_391887"].datablock)[:, 1]
    assert np.allclose(q1, 1.5 * q0)
    h0 = np.array(base["141134_395446"].datablock)[:, 1]
    h1 = np.array(new["141134_395446"].datablock)[:, 1]
    assert np.allclose(h1, h0 + 0.25)
    # Forcings that are not targeted stay unchanged.
    for name in base:
        if name not in ("140703_391887", "141134_395446"):
            assert np.allclose(np.array(new[name].datablock), np.array(base[name].datablock))
    assert (manifest.filter(like="edits_") > 0).all(axis=None)


def test_crest_offset(tmp_path):
    model = StructureModel(BASE / "structure.ini")
    from research.gsa import perturb
    before = {s.id: s.crestlevel for s in model.structure
              if isinstance(getattr(s, "crestlevel", None), (int, float))}
    edits = perturb.offset_crest(model, 0.1)
    assert edits == len(before) > 0
    for s in model.structure:
        if s.id in before:
            assert s.crestlevel == pytest.approx(before[s.id] + 0.1)


def test_morris_recovers_known_ranking():
    # Synthetic output where q_upstream dominates and h_downstream is inert.
    samples = sampling.morris(FACTORS, trajectories=20)
    y = 10 * samples["q_upstream"] + 1 * samples["n_main"]
    result = analysis.morris(FACTORS, samples, y)
    assert list(result.index[:2]) == ["q_upstream", "n_main"]
    assert result.loc["h_downstream", "mu_star"] == pytest.approx(0.0)


def test_sobol_recovers_known_ranking():
    samples = sampling.sobol(FACTORS, base_samples=256)
    y = 10 * samples["q_upstream"] + 1 * samples["n_main"]
    result = analysis.sobol(FACTORS, y)
    assert result.index[0] == "q_upstream"
    assert result.loc["q_upstream", "ST"] > 0.9
