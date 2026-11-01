"""Sensitivity indices from model outputs."""

from typing import List

import numpy as np
import pandas as pd
from SALib.analyze import morris as morris_analyze
from SALib.analyze import sobol as sobol_analyze

from .factors import Factor, to_salib_problem


def morris(factors: List[Factor], samples: pd.DataFrame, output: pd.Series,
           levels: int = 4, seed: int = 42) -> pd.DataFrame:
    """Morris mu*, mu and sigma for one scalar output (e.g. flooded area)."""
    problem = to_salib_problem(factors)
    result = morris_analyze.analyze(problem, samples.to_numpy(),
                                    output.loc[samples.index].to_numpy(),
                                    num_levels=levels, seed=seed)
    frame = pd.DataFrame({k: result[k] for k in ("mu", "mu_star", "sigma", "mu_star_conf")},
                         index=problem["names"])
    return frame.sort_values("mu_star", ascending=False)


def sobol(factors: List[Factor], output: pd.Series, second_order: bool = False,
          seed: int = 42) -> pd.DataFrame:
    """First-order (S1) and total-order (ST) Sobol indices for one scalar output."""
    problem = to_salib_problem(factors)
    result = sobol_analyze.analyze(problem, np.asarray(output),
                                   calc_second_order=second_order, seed=seed)
    frame = pd.DataFrame({k: result[k] for k in ("S1", "S1_conf", "ST", "ST_conf")},
                         index=problem["names"])
    return frame.sort_values("ST", ascending=False)
