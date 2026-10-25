"""Sampling designs for screening (Morris), variance-based (Sobol) GSA and Monte Carlo."""

from typing import List

import numpy as np
import pandas as pd
from SALib.sample import latin
from SALib.sample import morris as morris_sample
from SALib.sample import sobol as sobol_sample

from .factors import Factor, to_salib_problem


def morris(factors: List[Factor], trajectories: int = 20, levels: int = 4,
           seed: int = 42) -> pd.DataFrame:
    """Morris elementary-effects design with (k + 1) * trajectories runs."""
    problem = to_salib_problem(factors)
    x = morris_sample.sample(problem, trajectories, num_levels=levels, seed=seed)
    return _frame(problem, x)


def sobol(factors: List[Factor], base_samples: int = 256,
          second_order: bool = False, seed: int = 42) -> pd.DataFrame:
    """Saltelli/Sobol design with N * (k + 2) runs (N * (2k + 2) with second order)."""
    problem = to_salib_problem(factors)
    x = sobol_sample.sample(problem, base_samples,
                            calc_second_order=second_order, seed=seed)
    return _frame(problem, x)


def latin_hypercube(factors: List[Factor], samples: int = 100, seed: int = 42) -> pd.DataFrame:
    """Latin hypercube design with ``samples`` runs, e.g. for GLUE-style conditioning."""
    problem = to_salib_problem(factors)
    return _frame(problem, latin.sample(problem, samples, seed=seed))


def _frame(problem, x: np.ndarray) -> pd.DataFrame:
    frame = pd.DataFrame(x, columns=problem["names"])
    frame.index.name = "run"
    return frame
