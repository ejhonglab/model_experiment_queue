#!/usr/bin/env python3
"""E121 step 0: the onset race, read from the family's spikes and APL's trace (no new model code).

The family is the breadth-2 APL-effect fit on its eight wirings (E72): S = the hemibrain wiring (C2, E67's
rate0.068_aplboth_k2_cap3), P1, P2, M1; T = P3, P4, M2, M3. For each wiring the fitted model is rebuilt from its
rebuild archive (pyolfsysm.rebuild, data/model_outputs/<exp>/rebuild/<tag>/) at the default 0.5 s pulse, APL on and
APL off (acute), with olfsysm's spike recordings and APL trace (`inh_sims`) saved: one `fit_mb_model` call per
condition with `return_olfsysm_vars=True`, then the KC layer run again with `mp.kc.save_spike_recordings` and
`mp.kc.save_inh_sims` set (committed olfsysm, no option added). The rerun must reproduce the archive's saved spike
counts and responses exactly, else the run stops.

Per odor and wiring: APL's rise, the time from odor onset to 90 % of its plateau (`apl_rise`; the plateau is the
mean of the trace over the last 100 ms of the pulse; a second reading takes the trace's maximum over the pulse).
Per responding KC and odor (2 or more spikes) the time of its second spike, the response's completion under the
2-spike criterion. A response "completes before APL's rise" when its second spike is at or before APL's 90 % point
of that odor and wiring (APL on's trace in both conditions: APL off is the reference, read against the same clock).

The gate (E121's entry): fewer than 20 % of APL-on responses complete before APL's 90 % point: no race, and the
experiment closes without code.

Groups (E77_input_screen's): ester group 2h, IaA, pa, eb, ep, ms; alcohol group 1-8ol, 1-5ol, 1-6ol, 6al; the rest
(2-but, t2h, paa, -aPine, euc, B-myr) "other". Breadth: the number of the 16 odors the KC responds to in that
condition (bins 1, 2-3, 4-7, 8+).

Outputs (results/E121_onset_race/): step0_shares.csv (share complete before APL's 90 % point, overall, by group,
by breadth, by odor, per wiring and pooled, APL on and off, both plateau readings), step0_apl_rise.csv (per wiring
and odor: plateau, time to 50 % and 90 %, both readings), step0_gate.json, provenance.json. Per-wiring spike times
and APL traces under outputs/E121_onset_race/step0/.

Usage (from the repository root):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E121_onset_race.py step0
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E121_onset_race'
GATE_SHARE = 0.20
RISE_FRACTION = 0.9


def apl_rise(trace: np.ndarray, time_s: np.ndarray, stim_end: float = 0.5, plateau: str = 'late',
             fraction: float = RISE_FRACTION) -> Tuple[float, float]:
    """Returns (plateau, time from onset to `fraction` of it) of one APL trace."""
    raise NotImplementedError


def second_spike_times(spikes: pd.DataFrame) -> pd.DataFrame:
    """Returns the second spike's time per (kc, odor) with 2 or more spikes."""
    raise NotImplementedError


def shares(completion: pd.DataFrame, rise: pd.DataFrame) -> pd.DataFrame:
    """Returns the shares of responses complete before APL's rise, overall and by group, breadth, odor."""
    raise NotImplementedError


def run_wiring(label: str, archive: Path, out: Path) -> Dict[str, pd.DataFrame]:
    """Rebuilds one wiring's fit, APL on and off, and returns its spikes and APL traces."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('part', choices=['step0'])
    ap.parse_args()
    raise NotImplementedError


if __name__ == '__main__':
    main()
