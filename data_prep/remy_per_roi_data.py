#!/usr/bin/env python3
"""E119 part 0: Remy's per-ROI megamat recordings (PN boutons, KC claws) reduced to ROI x trial window means, each ROI's
pre-odor noise SD and a lock file, with the checks that they reproduce Remy's own per-fly distance matrices and the
decoded dashboard matrices.

Reads (read-only) Remy's analysis outputs on tensor, <REMY_ROOT>/data_by_acq/:
- kc_claws/volumetric_with_neuropil_correction/<acq>/trials/xrds_bc_trials.nc (4 flies) and
  pn_boutons/volumetric_without_neuropil_correction/<acq>/trials/xrds_bc_trials.nc (6 flies): trials (17 odors x 3) x
  ROIs x time (-5 to 19.8 s, 0.2 s frames), baseline-corrected (median of -5 to 0 s), variables F, Fneu, Fc, spks,
  F_zscore and Fc_zscore;
- <acq>/respvec/xrds_mean_peak.nc: Remy's own per-trial window means (0.5 to 2.5 s), checked against ours;
- <stage dir>/aggregate/xrds_stim_rdm_concat.nc: Remy's per-fly 17 x 17 correlation distances (1 - Pearson r across
  ROIs of the trial-mean window means), check 1;
- <shared>/data/remy_rdm_timeseries/ (decoded from Remy's dashboards): the per-fly per-frame matrices and their window
  means, check 2.

Writes <shared>/data/remy_per_roi/ (README.md there):
- window_means.parquet: one row per (stage, acq, roi, trial), the trial's odor and repeat, and per variable (VARIABLES)
  and window (WINDOWS) the mean over the window's frames (closed windows, as Remy's respvec and the dashboards'
  window_csvs.py select frames);
- noise_sd.parquet: per (stage, acq, roi) and variable, the SD of the frames within the 51 pre-odor windows (-5 to 0 s),
  each window's mean removed (the baseline-noise unit of the realism bands' data side, notes/response_criterion_decision.md);
- frame_rdms.npz: per stage and acq the per-frame trial-mean correlation distances (frames x 17 x 17) of the variable
  Remy's figures use (DIST_VAR), for check 2;
- lock.json: the files read (path, size, mtime, sha256), the flies, ROI counts, frame times, odor names and their
  trial order, the windows, and the checks' results;
- checks.csv, provenance.json.

Part 0b (the same fly, megamat1_redo, from Remy's suite2p output on matrix) has functions here (`cut_trials`,
`suite2p_traces`, `same_fly_checks`) but is NOT run until viki's user clears that source: `main` never reads
/local/matrix.

Usage (from the repository root):
    nice -n 10 ./venv/bin/python -u data_prep/remy_per_roi_data.py [--out DIR]
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

REMY_ROOT = Path('/local/tensor/Remy/odor_space_collab/analysis_outputs')
STAGE_DIRS = {'bouton': 'pn_boutons/volumetric_without_neuropil_correction',
              'claw': 'kc_claws/volumetric_with_neuropil_correction'}
VARIABLES = ('F', 'Fc', 'F_zscore', 'Fc_zscore', 'spks')
# the variable Remy's figures and the dashboards use per stage
DIST_VAR = {'bouton': 'F_zscore', 'claw': 'Fc_zscore'}
# the windows (s after odor onset), closed: Remy's respvec window, E90's STAGE_WINDOW and the dashboards' windows
WINDOWS = {'remy': (0.5, 2.5), 'e90': (0.5, 3.0), 'early': (0.5, 1.0), 'mid': (1.0, 2.0), 'late': (2.0, 3.0),
           'pre': (-5.0, -0.2)}
PRE_NOISE = (-5.0, 0.0)                 # the pre-odor frames whose SD (each trial's mean removed) is the noise unit
EXACT_TOL = 1e-6                        # check 1: Remy's stim RDMs and window means reproduced to this
DASH_R_MIN = 0.98                       # check 2: per fly, r over the 136 pairs of the dashboard's window means


def out_dir() -> Path:
    raise NotImplementedError


def odor_of_stim(stim: str) -> str:
    """'ms @ -3.0' -> 'ms'."""
    raise NotImplementedError


def acquisitions() -> pd.DataFrame:
    """One row per stage and acquisition: stage, acq (the folder name), path of xrds_bc_trials.nc."""
    raise NotImplementedError


def frame_mask(t: np.ndarray, window: Tuple[float, float]) -> np.ndarray:
    raise NotImplementedError


def window_means(x: np.ndarray, t: np.ndarray, windows: Dict[str, Tuple[float, float]] = WINDOWS) -> Dict[str, np.ndarray]:
    """{window: trials x ROIs} means of x (trials x ROIs x frames) over each window's frames."""
    raise NotImplementedError


def noise_sd(x: np.ndarray, t: np.ndarray, window: Tuple[float, float] = PRE_NOISE) -> np.ndarray:
    """Per ROI the SD of the pre-odor frames of all trials, each trial's pre-odor mean removed (ROIs)."""
    raise NotImplementedError


def trial_repeats(stims: Sequence[str]) -> np.ndarray:
    """The repeat (0, 1, 2) of every trial: its odor's occurrence number in the trial order."""
    raise NotImplementedError


def stim_rdm(wm: np.ndarray, odors_of_trials: Sequence[str], odors: Sequence[str]) -> np.ndarray:
    """1 - Pearson r across ROIs of the trial-mean window means (trials x ROIs) of every two odors (odors x odors)."""
    raise NotImplementedError


def frame_rdms(x: np.ndarray, odors_of_trials: Sequence[str], odors: Sequence[str]) -> np.ndarray:
    """Per frame, 1 - Pearson r across ROIs of the trial-mean responses (frames x odors x odors)."""
    raise NotImplementedError


def reduce_acq(stage: str, acq: str, path: Path) -> dict:
    """One acquisition's window means, noise SDs, frame RDMs and its lock entry."""
    raise NotImplementedError


def check_remy_rdms(stage: str, rdms: Dict[str, np.ndarray], odors: Sequence[str]) -> pd.DataFrame:
    """Check 1: our stim RDMs (Remy's window, each variable) against her aggregate xrds_stim_rdm_concat.nc."""
    raise NotImplementedError


def check_dashboard(stage: str, frames: Dict[str, np.ndarray], t: np.ndarray, odors: Sequence[str]) -> pd.DataFrame:
    """Check 2: our per-frame RDMs' window means against the decoded dashboard's, flies matched one to one."""
    raise NotImplementedError


# part 0b: functions only (the source on matrix is not read until cleared)

def suite2p_traces(folder: Path) -> dict:
    """F, Fneu, iscell of a suite2p 'combined' folder (F.npy, Fneu.npy, iscell.npy)."""
    raise NotImplementedError


def cut_trials(traces: np.ndarray, frame_times: np.ndarray, onsets: Sequence[float],
               t: np.ndarray = np.round(np.arange(-5.0, 19.81, 0.2), 6)) -> np.ndarray:
    """Cut traces (ROIs x frames) into trials (trials x ROIs x len(t)) around the odor onsets, by linear interpolation
    onto t, baseline-corrected by the median of t in [-5, 0]."""
    raise NotImplementedError


def same_fly_checks(claw_frames: np.ndarray, soma_frames: np.ndarray, t: np.ndarray, odors: Sequence[str]) -> dict:
    """Part 0b's gate: the decoded same-fly series reproduced (r >= 0.98) and the claws-against-somata agreement (1-2 s)
    recomputed (Remy's 0.86)."""
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
