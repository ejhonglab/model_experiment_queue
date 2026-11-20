#!/usr/bin/env python3
"""E114: a compressive map from the ORN recordings to the model's input rates. Every odor's positive input to an exponent
below 1 (or through a saturating map), on the family's eight wirings at fixed scales with both response rates re-held,
then the family fitted again at the exponent the gate prefers.

The entry is item E114 of experiment_priorities.md (proposed by viki 2026-11-19 13:58 at Tom's request of 13:49 relayed
by hal; planned by a Fable subagent, claude-fable-5-2, spawned from viki 13:50; claimed by viki 2026-11-19 19:04 on hal's
yes). It comes from handle S3 of the input side of the error common to all wirings (E77; outcome (a3)). "The family" is
C2's (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each of E72 part B's eight wirings at its own fit
(E83's `saved_fit`; S: real, P1, P2, M1; T: P3, P4, M2, M3), read as the mean matrix over the eight wirings (E72's
population, 'all8'), and over S and T beside. "The two ρ's" are ρ on odor-odor distances with APL on against the
control flies and with APL off against APL output blocked (TNT), on raw counts (cosine distance on spike counts), each
KC matrix the mean over all recordings or over one of E108's halves A to D (E108's `dist_targets`); the APL-effect
statistics are E108's `unit_scores` (raw counts and the calibrated readout, against all recordings and the halves).

Concentrations. paa's input is the lab's recording at 10^-4 for the flies' 10^-3; euc and -aPine stay at the recorded
10^-1.75 and 10^-1.5 for the flies' 10^-2; B-myr and 2-but are at the flies' concentrations. Every output involves paa,
euc and -aPine (the map acts on every odor). The holds stand: the model's hemibrain wiring and E72's re-draws of it only;
no other connectome table; no raw-only recording. Nothing is published to claude.ai.

The inputs (glomeruli x the 16 odors, the lab's 36 rows; inputs/<name>.csv with sha256.json; fixed before any run):
- pow<e>: E77's S3 rule (`E77_input_screen.orn_power`): each odor's positive deltas to the power e, rescaled so that
  their sum over glomeruli is the odor's own sum of positive deltas; zero and negative deltas unchanged. Part A:
  e = 0.6, 0.7, 0.8, 0.9 (EXPONENTS); the wrong-way control e = 1.15 (WRONG_WAY); pow1 is the identity (part 0).
- sat<c>: the saturating map x / (1 + x / c) on each positive delta (x in the model's input units, estimated spike
  deltas in Hz), c = 100, 200, 400 Hz (SAT_C), rescaled by the same rule (the odor's sum of positive deltas kept); zero
  and negative deltas unchanged. Part B.
- I0: the lab's input, `H.apl16_orn_deltas()` (the reference: the saved fits).

Parts (CLI `--part`; every simulation part is resumable, a point done when its row is in the part's table):
- inputs (no simulations): the inputs' files, sha256.json and inputs_summary.csv (per input and odor: the largest
  delta, the sum of positive deltas, the participation ratio over glomeruli).
- 0 (simulations): the identity. pow1 equals the lab's input bit for bit; on each wiring (`--wirings` restricts), the
  simulator built with pow1 and its breadth recomputed (E85's `build`, this experiment's cache) gives the wiring's
  breadth and reproduces the saved APL-on and APL-off runs to the last count and response, and the held reading of pow1
  (E85's `held_reading`) reproduces them with both factors exactly 1 (identity.json; 'all_exact' needs all eight).
  If anything differs the experiment stops there.
- 0c (one E10_joint_fit.py smoke run, 5 evaluations): the refit's path. E72's part-B command on the real wiring with
  `--orn-input results/E114_input_compression/inputs/pow1.csv` and `--init-from <the fit> --init-from-exact`: its first
  evaluation equals the fit's best history row and its fits.csv row on every shared column (E110's `exact_check`;
  check.json).
- A (simulations): each of pow0.6 to pow0.9 and pow1.15 on the eight wirings at fixed scales with both response rates
  re-held (the breadth recomputed for the input, the somatic threshold bisected to the fit's APL-off response count,
  then the tied APL scale to its APL-on count; E85's `held_reading`). Rows in runs_A.csv (absolute parameters on and
  off, the factors, realism under both response criteria, `H.code_commits()`); runs under outputs/.../A/.
- B (simulations): sat100, sat200, sat400 on the eight wirings, the same way (runs_B.csv; outputs/.../B/).
- score (no simulations): readings.csv (per input, unit and readout, E108's `unit_scores` and the runs' means; 'd_'
  columns the change from I0's same unit and readout), runs.csv, gate.csv and gate.json (`gate` per input of parts A,
  B and C), choice.json (part C's exponent, `choose_c`), outcome.json (`outcome_of`), beside tables when their inputs
  exist (e112.csv, pn_stage.csv, bouton_stage.csv, d1.csv, d2.csv), comparison_table.md, realism_table.md.
- Cstart (simulations only when the chosen exponent has no part-A runs, i.e. 0.75): the held readings of the chosen
  exponent on the eight wirings (runs_Cs.csv; outputs/.../Cs/), then starts.json (C_<wiring>: the held point's APL-on
  parameters, absolute).
- C (fits, one E10_joint_fit.py process each): the family fitted again with the chosen map on each of the eight
  wirings: E72's part-B command (`-k 2 --aplpn-max 3 --rate-target ctrl_frac --loss apl --apl-readout both --tag
  _cap3_C --maxfev 150`, the wiring's options as E110's `wiring_args`) with `--orn-input
  results/E114_input_compression/inputs/pow<e>.csv --init-params results/E114_input_compression/starts.json:C_<w>`,
  into results/ and outputs/ of E114_input_compression/fits/C_<w>/ (launch.json beside; logs under outputs/.../logs/).
- D (simulations; only when part C's class is "better", `--force` overrides): D1, the chosen map on the validation2
  ORN input (E107's `panel_input('validation2')`, the same rule per odor), T's four wirings at their refits'
  parameters with the breadth of the mapped megamat input on that wiring (the refit's own), nothing re-held, APL on
  and acutely off; read against Remy's validation2 somata as E109's `remy_readings` (E107's carrying set of pairs and
  the 31 of paa and B-myr, her copy under E107's lock), with E107's saved runs of the family as fitted beside
  (E109's `e107_reference`). D2, the map on top of E113's best input when E113 has one (`--e113-input <csv>`, glomeruli
  x abbreviated odors as E10's `--orn-input` reads it; without it D2 records that E113 has none): E113's input and the
  map on it, each on the eight wirings held as in part A.
- beside (simulations, no class): the PN stage of every input (E90's `al_stage` at the bouton base's parameters,
  `pn_statistics` against the PN dendrites, E90's al_transfer.csv); the bouton stage on the real wiring, APL on, at each
  input's held point (I0 at the fit, part C's input at its refit): E90's `read_stages` over its MODEL_WINDOW, the
  boutons of the 28 glomeruli the flies share, 1 - Pearson across boutons per megamat pair, group means over E90's
  GROUPS for the 'rate' and 'session_z' readings, beside the flies' boutons as measured and with the noise removed
  (E90's stages_flies.csv). These need no KC data.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive and kcq/apl_effect.py.

The gate (the entry's, fixed before any run; `gate`; each input's population of eight against I0's, each wiring paired
with itself; E108's constants):
- Gain: one ρ higher by GAIN (0.010) or more for the population against all recordings, on at least 6 of the 8 wirings
  (against all recordings) and for the population against each of halves A to D; the other ρ not lower by more than
  OTHER_SLACK (0.005) for the population against all recordings.
- The wrong-way control loses: on each ρ that gains, pow1.15's population change against all recordings is below zero.
  The same part-A control serves parts B and C (it is never refitted). A failure is named among the fails.
- The APL effect: the pattern part under the calibrated readout not higher by more than 0.003 against all recordings
  and each of halves A to D; the raw-count loss (all recordings) not higher by more than 0.005; the mean APL effect
  (calibrated) not lower by more than 0.01. The amplitude part, the shape part, the spread ratio, ρ and Pearson r on the
  per-pair APL effect go beside.
- Realism: every one of the input's eight runs inside the four bands (E83's `realism`, the bands of 2026-11-10; the
  flags under Yang's z rule beside).
- Response rates: response rate ρ with APL on and with APL off (the bands' units, the mean over the eight runs) not
  lower by more than 0.05.
- The block, the natural-source alignment with APL on, off and on minus off, and the structure against the wild-type
  somata are columns only: they enter no class (the entry: "The block is a column only").
- Classes (`class_of`): "better" (the gain, the control and every constraint hold); "a price" (the gain holds and a
  constraint or the control fails; each failure named with its size); "an equal" (no gain, every constraint holds and
  neither ρ lower than I0's by more than OTHER_SLACK); "no gain" (otherwise).

Choices the entry left open (fixed here before any run, with the reason):
- "Loses" for the wrong-way control: its population change below zero on each ρ the candidate gains on (the entry's
  word; E77's screen had 1.33 lower on both ρ's, so a sign test is the plainest reading of 1.15).
- Part C's exponent ("the largest gated gain"): among pow0.6 to pow0.9, those whose class is "better"; of them the one
  with the largest population change against all recordings on a ρ that gains; on a tie the exponent nearer 1. None
  "better": 0.75, as the entry says ("if none passes"), whose held points are run first (part Cstart) as the refit's
  start and enter no class. The saturating form is not refitted (the entry refits "the exponent").
- The refit's start: each wiring's held point of the chosen exponent (the fit's own parameters with the somatic
  threshold and the tied APL scale re-held for the new input; E110 part B's precedent), one start, 150 evaluations
  (E72's), KC-KC free. The bounds are E10's (relative to E3's base parameters, which do not depend on the input).
- Part C's gate: the refits' population against the family's saved fits (I0), with part A's pow1.15 as its control.
- Part D runs on part C's class "better" only ("only on a pass after the refit"); D1's wirings are T's at their refits;
  D2's reference is E113's input held the same way, and its pass is not gated (beside).
- Outcome letters (`outcome_of`, after part C; a provisional letter from part A before it):
  (a) part C's class is "better";
  (b) part C's class is "a price" with the calibrated pattern part among its fails, and part A's chosen exponent's also
      (the price at fixed scales and after the refit);
  (c) no exponent of part A has the population gain, but on the real wiring alone one ρ is higher by GAIN or more with
      the other not lower by more than OTHER_SLACK (the hemibrain wiring's own point, E77's);
  (d) no exponent of part A has the population gain, and the saturating form differs from the exponent: a sat<c> has
      the population gain, or the best population change of a ρ over the sat<c> differs from that over the exponents
      by GAIN or more;
  '-' none of these (reported with the classes).
- The reading E112 names: E112 had not ended when this was written. score gives, beside, the two ρ's of every input's
  population under E112's readings R0 to R3 (`E112_rho_like_for_like.run_readings` / `reading_scores`) with E112's
  index row state, and names none as carrying; that is E112's.
- E108's halves A to D are the halves (E77's `recording_halves(0)` and `(1)`).

Departures from the entry: none so far.

Provenance: the other checkouts must be clean; this repository clean in every file a process loads (E108's
`_commits_clean`; for the fits E110's `_fit_clean`, nothing modified under experiments/ or kcq/); every row records
`H.code_commits()`; every simulation part's provenance has `H.model_provenance(orn_deltas=...)` per point.

Outputs: <shared>/results/E114_input_compression/; spike counts and responses under outputs/E114_input_compression/
(A, B, Cs, D, beside, fits/), then data/model_outputs/E114_input_compression/.

Usage (from the repository root; one 16-thread lane on viki, under systemd-run --user; outputs/E114_input_compression/
chain.sh runs the parts in order):
    ./venv/bin/python -u experiments/E114_input_compression.py --part inputs
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part 0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part 0c
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part B
    ./venv/bin/python -u experiments/E114_input_compression.py --part score
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part Cstart
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part C
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part D
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part beside
    ./venv/bin/python -u experiments/E114_input_compression.py --part score
    ./venv/bin/python -u experiments/E114_input_compression.py --part variants
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E114_input_compression'
EXPONENTS = (0.6, 0.7, 0.8, 0.9)
WRONG_WAY = 1.15
FALLBACK_EXPONENT = 0.75
SAT_C = (100.0, 200.0, 400.0)


def power_map(orn: pd.DataFrame, e: float) -> pd.DataFrame:
    """pow<e> (E77's S3 rule, `E77_input_screen.orn_power`)."""
    raise NotImplementedError


def saturating_map(orn: pd.DataFrame, c: float) -> pd.DataFrame:
    """sat<c>: x / (1 + x / c) on each positive delta, rescaled to the odor's sum of positive deltas."""
    raise NotImplementedError


def input_of(name: str) -> pd.DataFrame:
    """An input by name ('I0', 'pow<e>', 'sat<c>')."""
    raise NotImplementedError


def write_inputs() -> dict:
    """Part inputs."""
    raise NotImplementedError


def identity(units: Optional[Sequence[str]] = None) -> dict:
    """Part 0."""
    raise NotImplementedError


def fit_path_check() -> dict:
    """Part 0c."""
    raise NotImplementedError


def part_sims(part: str, wirings: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Parts A, B and Cstart's held readings."""
    raise NotImplementedError


def readings() -> pd.DataFrame:
    """readings.csv and runs.csv."""
    raise NotImplementedError


def gate(rd: pd.DataFrame, runs: pd.DataFrame, name: str, ref: str = 'I0') -> dict:
    """The entry's gate for one input against I0."""
    raise NotImplementedError


def class_of(g: dict) -> str:
    """'better', 'a price', 'an equal' or 'no gain'."""
    raise NotImplementedError


def choose_c(gates: Dict[str, dict]) -> dict:
    """Part C's exponent."""
    raise NotImplementedError


def outcome_of(gates: Dict[str, dict], rd: pd.DataFrame) -> dict:
    """The entry's outcome letter."""
    raise NotImplementedError


def part_c(wirings: Optional[Sequence[str]] = None) -> List[str]:
    """Part C's fits."""
    raise NotImplementedError


def part_d(force: bool = False, e113_input: Optional[str] = None) -> dict:
    """Part D (only on a pass after the refit)."""
    raise NotImplementedError


def part_beside() -> dict:
    """The PN and bouton stages (no KC data)."""
    raise NotImplementedError


def score() -> dict:
    """Part score."""
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    """variants.csv and model_pair_distances_long.csv."""
    raise NotImplementedError


def main() -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
