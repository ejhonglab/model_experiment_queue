#!/usr/bin/env python3
"""E110: KC-KC inhibition held at zero in the best family: what the model's only APL-independent inhibition is needed for,
at fixed scales and with the other parameters fitted again.

The entry is item E110 of experiment_priorities.md (proposed by viki 2026-11-18 14:54; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-18 14:49; moved into the list 14:58 on Tom's words of 14:48, relayed by hal).
"The family" is C2's (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse) on E72 part B's eight wirings, each at
its own fit (E83's `saved_fit`; real, P1, P2, M1, P3, P4, M2, M3); "the population" is the mean matrix over the eight
wirings per condition (E72's, `E83.population`); "the two ρ's" are ρ on odor-odor distances with APL on against the
control flies and with APL off against APL output blocked (TNT), on raw counts (`E106_rho_gap_parts.two_rhos`, each KC
matrix the mean over recordings of the cosine distance). The family's saved fits are the reference and are not refitted.

Every fit is E72's / E67's part-B command, `E10_joint_fit.py -k 2 --aplpn-max 3 --rate-target ctrl_frac --loss apl
--apl-readout both --tag _cap3...` (150 evaluations unless named), with the lab's input and the wiring's options (`--wiring
claw_permutation_within_type --wiring-seed n` for Pn, E59's matched curveball matrix `--wiring-file ... --wiring-tag
curveball_within_type_matched-seed<n>` for Mn), run as its own process into its own folder: results and outputs of
`--experiment E110_no_kckc/fits/<name>` (fits.csv, history, provenance; spike counts and responses), so that hosts never
write one fits.csv over NFS.

Parts (CLI `--part`; `--wirings` takes a subset, so a part's points can be split over lanes and hosts; `--lane` names the
lane's tables, default the host's short name; a point is claimed by `mkdir` of results/E110_no_kckc/claims/<name>/ and
skipped when done or claimed, `--reclaim` taking over a claim left by a process that died):
- 0, the checks: (i) `--kckc-fixed` at the hemibrain wiring's fitted value (C2's kc_kc_inh) from C2's exact start
  (`--init-from-exact`) reproduces the fit at its first evaluation: every column of the first row of the check's history
  equals C2's history at its best evaluation, and the fit's row (fits.csv, `fit_<term>` and `loss`), exactly; (ii)
  `--kckc-fixed 0` runs (a smoke fit, `--smoke`: five evaluations), every evaluation at KC-KC 0. check.json. Run before
  any other part; each check is one E10_joint_fit.py process (results/E110_no_kckc/check/<name>/, its launch record
  launch.json).
- A (simulations): KC-KC at 0 at fixed scales on the eight wirings, both response rates held (`E77_input_screen.hold_rates`
  on `E83.build_simulator(w)` at the fit's parameters with wKCKC_scale 0: the somatic threshold bisected with APL off to
  the fit's APL-off response count, then the tied APL scales, APL->KC with KC->APL and bouton->APL, to its APL-on count;
  breadth weights and APL->bouton as fitted). Rows in runs_A_<lane>.csv (absolute parameters on and off, the factors,
  realism, `H.code_commits()`); runs under outputs/E110_no_kckc/A/ (tag A_<w>); each held point enters starts.json as
  A_<w>. analyses/yang_family_no_kckc.py (today's run for Yang) is not committed, so part A runs here, as the entry says;
  part score compares its held points with that run's beside (probe_check.csv), choosing nothing.
- B (fits): KC-KC held at 0 (`--kckc-fixed 0`), the other three parameters fitted (fixed_thr, APL->KC with KC->APL and
  bouton->APL tied, APL->bouton under its cap), started at part A's held point of the wiring (`--init-params
  results/E110_no_kckc/starts.json:A_<w>`, all five scales); on real, P1, P2 and M1 a second start (B2) from the fit's own
  parameters with KC-KC replaced by 0 (`--init-from <fit> --init-from-exact --kckc-fixed 0`), 200 evaluations. Of two
  starts the lower `fit_apl_cal_rmse` counts (E72's rule for the family's two starts); both are reported.
- C (fits, after B): the return fit on real, P1, P3 and M2: from part B's counted optimum with KC-KC started at 2 and
  free (starts.json:C_<w>), 150 evaluations.
- D (fits, only if part B's verdict is "needed"; `--force` overrides): half the fitted KC-KC held (`--kckc-fixed`) on
  real, P1, P2 and M1, from the fit's own parameters (exact start).
- E (fits): the breadth exponent at 1.5 and at 3 (`-k`) with KC-KC held at 0 on real, P1, P2 and M1, from part A's held
  point (as part B). E-rest: an exponent whose pattern part under the readout (the population of the four) is lower than
  exponent 2's (part B's counted runs on the same four) by E_ALT_GAIN or more is run on P3, P4, M2 and M3 and reported
  as the alternative; exponent 2 carries the verdict (`--force` runs E-rest anyway).
- score (no simulations): readings.csv (per kind of run, unit and readout: the two ρ's, against all recordings and E108's
  halves A to D, the APL-effect loss, its pattern part, amplitude part, shape part, spread ratio, ρ and Pearson r on the
  per-pair APL effect, the mean APL effect, the block, the natural-source alignment on, off and on minus off, the
  structure against wild-type somata (`E108_paa_extrapolation.unit_scores`); the runs' response rates, ratio, response
  rate ρ's, odors per KC, largest counts, KCs with 20 or more spikes to some odor, silent KCs, parameters; 'd_' columns
  the change from the family's same unit and readout), runs.csv, changes_B.csv (the verdict's measures: population, mean
  paired difference, each wiring), bands_B.csv, verdict.json, c_reading.csv, e_reading.csv, d_reading.csv,
  reference_check.json (the family's population against E108's I0), start_check.csv (each part-B fit's first
  evaluation against part A's held run), probe_check.csv, yardsticks_e106.csv, comparison_table.md,
  provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive and kcq/apl_effect.py.
- commands: prints the lanes' command lines and every fit's E10_joint_fit.py command; starts: (re)writes starts.json.

The verdict (fixed in the entry before any run; part B at exponent 2 against the family; `verdict`; the margins are
E108's, the 2026-11-17 15:40 answer):
- "not needed": the pattern part under the calibrated readout not higher by more than 0.003, the raw-count loss not
  higher by more than 0.005, the mean APL effect not lower by more than 0.01 under either reading, each ρ not lower by
  more than 0.005, response rate ρ with APL on and with APL off not lower by more than 0.05 (each on the population, as
  the mean paired difference over the eight wirings and on at least 6 of the 8 wirings), and every run inside the four
  bands;
- "needed": a margin fails on the population and on at least 6 of the 8 wirings (paired), or a band on at least 6 of the
  8 runs; the failing measures are named with their sizes, grouped as the main targets, the two ρ's and the APL-off
  measures (ρ with APL off, the APL-off response rate, response rate ρ with APL off, odors per KC with APL off);
- "better without": the pattern part lower by 0.003 or more, or one ρ higher by 0.010 or more, on the population and on
  at least 6 of the 8 wirings, with every "not needed" condition holding;
- "not determined": none of these, or part C's return fit ends within 0.002 of both the fit and part B on at least 3 of
  its 4 wirings;
- reported apart, pending Tom's answers of 2026-11-17 14:44 (not in the class): the block (0.010), the natural-source
  alignment with APL on, with APL off and their difference, and the structure against wild-type somata (0.05 each);
- E106's yardsticks beside (yardsticks_e106.csv). No threshold moves after a run is read.

Deviations and choices (fixed here before any run):
- The fits' folders: one per fit (`--experiment E110_no_kckc/fits/<name>`, <name> B_<w>, B2_<w>, C_<w>, D_<w>,
  E<k>_<w>); E10_joint_fit.py's tags add a suffix to E72's `_cap3` (`_cap3_B`, `_cap3_B2`, `_cap3_C`, `_cap3_D`, `_cap3_E`)
  and its own kckc<S>_ prefix. The archive keeps the layout: data/model_outputs/E110_no_kckc/A/ and
  data/model_outputs/E110_no_kckc/fits/<name>/; the scoring reads this host's outputs/E110_no_kckc/ first, then those.
- Which of part B's two starts counts: the lower `fit_apl_cal_rmse` (the calibrated APL-effect loss of the fit's own
  loss; E72's `new_fit_members` chose between the family's starts by it). Part C's 0.002 comparisons are on the same
  column: "returns" is KC-KC above 10 with `fit_apl_cal_rmse` within 0.002 of the family fit's; "flat" is within 0.002 of
  both the fit's and part B's counted fit's; "more evaluations explain a gain" is KC-KC at most 10 with C lower than part
  B by more than 0.002. The verdict's "not determined" from part C needs "flat" on 3 of its 4 wirings (E108's
  `wirings_needed`, the 6-of-8 share); until part C is complete the verdict says "part C pending".
- Part B's second start and part D start from the fit's own parameters by `--init-from-exact` with `--kckc-fixed`
  (E10_joint_fit.py's documented behaviour: the seed fit's parameters with its KC-KC scale replaced), the path part 0
  checks.
- Part E starts at part A's held point (as part B), so that it differs from part B only in the exponent and needs no
  part-B result; its rule reads the population of real, P1, P2 and M1 under the calibrated readout.
- "Not needed" needs each margin to hold on the population, on the mean paired difference and on 6 of the 8 wirings
  (the entry's "neither holds on 6 of 8 wirings" reads both classes on the wirings; the header names the population and
  the mean paired difference); "needed" and the gains of "better without" are read on the population and on 6 of 8
  wirings. A margin that holds on the population but on fewer than 6 wirings, without failing on 6, leaves the verdict
  "not determined". The response rate ρ's population value is the mean over the
  eight runs (E108). The bands are A.realism_flag's five checks (APL-on rate, ratio, response rate ρ with APL on,
  odors per KC on and off), each read from the run's numbers with kcq.apl_effect's constants.
- The pending items are reported apart and do not enter the class (the entry: "reported apart").
- The ρ's of the verdict are `E106.two_rhos` (raw counts, all recordings); readings.csv carries E108's `unit_scores`
  values beside them (the same matrices; reference_check.json records any difference).

Provenance: simulation parts and fit launches need the other checkouts clean and this repository clean except for
files no process of the part loads (another session's work in progress in the same checkout; `E108._commits_clean`,
`_fit_clean` refuses any modified file under experiments/ or kcq/ too, since E10_joint_fit.py loads modules lazily); each
row records `H.code_commits()` with 'queue_dirty_paths'; part A's provenance_A_<lane>.json has `H.model_provenance` per
point; each fit's folder has E10_joint_fit.py's own provenance and the launch record (launch.json: command, host, commits).

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): only the model's hemibrain wiring and E72's seven re-draws of
it; no other connectome table; no raw-only recording. Nothing is published to claude.ai.

Usage (from the repository root; threads = odors, 16, or fewer on a shared machine; under nice). Part 0 runs once (its
check.json is shared); then any number of lanes on any hosts, each running the same parts in this order (the claims
split the points; `--wirings` restricts a lane; `--wait` queues a lane behind starts another lane is still making):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E110_no_kckc.py --part 0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E110_no_kckc.py --part A --lane hal1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E110_no_kckc.py --part B --lane hal1 --wait
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E110_no_kckc.py --part E --lane hal1 --wait
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E110_no_kckc.py --part C --lane hal1 --wait
    ./venv/bin/python -u experiments/E110_no_kckc.py --part commands      # every fit's E10_joint_fit.py command
then, on every host that ran a part (the scoring reads data/model_outputs/E110_no_kckc/A/ and .../fits/<name>/):
    ./venv/bin/python archive_tools/archive_model_outputs.py E110_no_kckc --exclude 'check,distances_cache.pkl,fits/*/*base_acute_off_spike_counts.parquet'
    ./venv/bin/python archive_tools/migrate_model_output_params.py E110_no_kckc --apply
and, on any host:
    ./venv/bin/python -u experiments/E110_no_kckc.py --part score
    ./venv/bin/python -u experiments/E110_no_kckc.py --part variants      # before the archive's params.csv
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import shutil
import socket
import subprocess
import time
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E108_paa_extrapolation as E108

EXPERIMENT = 'E110_no_kckc'
WIRINGS = E83.ALL_WIRINGS
S_WIRINGS = E83.S_WIRINGS                          # real, P1, P2, M1
B2_WIRINGS = ('real', 'P1', 'P2', 'M1')
C_WIRINGS = ('real', 'P1', 'P3', 'M2')
D_WIRINGS = ('real', 'P1', 'P2', 'M1')
E_WIRINGS = ('real', 'P1', 'P2', 'M1')
E_REST = ('P3', 'P4', 'M2', 'M3')
E_EXPONENTS = (1.5, 3.0)
K = 2.0
MAXFEV, MAXFEV_B2 = 150, 200
C_KCKC_START = 2.0
C_RETURN_KCKC = 10.0
C_WITHIN = 0.002
E_ALT_GAIN = 0.003
CHOOSE_COL = 'fit_apl_cal_rmse'
UNITS = {'all8': WIRINGS, 'S4': S_WIRINGS, 'C4': C_WIRINGS, 'T4': E83.T_WIRINGS}
READOUTS = ('raw', 'cal')
PARAM_KEYS = E83.PARAM_KEYS                        # fixed_thr, the four APL scales, wKCKC_scale
# the verdict's margins (E108's: the 2026-11-17 15:40 answer)
PATTERN_SLACK, RAW_SLACK, MEAN_SLACK = E108.EQUAL_PATTERN, E108.EQUAL_RAW, E108.EQUAL_MEAN
RHO_SLACK, RATE_RHO_SLACK = E108.OTHER_SLACK, E108.RATE_RHO_SLACK
BETTER_PATTERN, BETTER_RHO = 0.003, E108.GAIN
N_WIRINGS = E108.GAIN_WIRINGS                      # 6 of the 8
BLOCK_SLACK, STRUCTURE_SLACK = E108.BLOCK_SLACK, E108.STRUCTURE_SLACK
# the family's population on raw counts as the entry quotes it (E108's I0)
ENTRY_REFERENCE = {'rho_on': 0.690, 'rho_off': 0.756}
FIT_EXPERIMENT_PREFIX = f'{EXPERIMENT}/fits'
CHECK_EXPERIMENT_PREFIX = f'{EXPERIMENT}/check'

Measure = namedtuple('Measure', 'key label col readout worse slack groups')
# worse: +1 when higher is worse, -1 when lower is worse
MEASURES = (
    Measure('pattern_cal', 'pattern part (readout)', 'all_pattern', 'cal', +1, PATTERN_SLACK, ('main',)),
    Measure('loss_raw', 'APL-effect loss (raw counts)', 'all_loss', 'raw', +1, RAW_SLACK, ('main',)),
    Measure('mean_raw', 'mean APL effect (raw counts)', 'mean_effect', 'raw', -1, MEAN_SLACK, ('main',)),
    Measure('mean_cal', 'mean APL effect (readout)', 'mean_effect', 'cal', -1, MEAN_SLACK, ('main',)),
    Measure('rho_on', 'ρ on odor-odor distances, APL on', 'rho_on', 'raw', -1, RHO_SLACK, ('rho',)),
    Measure('rho_off', 'ρ on odor-odor distances, APL off', 'rho_off', 'raw', -1, RHO_SLACK, ('rho', 'apl_off')),
    Measure('rate_rho_on', 'response rate ρ, APL on', 'run_rate_rho_on', 'raw', -1, RATE_RHO_SLACK, ('rates',)),
    Measure('rate_rho_off', 'response rate ρ, APL off', 'run_rate_rho_off', 'raw', -1, RATE_RHO_SLACK, ('apl_off',)),
)
MEASURE = {m.key: m for m in MEASURES}
Band = namedtuple('Band', 'key label groups')
BANDS = (
    Band('band_rate_on', 'APL-on response rate (band)', ('realism',)),
    Band('band_ratio', 'APL-off / APL-on response rate (band)', ('realism', 'apl_off')),
    Band('band_rate_rho_on', 'response rate ρ, APL on (band)', ('realism',)),
    Band('band_opk_on', 'odors per KC, APL on (band)', ('realism',)),
    Band('band_opk_off', 'odors per KC, APL off (band)', ('realism', 'apl_off')),
)
GROUPS = {'main': 'the main targets', 'rho': "the two ρ's", 'apl_off': 'the APL-off measures', 'rates': 'response rates',
          'realism': 'realism bands'}
# reported apart (pending Tom's answers of 2026-11-17 14:44): key, label, readings column, slack (lower is worse)
PENDING = (('block', 'the block', 'block', BLOCK_SLACK),
           ('vcf_on', 'natural-source alignment, APL on', 'vcf_on', STRUCTURE_SLACK),
           ('vcf_off', 'natural-source alignment, APL off', 'vcf_off', STRUCTURE_SLACK),
           ('vcf_diff', 'natural-source alignment, APL on minus off', 'vcf_diff', STRUCTURE_SLACK),
           ('struct_wt_on', 'structure vs wild-type somata', 'struct_wt_on', STRUCTURE_SLACK))

CLASSES = ('not needed', 'needed', 'better without', 'not determined')

_CACHE: Dict[object, object] = {}

FitSpec = namedtuple('FitSpec', 'name part wiring k kckc_fixed maxfev start suffix')
# start: ('params', <starts.json key>) or ('exact',) (the wiring's own fit, --init-from-exact)


# ---------------------------------------------------------------------------
# folders and small helpers

def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    return H.outputs_dir(f'{EXPERIMENT}/{part}' if part else EXPERIMENT)


def archive_dir() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M:%S %Z')


# ---------------------------------------------------------------------------
# the wirings and the fits' commands

def wiring_prefix(w: str) -> str:
    """The wiring's part of E10_joint_fit.py's tag: '' for the hemibrain wiring, E72's DRAWS prefix for a re-draw
    (wiring-claw_permutation_within_type-seed<n>_, wiring-curveball_within_type_matched-seed<n>_)."""
    from experiments import E72_wiring_population as E72
    if w == 'real':
        return ''
    return dict(E72.DRAWS)[w]


def wiring_args(w: str) -> List[str]:
    """E10_joint_fit.py's wiring options for a wiring label, as E72 part B ran them: none for 'real'; `--wiring
    claw_permutation_within_type --wiring-seed n` for Pn; E59's matched curveball matrix of this host's shared directory
    with its `--wiring-tag` for Mn (E83's `wiring_kws`)."""
    if w == 'real':
        return []
    pre = wiring_prefix(w)
    seed = int(pre.split('seed')[1].rstrip('_'))
    if w.startswith('P'):
        return ['--wiring', 'claw_permutation_within_type', '--wiring-seed', str(seed)]
    f = H.shared_dir() / 'data' / 'model_outputs' / 'E59_coinput_shuffle' / f'wPNKC_curveball_within_type-seed{seed}_matched.parquet'
    return ['--wiring-file', str(f), '--wiring-tag', f'curveball_within_type_matched-seed{seed}']


def e10_tag(w: str, k: float, kckc_fixed: Optional[float], suffix: str, rate: Optional[float] = None) -> str:
    """The tag E10_joint_fit.py gives a fit of this experiment (its prefix rules for the options used here: the wiring,
    rate<R>_, aplboth_, kckc<S>_, then k<k> and the --tag suffix). `rate`: the APL-on rate target (default the registry's
    ctrl_frac, `--rate-target ctrl_frac`)."""
    r = float(H.data_target('ctrl_frac')) if rate is None else float(rate)
    kc = '' if kckc_fixed is None else f'kckc{float(kckc_fixed):g}_'
    return f'{wiring_prefix(w)}rate{r:.3f}_aplboth_{kc}k{float(k):g}{suffix}'


def fit_tag(spec: FitSpec, rate: Optional[float] = None) -> str:
    return e10_tag(spec.wiring, spec.k, spec.kckc_fixed, f'_cap3_{spec.suffix}', rate)


def half_kckc(w: str) -> float:
    """Part D's KC-KC scale: half the wiring's fitted one."""
    return 0.5 * float(E83.saved_fit(w)['p']['wKCKC_scale'])


def fit_specs(part: str, wirings: Optional[Sequence[str]] = None, exponents: Sequence[float] = E_EXPONENTS
              ) -> List[FitSpec]:
    """The fits of a part ('B', 'C', 'D', 'E', 'E-rest'), in order; `wirings` keeps a subset."""
    out: List[FitSpec] = []
    if part == 'B':
        for w in WIRINGS:
            out.append(FitSpec(f'B_{w}', 'B', w, K, 0.0, MAXFEV, ('params', f'A_{w}'), 'B'))
            if w in B2_WIRINGS:
                out.append(FitSpec(f'B2_{w}', 'B2', w, K, 0.0, MAXFEV_B2, ('exact',), 'B2'))
    elif part == 'C':
        out = [FitSpec(f'C_{w}', 'C', w, K, None, MAXFEV, ('params', f'C_{w}'), 'C') for w in C_WIRINGS]
    elif part == 'D':
        out = [FitSpec(f'D_{w}', 'D', w, K, half_kckc(w), MAXFEV, ('exact',), 'D') for w in D_WIRINGS]
    elif part in ('E', 'E-rest'):
        ws = E_WIRINGS if part == 'E' else E_REST
        out = [FitSpec(f'E{float(k):g}_{w}', part, w, float(k), 0.0, MAXFEV, ('params', f'A_{w}'), 'E')
               for k in exponents for w in ws]
    else:
        raise ValueError(f'no fits in part {part!r}')
    return [s for s in out if wirings is None or s.wiring in wirings]


def starts_spec(key: str) -> str:
    """`--init-params`' FILE:KEY of a start in starts.json (FILE relative to the shared directory)."""
    return f'results/{EXPERIMENT}/starts.json:{key}'


def fit_command(spec: FitSpec, python: Optional[str] = None) -> List[str]:
    """The E10_joint_fit.py command of one fit (E72's part-B command plus this fit's options), into its own folders."""
    py = str(H.HERE / 'venv' / 'bin' / 'python') if python is None else python
    cmd = [py, '-u', 'experiments/E10_joint_fit.py', '-k', f'{float(spec.k):g}', '--aplpn-max', '3', '--rate-target',
           'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', f'_cap3_{spec.suffix}', '--maxfev',
           str(int(spec.maxfev)), '--experiment', f'{FIT_EXPERIMENT_PREFIX}/{spec.name}'] + wiring_args(spec.wiring)
    if spec.kckc_fixed is not None:
        cmd += ['--kckc-fixed', repr(float(spec.kckc_fixed))]
    if spec.start[0] == 'params':
        cmd += ['--init-params', starts_spec(spec.start[1])]
    elif spec.start[0] == 'exact':
        folder, tag = E83.members()[spec.wiring]
        cmd += ['--init-from', tag, '--init-from-experiment', folder, '--init-from-exact']
    else:
        raise ValueError(f'unknown start {spec.start}')
    return cmd


def fit_results_dir(name: str, check: bool = False) -> Path:
    """results/E110_no_kckc/fits/<name>/ (or check/<name>/), not created."""
    return H.shared_dir() / 'results' / (CHECK_EXPERIMENT_PREFIX if check else FIT_EXPERIMENT_PREFIX) / name


def fit_row(name: str, tag: str) -> Optional[pd.Series]:
    """The fit's row in its folder's fits.csv, or None."""
    p = fit_results_dir(name) / 'fits.csv'
    if not p.exists():
        return None
    f = pd.read_csv(p)
    r = f[f['tag'] == tag]
    return r.iloc[-1] if len(r) else None


def load_fit_outputs(name: str, tag: str) -> Optional[Dict[str, tuple]]:
    """A fit's saved runs ({'on', 'off'}: (responses, spike counts)): this host's outputs/E110_no_kckc/fits/<name>/, else
    the archive's data/model_outputs/E110_no_kckc/fits/<name>/; or None."""
    for d in (H.HERE / 'outputs' / FIT_EXPERIMENT_PREFIX / name, archive_dir() / 'fits' / name):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def load_a(w: str) -> Optional[Dict[str, tuple]]:
    """Part A's run of a wiring: this host's outputs/E110_no_kckc/A/, else the archive's A/."""
    for d in (H.HERE / 'outputs' / EXPERIMENT / 'A', archive_dir() / 'A'):
        o = E83._load_pair(d, f'A_{w}')
        if o is not None:
            return o
    return None


def params_of_row(r: Mapping) -> Dict[str, float]:
    """Absolute parameters of a fits.csv row (E72's `fitted_params`): the five scales and wKCKC_scale (kc_kc_inh)."""
    p = {k: float(r[k]) for k in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale')}
    p['wKCKC_scale'] = float(r['kc_kc_inh'])
    return p


# ---------------------------------------------------------------------------
# the starts file

def start_entry(params: Mapping[str, float], **meta) -> dict:
    """One starts.json entry: the absolute parameters (PARAM_KEYS) and string metadata, which E10_joint_fit.read_params
    ignores."""
    e = {k: float(params[k]) for k in PARAM_KEYS}
    e.update({k: str(v) for k, v in meta.items()})
    return e


def a_start(row: Mapping) -> dict:
    """Part A's held point (a runs_A row) as part B's and part E's start: its APL-on parameters, KC-KC 0."""
    p = json.loads(row['params_on']) if isinstance(row['params_on'], str) else dict(row['params_on'])
    if float(p['wKCKC_scale']) != 0.0:
        raise ValueError(f"part A's point {row['tag']} has KC-KC {p['wKCKC_scale']}")
    return start_entry(p, source=f"part A's held point {row['tag']} (results/{EXPERIMENT}/runs_A_*.csv)")


def c_start(row: Mapping, name: str) -> dict:
    """Part C's start: part B's counted optimum (a fits.csv row) with KC-KC set to C_KCKC_START."""
    p = params_of_row(row)
    p['wKCKC_scale'] = C_KCKC_START
    return start_entry(p, source=f"part B's counted optimum {name}:{row['tag']} with KC-KC {C_KCKC_START:g}")


class _Lock:
    """A lock by mkdir of `<path>.lock` (atomic on NFS); waits up to `wait` s, and takes over a lock older than `stale`
    s (a process that died holding it)."""

    def __init__(self, path: Path, wait: float = 120.0, stale: float = 120.0):
        self.d, self.wait, self.stale = path.with_name(path.name + '.lock'), wait, stale

    def __enter__(self):
        t0 = time.time()
        while True:
            try:
                self.d.mkdir(parents=True)
                return self
            except FileExistsError:
                try:
                    if time.time() - self.d.stat().st_mtime > self.stale:
                        shutil.rmtree(self.d, ignore_errors=True)
                        continue
                except FileNotFoundError:
                    continue
                if time.time() - t0 > self.wait:
                    raise TimeoutError(f'{self.d} held for more than {self.wait:.0f} s')
                time.sleep(0.5)

    def __exit__(self, *exc):
        shutil.rmtree(self.d, ignore_errors=True)
        return False


def merge_starts(path: Path, entries: Mapping[str, dict]) -> dict:
    """Adds `entries` to the starts file `path` (JSON {key: entry}) under a lock and returns the merged content. An
    existing key must have the same parameters (RuntimeError otherwise: a start never changes once a fit may have used
    it); keys are never removed. Written to a temporary file and moved into place (atomic on NFS), so a reader never sees
    half a file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with _Lock(path):
        cur = _read_json(path)
        changed = False
        for k, e in entries.items():
            if k in cur:
                diff = {p: (cur[k].get(p), e[p]) for p in PARAM_KEYS if float(cur[k].get(p, np.nan)) != float(e[p])}
                if diff:
                    raise RuntimeError(f'starts.json: {k} would change {diff} (a start is fixed once written)')
                continue
            cur[k] = dict(e)
            changed = True
        if changed or not path.exists():
            tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
            tmp.write_text(json.dumps(cur, indent=1) + '\n')
            os.replace(tmp, path)
    return cur


def starts_path() -> Path:
    return out_dir() / 'starts.json'


def write_starts() -> dict:
    """starts.json from what is saved: A_<w> for part A's held points, C_<w> for part B's counted optima of C's wirings."""
    entries = {}
    ra = runs_table_of('A')
    for _, r in ra.iterrows():
        entries[f"A_{r['wiring']}"] = a_start(r)
    for w in C_WIRINGS:
        b = counted_b(w)
        if b is not None:
            entries[f'C_{w}'] = c_start(b['row'], b['name'])
    return merge_starts(starts_path(), entries)


# ---------------------------------------------------------------------------
# claims and lane tables

def claim(name: str, lane: str, reclaim: bool = False) -> bool:
    """Claims a point by mkdir of results/E110_no_kckc/claims/<name>/ (atomic on NFS); False if another process holds it
    (unless `reclaim`). The claim records host, lane, pid and time."""
    d = out_dir() / 'claims' / name
    try:
        d.mkdir(parents=True)
    except FileExistsError:
        if not reclaim:
            return False
    _write_json(d / 'claim.json', dict(name=name, host=_host(), lane=lane, pid=os.getpid(), time=_now()))
    return True


def release(name: str) -> None:
    shutil.rmtree(out_dir() / 'claims' / name, ignore_errors=True)


def runs_table_of(kind: str) -> pd.DataFrame:
    """Every lane's rows of a kind ('A' or 'fits'): runs_<kind>_<lane>.csv, last row per tag."""
    ps = sorted(out_dir().glob(f'runs_{kind}_*.csv'))
    if not ps:
        return pd.DataFrame(columns=['tag', 'wiring', 'name'])
    return pd.concat([E108._read_rows(p) for p in ps], ignore_index=True, sort=False).drop_duplicates('tag', keep='last')


def _fit_clean() -> dict:
    """`E108._commits_clean()` (the other checkouts clean; this repository's modified files none that this process
    loads), and no modified file under experiments/ or kcq/ at all, since E10_joint_fit.py loads modules lazily."""
    cc = E108._commits_clean()
    bad = [p for p in str(cc.get('queue_dirty_paths') or '').split() if p.startswith(('experiments/', 'kcq/'))]
    if bad:
        raise RuntimeError(f'commit first: {bad} (an E10_joint_fit.py process may load them)')
    return cc


# ---------------------------------------------------------------------------
# part 0: the checks

def check_specs() -> Dict[str, Tuple[str, List[str]]]:
    """Part 0's two E10_joint_fit.py smoke runs on the hemibrain wiring: {name: (tag, command)}."""
    folder, tag = E83.members()['real']
    kc = float(E83.saved_fit('real')['p']['wKCKC_scale'])
    start = ['--init-from', tag, '--init-from-experiment', folder, '--init-from-exact']
    out = {}
    for name, kckc, suffix in (('exact_real', kc, 'check0'), ('smoke0_real', 0.0, 'smoke0')):
        cmd = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-max', '3',
               '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', f'_cap3_{suffix}', '--maxfev',
               str(MAXFEV), '--experiment', f'{CHECK_EXPERIMENT_PREFIX}/{name}'] + start + ['--kckc-fixed', repr(kckc), '--smoke']
        out[name] = (e10_tag('real', K, kckc, f'_cap3_{suffix}'), cmd)
    return out


def _child_env() -> dict:
    """The environment of an E10_joint_fit.py process: this one's, with OMP_THREAD_LIMIT set to OMP_NUM_THREADS when only
    the latter is set (importing pyolfsysm.mb_model resets OpenMP to every core otherwise; CLAUDE.md, "Gotchas")."""
    env = dict(os.environ)
    if 'OMP_NUM_THREADS' in env and 'OMP_THREAD_LIMIT' not in env:
        env['OMP_THREAD_LIMIT'] = env['OMP_NUM_THREADS']
    env.setdefault('PYTHONUNBUFFERED', '1')
    return env


def launch(cmd: List[str], res: Path, log: Path, cc: dict, name: str) -> dict:
    """Runs one E10_joint_fit.py process (log `log`), with its launch record res/launch.json (command, host, start, the
    commits, the queue's modified paths; then end, return code, seconds)."""
    rec = dict(name=name, command=cmd, cwd=str(H.HERE), host=_host(), started=_now(), commits=cc,
               queue_dirty_paths=str(cc.get('queue_dirty_paths') or '').split(), log=str(log.relative_to(H.HERE)),
               threads=dict(OMP_NUM_THREADS=os.environ.get('OMP_NUM_THREADS'), OMP_THREAD_LIMIT=_child_env().get('OMP_THREAD_LIMIT')),
               launcher='experiments/E110_no_kckc.py')
    res.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    _write_json(res / 'launch.json', rec)
    t0 = time.time()
    with open(log, 'w') as fh:
        r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT, env=_child_env())
    rec.update(finished=_now(), returncode=r.returncode, seconds=time.time() - t0)
    _write_json(res / 'launch.json', rec)
    return rec


def same_command(a: Sequence[str], b: Sequence[str]) -> bool:
    """Two command lines alike, argument by argument (as strings, or as numbers: '0' and '0.0' are the same value)."""
    if len(a) != len(b):
        return False
    for x, y in zip(a, b):
        if x == y:
            continue
        try:
            if float(x) != float(y):
                return False
        except ValueError:
            return False
    return True


def _num_same(a, b) -> bool:
    """Two cells equal exactly (numbers; NaN equals NaN; otherwise as strings)."""
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)
    return (np.isnan(fa) and np.isnan(fb)) or fa == fb


def exact_check(h0: Mapping, best: Mapping, row: Mapping) -> dict:
    """Check (i): `h0` (the check's first history row) against `best` (the fit's history at its best evaluation) on every
    column they share but 'eval' and 'time_s', and against the fit's row (`fit_<col>` and 'loss'). Returns the
    differences and 'exact'."""
    skip = ('eval', 'time_s', 'Unnamed: 0')
    keys = list(h0.keys())
    cols = [c for c in keys if c in best and c not in skip]
    d_hist = {c: [float(h0[c]), float(best[c])] for c in cols if not _num_same(h0[c], best[c])}
    rcols = [c for c in keys if c not in skip and (f'fit_{c}' in row or c == 'loss')]
    d_row = {c: [float(h0[c]), float(row['loss' if c == 'loss' else f'fit_{c}'])] for c in rcols
             if not _num_same(h0[c], row['loss' if c == 'loss' else f'fit_{c}'])}
    return dict(n_history_columns=len(cols), n_row_columns=len(rcols), history_differs=d_hist, row_differs=d_row,
                exact=bool(cols and rcols and not d_hist and not d_row))


def part_0(rerun: bool = False) -> dict:
    """Part 0 (the module docstring): runs each check unless its history and launch record exist (`rerun` runs it again),
    then reads them; check.json. Raises if a check fails."""
    specs = check_specs()
    res = dict(_read_json(out_dir() / 'check.json'))
    cc = None
    for name, (tag, cmd) in specs.items():
        d = fit_results_dir(name, check=True)
        hp = d / f'history_{tag}.csv'
        rec = _read_json(d / 'launch.json')
        if (rerun or not hp.exists() or rec.get('returncode') is None
                or not same_command(rec.get('command') or [], cmd)):
            cc = cc or _fit_clean()
            print(f'E110 0 {name}: ' + ' '.join(cmd), flush=True)
            rec = launch(cmd, d, local_dir('logs') / f'check_{name}.log', cc, name)
        res[name] = dict(tag=tag, command=' '.join(rec['command']), launch={k: rec.get(k) for k in
                         ('host', 'started', 'finished', 'returncode', 'seconds', 'log', 'launcher', 'threads')},
                         commits=rec.get('commits'), queue_dirty_paths=rec.get('queue_dirty_paths'))
        if rec.get('returncode') != 0 or not hp.exists():
            res[name]['passes'] = False
            continue
        h = pd.read_csv(hp)
        res[name].update(n_evaluations=int(len(h)), seconds_per_evaluation=float(h['time_s'].mean()) if 'time_s' in h else None,
                         kckc_values=sorted({float(x) for x in h['wKCKC_scale']}))
        if name == 'exact_real':
            folder, ftag = E83.members()['real']
            hist = pd.read_csv(H.results_dir(folder) / f'history_{ftag}.csv')
            best = hist.loc[hist['loss'].idxmin()]
            row = pd.read_csv(H.results_dir(folder) / 'fits.csv').set_index('tag').loc[ftag]
            ex = exact_check(h.iloc[0], best, row)
            res[name].update(reference=f'{folder}:{ftag} (history row {int(best["eval"])}, the best; fits.csv row)', **ex)
            res[name]['passes'] = ex['exact']
        else:
            res[name]['passes'] = bool(len(h) >= 2 and all(float(x) == 0.0 for x in h['wKCKC_scale']))
    res['all_pass'] = bool(all(res.get(n, {}).get('passes') for n in specs))
    res['checked'] = _now()
    res['evaluated_by'] = dict(script='experiments/E110_no_kckc.py', host=_host(), **H.code_commits())
    _write_json(out_dir() / 'check.json', res)
    print(json.dumps({n: {k: res[n].get(k) for k in ('passes', 'n_evaluations', 'seconds_per_evaluation', 'history_differs',
                                                      'row_differs')} for n in specs}, indent=1, default=str), flush=True)
    if not res['all_pass']:
        raise RuntimeError('E110 part 0 fails (check.json); nothing else runs until it passes')
    return res


def check_passed() -> bool:
    return bool(_read_json(out_dir() / 'check.json').get('all_pass'))


# ---------------------------------------------------------------------------
# part A: fixed scales, both response rates held

def run_extras(o: Dict[str, tuple]) -> dict:
    """Per condition: silent KCs (responding to no odor; two or more spikes), KCs with 20 or more spikes to some odor."""
    row = {}
    for c in ('on', 'off'):
        r = o[c][0].to_numpy(dtype=float) > 0.5
        s = o[c][1].to_numpy(dtype=float)
        row[f'silent_{c}'] = float((r.sum(axis=1) == 0).mean())
        row[f'n_kc_ge20_{c}'] = int((s.max(axis=1) >= 20).sum())
    return row


def save_a(w: str, o: Dict[str, tuple]) -> None:
    d = local_dir('A')
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'A_{w}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'A_{w}_{cc}_spike_counts.parquet')


def part_a(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, reclaim: bool = False) -> pd.DataFrame:
    """Part A (the module docstring), resumable: a wiring is done when its row is in any lane's runs_A table."""
    from experiments.E77_input_screen import hold_rates
    if not check_passed():
        raise RuntimeError('run --part 0 first')
    lane = lane or _host()
    cc = E108._commits_clean()
    path = out_dir() / f'runs_A_{lane}.csv'
    prov_path = out_dir() / f'provenance_A_{lane}.json'
    for w in wirings:
        tag = f'A_{w}'
        if tag in set(runs_table_of('A')['tag']):
            print(f'E110 A {w}: done', flush=True)
            continue
        if not claim(tag, lane, reclaim):
            print(f'E110 A {w}: claimed by another process, skipped', flush=True)
            continue
        try:
            E108._still_clean(cc)
            P = E83.wiring_parts(w)
            res = hold_rates(lambda: E83.build_simulator(w), target_off=P['target_off'], target_on=P['target_on'],
                             params=dict(P['p'], wKCKC_scale=0.0), label=f'E110 A {w}')
            o = {'on': res['on'], 'off': res['off']}
            save_a(w, o)
            row = dict(part='A', tag=tag, wiring=w, fit_folder=P['folder'], fit_tag=P['tag'],
                       fitted_wKCKC_scale=float(P['p']['wKCKC_scale']), thr_factor=res['thr_factor'],
                       apl_factor=res['apl_factor'], n_on=res['n_on'], n_off=res['n_off'], target_on=res['target_on'],
                       target_off=res['target_off'], at_bound_thr=res['at_bound_thr'], at_bound_apl=res['at_bound_apl'],
                       n_sims=res['n_sims'], seconds=res.get('seconds'), **{k: res['params'][k] for k in PARAM_KEYS},
                       params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
                       params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}),
                       K_run=E83.readout_K(o), breadth_sha256=res.get('breadth_sha256'), host=_host(), lane=lane)
            row.update(E83.realism(o))
            row.update(run_extras(o))
            row.update(cc)
            E83._append_row(path, row)
            prov = _read_json(prov_path) or dict(part='A', lane=lane, script='experiments/E110_no_kckc.py', points={})
            prov['points'][tag] = dict(H.model_provenance(orn_deltas=P['orn']), held=dict(
                change="wKCKC_scale = 0; fixed_thr re-tuned with APL off to the fit's APL-off response count, then the tied "
                       "APL scales (APL->KC, KC->APL, bouton->APL) to its APL-on count (E77_input_screen.hold_rates)",
                fit=f"{P['folder']}:{P['tag']}", kws=res.get('kws')))
            _write_json(prov_path, prov)
            merge_starts(starts_path(), {tag: a_start(row)})
            print(f"E110 A {w}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} "
                  f"on, {res['n_off']}/{res['target_off']} off); flag '{row.get('realism_flag') or ''}'; "
                  f"{res.get('seconds', 0):.0f}s", flush=True)
            del res, o
            gc.collect()
        except BaseException:
            release(tag)
            raise
    return runs_table_of('A')


# ---------------------------------------------------------------------------
# the fits (parts B, C, D, E)

def fit_done(spec: FitSpec) -> bool:
    return fit_row(spec.name, fit_tag(spec)) is not None


def counted_b(w: str) -> Optional[dict]:
    """Part B's counted fit of a wiring: of B_<w> and B2_<w> (where it exists) that done, the lower CHOOSE_COL
    ({'name', 'tag', 'row', 'n_starts'}); None when B_<w> is not done, or B2 is planned and not done."""
    cands = []
    for s in fit_specs('B', [w]):
        r = fit_row(s.name, fit_tag(s))
        if r is None:
            return None
        cands.append(dict(name=s.name, tag=fit_tag(s), row=r))
    best = choose_counted([c['row'] for c in cands])
    out = dict(cands[best])
    out['n_starts'] = len(cands)
    return out


def choose_counted(rows: Sequence[Mapping]) -> int:
    """The index of the counted start: the lowest CHOOSE_COL (the first on a tie)."""
    v = [float(r[CHOOSE_COL]) for r in rows]
    return int(np.argmin(v))


def _prerequisite(spec: FitSpec) -> Optional[str]:
    """Why a fit cannot start yet (None when it can)."""
    if spec.start[0] == 'params' and spec.start[1] not in _read_json(starts_path()):
        return f'start {spec.start[1]} not in starts.json (part {"A" if spec.start[1].startswith("A_") else "B"} first)'
    return None


WAIT_SECONDS = 300


def run_fits(specs: Sequence[FitSpec], lane: Optional[str] = None, reclaim: bool = False, dry_run: bool = False,
             wait: bool = False) -> List[str]:
    """Runs the fits one process after another (`fit_command`; logs under outputs/E110_no_kckc/logs/), skipping those
    done (their row in their folder's fits.csv) or claimed by another process; a row per fit in runs_fits_<lane>.csv.
    With `wait`, fits whose start is not saved yet (part A's or part B's, on another lane or host) are retried every
    WAIT_SECONDS until none is left waiting."""
    if not dry_run and not check_passed():
        raise RuntimeError('run --part 0 first')
    # the commits are taken once: a pull that changes a loaded file stops the lane (`E108._still_clean`)
    cc = None if dry_run else _fit_clean()
    while True:
        done, waiting = _run_fits_pass(specs, lane, reclaim, dry_run, cc)
        if not (wait and waiting) or dry_run:
            return done
        print(f'E110: waiting for the starts of {waiting} ({WAIT_SECONDS} s)', flush=True)
        time.sleep(WAIT_SECONDS)


def _run_fits_pass(specs: Sequence[FitSpec], lane: Optional[str], reclaim: bool, dry_run: bool, cc: Optional[dict]
                   ) -> Tuple[List[str], List[str]]:
    lane = lane or _host()
    path = out_dir() / f'runs_fits_{lane}.csv'
    done, waiting = [], []
    for spec in specs:
        tag = fit_tag(spec)
        if fit_done(spec):
            print(f'E110 {spec.name}: done ({tag})', flush=True)
            done.append(spec.name)
            continue
        if not dry_run:
            write_starts()
        why = _prerequisite(spec)
        cmd = fit_command(spec)
        if dry_run:
            print(f'E110 {spec.name}{" (" + why + ")" if why else ""}: ' + ' '.join(cmd), flush=True)
            continue
        if why:
            print(f'E110 {spec.name}: not yet, {why}', flush=True)
            waiting.append(spec.name)
            continue
        if not claim(spec.name, lane, reclaim):
            print(f'E110 {spec.name}: claimed by another process, skipped', flush=True)
            continue
        try:
            E108._still_clean(cc)
            print(f'E110 {spec.name}: ' + ' '.join(cmd), flush=True)
            rec = launch(cmd, fit_results_dir(spec.name), local_dir('logs') / f'{spec.name}.log', cc, spec.name)
            if rec['returncode']:
                raise RuntimeError(f"E110 {spec.name}: E10_joint_fit.py exited {rec['returncode']} (log {rec['log']})")
            r = fit_row(spec.name, tag)
            if r is None:
                raise RuntimeError(f"E110 {spec.name}: no row {tag} in its fits.csv (log {rec['log']})")
        except BaseException:
            release(spec.name)
            raise
        row = dict(name=spec.name, part=spec.part, wiring=spec.wiring, k=spec.k, kckc_fixed=spec.kckc_fixed,
                   maxfev=spec.maxfev, start=':'.join(spec.start), tag=tag, folder=f'{FIT_EXPERIMENT_PREFIX}/{spec.name}',
                   host=_host(), lane=lane, seconds=rec['seconds'], returncode=rec['returncode'], command=' '.join(cmd),
                   **params_of_row(r), loss=float(r['loss']), **{CHOOSE_COL: float(r[CHOOSE_COL])}, n_evals=int(r['n_evals']),
                   **cc)
        E83._append_row(path, row)
        print(f"E110 {spec.name}: loss {float(r['loss']):.4f}, {CHOOSE_COL} {float(r[CHOOSE_COL]):.4f}, KC-KC "
              f"{float(r['kc_kc_inh']):.3g}, {rec['seconds'] / 60:.0f} min", flush=True)
        done.append(spec.name)
    return done, waiting


# ---------------------------------------------------------------------------
# readings

def _dist(o: Dict[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (`E83.pair_distances`; 'cal' with K set on the run's APL-on counts), cached in this process
    and in outputs/E110_no_kckc/distances_cache.pkl by the sha256 of the run's spike counts and the readout."""
    import pickle
    from experiments.E77_input_screen import frame_sha256
    path = local_dir() / 'distances_cache.pkl'
    if 'dist_store' not in _CACHE:
        try:
            _CACHE['dist_store'] = pickle.loads(path.read_bytes()) if path.exists() else {}
        except Exception:
            _CACHE['dist_store'] = {}
    store = _CACHE['dist_store']
    key = (frame_sha256(o['on'][1]), frame_sha256(o['off'][1]), readout)
    if key not in store:
        store[key] = E83.pair_distances(o, readout)
        _CACHE['dist_dirty'] = True
    return store[key]


def _save_dist_cache() -> None:
    import pickle
    if _CACHE.get('dist_dirty'):
        try:
            (local_dir() / 'distances_cache.pkl').write_bytes(pickle.dumps(_CACHE['dist_store']))
        except OSError:
            pass
        _CACHE['dist_dirty'] = False


def scores_of(d: Dict[str, pd.Series]) -> dict:
    """The readings of one unit's distances ({'on', 'off'}): E108's `unit_scores` (against all recordings and halves A to
    D), the two ρ's by `E106.two_rhos` ('rho_on', 'rho_off'; all recordings) and the natural-source alignment's
    difference ('vcf_diff', APL on minus off)."""
    from experiments import E106_rho_gap_parts as E106
    row = E108.unit_scores(d)
    t = E106.two_rhos(E83._cs(d['on']), E83._cs(d['off']))
    row.update(rho_on=float(t['on']), rho_off=float(t['off']), vcf_diff=row['vcf_on'] - row['vcf_off'])
    return row


RUN_COLS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv', 'max_count_on',
            'max_count_off', 'silent_on', 'silent_off', 'n_kc_ge20_on', 'n_kc_ge20_off', 'fixed_thr', 'wAPLKC_scale',
            'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale')


def family_fit_row(w: str) -> pd.Series:
    folder, tag = E83.members()[w]
    return pd.read_csv(H.results_dir(folder) / 'fits.csv').set_index('tag').loc[tag]


def collect_runs() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """Every run there is: ({(kind, wiring): outputs}, runs table). Kinds: 'family' (the saved fits), 'A', 'B1' (part B's
    first start), 'B2' (its second), 'B' (the counted one), 'C', 'D', 'E<k>'. The table has one row per run: kind,
    wiring, tag, folder, the absolute parameters, the fit's loss and CHOOSE_COL, n_evals, realism (`E83.realism`),
    run_extras, the bands and 'flagged'."""
    outs: Dict[Tuple[str, str], Dict[str, tuple]] = {}
    rows = []

    def add(kind, w, o, **meta):
        outs[(kind, w)] = o
        rows.append(dict(kind=kind, wiring=w, **meta, **E83.realism(o), **run_extras(o)))

    for w in WIRINGS:
        F = E83.saved_fit(w)
        fr = family_fit_row(w)
        add('family', w, F['saved'], tag=F['tag'], folder=F['folder'], **F['p'], loss=float(fr['loss']),
            **{CHOOSE_COL: float(fr[CHOOSE_COL])}, n_evals=int(fr['n_evals']))
    ra = runs_table_of('A').set_index('tag') if len(runs_table_of('A')) else None
    for w in WIRINGS:
        o = load_a(w)
        if o is not None and ra is not None and f'A_{w}' in ra.index:
            r = ra.loc[f'A_{w}']
            add('A', w, o, tag=f'A_{w}', folder=f'{EXPERIMENT}/A', **{k: float(r[k]) for k in PARAM_KEYS},
                thr_factor=float(r['thr_factor']), apl_factor=float(r['apl_factor']))
    specs = fit_specs('B') + fit_specs('C') + fit_specs('E') + fit_specs('E-rest')
    try:
        specs += fit_specs('D')
    except Exception:
        pass
    for s in specs:
        tag = fit_tag(s)
        r = fit_row(s.name, tag)
        if r is None:
            continue
        o = load_fit_outputs(s.name, tag)
        if o is None:
            print(f'E110 score: {s.name} has a fits.csv row but no outputs here or in the archive', flush=True)
            continue
        kind = {'B': 'B1', 'B2': 'B2', 'C': 'C', 'D': 'D'}.get(s.part, f'E{float(s.k):g}')
        add(kind, s.wiring, o, tag=tag, folder=f'{FIT_EXPERIMENT_PREFIX}/{s.name}', name=s.name, **params_of_row(r),
            loss=float(r['loss']), **{CHOOSE_COL: float(r[CHOOSE_COL])}, n_evals=int(r['n_evals']))
    for w in WIRINGS:
        b = counted_b(w)
        if b is not None and (('B1' if b['name'].startswith('B_') else 'B2'), w) in outs:
            src = 'B1' if b['name'].startswith('B_') else 'B2'
            o = outs[(src, w)]
            meta = next(dict(r) for r in rows if r['kind'] == src and r['wiring'] == w)
            meta = {k: v for k, v in meta.items() if k not in ('kind', 'wiring')}
            meta.update(counted_from=src, n_starts=b['n_starts'])
            outs[('B', w)] = o
            rows.append(dict(kind='B', wiring=w, **meta))
    runs = pd.DataFrame(rows)
    bands = bands_table(runs)
    runs = pd.concat([runs.reset_index(drop=True), bands.reset_index(drop=True)], axis=1)
    runs['flagged'] = [E108._flagged(x) for x in runs['realism_flag']]
    return outs, runs


def bands_of(rate_on: float, rate_off: float, rate_rho_on: float, opk_on_tv: float, opk_off_tv: float,
             opk_bands: Optional[Tuple[float, float]] = None) -> Dict[str, bool]:
    """Inside (True) or outside each of the bands of kcq.apl_effect.realism_flag (units 'noise'), from a run's numbers."""
    from kcq import apl_effect as A
    ob = (A.opk_band('on', 'noise'), A.opk_band('off', 'noise')) if opk_bands is None else opk_bands
    ratio = rate_off / rate_on if rate_on > 0 else np.inf
    return dict(band_rate_on=bool(A.RATE_ON_BAND[0] <= rate_on <= A.RATE_ON_BAND[1]),
                band_ratio=bool(A.RATE_RATIO_BAND[0] <= ratio <= A.RATE_RATIO_BAND[1]),
                band_rate_rho_on=bool(not rate_rho_on < A.RATE_RHO_MIN),
                band_opk_on=bool(not opk_on_tv > ob[0]), band_opk_off=bool(not opk_off_tv > ob[1]))


def bands_table(runs: pd.DataFrame, opk_bands: Optional[Tuple[float, float]] = None) -> pd.DataFrame:
    """`bands_of` per row of a runs table (rate_on, rate_off, rate_rho_on, opk_on_tv, opk_off_tv)."""
    return pd.DataFrame([bands_of(float(r.rate_on), float(r.rate_off), float(r.rate_rho_on), float(r.opk_on_tv),
                                  float(r.opk_off_tv), opk_bands) for r in runs.itertuples()], index=runs.index)


def readings_table(outs: Mapping[Tuple[str, str], Dict[str, tuple]], runs: pd.DataFrame) -> pd.DataFrame:
    """Per kind, unit (each wiring with a run, and each population of UNITS whose wirings all have one) and readout:
    `scores_of`, the runs' means of RUN_COLS ('run_<col>'), 'n_flagged', and the change from the family's same unit and
    readout ('d_<col>')."""
    rows = []
    kinds = list(dict.fromkeys(k for k, _ in outs))
    for kind in kinds:
        have = [w for w in WIRINGS if (kind, w) in outs]
        sub = runs[runs['kind'] == kind].set_index('wiring')
        units = [(w, (w,)) for w in have] + [(u, ws) for u, ws in UNITS.items() if all(w in have for w in ws)]
        for ro in READOUTS:
            dd = {w: _dist(outs[(kind, w)], ro) for w in have}
            for unit, ws in units:
                single = len(ws) == 1 and unit in dd
                d = dd[unit] if single else E83.population([dd[w] for w in ws])
                s = sub.loc[list(ws)]
                rows.append(dict(kind=kind, unit=unit, unit_kind='wiring' if single else 'population', n_wirings=len(ws),
                                 readout=ro, n_flagged=int(s['flagged'].sum()),
                                 **{f'run_{c}': float(s[c].astype(float).mean()) for c in RUN_COLS if c in s.columns},
                                 **scores_of(d)))
    _save_dist_cache()
    rd = pd.DataFrame(rows)
    base = rd[rd['kind'] == 'family'].set_index(['unit', 'readout'])
    cols = [c for c in rd.columns if c not in ('kind', 'unit', 'unit_kind', 'n_wirings', 'readout')]
    for c in cols:
        rd[f'd_{c}'] = [v - base[c].get((u, r), np.nan) for u, r, v in zip(rd['unit'], rd['readout'], rd[c])]
    return rd


def _v(rd: pd.DataFrame, kind: str, unit: str, col: str, readout: str) -> float:
    r = rd[(rd['kind'] == kind) & (rd['unit'] == unit) & (rd['readout'] == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


def changes_table(rd: pd.DataFrame, kind: str = 'B', ref: str = 'family', units: Sequence[str] = WIRINGS,
                  pop: str = 'all8', measures: Sequence[Measure] = MEASURES) -> pd.DataFrame:
    """The verdict's measures, `kind` minus `ref`: one row per measure (index its key), columns 'pop' (the change of the
    population `pop`), 'mean_paired' (the mean over `units` of each wiring's paired change) and one per wiring."""
    rows = {}
    for m in measures:
        per = {w: _v(rd, kind, w, m.col, m.readout) - _v(rd, ref, w, m.col, m.readout) for w in units}
        rows[m.key] = dict(pop=_v(rd, kind, pop, m.col, m.readout) - _v(rd, ref, pop, m.col, m.readout),
                           mean_paired=float(np.mean(list(per.values()))), **per)
    return pd.DataFrame.from_dict(rows, orient='index')


def pending_table(rd: pd.DataFrame, kind: str = 'B', ref: str = 'family', units: Sequence[str] = WIRINGS,
                  pop: str = 'all8') -> pd.DataFrame:
    """The items reported apart (PENDING), on raw counts, `kind` minus `ref`, as `changes_table`, with 'fails' (the
    population's change lower than minus its slack) and 'slack'."""
    ms = [Measure(k, lab, col, 'raw', -1, sl, ('pending',)) for k, lab, col, sl in PENDING]
    t = changes_table(rd, kind, ref, units, pop, ms)
    t['slack'] = [m.slack for m in ms]
    t['label'] = [m.label for m in ms]
    t['fails'] = t['pop'] < -t['slack']
    return t


# ---------------------------------------------------------------------------
# the verdict

def worse_by(m: Measure, d: float) -> bool:
    """Whether a change `d` is worse than the measure's margin (NaN is not)."""
    if not np.isfinite(d):
        return False
    return bool(d > m.slack) if m.worse > 0 else bool(d < -m.slack)


def _gain(key: str, d: float) -> bool:
    if not np.isfinite(d):
        return False
    if key == 'pattern_cal':
        return bool(d <= -BETTER_PATTERN)
    return bool(d >= BETTER_RHO)


def wirings_needed(n: int) -> int:
    """N_WIRINGS of 8, the same share (rounded up) of another number of wirings (E108's rule)."""
    return E108.wirings_needed(n)


def verdict(ch: pd.DataFrame, bands: pd.DataFrame, c_flat: Optional[bool] = None) -> dict:
    """The verdict of the entry (the module docstring) from `changes_table`'s table `ch` (index the MEASURES' keys;
    columns 'pop', 'mean_paired' and one per wiring) and the runs' bands `bands` (index the wirings, columns the BANDS'
    keys, True inside). `c_flat`: part C's return fit flat (None: part C not complete).

    Returns dict: 'cls' (one of CLASSES), 'reason', 'not_needed_fails' / 'needed' (key -> dict(label, pop, mean_paired,
    n_wirings_worse) or, for bands, dict(label, n_runs_outside)), 'gains' (key -> dict(label, pop, n_wirings)), 'groups'
    (group name -> ['<label> <size>']), 'per_measure', 'n_needed', 'c_flat', 'c_pending'."""
    wirings = list(bands.index)
    n_need = wirings_needed(len(wirings))
    nn_fail, needed, per = {}, {}, {}
    for m in MEASURES:
        r = ch.loc[m.key]
        n_w = int(sum(worse_by(m, float(r[w])) for w in wirings))
        pw, mw = worse_by(m, float(r['pop'])), worse_by(m, float(r['mean_paired']))
        per[m.key] = dict(label=m.label, pop=float(r['pop']), mean_paired=float(r['mean_paired']), n_wirings_worse=n_w,
                          pop_worse=pw, mean_paired_worse=mw, slack=m.slack)
        # the margin holds on the population, on the mean paired difference and on n_need of the wirings
        if pw or mw or len(wirings) - n_w < n_need:
            nn_fail[m.key] = dict(label=m.label, pop=float(r['pop']), mean_paired=float(r['mean_paired']), n_wirings_worse=n_w)
        if pw and n_w >= n_need:
            needed[m.key] = dict(label=m.label, pop=float(r['pop']), mean_paired=float(r['mean_paired']), n_wirings_worse=n_w)
    for b in BANDS:
        n_out = int((~bands[b.key].astype(bool)).sum())
        per[b.key] = dict(label=b.label, n_runs_outside=n_out)
        if n_out:
            nn_fail[b.key] = dict(label=b.label, n_runs_outside=n_out)
        if n_out >= n_need:
            needed[b.key] = dict(label=b.label, n_runs_outside=n_out)
    gains = {}
    for key in ('pattern_cal', 'rho_on', 'rho_off'):
        r = ch.loc[key]
        n_w = int(sum(_gain(key, float(r[w])) for w in wirings))
        if _gain(key, float(r['pop'])) and n_w >= n_need:
            gains[key] = dict(label=MEASURE[key].label, pop=float(r['pop']), n_wirings=n_w)
    groups: Dict[str, List[str]] = {}
    src = needed if needed else nn_fail
    for key, v in src.items():
        gs = (MEASURE[key].groups if key in MEASURE else next(b.groups for b in BANDS if b.key == key))
        size = (f"{v['pop']:+.4f} (population), {v['mean_paired']:+.4f} (mean paired), {v['n_wirings_worse']} of "
                f"{len(wirings)} wirings" if 'pop' in v else f"{v['n_runs_outside']} of {len(wirings)} runs outside")
        for g in gs:
            groups.setdefault(GROUPS[g], []).append(f"{v['label']}: {size}")
    if c_flat:
        cls, why = 'not determined', "part C's return fit ends within 0.002 of both the fit and part B"
    elif gains and not nn_fail:
        cls, why = 'better without', 'a gain on the population and on enough wirings, every "not needed" condition holding'
    elif not nn_fail:
        cls, why = 'not needed', (f'every margin holds on the population, the mean paired difference and at least {n_need} of '
                                  f'{len(wirings)} wirings; every run inside the bands')
    elif needed:
        cls, why = 'needed', f'a margin fails on the population and on at least {n_need} of {len(wirings)} wirings'
    else:
        cls, why = 'not determined', f'a margin fails, but not on the population and on {n_need} of {len(wirings)} wirings'
    return dict(cls=cls, reason=why, not_needed_fails=nn_fail, needed=needed, gains=gains, groups=groups,
                per_measure=per, n_needed=n_need, c_flat=c_flat, c_pending=c_flat is None)


def c_reading(rows: pd.DataFrame) -> dict:
    """Part C's return fits: `rows` one per wiring of C_WIRINGS with 'wiring', 'c_kckc' (C's fitted KC-KC), 'c_loss' (C's
    CHOOSE_COL), 'fit_loss' (the family fit's) and 'b_loss' (part B's counted fit's); NaN rows for fits not done. Per
    wiring 'returns' (KC-KC above C_RETURN_KCKC and within C_WITHIN of the fit), 'flat' (within C_WITHIN of both the fit
    and part B), 'more_evaluations' (KC-KC at most C_RETURN_KCKC and lower than part B by more than C_WITHIN); 'complete',
    the counts and 'c_flat' (flat on `wirings_needed` of the wirings, None until complete)."""
    per = []
    for r in rows.itertuples():
        ok = all(np.isfinite(float(x)) for x in (r.c_kckc, r.c_loss, r.fit_loss, r.b_loss))
        d_fit, d_b = float(r.c_loss) - float(r.fit_loss), float(r.c_loss) - float(r.b_loss)
        per.append(dict(wiring=r.wiring, done=ok, c_kckc=float(r.c_kckc), d_fit=d_fit, d_b=d_b,
                        returns=bool(ok and r.c_kckc > C_RETURN_KCKC and abs(d_fit) <= C_WITHIN),
                        flat=bool(ok and abs(d_fit) <= C_WITHIN and abs(d_b) <= C_WITHIN),
                        more_evaluations=bool(ok and r.c_kckc <= C_RETURN_KCKC and d_b < -C_WITHIN)))
    p = pd.DataFrame(per)
    complete = bool(len(p) == len(C_WIRINGS) and p['done'].all())
    need = wirings_needed(len(C_WIRINGS))
    return dict(per=p, complete=complete, n_returns=int(p['returns'].sum()), n_flat=int(p['flat'].sum()),
                n_more_evaluations=int(p['more_evaluations'].sum()), n_needed=need,
                c_flat=(bool(p['flat'].sum() >= need) if complete else None))


def e_alternatives(pattern: Mapping[float, float], pattern_k2: float) -> List[float]:
    """The exponents whose pattern part (the population of E_WIRINGS under the readout) is lower than exponent 2's (part
    B's counted runs on the same wirings) by E_ALT_GAIN or more."""
    return [k for k, v in pattern.items() if np.isfinite(v) and np.isfinite(pattern_k2) and v <= pattern_k2 - E_ALT_GAIN]


# ---------------------------------------------------------------------------
# part score

def _fmt(x, spec='.3f') -> str:
    try:
        return format(float(x), spec) if np.isfinite(float(x)) else '-'
    except (TypeError, ValueError):
        return str(x)


def reference_check(rd: pd.DataFrame) -> dict:
    """The family's population (raw counts) against E108's I0 (results/E108_paa_extrapolation/readings.csv) and the
    entry's 0.690 / 0.756."""
    ours = {c: _v(rd, 'family', 'all8', c, 'raw') for c in ('rho_on', 'rho_off', 'rho_on_all', 'rho_off_all', 'all_pattern',
                                                             'mean_effect')}
    ours.update({f'{c}_cal': _v(rd, 'family', 'all8', c, 'cal') for c in ('all_pattern', 'mean_effect', 'all_loss')})
    out = dict(ours=ours, entry=dict(ENTRY_REFERENCE))
    p = H.results_dir('E108_paa_extrapolation') / 'readings.csv'
    if p.exists():
        r = pd.read_csv(p)
        r = r[(r['input'] == 'I0') & (r['unit'] == 'all8')].set_index('readout')
        e108 = dict(rho_on=float(r.loc['raw', 'rho_on_all']), rho_off=float(r.loc['raw', 'rho_off_all']),
                    all_pattern=float(r.loc['raw', 'all_pattern']), mean_effect=float(r.loc['raw', 'mean_effect']),
                    all_pattern_cal=float(r.loc['cal', 'all_pattern']), mean_effect_cal=float(r.loc['cal', 'mean_effect']),
                    all_loss_cal=float(r.loc['cal', 'all_loss']))
        out['e108_I0'] = e108
        out['diff_vs_e108'] = {k: ours.get(k, np.nan) - v for k, v in e108.items()}
    out['diff_vs_entry'] = {k: round(ours[k], 3) - v for k, v in ENTRY_REFERENCE.items()}
    out['two_rhos_vs_unit_scores'] = {c: ours[c] - ours[f'{c}_all'] for c in ('rho_on', 'rho_off')}
    return out


def start_check(runs: pd.DataFrame) -> pd.DataFrame:
    """Each part-B / part-E fit started at part A's held point: its first evaluation's response rates (history row 0)
    against part A's run of the wiring (equal when the fit's simulator at the start is part A's)."""
    rows = []
    a = runs[runs['kind'] == 'A'].set_index('wiring')
    for s in fit_specs('B') + fit_specs('E') + fit_specs('E-rest'):
        if s.start[0] != 'params' or s.wiring not in a.index:
            continue
        hp = fit_results_dir(s.name) / f'history_{fit_tag(s)}.csv'
        if not hp.exists():
            continue
        h0 = pd.read_csv(hp).iloc[0]
        ra = a.loc[s.wiring]
        # rates to 1e-12 (the history's and the realism's means sum in different orders); the start's parameters exactly
        rows.append(dict(name=s.name, k=s.k, rate_on_first=float(h0['rate_on']), rate_off_first=float(h0['rate_off']),
                         rate_on_A=float(ra['rate_on']), rate_off_A=float(ra['rate_off']),
                         thr_first=float(h0['fixed_thr']), thr_A=float(ra['fixed_thr']),
                         params_same=bool(all(float(h0[p]) == float(ra[p]) for p in PARAM_KEYS)),
                         same=bool(s.k == K and abs(float(h0['rate_on']) - float(ra['rate_on'])) < 1e-12
                                   and abs(float(h0['rate_off']) - float(ra['rate_off'])) < 1e-12)))
    return pd.DataFrame(rows)


def probe_check(runs: pd.DataFrame) -> pd.DataFrame:
    """Part A's held points beside today's uncommitted run for Yang (results/yang_family_no_kckc/held_<w>.json), where it
    has the wiring; chooses nothing."""
    rows = []
    d = H.shared_dir() / 'results' / 'yang_family_no_kckc'
    a = runs[runs['kind'] == 'A'].set_index('wiring')
    for w in WIRINGS:
        j = _read_json(d / f'held_{w}.json')
        if not j or w not in a.index:
            continue
        r = a.loc[w]
        rows.append(dict(wiring=w, thr_A=float(r['fixed_thr']), thr_probe=float(j['params_on']['fixed_thr']),
                         aplkc_A=float(r['wAPLKC_scale']), aplkc_probe=float(j['params_on']['wAPLKC_scale']),
                         rate_on_A=float(r['rate_on']), rate_on_probe=float(j['rate_on']), rate_off_A=float(r['rate_off']),
                         rate_off_probe=float(j['rate_off']),
                         same_params=bool(all(float(r[k]) == float(j['params_on'][k]) for k in PARAM_KEYS))))
    return pd.DataFrame(rows)


def c_rows(runs: pd.DataFrame) -> pd.DataFrame:
    """`c_reading`'s input from the runs table."""
    rows = []
    for w in C_WIRINGS:
        c = runs[(runs['kind'] == 'C') & (runs['wiring'] == w)]
        b = runs[(runs['kind'] == 'B') & (runs['wiring'] == w)]
        f = runs[(runs['kind'] == 'family') & (runs['wiring'] == w)]
        rows.append(dict(wiring=w, c_kckc=float(c['wKCKC_scale'].iloc[0]) if len(c) else np.nan,
                         c_loss=float(c[CHOOSE_COL].iloc[0]) if len(c) else np.nan,
                         fit_loss=float(f[CHOOSE_COL].iloc[0]), b_loss=float(b[CHOOSE_COL].iloc[0]) if len(b) else np.nan))
    return pd.DataFrame(rows)


def e_reading(rd: pd.DataFrame) -> Tuple[pd.DataFrame, List[float]]:
    """Part E on the population of E_WIRINGS (S4) under the readout and on raw counts, beside part B's counted runs (exponent
    2) and the family; and the alternatives (`e_alternatives`)."""
    rows = []
    for kind in ['family', 'B'] + [f'E{float(k):g}' for k in E_EXPONENTS]:
        for unit in ('S4', 'all8'):
            if not len(rd[(rd['kind'] == kind) & (rd['unit'] == unit)]):
                continue
            rows.append(dict(kind=kind, unit=unit, **{f'{c}_{ro}': _v(rd, kind, unit, c, ro) for c in
                             ('all_pattern', 'all_loss', 'mean_effect', 'all_amplitude', 'all_shape', 'all_sd_ratio',
                              'rho_on', 'rho_off', 'n_flagged') for ro in READOUTS}))
    t = pd.DataFrame(rows)
    pk2 = _v(rd, 'B', 'S4', 'all_pattern', 'cal')
    pat = {float(k): _v(rd, f'E{float(k):g}', 'S4', 'all_pattern', 'cal') for k in E_EXPONENTS}
    return t, e_alternatives(pat, pk2)


COMPARISON_ROWS = (
    ("ρ on odor-odor distances, APL on (raw)", 'rho_on', 'raw'), ("ρ on odor-odor distances, APL off (raw)", 'rho_off', 'raw'),
    ('APL-effect loss (readout)', 'all_loss', 'cal'), ('pattern part (readout)', 'all_pattern', 'cal'),
    ('amplitude part (readout)', 'all_amplitude', 'cal'), ('shape part (readout)', 'all_shape', 'cal'),
    ('spread ratio (readout)', 'all_sd_ratio', 'cal'), ('mean APL effect (readout)', 'mean_effect', 'cal'),
    ('ρ on the per-pair APL effect (readout)', 'all_rho_effect', 'cal'),
    ('Pearson r on the per-pair APL effect (readout)', 'all_pearson_effect', 'cal'),
    ('APL-effect loss (raw)', 'all_loss', 'raw'), ('pattern part (raw)', 'all_pattern', 'raw'),
    ('amplitude part (raw)', 'all_amplitude', 'raw'), ('shape part (raw)', 'all_shape', 'raw'),
    ('spread ratio (raw)', 'all_sd_ratio', 'raw'), ('mean APL effect (raw)', 'mean_effect', 'raw'),
    ('ρ on the per-pair APL effect (raw)', 'all_rho_effect', 'raw'), ('Pearson r on the per-pair APL effect (raw)', 'all_pearson_effect', 'raw'),
    ('response rate, APL on', 'run_rate_on', 'raw'), ('response rate, APL off', 'run_rate_off', 'raw'),
    ('APL-off / APL-on response rate', 'run_ratio', 'raw'), ('response rate ρ, APL on', 'run_rate_rho_on', 'raw'),
    ('response rate ρ, APL off', 'run_rate_rho_off', 'raw'), ('odors per KC, APL on (distance from the flies)', 'run_opk_on_tv', 'raw'),
    ('odors per KC, APL off (distance from the flies)', 'run_opk_off_tv', 'raw'), ('silent KCs, APL on', 'run_silent_on', 'raw'),
    ('silent KCs, APL off', 'run_silent_off', 'raw'), ('KCs with 20 or more spikes to some odor, APL off', 'run_n_kc_ge20_off', 'raw'),
    ('largest spike count, APL off', 'run_max_count_off', 'raw'), ('the block (raw)', 'block', 'raw'),
    ('natural-source alignment, APL on (raw)', 'vcf_on', 'raw'), ('natural-source alignment, APL off (raw)', 'vcf_off', 'raw'),
    ('natural-source alignment, APL on minus off (raw)', 'vcf_diff', 'raw'),
    ('structure vs wild-type somata (raw)', 'struct_wt_on', 'raw'), ('somatic threshold (Hz)', 'run_fixed_thr', 'raw'),
    ('APL->KC scale', 'run_wAPLKC_scale', 'raw'), ('APL->bouton scale', 'run_wAPLPN_scale', 'raw'),
    ('KC-KC scale', 'run_wKCKC_scale', 'raw'), ('runs outside the bands', 'n_flagged', 'raw'))


def comparison_table(rd: pd.DataFrame, v: dict, cr: dict, et: pd.DataFrame, alts: List[float]) -> str:
    """comparison_table.md: the population of the eight per kind, then the verdict with its measures per wiring."""
    L = ['# E110: KC-KC inhibition held at zero in the best family', '',
         'Each column the mean matrix over the eight wirings (the family: E72 part B\'s fits; A: KC-KC 0 at fixed scales, both '
         'response rates held; B: KC-KC 0, the other three parameters refitted, the counted start; B1: its first start, at '
         'part A\'s held point). Run statistics are means over the eight runs.', '']
    kinds = [k for k in ('family', 'A', 'B', 'B1') if len(rd[(rd['kind'] == k) & (rd['unit'] == 'all8')])]
    L += ['| statistic | ' + ' | '.join(kinds) + ' |', '|---|' + '---|' * len(kinds)]
    for lab, col, ro in COMPARISON_ROWS:
        sp = '.0f' if col in ('n_flagged', 'run_n_kc_ge20_off', 'run_max_count_off') else ('.4g' if 'scale' in col else '.3f')
        L.append(f'| {lab} | ' + ' | '.join(_fmt(_v(rd, k, 'all8', col, ro), sp) for k in kinds) + ' |')
    L += ['', f"## The verdict (part B against the family): {v['cls']}", '', v['reason'] + '.',
          '' if not v['c_pending'] else "Part C pending: the verdict is provisional until part C's return fits are scored.", '']
    for g, items in v['groups'].items():
        L += [f'- {g}: ' + '; '.join(items)]
    if v['gains']:
        L += ['- gains: ' + '; '.join(f"{x['label']} {x['pop']:+.4f} ({x['n_wirings']} wirings)" for x in v['gains'].values())]
    L += ['', '| measure | margin | population | mean paired | wirings worse |', '|---|---|---|---|---|']
    for key, p in v['per_measure'].items():
        if 'pop' in p:
            L.append(f"| {p['label']} | {p['slack']:g} | {p['pop']:+.4f} | {p['mean_paired']:+.4f} | {p['n_wirings_worse']} |")
        else:
            L.append(f"| {p['label']} | inside | | | {p['n_runs_outside']} runs outside |")
    if cr:
        L += ['', f"## Part C (return fits): complete {cr['complete']}, returns {cr['n_returns']}, flat {cr['n_flat']}, "
              f"more evaluations {cr['n_more_evaluations']} of {len(C_WIRINGS)}", '']
    if len(et):
        L += ['', f'## Part E (breadth exponent, KC-KC 0): alternatives {alts}', '',
              '| kind | unit | pattern (readout) | pattern (raw) | mean effect (readout) | ρ on | ρ off | outside bands |',
              '|---|---|---|---|---|---|---|---|']
        for r in et.itertuples():
            L.append(f'| {r.kind} | {r.unit} | {_fmt(r.all_pattern_cal)} | {_fmt(r.all_pattern_raw)} | {_fmt(r.mean_effect_cal)} | '
                     f'{_fmt(r.rho_on_raw)} | {_fmt(r.rho_off_raw)} | {_fmt(r.n_flagged_raw, ".0f")} |')
    return '\n'.join(L) + '\n'


def score() -> dict:
    """Part score (no simulations; the module docstring)."""
    out = out_dir()
    outs, runs = collect_runs()
    runs.to_csv(out / 'runs.csv', index=False)
    rd = readings_table(outs, runs)
    rd.to_csv(out / 'readings.csv', index=False)
    ref = reference_check(rd)
    _write_json(out / 'reference_check.json', ref)
    start_check(runs).to_csv(out / 'start_check.csv', index=False)
    probe_check(runs).to_csv(out / 'probe_check.csv', index=False)
    E108.e106_yardsticks().to_csv(out / 'yardsticks_e106.csv', index=False)
    cr = c_reading(c_rows(runs))
    cr['per'].to_csv(out / 'c_reading.csv', index=False)
    et, alts = e_reading(rd)
    et.to_csv(out / 'e_reading.csv', index=False)
    res = dict(reference_check=ref, e_alternatives=alts, c={k: v for k, v in cr.items() if k != 'per'})
    have_b = all(('B', w) in outs for w in WIRINGS)
    v = None
    if have_b:
        ch = changes_table(rd, 'B')
        ch.to_csv(out / 'changes_B.csv')
        b = runs[runs['kind'] == 'B'].set_index('wiring').loc[list(WIRINGS), [x.key for x in BANDS]]
        b.to_csv(out / 'bands_B.csv')
        v = verdict(ch, b, cr['c_flat'])
        pt = pending_table(rd, 'B')
        pt.to_csv(out / 'pending_B.csv')
        v['pending'] = {k: dict(label=r['label'], pop=float(r['pop']), mean_paired=float(r['mean_paired']),
                                fails=bool(r['fails'])) for k, r in pt.iterrows()}
        res['verdict'] = v
        # part D's reading (graded cost), when run
        if all(('D', w) in outs for w in D_WIRINGS):
            dt = pd.concat({'D_vs_family': changes_table(rd, 'D', 'family', D_WIRINGS, 'S4'),
                            'B_vs_family': changes_table(rd, 'B', 'family', D_WIRINGS, 'S4')})
            dt.to_csv(out / 'd_reading.csv')
    (out / 'comparison_table.md').write_text(comparison_table(rd, v, cr, et, alts) if v else
                                              '# E110\n\nPart B incomplete: no verdict yet.\n')
    _write_json(out / 'verdict.json', dict(res, generated=_now(), **H.code_commits()))
    _write_json(out / 'provenance_score.json', dict(part='score', script='experiments/E110_no_kckc.py', generated=_now(),
                                                    host=_host(), **H.code_commits()))
    print(json.dumps(dict(verdict=None if v is None else {k: v[k] for k in ('cls', 'reason', 'groups', 'gains', 'c_pending')},
                          reference_check=ref.get('diff_vs_e108'), e_alternatives=alts), indent=1, default=str), flush=True)
    return res


# ---------------------------------------------------------------------------
# the archive's tables

def write_variants() -> pd.DataFrame:
    """variants.csv (one row per part-A run and fit, absolute parameters) and model_pair_distances_long.csv (raw spike
    counts, the harness's primary model RDM) for the model-output archive and kcq/apl_effect.py."""
    out = out_dir()
    outs, runs = collect_runs()
    rows, parts = [], []
    for r in runs[~runs['kind'].isin(['family', 'B'])].itertuples():
        o = outs[(r.kind, r.wiring)]
        kind = {'A': "KC-KC 0 at fixed scales, both response rates held at the fit's (E77_input_screen.hold_rates)",
                'B1': "KC-KC held at 0, the other three refitted from part A's held point",
                'B2': "KC-KC held at 0, the other three refitted from the fit's own parameters (200 evaluations)",
                'C': "the return fit: from part B's counted optimum with KC-KC started at 2 and free",
                'D': "half the fitted KC-KC held, the other three refitted from the fit's own parameters"}.get(
                    r.kind, f'breadth exponent {r.kind[1:]} with KC-KC held at 0, refitted from part A\'s held point')
        rows.append(dict(tag=r.tag, variant=r.tag, part=r.kind, base='bouton', k=K if not str(r.kind).startswith('E')
                         else float(r.kind[1:]), wiring=r.wiring, kind=f"C2's family's {r.wiring}: {kind}",
                         ref_experiment=E83.members()[r.wiring][0], ref_fit=E83.members()[r.wiring][1],
                         fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale, wKCAPL_scale=r.wKCAPL_scale,
                         wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale,
                         kc_kc_inh=r.wKCKC_scale, rate_on=r.rate_on, rate_off=r.rate_off, max_count_on=r.max_count_on,
                         max_count_off=r.max_count_off, realism_flag=r.realism_flag, subdir=('A' if r.kind == 'A' else r.folder.split('/', 1)[1])))
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# commands

def lane_line(part: str, wirings: Sequence[str], lane: str, threads: int) -> str:
    return (f'OMP_NUM_THREADS={threads} OMP_THREAD_LIMIT={threads} nice -n 10 ./venv/bin/python -u experiments/E110_no_kckc.py '
            f'--part {part} --wirings {" ".join(wirings)} --lane {lane}')


def print_commands() -> None:
    """Every fit's E10_joint_fit.py command (as `run_fits` launches it) and example lane lines."""
    for part in ('B', 'C', 'E', 'E-rest'):
        for s in fit_specs(part):
            print(f'# {s.name} -> {fit_tag(s)}')
            print(' '.join(['./venv/bin/python' if i == 0 else x for i, x in enumerate(fit_command(s))]))
    print('# part D (only on a "needed" verdict):')
    for s in fit_specs('D'):
        print(' '.join(['./venv/bin/python' if i == 0 else x for i, x in enumerate(fit_command(s))]))
    print('# lanes (each lane runs the parts in this order; the claims split the points):')
    for lane, t in (('atlas', 8), ('hal1', 16), ('hal2', 16)):
        for part in ('A', 'B', 'E', 'C'):
            print(lane_line(part, WIRINGS if part in ('A', 'B') else (E_WIRINGS if part == 'E' else C_WIRINGS), lane, t)
                  + ('' if part == 'A' else ' --wait'))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'C', 'D', 'E', 'E-rest', 'score', 'variants', 'commands',
                                                      'starts'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=list(WIRINGS), help="a subset of the part's wirings (a lane)")
    ap.add_argument('--lane', default=None, help="the lane's name (its tables runs_<A|fits>_<lane>.csv; default the host)")
    ap.add_argument('--reclaim', action='store_true', help='take over claims (a process that died left them)')
    ap.add_argument('--dry-run', action='store_true', help='print the fits\' commands, run nothing')
    ap.add_argument('--force', action='store_true', help='part D without a "needed" verdict, E-rest without an alternative')
    ap.add_argument('--rerun', action='store_true', help='part 0: run the checks again')
    ap.add_argument('--wait', action='store_true',
                    help=f'fits: retry those whose start is not saved yet every {WAIT_SECONDS} s (a part-A or part-B lane elsewhere)')
    a = ap.parse_args()
    if a.part == '0':
        part_0(a.rerun)
    elif a.part == 'A':
        part_a(a.wirings or WIRINGS, a.lane, a.reclaim)
    elif a.part in ('B', 'C', 'E'):
        run_fits(fit_specs(a.part, a.wirings), a.lane, a.reclaim, a.dry_run, a.wait)
    elif a.part in ('D', 'E-rest'):
        v = _read_json(out_dir() / 'verdict.json')
        alts = list(v.get('e_alternatives') or [])
        due = (v.get('verdict', {}).get('cls') == 'needed') if a.part == 'D' else bool(alts)
        if not due and not a.force:
            raise SystemExit(f'E110 {a.part}: not due (verdict.json: {v.get("verdict", {}).get("cls")}, alternatives '
                             f'{alts}); --force runs it anyway')
        specs = (fit_specs('D', a.wirings) if a.part == 'D'
                 else fit_specs('E-rest', a.wirings, tuple(alts) if alts else E_EXPONENTS))
        run_fits(specs, a.lane, a.reclaim, a.dry_run, a.wait)
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'commands':
        print_commands()
    elif a.part == 'starts':
        print(json.dumps(write_starts(), indent=1))


if __name__ == '__main__':
    main()
