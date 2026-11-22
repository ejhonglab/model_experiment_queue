#!/usr/bin/env python3
"""E130: the best family on the male CNS's two calyces: a second fly's wiring, APL synapses and claw-to-claw synapses
among KCs in place of the hemibrain's.

The entry is item E130 of experiment_priorities.md (proposed by hal 2026-11-22 13:29, planned by a Fable subagent,
claude-fable-5-2; moved in under the idle-host rule). "The family" is C2's breadth-2 APL-effect fit (E67) on E72 part
B's eight wirings, each at its own fit (E83's `saved_fit`: real, P1, P2, M1, P3, P4, M2, M3). The calyces: the hemibrain
with multibouton_claws 'argmax' (each claw on the bouton with most of its synapses, olfsysm's default since 60cddc3; the
family's is 'first', E3's BOUTON pin) and the male CNS's left and right calyces (olfsysm's pn2kc_connections
'malecns-left' / 'malecns-right', 60cddc3), each through E10_joint_fit.py's `--connectome` (`J.connectome_kws`).

Parts (CLI `--part`):
- step0 (no simulation): the three calyces side by side from the model's tables (`J.bouton_base`'s wiring of each
  connectome, olfsysm's `connectome_APL_weights` counts, the KC-KC matrices below): KCs, claws and claws per KC by type,
  boutons; APL's counted synapses (the share of claws with none, APL->KC per KC within a claw count as E96's `_within`,
  the boutons' share of APL's output and input synapses in the calyx); the claw-to-claw matrix among the modelled KCs
  (synapses, the share within type, each type's incoming amount, the share of KCs with none, the share of model KCs in
  the table); claws per KC by type at E97's claw thresholds 3, 5 and 10 (`E97.counted_claws`). step0_*.csv, step0.json,
  step0_table.md.
- check (simulations, before step 1): (i) the bases (`tune_bases`: E3's tuning on each calyx, `J.connectome_base`;
  results/.../bases.json) and the same tuning on 'hemibrain-first' against E3's saved parameters (reported, the family's
  base stays E3's); (ii) the identity: this script's simulator on 'hemibrain-first' (E10's default) at C2's parameters
  reproduces C2's saved APL-on and APL-off runs exactly (0 count and 0 response differences); (iii) the male wirings'
  claws, KCs and boutons equal olfsysm's README (data/malecns_2026-10-25: 10,077 claws, 1,864 KCs, 471 boutons left;
  10,008, 1,884, 434 right); (iv) E10_joint_fit.py itself with `--connectome hemibrain-first` (the default, passed) from
  C2's exact start: its first evaluation equals C2's history at its best evaluation and C2's fits.csv row exactly
  (E110's `exact_check`); (v) 5-evaluation smoke fits of E10_joint_fit.py on 'malecns-left' (with its KC-KC matrix) and
  on 'hemibrain-argmax', the loader refusing nothing. check.json.
- 1 (simulations, resumable): six runs at fixed scales: each calyx with C2's factors on the calyx's own base and its own
  KC-KC matrix at C2's KC-KC scale ('S1_<calyx>') and at 0 ('S1_<calyx>-kckc0'), both response rates re-held
  (`E77_input_screen.hold_rates`: the threshold bisected with APL off, then APL->KC with KC->APL and bouton->APL tied, with
  APL on; APL->bouton and KC-KC as set). runs_1_<host>.csv, starts.json (the held points), provenance_1_<host>.json; runs
  under outputs/E130_malecns_family/step1/ (and the archive's step1/, written by archive_model_outputs.py).
- fit (`--targets`; E10_joint_fit.py processes, one after another): E72's part-B command (`-k 2 --aplpn-max 3
  --rate-target ctrl_frac --loss apl --apl-readout both --tag _cap3`, 150 evaluations) with `--connectome <calyx>
  --base-params results/E130_malecns_family/bases.json:<calyx>` and, on the male calyces, `--kckc-matrix` of the calyx's
  matrix; the targets 'hemibrain-argmax', 'malecns-left', 'malecns-right' (KC-KC free; start C2's factors on the
  calyx's base, `--init-from`, as E72's part B) and 'malecns-left-kckc0', 'malecns-right-kckc0' (`--kckc-fixed 0`; start
  the step-1 held point with KC-KC 0, `--init-params`, as E110's part B). A second start (`_cap3_s2`) when the first ends
  with fit_apl_cal_rmse above SECOND_START_ABOVE (0.066; E72's rule): the other of the two starts (the step-1 held point
  for the free targets, C2's factors with KC-KC 0 for the others). Each start claimed by mkdir of
  results/E130_malecns_family/claims/<name>/ and skipped when done or claimed (`--reclaim` takes over a dead claim); its
  folder results/ and outputs/E130_malecns_family/fits/<name>/; a row per fit in runs_fits_<lane>.csv.
- score (no simulation): readings.csv (each run and population, raw and calibrated: E110's `scores_of`, i.e. E108's
  `unit_scores` against all recordings and halves A to D with the amplitude part, shape part, spread ratio, ρ and Pearson
  r on the per-pair APL effect, and the two ρ's; the ester and alcohol groups' pattern errors; E92's line), runs.csv
  (realism, the four bands, response rates by KC type, parameters, bounds), populations.csv (the three real calyces
  against the 21 populations of three of the eight that contain the hemibrain wiring), verdict.json (the classes),
  comparison_table.md, provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive.
- commands: the lanes' command lines and every fit's E10_joint_fit.py command.

Choices fixed here before any run (the entry leaves them open):
- "The calyx's own base": E3's bouton-base tuning (fit_mb_model at target_sparsity 0.05 with APL on, E3's BOUTON kwargs)
  run on the calyx (`J.connectome_base`), once, on atlas, written to bases.json; every host's fits read it
  (`--base-params`). The family's base ('hemibrain-first') stays E3's saved parameters.
- "C2's factors" at fixed scales (step 1): all five scales as C2's over E3's base times the calyx's base (the threshold,
  APL->KC, KC->APL, APL->bouton, bouton->APL; on 'hemibrain-first' this is C2 exactly), KC-KC at C2's absolute scale
  (26.43; every calyx's matrix has a mean row sum of 1, olfsysm's normalization of the hemibrain's) or 0. The fits'
  start is E10's `--init-from` (E72's part B: the threshold, APL and APL->bouton factors on the base, bouton->APL the
  base's), so a family fit starts as the seven re-draws started.
- Step 1's held rates: C2's APL-on and APL-off response counts as fractions of its KCs x odors, times the calyx's
  (rounded), so each calyx is held at the family's two response rates.
- The male KC-KC matrices (`malecns_kckc`): Pratyush's KCclaw-to-KCclaw table of 2026-10-25, every synapse whose two
  KCs are modelled KCs of the calyx, counted per (post, pre) KC pair (claws ignored), autapses dropped, normalized to a
  mean row sum of 1 over the modelled KCs (olfsysm's `connectome_KCKC_weights` for the hemibrain). No ROI rule (the
  official calyx ROI misses part of the left calyx; E97 and the model tables use none); the hemibrain's table keeps
  presynapses in CA(R), and step 0 reports the male synapses outside CA(L) / CA(R) / CentralBrain-unspecified.
- The second start's column and the counted start: fit_apl_cal_rmse (the fit's calibrated APL-effect loss), the lower
  counting (E72's `new_fit_members`, E110's CHOOSE_COL).
- The population of the three real calyces: the family's hemibrain fit (C2, 'first'), the male left and the male right
  fits (counted starts), against the 21 populations {real, X, Y} of the eight; beside it (no class) the same with the
  'argmax' refit in place of C2.
- At a bound (E124's rule, BOUND_TOL): the threshold within [1/4, 4] x the base's, APL->KC within [1/30, 30] x the
  start's, KC-KC at 300; APL->bouton at its cap (3 x the base's) reported only.

The classes (fixed beforehand; the calibrated pattern part against all recordings and each of the four halves A to D;
WITHIN 0.003; a run outside a band gets no class and is listed for a "Realism flag:" post):
- the new default: "unchanged" if the 'argmax' refit is within 0.003 of the family's hemibrain fit (C2) on every
  target and inside the bands, else "changed";
- a male calyx: "inside" if its fit lies within the eight wirings' range on every target (the range per target),
  otherwise "outside (worse side)", "outside (better side)" or "outside (both sides)";
- the population: "real calyces average better" if the three real calyces as one population are below the lowest of
  the 21 and at least 0.003 below their median on every target, "worse" the mirror, otherwise "as the re-draws";
- KC-KC: "replicated" if on both male calyces the fit with KC-KC at 0 is higher by REPLICATED (0.005) or more on every
  target, "not replicated" if within 0.003 on both on every target, otherwise "mixed"; "not testable" on a calyx where
  fewer than half of the modelled KCs appear in the KC-KC table.
The amplitude part, shape part, spread ratio and Pearson r go beside every pattern part (comparison_table.md).

Holds (CLAUDE.md, "Connectome data"): the hemibrain (olfsysm's tables) and Pratyush's male CNS (olfsysm's model tables
of hong-data 8dc9cea, made from his run 2026-10-25-17-31-06, and his KCclaw-to-KCclaw table in data/Male-CNS_Tom/); no
other connectome. No olfsysm change. Nothing is published.

Concentrations: the model's input for paa, euc and -aPine is the ORN data within 1 log10 unit of the flies'.

Usage (from the repository root; 8 threads under nice on atlas):
    ./venv/bin/python -u experiments/E130_malecns_family.py --part step0
    OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 nice -n 10 ./venv/bin/python -u experiments/E130_malecns_family.py --part check
    OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 nice -n 10 ./venv/bin/python -u experiments/E130_malecns_family.py --part 1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E130_malecns_family.py --part fit \\
        --targets malecns-left hemibrain-argmax --lane hal-1
    ./venv/bin/python -u experiments/E130_malecns_family.py --part score && \\
        ./venv/bin/python -u experiments/E130_malecns_family.py --part variants
"""
from __future__ import annotations

import argparse
import gc
import itertools
import json
import os
import shutil
import socket
import time
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E10_joint_fit as J
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E130_malecns_family'
SCRIPT = 'experiments/E130_malecns_family.py'
FIRST = J.DEFAULT_CONNECTOME                           # 'hemibrain-first', the family's
CALYCES = ('hemibrain-argmax', 'malecns-left', 'malecns-right')
MALE = ('malecns-left', 'malecns-right')
STEP0_CALYCES = (FIRST,) + CALYCES
SHORT = {FIRST: 'hemibrain (first)', 'hemibrain-argmax': 'hemibrain (argmax)', 'malecns-left': 'male CNS left',
         'malecns-right': 'male CNS right'}
C2 = E83.C2                                            # ('E67_apl_effect_fits', 'rate0.068_aplboth_k2_cap3')
K = 2.0
MAXFEV = 150
SECOND_START_ABOVE = 0.066
CHOOSE_COL = 'fit_apl_cal_rmse'
WITHIN = 0.003
REPLICATED = 0.005
POP_BELOW = 0.003
BOUND_TOL = 0.01
THR_BOUNDS = (0.25, 4.0)
APLKC_BOUNDS = (1 / 30, 30.0)
KCKC_MAX = 300.0
TARGETS_READ = ('all', 'A', 'B', 'C', 'D')            # E108's effect targets: all recordings, halves A to D
READOUTS = ('raw', 'cal')
TYPES = ('g', "a'b'", 'ab')
TYPE_COL = {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}
PARAM5 = ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale')
PARAM_KEYS = PARAM5 + ('wKCKC_scale',)
# olfsysm's README of its male CNS tables (pyolfsysm/data/malecns_2026-10-25/README.md)
README_COUNTS = {'malecns-left': dict(claws=10077, kcs=1864, boutons=471),
                 'malecns-right': dict(claws=10008, kcs=1884, boutons=434)}
MALECNS_KCKC_FILE = 'KCclaw-to-KCclaw_Connectivity_Male-CNS_2026-10-25-17-31-06.parquet'
CLAW_THRESHOLDS = (3, 5, 10)
S1_RUNS = tuple(f'{c}{s}' for c in CALYCES for s in ('', '-kckc0'))
FIT_TARGETS = ('malecns-left', 'malecns-left-kckc0', 'hemibrain-argmax', 'malecns-right', 'malecns-right-kckc0')
LANES = {'hal-1': ('malecns-left', 'hemibrain-argmax'), 'hal-2': ('malecns-left-kckc0',),
         'atlas-1': ('malecns-right',), 'atlas-2': ('malecns-right-kckc0',)}
FIT_PREFIX = f'{EXPERIMENT}/fits'
_CACHE: dict = {}

FitSpec = namedtuple('FitSpec', 'name target connectome kckc0 start suffix maxfev')
# start: 'factors' (C2's, --init-from) or 'held' (the step-1 held point of the target, --init-params)


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E130_malecns_family/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def calyx_of(name: str) -> str:
    """The connectome of a step-1 run or fit target ('malecns-left-kckc0' -> 'malecns-left')."""
    c = name[:-len('-kckc0')] if name.endswith('-kckc0') else name
    if c not in STEP0_CALYCES:
        raise ValueError(f'{name!r}: not a calyx or a calyx with -kckc0')
    return c


def is_kckc0(name: str) -> bool:
    return name.endswith('-kckc0')


def c2_params() -> Dict[str, float]:
    """C2's fitted absolute parameters (E72's `fitted_params`; wKCKC_scale its kc_kc_inh)."""
    from experiments import E72_wiring_population as E72
    return E72.fitted_params(*C2)


def e3_base() -> Dict[str, float]:
    from experiments.E12_breadth_apl import e3_bouton_params
    return e3_bouton_params()


# ---------------------------------------------------------------------------
# the bases and the wirings

def bases_path() -> Path:
    return out_dir() / 'bases.json'


def base_params_spec(c: str) -> str:
    """`--base-params` FILE:KEY of a calyx's base (FILE relative to the shared directory)."""
    return f'results/{EXPERIMENT}/bases.json:{c}'


def tune_bases(retune_first: bool = True) -> dict:
    """Each calyx's own base (E3's tuning on it, `J.connectome_base`, cached under outputs/.../bases/) in bases.json,
    with E3's parameters under 'hemibrain-first'; with `retune_first`, the same tuning on 'hemibrain-first' recorded
    beside as 'hemibrain-first-retuned' (its relative differences from E3's saved parameters), not used by any run. A
    calyx already in bases.json is kept (a base never changes once a run may have used it)."""
    from pyolfsysm.mb_model import fit_mb_model
    cur = _read_json(bases_path())
    orn = H.apl16_orn_deltas()
    changed = False
    if FIRST not in cur:
        cur[FIRST] = dict(e3_base(), source="E3's saved bouton base (e3_bouton_params; the family's)")
        changed = True
    for c in CALYCES:
        if c in cur:
            continue
        t0 = time.time()
        p = J.connectome_base(orn, c, local_dir('bases'))
        rec = _read_json(local_dir('bases') / f'{J.connectome_tag(c)}base_params.json')
        cur[c] = dict(p, source=f"E3's tuning on {c} (E10_joint_fit.connectome_base)", rate_on=rec.get('rate_on'),
                      tuned_on=_host(), tuned=_now(), seconds=time.time() - t0, commits=H.code_commits())
        changed = True
        _write_json(bases_path(), cur)
    if retune_first and 'hemibrain-first-retuned' not in cur:
        t0 = time.time()
        r_on, _, _, params = fit_mb_model(orn, target_sparsity=0.05, **J.base_kws(None))
        p = {k: float(np.mean(params[k])) for k in PARAM5}
        e3 = e3_base()
        cur['hemibrain-first-retuned'] = dict(
            {f'retuned_{k}': v for k, v in p.items()}, rel_diff={k: p[k] / e3[k] - 1 for k in PARAM5},
            rate_on=float(r_on.to_numpy(dtype=float).mean()), seconds=time.time() - t0, tuned_on=_host(), tuned=_now(),
            note="E3's tuning re-run on the family's connectome, against E3's saved base (a check; not used by any run)")
        changed = True
    if changed:
        _write_json(bases_path(), cur)
    return cur


def calyx_base(c: str) -> Dict[str, float]:
    """A calyx's base parameters from bases.json (E3's for 'hemibrain-first')."""
    if c == FIRST:
        return e3_base()
    b = _read_json(bases_path())
    if c not in b:
        raise RuntimeError(f'{c}: no base in {bases_path()} yet (run --part check first)')
    return {k: float(b[c][k]) for k in PARAM5}


def calyx_parts(c: str) -> dict:
    """A calyx's pieces, once per process: 'orn', 'base', 'wnorm', 'W' (claw x bouton wiring), 'breadth'
    (`J.bouton_base` with the connectome and its base; for 'hemibrain-first' E10's default, as E83's 'real')."""
    if c not in _CACHE.setdefault('parts', {}):
        orn = H.apl16_orn_deltas()
        if c == FIRST:
            initial, wnorm, W, breadth = J.bouton_base(orn, H.outputs_dir('E83_cross_wiring_fit'), None)
        else:
            initial, wnorm, W, breadth = J.bouton_base(orn, local_dir('base'), None, connectome=c,
                                                       base_params=base_params_spec(c))
        _CACHE['parts'][c] = dict(orn=orn, base=initial, wnorm=wnorm, W=W, breadth=breadth)
    return _CACHE['parts'][c]


def model_kc_ids(W: pd.DataFrame) -> pd.Index:
    """The model's KC IDs in its KC order (the wiring's claws' KC IDs, first appearance)."""
    return pd.Index(pd.unique(W.index.get_level_values(H.KC_ID_LEVEL)), name=H.KC_ID_LEVEL)


def kc_types_of_wiring(W: pd.DataFrame) -> pd.Series:
    """Each model KC's type (the wiring's kc_type level), indexed by KC ID in the model's order."""
    f = W.index.to_frame(index=False)[[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL]].drop_duplicates(H.KC_ID_LEVEL)
    return pd.Series(f[H.KC_TYPE_LEVEL].astype(str).to_numpy(), index=pd.Index(f[H.KC_ID_LEVEL].to_numpy(), name=H.KC_ID_LEVEL))


# ---------------------------------------------------------------------------
# the KC-KC matrices

def malecns_kckc_table() -> pd.DataFrame:
    """Pratyush's male CNS KCclaw-to-KCclaw synapses (one row per synapse; both sides), the columns used here."""
    if 'kckc_table' not in _CACHE:
        from experiments.E97_glomerulus_claw_shares import malecns_dir
        _CACHE['kckc_table'] = pd.read_parquet(malecns_dir() / MALECNS_KCKC_FILE,
                                               columns=['bodyId_pre', 'bodyId_post', 'roi_pre', 'roi_post'])
    return _CACHE['kckc_table']


def kckc_counts(table: pd.DataFrame, kc_ids: Sequence[int]) -> Tuple[pd.DataFrame, dict]:
    """KC x KC synapse counts (post x pre, `kc_ids` order on both axes) of the synapses in `table` whose two KCs are in
    `kc_ids`, autapses dropped; and what was kept (synapses kept, autapses dropped, roi_pre counts of the kept)."""
    ids = pd.Index(kc_ids)
    if ids.duplicated().any():
        raise ValueError('kc_ids must not have duplicates')
    t = table[table['bodyId_pre'].isin(ids) & table['bodyId_post'].isin(ids)]
    auto = t['bodyId_pre'] == t['bodyId_post']
    info = dict(n_synapses=int((~auto).sum()), n_autapses_dropped=int(auto.sum()))
    if 'roi_pre' in t.columns:
        info['roi_pre'] = {str(k): int(v) for k, v in t.loc[~auto, 'roi_pre'].value_counts().items()}
    t = t[~auto]
    c = t.groupby(['bodyId_post', 'bodyId_pre']).size()
    M = np.zeros((len(ids), len(ids)))
    M[ids.get_indexer(c.index.get_level_values(0)), ids.get_indexer(c.index.get_level_values(1))] = c.to_numpy()
    return (pd.DataFrame(M, index=pd.Index(ids, name=f'{H.KC_ID_LEVEL}_post'), columns=pd.Index(ids, name=f'{H.KC_ID_LEVEL}_pre')),
            info)


def normalize_kckc(counts: pd.DataFrame) -> pd.DataFrame:
    """Counts divided by their mean row sum (over all the KCs, those without input included), as olfsysm's
    `connectome_KCKC_weights(normalize=True)`."""
    total = float(counts.to_numpy().sum())
    if not total > 0:
        raise ValueError('no KC-KC synapses among these KCs')
    return counts / (total / len(counts))


def malecns_kckc(connectome: str, kc_ids: Sequence[int]) -> pd.DataFrame:
    """The male calyx's KC x KC matrix (post x pre) among the model's KCs `kc_ids`, from Pratyush's KCclaw-to-KCclaw
    table, normalized as olfsysm's `connectome_KCKC_weights` (mean row sum 1)."""
    if connectome not in MALE:
        raise ValueError(f'{connectome=} not in {MALE}')
    counts, _ = kckc_counts(malecns_kckc_table(), kc_ids)
    return normalize_kckc(counts)


def kckc_raw(c: str, W: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """A calyx's KC-KC synapse counts among the model KCs of its wiring `W` (post x pre, the model's KC order) and the
    kept-synapse record: the male CNS's from Pratyush's table, the hemibrain's from olfsysm's (`connectome_KCKC_weights`)."""
    ids = model_kc_ids(W)
    if c in MALE:
        return kckc_counts(malecns_kckc_table(), ids)
    from pyolfsysm.mb_model import connectome_KCKC_weights
    M = connectome_KCKC_weights(ids, normalize=False)
    return M, dict(n_synapses=int(M.to_numpy().sum()), n_autapses_dropped=0, roi_pre={'CA(R)': int(M.to_numpy().sum())})


def matrix_path(c: str) -> Path:
    """The male calyx's matrix as E10_joint_fit.py's --kckc-matrix reads it (relative to the repository root)."""
    return Path('outputs') / EXPERIMENT / 'matrices' / f'kckc_{c}.parquet'


def write_matrix(c: str) -> Tuple[Path, pd.DataFrame, str]:
    """The male calyx's normalized matrix, written to outputs/.../matrices/kckc_<c>.parquet (columns as str, as E123
    writes them) unless there with the same values; its sha256 (E77's `frame_sha256`) checked against matrices.json in
    results (written the first time, so every host's matrix is the same). (path relative to the root, matrix, sha256)."""
    from experiments.E77_input_screen import frame_sha256
    if c not in MALE:
        raise ValueError(f'{c}: only the male calyces have a matrix of their own')
    if c in _CACHE.setdefault('mat', {}):
        return _CACHE['mat'][c]
    M = malecns_kckc(c, model_kc_ids(calyx_parts(c)['W']))
    sha = frame_sha256(M)
    rec = _read_json(out_dir() / 'matrices.json')
    if c in rec and rec[c]['sha256'] != sha:
        raise RuntimeError(f"{c}: this host's KC-KC matrix differs from matrices.json's ({sha} vs {rec[c]['sha256']})")
    if c not in rec:
        rec[c] = dict(sha256=sha, n_kcs=len(M), written_on=_host(), written=_now(), source=MALECNS_KCKC_FILE,
                      rule='synapses between modelled KCs of the calyx, per (post, pre) KC pair, autapses dropped, '
                           'normalized to a mean row sum of 1')
        _write_json(out_dir() / 'matrices.json', rec)
    p = H.HERE / matrix_path(c)
    p.parent.mkdir(parents=True, exist_ok=True)
    ok = False
    if p.exists():
        old = pd.read_parquet(p)
        old.columns = pd.Index([int(x) for x in old.columns], name=old.columns.name)
        ok = frame_sha256(old) == sha
    if not ok:
        Wr = M.copy()
        Wr.columns = Wr.columns.astype(str)
        Wr.to_parquet(p)
    _CACHE['mat'][c] = (matrix_path(c), M, sha)
    return _CACHE['mat'][c]


# ---------------------------------------------------------------------------
# step 0

def apl_counts(c: str, W: pd.DataFrame) -> Dict[str, pd.Series]:
    """APL's counted synapses on the wiring's claws and boutons (olfsysm's `connectome_APL_weights`, not normalized):
    'apl_kc', 'kc_apl' per claw, 'apl_pn', 'pn_apl' per bouton."""
    from pyolfsysm.mb_model import connectome_APL_weights
    conn = 'hemibrain' if c.startswith('hemibrain') else c
    a = connectome_APL_weights(conn, prat_claws=True, prat_boutons=True, wPNKC=W)
    return dict(zip(('apl_kc', 'kc_apl', 'apl_pn', 'pn_apl'), a[:4]))


def step0_wiring(c: str, W: pd.DataFrame) -> pd.DataFrame:
    """KCs, claws, claws per KC (mean, median) by type ('all' first) and the calyx's boutons and glomeruli."""
    f = W.index.to_frame(index=False)
    per_kc = f.groupby(H.KC_ID_LEVEL).agg(n_claws=(H.KC_ID_LEVEL, 'size'), kc_type=(H.KC_TYPE_LEVEL, 'first'))
    rows = []
    for t in ('all',) + TYPES:
        s = per_kc if t == 'all' else per_kc[per_kc['kc_type'] == t]
        rows.append(dict(calyx=c, kc_type=t, n_kcs=int(len(s)), n_claws=int(s['n_claws'].sum()),
                         claws_per_kc_mean=float(s['n_claws'].mean()), claws_per_kc_median=float(s['n_claws'].median()),
                         n_boutons=int(W.shape[1]) if t == 'all' else np.nan,
                         n_glomeruli=int(W.columns.get_level_values('glomerulus').nunique()) if t == 'all' else np.nan))
    return pd.DataFrame(rows)


def step0_apl(c: str, W: pd.DataFrame, a: Mapping[str, pd.Series]) -> dict:
    """APL's counted synapses: the share of claws with no APL->claw (and no claw->APL) synapse, APL->KC per KC by type,
    the evenness of the per-KC APL->KC totals among KCs with the same number of claws (E96's `_within`: the CV over the
    KCs of each claw count 2 to 10 with 10 KCs or more; min, median, max), the boutons' share of APL's output and input
    synapses (bouton / (bouton + claw))."""
    from experiments.E96_offclaw_apl import _within
    kc = W.index.get_level_values(H.KC_ID_LEVEL)
    tot = pd.Series(a['apl_kc'].to_numpy(dtype=float), index=kc).groupby(level=0).sum()
    n = pd.Series(1, index=kc).groupby(level=0).sum().reindex(tot.index)
    typ = kc_types_of_wiring(W).reindex(tot.index)
    w = _within(n.to_numpy(dtype=float), tot.to_numpy(dtype=float))
    s_kc, s_pn = float(a['apl_kc'].sum()), float(a['apl_pn'].sum())
    i_kc, i_pn = float(a['kc_apl'].sum()), float(a['pn_apl'].sum())
    row = dict(calyx=c, claws_no_apl_kc=float((a['apl_kc'].to_numpy() == 0).mean()),
               claws_no_kc_apl=float((a['kc_apl'].to_numpy() == 0).mean()),
               kcs_no_apl_kc=float((tot.to_numpy() == 0).mean()),
               apl_kc_synapses=s_kc, apl_pn_synapses=s_pn, kc_apl_synapses=i_kc, pn_apl_synapses=i_pn,
               bouton_share_apl_output=s_pn / (s_pn + s_kc), bouton_share_apl_input=i_pn / (i_pn + i_kc),
               apl_kc_per_claw=float(a['apl_kc'].mean()), apl_kc_per_kc=float(tot.mean()),
               apl_kc_cv_kcs=float(tot.std() / tot.mean()), **{f'evenness_{k}': v for k, v in w.items()
                                                               if k.startswith('within_cv')})
    for t in TYPES:
        row[f'apl_kc_per_kc_{TYPE_COL[t]}'] = float(tot[typ == t].mean())
    return row


def step0_kckc(c: str, M: pd.DataFrame, info: Mapping, types: pd.Series) -> dict:
    """The claw-to-claw matrix among the modelled KCs (counts, post x pre): synapses, the share within type, each
    type's incoming amount (mean row sum per post type over the mean row sum of all KCs, i.e. in the normalized
    matrix's units; and synapses per KC), the share of KCs with no input, with no output, and of the model's KCs in the
    table at all (in or out)."""
    A = M.to_numpy(dtype=float)
    t = types.reindex(M.index).to_numpy()
    total = A.sum()
    same = (t[:, None] == t[None, :])
    rs, cs = A.sum(axis=1), A.sum(axis=0)
    row = dict(calyx=c, n_kcs=int(len(A)), synapses=float(total), within_type_share=float(A[same].sum() / total),
               kcs_no_input=float((rs == 0).mean()), kcs_no_output=float((cs == 0).mean()),
               kcs_in_table=float(((rs + cs) > 0).mean()), synapses_per_kc=float(total / len(A)),
               n_autapses_dropped=int(info.get('n_autapses_dropped', 0)),
               roi_pre=json.dumps(info.get('roi_pre', {}), sort_keys=True))
    mean_rs = total / len(A)
    for tp in TYPES:
        m = t == tp
        row[f'incoming_{TYPE_COL[tp]}'] = float(rs[m].mean() / mean_rs) if m.any() else np.nan
        row[f'synapses_in_per_kc_{TYPE_COL[tp]}'] = float(rs[m].mean()) if m.any() else np.nan
    row['testable'] = bool(row['kcs_in_table'] >= 0.5)
    return row


def gamma_thresholds() -> pd.DataFrame:
    """Claws per KC by class at E97's claw thresholds 3, 5 and 10 (`E97.counted_claws`, each brain's ROI rule), per
    mushroom body; the γ / αβ and γ / all ratios beside."""
    from experiments import E97_glomerulus_claw_shares as E97
    rows = []
    for t in CLAW_THRESHOLDS:
        for brain in ('hemibrain', 'malecns'):
            cl = E97.counted_claws(brain, syn_min=t)
            per = cl.groupby(['mb', 'kc_class', 'kc_id']).size().rename('n').reset_index()
            for mb, g in per.groupby('mb'):
                r = dict(threshold=t, mb=mb, claws_per_kc_all=float(g['n'].mean()), n_kcs_all=int(len(g)))
                for k in TYPES:
                    gg = g[g['kc_class'] == k]
                    r[f'claws_per_kc_{TYPE_COL[k]}'] = float(gg['n'].mean()) if len(gg) else np.nan
                    r[f'n_kcs_{TYPE_COL[k]}'] = int(len(gg))
                r['g_over_ab'] = r['claws_per_kc_g'] / r['claws_per_kc_ab']
                r['g_over_all'] = r['claws_per_kc_g'] / r['claws_per_kc_all']
                rows.append(r)
    return pd.DataFrame(rows)


def step0_wiring_only(c: str) -> pd.DataFrame:
    """A calyx's claw x bouton wiring without any simulation (olfsysm's `connectome_wPNKC` with the base's kwargs), the
    same as `J.bouton_base`'s (checked in part check)."""
    from pyolfsysm.mb_model import connectome_wPNKC
    if c not in _CACHE.setdefault('W0', {}):
        kw = J.base_kws(None, connectome=c)
        conn = kw.get('pn2kc_connections', 'hemibrain')
        _CACHE['W0'][c] = connectome_wPNKC(conn, prat_claws=True, prat_boutons=True,
                                           multibouton_claws=kw['multibouton_claws'])
    return _CACHE['W0'][c]


def step0() -> Dict[str, pd.DataFrame]:
    """The three calyces side by side (no simulation); the family's 'first' beside the hemibrain's 'argmax'."""
    out = out_dir()
    wr, ap, kk = [], [], []
    for c in STEP0_CALYCES:
        W = step0_wiring_only(c)
        wr.append(step0_wiring(c, W))
        ap.append(step0_apl(c, W, apl_counts(c, W)))
        M, info = kckc_raw(c, W)
        kk.append(step0_kckc(c, M, info, kc_types_of_wiring(W)))
    res = dict(wiring=pd.concat(wr, ignore_index=True), apl=pd.DataFrame(ap), kckc=pd.DataFrame(kk),
               gamma=gamma_thresholds())
    for k, df in res.items():
        df.to_csv(out / f'step0_{k}.csv', index=False)
    (out / 'step0_table.md').write_text(step0_table(res))
    _write_json(out / 'step0.json', dict(generated=_now(), host=_host(), script=SCRIPT, commits=H.code_commits(),
                                         testable={r['calyx']: r['testable'] for r in kk}))
    print(step0_table(res), flush=True)
    return res


def _f(x, sp: str = '.3f') -> str:
    try:
        return '' if x is None or (isinstance(x, float) and not np.isfinite(x)) else format(float(x), sp)
    except (TypeError, ValueError):
        return str(x)


def step0_table(res: Mapping[str, pd.DataFrame]) -> str:
    """Step 0's tables as markdown (statistic down the first column, the calyces across)."""
    cs = list(STEP0_CALYCES)
    wr = res['wiring'].set_index(['calyx', 'kc_type'])
    ap = res['apl'].set_index('calyx')
    kk = res['kckc'].set_index('calyx')
    L = ['| statistic | ' + ' | '.join(SHORT[c] for c in cs) + ' |', '|---|' + '---:|' * len(cs)]

    def line(lab, vals):
        L.append(f'| {lab} | ' + ' | '.join(vals) + ' |')
    for t in ('all',) + TYPES:
        line(f'KCs, {t}', [_f(wr.loc[(c, t), 'n_kcs'], '.0f') for c in cs])
        line(f'claws, {t}', [_f(wr.loc[(c, t), 'n_claws'], '.0f') for c in cs])
        line(f'claws per KC, {t}', [_f(wr.loc[(c, t), 'claws_per_kc_mean'], '.2f') for c in cs])
    line('boutons', [_f(wr.loc[(c, 'all'), 'n_boutons'], '.0f') for c in cs])
    line('glomeruli', [_f(wr.loc[(c, 'all'), 'n_glomeruli'], '.0f') for c in cs])
    for col, lab, sp in (('claws_no_apl_kc', 'claws with no APL->claw synapse', '.3f'),
                         ('claws_no_kc_apl', 'claws with no claw->APL synapse', '.3f'),
                         ('apl_kc_per_kc', 'APL->KC synapses per KC', '.2f'),
                         ('apl_kc_per_kc_g', 'APL->KC synapses per KC, g', '.2f'),
                         ('apl_kc_per_kc_apbp', "APL->KC synapses per KC, a'b'", '.2f'),
                         ('apl_kc_per_kc_ab', 'APL->KC synapses per KC, ab', '.2f'),
                         ('evenness_within_cv_median', 'CV of APL->KC per KC within a claw count (median)', '.3f'),
                         ('evenness_within_cv_min', 'same, lowest claw count group', '.3f'),
                         ('evenness_within_cv_max', 'same, highest', '.3f'),
                         ('bouton_share_apl_output', "boutons' share of APL's output synapses", '.3f'),
                         ('bouton_share_apl_input', "boutons' share of APL's input synapses", '.3f')):
        line(lab, [_f(ap.loc[c, col], sp) for c in cs])
    for col, lab, sp in (('synapses', 'KC-KC synapses among the modelled KCs', '.0f'),
                         ('within_type_share', 'KC-KC share within type', '.3f'),
                         ('incoming_g', 'KC-KC incoming amount, g (x mean)', '.3f'),
                         ('incoming_apbp', "KC-KC incoming amount, a'b' (x mean)", '.3f'),
                         ('incoming_ab', 'KC-KC incoming amount, ab (x mean)', '.3f'),
                         ('kcs_no_input', 'KCs with no KC-KC input', '.3f'),
                         ('kcs_no_output', 'KCs with no KC-KC output', '.3f'),
                         ('kcs_in_table', 'modelled KCs in the KC-KC table', '.3f'),
                         ('n_autapses_dropped', 'autapses dropped', '.0f')):
        line(lab, [_f(kk.loc[c, col], sp) for c in cs])
    g = res['gamma']
    L += ['', '| threshold | mushroom body | claws per KC, all | g | ab | a\'b\' | g / ab |', '|---:|---|---:|---:|---:|---:|---:|']
    for r in g.itertuples():
        L.append(f'| {r.threshold} | {r.mb} | {r.claws_per_kc_all:.2f} | {r.claws_per_kc_g:.2f} | {r.claws_per_kc_ab:.2f} | '
                 f'{r.claws_per_kc_apbp:.2f} | {r.g_over_ab:.2f} |')
    return '\n'.join(L) + '\n'


# ---------------------------------------------------------------------------
# the simulator and step 1

def build_simulator(c: str, params: Mapping[str, float]) -> Tuple:
    """C2's family on a calyx as a fixed-parameter simulator at `params` (E83's `build_simulator` with the connectome's
    kwargs, the calyx's wiring, breadth and normalized APL weights, and on the male calyces its own KC-KC matrix as
    `_wKCKC`); (simulate, info) with info's 'kws' the scalar model kwargs."""
    from pyolfsysm.mb_model import make_fixed_param_kc_simulator
    from experiments.E77_input_screen import _scalar_kws
    P = calyx_parts(c)
    p = dict(params)
    wnorm = P['wnorm']
    wAPLKC = J.breadth_weights(wnorm['wAPLKC'], P['breadth'], K)
    extra = {}
    if c in MALE:
        extra['_wKCKC'] = write_matrix(c)[1]
    kws = J.base_kws(None, connectome=c)
    sim, _ = make_fixed_param_kc_simulator(P['orn'], initial_params=p, **kws, **extra, _wPNKC=P['W'], _wAPLKC=wAPLKC,
                                           _wKCAPL=wnorm['wKCAPL'], _wAPLPN=wnorm['wAPLPN'], _wPNAPL=wnorm['wPNAPL'],
                                           kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU)
    return sim, dict(params=p, kws=_scalar_kws(dict(kws, kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU)),
                     connectome=c, multibouton_claws=kws['multibouton_claws'], breadth=P['breadth'], wAPLKC=wAPLKC,
                     W=P['W'], orn=P['orn'])


def transfer_params(base: Mapping[str, float], kckc: bool, c2: Optional[Mapping[str, float]] = None,
                    e3: Optional[Mapping[str, float]] = None) -> Dict[str, float]:
    """C2's factors on a base: each of the five scales C2's over E3's base times `base`'s; KC-KC C2's absolute scale
    (`kckc`) or 0."""
    c2 = c2_params() if c2 is None else c2
    e3 = e3_base() if e3 is None else e3
    p = {k: float(base[k]) * float(c2[k]) / float(e3[k]) for k in PARAM5}
    p['wKCKC_scale'] = float(c2['wKCKC_scale']) if kckc else 0.0
    return p


def held_targets(n_kcs: int, n_odors: int, c2_on: int, c2_off: int, c2_size: int) -> Tuple[int, int]:
    """(APL-on, APL-off) response counts on a calyx of `n_kcs` x `n_odors`: C2's counts as fractions of its size."""
    size = int(n_kcs) * int(n_odors)
    return int(round(c2_on / c2_size * size)), int(round(c2_off / c2_size * size))


def s1_point(name: str) -> Tuple[Dict[str, float], int, int]:
    """A step-1 run's start (C2's factors on the calyx's base, `transfer_params`) and its held counts (`held_targets`)."""
    c = calyx_of(name)
    F = E83.saved_fit('real')
    P = calyx_parts(c)
    p = transfer_params(P['base'], not is_kckc0(name))
    n_kcs = int(P['breadth'].shape[0])
    on, off = held_targets(n_kcs, P['orn'].shape[1], F['target_on'], F['target_off'], F['size'])
    return p, on, off


def runs_table() -> pd.DataFrame:
    fs = [pd.read_csv(p) for p in sorted(out_dir().glob('runs_1_*.csv'))]
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last') if fs else pd.DataFrame(columns=['tag'])


def save_run(tag: str, o: Dict[str, tuple]) -> None:
    d = local_dir('step1')
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(tag: str) -> Optional[Dict[str, tuple]]:
    for d in (H.HERE / 'outputs' / EXPERIMENT / 'step1', archive_dir() / 'step1'):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def type_rates(o: Dict[str, tuple]) -> dict:
    from experiments import E122_feedback_neuron as E122
    return E122.type_rates(o)


def starts_path() -> Path:
    return out_dir() / 'starts.json'


def add_start(key: str, params: Mapping[str, float], source: str) -> None:
    from experiments import E110_no_kckc as E110
    E110.merge_starts(starts_path(), {key: E110.start_entry(params, source=source)})


def part_1(names: Sequence[str] = S1_RUNS) -> pd.DataFrame:
    """Step 1 (the module docstring), resumable: a run is done when its tag is in any runs_1_*.csv."""
    from experiments import E108_paa_extrapolation as E108
    from experiments.E77_input_screen import hold_rates
    if not _read_json(out_dir() / 'check.json').get('passes'):
        raise RuntimeError('run --part check first (check.json)')
    table = out_dir() / f'runs_1_{_host()}.csv'
    prov_path = out_dir() / f'provenance_1_{_host()}.json'
    cc = E108._commits_clean()
    for name in names:
        tag = f'S1_{name}'
        if tag in set(runs_table()['tag']):
            print(f'E130 {tag}: done', flush=True)
            continue
        E108._still_clean(cc)
        c = calyx_of(name)
        p, on, off = s1_point(name)
        t0 = time.time()
        res = hold_rates(lambda: build_simulator(c, p), target_off=off, target_on=on, params=dict(p), label=f'E130 {tag}')
        o = {'on': res['on'], 'off': res['off']}
        save_run(tag, o)
        row = dict(part='step1', tag=tag, name=name, connectome=c, kckc0=is_kckc0(name),
                   multibouton_claws=J.base_kws(None, connectome=c)['multibouton_claws'],
                   start_params=json.dumps(p), thr_factor=res['thr_factor'], apl_factor=res['apl_factor'],
                   n_on=res['n_on'], n_off=res['n_off'], target_on=on, target_off=off, at_bound_thr=res['at_bound_thr'],
                   at_bound_apl=res['at_bound_apl'], n_sims=res['n_sims'], seconds=time.time() - t0,
                   **{k: float(res['params'][k]) for k in PARAM_KEYS},
                   params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
                   params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}), host=_host())
        row.update(E83.realism(o))
        row.update(type_rates(o))
        row.update(cc)
        E83._append_row(table, row)
        prov = _read_json(prov_path) or dict(part='step1', script=SCRIPT, points={})
        P = calyx_parts(c)
        prov['points'][tag] = dict(H.model_provenance(orn_deltas=P['orn']), held=dict(
            connectome=c, multibouton_claws=row['multibouton_claws'], base=P['base'], base_source=base_params_spec(c),
            start=p, targets=dict(on=on, off=off),
            change=("C2's five scales as factors on the calyx's own base, KC-KC at C2's scale" + (' 0' if is_kckc0(name) else '')
                    + "; fixed_thr re-held with APL off, then the tied APL scales with APL on (E77_input_screen.hold_rates)"),
            kckc_matrix=(dict(path=str(write_matrix(c)[0]), sha256=write_matrix(c)[2]) if c in MALE else
                         "olfsysm's hemibrain KC->KC counts (connectome_KCKC_weights)"),
            kws=res.get('kws'), breadth_sha256=res.get('breadth_sha256')))
        _write_json(prov_path, prov)
        add_start(tag, res['params'], f'step 1 held point {tag} (results/{EXPERIMENT}/runs_1_{_host()}.csv)')
        print(f"E130 {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{on} on, "
              f"{res['n_off']}/{off} off); flag '{row.get('realism_flag') or ''}'; {row['seconds']:.0f}s", flush=True)
        del res, o
        gc.collect()
    return runs_table()


# ---------------------------------------------------------------------------
# check

def e10_default_explained(rec: Mapping) -> dict:
    """Check (iv)'s reading when the first evaluation differs from C2's saved history: it passes if every column that
    differs is one of the calibrated readout's (apl_cal_*, mean_dist_cal_*, term_apl_cal, and the loss they enter) and
    differs by at most 1e-12 relative, and if today's code on C2's saved spike counts (no simulation; E72's `distances`,
    `H.apl_effect_decomposition`) gives the check's values of apl_cal_rmse and mean_dist_cal_off exactly: the
    difference is then in the calibration's last digits on the same counts, not in the run (the counts-based columns,
    rates, ρ's and block, are equal). With no difference it passes as exact."""
    from experiments import E72_wiring_population as E72
    diffs = {k: v for k, v in rec.get('history_differs', {}).items()}
    diffs.update(rec.get('row_differs', {}))
    if rec.get('exact'):
        return dict(passes=True, differences_explained='none (exact)')
    cal = all(k.startswith(('apl_cal_', 'mean_dist_cal_', 'term_apl_cal')) or k == 'loss' for k in diffs)
    tiny = all(abs(a - b) <= 1e-12 * max(abs(a), abs(b)) for a, b in diffs.values())
    o = E83.saved_fit('real')['saved']
    d = E72.distances(o, 'cal')
    dec, _ = H.apl_effect_decomposition(d['on'], d['off'], H.yang_apl_effect())
    now = dict(apl_cal_rmse=float(dec['apl_effect_rmse']), mean_dist_cal_off=float(d['off'].mean()))
    # (the history's CSV keeps 16 significant digits)
    same_now = all(k not in diffs or diffs[k][0] == float(f'{now[k]:.16g}') for k in now)
    ok = bool(rec.get('n_history_columns') and cal and tiny and same_now)
    return dict(passes=ok, differences_explained=dict(only_calibrated_columns=bool(cal), within_1e_12=bool(tiny),
                recomputed_from_saved_counts=now, check_equals_recomputation=bool(same_now),
                note=("the first evaluation's calibrated-readout columns equal today's code on C2's saved counts; C2's "
                      "history has them from the code of its run (E67), which differs in the last digits")))


def part_check(smoke: bool = True) -> dict:
    """The checks (the module docstring), each kept once passed; check.json. Raises if one fails."""
    from experiments import E110_no_kckc as E110
    path = out_dir() / 'check.json'
    res = _read_json(path)
    # (i) the bases
    b = tune_bases()
    res['bases'] = {k: v for k, v in b.items()}
    _write_json(path, res)
    # (ii) the identity on 'hemibrain-first'
    if not res.get('identity', {}).get('exact'):
        F = E83.saved_fit('real')
        P = calyx_parts(FIRST)
        P83 = E83.wiring_parts('real')
        sim, _ = build_simulator(FIRST, F['p'])
        diff = {}
        for cnd, q in (('on', dict(F['p'])), ('off', dict(F['p'], wAPLKC_scale=0.0, wAPLPN_scale=0.0))):
            r, s = sim(q)
            rs, ss = F['saved'][cnd]
            diff[cnd] = dict(counts=int((s.to_numpy(dtype=float) != ss.to_numpy(dtype=float)).sum()),
                             responses=int((r.to_numpy(dtype=float) != rs.to_numpy(dtype=float)).sum()),
                             same_index=bool(s.index.equals(ss.index) and s.columns.equals(ss.columns)))
        res['identity'] = dict(differences=diff, wiring_same_as_E83=bool(P['W'].equals(P83['W'])),
                               breadth_same_as_E83=bool(P['breadth'].equals(P83['breadth'])), checked=_now(), host=_host(),
                               exact=bool(all(d['counts'] == 0 and d['responses'] == 0 and d['same_index'] for d in diff.values())
                                          and P['W'].equals(P83['W'])))
        del sim
        gc.collect()
        _write_json(path, res)
        print(f"E130 check (ii) identity: {res['identity']}", flush=True)
    # (iii) the male wirings against olfsysm's README; every calyx's model wiring against step 0's
    counts = {}
    for c in MALE:
        W = calyx_parts(c)['W']
        got = dict(claws=int(len(W)), kcs=int(model_kc_ids(W).size), boutons=int(W.shape[1]))
        counts[c] = dict(got=got, readme=README_COUNTS[c], same=got == README_COUNTS[c])
    res['readme_counts'] = counts
    res['step0_wiring_same'] = {c: bool(calyx_parts(c)['W'].equals(step0_wiring_only(c))) for c in STEP0_CALYCES}
    _write_json(path, res)
    print(f'E130 check (iii) README counts: {counts}', flush=True)
    # (iv) E10_joint_fit.py with --connectome hemibrain-first, exact start
    if smoke and 'differences_explained' not in res.get('e10_default', {}) and not res.get('e10_default', {}).get('exact'):
        name = 'e10_default_exact'
        cmd = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-max', '3',
               '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', '_cap3_check0', '--maxfev',
               str(MAXFEV), '--experiment', f'{EXPERIMENT}/check/{name}', '--connectome', FIRST, '--init-from', C2[1],
               '--init-from-experiment', C2[0], '--init-from-exact', '--smoke']
        cc = _fit_clean()
        rec = E110.launch(cmd, H.results_dir(f'{EXPERIMENT}/check/{name}'), local_dir('logs') / f'{name}.log', cc, name)
        tag = f"rate{float(H.data_target('ctrl_frac')):.3f}_aplboth_k2_cap3_check0"
        hp = H.results_dir(f'{EXPERIMENT}/check/{name}') / f'history_{tag}.csv'
        h0 = pd.read_csv(hp).iloc[0].to_dict() if hp.exists() else {}
        hc2 = pd.read_csv(H.results_dir(C2[0]) / f'history_{C2[1]}.csv')
        best = hc2.loc[hc2['loss'].idxmin()].to_dict()
        row = pd.read_csv(H.results_dir(C2[0]) / 'fits.csv').set_index('tag').loc[C2[1]].to_dict()
        ex = E110.exact_check(h0, best, row) if h0 else dict(exact=False, note='no history')
        res['e10_default'] = dict(ex, returncode=rec['returncode'], seconds=rec['seconds'], command=cmd, checked=_now())
        _write_json(path, res)
        print(f"E130 check (iv) E10 --connectome {FIRST}: exact {res['e10_default']['exact']}", flush=True)
    if smoke and res.get('e10_default') and 'differences_explained' not in res['e10_default']:
        res['e10_default'].update(e10_default_explained(res['e10_default']))
        _write_json(path, res)
        print(f"E130 check (iv) differences explained: {res['e10_default']['differences_explained']}", flush=True)
    # (v) smoke fits on the new connectomes
    if smoke:
        for c in ('malecns-left', 'hemibrain-argmax'):
            key = f'smoke_{c}'
            if res.get(key, {}).get('passes'):
                continue
            spec = FitSpec(f'smoke_{c}', c, c, False, 'factors', 'cap3_smoke', 5)
            cmd = fit_command(spec, experiment=f'{EXPERIMENT}/check/{spec.name}') + ['--smoke']
            if c in MALE:
                write_matrix(c)
            cc = _fit_clean()
            name = spec.name
            rec = E110.launch(cmd, H.results_dir(f'{EXPERIMENT}/check/{name}'), local_dir('logs') / f'{name}.log', cc, name)
            hp = H.results_dir(f'{EXPERIMENT}/check/{name}') / f'history_{e10_tag(spec)}.csv'
            h = pd.read_csv(hp) if hp.exists() else None
            res[key] = dict(returncode=rec['returncode'], seconds=rec['seconds'], history=str(hp), n_evals=None if h is None else len(h),
                            first_loss=None if h is None else float(h['loss'].iloc[0]), command=cmd, checked=_now(),
                            passes=bool(rec['returncode'] == 0 and h is not None and len(h) >= 2))
            _write_json(path, res)
            print(f'E130 check (v) smoke {c}: {res[key]}', flush=True)
    res['passes'] = bool(res.get('identity', {}).get('exact') and all(v['same'] for v in res['readme_counts'].values())
                         and (not smoke or (res.get('e10_default', {}).get('passes')
                                            and all(res.get(f'smoke_{c}', {}).get('passes') for c in ('malecns-left', 'hemibrain-argmax')))))
    res['checked'] = _now()
    _write_json(path, res)
    if not res['passes']:
        raise RuntimeError(f'E130 check failed ({path})')
    return res


# ---------------------------------------------------------------------------
# step 2: the fits

def fit_specs(targets: Sequence[str] = FIT_TARGETS) -> List[FitSpec]:
    """Each target's two starts (the second run only when the first ends above SECOND_START_ABOVE)."""
    out = []
    for t in targets:
        c, k0 = calyx_of(t), is_kckc0(t)
        first, second = ('held', 'factors') if k0 else ('factors', 'held')
        out.append(FitSpec(t, t, c, k0, first, 'cap3', MAXFEV))
        out.append(FitSpec(f'{t}_s2', t, c, k0, second, 'cap3_s2', MAXFEV))
    return out


def e10_tag(spec: FitSpec) -> str:
    """The tag E10_joint_fit.py gives the fit (its prefix rules: cx-<calyx>_, rate<R>_, aplboth_, kckc0_, kcm-<calyx>_,
    then k2 and the --tag suffix)."""
    r = float(H.data_target('ctrl_frac'))
    return (f"{J.connectome_tag(spec.connectome)}rate{r:.3f}_aplboth_{'kckc0_' if spec.kckc0 else ''}"
            f"{'kcm-' + spec.connectome + '_' if spec.connectome in MALE else ''}k{K:g}_{spec.suffix}")


def fit_command(spec: FitSpec, python: Optional[str] = None, experiment: Optional[str] = None) -> List[str]:
    """The E10_joint_fit.py command of one fit (E72's part-B command with the calyx's options), into its own folders
    (`experiment`, default fits/<name>)."""
    py = str(H.HERE / 'venv' / 'bin' / 'python') if python is None else python
    cmd = [py, '-u', 'experiments/E10_joint_fit.py', '-k', f'{K:g}', '--aplpn-max', '3', '--rate-target', 'ctrl_frac',
           '--loss', 'apl', '--apl-readout', 'both', '--tag', f'_{spec.suffix}', '--maxfev', str(int(spec.maxfev)),
           '--experiment', experiment or f'{FIT_PREFIX}/{spec.name}', '--connectome', spec.connectome,
           '--base-params', base_params_spec(spec.connectome)]
    if spec.connectome in MALE:
        cmd += ['--kckc-matrix', str(matrix_path(spec.connectome)), '--kckc-matrix-tag', spec.connectome]
    if spec.kckc0:
        cmd += ['--kckc-fixed', '0']
    if spec.start == 'factors':
        cmd += ['--init-from', C2[1], '--init-from-experiment', C2[0]]
    elif spec.start == 'held':
        cmd += ['--init-params', f'results/{EXPERIMENT}/starts.json:S1_{spec.target}']
    else:
        raise ValueError(f'unknown start {spec.start!r}')
    return cmd


def fit_results_dir(name: str) -> Path:
    return H.shared_dir() / 'results' / FIT_PREFIX / name


def fit_row(spec: FitSpec) -> Optional[pd.Series]:
    p = fit_results_dir(spec.name) / 'fits.csv'
    if not p.exists():
        return None
    f = pd.read_csv(p)
    r = f[f['tag'] == e10_tag(spec)]
    return r.iloc[-1] if len(r) else None


def load_fit_outputs(spec: FitSpec) -> Optional[Dict[str, tuple]]:
    for d in (H.HERE / 'outputs' / FIT_PREFIX / spec.name, archive_dir() / 'fits' / spec.name):
        o = E83._load_pair(d, e10_tag(spec))
        if o is not None:
            return o
    return None


def needs_second(row: Optional[Mapping]) -> bool:
    """Whether a target's first start calls for the second (E72's rule: fit_apl_cal_rmse above 0.066)."""
    return row is not None and float(row[CHOOSE_COL]) > SECOND_START_ABOVE


def counted(target: str) -> Optional[dict]:
    """The target's counted start: of the done starts the lower CHOOSE_COL; None while the first is not done or the
    second is due and not done. {'spec', 'row', 'n_starts', 'rows'}."""
    s1, s2 = [s for s in fit_specs([target])]
    r1 = fit_row(s1)
    if r1 is None:
        return None
    cands = [(s1, r1)]
    r2 = fit_row(s2)
    if needs_second(r1):
        if r2 is None:
            return None
        cands.append((s2, r2))
    elif r2 is not None:
        cands.append((s2, r2))
    i = int(np.argmin([float(r[CHOOSE_COL]) for _, r in cands]))
    return dict(spec=cands[i][0], row=cands[i][1], n_starts=len(cands), rows=[r for _, r in cands])


def claim(name: str, lane: str, reclaim: bool = False) -> bool:
    """Claims a fit by mkdir of results/E130_malecns_family/claims/<name>/ (atomic on NFS); False if held (unless
    `reclaim`)."""
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


def _fit_clean(wait_s: float = 6 * 3600, every: float = 300) -> dict:
    """E110's `_fit_clean` (no modified file under experiments/ or kcq/, the other checkouts clean), waiting up to
    `wait_s` for another session's edits to be committed."""
    from experiments import E110_no_kckc as E110
    t0 = time.time()
    while True:
        try:
            return E110._fit_clean()
        except RuntimeError as e:
            if time.time() - t0 > wait_s:
                raise
            print(f'E130: waiting ({e})', flush=True)
            time.sleep(every)


def prerequisite(spec: FitSpec) -> Optional[str]:
    """Why a fit cannot start yet (None when it can)."""
    if spec.connectome not in _read_json(bases_path()):
        return f'no base for {spec.connectome} in bases.json (atlas: --part check)'
    if spec.start == 'held' and f'S1_{spec.target}' not in _read_json(starts_path()):
        return f'no step-1 held point S1_{spec.target} in starts.json (atlas: --part 1)'
    return None


def run_fits(targets: Sequence[str], lane: Optional[str] = None, reclaim: bool = False, dry_run: bool = False,
             wait: bool = True, wait_every: float = 300) -> List[str]:
    """Step 2: the fits of `targets` one process after another (`fit_command`; logs under outputs/.../logs/), each
    target's first start then, when due, its second; skipping those done or claimed; waiting (`wait`) for missing
    prerequisites (bases, held points, matrices)."""
    from experiments import E110_no_kckc as E110
    from experiments import E108_paa_extrapolation as E108
    lane = lane or _host()
    table = out_dir() / f'runs_fits_{lane}.csv'
    done = []
    cc = None
    for t in targets:
        for spec in fit_specs([t]):
            if spec.suffix == 'cap3_s2':
                r1 = fit_row(fit_specs([t])[0])
                if not needs_second(r1):
                    print(f'E130 {spec.name}: not due (first start {CHOOSE_COL} '
                          f"{'missing' if r1 is None else format(float(r1[CHOOSE_COL]), '.4f')})", flush=True)
                    continue
            if fit_row(spec) is not None:
                print(f'E130 {spec.name}: done ({e10_tag(spec)})', flush=True)
                done.append(spec.name)
                continue
            cmd = fit_command(spec)
            if dry_run:
                why = prerequisite(spec)
                print(f'E130 {spec.name}{" (" + why + ")" if why else ""}: ' + ' '.join(cmd), flush=True)
                continue
            while True:
                why = prerequisite(spec)
                if why is None or not wait:
                    break
                print(f'E130 {spec.name}: waiting, {why} ({_now()})', flush=True)
                time.sleep(wait_every)
            if why:
                raise RuntimeError(f'E130 {spec.name}: {why}')
            if not claim(spec.name, lane, reclaim):
                print(f'E130 {spec.name}: claimed by another process, skipped', flush=True)
                continue
            try:
                if spec.connectome in MALE:
                    write_matrix(spec.connectome)
                cc = _fit_clean() if cc is None else cc
                E108._still_clean(cc)
                print(f'E130 {spec.name}: ' + ' '.join(cmd), flush=True)
                rec = E110.launch(cmd, fit_results_dir(spec.name), local_dir('logs') / f'{spec.name}.log', cc, spec.name)
                if rec['returncode']:
                    raise RuntimeError(f"E130 {spec.name}: E10_joint_fit.py exited {rec['returncode']} (log {rec['log']})")
                r = fit_row(spec)
                if r is None:
                    raise RuntimeError(f"E130 {spec.name}: no row {e10_tag(spec)} in its fits.csv (log {rec['log']})")
            except BaseException:
                release(spec.name)
                raise
            row = dict(name=spec.name, target=spec.target, connectome=spec.connectome, kckc0=spec.kckc0, start=spec.start,
                       tag=e10_tag(spec), folder=f'{FIT_PREFIX}/{spec.name}', host=_host(), lane=lane, seconds=rec['seconds'],
                       command=' '.join(cmd), **E110.params_of_row(r), loss=float(r['loss']),
                       **{CHOOSE_COL: float(r[CHOOSE_COL])}, n_evals=int(r['n_evals']), **cc)
            E83._append_row(table, row)
            print(f"E130 {spec.name}: loss {float(r['loss']):.4f}, {CHOOSE_COL} {float(r[CHOOSE_COL]):.4f}, KC-KC "
                  f"{float(r['kc_kc_inh']):.3g}, {rec['seconds'] / 60:.0f} min", flush=True)
            done.append(spec.name)
    return done


# ---------------------------------------------------------------------------
# readings

def _dist(o: Dict[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (E83's `pair_distances`), cached by the sha256 of its spike counts and the readout."""
    import pickle
    from experiments.E77_input_screen import frame_sha256
    path = local_dir() / 'distances_cache.pkl'
    if 'dist' not in _CACHE:
        try:
            _CACHE['dist'] = pickle.loads(path.read_bytes()) if path.exists() else {}
        except Exception:
            _CACHE['dist'] = {}
    st = _CACHE['dist']
    key = (frame_sha256(o['on'][1]), frame_sha256(o['off'][1]), readout)
    if key not in st:
        st[key] = E83.pair_distances(o, readout)
        _CACHE['dist_dirty'] = True
    return st[key]


def _save_dist_cache() -> None:
    import pickle
    if _CACHE.get('dist_dirty'):
        try:
            (local_dir() / 'distances_cache.pkl').write_bytes(pickle.dumps(_CACHE['dist']))
        except OSError:
            pass
        _CACHE['dist_dirty'] = False


def group_errors(d: Dict[str, pd.Series]) -> dict:
    """The ester and alcohol groups' mean pattern error (model minus data per-pair APL effect, means removed; E77's
    `kind_masks`) against every effect target, and E92's line: the groups' APL-on distances relative to the mean pair's
    and the residual off E92's line through E89's weightings (S's populations), on this unit's readout."""
    from experiments import E77_input_side_common_error as E77
    from experiments import E108_paa_extrapolation as E108
    me = E83._cs(H.apl_effect(d['on'], d['off']))
    m = E77.kind_masks(me.index)
    row = {}
    for t, de0 in E108.effect_targets().items():
        de = E83._cs(de0).reindex(me.index)
        err = (me - me.mean()) - (de - de.mean())
        row[f'{t}_err_ester'] = float(err[m['ester']].mean())
        row[f'{t}_err_alcohol'] = float(err[m['alcohol']].mean())
    return row


def line_place(d: Dict[str, pd.Series], readout: str) -> dict:
    """E92's line: (ester, alcohol) relative APL-on distances and the residual off the line `E92.line(readout)['S']`."""
    from experiments import E92_steps_like_for_like as E92
    est, alc = E92.relative_distances(d['on'])
    L = E92.line(readout)['S']
    return dict(line_ester_rel=est, line_alcohol_rel=alc,
                line_resid=float(est - (L['intercept'] + L['slope'] * alc)))


def unit_readings(d: Dict[str, pd.Series], readout: str) -> dict:
    from experiments import E110_no_kckc as E110
    row = E110.scores_of(d)
    row.update(group_errors(d))
    try:
        row.update(line_place(d, readout))
    except Exception as e:          # (E92's line.json absent on this host: the line computed reads E89's runs)
        row['line_error'] = repr(e)[:200]
    return row


def bound_flags(row: Mapping, start: Mapping[str, float], base: Mapping[str, float]) -> dict:
    """At a bound of its fit (BOUND_TOL, relative): the threshold vs the base's, APL->KC vs the start's, KC-KC at
    KCKC_MAX; APL->bouton at its cap (3 x the base's) beside."""
    thr = float(row['fixed_thr']) / float(base['fixed_thr'])
    apl = float(row['wAPLKC_scale']) / float(start['wAPLKC_scale'])
    aplpn = float(row['wAPLPN_scale']) / float(base['wAPLPN_scale'])
    at = lambda x, lo, hi: bool(x <= lo * (1 + BOUND_TOL) or x >= hi * (1 - BOUND_TOL))
    return dict(thr_factor_vs_base=thr, aplkc_factor_vs_start=apl, aplpn_factor_vs_base=aplpn,
                thr_at_bound=at(thr, *THR_BOUNDS), aplkc_at_bound=at(apl, *APLKC_BOUNDS),
                kckc_at_bound=bool(float(row.get('wKCKC_scale', 0.0)) >= KCKC_MAX * (1 - BOUND_TOL)),
                aplpn_at_cap=bool(aplpn >= 3 * (1 - BOUND_TOL)))


def fit_start(spec: FitSpec) -> Dict[str, float]:
    """The start a fit's provenance records (E10's 'initial')."""
    p = fit_results_dir(spec.name) / f'provenance_{e10_tag(spec)}.json'
    return {k: float(v) for k, v in json.loads(p.read_text())['initial'].items()}


def collect() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """Every run there is: ({(kind, unit): outputs}, runs table). Kinds: 'family' (the eight wirings' saved fits),
    'E110_B' (E110's counted part-B fits, KC-KC 0, the eight), 'S1' (units S1 run names), 'fit' (units the targets, the
    counted start), 'fit_start' (units <target>:<start> for every start done)."""
    from experiments import E110_no_kckc as E110
    outs, rows = {}, []

    def add(kind, unit, o, **meta):
        outs[(kind, unit)] = o
        r = dict(kind=kind, unit=unit, **meta, **E83.realism(o), **type_rates(o), **E110.run_extras(o))
        r.update(E110.bands_of(float(r['rate_on']), float(r['rate_off']), float(r['rate_rho_on']), float(r['opk_on_tv']),
                               float(r['opk_off_tv'])))
        r['in_bands'] = all(r[b] for b in ('band_rate_on', 'band_ratio', 'band_rate_rho_on', 'band_opk_on', 'band_opk_off'))
        rows.append(r)

    for w in E83.ALL_WIRINGS:
        F = E83.saved_fit(w)
        add('family', w, F['saved'], tag=F['tag'], folder=F['folder'], connectome=FIRST, **F['p'])
        b = E110.counted_b(w)
        if b is not None:
            o = E110.load_fit_outputs(b['name'], b['row']['tag'])
            if o is not None:
                add('E110_B', w, o, tag=b['row']['tag'], folder=f"E110_no_kckc/fits/{b['name']}", connectome=FIRST,
                    **E110.params_of_row(b['row']))
    for r in runs_table().to_dict('records'):
        o = load_run(r['tag'])
        if o is not None:
            add('S1', r['name'], o, tag=r['tag'], folder=f'{EXPERIMENT}/step1', connectome=r['connectome'],
                thr_factor=r['thr_factor'], apl_factor=r['apl_factor'], **{k: float(r[k]) for k in PARAM_KEYS})
    for t in FIT_TARGETS:
        for spec in fit_specs([t]):
            r = fit_row(spec)
            if r is None:
                continue
            o = load_fit_outputs(spec)
            if o is None:
                print(f'E130 score: {spec.name} has a fits.csv row but no outputs here or in the archive', flush=True)
                continue
            p = E110.params_of_row(r)
            meta = dict(tag=e10_tag(spec), folder=f'{FIT_PREFIX}/{spec.name}', connectome=spec.connectome, start=spec.start,
                        loss=float(r['loss']), **{CHOOSE_COL: float(r[CHOOSE_COL])}, n_evals=int(r['n_evals']), **p)
            try:
                meta.update(bound_flags(p, fit_start(spec), calyx_base(spec.connectome)))
            except (FileNotFoundError, KeyError, RuntimeError):
                pass
            add('fit_start', f'{t}:{spec.start}', o, target=t, **meta)
        cnt = counted(t)
        if cnt is not None and ('fit_start', f"{t}:{cnt['spec'].start}") in outs:
            o = outs[('fit_start', f"{t}:{cnt['spec'].start}")]
            meta = next(dict(x) for x in rows if x['kind'] == 'fit_start' and x['unit'] == f"{t}:{cnt['spec'].start}")
            meta = {k: v for k, v in meta.items() if k not in ('kind', 'unit')}
            outs[('fit', t)] = o
            rows.append(dict(meta, kind='fit', unit=t, n_starts=cnt['n_starts'],
                             all_starts_outside_bands=all(not x['in_bands'] for x in rows
                                                          if x['kind'] == 'fit_start' and x.get('target') == t)))
    return outs, pd.DataFrame(rows)


def triples() -> List[Tuple[str, str, str]]:
    """The 21 populations of three of the eight that contain the hemibrain wiring: (real, X, Y), X < Y of the seven."""
    draws = [w for w in E83.ALL_WIRINGS if w != 'real']
    return [('real', x, y) for x, y in itertools.combinations(draws, 2)]


def population_units(outs: Mapping[Tuple[str, str], Dict[str, tuple]]) -> Dict[str, List[Tuple[str, str]]]:
    """{population name: [(kind, unit)]}: 'family8'; 'tri_<X>_<Y>' (the 21); 'real3' (C2 with the two male fits),
    'real3_argmax' (the 'argmax' refit in place of C2), 'S1_real3' (C2 with the two male step-1 runs); 'E110_B8'."""
    pops = {'family8': [('family', w) for w in E83.ALL_WIRINGS]}
    for _, x, y in triples():
        pops[f'tri_{x}_{y}'] = [('family', 'real'), ('family', x), ('family', y)]
    pops['real3'] = [('family', 'real'), ('fit', 'malecns-left'), ('fit', 'malecns-right')]
    pops['real3_argmax'] = [('fit', 'hemibrain-argmax'), ('fit', 'malecns-left'), ('fit', 'malecns-right')]
    pops['S1_real3'] = [('family', 'real'), ('S1', 'malecns-left'), ('S1', 'malecns-right')]
    pops['E110_B8'] = [('E110_B', w) for w in E83.ALL_WIRINGS]
    return {k: v for k, v in pops.items() if all(m in outs for m in v)}


RUN_COLS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv', 'max_count_on',
            'max_count_off', 'fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale') + tuple(
    f'{k}_{TYPE_COL[t]}' for t in TYPES for k in ('rate_on', 'rate_off', 'ratio'))


def readings(outs, runs) -> pd.DataFrame:
    """Per run (kind, unit) and population, per readout: `unit_readings`, the runs' mean of RUN_COLS ('run_<col>'),
    the number of members outside the bands."""
    rows = []
    idx = runs.drop_duplicates(['kind', 'unit'], keep='last').set_index(['kind', 'unit'])
    singles = {f'{k}:{u}': [(k, u)] for k, u in outs}
    for ro in READOUTS:
        dd = {m: _dist(o, ro) for m, o in outs.items()}
        for name, members in list(singles.items()) + list(population_units(outs).items()):
            d = dd[members[0]] if len(members) == 1 else E83.population([dd[m] for m in members])
            s = idx.loc[members]
            rows.append(dict(unit=name, n_members=len(members), readout=ro,
                             n_outside_bands=int((~s['in_bands'].astype(bool)).sum()),
                             **{f'run_{c}': float(s[c].astype(float).mean()) for c in RUN_COLS if c in s.columns},
                             **unit_readings(d, ro)))
    _save_dist_cache()
    return pd.DataFrame(rows)


def _v(rd: pd.DataFrame, unit: str, col: str, readout: str = 'cal') -> float:
    r = rd[(rd['unit'] == unit) & (rd['readout'] == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


# ---------------------------------------------------------------------------
# the classes (pure functions of the readings; tested on made-up values)

def class_new_default(argmax: Mapping[str, float], family: Mapping[str, float], in_bands: bool) -> dict:
    """'unchanged' if the 'argmax' refit's pattern part is within WITHIN of the family's hemibrain fit's on every target
    and inside the bands, 'changed' otherwise; no class outside the bands."""
    d = {t: float(argmax[t]) - float(family[t]) for t in TARGETS_READ}
    if not in_bands:
        return dict(cls=None, reason='outside the bands', d=d)
    return dict(cls='unchanged' if all(abs(v) <= WITHIN for v in d.values()) else 'changed', d=d)


def class_inside(v: Mapping[str, float], ranges: Mapping[str, Tuple[float, float]], in_bands: bool) -> dict:
    """A male calyx's fit against the eight wirings' range per target: 'inside' on every target, else 'outside (worse
    side)' / '(better side)' / '(both sides)'; no class outside the bands."""
    above = [t for t in TARGETS_READ if float(v[t]) > ranges[t][1]]
    below = [t for t in TARGETS_READ if float(v[t]) < ranges[t][0]]
    det = dict(above=above, below=below, values={t: float(v[t]) for t in TARGETS_READ},
               ranges={t: list(ranges[t]) for t in TARGETS_READ})
    if not in_bands:
        return dict(det, cls=None, reason='outside the bands')
    if not above and not below:
        cls = 'inside'
    elif above and not below:
        cls = 'outside (worse side)'
    elif below and not above:
        cls = 'outside (better side)'
    else:
        cls = 'outside (both sides)'
    return dict(det, cls=cls)


def class_population(real3: Mapping[str, float], tri: Mapping[str, Sequence[float]], in_bands: bool) -> dict:
    """The three real calyces against the 21 three-wiring populations: 'real calyces average better' if below the lowest
    and at least POP_BELOW below the median on every target, 'worse' the mirror, otherwise 'as the re-draws'; no class
    if a member is outside the bands."""
    det = {}
    better, worse = True, True
    for t in TARGETS_READ:
        x = np.asarray(tri[t], dtype=float)
        v = float(real3[t])
        med, lo, hi = float(np.median(x)), float(x.min()), float(x.max())
        det[t] = dict(value=v, lowest=lo, median=med, highest=hi, rank=int((x < v).sum()), n=int(len(x)))
        better &= v < lo and v <= med - POP_BELOW
        worse &= v > hi and v >= med + POP_BELOW
    if not in_bands:
        return dict(cls=None, reason='a member outside the bands', targets=det)
    return dict(cls='real calyces average better' if better else ('worse' if worse else 'as the re-draws'), targets=det)


def class_kckc(free: Mapping[str, Mapping[str, float]], zero: Mapping[str, Mapping[str, float]],
               testable: Mapping[str, bool], in_bands: bool) -> dict:
    """KC-KC on the male calyces: 'replicated' if on both the fit with KC-KC at 0 is higher by REPLICATED or more on
    every target, 'not replicated' if within WITHIN on both on every target, 'mixed' otherwise; 'not testable' if a
    calyx has fewer than half of its modelled KCs in the table; no class outside the bands."""
    d = {c: {t: float(zero[c][t]) - float(free[c][t]) for t in TARGETS_READ} for c in MALE}
    if not all(testable.get(c, False) for c in MALE):
        return dict(cls='not testable', d=d, testable=dict(testable))
    if not in_bands:
        return dict(cls=None, reason='a fit outside the bands', d=d)
    rep = all(v >= REPLICATED for c in MALE for v in d[c].values())
    notrep = all(abs(v) <= WITHIN for c in MALE for v in d[c].values())
    return dict(cls='replicated' if rep else ('not replicated' if notrep else 'mixed'), d=d)


def pattern_by_target(rd: pd.DataFrame, unit: str) -> Dict[str, float]:
    return {t: _v(rd, unit, f'{t}_pattern') for t in TARGETS_READ}


def verdict(rd: pd.DataFrame, runs: pd.DataFrame, testable: Mapping[str, bool]) -> dict:
    """The classes of the module docstring from the readings (calibrated readout); incomplete parts say so."""
    ib = runs.drop_duplicates(['kind', 'unit'], keep='last').set_index(['kind', 'unit'])['in_bands']
    have = lambda k, u: (k, u) in ib.index
    v = dict(rules=dict(within=WITHIN, replicated=REPLICATED, pop_below=POP_BELOW, targets=list(TARGETS_READ),
                        second_start_above=SECOND_START_ABOVE, choose_col=CHOOSE_COL))
    fam = pattern_by_target(rd, 'family:real')
    v['new_default'] = (class_new_default(pattern_by_target(rd, 'fit:hemibrain-argmax'), fam, bool(ib[('fit', 'hemibrain-argmax')]))
                        if have('fit', 'hemibrain-argmax') else dict(cls='incomplete'))
    ranges = {t: (min(_v(rd, f'family:{w}', f'{t}_pattern') for w in E83.ALL_WIRINGS),
                  max(_v(rd, f'family:{w}', f'{t}_pattern') for w in E83.ALL_WIRINGS)) for t in TARGETS_READ}
    v['male'] = {c: (class_inside(pattern_by_target(rd, f'fit:{c}'), ranges, bool(ib[('fit', c)])) if have('fit', c)
                     else dict(cls='incomplete')) for c in MALE}
    tris = [f'tri_{x}_{y}' for _, x, y in triples()]
    if (rd['unit'] == 'real3').any():
        tri = {t: [_v(rd, u, f'{t}_pattern') for u in tris] for t in TARGETS_READ}
        v['population'] = class_population(pattern_by_target(rd, 'real3'), tri,
                                           bool(ib[('fit', 'malecns-left')] and ib[('fit', 'malecns-right')] and ib[('family', 'real')]))
        if (rd['unit'] == 'real3_argmax').any():
            v['population_argmax_beside'] = class_population(pattern_by_target(rd, 'real3_argmax'), tri, True)
    else:
        v['population'] = dict(cls='incomplete')
    if all(have('fit', c) and have('fit', f'{c}-kckc0') for c in MALE):
        v['kckc'] = class_kckc({c: pattern_by_target(rd, f'fit:{c}') for c in MALE},
                               {c: pattern_by_target(rd, f'fit:{c}-kckc0') for c in MALE}, testable,
                               all(bool(ib[('fit', u)]) for c in MALE for u in (c, f'{c}-kckc0')))
    else:
        v['kckc'] = dict(cls='incomplete')
    # the hemibrain's eight beside (E110's part B against the family): population and per wiring
    if (rd['unit'] == 'E110_B8').any():
        v['kckc_hemibrain_eight'] = dict(
            population={t: _v(rd, 'E110_B8', f'{t}_pattern') - _v(rd, 'family8', f'{t}_pattern') for t in TARGETS_READ},
            per_wiring={w: _v(rd, f'E110_B:{w}', 'all_pattern') - _v(rd, f'family:{w}', 'all_pattern') for w in E83.ALL_WIRINGS})
    v['realism_flags'] = [dict(kind=r['kind'], unit=r['unit'], flag=r.get('realism_flag'),
                               bands={b: bool(r[b]) for b in ('band_rate_on', 'band_ratio', 'band_rate_rho_on', 'band_opk_on', 'band_opk_off')})
                          for r in runs.to_dict('records') if r['kind'] in ('S1', 'fit', 'fit_start') and not r['in_bands']]
    v['at_bounds'] = [dict(unit=r['unit'], **{k: r.get(k) for k in ('thr_at_bound', 'aplkc_at_bound', 'kckc_at_bound', 'aplpn_at_cap')})
                      for r in runs.to_dict('records') if r['kind'] == 'fit_start'
                      and any(bool(r.get(k)) for k in ('thr_at_bound', 'aplkc_at_bound', 'kckc_at_bound'))]
    v['no_class_both_starts'] = sorted({r['unit'] for r in runs.to_dict('records') if r['kind'] == 'fit'
                                        and r.get('all_starts_outside_bands')})
    return v


# ---------------------------------------------------------------------------
# score and variants

STATS = (('pattern part, readout', 'all_pattern', 'cal', '.4f'), ('pattern part, half A / B / C / D', None, 'cal', '.3f'),
         ('amplitude part, readout', 'all_amplitude', 'cal', '.4f'), ('shape part, readout', 'all_shape', 'cal', '.4f'),
         ('spread ratio, readout', 'all_sd_ratio', 'cal', '.3f'),
         ('Pearson r on the per-pair APL effect, readout', 'all_pearson_effect', 'cal', '.3f'),
         ('ρ on the per-pair APL effect, readout', 'all_rho_effect', 'cal', '.3f'),
         ('mean APL effect, readout', 'mean_effect', 'cal', '+.3f'), ('APL-effect loss, raw', 'all_loss', 'raw', '.4f'),
         ('mean APL effect, raw', 'mean_effect', 'raw', '+.3f'), ('pattern part, raw', 'all_pattern', 'raw', '.4f'),
         ('ρ on odor-odor distances, APL on', 'rho_on', 'raw', '.3f'), ('ρ on odor-odor distances, APL off', 'rho_off', 'raw', '.3f'),
         ('ester group pattern error, readout', 'all_err_ester', 'cal', '+.3f'),
         ('alcohol group pattern error, readout', 'all_err_alcohol', 'cal', '+.3f'),
         ("residual off E92's line, raw", 'line_resid', 'raw', '+.3f'),
         ('response rate, APL on', 'run_rate_on', 'raw', '.4f'), ('response rate, APL off', 'run_rate_off', 'raw', '.4f'),
         ('APL-off / APL-on response rate', 'run_ratio', 'raw', '.2f'), ('response rate ρ, APL on', 'run_rate_rho_on', 'raw', '.3f'),
         ('odors per KC (TV), APL on', 'run_opk_on_tv', 'raw', '.3f'), ('odors per KC (TV), APL off', 'run_opk_off_tv', 'raw', '.3f'),
         ('γ response rate, APL on / off', None, 'raw', '.3f'), ("α′β′ response rate, APL on / off", None, 'raw', '.3f'),
         ('αβ response rate, APL on / off', None, 'raw', '.3f'), ('γ off / on ratio', 'run_ratio_g', 'raw', '.2f'),
         ('largest spike count, APL on / off', None, 'raw', '.0f'), ('somatic threshold', 'run_fixed_thr', 'raw', '.1f'),
         ('APL->KC scale', 'run_wAPLKC_scale', 'raw', '.4f'), ('APL->bouton scale', 'run_wAPLPN_scale', 'raw', '.4f'),
         ('KC-KC scale', 'run_wKCKC_scale', 'raw', '.2f'), ('runs outside the bands', 'n_outside_bands', 'raw', '.0f'))


def comparison_table(rd: pd.DataFrame, units: Sequence[str]) -> str:
    us = [u for u in units if (rd['unit'] == u).any()]
    L = ['| statistic | ' + ' | '.join(us) + ' |', '|---|' + '---:|' * len(us)]
    pairs = {'γ response rate, APL on / off': ('run_rate_on_g', 'run_rate_off_g'),
             "α′β′ response rate, APL on / off": ('run_rate_on_apbp', 'run_rate_off_apbp'),
             'αβ response rate, APL on / off': ('run_rate_on_ab', 'run_rate_off_ab'),
             'largest spike count, APL on / off': ('run_max_count_on', 'run_max_count_off')}
    for lab, col, ro, sp in STATS:
        if lab in pairs:
            a, b = pairs[lab]
            cells = [f'{_f(_v(rd, u, a, ro), sp)} / {_f(_v(rd, u, b, ro), sp)}' for u in us]
        elif col is None:
            cells = [' / '.join(_f(_v(rd, u, f'{h}_pattern', ro), sp) for h in TARGETS_READ[1:]) for u in us]
        else:
            cells = [_f(_v(rd, u, col, ro), sp) for u in us]
        L.append(f'| {lab} | ' + ' | '.join(cells) + ' |')
    return '\n'.join(L)


def runs_bands_table(runs: pd.DataFrame) -> str:
    cols = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'opk_on_tv', 'opk_off_tv', 'rate_on_g', 'ratio_g', 'max_count_on',
            'max_count_off', 'fixed_thr', 'wKCKC_scale')
    L = ['| kind | unit | ' + ' | '.join(cols) + ' | in bands | flag |', '|---|---|' + '---:|' * len(cols) + '---|---|']
    for r in runs.itertuples():
        L.append(f'| {r.kind} | {r.unit} | ' + ' | '.join(_f(getattr(r, c, np.nan)) for c in cols)
                 + f" | {'yes' if r.in_bands else 'no'} | {r.realism_flag if isinstance(r.realism_flag, str) else ''} |")
    return '\n'.join(L)


def score() -> dict:
    out = out_dir()
    outs, runs = collect()
    rd = readings(outs, runs)
    rd.to_csv(out / 'readings.csv', index=False)
    runs.to_csv(out / 'runs.csv', index=False)
    testable = _read_json(out / 'step0.json').get('testable', {})
    v = verdict(rd, runs, testable)
    v.update(generated=_now(), host=_host(), n_runs=int(len(runs)))
    _write_json(out / 'verdict.json', v)
    tris = [f'tri_{x}_{y}' for _, x, y in triples()]
    pops = rd[rd['unit'].isin(['family8', 'real3', 'real3_argmax', 'S1_real3', 'E110_B8'] + tris)]
    pops.to_csv(out / 'populations.csv', index=False)
    units = (['family:real'] + [f'S1:{n}' for n in S1_RUNS] + [f'fit:{t}' for t in FIT_TARGETS]
             + ['real3', 'real3_argmax', 'S1_real3', 'family8'])
    fam = ['family:' + w for w in E83.ALL_WIRINGS]
    md = ['# E130 comparison tables (generated)', '',
          'Classes: ' + '; '.join(f"{k} {json.dumps(v[k].get('cls') if isinstance(v[k], dict) and 'cls' in v[k] else {c: x.get('cls') for c, x in v[k].items()})}"
                                   for k in ('new_default', 'male', 'population', 'kckc')), '',
          '## The calyces, their populations and the family', '', comparison_table(rd, units), '',
          '## The eight wirings', '', comparison_table(rd, fam), '',
          '## Runs and bands', '', runs_bands_table(runs), '']
    (out / 'comparison_table.md').write_text('\n'.join(md))
    _write_json(out / 'provenance_score.json', dict(generated=_now(), host=_host(), script=SCRIPT, commits=H.code_commits()))
    print(json.dumps({k: v[k] for k in ('new_default', 'male', 'population', 'kckc')}, default=str)[:4000], flush=True)
    return v


def write_variants() -> pd.DataFrame:
    """variants.csv (every run of this experiment, absolute parameters, the connectome and its claws' boutons) and
    model_pair_distances_long.csv (raw spike counts) for the model-output archive and kcq/apl_effect.py."""
    out = out_dir()
    outs, runs = collect()
    rows, parts = [], []
    for r in runs[runs['kind'].isin(['S1', 'fit_start'])].to_dict('records'):
        o = outs[(r['kind'], r['unit'])]
        c = r['connectome']
        kws = J.base_kws(None, connectome=c)
        name = r['unit'] if r['kind'] == 'S1' else r['target']
        what = ("C2's factors on the calyx's base, both response rates held (step 1)" if r['kind'] == 'S1'
                else f"refitted as E72's part B from the {r['start']} start (step 2)")
        rows.append(dict(tag=r['tag'], variant=r['tag'], part='S1' if r['kind'] == 'S1' else 'fit', base=f'bouton_{c}', k=K,
                         connectome=c, pn2kc_connections=kws.get('pn2kc_connections', 'hemibrain'),
                         multibouton_claws=kws['multibouton_claws'], target=name, kckc0=is_kckc0(name),
                         kind=f"C2's family on {SHORT[c]}{', KC-KC 0' if is_kckc0(name) else ''}: {what}",
                         kckc_matrix=(f'matrices/kckc_{c}.parquet' if c in MALE else 'hemibrain (olfsysm)'),
                         base_params=json.dumps(calyx_base(c)),
                         **{k: r.get(k) for k in PARAM_KEYS}, kc_kc_inh=r.get('wKCKC_scale'),
                         rate_on=r['rate_on'], rate_off=r['rate_off'], max_count_on=r['max_count_on'],
                         max_count_off=r['max_count_off'], realism_flag=r['realism_flag'],
                         subdir='step1' if r['kind'] == 'S1' else str(r['folder']).split('/', 1)[1]))
        parts.append(pd.concat({(r['tag'], 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r['tag'], 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    # (for archive_tools/migrate_model_output_params.py: each variant's base kwargs by its `base` name)
    _write_json(out / 'provenance.json', dict(
        generated=_now(), host=_host(), script=SCRIPT, **H.code_commits(), kc_kc_tau=J.TAU, breadth_k=K,
        bases={f'bouton_{c}': J.base_kws(None, connectome=c) for c in CALYCES},
        base_params={c: calyx_base(c) for c in CALYCES},
        note=("C2's family (E67) on each calyx: E10_joint_fit.py's bouton base with --connectome <calyx> (its kwargs under "
              "bases), each calyx's own tuned base (base_params, bases.json), breadth-weighted APL->KC (exponent 2), the "
              "calyx's KC-KC matrix (matrices/; the hemibrain's from olfsysm); per-run provenance in provenance_1_<host>.json "
              "and the fits' provenance_<tag>.json")))
    return df


def print_commands() -> None:
    """The lanes' command lines (this script) and every fit's E10_joint_fit.py command."""
    pre = 'OMP_NUM_THREADS={t} OMP_THREAD_LIMIT={t} PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u ' + SCRIPT
    for lane, ts in LANES.items():
        print(f"{lane}: {pre.format(t=16 if lane.startswith('hal') else 8)} --part fit --targets {' '.join(ts)} --lane {lane}")
    for spec in fit_specs():
        print(f'{spec.name}: ' + ' '.join(fit_command(spec, python='./venv/bin/python')))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('step0', 'check', '1', 'fit', 'score', 'variants', 'commands'))
    ap.add_argument('--targets', nargs='+', default=list(FIT_TARGETS), choices=list(FIT_TARGETS))
    ap.add_argument('--lane', default=None, help='the lane (default the host): its runs_fits_<lane>.csv and claims')
    ap.add_argument('--reclaim', action='store_true', help="take over the targets' claims (a process that died)")
    ap.add_argument('--dry-run', action='store_true', help='with --part fit: print the commands only')
    ap.add_argument('--no-smoke', action='store_true', help='with --part check: skip the E10_joint_fit.py processes')
    a = ap.parse_args()
    if a.part == 'step0':
        step0()
    elif a.part == 'check':
        part_check(smoke=not a.no_smoke)
    elif a.part == '1':
        part_1()
    elif a.part == 'fit':
        run_fits(a.targets, a.lane, a.reclaim, a.dry_run)
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'commands':
        print_commands()


if __name__ == '__main__':
    main()
