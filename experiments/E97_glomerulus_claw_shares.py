#!/usr/bin/env python3
"""E97: the glomeruli's shares of the claws in other connectomes. Every wiring that the best family has been read on
gives each glomerulus the hemibrain's number of claws; whether the proportions of FlyWire's and the male CNS's mushroom
bodies, put onto the model's wirings by moving claws, move the per-pair APL effect; and whether the error common to all
wirings is in part one fly's proportions.

The entry is item E97 of experiment_priorities.md (proposed by hal 2026-11-15 11:38; planned by a Fable subagent,
claude-fable-5-2, spawned from hal 2026-11-15 11:13). Built on E83 (the eight wirings at their own fits, the readings,
the populations), E85 and E93 (the held reading: the somatic threshold and the tied APL scale bisected so that both
response rates are the fit's), E89 (the targets, the co-response readings), E90 / E92 (the groups, the line, the six
near pairs), E86 (part E's searched factors) and E10_joint_fit.py (`--wiring-file --wiring-tag`, part D).

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B.
S: real, P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the recordings E77's, C / D the second split's (E89's
`targets89`). The input is the lab's except in B6 (E93's fuller input F).

Definitions (the entry's, fixed before any run):
- The derived tables (read only; Yang's own Claude's, single-linkage clusters of 2 um of each KC's PN synapses in the
  calyx): `<dataset>_KC_claws_derived.csv.gz` with the official ROIs and, in with_outside_ROI_synapses/, with the
  synapses outside the official ROIs added. A claw counts when it has SYN_MIN (5) or more PN synapses and its main PN's
  type names one glomerulus (no '+', not M); glomerulus names are mapped to the model's (RENAMES: the loader merges the
  hemibrain's VC3l and VC3m into VC3; FlyWire's VC3 and VC5 are those two, its VM6 the hemibrain's VC5). A KC's
  mushroom body is the side of the APL it makes most of its synapses with (`kc_sides`).
- A mushroom body's share for glomerulus g: its counted claws of g over its counted claws on the model's 54 glomeruli.
  Its ratio r_g: that share over the derived hemibrain's share; r_g = 1 (and listed) for a glomerulus with fewer than
  SMALL_CLAWS (30) counted claws in either.
- The five mushroom bodies: the hemibrain's, FlyWire's two and the male CNS's two (BANC's two for the record only). The
  mean proportions: the geometric mean of r_g over the three brains, the hemibrain's being 1 and the two sides of a brain
  averaged first (on log2).
- A wiring at proportions r (`move_claws`): within each KC type, the target claw count of g is the wiring's own count
  times r_g, rescaled so that the type's total is unchanged (largest-remainder rounding). Glomeruli above target lose
  claws chosen at random; each lost claw goes to a glomerulus below target, onto one of its boutons drawn in proportion
  to the boutons' claw counts in that KC type. Every KC keeps its claws (rows), its type and its rows' index; the matrix
  keeps the base's index and columns, as `E10_joint_fit.py --wiring-file` requires. A moved claw keeps its own APL->KC
  and KC->APL weights (indexed by claw) and takes the APL->bouton weight of its new bouton.
- A held reading (`held97`): the built wiring's breadth recomputed (`E10_joint_fit.bouton_base(wiring_matrix=...)`,
  its APL-off breadth run cached under outputs/E97_glomerulus_claw_shares/breadth/), the wiring's own fit's parameters,
  the somatic threshold bisected to the fit's APL-off response count, then the tied APL scale to its APL-on count
  (E77's `hold_rates`, as E85's and E93's held reading).

Parts (CLI `--part`):
- 0, the identities, before any scored point (identity.json; parts B and D refuse to run until it passes):
  - unchanged_<wiring>: with every r_g = 1 the builder returns each of the eight wirings unchanged; each wiring's matrix
    and claw table are cached under outputs/.../W/ (part A reads the real wiring's);
  - held_real: the held reading of the real wiring through this experiment's path (`build97` with the unchanged
    matrix, its breadth recomputed) reproduces its saved run with 0 count differences;
  - built (needs part A's ratios.csv, so run again after part A): every wiring that part B or C2 builds has each KC's
    claw count and type the base's, the per-type totals unchanged, the glomerulus totals per type equal to the targets,
    every moved claw on a bouton of its new glomerulus and every other row unchanged (`check_built`; no simulation);
    each built matrix's sha256 is recorded, and part B refuses a matrix that differs.
- A (no simulations; part 0's cache): A1 the counts (a1_counts.csv: per dataset, ROI version, mushroom body, threshold,
  KC class and glomerulus the claws, the boutons and PNs that carry them; a1_model.csv the model's wiring;
  a1_hemibrain_vs_model.csv and A1's rule; a1_roi.csv the two ROI versions per mushroom body; a1_flywire_pratyush.csv,
  FlyWire's proportions a second way, from olfsysm's FlyWire PN->KC table; a1_sides.csv the KCs' sides); A2 the ratios
  (ratios.csv: per glomerulus the counts, shares and ratios of the seven mushroom bodies, the mean proportions, B3's
  strata; a2_spread.csv against the counting noise; a2_agreement.csv, two sides of a brain against two brains;
  a2_departures.csv, where the hemibrain is furthest from the other brains); A3 the searched inputs against the mean
  proportions (a3_searched.csv per glomerulus, a3_rule.csv); A4 the claw-weighted input distances by group
  (a4_distances.csv); a_summary.json. Raises after writing them if A1's rule stops the experiment.
- B (`--arm` B1..B6, default all; `--wirings`, default all8; resumable, runs_B.csv; runs under
  outputs/E97_glomerulus_claw_shares/runs/ with the wiring matrix, the moved claws and the targets per type and
  glomerulus): held readings, 116 points: B1 each of the four other mushroom bodies' proportions on the eight (32); B2
  the mean proportions on the eight (8); B3 the controls on the eight (the mean proportions permuted within strata, two
  seeds, 16; inverted, 8; each mushroom body's permuted, one seed, 32); B4 a second draw of the moved claws at the mean
  proportions on S (4); B5 the mean proportions at half and twice their log ratios on S (8); B6 the mean proportions with
  E93's fuller input on the eight (8).
- score (`--no-co`, `--no-zcell`): the readings, part C's verdict and the tables.
- C2: only on a verdict of "a part" or "moves", the second draw on T (4 points, with B4), then `score`.
- D (`--pop S|T`, `--wirings` for a lane, `--dry-run`, `--force`): only on "a part" or "moves", C2's family fitted again
  at the mean proportions (E10_joint_fit.py with B2's matrix as --wiring-file).
- Dscore, E (documented stubs), variants, table.

Deviations from the entry, forced by the code or the data, or readings the entry leaves open:
- Names: the model has 54 glomeruli because olfsysm's loader merges the hemibrain's VC3l and VC3m into VC3; the derived
  tables are read the same way (the hemibrain's VC3l and VC3m, FlyWire's VC3 and VC5, are VC3; FlyWire's VM6 is the
  model's VC5). The shares are over the claws on the model's 54 glomeruli (the planning probe's over all uniglomerular
  labels; the builder's rescaling within each type makes the targets independent of that choice).
- Sides: a KC's mushroom body is the majority side of the APL it makes synapses with; the 30 to 52 BANC KCs per side
  without an APL synapse take the side of the largest gap in the KCs' mean claw x (both rules and the KC's own side
  label agree on every KC that has all three; a1_sides.csv).
- ROI versions: the rule "a mushroom body whose two ROI versions differ by more than that enters with the outside-ROI
  version" is read per mushroom body, as the SD over the glomeruli with A1_RULE_MIN_CLAWS model claws of log2 of the
  share with the outside synapses over the official one (a1_roi.csv); the male CNS's left side enters outside, its right
  side official.
- The mean proportions average a brain's two sides on log2 (a geometric mean, as the outer one). The small rule is
  applied per mushroom body; the mean proportions are 1 by rule only where all four other bodies are small (DA3, DA4m,
  VP4), and B3's permutations keep those (and, for one body's permutation, that body's small glomeruli) in place.
- B3's strata: ORN input (the lab's 36) or not, then thirds of the real wiring's claw total within that class.
- The moved claw's bouton is drawn in proportion to the boutons' claw counts within the claw's KC type in the base
  wiring (all types where a glomerulus has none in that type, which does not happen at these targets).
- A3: glomeruli whose mean-proportion ratio is 1 by rule are left out (none of the 36); for the nine glomeruli the
  permutations are of the labels over the 36, read on the nine.
- A4 is over the 66 megamat pairs (E93's `weighted_distances`; the groups are the megamat groups).
- Part C: "a part" requires the controls' conditions (both permuted seeds below half the gain, the inverted below
  0.0015); the other conditions of "moves" are in verdict.json beside it. (f)'s "another raises it" is read as any
  other body's pattern part above the reference's on the eight. R2 "not above the reference" is R2's pattern part at or
  below the eight wirings' as they are. The bootstrap resamples each genotype's recordings with replacement (BOOT_N).
- The brain term is given beside the move's own draw term (B2 against B4 on S) and the draw term between the eight
  wirings as they are, computed the same way here, as well as beside E72's and E85's 0.006.
- The co-response readings (E89's) are of the reference, B2 and the first permuted and the inverted controls only.
- Part D starts each fit from B2's held point (`--init-params`, as E93's part E): E10_joint_fit.py's exact start checks
  the wiring arguments of the base, which a moved wiring does not have.
- The builder is in this module, not in pyolfsysm.wpnkc_shuffles (the entry makes that optional).

Outputs: results/E97_glomerulus_claw_shares/ (shared); spike counts, responses, wiring matrices, moved-claw lists and
targets under outputs/E97_glomerulus_claw_shares/ (runs/, W/, breadth/).

Usage (from the repository root; 16-thread lanes; check `uptime` first; long runs with -u):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E97_glomerulus_claw_shares --part 0
    ./venv/bin/python -u -m experiments.E97_glomerulus_claw_shares --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E97_glomerulus_claw_shares --part 0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E97_glomerulus_claw_shares \\
        --part B --wirings S > <log> 2>&1
    ./venv/bin/python -u -m experiments.E97_glomerulus_claw_shares --part score
"""
from __future__ import annotations

import argparse
import json
import socket
import time
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E97_glomerulus_claw_shares'
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS

# ---------------------------------------------------------------------------
# the derived tables (read only)

DERIVED_DIR = Path('/local/tensor/Tom/claude_APL_modeling_YZ4TC/Claude_APL_analysis/apl_synapses')
ROI_VERSIONS = {'official': DERIVED_DIR, 'outside': DERIVED_DIR / 'with_outside_ROI_synapses'}
DATASETS = ('hemibrain', 'flywire', 'malecns', 'banc')
# the model's names (olfsysm's loader merges the hemibrain's VC3l and VC3m into VC3, `glomerulus_renames`); FlyWire's
# names (FlyWire, the male CNS, BANC) of the hemibrain's glomeruli as the derived tables' README and olfsysm's FlyWire
# loader give them: VC3 = the hemibrain's VC3l, VC5 = its VC3m, VM6 = its VC5
RENAMES = {'hemibrain': {'VC3l': 'VC3', 'VC3m': 'VC3'},
           'flywire': {'VC5': 'VC3', 'VM6': 'VC5'}, 'malecns': {'VC5': 'VC3', 'VM6': 'VC5'},
           'banc': {'VC5': 'VC3', 'VM6': 'VC5'}}
# the mushroom bodies: the hemibrain's (right side only), and both sides of the other three
MB_SIDES = {'hemibrain': ('R',), 'flywire': ('L', 'R'), 'malecns': ('L', 'R'), 'banc': ('L', 'R')}
FIVE_MBS = ('hemibrain', 'flywire_L', 'flywire_R', 'malecns_L', 'malecns_R')
OTHER_MBS = FIVE_MBS[1:]
BANC_MBS = ('banc_L', 'banc_R')
ALL_MBS = FIVE_MBS + BANC_MBS
BRAINS = {'flywire': ('flywire_L', 'flywire_R'), 'malecns': ('malecns_L', 'malecns_R')}
SYN_THRESHOLDS = (3, 5, 10)
SYN_MIN = 5                                         # a counted claw has this many PN synapses or more
SMALL_CLAWS = 30                                    # fewer counted claws in either mushroom body: r_g = 1
KC_CLASSES = ('g', 'ab', "a'b'")                    # the model's KC types; derived types map by prefix (`kc_class`)
FLYWIRE_PRATYUSH = 'from_pratyush/2024-09-13'       # olfsysm's FlyWire PN->KC table (A1's second way)
PRATYUSH_MIN_WEIGHT = 4                             # the loader's threshold on a FlyWire PN->KC connection

# ---------------------------------------------------------------------------
# part A's rules (fixed in the entry)

A1_RULE_SD = 0.15                                   # SD of log2(derived hemibrain share / model share) at most this
A1_RULE_MIN_CLAWS = 100                             # over the glomeruli with this many claws or more (the model's)
A3_NULL = 10_000                                    # permutations of the glomerulus labels
A3_QUANTILES = (0.05, 0.95)                         # "against" below the null's 5th percentile, "with" above its 95th
A3_SEED = 0
# the nine glomeruli of the eleven recurring single-glomerulus changes of the reduced searched input (E86's part E
# entries that meet its rule at the read size, `meets_read`; checked against partE_entries.csv)
A3_NINE = ('DC2', 'DM1', 'DM2', 'DM3', 'DM4', 'DM5', 'DP1m', 'VL2p', 'VM7d')
NOISE_DRAWS = 1000                                  # A2: multinomial draws for the counting noise
NOISE_SEED = 0

# ---------------------------------------------------------------------------
# part B (fixed before any run)

DRAW_SEEDS = (1, 2)                                 # the draws of the moved claws (B4 and C2 the second)
PERM_SEEDS = (1, 2)                                 # B3: the mean proportions permuted; each body's with the first
N_TERCILES = 3                                      # B3's permutation strata: ORN input or not x thirds of the claw total
B5_SCALES = (0.5, 2.0)                              # B5: the mean proportions' log ratios times these
B6_INPUT = 'F'                                      # E93's fuller input
ARMS = ('B1', 'B2', 'B3', 'B4', 'B5', 'B6')
ARMS_ON = {'B1': ALL_WIRINGS, 'B2': ALL_WIRINGS, 'B3': ALL_WIRINGS, 'B4': S_WIRINGS, 'B5': S_WIRINGS,
           'B6': ALL_WIRINGS, 'C2': T_WIRINGS}
ARM_ORDER = ('B2', 'B1', 'B3', 'B4', 'B5', 'B6')    # a lane's order (the main reading first)

# ---------------------------------------------------------------------------
# the references (C2's family at its fits; pattern part under the calibrated readout against all recordings, E83)

REF_S, REF_T, REF_8 = 0.0500, 0.0486, 0.0482
DRAW_TERM, INPUT_FLY_TERM = 0.006, 0.006            # E72's draw term (one wiring), E85's input-fly term
FLIES_MEAN = 0.178                                  # the flies' mean APL effect
BOOT_N = 2000                                       # paired bootstrap over the recordings
BOOT_SEED = 0

# ---------------------------------------------------------------------------
# part C, the verdict (fixed in the entry; thresholds on the pattern part)

NO_MOVE = 0.0015                                    # R1 and R2 within this of the reference
MOVE_8, MOVE_ST = 0.003, 0.0015                     # "moves": R1 this much lower on the eight, and on S and T each
PART_LO = 0.0015                                    # "a part": R1 0.0015 to 0.003 lower, the controls as for "moves"
CONTROL_RATIO = 0.5                                 # the permuted mean proportions gain less than half as much
INVERTED_MAX = 0.0015                               # the inverted ones gain less than this
SIX_KEEP = 0.5                                      # against the six typical recordings at least half of the gain kept
RAW_SLACK, MEAN_SLACK = 0.005, 0.01                 # the raw-count loss not above the reference's by more; the mean
WORSE = 0.0015                                      # R1 this much higher: "worse"
AVERAGING = 0.0015                                  # R2 gains this and its permuted control at least half as much
DISAGREE = 0.003                                    # (f): one mushroom body lowers the eight's pattern part this much
B6_DIFF = 0.0015                                    # (h): B6 differs from B2 on the eight by this much

# ---------------------------------------------------------------------------
# part D (only on "a part" or "moves")

D_MAXFEV = 200
D_HOLD = 0.5                                        # S's population keeps at least half of R1's gain on S
D_PASS_GAIN, D_WIRINGS_LOWER = 0.003, 6             # the eight by E91's rule for a fitted family
D_RAW_REF = 0.073

_CACHE: Dict[object, object] = {}


def out_dir() -> Path:
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _commits_clean() -> dict:
    """`H.code_commits()` at the start of a part, which must be clean (E89's `_commits_clean`)."""
    from experiments import E89_apl_breadth_weighting as E89
    return E89._commits_clean()


def _wirings(arg: str) -> Tuple[str, ...]:
    return {'S': S_WIRINGS, 'T': T_WIRINGS, 'all8': ALL_WIRINGS}.get(arg) or tuple(x for x in arg.split(',') if x)


def _rms(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x ** 2)))


def _sha(df: pd.DataFrame) -> str:
    from experiments.E77_input_screen import frame_sha256
    return frame_sha256(df)


# ---------------------------------------------------------------------------
# the derived tables: claws, boutons, sides, names

def glomerulus_name(pn_type, dataset: str) -> Optional[str]:
    """The model's (hemibrain) name of a derived table's PN type (e.g. 'VC3_adPN' in FlyWire -> 'VC3'; 'VM6_adPN' in
    FlyWire -> 'VC5'; 'VC3l_adPN' in the hemibrain -> 'VC3'), or None when the type names more than one glomerulus
    ('+'), is multiglomerular (M) or untyped. Names that are not among the model's glomeruli are returned as they are
    (and left out of the shares by the glomeruli they are read over)."""
    if pn_type is None or (isinstance(pn_type, float) and np.isnan(pn_type)):
        return None
    t = str(pn_type).strip()
    if not t or t.lower() == 'nan':
        return None
    g = t.split('_')[0]
    if not g or '+' in g or g == 'M':
        return None
    return RENAMES.get(dataset, {}).get(g, g)


def kc_class(kc_type) -> str:
    """A KC type's class: 'g' (KCg*, or the model's 'g'), 'ab' (KCab*), "a'b'" (KCa'b'*, KCapbp*), 'unknown' (KC
    alone, untyped)."""
    t = str(kc_type)
    if t in KC_CLASSES:
        return t
    if t.startswith('KCg'):
        return 'g'
    if t.startswith("KCa'b'") or t.startswith('KCapbp'):
        return "a'b'"
    if t.startswith('KCab'):
        return 'ab'
    return 'unknown'


def _side(x) -> Optional[str]:
    s = str(x).strip().lower()
    if s in ('l', 'left', 'apl_l'):
        return 'L'
    if s in ('r', 'right', 'apl_r'):
        return 'R'
    return None


def _ids(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors='coerce').astype('Int64')


def _read_claws(dataset: str, roi: str) -> pd.DataFrame:
    key = ('claws_raw', dataset, roi)
    if key not in _CACHE:
        df = pd.read_csv(ROI_VERSIONS[roi] / f'{dataset}_KC_claws_derived.csv.gz')
        df['kc_id'] = _ids(df['kc_id'])
        df['main_PN_id'] = _ids(df['main_PN_id'])
        _CACHE[key] = df
    return _CACHE[key]


def kc_sides(dataset: str, roi: str = 'official') -> pd.DataFrame:
    """Per KC of a dataset's claw table: 'side' ('L' / 'R'), the majority side of the APLs it makes synapses with
    (`<dataset>_APL_MB_synapses.csv.gz` of the same ROI version), and beside it 'apl_side' (that majority, None without
    an APL synapse or on a tie), 'own_side' (the KC's own side label, partner_side, majority), 'gap_side' (the side by the
    largest gap of the KCs' mean claw x, labelled by its majority agreement with the APL sides), 'n_apl_syn'. A KC
    without an APL side takes its gap side. The hemibrain's KCs are all 'R' (its one APL)."""
    key = ('sides', dataset, roi)
    if key in _CACHE:
        return _CACHE[key]
    claws = _read_claws(dataset, roi)
    kx = claws.groupby('kc_id')['centroid_x_nm'].mean().sort_values()
    out = pd.DataFrame(index=kx.index)
    if dataset == 'hemibrain':
        out['apl_side'] = 'R'
        out['own_side'] = 'R'
        out['gap_side'] = 'R'
        out['n_apl_syn'] = np.nan
        out['side'] = 'R'
        _CACHE[key] = out
        return out
    cols = {'apl_side', 'apl_instance', 'partner_id', 'partner_side'}
    syn = pd.read_csv(ROI_VERSIONS[roi] / f'{dataset}_APL_MB_synapses.csv.gz', usecols=lambda c: c in cols, dtype=str)
    syn['partner_id'] = _ids(syn['partner_id'])
    syn = syn[syn.partner_id.isin(kx.index)]
    aside = syn['apl_side'].map(_side) if 'apl_side' in syn else pd.Series(None, index=syn.index, dtype=object)
    if 'apl_instance' in syn:
        aside = aside.where(aside.notna(), syn['apl_instance'].map(_side))
    own = syn['partner_side'].map(_side) if 'partner_side' in syn else pd.Series(None, index=syn.index, dtype=object)
    syn = syn.assign(a=aside, o=own)

    def majority(col: str) -> pd.Series:
        s = syn.dropna(subset=[col])
        if not len(s):
            return pd.Series(dtype=object)
        c = s.groupby(['partner_id', col]).size().unstack(fill_value=0)
        win = c.idxmax(axis=1)
        if c.shape[1] > 1:
            srt = np.sort(c.to_numpy(), axis=1)
            win = win.where(srt[:, -1] > srt[:, -2])
        return win
    out['apl_side'] = majority('a').reindex(out.index)
    out['own_side'] = majority('o').reindex(out.index)
    out['n_apl_syn'] = syn.groupby('partner_id').size().reindex(out.index).fillna(0).astype(int)
    x = kx.to_numpy()
    i = int(np.argmax(np.diff(x)))
    lo = pd.Series(np.arange(len(x)) <= i, index=kx.index)
    agree_lo_L = float((out['apl_side'][lo] == 'L').mean() + (out['apl_side'][~lo] == 'R').mean())
    agree_lo_R = float((out['apl_side'][lo] == 'R').mean() + (out['apl_side'][~lo] == 'L').mean())
    lo_side, hi_side = ('L', 'R') if agree_lo_L >= agree_lo_R else ('R', 'L')
    out['gap_side'] = np.where(lo, lo_side, hi_side)
    out['side'] = out['apl_side'].where(out['apl_side'].notna(), out['gap_side'])
    _CACHE[key] = out
    return out


def derived_claws(dataset: str, roi: str = 'official') -> pd.DataFrame:
    """A dataset's derived claws with the columns kc_id, claw_id, kc_type, kc_class, n_PN_synapses, main_PN_id,
    main_bouton_id, main_PN_type, glomerulus (the model's name, or None), side, mb ('hemibrain', or '<dataset>_<side>')."""
    key = ('claws', dataset, roi)
    if key not in _CACHE:
        df = _read_claws(dataset, roi).copy()
        df['glomerulus'] = [glomerulus_name(t, dataset) for t in df['main_PN_type']]
        df['kc_class'] = [kc_class(t) for t in df['kc_type']]
        sides = kc_sides(dataset, roi)['side']
        df['side'] = df['kc_id'].map(sides)
        df['mb'] = 'hemibrain' if dataset == 'hemibrain' else dataset + '_' + df['side'].astype(str)
        _CACHE[key] = df
    return _CACHE[key]


def derived_boutons(dataset: str, roi: str = 'official') -> pd.DataFrame:
    """A dataset's derived boutons with bouton_id, pn_id, n_tbars, glomerulus (the model's name, or None) and mb (that
    of the KCs whose counted claws have the bouton's PN as their main PN, the majority; None for a PN without such a
    claw)."""
    b = pd.read_csv(ROI_VERSIONS[roi] / f'{dataset}_PN_boutons_derived.csv.gz')
    b['pn_id'] = _ids(b['pn_id'])
    b['glomerulus'] = [glomerulus_name(t, dataset) for t in b['pn_type']]
    c = derived_claws(dataset, roi)
    c = c[c.n_PN_synapses >= SYN_MIN].dropna(subset=['side'])
    pn_mb = c.groupby(['main_PN_id', 'mb']).size().unstack(fill_value=0)
    b['mb'] = b['pn_id'].map(pn_mb.idxmax(axis=1)) if len(pn_mb) else None
    return b


def mb_counts(claws: pd.DataFrame, syn_min: int = SYN_MIN, by_class: bool = False) -> pd.DataFrame:
    """Per mushroom body (and KC class, with `by_class`) and glomerulus: the claws with `syn_min` or more PN synapses
    whose main PN names one glomerulus ('n_claws'), the boutons and PNs that carry them as main bouton and PN
    ('n_boutons', 'n_pns') and their KCs ('n_kcs')."""
    c = claws[(claws['n_PN_synapses'] >= syn_min) & claws['glomerulus'].notna()]
    keys = ['mb'] + (['kc_class'] if by_class else []) + ['glomerulus']
    return (c.groupby(keys).agg(n_claws=('claw_id', 'size'), n_boutons=('main_bouton_id', 'nunique'),
                                n_pns=('main_PN_id', 'nunique'), n_kcs=('kc_id', 'nunique')).reset_index())


def shares(counts: pd.Series, glomeruli: Sequence[str]) -> pd.Series:
    """A mushroom body's share per glomerulus: its counted claws of g over its counted claws on `glomeruli`."""
    s = pd.Series(counts, dtype=float).reindex(list(glomeruli)).fillna(0.0)
    return s / s.sum()


def ratios(counts: pd.DataFrame, glomeruli: Sequence[str], ref: str = 'hemibrain', small: int = SMALL_CLAWS
           ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Each mushroom body's ratio r_g to `ref` (counts: glomeruli x mushroom bodies; shares over `glomeruli`), with
    r_g = 1 where either has fewer than `small` counted claws. Returns (r, small), both glomeruli x mushroom bodies."""
    C = counts.reindex(list(glomeruli)).fillna(0.0).astype(float)
    S = C / C.sum()
    with np.errstate(divide='ignore', invalid='ignore'):
        r = S.div(S[ref], axis=0)
    sm = C.lt(small) | pd.DataFrame(np.repeat((C[ref] < small).to_numpy()[:, None], C.shape[1], axis=1),
                                    index=C.index, columns=C.columns)
    r = r.where(~sm, 1.0)
    return r, sm


def mean_proportions(r: pd.DataFrame, brains: Mapping[str, Sequence[str]] = BRAINS) -> pd.Series:
    """The geometric mean of r_g over the three brains: the hemibrain's 1 and each other brain's sides averaged on log2
    first."""
    logs = [pd.Series(0.0, index=r.index)]
    for _, mbs in brains.items():
        logs.append(np.log2(r[list(mbs)].astype(float)).mean(axis=1))
    return pd.Series(2.0 ** (sum(logs) / len(logs)), index=r.index, name='r_mean')


def pratyush_flywire() -> pd.DataFrame:
    """FlyWire's PN->KC connections from olfsysm's tables (pyolfsysm/data/from_pratyush/2024-09-13, both files), with the
    loader's glomerulus names (from source_hemibrain_type: '+' and M dropped, VC3l and VC3m merged into VC3), the KC's side
    (target_side) and KC class; connections of PRATYUSH_MIN_WEIGHT or more synapses (the loader's threshold)."""
    import pyolfsysm
    d = Path(pyolfsysm.__file__).resolve().parent / 'data' / FLYWIRE_PRATYUSH
    df = pd.concat([pd.read_csv(d / f'FlyWire_PNKC_{s}.csv') for s in ('Left', 'Right')], ignore_index=True)
    df['glomerulus'] = [glomerulus_name(t, 'hemibrain') for t in df['source_hemibrain_type']]
    df['side'] = df['target_side'].map(_side)
    df['mb'] = 'flywire_' + df['side'].astype(str)
    df['kc_class'] = [kc_class(t) for t in df['target_hemibrain_type']]
    return df[df['weight'] >= PRATYUSH_MIN_WEIGHT].copy()


# ---------------------------------------------------------------------------
# the model's wirings, the builder and the controls

def claw_table(W: pd.DataFrame) -> pd.DataFrame:
    """Per row of a claw x bouton wiring (one row per claw, 0 / 1 with at most one 1): kc_id, claw_id, kc_type, kc_class,
    col (the bouton's column position, -1 for a row without one) and glomerulus (None for a row without a bouton)."""
    v = W.to_numpy()
    has = v.astype(bool).any(axis=1)
    col = np.where(has, v.astype(float).argmax(axis=1), -1)
    gl = np.asarray([str(x) for x in W.columns.get_level_values('glomerulus')], dtype=object)
    t = pd.DataFrame({H.KC_ID_LEVEL: np.asarray(W.index.get_level_values(H.KC_ID_LEVEL)),
                      'claw_id': np.asarray(W.index.get_level_values('claw_id')),
                      H.KC_TYPE_LEVEL: np.asarray(W.index.get_level_values(H.KC_TYPE_LEVEL)).astype(str),
                      'col': col.astype(int)})
    t['glomerulus'] = np.where(has, gl[np.maximum(col, 0)], None)
    t['kc_class'] = [kc_class(x) for x in t[H.KC_TYPE_LEVEL]]
    return t


def _check_format(W: pd.DataFrame) -> np.ndarray:
    v = W.to_numpy()
    if W.index.has_duplicates or W.columns.has_duplicates:
        raise ValueError('the wiring has duplicate rows or columns')
    if not (np.isin(v, (0, 1)).all() and (v.astype(int).sum(axis=1) <= 1).all()):
        raise ValueError('a claw x bouton wiring is 0 / 1 with at most one 1 per row')
    for lev in (H.KC_ID_LEVEL, 'claw_id', H.KC_TYPE_LEVEL):
        if lev not in W.index.names:
            raise ValueError(f'the wiring index has no {lev!r} level')
    if 'glomerulus' not in W.columns.names:
        raise ValueError("the wiring's columns have no 'glomerulus' level")
    return v


def type_glomerulus_counts(W: pd.DataFrame) -> pd.DataFrame:
    """KC type x glomerulus claw counts of a wiring (every glomerulus of its columns)."""
    t = claw_table(W)
    t = t[t.col >= 0]
    gl = list(dict.fromkeys(str(x) for x in W.columns.get_level_values('glomerulus')))
    c = pd.crosstab(t[H.KC_TYPE_LEVEL], t['glomerulus'])
    return c.reindex(columns=gl, fill_value=0).astype(np.int64)


def type_targets(counts: pd.DataFrame, r: Mapping[str, float]) -> pd.DataFrame:
    """Per KC type (rows) and glomerulus (columns), the target claw counts: each count times r_g (1 where `r` has none),
    rescaled to the type's total and rounded by largest remainder (ties by glomerulus name; only glomeruli with claws in
    the type receive the remainder), so that every type's total is unchanged."""
    rr = pd.Series(dict(r), dtype=float).reindex(counts.columns).fillna(1.0).to_numpy()
    if not (np.isfinite(rr).all() and (rr > 0).all()):
        raise ValueError('ratios must be finite and positive')
    out = pd.DataFrame(0, index=counts.index, columns=counts.columns, dtype=np.int64)
    names = np.asarray([str(c) for c in counts.columns])
    for t in counts.index:
        n = counts.loc[t].to_numpy(dtype=float)
        tot = int(round(n.sum()))
        if tot == 0:
            continue
        x = n * rr
        y = x * (tot / x.sum())
        f = np.floor(y + 1e-9)
        rem = np.clip(y - f, 0.0, None)
        left = int(tot - f.sum())
        cand = np.flatnonzero(n > 0)
        order = sorted(cand, key=lambda i: (-rem[i], names[i]))
        for i in order[:left]:
            f[i] += 1
        if int(f.sum()) != tot:
            raise RuntimeError(f'type {t}: targets sum to {f.sum()}, not {tot}')
        out.loc[t] = f.astype(np.int64)
    return out


def move_claws(W: pd.DataFrame, r: Mapping[str, float], seed) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The wiring `W` (claw x bouton, one row per claw) at proportions `r` (module docstring), the draws from
    `np.random.default_rng(seed)`.

    Returns:
        (the built matrix, with W's index, columns and dtypes; the moved claws, one row each: row (position), kc_id,
        claw_id, kc_type, from_glomerulus, from_col, from_bouton, to_glomerulus, to_col, to_bouton; the targets, one row
        per KC type and glomerulus: kc_type, glomerulus, r, base, target)."""
    v = _check_format(W)
    t = claw_table(W)
    has = t.col.to_numpy() >= 0
    types = t[H.KC_TYPE_LEVEL].to_numpy()
    glom_row = t.glomerulus.to_numpy()
    col = t.col.to_numpy()
    gl_col = np.asarray([str(x) for x in W.columns.get_level_values('glomerulus')])
    counts = type_glomerulus_counts(W)
    targets = type_targets(counts, r)
    rng = np.random.default_rng(seed)
    new_col = col.copy()
    moved = []
    bc_all = np.bincount(col[has], minlength=W.shape[1]).astype(float)
    for ty in sorted(counts.index):
        diff = (targets.loc[ty] - counts.loc[ty])
        gnames = sorted(diff.index)
        lost: List[int] = []
        for g in gnames:
            if diff[g] < 0:
                rows = np.flatnonzero(has & (types == ty) & (glom_row == g))
                lost.extend(np.sort(rng.choice(rows, size=int(-diff[g]), replace=False)).tolist())
        slots = [g for g in gnames for _ in range(int(max(diff[g], 0)))]
        if len(slots) != len(lost):
            raise RuntimeError(f'type {ty}: {len(lost)} claws lost, {len(slots)} places')
        if not lost:
            continue
        bc = np.bincount(col[has & (types == ty)], minlength=W.shape[1]).astype(float)
        order = rng.permutation(len(lost))
        for i, k in enumerate(order):
            row, g = lost[k], slots[i]
            cand = np.flatnonzero(gl_col == g)
            w = bc[cand]
            if w.sum() <= 0:
                w = bc_all[cand]
            if w.sum() <= 0:
                w = np.ones(len(cand))
            c = int(cand[rng.choice(len(cand), p=w / w.sum())])
            new_col[row] = c
            moved.append(dict(row=int(row), kc_id=t[H.KC_ID_LEVEL].iat[row], claw_id=t['claw_id'].iat[row], kc_type=ty,
                              from_glomerulus=glom_row[row], from_col=int(col[row]), from_bouton=str(W.columns[col[row]]),
                              to_glomerulus=g, to_col=c, to_bouton=str(W.columns[c])))
    nv = v.copy()
    if moved:
        rows = np.asarray([m['row'] for m in moved])
        nv[rows, :] = 0
        nv[rows, new_col[rows]] = 1
    W1 = pd.DataFrame(nv, index=W.index.copy(), columns=W.columns.copy())
    if W.dtypes.nunique() > 1:
        W1 = W1.astype(W.dtypes.to_dict())
    tg = targets.stack().rename('target').reset_index()
    tg.columns = [H.KC_TYPE_LEVEL, 'glomerulus', 'target']
    base = counts.stack().rename('base').reset_index()
    base.columns = [H.KC_TYPE_LEVEL, 'glomerulus', 'base']
    tg = base.merge(tg, on=[H.KC_TYPE_LEVEL, 'glomerulus'])
    tg.insert(2, 'r', [float(dict(r).get(g, 1.0)) for g in tg.glomerulus])
    mv = pd.DataFrame(moved, columns=['row', H.KC_ID_LEVEL, 'claw_id', H.KC_TYPE_LEVEL, 'from_glomerulus', 'from_col',
                                      'from_bouton', 'to_glomerulus', 'to_col', 'to_bouton'])
    return W1, mv, tg


def check_built(W0: pd.DataFrame, W1: pd.DataFrame, moved: pd.DataFrame, targets: pd.DataFrame) -> dict:
    """Part 0's invariants of a built wiring: the same index and columns; each row's claw count (0 or 1), so each KC's
    claw count, and its type the base's; the per-type totals unchanged; the per-type glomerulus counts equal to the
    targets; every moved claw on a bouton of its new glomerulus, from a bouton of its old one; every other row unchanged.
    Returns the checks with 'passes'."""
    res = dict(same_index=bool(W0.index.equals(W1.index)), same_columns=bool(W0.columns.equals(W1.columns)))
    if not (res['same_index'] and res['same_columns']):
        res['passes'] = False
        return res
    v0, v1 = W0.to_numpy().astype(int), W1.to_numpy().astype(int)
    res['row_claws_same'] = bool(np.array_equal(v0.sum(axis=1), v1.sum(axis=1)))
    kc = np.asarray(W0.index.get_level_values(H.KC_ID_LEVEL))
    res['kc_claws_same'] = bool(pd.Series(v0.sum(axis=1)).groupby(kc).sum().equals(pd.Series(v1.sum(axis=1)).groupby(kc).sum()))
    c0, c1 = type_glomerulus_counts(W0), type_glomerulus_counts(W1)
    res['type_totals_same'] = bool(c0.sum(axis=1).equals(c1.sum(axis=1)))
    tg = targets.pivot(index=H.KC_TYPE_LEVEL, columns='glomerulus', values='target').reindex(
        index=c1.index, columns=c1.columns).fillna(0).astype(np.int64)
    res['targets_met'] = bool(np.array_equal(tg.to_numpy(), c1.to_numpy()))
    mrows = moved['row'].to_numpy(dtype=int) if len(moved) else np.zeros(0, dtype=int)
    other = np.setdiff1d(np.arange(v0.shape[0]), mrows)
    res['others_unchanged'] = bool(np.array_equal(v0[other], v1[other]))
    gl = np.asarray([str(x) for x in W0.columns.get_level_values('glomerulus')])
    if len(moved):
        new = v1[mrows].argmax(axis=1)
        old = v0[mrows].argmax(axis=1)
        res['moved_on_new_glomerulus'] = bool((gl[new] == moved['to_glomerulus'].astype(str).to_numpy()).all()
                                              and (new == moved['to_col'].to_numpy()).all())
        res['moved_from_old_glomerulus'] = bool((gl[old] == moved['from_glomerulus'].astype(str).to_numpy()).all())
        res['moved_changed_glomerulus'] = bool((moved['to_glomerulus'] != moved['from_glomerulus']).all())
        res['moved_once'] = bool(not pd.Series(mrows).duplicated().any())
    else:
        res.update(moved_on_new_glomerulus=True, moved_from_old_glomerulus=True, moved_changed_glomerulus=True,
                   moved_once=True)
    res['n_moved'] = int(len(moved))
    res['passes'] = bool(all(v for k, v in res.items() if k != 'n_moved'))
    return res


def strata(claw_totals: pd.Series, with_input: Sequence[str], n: int = N_TERCILES) -> pd.Series:
    """B3's strata per glomerulus: 'input' / 'noinput' (ORN input in the lab's input or not) and the third of the claw
    total within that class ('_1' the fewest claws; ties by name)."""
    wi = set(with_input)
    out = {}
    for cls, gs in (('input', [g for g in claw_totals.index if g in wi]),
                    ('noinput', [g for g in claw_totals.index if g not in wi])):
        order = sorted(gs, key=lambda g: (float(claw_totals[g]), str(g)))
        for k, part in enumerate(np.array_split(np.asarray(order, dtype=object), n)):
            for g in part:
                out[g] = f'{cls}_{k + 1}'
    return pd.Series(out).reindex(claw_totals.index)


def permuted_ratios(r: pd.Series, strata_: pd.Series, seed: int, fixed: Sequence[str] = ()) -> pd.Series:
    """`r` permuted among the glomeruli within each stratum (one generator from `seed`, strata in sorted order,
    glomeruli by name); the glomeruli in `fixed` (the small ones, r_g = 1) keep theirs."""
    rng = np.random.default_rng(seed)
    out = r.astype(float).copy()
    fx = set(fixed)
    st = strata_.reindex(r.index)
    for s in sorted(st.dropna().unique()):
        gs = sorted(g for g in st.index[st == s] if g not in fx)
        if len(gs) > 1:
            out[gs] = r[gs].to_numpy(dtype=float)[rng.permutation(len(gs))]
    return out


def scaled_ratios(r: pd.Series, scale: float) -> pd.Series:
    """r ** scale (B5's half and double of the log ratios; -1 the inverted ratios of B3)."""
    return np.power(r.astype(float), float(scale))


# ---------------------------------------------------------------------------
# the variants

def arm_variants(arm: str) -> List[str]:
    """The variants of an arm (B1..B6, C2), in run order."""
    return {'B1': [f'B1_{m}' for m in OTHER_MBS],
            'B2': ['B2_mean'],
            'B3': [f'B3_permmean_s{s}' for s in PERM_SEEDS] + ['B3_invmean'] + [f'B3_perm_{m}_s{PERM_SEEDS[0]}' for m in OTHER_MBS],
            'B4': [f'B4_mean_d{DRAW_SEEDS[1]}'],
            'B5': [f'B5_mean_x{s:g}' for s in B5_SCALES],
            'B6': [f'B6_mean_{B6_INPUT}'],
            'C2': [f'C2_mean_d{DRAW_SEEDS[1]}']}[arm]


def parse_variant(v: str) -> dict:
    """A variant's fields: arm, proportions ('mean' or a mushroom body; 'hemibrain' for 'ref', the wirings as they are),
    control ('perm', 'inv' or None), seed (the permutation's), scale (of the log ratios), draw (the moved claws' draw: 1,
    or 2 for B4 / C2), input ('lab' or 'F')."""
    import re
    out = dict(variant=v, arm=None, proportions=None, control=None, seed=None, scale=1.0, draw=1, input='lab')
    if v == 'ref':
        out.update(arm='ref', proportions='hemibrain', draw=0)
        return out
    m = re.match(r'^(B[1-6]|C2)_(.+)$', v)
    if not m:
        raise ValueError(f'variant {v!r}: not one of this experiment\'s')
    arm, rest = m.groups()
    out['arm'] = arm
    if v not in arm_variants(arm):
        raise ValueError(f'variant {v!r}: not one of arm {arm}\'s')
    if arm == 'B1':
        out['proportions'] = rest
    elif arm == 'B2':
        out['proportions'] = 'mean'
    elif arm == 'B3':
        if rest == 'invmean':
            out.update(proportions='mean', control='inv', scale=-1.0)
        elif rest.startswith('permmean_s'):
            out.update(proportions='mean', control='perm', seed=int(rest.split('_s')[1]))
        else:
            mb, s = re.match(r'^perm_(.+)_s(\d+)$', rest).groups()
            out.update(proportions=mb, control='perm', seed=int(s))
    elif arm in ('B4', 'C2'):
        out.update(proportions='mean', draw=int(rest.split('_d')[1]))
    elif arm == 'B5':
        out.update(proportions='mean', scale=float(rest.split('_x')[1]))
    elif arm == 'B6':
        out.update(proportions='mean', input=rest.split('_')[1])
    return out


def run_tag(variant: str, wiring: str) -> str:
    return f'{variant}_{wiring}'


def plan(arms: Sequence[str] = ARMS, wirings: Sequence[str] = ALL_WIRINGS) -> List[Tuple[str, str]]:
    """[(variant, wiring)] of `arms` on `wirings` (each arm on its own wirings, ARMS_ON), wiring by wiring in
    ALL_WIRINGS order, each wiring's arms in ARM_ORDER (C2 last)."""
    out = []
    order = [a for a in ARM_ORDER if a in arms] + [a for a in arms if a not in ARM_ORDER]
    for w in ALL_WIRINGS:
        if w not in wirings:
            continue
        for a in order:
            if w in ARMS_ON[a]:
                out += [(v, w) for v in arm_variants(a)]
    return out


def ratios_table(path: Optional[Path] = None) -> pd.DataFrame:
    """Part A's ratios.csv (glomerulus as index)."""
    p = out_dir() / 'ratios.csv' if path is None else path
    if not p.exists():
        raise RuntimeError(f'{p} does not exist: run --part A first')
    key = ('ratios', str(p), p.stat().st_mtime)
    if key not in _CACHE:
        _CACHE[key] = pd.read_csv(p).set_index('glomerulus')
    return _CACHE[key]


def variant_ratios(variant: str, rt: Optional[pd.DataFrame] = None) -> pd.Series:
    """The ratios r_g of a variant over the model's glomeruli (ratios.csv of part A): a mushroom body's ('r_<mb>'), the
    mean proportions ('r_mean'), permuted within B3's strata (the small glomeruli fixed), inverted, or scaled."""
    rt = ratios_table() if rt is None else rt
    pv = parse_variant(variant)
    if pv['proportions'] == 'hemibrain':
        return pd.Series(1.0, index=rt.index)
    col = 'r_mean' if pv['proportions'] == 'mean' else f"r_{pv['proportions']}"
    r = rt[col].astype(float)
    if pv['control'] == 'perm':
        fixed = list(rt.index[rt[f"fixed_{pv['proportions']}"].astype(bool)])
        r = permuted_ratios(r, rt['stratum'], pv['seed'], fixed)
    if pv['scale'] != 1.0:
        r = scaled_ratios(r, pv['scale'])
    return r.rename(variant)


def draw_seed(variant: str, wiring: str) -> List[int]:
    """The moved claws' draw: `np.random.default_rng([DRAW_SEEDS[draw - 1], position of the wiring in ALL_WIRINGS])`."""
    return [int(DRAW_SEEDS[max(parse_variant(variant)['draw'], 1) - 1]), int(ALL_WIRINGS.index(wiring))]


def wiring_tag_of(W: pd.DataFrame) -> str:
    """The name of a built matrix for E10_joint_fit.bouton_base's cache: 'E97-' and the first 12 hex digits of its
    sha256 (the same matrix, the same cached APL-off breadth run)."""
    return f'E97-{_sha(W)[:12]}'


def built_wiring(variant: str, wiring: str, rt: Optional[pd.DataFrame] = None
                 ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    """(matrix, moved claws, targets, wiring tag) of a variant on a wiring: `move_claws` of the wiring's own matrix
    (`base_wiring`) at the variant's ratios and draw."""
    W0 = base_wiring(wiring)
    W1, mv, tg = move_claws(W0, variant_ratios(variant, rt).to_dict(), draw_seed(variant, wiring))
    return W1, mv, tg, wiring_tag_of(W1)


# ---------------------------------------------------------------------------
# the base wirings (cached) and the held reading

def _w_dir() -> Path:
    return local_dir('W')


def base_wiring(label: str, simulate: bool = True) -> pd.DataFrame:
    """A wiring's claw x bouton matrix: part 0's cache (outputs/.../W/<label>_W.parquet), else E83's `wiring_parts`
    (one simulation of the bouton base) and cached, with its claw table (<label>_claws.parquet)."""
    key = ('W', label)
    if key in _CACHE:
        return _CACHE[key]
    p = _w_dir() / f'{label}_W.parquet'
    if p.exists():
        W = pd.read_parquet(p)
    elif simulate:
        W = E83.wiring_parts(label)['W']
        W.to_parquet(p)
        claw_table(W).to_parquet(_w_dir() / f'{label}_claws.parquet')
        if not pd.read_parquet(p).equals(W):
            raise RuntimeError(f'{label}: the cached matrix does not read back equal')
    else:
        raise RuntimeError(f'{label}: no cached matrix (run --part 0 first: it caches them under {_w_dir()})')
    _CACHE[key] = W
    return W


def model_claws(label: str = 'real') -> pd.DataFrame:
    """A wiring's claw table from part 0's cache (no simulation)."""
    p = _w_dir() / f'{label}_claws.parquet'
    if not p.exists():
        return claw_table(base_wiring(label, simulate=False))
    return pd.read_parquet(p)


def build97(label: str, W: pd.DataFrame, wiring_tag: str, orn: Optional[pd.DataFrame] = None,
            input_tag: Optional[str] = None):
    """C2's family on wiring `label`'s fit with the claw x bouton matrix `W` (the wiring's index and columns) and its
    breadth recomputed (`E10_joint_fit.bouton_base(wiring_matrix=W, wiring_tag=...)`, the APL-off breadth run cached
    under outputs/.../breadth/; with another input `orn`, `input_tag` names it there too); the normalized APL weights
    must be the wiring's. The simulator as E83's `build_simulator` builds it. Returns (simulate, info)."""
    from experiments import E10_joint_fit as J
    from experiments.E77_input_screen import _scalar_kws
    from pyolfsysm.mb_model import make_fixed_param_kc_simulator
    P = E83.wiring_parts(label)
    if not (W.index.equals(P['W'].index) and W.columns.equals(P['W'].columns)):
        raise ValueError(f'{label}: the matrix does not have the wiring\'s index and columns')
    kw = dict(wiring_matrix=W, wiring_tag=wiring_tag)
    if orn is None:
        orn = P['orn']
    else:
        if input_tag is None:
            raise ValueError('another ORN input needs input_tag (it names the cached breadth run)')
        kw['input_tag'] = input_tag
    _, wnorm, W2, breadth = J.bouton_base(orn, local_dir('breadth'), None, **kw)
    if not W2.equals(W):
        raise RuntimeError(f'{label}: bouton_base returned another matrix')
    for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL'):
        if not P['wnorm'][k].equals(wnorm[k]):
            raise RuntimeError(f'{label}: {k} built with {wiring_tag} differs from the wiring\'s')
    p = dict(P['p'])
    wAPLKC = J.breadth_weights(wnorm['wAPLKC'], breadth, E83.BREADTH_K)
    sim, _ = make_fixed_param_kc_simulator(orn, initial_params=p, **J.base_kws(None), _wPNKC=W, _wAPLKC=wAPLKC,
                                           _wKCAPL=wnorm['wKCAPL'], _wAPLPN=wnorm['wAPLPN'], _wPNAPL=wnorm['wPNAPL'],
                                           kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU)
    cells = dict(zip(sim.__code__.co_freevars, sim.__closure__))
    kws = _scalar_kws(dict(J.base_kws(None), kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU))
    return sim, dict(label=label, params=p, kws=kws, breadth=breadth, wAPLKC=wAPLKC, wAPLPN=wnorm['wAPLPN'], W=W,
                     orn=orn, cells=cells)


def held97(label: str, W: pd.DataFrame, wiring_tag: str, orn: Optional[pd.DataFrame] = None,
           input_tag: Optional[str] = None) -> dict:
    """The held reading of `build97`'s simulator: E77's `hold_rates` at the wiring's fit's response counts and parameters
    (the somatic threshold bisected with APL off, then the tied APL scale with APL on). Returns hold_rates' dict with
    'wiring_tag', 'input_tag' and 'breadth_equal_fit' (the recomputed breadth equal to the wiring's fit's)."""
    from experiments import E77_input_screen as S
    P = E83.saved_fit(label)
    keep = {}

    def build():
        sim, info = build97(label, W, wiring_tag, orn, input_tag)
        b0 = E83.wiring_parts(label)['breadth']
        keep['breadth_equal_fit'] = bool(info['breadth'].index.equals(b0.index)
                                         and np.array_equal(info['breadth'].to_numpy(dtype=float), b0.to_numpy(dtype=float)))
        return sim, info
    res = S.hold_rates(build, target_off=P['target_off'], target_on=P['target_on'], params=P['p'],
                       label=f'{label} {wiring_tag}')
    res.update(wiring_tag=wiring_tag, input_tag=input_tag, **keep)
    return res


# ---------------------------------------------------------------------------
# part 0: the identities

def identity_passed() -> bool:
    return bool(_read_json(out_dir() / 'identity.json').get('passes'))


def built_records() -> dict:
    return _read_json(out_dir() / 'identity.json').get('built', {})


def part_0(wirings: Sequence[str] = ALL_WIRINGS, held: bool = True) -> dict:
    """The identities (module docstring, part 0); identity.json (resumable per key); raises if one fails."""
    import gc
    path = out_dir() / 'identity.json'
    res = _read_json(path)
    res['commits'] = _commits_clean()
    res['host'] = _host()
    for w in wirings:
        key = f'unchanged_{w}'
        if res.get(key, {}).get('passes'):
            continue
        t0 = time.time()
        W = base_wiring(w)
        W1, mv, tg = move_claws(W, {}, draw_seed('B2_mean', w))
        t = claw_table(W)
        chk = check_built(W, W1, mv, tg)
        res[key] = dict(passes=bool(W1.equals(W) and len(mv) == 0 and chk['passes']), n_rows=int(W.shape[0]),
                        n_columns=int(W.shape[1]), n_claws=int((t.col >= 0).sum()), n_kcs=int(t[H.KC_ID_LEVEL].nunique()),
                        n_glomeruli=int(t.glomerulus.dropna().nunique()), types=sorted(t[H.KC_TYPE_LEVEL].unique()),
                        column_levels=[str(x) for x in W.columns.names], index_levels=[str(x) for x in W.index.names],
                        dtype=str(W.dtypes.iloc[0]), sha256=_sha(W), seconds=round(time.time() - t0, 1))
        print(f'E97 0 {key}: {res[key]}', flush=True)
        _write_json(path, res)
    if held and 'real' in wirings and not res.get('held_real', {}).get('passes'):
        t0 = time.time()
        W = base_wiring('real')
        r = held97('real', W, wiring_tag_of(W))
        ns, nr = E83._same({'on': r['on'], 'off': r['off']}, E83.saved_fit('real')['saved'])
        res['held_real'] = dict(count_diffs=ns, response_diffs=nr, thr_factor=r['thr_factor'], apl_factor=r['apl_factor'],
                                n_on=r['n_on'], n_off=r['n_off'], target_on=r['target_on'], target_off=r['target_off'],
                                n_sims=r['n_sims'], breadth_equal_fit=r.get('breadth_equal_fit'),
                                wiring_tag=r['wiring_tag'], seconds=round(time.time() - t0, 1),
                                passes=bool(ns == 0 and nr == 0))
        print(f"E97 0 held_real: {res['held_real']}", flush=True)
        _write_json(path, res)
        del r
        gc.collect()
    if (out_dir() / 'ratios.csv').exists():
        rt = ratios_table()
        built = res.setdefault('built', {})
        for v, w in plan(ARMS + ('C2',), wirings):
            key = run_tag(v, w)
            if built.get(key, {}).get('passes'):
                continue
            W0 = base_wiring(w)
            W1, mv, tg, wtag = built_wiring(v, w, rt)
            chk = check_built(W0, W1, mv, tg)
            n = int((W0.to_numpy().astype(int).sum(axis=1) > 0).sum())
            built[key] = dict(chk, sha256=_sha(W1), wiring_tag=wtag, moved_share=len(mv) / n,
                              ratios_sha256=_sha(variant_ratios(v, rt).to_frame()))
            if not chk['passes']:
                print(f'E97 0 built {key}: FAILS {chk}', flush=True)
        _write_json(path, res)
        nb = len(built)
        print(f"E97 0 built: {sum(b['passes'] for b in built.values())} of {nb} pass; moved share "
              f"{min(b['moved_share'] for b in built.values()):.3f} to {max(b['moved_share'] for b in built.values()):.3f}",
              flush=True)
    else:
        print('E97 0: no ratios.csv yet (run --part A, then --part 0 again for the built wirings)', flush=True)
    keys = [f'unchanged_{w}' for w in ALL_WIRINGS] + ['held_real']
    failed = [k for k in keys if k in res and not res[k].get('passes')]
    failed += [f'built:{k}' for k, b in res.get('built', {}).items() if not b.get('passes')]
    n_plan = len(plan(ARMS + ('C2',), ALL_WIRINGS))
    complete = bool(all(k in res for k in keys) and len(res.get('built', {})) >= n_plan)
    res.update(failed=failed, complete=complete, passes=bool(complete and not failed),
               passes_so_far=bool(not failed))
    _write_json(path, res)
    if failed:
        raise RuntimeError(f'E97 part 0: the identity fails for {failed} (identity.json); the experiment stops here')
    print(f"E97 0: {'passes' if res['passes'] else 'so far exact (not all keys yet)'}", flush=True)
    return res


# ---------------------------------------------------------------------------
# part A (no simulations)

def model_glomeruli(claws: Optional[pd.DataFrame] = None) -> List[str]:
    """The model's glomeruli (the real wiring's columns, in their order; part 0's cache)."""
    p = _w_dir() / 'real_W.parquet'
    if ('W', 'real') in _CACHE or p.exists():
        W = base_wiring('real', simulate=False)
        return list(dict.fromkeys(str(x) for x in W.columns.get_level_values('glomerulus')))
    c = model_claws('real') if claws is None else claws
    return sorted(c.glomerulus.dropna().unique())


def a1_rule(model_counts: pd.Series, derived_counts: pd.Series, min_claws: int = A1_RULE_MIN_CLAWS,
            max_sd: float = A1_RULE_SD) -> dict:
    """A1's rule: the SD over the glomeruli with `min_claws` or more claws in the model of log2(derived share / model
    share), each share over all the model's glomeruli; 'stops' if above `max_sd` (or a glomerulus has none)."""
    gl = list(model_counts.index)
    sm, sd_ = shares(model_counts, gl), shares(derived_counts, gl)
    keep = model_counts.reindex(gl) >= min_claws
    with np.errstate(divide='ignore'):
        lr = np.log2(sd_[keep] / sm[keep])
    finite = np.isfinite(lr)
    sd = float(lr[finite].std(ddof=1)) if finite.sum() > 1 else np.nan
    worst = lr[finite].abs().sort_values(ascending=False).head(5)
    return dict(sd=sd, n_glomeruli=int(keep.sum()), n_infinite=int((~finite).sum()), max_sd=max_sd, min_claws=min_claws,
                mean_log2=float(lr[finite].mean()), largest={g: float(lr[g]) for g in worst.index},
                stops=bool(not np.isfinite(sd) or sd > max_sd or (~finite).any()))


def a3_rule(agreement: float, null: np.ndarray, q: Tuple[float, float] = A3_QUANTILES) -> str:
    """A3's rule (as E93's A2): 'with' above the null's upper quantile, 'against' below its lower, otherwise 'neither'."""
    lo, hi = np.nanquantile(np.asarray(null, dtype=float), q)
    if agreement > hi:
        return 'with'
    if agreement < lo:
        return 'against'
    return 'neither'


def sign_agreement(x: np.ndarray, y: np.ndarray) -> float:
    """The share of glomeruli whose two values have the same sign (zeros count as disagreement; E93's)."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.size == 0:
        return np.nan
    return float(np.mean((np.sign(x) == np.sign(y)) & (np.sign(x) != 0)))


def _spearman(x, y) -> float:
    from scipy.stats import spearmanr
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.size < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(spearmanr(x, y)[0])


def counting_noise_sd(n_a, n_b, n_draws: int = NOISE_DRAWS, seed: int = NOISE_SEED) -> dict:
    """The SD over glomeruli of log2(share_a / share_b) expected from counting alone: analytic (sqrt of the mean over
    glomeruli of (1 / n_a + 1 / n_b), over ln 2) and from multinomial draws of both bodies' counts at their pooled shares
    and totals (the mean SD over `n_draws`). Glomeruli with no claws in either are left out."""
    a, b = np.asarray(n_a, dtype=float), np.asarray(n_b, dtype=float)
    ok = (a > 0) & (b > 0)
    a, b = a[ok], b[ok]
    analytic = float(np.sqrt(np.mean(1.0 / a + 1.0 / b)) / np.log(2))
    p = (a + b) / (a + b).sum()
    rng = np.random.default_rng(seed)
    sds = []
    for _ in range(n_draws):
        xa = rng.multinomial(int(a.sum()), p).astype(float)
        xb = rng.multinomial(int(b.sum()), p).astype(float)
        m = (xa > 0) & (xb > 0)
        sds.append(np.std(np.log2((xa[m] / xa.sum()) / (xb[m] / xb.sum())), ddof=1))
    return dict(analytic=analytic, simulated=float(np.mean(sds)), n_glomeruli=int(ok.sum()))


def a1_counts(out: Path) -> dict:
    """A1 (module docstring): writes a1_counts.csv, a1_model.csv, a1_hemibrain_vs_model.csv, a1_roi.csv, a1_sides.csv,
    a1_flywire_pratyush.csv; returns the rule, the ROI choice per mushroom body, the chosen counts (glomeruli x
    mushroom bodies at SYN_MIN, each body in the version it enters with) and the model's counts."""
    mc = model_claws('real')
    gl = model_glomeruli(mc)
    mcl = mc[mc.col >= 0]
    model = mcl.groupby('glomerulus').agg(n_claws=('claw_id', 'size'), n_boutons=('col', 'nunique'),
                                          n_kcs=(H.KC_ID_LEVEL, 'nunique')).reindex(gl).fillna(0)
    for cls in KC_CLASSES:
        model[f'n_claws_{cls}'] = mcl[mcl.kc_class == cls].groupby('glomerulus').size().reindex(gl).fillna(0)
    model.index.name = 'glomerulus'
    model.reset_index().to_csv(out / 'a1_model.csv', index=False)
    rows, side_rows = [], []
    counts5: Dict[Tuple[str, str], pd.Series] = {}
    for ds in DATASETS:
        for roi in ROI_VERSIONS:
            c = derived_claws(ds, roi)
            sd = kc_sides(ds, roi)
            for s in MB_SIDES[ds]:
                k = sd[sd.side == s]
                both = k.apl_side.notna() & k.own_side.notna()
                side_rows.append(dict(dataset=ds, roi=roi, side=s, n_kcs=int(len(k)), by_apl=int(k.apl_side.notna().sum()),
                                      agree_own=float((k.own_side == k.apl_side)[both].mean()) if both.any() else np.nan,
                                      agree_gap=float((k.gap_side == k.apl_side)[k.apl_side.notna()].mean())
                                      if k.apl_side.notna().any() else np.nan))
            for thr in SYN_THRESHOLDS:
                for by_class in (False, True):
                    t = mb_counts(c, thr, by_class)
                    if not by_class:
                        t.insert(1, 'kc_class', 'all')
                        if thr == SYN_MIN:
                            for mb, g in t.groupby('mb'):
                                counts5[(mb, roi)] = g.set_index('glomerulus')['n_claws']
                    t.insert(0, 'threshold', thr)
                    t.insert(0, 'roi', roi)
                    t.insert(0, 'dataset', ds)
                    rows.append(t)
            # the boutons and PNs of the bouton table, by the side of their PN's counted claws
            b = derived_boutons(ds, roi)
            bt = b.dropna(subset=['glomerulus', 'mb']).groupby(['mb', 'glomerulus']).agg(
                n_boutons=('bouton_id', 'size'), n_pns=('pn_id', 'nunique')).reset_index()
            bt.insert(0, 'kc_class', 'bouton_table')
            bt.insert(0, 'threshold', SYN_MIN)
            bt.insert(0, 'roi', roi)
            bt.insert(0, 'dataset', ds)
            rows.append(bt)
    a1 = pd.concat(rows, ignore_index=True, sort=False)
    a1['in_model'] = a1.glomerulus.isin(gl)
    a1.to_csv(out / 'a1_counts.csv', index=False)
    pd.DataFrame(side_rows).to_csv(out / 'a1_sides.csv', index=False)
    # the derived hemibrain against the model's wiring
    hv = pd.DataFrame(index=pd.Index(gl, name='glomerulus'))
    hv['model_claws'] = model['n_claws']
    hv['model_share'] = shares(model['n_claws'], gl)
    rules = {}
    for roi in ROI_VERSIONS:
        for thr in SYN_THRESHOLDS:
            t = a1[(a1.dataset == 'hemibrain') & (a1.roi == roi) & (a1.threshold == thr) & (a1.kc_class == 'all')]
            n = t.set_index('glomerulus')['n_claws'].reindex(gl).fillna(0)
            hv[f'derived_claws_{roi}_thr{thr}'] = n
            hv[f'derived_share_{roi}_thr{thr}'] = shares(n, gl)
            with np.errstate(divide='ignore'):
                hv[f'log2_ratio_{roi}_thr{thr}'] = np.log2(hv[f'derived_share_{roi}_thr{thr}'] / hv['model_share'])
            rules[f'{roi}_thr{thr}'] = a1_rule(model['n_claws'], n)
    for cls in KC_CLASSES:
        t = a1[(a1.dataset == 'hemibrain') & (a1.roi == 'official') & (a1.threshold == SYN_MIN) & (a1.kc_class == cls)]
        hv[f'derived_claws_official_thr{SYN_MIN}_{cls}'] = t.set_index('glomerulus')['n_claws'].reindex(gl).fillna(0)
        hv[f'model_claws_{cls}'] = model[f'n_claws_{cls}']
    hv.reset_index().to_csv(out / 'a1_hemibrain_vs_model.csv', index=False)
    # the two ROI versions per mushroom body: the version a body enters with
    keep = model['n_claws'] >= A1_RULE_MIN_CLAWS
    roi_rows, choice = [], {}
    for mb in ALL_MBS:
        a, b = counts5.get((mb, 'official')), counts5.get((mb, 'outside'))
        if a is None or b is None:
            continue
        sa, sb = shares(a, gl), shares(b, gl)
        with np.errstate(divide='ignore', invalid='ignore'):
            lr = np.log2(sb[keep] / sa[keep])
        lr = lr[np.isfinite(lr)]
        sd = float(lr.std(ddof=1))
        choice[mb] = 'outside' if sd > A1_RULE_SD else 'official'
        roi_rows.append(dict(mb=mb, n_official=float(a.reindex(gl).fillna(0).sum()),
                             n_outside=float(b.reindex(gl).fillna(0).sum()), sd_log2_outside_over_official=sd,
                             max_abs_log2=float(lr.abs().max()), glomerulus_max=str(lr.abs().idxmax()),
                             enters_with=choice[mb]))
    pd.DataFrame(roi_rows).to_csv(out / 'a1_roi.csv', index=False)
    chosen = pd.DataFrame({mb: counts5[(mb, choice[mb])].reindex(gl).fillna(0) for mb in ALL_MBS if mb in choice})
    chosen.index.name = 'glomerulus'
    rule = dict(rules[f"{choice['hemibrain']}_thr{SYN_MIN}"], roi=choice['hemibrain'], threshold=SYN_MIN)
    # FlyWire a second way (Pratyush's table)
    pf = pratyush_flywire()
    fr = []
    for mb, g in pf.dropna(subset=['glomerulus', 'side']).groupby('mb'):
        n_conn = g[g.weight >= SYN_MIN].groupby('glomerulus').size().reindex(gl).fillna(0)
        n_syn = g.groupby('glomerulus')['weight'].sum().reindex(gl).fillna(0)
        d = chosen[mb] if mb in chosen else pd.Series(0.0, index=gl)
        fr.append(pd.DataFrame(dict(mb=mb, glomerulus=gl, pratyush_connections_ge5=n_conn.to_numpy(),
                                    pratyush_synapses=n_syn.to_numpy(), derived_claws=d.reindex(gl).to_numpy())))
    fr = pd.concat(fr, ignore_index=True)
    summ = {}
    hb = shares(chosen['hemibrain'], gl)
    for mb, g in fr.groupby('mb'):
        g = g.set_index('glomerulus')
        sh_c, sh_s, sh_d = shares(g.pratyush_connections_ge5, gl), shares(g.pratyush_synapses, gl), shares(g.derived_claws, gl)
        m = (g.derived_claws >= SMALL_CLAWS) & (g.pratyush_connections_ge5 >= SMALL_CLAWS) & (chosen['hemibrain'] >= SMALL_CLAWS)
        lc, ld = np.log2(sh_c[m] / hb[m]), np.log2(sh_d[m] / hb[m])
        summ[mb] = dict(n_glomeruli=int(m.sum()), r_log2_ratio_conn_vs_derived=float(np.corrcoef(lc, ld)[0, 1]),
                        sd_log2_conn_over_derived=float(np.log2(sh_c[m] / sh_d[m]).std(ddof=1)),
                        sd_log2_syn_over_derived=float(np.log2(sh_s[m] / sh_d[m]).std(ddof=1)),
                        r_log2_share_conn_vs_derived=float(np.corrcoef(np.log2(sh_c[m]), np.log2(sh_d[m]))[0, 1]))
        idx = fr.mb == mb
        fr.loc[idx, 'share_conn'] = fr.loc[idx, 'glomerulus'].map(sh_c).to_numpy()
        fr.loc[idx, 'share_synapses'] = fr.loc[idx, 'glomerulus'].map(sh_s).to_numpy()
        fr.loc[idx, 'share_derived'] = fr.loc[idx, 'glomerulus'].map(sh_d).to_numpy()
    fr.to_csv(out / 'a1_flywire_pratyush.csv', index=False)
    return dict(rule=rule, rules=rules, roi_choice=choice, chosen=chosen, model=model, flywire_pratyush=summ)


def a2_ratios(out: Path, chosen: pd.DataFrame, model: pd.DataFrame) -> dict:
    """A2 (module docstring): ratios.csv, a2_spread.csv, a2_agreement.csv, a2_departures.csv."""
    gl = list(chosen.index)
    orn = H.apl16_orn_deltas()
    with_input = [g for g in gl if g in set(orn.index)]
    r, sm = ratios(chosen, gl)
    rm = mean_proportions(r)
    rt = pd.DataFrame(index=pd.Index(gl, name='glomerulus'))
    rt['has_input'] = [g in set(with_input) for g in gl]
    rt['model_claws'] = model['n_claws'].reindex(gl).to_numpy()
    rt['stratum'] = strata(rt['model_claws'], with_input)
    for mb in chosen.columns:
        rt[f'n_{mb}'] = chosen[mb]
        rt[f'share_{mb}'] = shares(chosen[mb], gl)
        rt[f'r_{mb}'] = r[mb]
        rt[f'small_{mb}'] = sm[mb]
        rt[f'fixed_{mb}'] = sm[mb]
    rt['r_mean'] = rm
    # the mean proportions are 1 by the small rule only where every body of the three brains has r_g = 1
    rt['fixed_mean'] = sm[[m for mbs in BRAINS.values() for m in mbs]].all(axis=1)
    rt['log2_r_mean'] = np.log2(rt['r_mean'])
    rt.reset_index().to_csv(out / 'ratios.csv', index=False)
    # spread against the counting noise
    sp = []
    for mb in list(OTHER_MBS) + list(BANC_MBS) + ['mean']:
        for subset, gs in (('all54', gl), ('input36', with_input)):
            col = 'r_mean' if mb == 'mean' else f'r_{mb}'
            small = rt['fixed_mean'] if mb == 'mean' else rt[f'small_{mb}']
            g_ok = [g for g in gs if not small[g]]
            lr = np.log2(rt.loc[g_ok, col].astype(float))
            row = dict(mb=mb, subset=subset, n_glomeruli=len(g_ok), n_small=len(gs) - len(g_ok),
                       sd_log2=float(lr.std(ddof=1)), q05=float(2 ** lr.quantile(0.05)), median=float(2 ** lr.median()),
                       q95=float(2 ** lr.quantile(0.95)), min=float(2 ** lr.min()), max=float(2 ** lr.max()))
            if mb != 'mean':
                cn = counting_noise_sd(rt.loc[g_ok, f'n_{mb}'], rt.loc[g_ok, 'n_hemibrain'])
                row.update(noise_sd_analytic=cn['analytic'], noise_sd_simulated=cn['simulated'],
                           excess_sd=float(np.sqrt(max(row['sd_log2'] ** 2 - cn['simulated'] ** 2, 0.0))))
            sp.append(row)
    spread = pd.DataFrame(sp)
    spread.to_csv(out / 'a2_spread.csv', index=False)
    # two sides of a brain against two brains
    ag = []
    for a, b in combinations([m for m in ALL_MBS if m in chosen.columns], 2):
        ok = (chosen[a] >= SMALL_CLAWS) & (chosen[b] >= SMALL_CLAWS)
        la, lb = np.log2(shares(chosen[a], gl)[ok]), np.log2(shares(chosen[b], gl)[ok])
        cn = counting_noise_sd(chosen.loc[ok, a], chosen.loc[ok, b])
        ag.append(dict(mb_a=a, mb_b=b, kind='same brain' if a.split('_')[0] == b.split('_')[0] else 'two brains',
                       banc=bool('banc' in a or 'banc' in b), n_glomeruli=int(ok.sum()),
                       r_log2_share=float(np.corrcoef(la, lb)[0, 1]), sd_log2_ratio=float((la - lb).std(ddof=1)),
                       noise_sd=cn['simulated']))
    agree = pd.DataFrame(ag)
    agree.to_csv(out / 'a2_agreement.csv', index=False)
    # where the hemibrain is furthest from the other brains
    lr4 = np.log2(rt[[f'r_{m}' for m in OTHER_MBS]].astype(float))
    dep = pd.DataFrame(index=rt.index)
    dep['has_input'] = rt['has_input']
    dep['model_claws'] = rt['model_claws']
    dep['n_hemibrain'] = rt['n_hemibrain']
    dep['mean_log2_r4'] = lr4.mean(axis=1)
    dep['sd_log2_r4'] = lr4.std(axis=1, ddof=1)
    n4 = rt[[f'n_{m}' for m in OTHER_MBS]].clip(lower=1)
    se = np.sqrt((1.0 / n4).mean(axis=1) / 4 + 1.0 / rt['n_hemibrain'].clip(lower=1)) / np.log(2)
    dep['counting_se'] = se
    dep['z'] = dep['mean_log2_r4'] / se
    dep['log2_r_mean'] = rt['log2_r_mean']
    dep['all_four_same_sign'] = (np.sign(lr4).nunique(axis=1) == 1) & (lr4 != 0).all(axis=1)
    sh5 = rt[[f'share_{m}' for m in FIVE_MBS]]
    dep['hemibrain_rank_of_5'] = sh5.rank(axis=1, method='min')['share_hemibrain']
    dep['small_any'] = rt[[f'small_{m}' for m in OTHER_MBS]].any(axis=1)
    dep = dep.sort_values('z', key=np.abs, ascending=False)
    dep.reset_index().to_csv(out / 'a2_departures.csv', index=False)
    lr_mean = np.log2(rt.loc[~rt.fixed_mean, 'r_mean'])
    nb = agree[~agree.banc]
    return dict(spread=spread, agreement=agree, departures=dep, ratios=rt,
                summary=dict(n_glomeruli=len(gl), n_with_input=len(with_input), n_fixed_mean=int(rt.fixed_mean.sum()),
                             mean_sd_log2=float(lr_mean.std(ddof=1)),
                             mean_q05=float(2 ** lr_mean.quantile(0.05)), mean_q95=float(2 ** lr_mean.quantile(0.95)),
                             within_brain_sd=float(nb[nb.kind == 'same brain'].sd_log2_ratio.mean()),
                             between_brain_sd=float(nb[nb.kind == 'two brains'].sd_log2_ratio.mean()),
                             within_brain_r=float(nb[nb.kind == 'same brain'].r_log2_share.mean()),
                             between_brain_r=float(nb[nb.kind == 'two brains'].r_log2_share.mean())))


def a3_searched(out: Path, rt: pd.DataFrame, n_null: int = A3_NULL, seed: int = A3_SEED) -> dict:
    """A3 (module docstring): per glomerulus the six searches' mean log2 factor over the ten odors (E86's
    partE_entries.csv, `mean_log2`) against log2 of the mean proportions' ratio; sign agreement and ρ over the glomeruli
    with input (the small ones, r_g = 1 by rule, left out) and over A3_NINE, each against `n_null` permutations of the
    glomerulus labels over the glomeruli with input (read on the nine for the nine). a3_searched.csv, a3_rule.csv."""
    e = pd.read_csv(H.results_dir('E86_reduced_input_candidate') / 'partE_entries.csv')
    nine = sorted(e.loc[e.meets_read.astype(bool), 'glomerulus'].unique())
    if tuple(nine) != A3_NINE:
        raise RuntimeError(f"E86's recurring changes are on {nine}, not A3_NINE {A3_NINE}")
    x = e.groupby('glomerulus')['mean_log2'].mean()
    t = pd.DataFrame({'searched_mean_log2': x})
    t['n_odors'] = int(e.odor.nunique())
    t['log2_r_mean'] = rt['log2_r_mean'].reindex(t.index)
    t['fixed_mean'] = rt['fixed_mean'].reindex(t.index).fillna(True).astype(bool)
    t['in_nine'] = t.index.isin(A3_NINE)
    t['same_sign'] = np.sign(t.searched_mean_log2) == np.sign(t.log2_r_mean)
    t.index.name = 'glomerulus'
    t.reset_index().to_csv(out / 'a3_searched.csv', index=False)
    use = t[~t.fixed_mean & t.log2_r_mean.notna()]
    xs, ys = use.searched_mean_log2.to_numpy(), use.log2_r_mean.to_numpy()
    nine_m = use.in_nine.to_numpy()
    rng = np.random.default_rng(seed)
    null = {k: [] for k in (('input', 'sign'), ('input', 'rho'), ('nine', 'sign'), ('nine', 'rho'))}
    for _ in range(n_null):
        yp = ys[rng.permutation(len(ys))]
        null[('input', 'sign')].append(sign_agreement(xs, yp))
        null[('input', 'rho')].append(_spearman(xs, yp))
        null[('nine', 'sign')].append(sign_agreement(xs[nine_m], yp[nine_m]))
        null[('nine', 'rho')].append(_spearman(xs[nine_m], yp[nine_m]))
    rows = []
    for subset, m in (('input', np.ones(len(xs), dtype=bool)), ('nine', nine_m)):
        sa, rho = sign_agreement(xs[m], ys[m]), _spearman(xs[m], ys[m])
        ns, nr = np.asarray(null[(subset, 'sign')]), np.asarray(null[(subset, 'rho')])
        rows.append(dict(subset=subset, n_glomeruli=int(m.sum()), sign_agreement=sa,
                         null_p05=float(np.nanquantile(ns, 0.05)), null_median=float(np.nanmedian(ns)),
                         null_p95=float(np.nanquantile(ns, 0.95)), rho=rho, rho_null_p05=float(np.nanquantile(nr, 0.05)),
                         rho_null_p95=float(np.nanquantile(nr, 0.95)), rule=a3_rule(sa, ns), n_null=n_null,
                         n_left_out_small=int(t.fixed_mean.sum())))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'a3_rule.csv', index=False)
    return dict(table=t, rule=df)


def a4_distances(out: Path, rt: pd.DataFrame) -> pd.DataFrame:
    """A4 (module docstring): with each mushroom body's proportions as weights (the model's real-wiring claws per
    glomerulus times r_g), the claw-weighted distances of the lab's ORN input over the 66 megamat pairs (E93's
    `weighted_distances`), by group (E92's groups), against the hemibrain's (the model's claws as they are);
    a4_distances.csv."""
    from experiments import E92_steps_like_for_like as E92
    from experiments import E93_missing_glomeruli_input as E93
    orn = H.apl16_orn_deltas()
    base = rt['model_claws'].astype(float)
    weights = {'hemibrain': base}
    for mb in list(OTHER_MBS) + list(BANC_MBS):
        weights[mb] = base * rt[f'r_{mb}']
    weights['mean'] = base * rt['r_mean']
    d0 = E93.weighted_distances(orn, weights['hemibrain'])
    m = E92.group_masks(d0.index)
    rows = []
    for name, w in weights.items():
        d = E93.weighted_distances(orn, w).reindex(d0.index)
        for g in E92.TABLE_GROUPS:
            if g not in m:
                continue
            v, v0 = float(d[m[g]].mean()), float(d0[m[g]].mean())
            rows.append(dict(proportions=name, group=g, n_pairs=int(m[g].sum()), mean_distance=v, hemibrain=v0,
                             difference=v - v0, rho_66_with_hemibrain=_spearman(d.to_numpy(), d0.to_numpy())))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'a4_distances.csv', index=False)
    return df


def part_a() -> dict:
    """Part A (module docstring): A1 to A4; a_summary.json; raises after writing if A1's rule stops the experiment."""
    out = out_dir()
    t0 = time.time()
    a1 = a1_counts(out)
    a2 = a2_ratios(out, a1['chosen'], a1['model'])
    a3 = a3_searched(out, a2['ratios'])
    a4 = a4_distances(out, a2['ratios'])
    rt = a2['ratios']
    summ = dict(commits=dict(H.code_commits()), host=_host(), seconds=round(time.time() - t0, 1),
                a1_rule=a1['rule'], a1_rules=a1['rules'], roi_choice=a1['roi_choice'],
                flywire_pratyush=a1['flywire_pratyush'], a2=a2['summary'],
                a3={r.subset: dict(sign_agreement=r.sign_agreement, null_p05=r.null_p05, null_p95=r.null_p95, rho=r.rho,
                                   rule=r.rule, n=r.n_glomeruli) for r in a3['rule'].itertuples()},
                a4={r.group: r.difference for r in a4[a4.proportions == 'mean'].itertuples()},
                small_mean=sorted(rt.index[rt.fixed_mean.astype(bool)]),
                stops=bool(a1['rule']['stops']))
    _write_json(out / 'a_summary.json', summ)
    print(json.dumps({k: summ[k] for k in ('a1_rule', 'roi_choice', 'a2', 'a3', 'a4', 'stops')}, indent=1, default=float),
          flush=True)
    if summ['stops']:
        raise RuntimeError(f"E97 part A: A1's rule stops the experiment (SD {a1['rule']['sd']:.3f} > {A1_RULE_SD}; "
                           'outcome (i))')
    return summ




# ---------------------------------------------------------------------------
# part B: the held readings

def _done(path: Path) -> set:
    return set(pd.read_csv(path).tag) if path.exists() else set()


def runs_b() -> pd.DataFrame:
    p = out_dir() / 'runs_B.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=['tag', 'variant', 'wiring'])


def _run_dirs(sub: str = 'runs') -> List[Path]:
    arc = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    return [H.outputs_dir(EXPERIMENT) / sub, arc / sub, arc / 'from_viki' / sub, arc / 'from_hal' / sub]


def save_run(tag: str, o: Dict[str, tuple], W: pd.DataFrame, moved: pd.DataFrame, targets: pd.DataFrame,
             sub: str = 'runs') -> Dict[str, str]:
    """A run's outputs (APL on and off: responses, spike counts), its wiring matrix, its moved claws and its targets per
    type and glomerulus, under outputs/E97_glomerulus_claw_shares/<sub>/. Returns the files' names."""
    d = local_dir(sub)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')
    W.to_parquet(d / f'{tag}_wPNKC.parquet')
    moved.to_csv(d / f'{tag}_moved.csv', index=False)
    targets.to_csv(d / f'{tag}_targets.csv', index=False)
    return dict(wiring_file=f'{tag}_wPNKC.parquet', moved_file=f'{tag}_moved.csv', targets_file=f'{tag}_targets.csv')


def load_run(tag: str, sub: str = 'runs') -> Optional[Dict[str, tuple]]:
    """A saved run ({'on', 'off'}): this host's outputs, else the model-output archive (with viki's or hal's lane)."""
    for d in _run_dirs(sub):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def wiring_file(tag: str, sub: str = 'runs') -> Optional[Path]:
    for d in _run_dirs(sub):
        p = d / f'{tag}_wPNKC.parquet'
        if p.exists():
            return p
    return None


def fuller_input() -> Tuple[pd.DataFrame, str, dict]:
    """E93's fuller input F (`E93.input_of('F')`), its cache tag, and the record of E93's saved file (results/E93_.../
    inputs/F.csv, its sha256), checked against the frame."""
    from experiments import E93_missing_glomeruli_input as E93
    orn = E93.input_of(B6_INPUT)
    rec = _read_json(H.results_dir(E93.EXPERIMENT) / 'inputs' / 'sha256.json').get(B6_INPUT, {})
    sha = _sha(orn)
    if rec and rec.get('frame_sha256') != sha:
        raise RuntimeError(f"E93's input {B6_INPUT} differs from its saved record ({rec.get('frame_sha256')} vs {sha})")
    return orn, E83.input_tag_for(orn, f'E97_{B6_INPUT}'), dict(rec, frame_sha256=sha)


def zero_input_share(W: pd.DataFrame, orn_glomeruli: Sequence[str]) -> float:
    """The share of a wiring's claws on glomeruli without ORN input (17 % in the model's wiring with the lab's input)."""
    t = claw_table(W)
    t = t[t.col >= 0]
    return float((~t.glomerulus.isin(set(orn_glomeruli))).mean())


def _provenance_add(part: str, tag: str, prov: dict) -> None:
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script='experiments/E97_glomerulus_claw_shares.py', points={})
    d.setdefault('points', {})[tag] = prov
    _write_json(path, d)


def run_point(variant: str, wiring: str, rt: Optional[pd.DataFrame] = None, cc: Optional[dict] = None,
              built: Optional[dict] = None, part: str = 'B') -> dict:
    """One point: the variant's wiring built (its sha256 checked against part 0's record and its invariants again),
    the held reading (`held97`, with E93's fuller input for B6), the run saved with its matrix, moved claws and targets
    (`save_run`), its row appended to runs_B.csv (absolute parameters, factors, realism, the moved share, the share of
    claws without input, the ratios, the commits) and its provenance to provenance_B.json. Returns the row."""
    import gc
    rt = ratios_table() if rt is None else rt
    built = built_records() if built is None else built
    cc = H.code_commits() if cc is None else cc
    tag = run_tag(variant, wiring)
    pv = parse_variant(variant)
    W0 = base_wiring(wiring)
    W1, mv, tg, wtag = built_wiring(variant, wiring, rt)
    sha = _sha(W1)
    rec = built.get(tag)
    if rec is None or rec.get('sha256') != sha:
        raise RuntimeError(f"{tag}: the built matrix ({sha[:12]}) is not part 0's record ({(rec or {}).get('sha256', '')[:12]}); "
                           'run --part 0 again')
    chk = check_built(W0, W1, mv, tg)
    if not chk['passes']:
        raise RuntimeError(f'{tag}: the built wiring fails its invariants: {chk}')
    if pv['input'] == 'lab':
        orn, itag, irec = None, None, dict(file='lab', frame_sha256=_sha(H.apl16_orn_deltas()))
        gl_in = list(H.apl16_orn_deltas().index)
    else:
        orn, itag, irec = fuller_input()
        gl_in = list(orn.index)
    res = held97(wiring, W1, wtag, orn, itag)
    o = {'on': res['on'], 'off': res['off']}
    files = save_run(tag, o, W1, mv, tg)
    r = variant_ratios(variant, rt)
    n_claws = int((W0.to_numpy().astype(int).sum(axis=1) > 0).sum())
    row = E83._run_row(part, variant, wiring, tag, res, o, variant=variant, arm=pv['arm'], proportions=pv['proportions'],
                       control=pv['control'] or '', perm_seed=pv['seed'], scale=pv['scale'], draw=pv['draw'],
                       draw_seed=json.dumps(draw_seed(variant, wiring)), input=pv['input'], input_tag=itag or '',
                       input_record=json.dumps(irec), wiring_tag=wtag, wiring_sha256=sha, n_moved=int(len(mv)),
                       moved_share=len(mv) / n_claws, zero_input_share=zero_input_share(W1, gl_in),
                       zero_input_share_base=zero_input_share(W0, gl_in), breadth_equal_fit=res.get('breadth_equal_fit'),
                       ratios=json.dumps({g: float(v) for g, v in r.items()}),
                       params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
                       params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}), **files)
    row.update({f'process_{k}': v for k, v in cc.items()})
    E83._append_row(out_dir() / 'runs_B.csv', row)
    _provenance_add('B', tag, dict(H.model_provenance(orn_deltas=H.apl16_orn_deltas() if orn is None else orn, _wPNKC=W1),
                                   wiring_sha256=sha, ratios_sha256=_sha(r.to_frame()), host=_host()))
    tgts = targets()
    p0 = _rms(E83.pattern_error(E83.pair_distances(E83.saved_fit(wiring)['saved'], 'cal'), tgts['z8']).to_numpy())
    p1 = _rms(E83.pattern_error(E83.pair_distances(o, 'cal'), tgts['z8']).to_numpy())
    print(f"E97 {part} {tag}: {len(mv)} claws moved ({row['moved_share']:.3f}); thr x{res['thr_factor']:.4f}, "
          f"APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} on, {res['n_off']}/{res['target_off']} off); "
          f"pattern part (readout, all recordings) {p0:.4f} -> {p1:.4f}; flag '{row.get('realism_flag') or ''}'; "
          f"{res['seconds']:.0f}s", flush=True)
    del res, o
    gc.collect()
    return row


def part_b(arms: Sequence[str] = ARMS, wirings: Sequence[str] = ALL_WIRINGS, dry_run: bool = False,
           max_points: Optional[int] = None) -> pd.DataFrame:
    """Part B's points of `arms` on `wirings` (each arm on its own wirings, ARMS_ON; `plan`), resumable (runs_B.csv);
    needs part 0 to have passed (its built records). `max_points` stops after that many new points (the smoke test)."""
    todo = [(v, w) for v, w in plan(arms, wirings) if run_tag(v, w) not in _done(out_dir() / 'runs_B.csv')]
    print(f'E97 B: {len(todo)} points to run: ' + ', '.join(run_tag(v, w) for v, w in todo), flush=True)
    if dry_run:
        return runs_b()
    if not identity_passed():
        raise RuntimeError('run --part 0 (after --part A) first: the identities must hold before any scored point')
    cc = _commits_clean()
    rt, built = ratios_table(), built_records()
    for i, (v, w) in enumerate(todo):
        if max_points is not None and i >= max_points:
            break
        run_point(v, w, rt, cc, built, part='B' if not v.startswith('C2') else 'C2')
    return runs_b()


# ---------------------------------------------------------------------------
# scoring: the runs, the distances, the targets

READOUTS = ('cal', 'raw', 'zcell', 'zcell_lin')
REF = 'ref'                                         # the wirings as they are (their saved fits)
E93F = 'E93_F'                                      # E93's held points with the fuller input on the unmoved wirings
R2 = 'R2_mbs40'                                     # the 40 model flies: the eight wirings by the five mushroom bodies
R2_PERM = 'R2_perm40'                               # their permuted control: the eight as they are and B3's 32
D2 = 'draw2'                                        # the mean proportions with the second draw (B4 on S, C2 on T)
UNITS = (('S', S_WIRINGS), ('T', T_WIRINGS), ('all8', ALL_WIRINGS))


def targets() -> Dict[str, pd.Series]:
    """E89's `targets89`: halves A to D, all recordings (z8, the main target), the six typical (z6), the noise-free
    effect, baseline-noise units."""
    from experiments import E89_apl_breadth_weighting as E89
    return E89.targets89()


def runs_table() -> pd.DataFrame:
    """Every run of the readings: the eight saved fits ('ref', with their realism), E93's fuller-input points on the
    unmoved wirings ('E93_F', where saved) and part B's / C2's (runs_B.csv), with 'flagged' and 'at_bound'."""
    from experiments import E93_missing_glomeruli_input as E93
    rows = []
    for w in ALL_WIRINGS:
        P = E83.saved_fit(w)
        rows.append(dict(tag=f'{REF}_{w}', variant=REF, wiring=w, source='saved', at_bound_thr=False, at_bound_apl=False,
                         thr_factor=1.0, apl_factor=1.0, fit_folder=P['folder'], fit_tag=P['tag'], **E83.realism(P['saved'])))
    try:
        e93 = pd.read_csv(H.results_dir(E93.EXPERIMENT) / 'runs_B.csv')
        e93 = e93[e93.input == B6_INPUT] if 'input' in e93 else e93.iloc[:0]
        for r in e93.drop_duplicates('wiring', keep='last').itertuples():
            rows.append(dict(tag=r.tag, variant=E93F, wiring=r.wiring, source='E93', at_bound_thr=r.at_bound_thr,
                             at_bound_apl=r.at_bound_apl, thr_factor=r.thr_factor, apl_factor=r.apl_factor,
                             realism_flag=r.realism_flag, rate_on=r.rate_on, rate_off=r.rate_off))
    except OSError:
        pass
    df = pd.DataFrame(rows)
    rb = runs_b()
    if len(rb):
        rb = rb.assign(source='B')
        df = pd.concat([df, rb], ignore_index=True, sort=False)
    df = df.drop_duplicates(['variant', 'wiring'], keep='last').reset_index(drop=True)
    df['flagged'] = df.realism_flag.fillna('').astype(str).str.strip() != ''
    df['at_bound'] = df.at_bound_thr.fillna(False).astype(bool) | df.at_bound_apl.fillna(False).astype(bool)
    return df


def run_outputs(variant: str, wiring: str) -> Dict[str, tuple]:
    """A run of the readings: 'ref' the wiring's saved fit, 'E93_F' E93's saved point, else this experiment's."""
    if variant == REF:
        return E83.saved_fit(wiring)['saved']
    if variant == E93F:
        from experiments import E93_missing_glomeruli_input as E93
        o = E93.load_run(E93.run_tag(B6_INPUT, wiring))
    else:
        o = load_run(run_tag(variant, wiring))
    if o is None:
        raise RuntimeError(f'no saved run {variant} on {wiring}')
    return o


def distances(o: Dict[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances: 'cal' (E83's calibrated readout, K set on the run's APL-on counts), 'raw', 'zcell' /
    'zcell_lin' (E73's per-KC reading under the saturating and the linear map, E95's `zcell_distances_map`)."""
    if readout in ('zcell', 'zcell_lin'):
        from experiments import E95_static_apl as E95
        return E95.zcell_distances_map(o, 'saturating' if readout == 'zcell' else 'linear')
    return E83.pair_distances(o, readout)


def _dist(variant: str, wiring: str, readout: str) -> Dict[str, pd.Series]:
    """`distances` of a run, cached in this process and in outputs/.../distances_cache.pkl (by the sha256 of the run's
    spike counts and the readout)."""
    import pickle
    mem = _CACHE.setdefault('dist_mem', {})
    key = (variant, wiring, readout)
    if key in mem:
        return mem[key]
    path = local_dir() / 'distances_cache.pkl'
    if 'dist_store' not in _CACHE:
        _CACHE['dist_store'] = pickle.loads(path.read_bytes()) if path.exists() else {}
    store = _CACHE['dist_store']
    o = run_outputs(variant, wiring)
    skey = (_sha(o['on'][1]), _sha(o['off'][1]), readout)
    if skey not in store:
        store[skey] = distances(o, readout)
        _CACHE['dist_dirty'] = True
    mem[key] = store[skey]
    return mem[key]


def _flush_dist_cache() -> None:
    import pickle
    if _CACHE.get('dist_dirty'):
        try:
            (local_dir() / 'distances_cache.pkl').write_bytes(pickle.dumps(_CACHE['dist_store']))
        except OSError:
            pass
        _CACHE['dist_dirty'] = False


def populations(runs: pd.DataFrame) -> List[Tuple[str, str, List[Tuple[str, str]]]]:
    """[(variant, unit, [(variant, wiring), ...])] of every reading: each run alone; per variant S, T and the eight where
    it has all their runs; R2's 40 model flies (the eight as they are and B1's four bodies) and their permuted control
    (the eight as they are and B3's four permuted bodies) on the eight and on S and T (20 each); the second draw of the
    mean proportions (B4 on S, C2 on T)."""
    have = set(zip(runs.variant, runs.wiring))
    out = []
    for v, g in runs.groupby('variant', sort=False):
        for w in g.wiring:
            out.append((v, w, [(v, w)]))
        for unit, ws in UNITS:
            if all((v, w) in have for w in ws):
                out.append((v, unit, [(v, w) for w in ws]))
    for name, vs in ((R2, [REF] + [f'B1_{m}' for m in OTHER_MBS]),
                     (R2_PERM, [REF] + [f'B3_perm_{m}_s{PERM_SEEDS[0]}' for m in OTHER_MBS])):
        for unit, ws in UNITS:
            mem = [(v, w) for v in vs for w in ws]
            if all(m in have for m in mem):
                out.append((name, unit, mem))
    d2 = {w: (f'B4_mean_d{DRAW_SEEDS[1]}' if w in S_WIRINGS else f'C2_mean_d{DRAW_SEEDS[1]}') for w in ALL_WIRINGS}
    for unit, ws in UNITS:
        mem = [(d2[w], w) for w in ws]
        if all(m in have for m in mem):
            out.append((D2, unit, mem))
    return out


def _effect(d: Dict[str, pd.Series]) -> pd.Series:
    return E83._cs(H.apl_effect(d['on'], d['off']))


def score_row(d: Dict[str, pd.Series]) -> dict:
    """E92's `score_row92` (E83's `effect_scores` against every target of `targets89`, the amplitude and shape parts, the
    spread ratio and Pearson r on the per-pair APL effect, checked to add up) plus, against all recordings, the mean
    and RMS pattern error over each of E92's groups (the six near pairs among them) and the mean distances."""
    from experiments import E92_gain as E92G
    from experiments import E92_steps_like_for_like as E92
    row = E92G.score_row92(d)
    me = _effect(d)
    de = E83._cs(targets()['z8']).reindex(me.index)
    err = (me - me.mean()) - (de - de.mean())
    m = E92.group_masks(me.index)
    for g in E92.TABLE_GROUPS:
        if g in m and m[g].any():
            row[f'z8_err_{g}_mean'] = float(err[m[g]].mean())
            row[f'z8_err_{g}_rms'] = _rms(err[m[g]].to_numpy())
            row[f'effect_{g}'] = float(me[m[g]].mean())
    row['mean_on'] = float(E83._cs(d['on']).mean())
    row['mean_off'] = float(E83._cs(d['off']).mean())
    return row


def boot_targets(n: int = BOOT_N, seed: int = BOOT_SEED) -> Tuple[pd.MultiIndex, np.ndarray]:
    """The per-pair APL effect of `n` resamples of the recordings (each genotype's recordings with replacement; Yang's z,
    all recordings), pairs in E83's canonical order, checked against the main target at the full sample."""
    key = ('boot', n, seed)
    if key in _CACHE:
        return _CACHE[key]
    mats = {}
    for g in ('TNT3c', 'TNT'):
        x = H.yang_recording_rdms_values(g, 'graded')
        x = x.copy()
        x.columns = pd.MultiIndex.from_tuples([tuple(sorted(map(str, c))) for c in x.columns], names=['odor_a', 'odor_b'])
        mats[g] = x.sort_index(axis=1)
    if not mats['TNT3c'].columns.equals(mats['TNT'].columns):
        raise RuntimeError('the two genotypes\' recordings have different pairs')
    pairs = mats['TNT3c'].columns
    full = mats['TNT3c'].to_numpy().mean(axis=0) - mats['TNT'].to_numpy().mean(axis=0)
    z8 = E83._cs(targets()['z8']).reindex(pairs).to_numpy()
    if not np.allclose(full, z8, atol=1e-10):
        raise RuntimeError('the resampling\'s full sample is not the main target')
    rng = np.random.default_rng(seed)
    C, T = mats['TNT3c'].to_numpy(), mats['TNT'].to_numpy()
    B = np.empty((n, len(pairs)))
    for b in range(n):
        B[b] = C[rng.integers(0, len(C), len(C))].mean(axis=0) - T[rng.integers(0, len(T), len(T))].mean(axis=0)
    _CACHE[key] = (pairs, B)
    return pairs, B


def boot_patterns(me: pd.Series, pairs: pd.MultiIndex, B: np.ndarray) -> np.ndarray:
    """The pattern part of a per-pair effect against each resampled target (rows of B)."""
    m = E83._cs(me).reindex(pairs).to_numpy()
    m = m - m.mean()
    Bc = B - B.mean(axis=1, keepdims=True)
    return np.sqrt(np.mean((m[None, :] - Bc) ** 2, axis=1))


def readings(readouts: Sequence[str] = READOUTS) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """readings.csv (per variant, unit and readout `score_row`, the runs outside the bands and at a bracket's bound, the
    change against 'ref' on the same unit and readout as 'd_<column>'), pairs.csv (each population's per-pair distances
    and effect under 'cal' and 'raw') and boot.csv (per population under the readout against all recordings: the
    difference of its pattern part from the reference's on the same unit, and a few other pairs, with the 90 % interval
    of a paired bootstrap over the recordings)."""
    out = out_dir()
    runs = runs_table()
    meta = {(r.variant, r.wiring): r for r in runs.itertuples()}
    rows, prs, eff = [], [], {}
    for v, unit, mem in populations(runs):
        for ro in readouts:
            ds = [_dist(mv, mw, ro) for mv, mw in mem]
            d = ds[0] if len(ds) == 1 else E83.population(ds)
            rows.append(dict(variant=v, unit=unit, kind='wiring' if len(mem) == 1 else 'population', n_runs=len(mem),
                             readout=ro, n_flagged=int(sum(meta[m].flagged for m in mem)),
                             n_at_bound=int(sum(meta[m].at_bound for m in mem)), **score_row(d)))
            if ro in ('cal', 'raw'):
                me = _effect(d)
                eff[(v, unit, ro)] = me
                if len(mem) > 1:
                    prs.append(pd.DataFrame(dict(variant=v, unit=unit, readout=ro, odor_a=me.index.get_level_values(0),
                                                 odor_b=me.index.get_level_values(1),
                                                 on=E83._cs(d['on']).reindex(me.index).to_numpy(),
                                                 off=E83._cs(d['off']).reindex(me.index).to_numpy(), effect=me.to_numpy())))
    _flush_dist_cache()
    rd = pd.DataFrame(rows)
    base = rd[rd.variant == REF].set_index(['unit', 'readout'])
    cols = [c for c in rd.columns if c.endswith(('_loss', '_pattern', '_rho', '_amplitude', '_shape', '_sd_ratio', '_pearson'))
            or c in ('mean_effect', 'mean_on', 'mean_off')]
    for c in cols:
        rd[f'd_{c}'] = [v - base[c].get((u, ro), np.nan) if (u, ro) in base.index else np.nan
                        for u, ro, v in zip(rd.unit, rd.readout, rd[c])]
    rd.to_csv(out / 'readings.csv', index=False)
    if prs:
        pd.concat(prs, ignore_index=True).to_csv(out / 'pairs.csv', index=False)
    # the paired bootstrap over the recordings
    pairs, B = boot_targets()
    bp = {k: boot_patterns(me, pairs, B) for k, me in eff.items() if k[2] == 'cal'}
    brows = []

    def add(a, b):
        if a in bp and b in bp:
            dlt = bp[a] - bp[b]
            brows.append(dict(variant=a[0], unit=a[1], ref_variant=b[0], ref_unit=b[1], readout='cal',
                              diff_median=float(np.median(dlt)), diff_lo=float(np.quantile(dlt, 0.05)),
                              diff_hi=float(np.quantile(dlt, 0.95)), share_below_zero=float((dlt < 0).mean()), n_boot=len(dlt)))
    for (v, unit, ro) in bp:
        if v != REF and unit in ('S', 'T', 'all8'):
            add((v, unit, ro), (REF, unit, 'cal'))
    for unit in ('S', 'T', 'all8'):
        add((f'B6_mean_{B6_INPUT}', unit, 'cal'), (E93F, unit, 'cal'))
        add((f'B6_mean_{B6_INPUT}', unit, 'cal'), ('B2_mean', unit, 'cal'))
        add((R2_PERM, unit, 'cal'), (R2, unit, 'cal'))
        add((D2, unit, 'cal'), ('B2_mean', unit, 'cal'))
    add((f'B4_mean_d{DRAW_SEEDS[1]}', 'S', 'cal'), ('B2_mean', 'S', 'cal'))
    bt = pd.DataFrame(brows)
    bt.to_csv(out / 'boot.csv', index=False)
    return rd, pd.concat(prs, ignore_index=True) if prs else pd.DataFrame(), bt


def brain_term(readout: str = 'cal') -> pd.DataFrame:
    """brain_term.csv: on each wiring, the root mean square difference of the per-pair APL effect between two mushroom
    bodies' proportions (the hemibrain's the wiring as it is; B1's four), raw and with both means removed; beside it the
    same between the move's two draws at the mean proportions (B2 and B4 on S) and between the eight wirings as they are
    (the draw term computed the same way); the summary rows 'brain_term', 'move_draw_term', 'wiring_draw_term' are the
    root mean square over their pairs."""
    runs = runs_table()
    have = set(zip(runs.variant, runs.wiring))
    rows = []

    def diff(a, b):
        ea, eb = _effect(_dist(*a, readout)), _effect(_dist(*b, readout))
        x = (ea - eb.reindex(ea.index)).to_numpy()
        return _rms(x), _rms(x - x.mean())
    mbs = ['hemibrain'] + list(OTHER_MBS)
    vof = {'hemibrain': REF, **{m: f'B1_{m}' for m in OTHER_MBS}}
    for w in ALL_WIRINGS:
        for a, b in combinations(mbs, 2):
            if (vof[a], w) in have and (vof[b], w) in have:
                raw, cen = diff((vof[a], w), (vof[b], w))
                rows.append(dict(kind='brain', wiring=w, a=a, b=b, rms=raw, rms_centered=cen))
        d2 = f'B4_mean_d{DRAW_SEEDS[1]}' if w in S_WIRINGS else f'C2_mean_d{DRAW_SEEDS[1]}'
        if ('B2_mean', w) in have and (d2, w) in have:
            raw, cen = diff(('B2_mean', w), (d2, w))
            rows.append(dict(kind='move_draw', wiring=w, a='B2_mean', b=d2, rms=raw, rms_centered=cen))
    for a, b in combinations(ALL_WIRINGS, 2):
        raw, cen = diff((REF, a), (REF, b))
        rows.append(dict(kind='wiring_draw', wiring=f'{a}-{b}', a=a, b=b, rms=raw, rms_centered=cen))
    df = pd.DataFrame(rows)
    _flush_dist_cache()
    summ = []
    for k in ('brain', 'move_draw', 'wiring_draw'):
        g = df[df.kind == k]
        if len(g):
            summ.append(dict(kind=f'{k}_term', wiring='', a='', b='', rms=_rms(g.rms), rms_centered=_rms(g.rms_centered),
                             n=len(g)))
    df = pd.concat([df, pd.DataFrame(summ)], ignore_index=True)
    df['readout'] = readout
    df.to_csv(out_dir() / 'brain_term.csv', index=False)
    return df


def line_reading(pr: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """line.csv: each population's ester and alcohol groups' mean APL-on distances relative to the mean over the 120
    pairs (E92's `relative_distances`), placed against E92's part A line (`E92.line(readout)['S']`) by E92's
    `place_on_line`; on T and the eight the residual is taken relative to the reference's on the same unit (E93's
    `line_reading`)."""
    from experiments import E92_steps_like_for_like as E92
    pr = pd.read_csv(out_dir() / 'pairs.csv') if pr is None else pr
    pts = {}
    for (v, unit, ro), g in pr.groupby(['variant', 'unit', 'readout'], sort=False):
        on = pd.Series(g.on.to_numpy(), index=pd.MultiIndex.from_arrays([g.odor_a, g.odor_b]))
        pts[(v, unit, ro)] = E92.relative_distances(on)
    rows = []
    for (v, unit, ro), (est, alc) in pts.items():
        L = E92.line(ro)['S']
        ref = 0.0
        if unit != 'S' and (REF, unit, ro) in pts:
            e0, a0 = pts[(REF, unit, ro)]
            ref = float(e0 - (L['intercept'] + L['slope'] * a0))
        rows.append(dict(variant=v, unit=unit, readout=ro, ester_rel=est, alcohol_rel=alc, slope=L['slope'],
                         intercept=L['intercept'], ref_resid=ref, **E92.place_on_line(est, alc, L, ref)))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'line.csv', index=False)
    return df


CO_VARIANTS = (REF, 'B2_mean', f'B3_permmean_s{PERM_SEEDS[0]}', 'B3_invmean')
CO_CANDIDATE = 'gm'


def _co_job97(job: tuple) -> dict:
    """One run read by E88's readings 1 and 1b and by the geometric mean of E88's two estimates at L1 and L3 (E93's
    `_co_job93` on this experiment's runs)."""
    import os
    import warnings
    from experiments import E89_apl_breadth_weighting as E89
    variant, wiring, K, cand, rel, gap = job
    warnings.simplefilter('ignore')
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(k, '1')
    o = run_outputs(variant, wiring)
    mr = E89.ModelReading89(o, label=f'{variant}_{wiring}', K=K, rel_target=rel, gap_target=gap)
    res = dict(variant=variant, wiring=wiring, readings={})
    for rdn in ('r1', 'r1b'):
        F = mr.mean_factors(rdn)
        res['readings'][rdn] = {c: dict(s=F[c].s, O=F[c].O, R=F[c].R) for c in E89.CONDS}
    for level in ('L1', 'L3'):
        v = {c: mr.level_values(level, c, None)[0] for c in E89.CONDS}
        F = {c: E89.cand_factors(v[c], cand) for c in E89.CONDS}
        res['readings'][f'cand_{level}'] = {c: dict(s=F[c].s, O=F[c].O, R=F[c].R) for c in E89.CONDS}
    return res


def coresponse97(n_workers: int = 4, variants: Sequence[str] = CO_VARIANTS) -> pd.DataFrame:
    """coresponse.csv: per variant (CO_VARIANTS) and population (S, T) the changes of L (co-response) and R
    (recruitment) under APL of the four nearest ester pairs, the ester and the alcohol groups (E88's `log_changes`) in
    readings 1 and 1b and with the geometric mean of E88's two estimates at L1 and L3, against the flies' (90 %
    intervals over resampled recordings), as E93's `coresponse93`. Jobs cached in outputs/.../coresponse_cache.pkl."""
    import multiprocessing
    import pickle
    from concurrent.futures import ProcessPoolExecutor
    from experiments import E88_responder_sets as E88
    from experiments import E89_apl_breadth_weighting as E89
    runs = runs_table()
    runs = runs[runs.variant.isin(list(variants))]
    rel = E88.flies_reliability('z')
    gap = {c: v['thresholded'] for c, v in E89.flies_gap_target('z').items()}
    cache = local_dir() / 'coresponse_cache.pkl'
    done = pickle.loads(cache.read_bytes()) if cache.exists() else {}
    units = {n: [(u, wu) for u, wu in UNITS[:2] if set(wu) <= set(g.wiring)] for n, g in runs.groupby('variant')}
    jobs = [(r.variant, r.wiring, E83.saved_fit(r.wiring)['K'], CO_CANDIDATE, rel, gap) for r in runs.itertuples()
            if (r.variant, r.wiring) not in done and any(r.wiring in wu for _, wu in units[r.variant])]
    if jobs:
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=multiprocessing.get_context('spawn')) as ex:
            for res, j in zip(ex.map(_co_job97, jobs), jobs):
                done[(j[0], j[1])] = res
                cache.write_bytes(pickle.dumps(done))
                print(f'E97 coresponse {j[0]} {j[1]}', flush=True)
    fl = {}
    fv = E89.flies_values({})
    for subset in E88.SUBSETS:
        tnt = E88.flies_subset(subset)
        rf = E88.flies_factors('z', 'r1', tnt)
        F = {c: rf[c].mean() for c in E89.CONDS}
        fl[('r1', subset)] = (E88.log_changes(F, ('L', 'R')), [E88.log_changes(b, ('L', 'R')) for b in E88.flies_boot(rf)])
        fl[('cand', subset)] = (E88.log_changes(E89.flies_cand(fv, CO_CANDIDATE, None, tnt), ('L', 'R')),
                                [E88.log_changes(b, ('L', 'R')) for b in E89.flies_cand_boot(fv, CO_CANDIDATE, None, tnt)])
    rows = []
    for name, us in units.items():
        for unit, wu in us:
            for rdn in ('r1', 'r1b', 'cand_L1', 'cand_L3'):
                Fs = [done[(name, w)]['readings'][rdn] for w in wu]
                Fp = {c: E88.Factors(*(E88._mean_rows(np.stack([f[c][k] for f in Fs])) for k in ('s', 'O', 'R')))
                      for c in E89.CONDS}
                ch = E88.log_changes(Fp, ('L', 'R'))
                for subset in E88.SUBSETS:
                    fch, fb = fl[('r1' if rdn.startswith('r1') else 'cand', subset)]
                    for f in ('L', 'R'):
                        for grp in ('nearest4', 'ester', 'alcohol'):
                            k = f'chg_{f}_{grp}'
                            lo, hi = np.nanpercentile([b[k] for b in fb], E88.INTERVAL)
                            x = ch[k]
                            differs = bool(not (lo <= x <= hi) and abs(x - fch[k]) >= np.log(E89.B1_FACTOR))
                            rows.append(dict(variant=name, unit=unit, reading=rdn,
                                             candidate=CO_CANDIDATE if rdn.startswith('cand') else '', subset=subset,
                                             factor=f, group=grp, ratio=float(np.exp(x)), flies_ratio=float(np.exp(fch[k])),
                                             flies_lo=float(np.exp(lo)), flies_hi=float(np.exp(hi)), differs=differs))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'coresponse.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# part C: the verdict

def _val(rd: pd.DataFrame, variant: str, unit: str, col: str, readout: str = 'cal') -> float:
    r = rd[(rd.variant == variant) & (rd.unit == unit) & (rd.readout == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


def _gain(rd: pd.DataFrame, variant: str, unit: str = 'all8', target: str = 'z8', readout: str = 'cal') -> float:
    """The reference's pattern part minus the variant's on the same unit (positive: the variant is lower, better)."""
    return _val(rd, REF, unit, f'{target}_pattern', readout) - _val(rd, variant, unit, f'{target}_pattern', readout)


def _boot_hi(bt: pd.DataFrame, variant: str, unit: str, ref: str = REF) -> float:
    r = bt[(bt.variant == variant) & (bt.unit == unit) & (bt.ref_variant == ref)] if len(bt) else bt
    return float(r.diff_hi.iloc[0]) if len(r) else np.nan


def classify(gain8: float, gain_s: float, gain_t: float, controls_ok: bool, others_ok: bool, r2_gain: float,
             perm40_gain: float) -> dict:
    """Part C's categories from the gains (reference minus variant, on the pattern part; positive is lower): 'moves'
    (MOVE_8 or more on the eight, MOVE_ST on S and T, the controls and every other condition), 'a part' (PART_LO to
    MOVE_8 on the eight with the controls' conditions), 'worse' (WORSE or more higher), 'no move' (R1 and R2 within
    NO_MOVE), 'averaging' (R2 gains AVERAGING and its permuted control at least CONTROL_RATIO of that); 'outcome' the
    first that holds in that order, else 'none of the entry's'."""
    flags = dict(
        moves=bool(gain8 >= MOVE_8 and gain_s >= MOVE_ST and gain_t >= MOVE_ST and controls_ok and others_ok),
        a_part=bool(PART_LO <= gain8 < MOVE_8 and controls_ok),
        worse=bool(gain8 <= -WORSE),
        no_move=bool(abs(gain8) < NO_MOVE and (not np.isfinite(r2_gain) or abs(r2_gain) < NO_MOVE)),
        averaging=bool(np.isfinite(r2_gain) and r2_gain >= AVERAGING and np.isfinite(perm40_gain)
                       and perm40_gain >= CONTROL_RATIO * r2_gain))
    names = {'moves': 'moves', 'a_part': 'a part', 'worse': 'worse', 'no_move': 'no move', 'averaging': 'averaging'}
    outcome = next((names[k] for k in ('moves', 'a_part', 'worse', 'no_move', 'averaging') if flags[k]),
                   "none of the entry's")
    return dict(flags, outcome=outcome)


def verdict97(rd: pd.DataFrame, bt: pd.DataFrame, runs: Optional[pd.DataFrame] = None) -> dict:
    """Part C's verdict, fixed in the entry, on the pattern part under the calibrated readout against all recordings
    (R1: B2_mean on the eight, S and T; R2: the 40 model flies), with the controls (B3's permuted mean proportions on
    both seeds, the inverted ones), the six typical recordings, the bands, the raw-count loss, the mean APL effect and
    the bootstrap interval; outcome (f) (one body lowers the eight by DISAGREE or more, another raises it), (h) (B6
    against B2 on the eight); and, where the second draw exists (B4 on S, C2 on T), the same reading on it ('a part' and
    'moves' must hold there). verdict.json."""
    v = 'B2_mean'
    g8, gs, gt = _gain(rd, v, 'all8'), _gain(rd, v, 'S'), _gain(rd, v, 'T')
    perm = {s: _gain(rd, f'B3_permmean_s{s}', 'all8') for s in PERM_SEEDS}
    inv = _gain(rd, 'B3_invmean', 'all8')
    r2, r2p = _gain(rd, R2, 'all8'), _gain(rd, R2_PERM, 'all8')
    six = _gain(rd, v, 'all8', 'z6')
    n_flag = _val(rd, v, 'all8', 'n_flagged')
    raw = _val(rd, v, 'all8', 'z8_loss', 'raw') - _val(rd, REF, 'all8', 'z8_loss', 'raw')
    mean = _val(rd, v, 'all8', 'mean_effect')
    hi = _boot_hi(bt, v, 'all8')

    def conds(gain8):
        c = dict(controls_perm=bool(all(np.isfinite(p) and p < CONTROL_RATIO * gain8 for p in perm.values())),
                 controls_inv=bool(np.isfinite(inv) and inv < INVERTED_MAX),
                 r2_not_above=bool(np.isfinite(r2) and r2 >= 0),
                 six_keeps=bool(np.isfinite(six) and six >= SIX_KEEP * gain8),
                 bands=bool(n_flag == 0), raw_loss=bool(raw <= RAW_SLACK), mean_effect=bool(abs(mean - FLIES_MEAN) <= MEAN_SLACK),
                 interval_below_zero=bool(np.isfinite(hi) and hi < 0))
        return c, bool(c['controls_perm'] and c['controls_inv']), bool(all(c[k] for k in (
            'r2_not_above', 'six_keeps', 'bands', 'raw_loss', 'mean_effect', 'interval_below_zero')))
    c1, ctl, oth = conds(g8)
    cls = classify(g8, gs, gt, ctl, oth, r2, r2p)
    res = dict(variant=v, readout='cal', target='z8', gain8=g8, gainS=gs, gainT=gt,
               ref8=_val(rd, REF, 'all8', 'z8_pattern'), refS=_val(rd, REF, 'S', 'z8_pattern'),
               refT=_val(rd, REF, 'T', 'z8_pattern'), R1_8=_val(rd, v, 'all8', 'z8_pattern'),
               R1_S=_val(rd, v, 'S', 'z8_pattern'), R1_T=_val(rd, v, 'T', 'z8_pattern'),
               R2=_val(rd, R2, 'all8', 'z8_pattern'), R2_gain=r2, R2_perm=_val(rd, R2_PERM, 'all8', 'z8_pattern'),
               R2_perm_gain=r2p, perm_gain=perm, inv_gain=inv, six_gain=six, n_flagged=n_flag, raw_loss_change=raw,
               mean_effect=mean, boot_hi=hi, conditions=c1, **cls)
    # the parts of the pattern part beside it (the reading the verdict is drawn on)
    for k in ('amplitude', 'shape', 'sd_ratio', 'pearson', 'rho', 'loss'):
        res[f'R1_8_{k}'] = _val(rd, v, 'all8', f'z8_{k}')
        res[f'ref8_{k}'] = _val(rd, REF, 'all8', f'z8_{k}')
    # (f): the bodies alone on the eight
    mbg = {m: _gain(rd, f'B1_{m}', 'all8') for m in OTHER_MBS}
    res['mb_gain'] = mbg
    fin = [x for x in mbg.values() if np.isfinite(x)]
    res['disagree'] = bool(fin and max(fin) >= DISAGREE and min(fin) < 0)
    # (h): B6 against B2 on the eight
    b6 = _val(rd, f'B6_mean_{B6_INPUT}', 'all8', 'z8_pattern')
    res['B6_8'] = b6
    res['B6_minus_B2'] = b6 - res['R1_8'] if np.isfinite(b6) else np.nan
    res['B6_minus_E93F'] = b6 - _val(rd, E93F, 'all8', 'z8_pattern') if np.isfinite(b6) else np.nan
    res['B6_differs'] = bool(np.isfinite(res['B6_minus_B2']) and abs(res['B6_minus_B2']) >= B6_DIFF)
    # B5: graded with the size of the change
    res['B5'] = {f'x{s:g}': _gain(rd, f'B5_mean_x{s:g}', 'S') for s in B5_SCALES}
    res['B5']['x1'] = gs
    # the second draw
    if np.isfinite(_val(rd, D2, 'all8', 'z8_pattern')):
        d8, ds_, dt = _gain(rd, D2, 'all8'), _gain(rd, D2, 'S'), _gain(rd, D2, 'T')
        c2, ctl2, oth2 = conds(d8)
        cls2 = classify(d8, ds_, dt, ctl2, oth2, r2, r2p)
        res['draw2'] = dict(gain8=d8, gainS=ds_, gainT=dt, conditions=c2, **cls2)
        res['holds_on_draw2'] = bool(cls2['outcome'] == res['outcome'] or (res['outcome'] == 'a part' and cls2['moves']))
    else:
        res['draw2'] = None
        res['holds_on_draw2'] = None
    res['needs_draw2'] = bool(res['outcome'] in ('a part', 'moves'))
    res['final'] = (res['outcome'] if not res['needs_draw2'] else
                    (res['outcome'] if res['holds_on_draw2'] else ('pending: C2' if res['holds_on_draw2'] is None
                                                                   else f"{res['outcome']} on the first draw only")))
    _write_json(out_dir() / 'verdict.json', res)
    return res


def score(coresponse: bool = True, zcell: bool = True, n_workers: int = 4) -> dict:
    """The readings (`readings`, `brain_term`, `line_reading`, `coresponse97`), part C's verdict and part D's plan,
    the archive tables and comparison_table.md."""
    ro = READOUTS if zcell else ('cal', 'raw')
    rd, pr, bt = readings(ro)
    brain_term()
    if len(pr):
        line_reading(pr)
    if coresponse:
        coresponse97(n_workers)
    v = verdict97(rd, bt)
    print(json.dumps({k: v[k] for k in ('outcome', 'final', 'gain8', 'gainS', 'gainT', 'R1_8', 'ref8', 'R2', 'perm_gain',
                                        'inv_gain', 'mb_gain', 'disagree', 'B6_minus_B2')}, indent=1, default=float),
          flush=True)
    write_variants()
    comparison_table()
    return v


def part_c2(dry_run: bool = False, force: bool = False) -> pd.DataFrame:
    """The second draw on T (C2_mean_d2), only on a verdict of "a part" or "moves" (verdict.json), then `score` is run
    again."""
    v = _read_json(out_dir() / 'verdict.json')
    ok = v.get('outcome') in ('a part', 'moves')
    print(f"E97 C2: verdict {v.get('outcome')!r}: {'run' if ok else 'not run'}", flush=True)
    if not ok and not force:
        return runs_b()
    return part_b(('C2',), T_WIRINGS, dry_run=dry_run)


# ---------------------------------------------------------------------------
# part D: the family fitted again at the mean proportions (only on "a part" or "moves")

D_VARIANT = 'B2_mean'


def d_member(label: str) -> dict:
    """A wiring's fit of C2's family (E91's `d_member`)."""
    from experiments import E91_summed_drive_inhibition as E91
    return E91.d_member(label)


def d_wiring_tag(label: str) -> str:
    return f'e97mean-{label}'


def d_tag(label: str) -> str:
    """The tag E10_joint_fit.py gives a part D fit (its prefix rules: wiring-<tag>_, rate0.068_, aplboth_, then
    k2_cap3)."""
    return f"wiring-{d_wiring_tag(label)}_rate{H.data_target('ctrl_frac'):.3f}_aplboth_k2_cap3"


def d_init_params() -> Path:
    """D/init_params_mean.json in this experiment's results: per wiring the absolute parameters of B2's held point (the
    wiring's own fit with the somatic threshold and the tied APL scale re-held at the mean proportions), the start of each
    refit (an exact start at the fit's own parameters is refused by E10_joint_fit.py on another wiring's base)."""
    rb = runs_b()
    d = {r.wiring: {k: float(x) for k, x in json.loads(r.params_on).items()} for r in rb[rb.variant == D_VARIANT].itertuples()}
    path = out_dir() / 'D' / 'init_params_mean.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    _write_json(path, d)
    return path


def d_command(label: str, maxfev: int = D_MAXFEV) -> List[str]:
    """E10_joint_fit.py's command for wiring `label`'s fit at the mean proportions: C2's family (`-k 2 --aplpn-max 3
    --rate-target ctrl_frac --loss apl --apl-readout both --tag _cap3`), B2's matrix as --wiring-file (it holds the
    wiring's own re-drawing, so no --wiring), from B2's held point (`--init-params`, `d_init_params`)."""
    wf = wiring_file(run_tag(D_VARIANT, label))
    if wf is None:
        raise RuntimeError(f'no matrix of {D_VARIANT} on {label} (run part B first)')
    init = d_init_params()
    if label not in json.loads(init.read_text()):
        raise RuntimeError(f'no held point of {D_VARIANT} on {label} in runs_B.csv')
    return [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-max', '3',
            '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', '_cap3', '--maxfev',
            str(int(maxfev)), '--experiment', EXPERIMENT, '--init-params', f'{init}:{label}', '--wiring-file', str(wf),
            '--wiring-tag', d_wiring_tag(label)]


def _fits() -> pd.DataFrame:
    p = out_dir() / 'fits.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=['tag'])


def d_allowed(pop: str = 'S') -> Tuple[bool, str]:
    """The entry's rules for part D: S only on "a part" or "moves" that held on the second draw; T only if S's fitted
    population holds (`score_d`)."""
    v = _read_json(out_dir() / 'verdict.json')
    if not v:
        return False, 'no verdict (run --part score)'
    if pop == 'S':
        ok = v.get('outcome') in ('a part', 'moves') and bool(v.get('holds_on_draw2'))
        return ok, f"verdict {v.get('outcome')!r}, second draw {'holds' if v.get('holds_on_draw2') else 'not read / fails'}"
    sd = score_d(write=False)
    return bool(sd.get('S_holds')), f"S's fitted population {'holds' if sd.get('S_holds') else 'does not hold'}"


def part_d(pop: str = 'S', wirings: Optional[Sequence[str]] = None, dry_run: bool = False, force: bool = False
           ) -> List[str]:
    """Part D's fits, one process each (`d_command`; logs under outputs/.../D/), skipping those in fits.csv; refuses
    unless `d_allowed` (or `force`)."""
    import subprocess
    ok, why = d_allowed(pop)
    print(f'E97 D {pop}: {"allowed" if ok else "not allowed"} ({why})', flush=True)
    if not ok and not force:
        return []
    ws = [w for w in (S_WIRINGS if pop == 'S' else T_WIRINGS) if wirings is None or w in wirings]
    logs = local_dir('D')
    tags = []
    for w in ws:
        tag = d_tag(w)
        tags.append(tag)
        if tag in set(_fits().tag):
            print(f'E97 D {tag}: done', flush=True)
            continue
        cmd = d_command(w)
        print('E97 D:', ' '.join(cmd), flush=True)
        if dry_run:
            continue
        _commits_clean()
        with open(logs / f'{tag}.log', 'w') as fh:
            r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f'E97 D {tag}: E10_joint_fit.py exited {r.returncode} (log {logs / (tag + ".log")})')
    return tags


def score_d(write: bool = True) -> dict:
    """Part D's scores: per wiring and population (S, T, the eight) the pattern part under the readout against all
    recordings, refitted against the saved fits; S holds if its fitted population keeps at least D_HOLD of R1's gain on S
    against REF_S with every fit inside the bands; the eight by E91's rule (D_PASS_GAIN below REF_8, at least
    D_WIRINGS_LOWER lower, every fit inside the bands, the raw-count loss not above D_RAW_REF by more than RAW_SLACK, the
    mean APL effect within MEAN_SLACK of the flies'). d_scores.csv, d_scores.json."""
    from experiments import E72_wiring_population as E72
    fits = _fits()
    tg = targets()
    got, rows, res = {}, [], {}
    for w in ALL_WIRINGS:
        tag = d_tag(w)
        if tag in set(fits.tag):
            o = E72.load_fit(EXPERIMENT, tag)
            if o is not None:
                got[w] = o
    for w, o in got.items():
        sc = E83.effect_scores(E83.pair_distances(o, 'cal'), tg)
        sc0 = E83.effect_scores(E83.pair_distances(E83.saved_fit(w)['saved'], 'cal'), tg)
        rows.append(dict(unit=w, pattern=sc['z8_pattern'], pattern_ref=sc0['z8_pattern'], mean_effect=sc['mean_effect'],
                         raw_loss=E83.effect_scores(E83.pair_distances(o, 'raw'), tg)['z8_loss'],
                         lower=bool(sc['z8_pattern'] < sc0['z8_pattern']), **E83.realism(o)))
    by = {r['unit']: r for r in rows}
    for unit, ws in UNITS:
        if not all(w in got for w in ws):
            continue
        dp = {ro: E83.population([E83.pair_distances(got[w], ro) for w in ws]) for ro in ('cal', 'raw')}
        sc = E83.effect_scores(dp['cal'], tg)
        r = dict(unit=unit, pattern=sc['z8_pattern'], mean_effect=sc['mean_effect'],
                 raw_loss=E83.effect_scores(dp['raw'], tg)['z8_loss'], n_lower=int(sum(by[w]['lower'] for w in ws)),
                 n_flagged=int(sum(bool(str(by[w].get('realism_flag') or '').strip()) for w in ws)))
        rows.append(r)
        res[unit] = r
    v = _read_json(out_dir() / 'verdict.json')
    if 'S' in res and v:
        gain_r1 = float(v.get('gainS', np.nan))
        res['S_holds'] = bool(np.isfinite(gain_r1) and REF_S - res['S']['pattern'] >= D_HOLD * gain_r1
                              and res['S']['n_flagged'] == 0)
    if 'all8' in res:
        e = res['all8']
        res['eight_passes'] = bool(e['pattern'] <= REF_8 - D_PASS_GAIN and e['n_lower'] >= D_WIRINGS_LOWER
                                   and e['n_flagged'] == 0 and e['raw_loss'] <= D_RAW_REF + RAW_SLACK
                                   and abs(e['mean_effect'] - FLIES_MEAN) <= MEAN_SLACK)
    if write:
        pd.DataFrame(rows).to_csv(out_dir() / 'd_scores.csv', index=False)
        _write_json(out_dir() / 'd_scores.json', res)
    return res


def part_e() -> None:
    """Part E (only on a pass by D): the firing-ceiling family at the flies' 2 s pulse with both response rates re-held
    and S through the trial noise as E72's part E ran it, and the mean proportions in E95's refit of the family at the
    connectome's bouton share (S). Not written: the entry runs it only on a pass by D."""
    raise NotImplementedError("part E runs only on a pass by part D (d_scores.json's 'eight_passes'); not written")


# ---------------------------------------------------------------------------
# the archive tables and the comparison table

def write_variants() -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike
    counts, the harness's primary model RDM): one row per part B / C2 run with its absolute parameters, the variant's
    proportions, its wiring file and sha256, the factors and realism; part D's fits (fits.csv rows, whose distances
    E10_joint_fit.py wrote) are kept."""
    out = out_dir()
    rb = runs_b()
    rows, parts = [], []
    for _, r in rb.drop_duplicates('tag', keep='last').iterrows():
        pv = parse_variant(r.variant)
        what = {None: '', 'perm': ' permuted within strata (seed %s)' % pv['seed'], 'inv': ' inverted'}[pv['control']]
        rows.append(dict(tag=r.tag, variant=r.tag, part=r.part, arm=r.arm, base='bouton', k=E83.BREADTH_K,
                         breadth_k=E83.BREADTH_K,
                         kind=(f"C2's family on {r.wiring} with the claws moved between glomeruli to the "
                               f"{pv['proportions']} proportions{what}" + (f", log ratios x{pv['scale']:g}" if pv['scale'] not in (1.0, -1.0) else '')
                               + f" (draw {pv['draw']}); input {pv['input']}; breadth recomputed for the wiring; the somatic "
                               "threshold and the tied APL scale re-held at the fit's APL-off and APL-on counts"),
                         config=r.variant, wiring=r.wiring, ref_experiment=r.fit_folder, ref_fit=r.fit_tag,
                         wiring_file=f'runs/{r.wiring_file}', wiring_sha256=r.wiring_sha256, n_moved=r.n_moved,
                         moved_share=r.moved_share, fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale,
                         wKCAPL_scale=r.wKCAPL_scale, wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale,
                         wKCKC_scale=r.wKCKC_scale, kc_kc_inh=r.wKCKC_scale, thr_factor=r.thr_factor,
                         apl_factor=r.apl_factor, rate_on=r.rate_on, rate_off=r.rate_off, max_count_on=r.max_count_on,
                         max_count_off=r.max_count_off, realism_flag=r.realism_flag, subdir='runs'))
        o = load_run(r.tag)
        if o is None:
            continue
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    long_path = out / 'model_pair_distances_long.csv'
    keep = []
    if long_path.exists():
        old = pd.read_csv(long_path)
        keep = [old[old.tag.isin(set(_fits().tag))]]
    if keep or parts:
        pd.concat(keep + parts, ignore_index=True).to_csv(long_path, index=False)
    return df


def _f(x, spec: str = '.4f') -> str:
    try:
        return format(float(x), spec) if np.isfinite(float(x)) else '-'
    except (TypeError, ValueError):
        return '-'


VARIANT_NAMES = {REF: 'the eight wirings as they are (the hemibrain\'s proportions)',
                 'B2_mean': 'the mean proportions of the three brains', D2: 'the mean proportions, second draw',
                 R2: 'the 40 model flies (eight wirings by five mushroom bodies)', R2_PERM: 'their permuted control',
                 'B3_permmean_s1': 'the mean proportions permuted (seed 1)', 'B3_permmean_s2': 'the mean proportions permuted (seed 2)',
                 'B3_invmean': 'the mean proportions inverted', f'B6_mean_{B6_INPUT}': "the mean proportions with E93's fuller input",
                 E93F: "E93's fuller input on the wirings as they are"}


def comparison_table() -> str:
    """comparison_table.md: part A's key numbers and, where part B has run, the populations under the calibrated readout
    against all recordings (pattern part with its amplitude and shape parts, spread ratio, Pearson r and ρ on the per-pair
    APL effect, the mean APL effect, the raw-count loss, the response rates and the runs outside the bands) and the
    verdict."""
    out = out_dir()
    lines = ['# E97: the glomeruli\'s shares of the claws in other connectomes', '',
             'What stays the hemibrain\'s or the model\'s: only the glomeruli\'s claw totals per KC type are another '
             'brain\'s; which KC samples which glomerulus, the KCs, their types and each KC\'s number of claws are the '
             'model\'s wiring\'s; a moved claw keeps its own APL weights and takes the APL->bouton weight of its new bouton.', '']
    a = _read_json(out / 'a_summary.json')
    if a:
        r = a['a1_rule']
        lines += ['## Part A', '',
                  f"A1's rule: SD of log2(derived hemibrain share / model share) over {r['n_glomeruli']} glomeruli with "
                  f"{A1_RULE_MIN_CLAWS} or more claws = {_f(r['sd'], '.3f')} (stops above {A1_RULE_SD}): "
                  f"{'stops' if r['stops'] else 'passes'}. ROI versions: " +
                  ', '.join(f'{k} {v}' for k, v in a['roi_choice'].items()) + '.', '']
        s = a['a2']
        lines += [f"A2: the mean proportions' ratios, SD of log2 {_f(s['mean_sd_log2'], '.3f')}, 5th to 95th percentile "
                  f"{_f(s['mean_q05'], '.2f')} to {_f(s['mean_q95'], '.2f')}; two sides of a brain SD {_f(s['within_brain_sd'], '.3f')}, "
                  f"two brains {_f(s['between_brain_sd'], '.3f')}.", '']
        lines += ['| A3 subset | glomeruli | sign agreement | null 5th / 95th | ρ | rule |', '|---|---:|---:|---:|---:|---|']
        for k, x in a['a3'].items():
            lines.append(f"| {k} | {x['n']} | {_f(x['sign_agreement'], '.2f')} | {_f(x['null_p05'], '.2f')} / "
                         f"{_f(x['null_p95'], '.2f')} | {_f(x['rho'], '.2f')} | {x['rule']} |")
        lines.append('')
    rp = out / 'readings.csv'
    if rp.exists():
        rd = pd.read_csv(rp)
        rd = rd[(rd.readout == 'cal') & rd.unit.isin(['S', 'T', 'all8'])]
        lines += ['## Part B: populations, calibrated readout, all recordings', '',
                  '| variant | unit | pattern part | amplitude | shape | spread ratio | Pearson r | ρ | mean APL effect | runs outside the bands |',
                  '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        for r in rd.itertuples():
            lines.append(f"| {VARIANT_NAMES.get(r.variant, r.variant)} | {r.unit} | {_f(r.z8_pattern)} | {_f(r.z8_amplitude)} | "
                         f"{_f(r.z8_shape)} | {_f(r.z8_sd_ratio, '.2f')} | {_f(r.z8_pearson, '.2f')} | {_f(r.z8_rho, '.2f')} | "
                         f"{_f(r.mean_effect, '+.3f')} | {int(r.n_flagged)} |")
        lines.append('')
        lines.append('KC data: the floor of the pattern part 0.029, the flies\' mean APL effect +0.178.')
        lines.append('')
    v = _read_json(out / 'verdict.json')
    if v:
        lines += ['## Part C: the verdict', '', f"Outcome: {v['outcome']} (final: {v['final']}).", '',
                  f"R1 on the eight {_f(v['R1_8'])} against {_f(v['ref8'])}; on S {_f(v['R1_S'])} against {_f(v['refS'])}; "
                  f"on T {_f(v['R1_T'])} against {_f(v['refT'])}; R2 {_f(v['R2'])}.", '']
    text = '\n'.join(lines) + '\n'
    (out / 'comparison_table.md').write_text(text)
    return text


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the part's run record (commits, host, time, settings), points kept."""
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script='experiments/E97_glomerulus_claw_shares.py', points={})
    d.setdefault('runs', []).append(dict(commits=dict(H.code_commits()), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                                         **extra))
    _write_json(path, d)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True,
                    choices=['0', 'A', 'B', 'score', 'C2', 'D', 'Dscore', 'E', 'variants', 'table'])
    ap.add_argument('--arm', default='all', help="part B's arms, comma-separated (B1..B6) or 'all'")
    ap.add_argument('--wirings', default='all8', help="S, T, all8 or a comma-separated list of wiring labels")
    ap.add_argument('--pop', default='S', choices=['S', 'T'], help="part D's population")
    ap.add_argument('--max-points', type=int, default=None, help='part B: stop after this many new points')
    ap.add_argument('--dry-run', action='store_true', help='print the plan (and the commands) without running')
    ap.add_argument('--force', action='store_true', help='part C2 / D without their condition')
    ap.add_argument('--no-held', action='store_true', help='part 0 without the held reading of the real wiring')
    ap.add_argument('--no-co', action='store_true', help='score without the co-response readings')
    ap.add_argument('--no-zcell', action='store_true', help="score without E73's per-KC reading")
    ap.add_argument('--workers', type=int, default=4, help='score: processes for the co-response readings')
    a = ap.parse_args()
    ws = _wirings(a.wirings)
    arms = ARMS if a.arm == 'all' else tuple(x for x in a.arm.split(',') if x)
    if a.part == '0':
        provenance('0', wirings=list(ws))
        part_0(ws, held=not a.no_held)
    elif a.part == 'A':
        provenance('A')
        part_a()
    elif a.part == 'B':
        if not a.dry_run:
            provenance('B', arms=list(arms), wirings=list(ws), max_points=a.max_points)
        part_b(arms, ws, dry_run=a.dry_run, max_points=a.max_points)
    elif a.part == 'score':
        score(coresponse=not a.no_co, zcell=not a.no_zcell, n_workers=a.workers)
    elif a.part == 'C2':
        part_c2(dry_run=a.dry_run, force=a.force)
    elif a.part == 'D':
        part_d(a.pop, None if a.wirings == 'all8' else ws, dry_run=a.dry_run, force=a.force)
    elif a.part == 'Dscore':
        print(json.dumps(score_d(), indent=1, default=float), flush=True)
    elif a.part == 'E':
        part_e()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'table':
        print(comparison_table(), flush=True)


if __name__ == '__main__':
    main()
