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
  type names one glomerulus (no '+', not M); glomerulus names are mapped to the model's (hemibrain) names (NAME_MAP).
  A KC's mushroom body is the side of the APL it makes most of its synapses with (`kc_sides`).
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
  its APL-off breadth run cached under outputs/E97_glomerulus_claw_shares/), the wiring's own fit's parameters, the
  somatic threshold bisected to the fit's APL-off response count, then the tied APL scale to its APL-on count (E77's
  `hold_rates`, as E85's and E93's held reading).

Parts (CLI `--part`):
- 0, the identities, before any scored point (identity.json; parts B and D refuse to run until it passes):
  - unchanged_<wiring>: with every r_g = 1 the builder returns each of the eight wirings unchanged;
  - held_real: the held reading of the real wiring through this experiment's path (`build97` with the unchanged
    matrix) reproduces its saved run with 0 count differences;
  - built_<variant>_<wiring>: every wiring that part B or C builds has each KC's claw count and type the base's, the
    per-type totals unchanged, the glomerulus totals equal to the targets and every moved claw on a bouton of its new
    glomerulus (`check_built`; no simulation). Each wiring's claw table is cached under outputs/.../W/ for part A.
- A (no simulations): A1 the counts (a1_counts.csv: per mushroom body, ROI version, threshold, KC type class and
  glomerulus the claws, boutons and PNs; a1_hemibrain_vs_model.csv and A1's rule, a1_rule.json; a1_flywire_pratyush.csv,
  FlyWire's proportions a second way from olfsysm's FlyWire PN->KC table; a1_sides.csv); A2 the ratios (ratios.csv,
  the per-glomerulus ratios of the four other mushroom bodies, BANC's two and the mean proportions; a2_spread.csv
  against the counting noise; a2_agreement.csv, two sides of a brain against two brains; a2_departures.csv); A3 the
  searched inputs against the mean proportions (a3_searched.csv, a3_rule.json); A4 the claw-weighted input distances by
  group (a4_distances.csv); a_summary.json.
- B (`--arm` B1..B6, default all; `--wirings`, default all8; resumable, runs_B.csv; runs under
  outputs/E97_glomerulus_claw_shares/runs/ with the wiring matrix and the moved claws): held readings, 116 points:
  B1 each of the four other mushroom bodies' proportions on the eight (32); B2 the mean proportions on the eight (8);
  B3 the controls on the eight (the mean proportions permuted within classes, two seeds, 16; inverted, 8; each mushroom
  body's permuted, one seed, 32); B4 a second draw of the moved claws at the mean proportions on S (4); B5 the mean
  proportions at half and twice their log ratios on S (8); B6 the mean proportions with E93's fuller input on the eight
  (8).
- score (`--no-co`, `--no-zcell`): readings.csv per variant and unit (each wiring; S, T and the eight where all have a
  run; the 40 model flies of R2 and their permuted control) and readout ('cal', 'raw', 'zcell', 'zcell_lin') against
  every target of `targets89`, with the amplitude and shape parts, the spread ratio and Pearson r (spread.csv), the
  paired bootstrap over the recordings of each difference from its reference (boot.csv), the brain term (brain_term.csv),
  the groups against E92's line (line.csv), the six near pairs, the alcohol pairs' co-response change (coresponse.csv);
  then part C's verdict (verdict.json), part D's plan, comparison_table.md and the figure.
- C2 (`--part C2`): only on a verdict of "a part" or "moves", the second draw on T (4 points, with B4), then `score`.
- D (`--pop S|T`, `--wirings` for a lane, `--dry-run`, `--force`): only on "a part" or "moves", C2's family fitted again
  at the mean proportions (`E10_joint_fit.py -k 2 --wiring-file <B2's matrix> --wiring-tag <tag> --rate-target ctrl_frac
  --loss apl --apl-readout both`, D_MAXFEV evaluations, from each wiring's own fit); T only on S's hold.
- Dscore: S's hold (D_HOLD of R1's gain on S, every fit inside the bands), the eight against REF_8 by E91's rule.
- E: documented stubs (`part_e`): runs only on a pass by D.
- variants: variants.csv and model_pair_distances_long.csv for the archive (`archive_tools/archive_model_outputs.py
  E97_glomerulus_claw_shares`, then `archive_tools/migrate_model_output_params.py --apply`).
- table: comparison_table.md and the figure from what exists.

Outputs: results/E97_glomerulus_claw_shares/ (shared); spike counts, responses, wiring matrices, moved-claw lists and
targets under outputs/E97_glomerulus_claw_shares/ (runs/, W/, the breadth cache).

Usage (from the repository root; 16-thread lanes; check `uptime` first; long runs with -u):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E97_glomerulus_claw_shares --part 0
    ./venv/bin/python -u -m experiments.E97_glomerulus_claw_shares --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E97_glomerulus_claw_shares \
        --part B --wirings S > <log> 2>&1
    ./venv/bin/python -u -m experiments.E97_glomerulus_claw_shares --part score
"""
from __future__ import annotations

import argparse
import json
import socket
import time
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
# FlyWire's names of the hemibrain's (the model's) glomeruli, as the derived tables' README gives them
NAME_MAP = {'VC3': 'VC3l', 'VC5': 'VC3m', 'VM6': 'VC5'}
NAME_MAPPED = ('flywire', 'malecns', 'banc')
# the mushroom bodies: the hemibrain's (right side only), and both sides of the other three
MB_SIDES = {'hemibrain': ('R',), 'flywire': ('L', 'R'), 'malecns': ('L', 'R'), 'banc': ('L', 'R')}
FIVE_MBS = ('hemibrain', 'flywire_L', 'flywire_R', 'malecns_L', 'malecns_R')
OTHER_MBS = FIVE_MBS[1:]
BANC_MBS = ('banc_L', 'banc_R')
BRAINS = {'flywire': ('flywire_L', 'flywire_R'), 'malecns': ('malecns_L', 'malecns_R')}
# the ROI version each brain enters with (A1's rule may switch the male CNS to 'outside'; part A records the choice)
DEFAULT_ROI = {'hemibrain': 'official', 'flywire': 'official', 'malecns': 'official', 'banc': 'official'}
SYN_THRESHOLDS = (3, 5, 10)
SYN_MIN = 5                                         # a counted claw has this many PN synapses or more
SMALL_CLAWS = 30                                    # fewer counted claws in either mushroom body: r_g = 1
KC_CLASSES = ('g', 'ab', "a'b'")                    # the model's KC types; derived types map by prefix (`kc_class`)
FLYWIRE_PRATYUSH = 'from_pratyush/2024-09-13'       # olfsysm's FlyWire PN->KC table (A1's second way)

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

# ---------------------------------------------------------------------------
# part B (fixed before any run)

DRAW_SEEDS = (1, 2)                                 # the draws of the moved claws (B4 and C2 the second)
PERM_SEEDS = (1, 2)                                 # B3: the mean proportions permuted; each body's with the first
N_TERCILES = 3                                      # B3's permutation strata: ORN input or not x thirds of the claw total
B5_SCALES = (0.5, 2.0)                              # B5: the mean proportions' log ratios times these
B6_INPUT = 'F'                                      # E93's fuller input
ARMS_ON = {'B1': 'all8', 'B2': 'all8', 'B3': 'all8', 'B4': 'S', 'B5': 'S', 'B6': 'all8', 'C2': 'T'}

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

# ---------------------------------------------------------------------------
# part D (only on "a part" or "moves")

D_MAXFEV = 200
D_HOLD = 0.5                                        # S's population keeps at least half of R1's gain on S
D_PASS_GAIN, D_WIRINGS_LOWER = 0.003, 6             # the eight by E91's rule for a fitted family
D_RAW_REF = 0.073
D_TAG = '_e97mean'

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


# ---------------------------------------------------------------------------
# the derived tables: claws, boutons, sides, names

def glomerulus_name(pn_type: str, dataset: str) -> Optional[str]:
    """The model's (hemibrain) name of a derived table's PN type (e.g. 'VC3_adPN' in FlyWire -> 'VC3l'), or None when
    the type names more than one glomerulus ('+'), is multiglomerular (M) or untyped."""
    raise NotImplementedError


def kc_class(kc_type: str) -> str:
    """A derived KC type's class: 'g' (KCg*), 'ab' (KCab*), "a'b'" (KCa'b'*, KCapbp*), 'unknown' (KC alone, other)."""
    raise NotImplementedError


def kc_sides(dataset: str, roi: str = 'official') -> pd.DataFrame:
    """Per KC of a dataset's claw table: its side ('L' / 'R'), the majority side of the APLs it makes synapses with
    (`<dataset>_APL_MB_synapses.csv.gz`), with the KC's own side label (partner_side) and the side by the largest gap
    of the KCs' mean claw x beside, for A1's check; the hemibrain's KCs are all 'R'."""
    raise NotImplementedError


def derived_claws(dataset: str, roi: str = 'official') -> pd.DataFrame:
    """A dataset's derived claws with the columns kc_id, claw_id, kc_type, kc_class, n_PN_synapses, main_PN_id,
    main_bouton_id, glomerulus (the model's name or None), side, mb ('<dataset>' for the hemibrain, '<dataset>_<side>'
    otherwise)."""
    raise NotImplementedError


def derived_boutons(dataset: str, roi: str = 'official') -> pd.DataFrame:
    """A dataset's derived boutons with bouton_id, pn_id, n_tbars, glomerulus (the model's name or None), side, mb (the
    side of the KCs its claws are on, else by the bouton's x against the claws' sides)."""
    raise NotImplementedError


def mb_counts(claws: pd.DataFrame, boutons: Optional[pd.DataFrame] = None, syn_min: int = SYN_MIN,
              by_class: bool = False) -> pd.DataFrame:
    """Per mushroom body (and KC class, with `by_class`) and glomerulus: the claws with `syn_min` or more PN synapses,
    the boutons of those claws' main PNs that carry them, and the PNs."""
    raise NotImplementedError


def shares(counts: pd.Series, glomeruli: Sequence[str]) -> pd.Series:
    """A mushroom body's share per glomerulus: its counted claws of g over its counted claws on `glomeruli`."""
    raise NotImplementedError


def ratios(counts: pd.DataFrame, glomeruli: Sequence[str], ref: str = 'hemibrain', small: int = SMALL_CLAWS
           ) -> pd.DataFrame:
    """Each mushroom body's ratio r_g to `ref` (glomeruli x mushroom bodies, counts as glomeruli x mushroom bodies), with
    r_g = 1 where either has fewer than `small` counted claws; a column 'small' lists those."""
    raise NotImplementedError


def mean_proportions(r: pd.DataFrame, brains: Mapping[str, Sequence[str]] = BRAINS) -> pd.Series:
    """The geometric mean of r_g over the three brains (the hemibrain 1; each brain's sides averaged on log2 first)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the model's wirings, the builder and the controls

def claw_table(W: pd.DataFrame) -> pd.DataFrame:
    """Per row of a claw x bouton wiring: kc_id, claw_id, kc_type, the column (bouton) index, its glomerulus (None for
    a row without a bouton)."""
    raise NotImplementedError


def type_targets(counts: pd.DataFrame, r: Mapping[str, float]) -> pd.DataFrame:
    """Per KC type (rows) and glomerulus (columns), the target claw counts: each count times r_g, rescaled to the type's
    total and rounded by largest remainder (ties by glomerulus name), so that every type's total is unchanged."""
    raise NotImplementedError


def move_claws(W: pd.DataFrame, r: Mapping[str, float], seed: int) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The wiring `W` at proportions `r` (module docstring). Returns (the built matrix with W's index and columns, the
    moved claws (row, kc_id, claw_id, kc_type, from_glomerulus, from_bouton, to_glomerulus, to_bouton), the targets per
    type and glomerulus with the base's counts)."""
    raise NotImplementedError


def check_built(W0: pd.DataFrame, W1: pd.DataFrame, moved: pd.DataFrame, targets: pd.DataFrame) -> dict:
    """Part 0's invariants of a built wiring: each KC's claw count and type the base's, the per-type totals unchanged,
    the glomerulus totals per type equal to the targets, every moved claw on a bouton of its new glomerulus, every row
    not in `moved` unchanged; 'passes'."""
    raise NotImplementedError


def permuted_ratios(r: pd.Series, strata: pd.Series, seed: int, fixed: Sequence[str] = ()) -> pd.Series:
    """`r` permuted among the glomeruli within each stratum (one generator from `seed`); the glomeruli in `fixed` (the
    small ones, r_g = 1) keep theirs."""
    raise NotImplementedError


def strata(claw_totals: pd.Series, with_input: Sequence[str], n: int = N_TERCILES) -> pd.Series:
    """B3's strata per glomerulus: with ORN input or not, and the third of the claw total within that class."""
    raise NotImplementedError


def scaled_ratios(r: pd.Series, scale: float) -> pd.Series:
    """r ** scale (B5's half and double of the log ratios; -1 the inverted ratios of B3)."""
    raise NotImplementedError


def variant_ratios(variant: str) -> pd.Series:
    """The ratios of a variant (ratios.csv of part A): B1_<mb>, B2_mean, B3_permmean_s<k>, B3_invmean,
    B3_perm_<mb>_s1, B4_mean_d2, B5_mean_x<s>, B6_mean_F, C2_mean_d2."""
    raise NotImplementedError


def parse_variant(v: str) -> dict:
    """A variant's fields: arm, proportions ('mean' or a mushroom body), control ('perm', 'inv' or None), seed (of the
    permutation), scale, draw (the moved claws' draw: 1, or 2 for B4 / C2), input ('lab' or 'F')."""
    raise NotImplementedError


def arm_variants(arm: str) -> List[str]:
    """The variants of an arm (B1..B6, C2)."""
    raise NotImplementedError


def run_tag(variant: str, wiring: str) -> str:
    return f'{variant}_{wiring}'


def built_wiring(variant: str, wiring: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    """(matrix, moved claws, targets, wiring tag) of a variant on a wiring: `move_claws` of the wiring's own matrix at
    the variant's ratios and draw (seed from DRAW_SEEDS and the wiring's position in ALL_WIRINGS); the tag names the
    matrix by its sha256 (it names the cached APL-off breadth run)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the held reading

def build97(label: str, W: pd.DataFrame, wiring_tag: str, orn: Optional[pd.DataFrame] = None):
    """C2's family on wiring `label`'s fit with the claw x bouton matrix `W` (same index and columns as the wiring's) and
    its breadth recomputed (`E10_joint_fit.bouton_base(wiring_matrix=W, wiring_tag=...)`, cached under this
    experiment's outputs); the normalized APL weights must be the wiring's. Returns (simulate, info) as E83's
    `build_simulator`."""
    raise NotImplementedError


def held97(label: str, W: pd.DataFrame, wiring_tag: str, orn: Optional[pd.DataFrame] = None) -> dict:
    """The held reading of `build97`'s simulator (E77's `hold_rates` at the wiring's fit's targets and parameters)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# parts

def identity_passed() -> bool:
    raise NotImplementedError


def part_0(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """The identities (module docstring, part 0); identity.json (resumable per key); raises if one fails."""
    raise NotImplementedError


def a1_rule(model_counts: pd.Series, derived_counts: pd.Series, min_claws: int = A1_RULE_MIN_CLAWS,
            max_sd: float = A1_RULE_SD) -> dict:
    """A1's rule: the SD over the glomeruli with `min_claws` or more claws (the model's) of log2(derived share / model
    share); 'stops' if above `max_sd`."""
    raise NotImplementedError


def a3_rule(agreement: float, null: np.ndarray, q: Tuple[float, float] = A3_QUANTILES) -> str:
    """A3's rule (as E93's A2): 'with' above the null's upper quantile, 'against' below its lower, otherwise 'neither'."""
    raise NotImplementedError


def part_a() -> dict:
    """Part A (module docstring): A1 to A4; a_summary.json; raises after writing if A1's rule stops the experiment."""
    raise NotImplementedError


def part_b(arms: Sequence[str] = ('B1', 'B2', 'B3', 'B4', 'B5', 'B6'), wirings: Sequence[str] = ALL_WIRINGS,
           dry_run: bool = False) -> pd.DataFrame:
    """Part B's points of `arms` on `wirings` (each arm on its own wirings, ARMS_ON); resumable (runs_B.csv)."""
    raise NotImplementedError


def score(coresponse: bool = True, zcell: bool = True) -> dict:
    """The readings, part C's verdict and the tables (module docstring, score)."""
    raise NotImplementedError


def verdict97(rd: pd.DataFrame, runs: pd.DataFrame) -> dict:
    """Part C's verdict, fixed in the entry (no move, moves, a part, worse, averaging; outcome (f)'s disagreement)."""
    raise NotImplementedError


def part_c2(dry_run: bool = False) -> pd.DataFrame:
    """The second draw on T, only on a verdict of "a part" or "moves"."""
    raise NotImplementedError


def d_command(label: str, maxfev: int = D_MAXFEV) -> List[str]:
    """E10_joint_fit.py's command for wiring `label`'s fit at the mean proportions (B2's matrix, from the wiring's own
    fit)."""
    raise NotImplementedError


def part_d(pop: str = 'S', wirings: Optional[Sequence[str]] = None, dry_run: bool = False, force: bool = False
           ) -> List[str]:
    """Part D's fits on their condition."""
    raise NotImplementedError


def score_d() -> dict:
    raise NotImplementedError


def part_e() -> None:
    """Part E (only on a pass by D): the firing-ceiling family at the flies' 2 s pulse, S through the trial noise, and
    the mean proportions in E95's refit at the connectome's bouton share. Not written."""
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    raise NotImplementedError


def comparison_table() -> str:
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True,
                    choices=['0', 'A', 'B', 'score', 'C2', 'D', 'Dscore', 'E', 'variants', 'table'])
    ap.add_argument('--arm', default='all', help="part B's arms, comma-separated (B1..B6) or 'all'")
    ap.add_argument('--wirings', default='all8', help="S, T, all8 or a comma-separated list of wiring labels")
    ap.add_argument('--pop', default='S', choices=['S', 'T'], help="part D's population")
    ap.add_argument('--dry-run', action='store_true', help='print the plan (and the commands) without running')
    ap.add_argument('--force', action='store_true', help="part D without its condition")
    ap.add_argument('--no-co', action='store_true', help="score without the co-response readings")
    ap.add_argument('--no-zcell', action='store_true', help="score without E73's per-KC reading")
    a = ap.parse_args()
    raise NotImplementedError(a.part)


if __name__ == '__main__':
    main()
