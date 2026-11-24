"""The model read as each stage is recorded (E134; the main goal adopted on 2026-11-24, CLAUDE.md "The main goal for
new experiments", "Like for like"): a run's units sampled and pooled as the recordings were made, its values through one
monotone curve per stage set on pooled values only, and through the flies' measurement noise, one population per
(model wiring, stage, panel, condition).

STUB (phase 1, 2026-11-24): the interface only; bodies raise NotImplementedError except trivial helpers. Built on
experiments/E90_model_stages.py (the saved stage matrices, `read_stages`, `stage_units`, `claw_signal`, `subset_mask`,
`unit_meta`) and experiments/E119_stages_per_roi.py (the settled pre-odor level, `model_stage`; the pooling into ROIs,
`pool_units`; the noisy trials, `noisy_trials`), which phase 2 reuses rather than copies.

The interfaces it follows (the stubs of 28a0f73 and cbff04e, 2026-11-24 10:08 and 10:09; their names and types are
followed exactly):
- viki's data_prep/stage_population_data.py: the stage and panel names (`STAGES`, `PANELS`), the cohorts ('wild_type',
  'same_fly', 'control', 'tnt'), its `Population` (one stage's recorded population on one panel: per fly a units x
  odors table of trial means in baseline-noise units, `tables`; single trials, `trials`, columns (odor, repeat); the
  stage's `criterion`), `panel_odors`, `response_criterion`, `load_population`, `require_lock` and the lock's committed
  copy experiments/E134_lock.json. Here that class is `D.Population` (the flies' side).
- hal's kcq/population_fit.py: its `Population` (one fly's or one wiring's population: `values`, `responds`, `name`,
  `named_units`; `Population.from_criterion`), which every reading here returns (`P.Population`); its views, scores and
  the energy guard read them. A model condition 'on' / 'off' (as hal's E134_archive_score.py and E90 name them; E90's
  files say 'on' / 'acute_off') is read against the cohort `cohort_of(panel, condition)` gives.

What a reading does, stage by stage (defaults in `DEFAULT_READINGS`, the alternatives in `ALTERNATIVE_READINGS`; which
reading is the stage's own and which its second is the lock's to fix, `reading_definitions` is what the lock records):
- 'orn': the model's ORN input itself (spike deltas, glomeruli x odors) on the glomeruli the recorded flies have (the
  recordings are of axon terminals; viki's lock says so). Named units.
- 'pn_dendrite': the model's PN rates (window mean minus the settled pre-odor mean, E119's level) on the driver line's
  glomeruli (GH146's, those the recorded flies have; `driver_glomeruli`). Named units.
- 'pn_bouton': the model's boutons (after APL's inhibition at the bouton), those of the shared glomeruli (E90's and
  E119's subset), pooled into regions as the recordings' ROIs (E119's verdict pool: none), the glomerulus labels then
  dropped: a bouton region has no glomerulus label.
- 'kc_claw': the model's claws (the net drive, E90's `claw_inputs`, with a share of the KC's spike count where a
  reading says so, E90's `claw_signal`), pooled into claw regions of about 14 claws (`CLAW_REGION_SIZE`; E119's 'g14':
  claws of one glomerulus, as a claw ROI surrounds one PN bouton), not assigned to a KC.
- 'kc_soma': the spike counts per KC and odor, every KC of the run.
Then, for every stage but where a reading says 'none': one monotone curve (`fit_stage_curve`), set on the run's pooled
values (every unit and odor of the condition it is set on, as one distribution) against the flies' pooled values of the
same stage and panel (`pooled_fly_values`), never on anything about odors: the function sees two arrays of numbers and
no odor label. By default it is set on the APL-on (wild-type) condition and applied unchanged to APL off, so that the
change APL makes is not absorbed by the curve. Then the flies' measurement noise (`fly_noise`, `add_measurement_noise`):
independent noise on each of the flies' number of trials, the trial mean kept, at the flies' own trial-to-trial residual
in baseline-noise units. A unit responds where its value is above the flies' criterion for the stage
(`D.response_criterion`); the 'raw' soma reading keeps the model's own criterion, two or more spikes (the planner's probe,
hal's archive 'raw').

The panels a saved run is read on (`run_panels`): E90's runs are of the APL16 input (16 odors), so their PN, bouton and
claw stages are read on the panel's 12 megamat odors as panel 'megamat' (the recorded flies have all 17; the comparison
is on the odors both sides have), and their somata as 'apl16' (16 odors) and as 'megamat' (the 12). E107's runs are
somata only, on 'megamat' (17) and 'validation2' (22). E134's own step-3 runs (experiments/E134_population_baseline.py)
save every stage on their panel's input.

Seeds: the pooling's POOL_SEED + the wiring's index and the noise's NOISE_SEED + the wiring's index (E119's), the same
for both conditions; the index is `seed_index` (the family's wirings in E83's ALL_WIRINGS order).

Provenance: every population read carries `reading_provenance`'s dict (the run's files with their sha256, the settled
levels' file, the reading, the curve's knots, the noise, the pooling and its seeds, the flies' files and the lock's
sha256, the code's commits).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from kcq import population_fit as P

if TYPE_CHECKING:
    from data_prep import stage_population_data as D

# the names viki's data module uses (data_prep/stage_population_data.py; the test checks they are the same)
STAGES: Tuple[str, ...] = ('orn', 'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'mixtures', 'apl16')
NAMED_UNITS: Dict[str, bool] = {'orn': True, 'pn_dendrite': True, 'pn_bouton': False, 'kc_claw': False,
                                'kc_soma': False}
# the model's conditions (hal's and E90's names) and E90's file names for them
CONDITIONS: Tuple[str, ...] = ('on', 'off')
CONDITION_FILE: Dict[str, str] = {'on': 'on', 'off': 'acute_off'}
# E90's names of the stages it saved (<tag>_<cond>_<name>_<win|pre>.parquet); the somata are its spike counts
E90_STAGE: Dict[str, str] = {'pn_dendrite': 'pn', 'pn_bouton': 'bouton', 'kc_claw': 'claw'}

# the windows (s): E90's stage window for the 0.5 s pulse and its pre-odor window, and the settled pre-odor window (E98's
# SETTLED_PRE_WINDOW; E90_settled_preodor's levels replace E90's pre-odor means of its saved runs, as E119 reads them)
MODEL_WINDOW: Tuple[float, float] = (0.0, 0.5)
PRE_WINDOW: Tuple[float, float] = (-0.5, 0.0)
SETTLED_PRE_WINDOW: Tuple[float, float] = (-0.25, 0.0)
SETTLED_LEVELS = ('results', 'E90_model_stages', 'settled_preodor', 'preodor_levels.parquet')

# pooling into the recordings' ROIs, as labels (E119's): '1' none; 'g<n>' groups of n units of one glomerulus (a claw ROI
# surrounds one PN bouton, so its claws share their PN's input); 'b<n>' groups of n claws of one bouton (a
# microglomerulus; needs the wiring's claw x bouton matrix); 's<n>' groups of n neighbouring claws by position within a
# glomerulus (the claws' coordinates in E90's claw index); 'r<n>' random groups of n
CLAW_REGION_SIZE = 14
CLAW_POOLS: Tuple[str, ...] = ('g14', 'b14', 's14', 'r14', '1')
BOUTON_POOLS: Tuple[str, ...] = ('1', 'g3', 'r3')
# the units a stage keeps before pooling (E90's VERDICT_SUBSETS, E119's SUBSET): boutons of the 28 glomeruli the PN
# recordings and the model's input share ('shared'), or of the driver line's glomeruli ('driver'), or all; claws all
SUBSETS: Dict[str, Tuple[str, ...]] = {'pn_bouton': ('shared', 'driver', 'all'), 'kc_claw': ('all', 'shared')}
# E90's versions of what a claw's calcium reports (the spike signal's share of the amplitude)
CLAW_SPIKE_SHARES: Tuple[float, ...] = (0.0, 0.25, 0.5)
# which glomeruli a named stage is read on: those every recorded fly of the stage has, those any fly has, or the 28 the
# GH146 recordings, the pebbled recordings and the model's input share (E90's `shared_glomeruli`)
GLOMERULI_RULES: Tuple[str, ...] = ('every_fly', 'any_fly', 'shared28')

# the monotone curve per stage: 'quantile' maps the pooled values rank for rank onto the flies' pooled distribution
# (ties, the many zeros of the somata, map to the flies' mean over the quantile range they hold, so equal model values
# stay equal); 'affine' matches the pooled mean and SD; 'saturating' is pyolfsysm.readout's x / (x + K) form with its
# scale and K set on the pooled distribution; 'none' leaves the model's units (only with the model's own criterion)
CURVES: Tuple[str, ...] = ('quantile', 'affine', 'saturating', 'none')
# where a curve is set: on the APL-on (wild-type) condition, applied unchanged to APL off ('on'), or on each condition
# against its own cohort ('own'); per wiring ('wiring', each wiring a fly with its own indicator scale) or on the
# family's pooled values ('family')
CURVE_SET_ON: Tuple[str, ...] = ('on', 'own')
CURVE_SCOPES: Tuple[str, ...] = ('wiring', 'family')
N_CURVE_KNOTS = 201
# the flies' weight in the pooled distribution: each fly counts equally ('fly'), or each unit and odor ('entry')
FLY_WEIGHTING: Tuple[str, ...] = ('fly', 'entry')
# the measurement noise: 'residual' (the flies' trial-to-trial residual SD in baseline-noise units, pooled over flies,
# units and odors, with the n / (n - 1) correction); 'residual_binned' (the same, as a function of the trial mean, in
# bins of the flies' pooled values); 'reliability' (E90's and E119's: per odor, the SD that gives the flies' trial-to-
# trial reliability across units); 'none'
NOISE_MODELS: Tuple[str, ...] = ('residual', 'residual_binned', 'reliability', 'none')
N_NOISE_BINS = 10
# the response criterion: the flies' for the stage, applied after the curve ('flies'), or the model's own, two or more
# spikes, on raw counts ('model_spikes'; the KC somata only)
CRITERIA: Tuple[str, ...] = ('flies', 'model_spikes')
MODEL_SPIKES = 2.0
POOL_SEED = 0
NOISE_SEED = 0


@dataclass(frozen=True)
class Reading:
    """How one stage of a model is read, fixed before any model is read on it (the lock records it).

    Attributes:
        name: the reading's label in tables and file names.
        stage: one of STAGES.
        curve: one of CURVES.
        curve_set_on: one of CURVE_SET_ON.
        curve_scope: one of CURVE_SCOPES.
        fly_weighting: one of FLY_WEIGHTING (the flies' pooled distribution the curve is set against).
        noise_aware_curve: set the curve so that the model's values after the curve and the noise match the flies'
            pooled distribution (fixed-point iterations at a fixed seed), not the values before the noise.
        noise: one of NOISE_MODELS.
        criterion: one of CRITERIA.
        pool: the pooling label (CLAW_POOLS for 'kc_claw', BOUTON_POOLS for 'pn_bouton'; None elsewhere).
        subset: the units kept before pooling (SUBSETS; None where a stage has no subset).
        claw_spike_share: the KC spike signal's share of a claw's signal (CLAW_SPIKE_SHARES; 'kc_claw' only).
        glomeruli: the rule for a named stage's glomeruli (GLOMERULI_RULES; 'orn' and 'pn_dendrite' only).
        settled_pre: the pre-odor mean at the settled level (SETTLED_PRE_WINDOW) rather than E90's PRE_WINDOW.
    """
    name: str
    stage: str
    curve: str = 'quantile'
    curve_set_on: str = 'on'
    curve_scope: str = 'wiring'
    fly_weighting: str = 'fly'
    noise_aware_curve: bool = False
    noise: str = 'residual'
    criterion: str = 'flies'
    pool: Optional[str] = None
    subset: Optional[str] = None
    claw_spike_share: float = 0.0
    glomeruli: Optional[str] = None
    settled_pre: bool = True


# PROVISIONAL until viki's lock fixes them: each stage's reading and its alternatives (the first alternative is the
# stage's second reading, the one a "not readable" verdict compares with)
DEFAULT_READINGS: Dict[str, Reading] = {
    'orn': Reading('curve_noise', 'orn', glomeruli='every_fly'),
    'pn_dendrite': Reading('curve_noise', 'pn_dendrite', glomeruli='every_fly'),
    'pn_bouton': Reading('curve_noise', 'pn_bouton', pool='1', subset='shared'),
    'kc_claw': Reading('curve_noise', 'kc_claw', pool='g14', subset='all', claw_spike_share=0.0),
    'kc_soma': Reading('curve_noise', 'kc_soma'),
}
ALTERNATIVE_READINGS: Dict[str, Tuple[Reading, ...]] = {
    'orn': (Reading('affine_noise', 'orn', curve='affine', glomeruli='every_fly'),
            Reading('curve_nonoise', 'orn', noise='none', glomeruli='every_fly')),
    'pn_dendrite': (Reading('affine_noise', 'pn_dendrite', curve='affine', glomeruli='every_fly'),
                    Reading('curve_noise_shared28', 'pn_dendrite', glomeruli='shared28')),
    'pn_bouton': (Reading('curve_noise_g3', 'pn_bouton', pool='g3', subset='shared'),
                  Reading('affine_noise', 'pn_bouton', curve='affine', pool='1', subset='shared')),
    'kc_claw': (Reading('curve_noise_b14', 'kc_claw', pool='b14', subset='all'),
                Reading('curve_noise_share25', 'kc_claw', pool='g14', subset='all', claw_spike_share=0.25),
                Reading('affine_noise', 'kc_claw', curve='affine', pool='g14', subset='all')),
    'kc_soma': (Reading('raw', 'kc_soma', curve='none', noise='none', criterion='model_spikes'),
                Reading('affine_noise', 'kc_soma', curve='affine'),
                Reading('curve_noise_own', 'kc_soma', curve_set_on='own')),
}


@dataclass
class SavedRun:
    """One saved model run (both conditions) as this module reads it.

    Attributes:
        source: where it was saved: 'E90:B' (E90's part B, the family on the APL16 input), 'E90:D1' (the firing-ceiling
            fits), 'E90:B5' (the 2 s pulse), 'E107:B' (somata only), 'E134:step3' (every stage, on its panel's input).
        tag: the run's tag in its folder (E90's wiring label or D1 model; E107's / E134's '<panel>_<reading>_<model>').
        model: the model ('C2' for the family; a beside model's name, experiments/E134_population_baseline.py
            BESIDE_MODELS).
        wiring: the wiring label (E83's ALL_WIRINGS for the family; the model's own otherwise).
        input_panel: the panel of the run's ORN input ('apl16' for E90's runs).
        conditions: 'on' / 'off' -> the saved matrices, E90's layout: 'pn_win', 'pn_pre', 'bouton_win', 'bouton_pre',
            'claw_win', 'claw_pre' (units x odors; E134's runs also 'pn_pre_settled', 'bouton_pre_settled',
            'claw_pre_settled'), 'apl' (APL's window and pre-odor means per odor), 'spike_counts', 'responses' (KCs x
            odors). A somata-only run has 'spike_counts' and 'responses' only.
        meta: the units' metadata, E90's units.json: 'bouton_glomerulus', 'claw_glomerulus', 'kc_of_claw'; and
            'claw_bouton' (each claw's bouton, from the wiring's claw x bouton matrix) where known.
        orn: the run's ORN input (glomeruli x odors, odors abbreviated), for the 'orn' stage; None if not known.
        files: every file read: dicts of path, bytes, sha256.
        settled: whether the pre-odor means are at the settled level.
    """
    source: str
    tag: str
    model: str
    wiring: str
    input_panel: str
    conditions: Dict[str, Dict[str, pd.DataFrame]]
    meta: Dict[str, np.ndarray] = field(default_factory=dict)
    orn: Optional[pd.DataFrame] = None
    files: List[dict] = field(default_factory=list)
    settled: bool = True

    def stages(self) -> Tuple[str, ...]:
        """The stages this run can be read at (every stage for a run with its stage matrices; 'kc_soma' alone for a
        somata-only run; 'orn' only where `orn` is known)."""
        c = self.conditions.get('on', {})
        out = []
        if self.orn is not None:
            out.append('orn')
        out += [s for s, k in E90_STAGE.items() if f'{k}_win' in c]
        if 'spike_counts' in c:
            out.append('kc_soma')
        return tuple(s for s in STAGES if s in out)


@dataclass
class StageCurve:
    """One stage's monotone curve from the model's values to the flies' units.

    Attributes:
        kind: one of CURVES.
        x: the curve's knots in the model's values (increasing).
        y: its values at the knots (non-decreasing); linear between knots, flat beyond the ends.
        params: the parametric curves' parameters ('affine': scale, offset; 'saturating': scale, K).
        set_on: what it was set on: the condition, the scope, the number of model and fly values, the flies' weighting,
            the sha256 of both pooled arrays.
    """
    kind: str
    x: np.ndarray
    y: np.ndarray
    params: Dict[str, float] = field(default_factory=dict)
    set_on: Dict[str, object] = field(default_factory=dict)

    def __call__(self, values: np.ndarray) -> np.ndarray:
        """The curve applied to an array of model values (any shape)."""
        raise NotImplementedError


@dataclass
class NoiseModel:
    """The flies' measurement noise of one stage and panel, in their units (baseline-noise SDs).

    Attributes:
        kind: one of NOISE_MODELS.
        n_trials: trials per odor of the recordings (the model's trials are averaged over as many).
        sd: the single-trial noise SD ('residual').
        bin_edges, bin_sd: the single-trial noise SD as a function of the trial mean ('residual_binned').
        reliability: per odor, the flies' trial-to-trial reliability across units ('reliability').
        from_flies: the flies, files and counts it was estimated from.
    """
    kind: str
    n_trials: int
    sd: float = np.nan
    bin_edges: Optional[np.ndarray] = None
    bin_sd: Optional[np.ndarray] = None
    reliability: Optional[pd.Series] = None
    from_flies: Dict[str, object] = field(default_factory=dict)


@dataclass
class ModelReading:
    """One model wiring's populations, read: per (stage, panel, condition) hal's `P.Population` and its provenance.

    Attributes:
        model, wiring: the model and wiring read.
        populations: (stage, panel, condition) -> the population (values in the flies' units after the curve and the
            noise; `responds` at the reading's criterion).
        before_curve: the same keys -> the model's pooled units x odors values before the curve (for the curve's
            checks and the 'raw' comparisons).
        curves: (stage, panel) -> the curve used (shared by both conditions where it is set on APL on).
        provenance: the same keys as `populations` -> `reading_provenance`'s dict.
    """
    model: str
    wiring: str
    populations: Dict[Tuple[str, str, str], P.Population] = field(default_factory=dict)
    before_curve: Dict[Tuple[str, str, str], pd.DataFrame] = field(default_factory=dict)
    curves: Dict[Tuple[str, str], StageCurve] = field(default_factory=dict)
    provenance: Dict[Tuple[str, str, str], dict] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# small helpers

def cohort_of(panel: str, condition: str) -> Optional[str]:
    """The flies' cohort (viki's `POPULATIONS`) a model condition is read against: 'control' / 'tnt' on the APL16
    panel, 'wild_type' with APL on elsewhere, and None with APL off outside the APL16 panel (no flies with APL output
    blocked there; the populations are still read, for the APL-effect statistics)."""
    if panel not in PANELS or condition not in CONDITIONS:
        raise ValueError(f'{panel=}, {condition=}')
    if panel == 'apl16':
        return 'control' if condition == 'on' else 'tnt'
    return 'wild_type' if condition == 'on' else None


def readings_of(stage: str) -> Tuple[Reading, ...]:
    """The stage's default reading, then its alternatives."""
    return (DEFAULT_READINGS[stage],) + tuple(ALTERNATIVE_READINGS.get(stage, ()))


def reading_by_name(stage: str, name: str) -> Reading:
    """The stage's reading called `name` (the default or an alternative); KeyError if none."""
    for r in readings_of(stage):
        if r.name == name:
            return r
    raise KeyError(f'{stage}: no reading {name!r}')


def reading_definitions() -> dict:
    """Every stage's readings as plain JSON-able dicts, the form the lock records them in: {stage: {'default': {...},
    'alternatives': [{...}, ...]}}, with the module's constants that a reading depends on (windows, region size,
    pools, seeds)."""
    return dict(
        stages={s: dict(default=asdict(DEFAULT_READINGS[s]),
                        alternatives=[asdict(r) for r in ALTERNATIVE_READINGS.get(s, ())]) for s in STAGES},
        constants=dict(MODEL_WINDOW=list(MODEL_WINDOW), PRE_WINDOW=list(PRE_WINDOW),
                       SETTLED_PRE_WINDOW=list(SETTLED_PRE_WINDOW), CLAW_REGION_SIZE=CLAW_REGION_SIZE,
                       N_CURVE_KNOTS=N_CURVE_KNOTS, N_NOISE_BINS=N_NOISE_BINS, MODEL_SPIKES=MODEL_SPIKES,
                       POOL_SEED=POOL_SEED, NOISE_SEED=NOISE_SEED))


def pooled_values(tables: Sequence[pd.DataFrame], odors: Optional[Sequence[str]] = None) -> np.ndarray:
    """Every finite entry of the tables (units x odors; restricted to `odors` where given), as one flat array: the
    pooled distribution, with no odor label left."""
    out = []
    for t in tables:
        x = (t[list(odors)] if odors is not None else t).to_numpy(dtype=float).ravel()
        out.append(x[np.isfinite(x)])
    return np.concatenate(out) if out else np.zeros(0)


# ---------------------------------------------------------------------------
# the saved runs

def load_e90_run(tag: str, part: str = 'B', settled: bool = True) -> SavedRun:
    """A run E90 saved with its stage matrices (local outputs/E90_model_stages/<part>/, else the archive
    data/model_outputs/E90_model_stages/<part>/), both conditions (E90's `load_stage`), with the units' metadata
    (<tag>_units.json, or the D1 model's wiring's) and the lab's input as its ORN stage.

    Args:
        tag: the wiring (part 'B': E83's ALL_WIRINGS), the firing-ceiling fit (part 'D1': E90's D1_MODELS) or
            'real_pulse2' (part 'B5').
        part: E90's part folder.
        settled: replace the pre-odor means by the settled levels (E90_settled_preodor's preodor_levels.parquet, as
            E119's `model_stage` does, with its check that the levels' E90-window column equals the saved pre-odor
            means).
    """
    raise NotImplementedError


def load_e107_run(panel: str, reading: str, wiring: str) -> SavedRun:
    """E107's saved spike counts and responses of one wiring on 'megamat' or 'validation2', 'fitted' or 'held'
    (outputs/E107_validation2_panel/B/<panel>_<reading>_<wiring>_<on|acute_off>_<responses|spike_counts>.parquet,
    else the archive), as a somata-only run, with its input (E107's `panel_input`) as its ORN stage."""
    raise NotImplementedError


def load_step3_run(tag: str) -> SavedRun:
    """A run of E134's step 3 (experiments/E134_population_baseline.py `load_run`): every stage, the settled pre-odor
    means saved beside E90's, on its panel's input."""
    raise NotImplementedError


def run_panels(run: SavedRun, stage: str) -> List[Tuple[str, Tuple[str, ...]]]:
    """The (panel, odors) a run's stage is read on: E90's APL16 runs give the PN, bouton and claw stages and the ORN
    input on 'megamat' (its 12 odors, H.MEGAMAT12) and the somata on 'apl16' (16) and 'megamat' (12); every other
    run its input's panel with all of its odors (viki's `panel_odors`, in panel order)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the units, sampled and pooled as recorded

def stage_changes(run: SavedRun, stage: str, condition: str, odors: Sequence[str],
                  reading: Optional[Reading] = None) -> pd.DataFrame:
    """A run's units x odors before sampling and pooling: 'orn' the input deltas; 'pn_dendrite', 'pn_bouton',
    'kc_claw' the window mean minus the pre-odor mean (E90's `stage_units`; at the settled level where the reading says
    so), the claws as E90's `claw_signal` at the reading's spike share (scales from the run's APL-on condition);
    'kc_soma' the spike counts. Odors abbreviated, restricted to `odors` in that order; the units' own index kept."""
    raise NotImplementedError


def driver_glomeruli(stage: str = 'pn_dendrite', panel: str = 'megamat', rule: str = 'every_fly') -> List[str]:
    """The glomeruli a named stage is read on: those the recorded flies of the stage and panel have (viki's
    `load_population`), every fly or any fly (`rule`), or the 28 shared ones (E90's `shared_glomeruli`), each only if
    the model has it."""
    raise NotImplementedError


def sample_driver_glomeruli(values: pd.DataFrame, glomeruli: Sequence[str]) -> pd.DataFrame:
    """The named units of `values` (glomeruli x odors) that are in `glomeruli`, in that order; KeyError if the model
    lacks one (the caller intersects first, `driver_glomeruli`)."""
    missing = [g for g in glomeruli if g not in values.index]
    if missing:
        raise KeyError(f'glomeruli not in the model: {missing}')
    return values.loc[list(glomeruli)]


def pool_boutons(values: pd.DataFrame, bouton_glomerulus: np.ndarray, pool: str = '1', subset: str = 'shared',
                 seed: int = POOL_SEED, subset_glomeruli: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """The model's boutons as the recordings' bouton regions: those of `subset` ('shared': the 28 glomeruli of E90's
    `shared_glomeruli`; 'driver': `subset_glomeruli`; 'all'), pooled by `pool` (BOUTON_POOLS; E119's `pool_units`),
    silent regions dropped (E90's `_drop_silent`), indexed 0..n-1 with no glomerulus label."""
    raise NotImplementedError


def pool_claws(values: pd.DataFrame, meta: Mapping[str, np.ndarray], pool: str = 'g14', subset: str = 'all',
               seed: int = POOL_SEED) -> pd.DataFrame:
    """The model's claws as claw regions of about CLAW_REGION_SIZE claws (CLAW_POOLS): 'g<n>' within each glomerulus
    (`meta['claw_glomerulus']`; E119's `pool_units`, a remainder of at least n / 2 kept as a smaller region), 'b<n>'
    within each bouton (`meta['claw_bouton']`), 's<n>' neighbouring claws by position within a glomerulus (the claw
    index's claw_x, claw_y, claw_z), 'r<n>' at random, '1' none; silent regions dropped, indexed 0..n-1, not assigned
    to a KC."""
    raise NotImplementedError


def sample_units(run: SavedRun, stage: str, condition: str, odors: Sequence[str], reading: Reading,
                 seed_index: int = 0) -> pd.DataFrame:
    """`stage_changes` then the stage's sampling or pooling under `reading` (`sample_driver_glomeruli`, `pool_boutons`,
    `pool_claws`; the somata as they are), with seed POOL_SEED + `seed_index` (the same in both conditions, so the
    same claws form the same regions)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the curve and the noise

def pooled_fly_values(flies: 'D.Population', odors: Optional[Sequence[str]] = None,
                      weighting: str = 'fly') -> Tuple[np.ndarray, np.ndarray]:
    """The flies' pooled distribution of one stage and panel: every fly's trial means (viki's `Population.tables`) on
    `odors` (default the panel's), flattened, and their weights ('fly': each fly sums to 1 / n_flies; 'entry': 1 each).
    Returns (values, weights); no odor label is kept."""
    raise NotImplementedError


def fit_stage_curve(model_values: np.ndarray, fly_values: np.ndarray, fly_weights: Optional[np.ndarray] = None,
                    kind: str = 'quantile', n_knots: int = N_CURVE_KNOTS, noise: Optional[NoiseModel] = None,
                    seed: int = NOISE_SEED) -> StageCurve:
    """The stage's monotone curve from the model's pooled values to the flies' pooled distribution.

    It is set on the two pooled distributions only: the arguments are flat arrays, with no odor, unit or pair label, so
    nothing about odors can enter it.

    Args:
        model_values: the model's pooled values (every unit and odor of the condition the curve is set on; for scope
            'family', of every wiring).
        fly_values: the flies' pooled values (`pooled_fly_values`).
        fly_weights: their weights (default equal).
        kind: one of CURVES ('none' returns the identity).
        n_knots: the knots of a 'quantile' curve (at evenly spaced quantiles of the model's values; ties map to the
            flies' weighted mean over the quantile range they hold).
        noise: when given, the curve is set so that the model's values after the curve and this noise (averaged over
            its n_trials) match the flies' pooled distribution (fixed-point iterations at `seed`).
        seed: the noise draws of a noise-aware fit.
    """
    raise NotImplementedError


def apply_curve(values: pd.DataFrame, curve: StageCurve) -> pd.DataFrame:
    """The curve applied entry by entry to a units x odors table (index and columns kept)."""
    return pd.DataFrame(curve(values.to_numpy(dtype=float)), index=values.index, columns=values.columns)


def fly_noise(flies: 'D.Population', kind: str = 'residual', odors: Optional[Sequence[str]] = None) -> NoiseModel:
    """The flies' measurement noise of one stage and panel from their single trials (viki's `Population.trials`,
    loaded with `load_population(..., trials=True)`): the residual of each trial from its unit's trial mean for the
    odor, pooled over flies, units and odors with the n / (n - 1) correction ('residual'), binned by the trial mean
    ('residual_binned'), or E119's per-odor reliability ('reliability'); n_trials the flies' median number of trials
    per odor. Raises ValueError where the population has no trials (the reading then needs noise 'none' and says so)."""
    raise NotImplementedError


def add_measurement_noise(values: pd.DataFrame, noise: NoiseModel, rng: np.random.Generator,
                          return_trials: bool = False):
    """`values` (units x odors, in the flies' units) through the flies' measurement noise: `noise.n_trials` trials,
    each the value plus independent noise of the model's kind, and their mean (units x odors). With `return_trials`
    also the trials (trials x units x odors), for readings that remove the noise as the flies' are (E119's
    `noise_removed`). 'none' returns `values` unchanged."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# readings

def fly_side(stage: str, panel: str, condition: str, odors: Optional[Sequence[str]] = None,
             trials: bool = True) -> Optional['D.Population']:
    """The flies a model condition is read against (viki's `load_population(stage, panel, cohort_of(panel,
    condition), trials=trials)`), None where there are none (APL off outside the APL16 panel; a population viki's
    `missing_populations` lists)."""
    raise NotImplementedError


def read_stage(run: SavedRun, stage: str, panel: str, odors: Sequence[str], reading: Reading,
               flies: Mapping[str, Optional['D.Population']], seed_index: int = 0,
               curve: Optional[StageCurve] = None) -> Tuple[Dict[str, P.Population], StageCurve, Dict[str, dict]]:
    """One stage of one run read on one panel, both conditions.

    Args:
        run: the saved run.
        stage, panel, odors: what is read (`run_panels`).
        reading: the stage's reading.
        flies: condition -> the flies it is read against (`fly_side`; the curve is set against the cohort of
            `reading.curve_set_on`, the noise taken from the same flies).
        seed_index: the wiring's index (pooling and noise seeds).
        curve: a curve already set (scope 'family'); None sets the run's own.

    Returns:
        ({'on', 'off'}: hal's `P.Population`, its responders at the reading's criterion), the curve, ({'on', 'off'}:
        `reading_provenance`).
    """
    raise NotImplementedError


def read_run(run: SavedRun, stages: Optional[Sequence[str]] = None, readings: Optional[Mapping[str, Reading]] = None,
             seed_index: int = 0, lock: Optional[dict] = None) -> ModelReading:
    """Every stage of a run (default all it saved, `SavedRun.stages`) on every panel it is read on (`run_panels`),
    both conditions, under `readings` (default `DEFAULT_READINGS`); the flies from viki's module (`fly_side`)."""
    raise NotImplementedError


def read_family(runs: Sequence[SavedRun], stages: Optional[Sequence[str]] = None,
                readings: Optional[Mapping[str, Reading]] = None, lock: Optional[dict] = None) -> Dict[str, ModelReading]:
    """`read_run` for each wiring of a family ({wiring: ModelReading}, seed index the wiring's position in `runs`);
    a reading of scope 'family' sets one curve per stage and panel on the wirings' pooled values first."""
    raise NotImplementedError


class SomaReading:
    """A KC-soma reading for spike counts alone (hal's experiments/E134_archive_score.py `Reading`: counts -> one
    `P.Population`), so that every archived variant is read as the family's runs are.

    Called with one condition's counts (KCs x odors, odors abbreviated): `reading(counts)` reads the APL-on condition;
    `reading(counts, condition='off', curve_counts=counts_on)` the APL-off condition with the curve set on the APL-on
    counts (a reading with `curve_set_on` 'on' refuses APL off without them, ValueError).

    Args:
        name: a 'kc_soma' reading's name (`readings_of('kc_soma')`).
        panel: the panel of the counts ('apl16' for the archive).
        seed_index: the noise seed's index.
    """

    def __init__(self, name: str, panel: str = 'apl16', seed_index: int = 0):
        self.reading = reading_by_name('kc_soma', name)
        self.panel = panel
        self.seed_index = seed_index

    def __call__(self, counts: pd.DataFrame, condition: str = 'on',
                 curve_counts: Optional[pd.DataFrame] = None) -> P.Population:
        raise NotImplementedError


def soma_readings() -> Dict[str, Callable[..., P.Population]]:
    """{name: `SomaReading(name)`} for every KC-soma reading, the readings hal's archive scoring may name."""
    return {r.name: SomaReading(r.name) for r in readings_of('kc_soma')}


def reading_provenance(run: SavedRun, stage: str, panel: str, condition: str, reading: Reading,
                       curve: Optional[StageCurve], noise: Optional[NoiseModel],
                       flies: Optional['D.Population'], n_units_in: int, n_units_out: int, seed_index: int,
                       lock: Optional[dict] = None) -> dict:
    """What a population read records: the run (source, tag, model, wiring, input panel, files with sha256, whether
    settled and the levels file's sha256), the stage, panel, odors and condition, the reading (`asdict`), the units
    before and after sampling or pooling, the curve (kind, knots or parameters, what it was set on), the noise (kind,
    SD or bins, n_trials, seed), the flies (cohort, flies, criterion, files from viki's lock), the lock's sha256, and
    `H.code_commits()`."""
    raise NotImplementedError
