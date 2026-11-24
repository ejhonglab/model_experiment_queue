"""The model read as each stage is recorded (E134; the main goal adopted on 2026-11-24, CLAUDE.md "The main goal for
new experiments", "Like for like"): a run's units sampled and pooled as the recordings were made, its values through one
monotone curve per stage set on pooled values only, and through the flies' measurement noise, one population per
(model wiring, stage, panel, condition).

Built on experiments/E90_model_stages.py (the saved stage matrices and their layout, `stage_units`, `claw_signal`,
`claw_inputs`, `shared_glomeruli`, `_drop_silent`) and experiments/E119_stages_per_roi.py (the settled pre-odor level,
as its `model_stage` reads it; the pooling into ROIs, `pool_units`), which it calls rather than copies.

The interfaces it follows (viki's and hal's modules as of 20f6603, 2026-11-24: viki's data module 4cf2c9f, its lock
607da67):
- viki's data_prep/stage_population_data.py: the stage and panel names (`STAGES`, `PANELS`: the kiwi and control
  mixtures are two panels, 'kmix' and 'cmix'), the cohorts ('wild_type', 'same_fly', 'control', 'tnt'), its
  `Population` (one stage's recorded population on one panel: per fly a units x odors table of trial means in the
  stage's units, `tables`; single trials, `trials`, columns (odor, repeat), where `has_trials`; the stage's
  `criterion`), `panel_odors`, `load_population`, `breadth_bins`, `require_lock` and the lock's committed copy
  experiments/E134_lock.json, which records `reading_definitions()` as its 'readings'. Here that class is
  `D.Population` (the flies' side). Every reading of real flies here (`fly_side`) refuses before the lock is committed.
- hal's kcq/population_fit.py: its `Population` (one fly's or one wiring's population: `values`, `responds`, `name`,
  `named_units`, `graded`; `Population.from_criterion`), which every reading here returns (`P.Population`); its views,
  scores and the energy guard read them. A model condition 'on' / 'off' (as hal's E134_archive_score.py and E90 name
  them; E90's files say 'on' / 'acute_off') is read against the cohort `cohort_of(panel, condition)` gives.

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
  claws of one glomerulus, as a claw ROI surrounds one PN bouton; 'b14' the claws of one bouton), not assigned to a KC.
- 'kc_soma': the spike counts per KC and odor, every KC of the run.
Regions (boutons, claws) that are silent through every odor of the condition are dropped, as E119 drops them (the
recordings' ROIs are of responsive structure).

Then, for every stage but where a reading says 'none': one monotone curve (`fit_stage_curve`), set on the run's pooled
values (every unit and odor of the condition it is set on, as one distribution) against the flies' pooled values of the
same stage and panel (`pooled_fly_values`, on the odors the run is read on), never on anything about odors: the function
sees two arrays of numbers and no odor label. By default it is set on the APL-on (wild-type) condition and applied
unchanged to APL off, so that the change APL makes is not absorbed by the curve. A curve can be set noise-aware
(`Reading.noise_aware_curve`: the model's values after the curve and after the noise then match the flies' pooled
distribution, which already holds the flies' noise); no reading the lock records is. Then the flies' measurement noise
(`fly_noise`, `add_measurement_noise`): independent noise on each of the flies' number of trials, the trial mean kept,
at the flies' own trial-to-trial residual in their units. Where the flies have no single trials (the mixtures' KC
somata) no noise is added and the provenance says so (NO_TRIALS_RULE). A unit responds where its value is above the
flies' criterion for the stage (`D.Population.criterion`, the lock's); the 'raw' soma reading keeps the model's own
criterion, two or more spikes (the planner's probe, hal's archive 'raw').

The panels a saved run is read on (`run_panels`): E90's runs are of the APL16 input (16 odors), so their PN, bouton and
claw stages and the ORN input are read on the panel's 12 megamat odors as panel 'megamat' (the recorded flies have all
17; the comparison is on the odors both sides have), and their somata as 'apl16' (16 odors) and as 'megamat' (the 12).
E107's runs are somata only, on 'megamat' (17) and 'validation2' (22). E134's own step-3 runs
(experiments/E134_population_baseline.py, saved by `save_step3_run`) hold every stage their model has, on their panel's
input (INPUT_PANELS); a mixture run ('mixtures') is read on both mixture panels, each on its own stimuli.

Seeds: the pooling's POOL_SEED + the wiring's index and the noise's NOISE_SEED + the wiring's index (E119's), the same
for both conditions (a new generator per condition, so both get the same draws); the index is `seed_index` (the family's
wirings in E83's ALL_WIRINGS order).

Provenance: every population read carries `reading_provenance`'s dict (the run's files with their sha256, the settled
levels' file, the reading, the curve's knots, the noise, the pooling and its seeds, the flies' files and the lock's
sha256, the code's commits).
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from kcq import population_fit as P

if TYPE_CHECKING:
    from data_prep import stage_population_data as D

# the names viki's data module uses (data_prep/stage_population_data.py; the test checks they are the same): the
# recorded panels, the kiwi and control mixtures two of them
STAGES: Tuple[str, ...] = ('orn', 'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'kmix', 'cmix', 'apl16')
MIXTURE_PANELS: Tuple[str, ...] = ('kmix', 'cmix')
# the panels of the models' inputs: one mixture input (`natmix_orn_deltas`, both mixture panels' stimuli) read on both
INPUT_PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'mixtures', 'apl16')
NAMED_UNITS: Dict[str, bool] = {'orn': True, 'pn_dendrite': True, 'pn_bouton': False, 'kc_claw': False,
                                'kc_soma': False}
# the model's conditions (hal's and E90's names) and E90's file names for them
CONDITIONS: Tuple[str, ...] = ('on', 'off')
CONDITION_FILE: Dict[str, str] = {'on': 'on', 'off': 'acute_off'}
# E90's names of the stages it saved (<tag>_<cond>_<name>_<win|pre>.parquet); the somata are its spike counts
E90_STAGE: Dict[str, str] = {'pn_dendrite': 'pn', 'pn_bouton': 'bouton', 'kc_claw': 'claw'}
STAGE_FILES: Tuple[str, ...] = ('pn', 'bouton', 'claw')

# the windows (s): E90's stage window for the 0.5 s pulse and its pre-odor window, and the settled pre-odor window (E98's
# SETTLED_PRE_WINDOW; E90_settled_preodor's levels replace E90's pre-odor means of its saved runs, as E119 reads them)
MODEL_WINDOW: Tuple[float, float] = (0.0, 0.5)
PRE_WINDOW: Tuple[float, float] = (-0.5, 0.0)
SETTLED_PRE_WINDOW: Tuple[float, float] = (-0.25, 0.0)
SETTLED_LEVELS = ('results', 'E90_model_stages', 'settled_preodor', 'preodor_levels.parquet')
LEVEL_TOL = 1e-6

# where the saved runs are (the local outputs, then the model-output archive)
E90_EXP = 'E90_model_stages'
E107_EXP = 'E107_validation2_panel'
STEP3_EXP = 'E134_population_baseline'
STEP3_PART = 'step3'
# E90's D1 fits and their wirings (results/E90_model_stages/d1.json)
E90_D1_WIRING = {'fit1': 'real', 'P1_fc': 'P1', 'P2_fc': 'P2', 'M1_fc': 'M1'}

# pooling into the recordings' ROIs, as labels (E119's): '1' none; 'g<n>' groups of n units of one glomerulus (a claw ROI
# surrounds one PN bouton, so its claws share their PN's input); 'b<n>' groups of n claws of one bouton (a
# microglomerulus; needs each claw's bouton, `claw_boutons`); 'r<n>' random groups of n. Within a group label a remainder
# of at least n / 2 units is kept as a smaller region (E119's `pool_units`)
CLAW_REGION_SIZE = 14
CLAW_POOLS: Tuple[str, ...] = ('g14', 'b14', 'r14', '1')
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
# stay equal; above the largest model value it continues with the slope of its top tenth); 'affine' matches the pooled
# mean and SD; 'saturating' is scale x / (|x| + K), pyolfsysm.readout's form, its scale and K set by least squares on
# the two pooled distributions' quantiles; 'none' leaves the model's units (only with the model's own criterion)
CURVES: Tuple[str, ...] = ('quantile', 'affine', 'saturating', 'none')
# where a curve is set: on the APL-on (wild-type) condition, applied unchanged to APL off ('on'), or on each condition
# against its own cohort ('own'; APL off outside the APL16 panel, which has no such flies, takes the APL-on curve); per
# wiring ('wiring', each wiring a fly with its own indicator scale) or on the family's pooled values ('family')
CURVE_SET_ON: Tuple[str, ...] = ('on', 'own')
CURVE_SCOPES: Tuple[str, ...] = ('wiring', 'family')
N_CURVE_KNOTS = 201
# (not in the lock's readings: implementation details fixed here before any model was read, `implementation_constants`)
CURVE_TOP_SHARE = 0.1
N_NOISE_AWARE_ITERS = 12
# the flies' weight in the pooled distribution: each fly counts equally ('fly'), or each unit and odor ('entry')
FLY_WEIGHTING: Tuple[str, ...] = ('fly', 'entry')
# the measurement noise: 'residual' (the flies' trial-to-trial residual SD in baseline-noise units, pooled over flies,
# units and odors, with the n / (n - 1) correction); 'residual_binned' (the same, as a function of the trial mean, in
# bins of the flies' pooled trial means); 'reliability' (E90's and E119's: per odor, the SD that gives the flies' trial-to-
# trial reliability across units); 'none'
NOISE_MODELS: Tuple[str, ...] = ('residual', 'residual_binned', 'reliability', 'none')
N_NOISE_BINS = 10
# the response criterion: the flies' for the stage, applied after the curve ('flies'), or the model's own, two or more
# spikes, on raw counts ('model_spikes'; the KC somata only)
CRITERIA: Tuple[str, ...] = ('flies', 'model_spikes')
MODEL_SPIKES = 2.0
POOL_SEED = 0
NOISE_SEED = 0
# where the flies have no single trials (viki's `has_trials` False: the mixtures' KC somata, olfdb's trial means only),
# their measurement noise cannot be read from them: no noise is added there, and the provenance says so
NO_TRIALS_RULE = 'no noise where the flies have no single trials'
# viki's committed copy of the lock (data_prep/stage_population_data.py LOCK_COPY)
LOCK_COPY = H.HERE / 'experiments' / 'E134_lock.json'


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
            pooled distribution (fixed-point iterations at a fixed seed), not the values before the noise. Only with
            noise 'residual' or 'residual_binned' and curve 'quantile' or 'affine'.
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

    def __post_init__(self):
        bad = []
        if self.stage not in STAGES:
            bad.append(f'stage {self.stage!r}')
        for name, ok in (('curve', CURVES), ('curve_set_on', CURVE_SET_ON), ('curve_scope', CURVE_SCOPES),
                         ('fly_weighting', FLY_WEIGHTING), ('noise', NOISE_MODELS), ('criterion', CRITERIA)):
            if getattr(self, name) not in ok:
                bad.append(f'{name} {getattr(self, name)!r}')
        if self.criterion == 'model_spikes' and self.stage != 'kc_soma':
            bad.append("criterion 'model_spikes' outside the KC somata")
        if self.curve == 'none' and self.criterion == 'flies':
            bad.append("curve 'none' with the flies' criterion (the model's units are not the flies')")
        if self.stage == 'kc_claw' and self.pool not in CLAW_POOLS:
            bad.append(f'claw pool {self.pool!r}')
        if self.stage == 'pn_bouton' and self.pool not in BOUTON_POOLS:
            bad.append(f'bouton pool {self.pool!r}')
        if self.stage in SUBSETS and self.subset not in SUBSETS[self.stage]:
            bad.append(f'subset {self.subset!r}')
        if NAMED_UNITS[self.stage] and self.glomeruli not in GLOMERULI_RULES:
            bad.append(f'glomeruli {self.glomeruli!r}')
        if bad:
            raise ValueError(f'reading {self.name!r}: ' + ', '.join(bad))

    @property
    def needs_flies(self) -> bool:
        """Whether reading it needs the flies (for the curve, the noise or the criterion)."""
        return self.curve != 'none' or self.noise != 'none' or self.criterion == 'flies'


# each stage's reading and its alternatives, as viki's lock records them (607da67, experiments/E134_lock.json
# 'readings'; `lock_readings_agree` checks it): the first alternative is the stage's second reading, the one a "not
# readable" verdict compares with. They are not changed after the lock
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
            fits), 'E90:B5' (the 2 s pulse), 'E107:B' (somata only), 'E134:step3' (every stage the model has, on its
            panel's input).
        tag: the run's tag in its folder (E90's wiring label or D1 model; E107's / E134's '<panel>_<reading>_<model>').
        model: the model ('C2' for the family, 'fc' for E90's D1 fits; a beside model's name,
            experiments/E134_population_baseline.py BESIDE_MODELS).
        wiring: the wiring label (E83's ALL_WIRINGS for the family; the model's own otherwise).
        input_panel: the panel of the run's ORN input ('apl16' for E90's runs).
        conditions: 'on' / 'off' -> the saved matrices, E90's layout: 'pn_win', 'pn_pre', 'bouton_win', 'bouton_pre',
            'claw_win', 'claw_pre' (units x odors, odors abbreviated; the pre-odor means at the settled level where
            `settled`, E90's then kept as '<stage>_pre_e90'), 'apl' (APL's window and pre-odor means per odor),
            'spike_counts', 'responses' (KCs x odors, as saved). A somata-only run has 'spike_counts' and 'responses'
            only.
        meta: the units' metadata, E90's units.json: 'bouton_glomerulus', 'claw_glomerulus', 'kc_of_claw'; and
            'claw_bouton' (each claw's bouton, `claw_boutons`) where known.
        orn: the run's ORN input (glomeruli x odors, odors abbreviated), for the 'orn' stage; None if not known.
        files: every file read: dicts of path, bytes, sha256.
        settled: whether the pre-odor means are at the settled level.
        info: anything else the run's folder records about it (step 3's <tag>_run.json).
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
    info: Dict[str, object] = field(default_factory=dict)

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
        x: the curve's knots in the model's values (increasing; 'quantile' only).
        y: its values at the knots (non-decreasing); linear between knots, flat below the first, and above the last
            continued with slope `params['upper_slope']`.
        params: the parametric curves' parameters ('affine': scale, offset; 'saturating': scale, K; 'quantile':
            upper_slope).
        set_on: what it was set on: the number of model and fly values, the flies' weighting, the sha256 of both
            pooled arrays, the noise of a noise-aware fit.
    """
    kind: str
    x: np.ndarray = field(default_factory=lambda: np.zeros(0))
    y: np.ndarray = field(default_factory=lambda: np.zeros(0))
    params: Dict[str, float] = field(default_factory=dict)
    set_on: Dict[str, object] = field(default_factory=dict)

    def __call__(self, values) -> np.ndarray:
        """The curve applied to an array of model values (any shape; NaN stays NaN)."""
        x = np.asarray(values, dtype=float)
        if self.kind == 'none':
            return x.copy()
        if self.kind == 'affine':
            return self.params['scale'] * x + self.params['offset']
        if self.kind == 'saturating':
            return self.params['scale'] * x / (np.abs(x) + self.params['K'])
        if self.kind == 'quantile':
            y = np.interp(x, self.x, self.y)
            hi = x > self.x[-1]
            y[hi] = self.y[-1] + float(self.params.get('upper_slope', 0.0)) * (x[hi] - self.x[-1])
            return np.where(np.isnan(x), np.nan, y)
        raise ValueError(self.kind)

    def as_dict(self) -> dict:
        """JSON-able form (knots as lists)."""
        return dict(kind=self.kind, x=[float(v) for v in self.x], y=[float(v) for v in self.y],
                    params={k: float(v) for k, v in self.params.items()}, set_on=dict(self.set_on))


@dataclass
class NoiseModel:
    """The flies' measurement noise of one stage and panel, in their units (baseline-noise SDs).

    Attributes:
        kind: one of NOISE_MODELS.
        n_trials: trials per odor of the recordings (the model's trials are averaged over as many).
        sd: the single-trial noise SD ('residual').
        bin_edges, bin_sd: the single-trial noise SD as a function of the trial mean ('residual_binned'): bin k holds
            the trial means from bin_edges[k] to bin_edges[k + 1] (the outer bins open).
        reliability: per odor, the flies' trial-to-trial reliability across units ('reliability').
        from_flies: the flies and counts it was estimated from.
    """
    kind: str
    n_trials: int
    sd: float = np.nan
    bin_edges: Optional[np.ndarray] = None
    bin_sd: Optional[np.ndarray] = None
    reliability: Optional[pd.Series] = None
    from_flies: Dict[str, object] = field(default_factory=dict)

    def trial_sd(self, values: np.ndarray) -> np.ndarray:
        """The single-trial SD at each value ('residual' and 'residual_binned'; the binned SD looked up at the value)."""
        v = np.asarray(values, dtype=float)
        if self.kind == 'residual':
            return np.full(v.shape, float(self.sd))
        if self.kind == 'residual_binned':
            k = np.clip(np.searchsorted(self.bin_edges[1:-1], v, side='right'), 0, len(self.bin_sd) - 1)
            return np.asarray(self.bin_sd, dtype=float)[k]
        raise ValueError(f'{self.kind}: no SD by value')

    def as_dict(self) -> dict:
        return dict(kind=self.kind, n_trials=int(self.n_trials), sd=float(self.sd),
                    bin_edges=None if self.bin_edges is None else [float(v) for v in self.bin_edges],
                    bin_sd=None if self.bin_sd is None else [float(v) for v in self.bin_sd],
                    reliability=None if self.reliability is None else {str(k): float(v) for k, v in
                                                                        self.reliability.items()},
                    from_flies=dict(self.from_flies))


@dataclass
class ModelReading:
    """One model wiring's populations, read: per (stage, panel, condition) hal's `P.Population` and its provenance.

    Attributes:
        model, wiring: the model and wiring read.
        populations: (stage, panel, condition) -> the population (values in the flies' units after the curve and the
            noise; `responds` at the reading's criterion).
        before_curve: the same keys -> the model's sampled and pooled units x odors values before the curve.
        curves: (stage, panel, condition) -> the curve used (the same object for both conditions where it is set on APL
            on).
        provenance: the same keys as `populations` -> `reading_provenance`'s dict.
        skipped: (stage, panel) -> why that stage and panel were not read (no flies to read it against).
    """
    model: str
    wiring: str
    populations: Dict[Tuple[str, str, str], P.Population] = field(default_factory=dict)
    before_curve: Dict[Tuple[str, str, str], pd.DataFrame] = field(default_factory=dict)
    curves: Dict[Tuple[str, str, str], StageCurve] = field(default_factory=dict)
    provenance: Dict[Tuple[str, str, str], dict] = field(default_factory=dict)
    skipped: Dict[Tuple[str, str], str] = field(default_factory=dict)


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


def implementation_constants() -> dict:
    """The constants of the readings' implementation that the lock (written from `reading_definitions` before the
    implementation) does not record, fixed before any model was read: the slope of a quantile curve above its last
    knot, the noise-aware fit's iterations (no locked reading is noise-aware), the rule for flies without single
    trials."""
    return dict(CURVE_TOP_SHARE=CURVE_TOP_SHARE, N_NOISE_AWARE_ITERS=N_NOISE_AWARE_ITERS,
                no_trials_rule=NO_TRIALS_RULE)


def pooled_values(tables: Sequence[pd.DataFrame], odors: Optional[Sequence[str]] = None) -> np.ndarray:
    """Every finite entry of the tables (units x odors; restricted to `odors` where given), as one flat array: the
    pooled distribution, with no odor label left."""
    out = []
    for t in tables:
        x = (t[list(odors)] if odors is not None else t).to_numpy(dtype=float).ravel()
        out.append(x[np.isfinite(x)])
    return np.concatenate(out) if out else np.zeros(0)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def file_record(path: Path) -> dict:
    """{'path', 'bytes', 'sha256'} of a file read."""
    path = Path(path)
    return dict(path=str(path), bytes=int(path.stat().st_size), sha256=_sha256_file(path))


def array_sha256(x) -> str:
    """sha256 of an array's float64 bytes (C order)."""
    return hashlib.sha256(np.ascontiguousarray(np.asarray(x, dtype='<f8')).tobytes()).hexdigest()


def _abbrev(df: pd.DataFrame) -> pd.DataFrame:
    """A copy with odor columns abbreviated (H.odor_abbrev: '2h @ -3' -> '2h'; a (panel, odor) pair -> its odor's)."""
    df = df.copy()
    df.columns = [H.odor_abbrev(c[-1] if isinstance(c, tuple) else c) for c in df.columns]
    return df


def _archive_dirs(experiment: str, part: str) -> List[Path]:
    """The local outputs folder of an experiment's part, then its model-output archive folder (when the shared directory
    is mounted)."""
    dirs = [H.HERE / 'outputs' / experiment / part]
    try:
        dirs.append(H.shared_dir() / 'data' / 'model_outputs' / experiment / part)
    except FileNotFoundError:
        pass
    return dirs


def _committed(path: Path) -> bool:
    """Whether `path` is tracked by this repository's git and has no uncommitted change."""
    try:
        r = subprocess.run(['git', '-C', str(H.HERE), 'ls-files', '--error-unmatch', str(path)], capture_output=True)
        if r.returncode != 0:
            return False
        r = subprocess.run(['git', '-C', str(H.HERE), 'diff', '--quiet', 'HEAD', '--', str(path)], capture_output=True)
        return r.returncode == 0
    except OSError:
        return False


def require_lock() -> dict:
    """The committed lock, before any model is read against the flies: experiments/E134_lock.json exists, is tracked
    with no uncommitted change, and passes viki's check (`stage_population_data.require_lock`). Raises PermissionError
    otherwise.

    It does not yet check that the lock's readings of the model are `reading_definitions()`: viki's lock format is not
    known yet (phase 2 of 2026-11-24 was written before it); `lock_readings_agree` does it once the lock holds them.
    """
    if not LOCK_COPY.exists():
        raise PermissionError(f'{LOCK_COPY.name}: no lock yet (viki writes it, data_prep/stage_population_data.py '
                              '--lock); no model is read before it is committed')
    if not _committed(LOCK_COPY):
        raise PermissionError(f'{LOCK_COPY.name} is not committed (or has uncommitted changes)')
    from data_prep import stage_population_data as D
    try:
        return D.require_lock()
    except NotImplementedError:
        raise PermissionError("viki's lock check (stage_population_data.require_lock) is not implemented yet") from None
    except RuntimeError as e:
        raise PermissionError(f'the lock check failed: {e}') from None


def lock_readings_agree(lock: Mapping) -> Optional[bool]:
    """Whether the lock's readings of the model (its 'readings', viki's `model_readings()`) are
    `reading_definitions()`; None when the lock holds none."""
    got = lock.get('readings') if isinstance(lock, Mapping) else None
    if got is None:
        return None
    canon = lambda o: json.dumps(o, sort_keys=True, default=str)
    return canon(got) == canon(reading_definitions())


# ---------------------------------------------------------------------------
# the saved runs

def _find_dir(experiment: str, part: str, probe: str) -> Path:
    for d in _archive_dirs(experiment, part):
        if (d / probe).exists():
            return d
    raise FileNotFoundError(f'{probe} not in {[str(d) for d in _archive_dirs(experiment, part)]}')


def _read_stage_files(d: Path, prefix: str, with_settled: bool = False) -> Tuple[Dict[str, pd.DataFrame], List[dict]]:
    """E90's layout of one condition (`_save_stage`): <prefix>_<pn|bouton|claw>_<win|pre>.parquet, <prefix>_apl.csv,
    <prefix>_responses.parquet, <prefix>_spike_counts.parquet; with `with_settled` also <prefix>_<stage>_pre_settled
    (step 3's), where present. Stage matrices absent from the folder are left out."""
    out, files = {}, []
    for k in STAGE_FILES:
        for w in ('win', 'pre') + (('pre_settled',) if with_settled else ()):
            f = d / f'{prefix}_{k}_{w}.parquet'
            if f.exists():
                out[f'{k}_{w}'] = pd.read_parquet(f)
                files.append(file_record(f))
    f = d / f'{prefix}_apl.csv'
    if f.exists():
        out['apl'] = pd.read_csv(f, index_col=0)
        files.append(file_record(f))
    for k in ('responses', 'spike_counts'):
        f = d / f'{prefix}_{k}.parquet'
        out[k] = pd.read_parquet(f)
        files.append(file_record(f))
    for k in list(out):
        if k.split('_')[0] in STAGE_FILES and k.endswith(('_win', '_pre', '_pre_settled')):
            out[k] = _abbrev(out[k])
    return out, files


def settled_levels() -> Tuple[pd.DataFrame, dict]:
    """E90_settled_preodor's preodor_levels.parquet (every unit's pre-odor level over E90's window and the settled one,
    per part, tag, condition and stage) and its file record."""
    p = H.shared_dir().joinpath(*SETTLED_LEVELS)
    return pd.read_parquet(p), file_record(p)


def apply_settled(st: Dict[str, pd.DataFrame], levels: pd.DataFrame, part: str, tag: str, cc: str) -> None:
    """Replaces a saved condition's pre-odor means by the settled levels, in place, as E119's `model_stage` does: for
    each stage the levels' E90-window column must equal the saved pre-odor means (LEVEL_TOL, relative to the largest
    value), and every odor then gets the unit's settled level; E90's means are kept as '<stage>_pre_e90'."""
    lv = levels[(levels.part == part) & (levels.tag == tag) & (levels.cond == cc)]
    for k in STAGE_FILES:
        if f'{k}_pre' not in st:
            continue
        x = st[f'{k}_pre']
        lk = lv[lv.stage == k].sort_values('unit')
        if len(lk) != x.shape[0]:
            raise RuntimeError(f'{part} {tag} {cc} {k}: {len(lk)} settled levels for {x.shape[0]} units')
        scale = max(1.0, float(np.abs(x.to_numpy(dtype=float)).max()))
        if np.max(np.abs(lk.e90.to_numpy() - x.iloc[:, 0].to_numpy(dtype=float))) > LEVEL_TOL * scale:
            raise RuntimeError(f"{part} {tag} {cc} {k}: the levels' E90 window is not the saved pre-odor mean")
        st[f'{k}_pre_e90'] = x
        st[f'{k}_pre'] = pd.DataFrame(np.repeat(lk.settled.to_numpy(dtype=float)[:, None], x.shape[1], axis=1),
                                      index=x.index, columns=x.columns)


def _units_meta(d: Path, wiring: str) -> Tuple[Dict[str, np.ndarray], dict]:
    f = d / f'{wiring}_units.json'
    if not f.exists():
        f = _find_dir(E90_EXP, 'B', f'{wiring}_units.json') / f'{wiring}_units.json'
    return {k: np.asarray(v) for k, v in json.loads(f.read_text()).items()}, file_record(f)


def check_meta(run: SavedRun) -> None:
    """The units' metadata are the run's units: one bouton glomerulus per bouton, one claw glomerulus and KC per claw,
    the claws' KCs those of the claw matrices' index (E90's check); RuntimeError otherwise."""
    m = run.meta
    for cond, st in run.conditions.items():
        if 'claw_win' in st:
            kc = np.asarray(st['claw_win'].index.get_level_values(H.KC_ID_LEVEL)).astype(np.int64)
            if not np.array_equal(kc, np.asarray(m['kc_of_claw']).astype(np.int64)):
                raise RuntimeError(f"{run.tag} {cond}: the claws are not those of the units' metadata")
            if len(m['claw_glomerulus']) != len(kc) or ('claw_bouton' in m and len(m['claw_bouton']) != len(kc)):
                raise RuntimeError(f"{run.tag} {cond}: the claws' glomeruli / boutons do not match the claws")
        if 'bouton_win' in st and len(st['bouton_win']) != len(m['bouton_glomerulus']):
            raise RuntimeError(f"{run.tag} {cond}: the boutons are not those of the units' metadata")


def claw_boutons(W: pd.DataFrame) -> np.ndarray:
    """Each claw's bouton (the column of its largest PN->claw weight; the family's claws have one bouton each) as
    '<glomerulus>/<pn_id>/<bouton_id>' labels, in the wiring's claw order."""
    cols = W.columns
    j = np.asarray(W.to_numpy(dtype=float).argmax(axis=1))
    lab = np.asarray(['/'.join(str(v) for v in (c if isinstance(c, tuple) else (c,))) for c in cols])
    return lab[j]


def load_e90_run(tag: str, part: str = 'B', settled: bool = True, with_claw_bouton: bool = False) -> SavedRun:
    """A run E90 saved with its stage matrices (local outputs/E90_model_stages/<part>/, else the archive
    data/model_outputs/E90_model_stages/<part>/), both conditions, with the units' metadata (the wiring's
    <wiring>_units.json in part B's folder) and the lab's input (H.apl16_orn_deltas, odors abbreviated) as its ORN stage.

    Args:
        tag: the wiring (part 'B': E83's ALL_WIRINGS), the firing-ceiling fit (part 'D1': E90's D1_MODELS) or
            'real_pulse2' (part 'B5').
        part: E90's part folder.
        settled: replace the pre-odor means by the settled levels (`apply_settled`, with its check).
        with_claw_bouton: add each claw's bouton to the metadata (`claw_boutons` of E83's `wiring_parts(...)['W']`,
            which may simulate the bouton base's APL-off run once where its cache is missing), for the 'b14' pool.
    """
    d = _find_dir(E90_EXP, part, f'{tag}_on_claw_win.parquet')
    conds, files = {}, []
    for c in CONDITIONS:
        conds[c], f = _read_stage_files(d, f'{tag}_{CONDITION_FILE[c]}')
        files += f
    if part == 'B':
        wiring, model = tag, 'C2'
    elif part == 'D1':
        wiring, model = E90_D1_WIRING[tag], 'fc'
    elif part == 'B5':
        wiring, model = 'real', 'C2'
    else:
        raise ValueError(f'E90 part {part!r}')
    meta, mf = _units_meta(d.parent / 'B', wiring)
    files.append(mf)
    if with_claw_bouton:
        from experiments import E83_cross_wiring_fit as E83
        meta['claw_bouton'] = claw_boutons(E83.wiring_parts(wiring)['W'])
    if settled:
        lv, lf = settled_levels()
        for c in CONDITIONS:
            apply_settled(conds[c], lv, part, tag, CONDITION_FILE[c])
        files.append(lf)
    run = SavedRun(source=f'E90:{part}', tag=tag, model=model, wiring=wiring, input_panel='apl16', conditions=conds,
                   meta=meta, orn=_abbrev(H.apl16_orn_deltas()), files=files, settled=settled)
    check_meta(run)
    return run


def load_e107_run(panel: str, reading: str, wiring: str) -> SavedRun:
    """E107's saved spike counts and responses of one wiring on 'megamat' or 'validation2', 'fitted' or 'held'
    (outputs/E107_validation2_panel/B/<panel>_<reading>_<wiring>_<on|acute_off>_<responses|spike_counts>.parquet,
    else the archive), as a somata-only run, with its input (E107's `panel_input`, odors abbreviated) as its ORN
    stage."""
    from experiments import E107_validation2_panel as E107
    tag = f'{panel}_{reading}_{wiring}'
    d = _find_dir(E107_EXP, 'B', f'{tag}_on_spike_counts.parquet')
    conds, files = {}, []
    for c in CONDITIONS:
        st = {}
        for k in ('responses', 'spike_counts'):
            f = d / f'{tag}_{CONDITION_FILE[c]}_{k}.parquet'
            st[k] = pd.read_parquet(f)
            files.append(file_record(f))
        conds[c] = st
    return SavedRun(source='E107:B', tag=tag, model='C2', wiring=wiring, input_panel=panel, conditions=conds,
                    orn=_abbrev(E107.panel_input(panel)), files=files, settled=True, info=dict(reading=reading))


def step3_dirs() -> List[Path]:
    """Where step 3's runs are: outputs/E134_population_baseline/step3/, then the archive's."""
    return _archive_dirs(STEP3_EXP, STEP3_PART)


def save_step3_run(tag: str, conditions: Mapping[str, Mapping[str, pd.DataFrame]], meta: Mapping[str, np.ndarray],
                   orn: pd.DataFrame, info: Mapping[str, object], dirs: Sequence[Path]) -> List[Path]:
    """Writes one step-3 run in this module's layout, into each of `dirs` (atomically, file by file): per condition
    c ('on' / 'acute_off') <tag>_<c>_<stage>_<win|pre|pre_settled>.parquet for the stages it has, <tag>_<c>_apl.csv,
    <tag>_<c>_responses.parquet and <tag>_<c>_spike_counts.parquet; <tag>_units.json (the units' metadata),
    <tag>_orn.parquet (the input, odors abbreviated) and <tag>_run.json (`info`: model, wiring, panel, reading, ...,
    written last, so its presence says the run is complete)."""
    written = []

    def _put(p: Path, write) -> None:
        tmp = p.with_name(f'.{p.name}.{os.getpid()}.tmp')
        write(tmp)
        os.replace(tmp, p)
        written.append(p)

    for d in dirs:
        d = Path(d)
        d.mkdir(parents=True, exist_ok=True)
        for c, st in conditions.items():
            cc = CONDITION_FILE.get(c, c)
            for k, x in st.items():
                if k == 'apl':
                    _put(d / f'{tag}_{cc}_apl.csv', lambda p, x=x: x.to_csv(p))
                elif k in ('responses', 'spike_counts') or (k.split('_')[0] in STAGE_FILES
                                                           and k.endswith(('_win', '_pre', '_pre_settled'))):
                    y = x.copy()
                    y.columns = [str(col) for col in y.columns]
                    if isinstance(y.index, pd.MultiIndex):
                        y.index = y.index.set_names([str(n) for n in y.index.names])
                    _put(d / f'{tag}_{cc}_{k}.parquet', lambda p, y=y: y.to_parquet(p))
        _put(d / f'{tag}_units.json',
             lambda p: p.write_text(json.dumps({k: [v.item() if hasattr(v, 'item') else v for v in np.asarray(a)]
                                                for k, a in meta.items()}) + '\n'))
        o = _abbrev(orn)
        _put(d / f'{tag}_orn.parquet', lambda p, o=o: o.to_parquet(p))
        _put(d / f'{tag}_run.json',
             lambda p: p.write_text(json.dumps(dict(info), indent=1, default=lambda v: v.item() if hasattr(v, 'item')
                                               else str(v)) + '\n'))
    return written


def load_step3_run(tag: str, settled: bool = True, dirs: Optional[Sequence[Path]] = None) -> SavedRun:
    """A run of E134's step 3 (`save_step3_run`'s layout; locally, else the archive): every stage it saved, the
    settled pre-odor means ('<stage>_pre_settled') in place of E90's window where `settled` (E90's kept as
    '<stage>_pre_e90'). FileNotFoundError if its <tag>_run.json is in none of the folders."""
    dirs = list(dirs) if dirs is not None else step3_dirs()
    d = next((Path(x) for x in dirs if (Path(x) / f'{tag}_run.json').exists()), None)
    if d is None:
        raise FileNotFoundError(f'{tag}_run.json not in {[str(x) for x in dirs]}')
    info = json.loads((d / f'{tag}_run.json').read_text())
    conds, files = {}, [file_record(d / f'{tag}_run.json')]
    for c in CONDITIONS:
        st, f = _read_stage_files(d, f'{tag}_{CONDITION_FILE[c]}', with_settled=True)
        files += f
        if settled:
            for k in STAGE_FILES:
                if f'{k}_pre_settled' in st:
                    st[f'{k}_pre_e90'] = st[f'{k}_pre']
                    st[f'{k}_pre'] = st.pop(f'{k}_pre_settled')
                elif f'{k}_pre' in st:
                    raise RuntimeError(f'{tag} {c} {k}: no settled pre-odor means saved')
        conds[c] = st
    meta = {k: np.asarray(v) for k, v in json.loads((d / f'{tag}_units.json').read_text()).items()}
    files.append(file_record(d / f'{tag}_units.json'))
    orn = pd.read_parquet(d / f'{tag}_orn.parquet')
    files.append(file_record(d / f'{tag}_orn.parquet'))
    run = SavedRun(source='E134:step3', tag=tag, model=str(info['model']), wiring=str(info['wiring']),
                   input_panel=str(info['panel']), conditions=conds, meta=meta, orn=orn, files=files, settled=settled,
                   info=info)
    check_meta(run)
    return run


def run_odors(run: SavedRun) -> Tuple[str, ...]:
    """The run's odors (its spike counts' columns, abbreviated), in its order."""
    return tuple(_abbrev(run.conditions['on']['spike_counts'].iloc[:0]).columns)


def run_panels(run: SavedRun, stage: str) -> List[Tuple[str, Tuple[str, ...]]]:
    """The (panel, odors) a run's stage is read on: E90's APL16 runs give the PN, bouton and claw stages and the ORN
    input on 'megamat' (its 12 odors, H.MEGAMAT12) and the somata on 'apl16' (16) and 'megamat' (12); a mixture run
    ('mixtures') each mixture panel ('kmix', 'cmix') on its stimuli; every other run its input's panel. The odors are
    those of the run that are the panel's (viki's `panel_odors`), in the panel's order (the comparison with the flies
    is on the odors both sides have)."""
    if run.input_panel == 'apl16':
        if stage == 'kc_soma':
            return [('apl16', tuple(H.APL16_ODORS)), ('megamat', tuple(H.MEGAMAT12))]
        return [('megamat', tuple(H.MEGAMAT12))]
    from data_prep import stage_population_data as D
    have = set(run_odors(run))
    panels = MIXTURE_PANELS if run.input_panel == 'mixtures' else (run.input_panel,)
    out = []
    for p in panels:
        od = tuple(o for o in D.panel_odors(p) if o in have)
        if od:
            out.append((p, od))
    return out


# ---------------------------------------------------------------------------
# the units, sampled and pooled as recorded

def stage_changes(run: SavedRun, stage: str, condition: str, odors: Sequence[str],
                  reading: Optional[Reading] = None) -> pd.DataFrame:
    """A run's units x odors before sampling and pooling: 'orn' the input deltas; 'pn_dendrite', 'pn_bouton',
    'kc_claw' the window mean minus the pre-odor mean (E90's `stage_units`; at the settled level where the run was
    loaded so), the claws as E90's `claw_signal` at the reading's spike share (scales from the run's APL-on condition);
    'kc_soma' the spike counts. Odors abbreviated, restricted to `odors` in that order; the units' own index kept."""
    odors = list(odors)
    if stage == 'orn':
        if run.orn is None:
            raise KeyError(f'{run.tag}: no ORN input')
        return run.orn[odors].astype(float)
    st = run.conditions[condition]
    if stage == 'kc_soma':
        return _abbrev(st['spike_counts'].astype(float))[odors]
    from experiments import E90_model_stages as E90
    k = E90_STAGE[stage]
    X = E90.stage_units(st, k)
    share = float(reading.claw_spike_share) if (reading is not None and stage == 'kc_claw') else 0.0
    if share:
        on = run.conditions['on']
        scales = (E90._rms(E90.stage_units(on, 'claw').to_numpy()), E90._rms(on['spike_counts'].to_numpy(dtype=float)))
        X = E90.claw_signal(X, st['spike_counts'], np.asarray(run.meta['kc_of_claw']), share, scales)
    return X[odors].astype(float)


def driver_glomeruli(stage: str = 'pn_dendrite', panel: str = 'megamat', rule: str = 'every_fly',
                     flies: Optional['D.Population'] = None,
                     model_glomeruli: Optional[Sequence[str]] = None) -> List[str]:
    """The glomeruli a named stage is read on: those the recorded flies of the stage and panel have (`flies`, viki's
    population; a fly has a glomerulus where its row has a finite value), every fly or any fly (`rule`), or the 28
    shared ones (E90's `shared_glomeruli`), each only if the model has it (`model_glomeruli`), sorted."""
    if rule == 'shared28':
        from experiments import E90_model_stages as E90
        gl = set(E90.shared_glomeruli())
    elif rule in ('every_fly', 'any_fly'):
        if flies is None:
            flies = fly_side(stage, panel, 'on', trials=False)
            if flies is None:
                raise ValueError(f'{stage} {panel}: no recorded flies to take the glomeruli from')
        sets = [set(t.index[np.isfinite(t.to_numpy(dtype=float)).any(axis=1)].astype(str)) for t in flies.tables.values()]
        gl = set.intersection(*sets) if rule == 'every_fly' else set.union(*sets)
    else:
        raise ValueError(f'{rule=}')
    if model_glomeruli is not None:
        gl &= set(str(g) for g in model_glomeruli)
    return sorted(gl)


def sample_driver_glomeruli(values: pd.DataFrame, glomeruli: Sequence[str]) -> pd.DataFrame:
    """The named units of `values` (glomeruli x odors) that are in `glomeruli`, in that order; KeyError if the model
    lacks one (the caller intersects first, `driver_glomeruli`)."""
    missing = [g for g in glomeruli if g not in values.index]
    if missing:
        raise KeyError(f'glomeruli not in the model: {missing}')
    return values.loc[list(glomeruli)]


def _pool(U: pd.DataFrame, groups: Optional[np.ndarray], pool: str, seed: int) -> pd.DataFrame:
    """E119's `pool_units` ('b<n>' is its 'g<n>' with the claws' boutons as the groups), silent regions dropped (E90's
    `_drop_silent`), the regions indexed 0..n-1 ('region')."""
    from experiments import E90_model_stages as E90
    from experiments import E119_stages_per_roi as E119
    pool = str(pool)
    lab = 'g' + pool[1:] if pool.startswith('b') else pool
    V = E90._drop_silent(E119.pool_units(U, lab, seed, None if lab in ('1',) or lab.startswith('r') else groups))
    V = V.reset_index(drop=True)
    V.index.name = 'region'
    return V


def pool_boutons(values: pd.DataFrame, bouton_glomerulus: np.ndarray, pool: str = '1', subset: str = 'shared',
                 seed: int = POOL_SEED, subset_glomeruli: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """The model's boutons as the recordings' bouton regions: those of `subset` ('shared': `subset_glomeruli`, by default
    the 28 glomeruli of E90's `shared_glomeruli`; 'driver': `subset_glomeruli`, required; 'all'), pooled by `pool`
    (BOUTON_POOLS; E119's `pool_units`, 'g<n>' within a glomerulus), silent regions dropped, indexed 0..n-1 with no
    glomerulus label."""
    g = np.asarray(bouton_glomerulus).astype(str)
    if len(g) != len(values):
        raise ValueError(f'{len(g)} bouton glomeruli for {len(values)} boutons')
    if subset == 'all':
        m = np.ones(len(g), dtype=bool)
    elif subset in ('shared', 'driver'):
        if subset_glomeruli is None:
            if subset == 'driver':
                raise ValueError("subset 'driver' needs the driver line's glomeruli")
            from experiments import E90_model_stages as E90
            subset_glomeruli = E90.shared_glomeruli()
        m = np.isin(g, [str(x) for x in subset_glomeruli])
    else:
        raise ValueError(f'{subset=}')
    return _pool(values.iloc[np.flatnonzero(m)], g[m], pool, seed)


def pool_claws(values: pd.DataFrame, meta: Mapping[str, np.ndarray], pool: str = 'g14', subset: str = 'all',
               seed: int = POOL_SEED, subset_glomeruli: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """The model's claws as claw regions of about CLAW_REGION_SIZE claws (CLAW_POOLS): 'g<n>' within each glomerulus
    (`meta['claw_glomerulus']`; E119's `pool_units`, a remainder of at least n / 2 kept as a smaller region), 'b<n>'
    within each bouton (`meta['claw_bouton']`, `claw_boutons`), 'r<n>' at random, '1' none; of `subset` ('all', or
    'shared': the claws of `subset_glomeruli`, by default E90's 28); silent regions dropped, indexed 0..n-1, not
    assigned to a KC."""
    cg = np.asarray(meta['claw_glomerulus']).astype(str)
    if len(cg) != len(values):
        raise ValueError(f'{len(cg)} claw glomeruli for {len(values)} claws')
    if subset == 'all':
        m = np.ones(len(cg), dtype=bool)
    elif subset == 'shared':
        if subset_glomeruli is None:
            from experiments import E90_model_stages as E90
            subset_glomeruli = E90.shared_glomeruli()
        m = np.isin(cg, [str(x) for x in subset_glomeruli])
    else:
        raise ValueError(f'{subset=}')
    pool = str(pool)
    groups = None
    if pool.startswith('g'):
        groups = cg[m]
    elif pool.startswith('b'):
        if 'claw_bouton' not in meta:
            raise ValueError("pool 'b<n>' needs each claw's bouton (load the run with with_claw_bouton=True)")
        groups = np.asarray(meta['claw_bouton']).astype(str)[m]
    return _pool(values.iloc[np.flatnonzero(m)], groups, pool, seed)


def sample_units(run: SavedRun, stage: str, condition: str, odors: Sequence[str], reading: Reading,
                 seed_index: int = 0, glomeruli: Optional[Sequence[str]] = None,
                 subset_glomeruli: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """`stage_changes` then the stage's sampling or pooling under `reading` (`sample_driver_glomeruli` on `glomeruli`,
    required for the named stages, intersected with the model's; `pool_boutons`, `pool_claws` with `subset_glomeruli`;
    the somata as they are), with seed POOL_SEED + `seed_index` (the same in both conditions, so the same units form
    the same regions)."""
    X = stage_changes(run, stage, condition, odors, reading)
    seed = POOL_SEED + int(seed_index)
    if NAMED_UNITS[stage]:
        if glomeruli is None:
            raise ValueError(f"{stage}: the glomeruli to read it on are needed (`driver_glomeruli`)")
        have = set(X.index.astype(str))
        X.index = X.index.astype(str)
        return sample_driver_glomeruli(X, [g for g in glomeruli if g in have])
    if stage == 'pn_bouton':
        return pool_boutons(X, run.meta['bouton_glomerulus'], reading.pool, reading.subset, seed, subset_glomeruli)
    if stage == 'kc_claw':
        return pool_claws(X, run.meta, reading.pool, reading.subset, seed, subset_glomeruli)
    return X


# ---------------------------------------------------------------------------
# the curve and the noise

def pooled_fly_values(flies: 'D.Population', odors: Optional[Sequence[str]] = None,
                      weighting: str = 'fly') -> Tuple[np.ndarray, np.ndarray]:
    """The flies' pooled distribution of one stage and panel: every fly's trial means (viki's `Population.tables`) on
    `odors` (default all of each table's columns; odors a fly lacks are left out), flattened, and their weights ('fly':
    each fly sums to 1 / n_flies; 'entry': 1 each). Returns (values, weights); no odor label is kept."""
    vals, wts = [], []
    flies_kept = [f for f in flies.flies if f in flies.tables]
    for f in flies_kept:
        t = flies.tables[f]
        cols = list(t.columns) if odors is None else [o for o in odors if o in t.columns]
        x = t[cols].to_numpy(dtype=float).ravel()
        x = x[np.isfinite(x)]
        if not len(x):
            continue
        vals.append(x)
        wts.append(np.full(len(x), 1.0 / len(x) if weighting == 'fly' else 1.0))
    if weighting not in FLY_WEIGHTING:
        raise ValueError(f'{weighting=}')
    if not vals:
        raise ValueError(f'{flies.stage} {flies.panel}: no finite trial means on {odors}')
    w = np.concatenate(wts)
    return np.concatenate(vals), w / w.sum()


def _weighted_quantile_integral(values: np.ndarray, weights: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """(c, G): the cumulative weights c (from 0) of the sorted values and the integral G of their quantile function at
    those points, so that the mean of the quantile function over [a, b] is (G(b) - G(a)) / (b - a), G linear
    between the points."""
    o = np.argsort(values, kind='mergesort')
    v, w = np.asarray(values, dtype=float)[o], np.asarray(weights, dtype=float)[o]
    keep = w > 0
    v, w = v[keep], w[keep] / w[keep].sum()
    c = np.concatenate([[0.0], np.cumsum(w)])
    G = np.concatenate([[0.0], np.cumsum(v * w)])
    c[-1] = 1.0
    return c, G


def _interval_means(c: np.ndarray, G: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """The mean of a quantile function (`_weighted_quantile_integral`'s c, G) over each [lo, hi]."""
    a, b = np.interp(lo, c, G), np.interp(hi, c, G)
    width = np.asarray(hi, dtype=float) - np.asarray(lo, dtype=float)
    mid = 0.5 * (np.asarray(lo) + np.asarray(hi))
    # (a zero-width interval: the quantile function at that point)
    q = np.interp(mid, c[1:], np.diff(G) / np.diff(c)) if len(c) > 1 else np.zeros_like(mid)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(width > 0, (b - a) / np.where(width > 0, width, 1.0), q)


def _quantile_knots(model_values: np.ndarray, n_knots: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The knots of a 'quantile' curve: distinct model values at evenly spaced quantiles (values of the sample), and the
    share of the model's values below and at or below each (the quantile range each knot holds)."""
    ms = np.sort(model_values)
    n = len(ms)
    xk = np.unique(np.quantile(ms, np.linspace(0.0, 1.0, int(n_knots)), method='inverted_cdf'))
    lo = np.searchsorted(ms, xk, side='left') / n
    hi = np.searchsorted(ms, xk, side='right') / n
    return xk, lo, hi


def _upper_slope(x: np.ndarray, y: np.ndarray, p_mid: np.ndarray) -> float:
    """The slope of the curve's top CURVE_TOP_SHARE (from the knot at that quantile to the last), at least 0."""
    if len(x) < 2:
        return 0.0
    j = int(min(np.searchsorted(p_mid, 1.0 - CURVE_TOP_SHARE), len(x) - 2))
    dx = float(x[-1] - x[j])
    return max(float(y[-1] - y[j]) / dx, 0.0) if dx > 0 else 0.0


def _noise_draws(n: int, noise: NoiseModel, seed: int) -> np.ndarray:
    """Standard normal draws for the trial-mean noise of a noise-aware fit (fixed seed, one per model value)."""
    return np.random.default_rng(seed).standard_normal(n)


def fit_stage_curve(model_values: np.ndarray, fly_values: np.ndarray, fly_weights: Optional[np.ndarray] = None,
                    kind: str = 'quantile', n_knots: int = N_CURVE_KNOTS, noise: Optional[NoiseModel] = None,
                    seed: int = NOISE_SEED) -> StageCurve:
    """The stage's monotone curve from the model's pooled values to the flies' pooled distribution.

    It is set on the two pooled distributions only: the arguments are flat arrays, with no odor, unit or pair label, so
    nothing about odors can enter it (a permutation of either array gives the same curve).

    Args:
        model_values: the model's pooled values (every unit and odor of the condition the curve is set on; for scope
            'family', of every wiring).
        fly_values: the flies' pooled values (`pooled_fly_values`).
        fly_weights: their weights (default equal).
        kind: one of CURVES ('none' returns the identity).
        n_knots: the knots of a 'quantile' curve (at evenly spaced quantiles of the model's values; ties map to the
            flies' weighted mean over the quantile range they hold).
        noise: when given ('residual' or 'residual_binned'), the curve is set so that the model's values after the
            curve and this noise (averaged over its n_trials) match the flies' pooled distribution: for 'quantile'
            N_NOISE_AWARE_ITERS fixed-point steps on the knots' interval means (monotone after each step), for 'affine'
            the flies' variance less the noise's.
        seed: the noise draws of a noise-aware fit.
    """
    m = np.asarray(model_values, dtype=float).ravel()
    m = m[np.isfinite(m)]
    f = np.asarray(fly_values, dtype=float).ravel()
    w = np.ones_like(f) if fly_weights is None else np.asarray(fly_weights, dtype=float).ravel()
    ok = np.isfinite(f) & np.isfinite(w)
    f, w = f[ok], w[ok]
    if kind not in CURVES:
        raise ValueError(f'{kind=}')
    set_on = dict(n_model=int(len(m)), n_flies_values=int(len(f)), model_sha256=array_sha256(np.sort(m)),
                  flies_sha256=array_sha256(np.sort(f)), noise=None if noise is None else noise.as_dict())
    if kind == 'none':
        return StageCurve('none', set_on=set_on)
    if not len(m) or not len(f) or w.sum() <= 0:
        raise ValueError('a curve needs model and fly values')
    if noise is not None and noise.kind not in ('residual', 'residual_binned'):
        raise ValueError(f"a noise-aware curve takes 'residual' or 'residual_binned' noise, not {noise.kind!r}")
    w = w / w.sum()
    if kind == 'affine':
        mf = float(np.sum(w * f))
        vf = float(np.sum(w * (f - mf) ** 2))
        mm, sm = float(m.mean()), float(m.std())
        if noise is not None:
            vf = max(vf - float(np.mean((noise.trial_sd(np.full(1, mf)) ** 2))) / noise.n_trials, 0.0)
        scale = float(np.sqrt(vf) / sm) if sm > 0 else 0.0
        return StageCurve('affine', params=dict(scale=scale, offset=mf - scale * mm), set_on=set_on)
    if kind == 'saturating':
        if noise is not None:
            raise ValueError('a saturating curve is not fitted noise-aware')
        from scipy.optimize import minimize_scalar
        p = np.linspace(0.005, 0.995, 199)
        qm = np.quantile(m, p)
        c, G = _weighted_quantile_integral(f, w)
        qf = _interval_means(c, G, p, p)
        span = float(np.max(np.abs(qm))) or 1.0

        def sse(logK):
            g = qm / (np.abs(qm) + np.exp(logK))
            den = float(g @ g)
            a = float(g @ qf) / den if den > 0 else 0.0
            return float(np.sum((a * g - qf) ** 2)), a

        r = minimize_scalar(lambda t: sse(t)[0], bounds=(np.log(span * 1e-3), np.log(span * 1e3)), method='bounded')
        K = float(np.exp(r.x))
        return StageCurve('saturating', params=dict(scale=sse(r.x)[1], K=K), set_on=set_on)
    xk, lo, hi = _quantile_knots(m, n_knots)
    c, G = _weighted_quantile_integral(f, w)
    y = _interval_means(c, G, lo, hi)
    p_mid = 0.5 * (lo + hi)
    slope = _upper_slope(xk, y, p_mid)
    if noise is not None:
        # each knot compared over the quantile range from the midpoints to its neighbours' (at least its own range), so
        # that a comparison holds several of the model's values; monotone again by an isotonic fit after each step
        from sklearn.isotonic import IsotonicRegression
        mids = 0.5 * (p_mid[1:] + p_mid[:-1])
        L = np.maximum(np.minimum(lo, np.concatenate([[0.0], mids])), np.concatenate([[0.0], hi[:-1]]))
        U = np.minimum(np.maximum(hi, np.concatenate([mids, [1.0]])), np.concatenate([lo[1:], [1.0]]))
        target = _interval_means(c, G, L, U)
        e = _noise_draws(len(m), noise, seed)
        ms = np.sort(m)
        iso = IsotonicRegression(increasing=True)
        for _ in range(N_NOISE_AWARE_ITERS):
            cur = StageCurve('quantile', x=xk, y=y, params=dict(upper_slope=slope))
            v = cur(ms)
            z = v + e * noise.trial_sd(v) / np.sqrt(noise.n_trials)
            cz, Gz = _weighted_quantile_integral(z, np.ones_like(z))
            y = iso.fit_transform(np.arange(len(y)), y + (target - _interval_means(cz, Gz, L, U)),
                                  sample_weight=U - L)
            slope = _upper_slope(xk, y, p_mid)
    return StageCurve('quantile', x=xk, y=y, params=dict(upper_slope=slope), set_on=set_on)


def apply_curve(values: pd.DataFrame, curve: StageCurve) -> pd.DataFrame:
    """The curve applied entry by entry to a units x odors table (index and columns kept)."""
    return pd.DataFrame(curve(values.to_numpy(dtype=float)), index=values.index, columns=values.columns)


def _trial_blocks(trials: pd.DataFrame, odors: Optional[Sequence[str]]) -> Dict[str, np.ndarray]:
    """{odor: units x repeats} of one fly's single-trial table (columns (odor, repeat))."""
    have = list(dict.fromkeys(trials.columns.get_level_values(0)))
    use = have if odors is None else [o for o in odors if o in have]
    return {o: trials[o].to_numpy(dtype=float) for o in use}


def fly_noise(flies: 'D.Population', kind: str = 'residual', odors: Optional[Sequence[str]] = None) -> NoiseModel:
    """The flies' measurement noise of one stage and panel from their single trials (viki's `Population.trials`,
    loaded with `load_population(..., trials=True)`): the residual of each trial from its unit's trial mean for the
    odor, pooled over flies, units and odors with the n / (n - 1) correction ('residual'), binned by the trial mean in
    N_NOISE_BINS bins of equal count ('residual_binned'), or E119's per-odor reliability, the mean over pairs of repeats
    of the Pearson r across units, averaged over flies ('reliability'); n_trials the flies' median number of trials per
    unit and odor. 'none' gives a model with no noise. Raises ValueError where the population has no trials (the
    reading then needs noise 'none' and says so)."""
    if kind == 'none':
        return NoiseModel('none', n_trials=1)
    if kind not in NOISE_MODELS:
        raise ValueError(f'{kind=}')
    if not flies.trials:
        raise ValueError(f'{flies.stage} {flies.panel} {flies.cohort}: no single trials on hand for the noise')
    ss_all, dof_all, means_all, res2_all, dofs_all, n_rep = 0.0, 0, [], [], [], []
    rel: Dict[str, List[float]] = {}
    for fly in flies.flies:
        if fly not in flies.trials:
            continue
        for o, A in _trial_blocks(flies.trials[fly], odors).items():
            fin = np.isfinite(A)
            n = fin.sum(axis=1)
            rows = n >= 2
            if not rows.any():
                continue
            A2 = np.where(fin, A, np.nan)[rows]
            mu = np.nanmean(A2, axis=1)
            r2 = np.nansum((A2 - mu[:, None]) ** 2, axis=1)
            d = n[rows] - 1
            ss_all += float(r2.sum())
            dof_all += int(d.sum())
            means_all.append(mu)
            res2_all.append(r2)
            dofs_all.append(d)
            n_rep.append(n[rows])
            if kind == 'reliability':
                cs = [np.corrcoef(A[:, i], A[:, j])[0, 1] for i in range(A.shape[1]) for j in range(i + 1, A.shape[1])
                      if (np.isfinite(A[:, i]) & np.isfinite(A[:, j])).sum() > 2]
                if cs:
                    rel.setdefault(o, []).append(float(np.nanmean(cs)))
    if dof_all == 0:
        raise ValueError(f'{flies.stage} {flies.panel}: no unit with two or more trials of an odor')
    n_trials = int(round(float(np.median(np.concatenate(n_rep)))))
    src = dict(stage=flies.stage, panel=flies.panel, cohort=flies.cohort, flies=list(flies.flies), dof=int(dof_all))
    if kind == 'residual':
        return NoiseModel('residual', n_trials=n_trials, sd=float(np.sqrt(ss_all / dof_all)), from_flies=src)
    if kind == 'residual_binned':
        mu, r2, d = np.concatenate(means_all), np.concatenate(res2_all), np.concatenate(dofs_all)
        edges = np.unique(np.quantile(mu, np.linspace(0, 1, N_NOISE_BINS + 1)))
        k = np.clip(np.searchsorted(edges[1:-1], mu, side='right'), 0, len(edges) - 2)
        sd = np.array([np.sqrt(r2[k == b].sum() / max(d[k == b].sum(), 1)) for b in range(len(edges) - 1)])
        return NoiseModel('residual_binned', n_trials=n_trials, sd=float(np.sqrt(ss_all / dof_all)), bin_edges=edges,
                          bin_sd=sd, from_flies=src)
    return NoiseModel('reliability', n_trials=n_trials,
                      reliability=pd.Series({o: float(np.mean(v)) for o, v in rel.items()}), from_flies=src)


def add_measurement_noise(values: pd.DataFrame, noise: Optional[NoiseModel], rng: np.random.Generator,
                          return_trials: bool = False):
    """`values` (units x odors, in the flies' units) through the flies' measurement noise: `noise.n_trials` trials,
    each the value plus independent Gaussian noise of the model's kind ('residual': SD `sd`; 'residual_binned': the SD
    of the value's bin; 'reliability': per odor the SD giving that odor's reliability across the units,
    sigma^2 = variance over units x (1 - r) / r, E119's `noisy_trials`), and their mean (units x odors). With
    `return_trials` also the trials (trials x units x odors), for readings that remove the noise as the flies' are.
    None or 'none' returns `values` unchanged (one trial)."""
    X = values.to_numpy(dtype=float)
    if noise is None or noise.kind == 'none':
        out = values.copy()
        return (out, X[None]) if return_trials else out
    n = int(noise.n_trials)
    if noise.kind in ('residual', 'residual_binned'):
        sd = noise.trial_sd(X)
    elif noise.kind == 'reliability':
        r = noise.reliability.reindex([str(c) for c in values.columns]).to_numpy(dtype=float)
        with np.errstate(invalid='ignore', divide='ignore'):
            sd = np.broadcast_to(np.sqrt(np.nanvar(X, axis=0) * (1.0 - r) / r), X.shape)
        sd = np.where(np.isfinite(sd), sd, 0.0)
    else:
        raise ValueError(noise.kind)
    trials = X[None] + sd[None] * rng.standard_normal((n,) + X.shape)
    out = pd.DataFrame(trials.mean(axis=0), index=values.index, columns=values.columns)
    return (out, trials) if return_trials else out


# ---------------------------------------------------------------------------
# readings

def fly_side(stage: str, panel: str, condition: str, odors: Optional[Sequence[str]] = None,
             trials: bool = True, check_lock: bool = True) -> Optional['D.Population']:
    """The flies a model condition is read against (viki's `load_population(stage, panel, cohort_of(panel,
    condition), trials=trials)`), None where there are none (APL off outside the APL16 panel; a population whose files
    are missing). With `check_lock` (the default) it refuses before the lock is committed (`require_lock`)."""
    cohort = cohort_of(panel, condition)
    if cohort is None:
        return None
    if check_lock:
        require_lock()
    from data_prep import stage_population_data as D
    if (stage, panel, cohort) not in D.POPULATIONS:
        return None
    try:
        return D.load_population(stage, panel, cohort, trials=trials)
    except FileNotFoundError:
        return None


class FlyCache:
    """A mapping (stage, panel, condition) -> the flies (`fly_side`), loaded once each; `fixed` entries (planted
    populations, in tests) are returned as given, and with `fallback` False nothing else is loaded (None)."""

    def __init__(self, fixed: Optional[Mapping[Tuple[str, str, str], Optional['D.Population']]] = None,
                 check_lock: bool = True, fallback: bool = True):
        self.fixed = dict(fixed or {})
        self.check_lock = check_lock
        self.fallback = fallback
        self._got: Dict[Tuple[str, str, str], Optional['D.Population']] = {}

    def get(self, key: Tuple[str, str, str], default=None):
        if key in self.fixed:
            return self.fixed[key] if self.fixed[key] is not None else default
        if not self.fallback:
            return default
        if key not in self._got:
            stage, panel, condition = key
            self._got[key] = fly_side(stage, panel, condition, check_lock=self.check_lock)
        return self._got[key] if self._got[key] is not None else default


def _noise_cache_key(fl: 'D.Population', kind: str, odors: Sequence[str]) -> tuple:
    return (id(fl), kind, tuple(odors))


_NOISE_CACHE: Dict[tuple, NoiseModel] = {}


def _fly_noise_cached(fl: 'D.Population', kind: str, odors: Sequence[str]) -> NoiseModel:
    key = _noise_cache_key(fl, kind, odors)
    if key not in _NOISE_CACHE:
        _NOISE_CACHE[key] = fly_noise(fl, kind, odors)
    return _NOISE_CACHE[key]


def has_trials(flies: Optional['D.Population']) -> bool:
    """Whether the flies' single trials are on hand (loaded, and viki's `has_trials`)."""
    return flies is not None and bool(flies.trials) and bool(getattr(flies, 'has_trials', True))


def set_curve(model_values: np.ndarray, flies: 'D.Population', odors: Sequence[str], reading: Reading,
              seed_index: int = 0) -> StageCurve:
    """A reading's curve set on the model's pooled values against these flies' (on the odors both have)."""
    fo = [o for o in odors if any(o in t.columns for t in flies.tables.values())]
    fv, fw = pooled_fly_values(flies, fo, reading.fly_weighting)
    nm = (_fly_noise_cached(flies, reading.noise, fo) if reading.noise_aware_curve and has_trials(flies)
          and reading.noise in ('residual', 'residual_binned') else None)
    return fit_stage_curve(model_values, fv, fw, reading.curve, noise=nm, seed=NOISE_SEED + int(seed_index))


def finish_population(U: pd.DataFrame, stage: str, condition: str, reading: Reading, curve: StageCurve,
                      flies: Optional['D.Population'], odors: Sequence[str], seed_index: int,
                      name: str) -> Tuple[P.Population, Optional[NoiseModel], float]:
    """The units x odors values `U` through the curve, the noise (the flies' of `flies`, NOISE_SEED + `seed_index`; none
    where the flies have no single trials, NO_TRIALS_RULE) and the criterion (the flies' or the model's own): (hal's
    population, the noise used, the criterion)."""
    V = apply_curve(U, curve)
    nm = None
    if reading.noise != 'none':
        if flies is None:
            raise ValueError(f'{stage} {condition}: no flies for the noise')
        if has_trials(flies):
            nm = _fly_noise_cached(flies, reading.noise, [o for o in odors if o in U.columns])
            V = add_measurement_noise(V, nm, np.random.default_rng(NOISE_SEED + int(seed_index)))
    if reading.criterion == 'model_spikes':
        crit, ge = MODEL_SPIKES, True
    else:
        if flies is None:
            raise ValueError(f'{stage} {condition}: no flies for the criterion')
        crit, ge = float(flies.criterion), False
    pop = P.Population.from_criterion(V, crit, at_or_above=ge, name=name, named_units=NAMED_UNITS[stage], graded=True)
    return pop, nm, crit


def read_stage(run: SavedRun, stage: str, panel: str, odors: Sequence[str], reading: Reading,
               flies: Mapping[str, Optional['D.Population']], seed_index: int = 0,
               curve: Optional[StageCurve] = None, glomeruli: Optional[Sequence[str]] = None,
               subset_glomeruli: Optional[Sequence[str]] = None,
               lock: Optional[dict] = None) -> Tuple[Dict[str, P.Population], Dict[str, StageCurve], Dict[str, dict],
                                                     Dict[str, pd.DataFrame]]:
    """One stage of one run read on one panel, both conditions.

    Args:
        run: the saved run.
        stage, panel, odors: what is read (`run_panels`).
        reading: the stage's reading.
        flies: condition -> the flies it is read against (`fly_side`; None where there are none). The curve is set
            against the flies of `reading.curve_set_on`'s condition ('on', or each condition's own, falling back to
            APL on's where a condition has none); the noise and the criterion are each condition's flies', falling back
            to APL on's.
        seed_index: the wiring's index (pooling and noise seeds).
        curve: a curve already set (scope 'family'), used for every condition.
        glomeruli: the named stages' glomeruli (default `driver_glomeruli` of the reading's rule on the APL-on flies).
        subset_glomeruli: the bouton / claw subset's glomeruli ('shared': E90's 28 by default; 'driver': required).
        lock: the lock, for the provenance.

    Returns:
        ({'on', 'off'}: hal's `P.Population`), ({'on', 'off'}: the curve), ({'on', 'off'}: `reading_provenance`),
        ({'on', 'off'}: the values before the curve).
    """
    conds = [c for c in CONDITIONS if c in run.conditions]
    fl_on = flies.get('on')
    if reading.needs_flies and fl_on is None:
        raise ValueError(f'{run.tag} {stage} {panel}: no APL-on flies to read it against')
    if NAMED_UNITS[stage] and glomeruli is None:
        X0 = stage_changes(run, stage, 'on', odors, reading)
        glomeruli = driver_glomeruli(stage, panel, reading.glomeruli, flies=fl_on, model_glomeruli=X0.index.astype(str))
    if stage == 'pn_bouton' and reading.subset == 'driver' and subset_glomeruli is None:
        subset_glomeruli = driver_glomeruli('pn_dendrite', panel, 'every_fly',
                                            model_glomeruli=run.meta['bouton_glomerulus'])
    U = {c: sample_units(run, stage, c, odors, reading, seed_index, glomeruli, subset_glomeruli) for c in conds}
    curves: Dict[str, StageCurve] = {}
    for c in conds:
        if reading.curve == 'none':
            curves[c] = StageCurve('none')
        elif curve is not None:
            curves[c] = curve
        elif reading.curve_set_on == 'on' and c != 'on':
            curves[c] = curves['on']
        else:
            ref = (flies.get(c) or fl_on) if reading.curve_set_on == 'own' else fl_on
            curves[c] = set_curve(pooled_values([U[c]]), ref, odors, reading, seed_index)
    pops, prov = {}, {}
    name = f'{run.model}/{run.wiring}'
    for c in conds:
        fc = flies.get(c) or fl_on
        pops[c], nm, crit = finish_population(U[c], stage, c, reading, curves[c], fc, odors, seed_index, name)
        prov[c] = reading_provenance(run, stage, panel, c, reading, curves[c], nm, fc, odors,
                                     n_units_in=int(stage_changes(run, stage, c, odors, reading).shape[0]),
                                     n_units_out=int(U[c].shape[0]), seed_index=seed_index, criterion=crit,
                                     lock=lock, glomeruli=glomeruli)
    return pops, curves, prov, U


def read_run(run: SavedRun, stages: Optional[Sequence[str]] = None, readings: Optional[Mapping[str, Reading]] = None,
             seed_index: int = 0, lock: Optional[dict] = None, flies: Optional[FlyCache] = None,
             curves: Optional[Mapping[Tuple[str, str], StageCurve]] = None,
             subset_glomeruli: Optional[Sequence[str]] = None) -> ModelReading:
    """Every stage of a run (default all it saved, `SavedRun.stages`) on every panel it is read on (`run_panels`),
    both conditions, under `readings` (default `DEFAULT_READINGS`); the flies from `flies` (default a `FlyCache` of
    viki's module, which refuses before the lock). A stage and panel with no flies to read it against, where the
    reading needs them, is skipped and named in `ModelReading.skipped`. `curves` gives curves set for the family
    ((stage, panel) -> curve)."""
    readings = dict(DEFAULT_READINGS, **(readings or {}))
    flies = flies if flies is not None else FlyCache()
    out = ModelReading(model=run.model, wiring=run.wiring)
    for stage in (stages or run.stages()):
        if stage not in run.stages():
            out.skipped[(stage, '*')] = f'the run did not save {stage}'
            continue
        rd = readings[stage]
        for panel, odors in run_panels(run, stage):
            fl = {c: flies.get((stage, panel, c)) for c in CONDITIONS}
            if rd.needs_flies and fl['on'] is None:
                out.skipped[(stage, panel)] = f'no recorded {stage} on {panel} to read it against'
                continue
            pops, cv, prov, U = read_stage(run, stage, panel, odors, rd, fl, seed_index,
                                           curve=(curves or {}).get((stage, panel)), lock=lock,
                                           subset_glomeruli=subset_glomeruli)
            for c in pops:
                k = (stage, panel, c)
                out.populations[k], out.curves[k], out.provenance[k], out.before_curve[k] = pops[c], cv[c], prov[c], U[c]
    return out


def read_family(runs: Sequence[SavedRun], stages: Optional[Sequence[str]] = None,
                readings: Optional[Mapping[str, Reading]] = None, lock: Optional[dict] = None,
                flies: Optional[FlyCache] = None,
                subset_glomeruli: Optional[Sequence[str]] = None) -> Dict[str, ModelReading]:
    """`read_run` for each wiring of a family ({wiring: ModelReading}, seed index the wiring's position in `runs`);
    a reading of scope 'family' sets one curve per stage and panel on the wirings' pooled values first (APL on, against
    the APL-on flies)."""
    readings = dict(DEFAULT_READINGS, **(readings or {}))
    flies = flies if flies is not None else FlyCache()
    fam_curves: Dict[Tuple[str, str], StageCurve] = {}
    for stage in (stages or runs[0].stages()):
        rd = readings[stage]
        if rd.curve_scope != 'family' or rd.curve == 'none':
            continue
        for panel, odors in run_panels(runs[0], stage):
            fl = flies.get((stage, panel, 'on'))
            if fl is None:
                continue
            gl = None
            if NAMED_UNITS[stage]:
                gl = driver_glomeruli(stage, panel, rd.glomeruli, flies=fl,
                                      model_glomeruli=stage_changes(runs[0], stage, 'on', odors).index.astype(str))
            vals = [pooled_values([sample_units(r, stage, 'on', odors, rd, i, gl, subset_glomeruli)])
                    for i, r in enumerate(runs)]
            fam_curves[(stage, panel)] = set_curve(np.concatenate(vals), fl, odors, rd, 0)
    return {r.wiring: read_run(r, stages, readings, i, lock, flies, fam_curves, subset_glomeruli)
            for i, r in enumerate(runs)}


class SomaReading:
    """A KC-soma reading for spike counts alone (hal's experiments/E134_archive_score.py `Reading`: counts -> one
    `P.Population`), so that every archived variant is read as the family's runs are.

    Called with one condition's counts (KCs x odors): `reading(counts)` reads the APL-on condition;
    `reading(counts, condition='off', curve_counts=counts_on)` the APL-off condition with the curve set on the APL-on
    counts (a reading with `curve_set_on` 'on' refuses APL off without them, ValueError).

    Args:
        name: a 'kc_soma' reading's name (`readings_of('kc_soma')`).
        panel: the panel of the counts ('apl16' for the archive).
        seed_index: the noise seed's index.
        flies: condition -> the flies (planted ones in tests); default viki's (`fly_side`, after the lock).
    """

    def __init__(self, name: str, panel: str = 'apl16', seed_index: int = 0,
                 flies: Optional[Mapping[str, Optional['D.Population']]] = None):
        self.reading = reading_by_name('kc_soma', name)
        self.panel = panel
        self.seed_index = seed_index
        self._flies = (FlyCache({('kc_soma', panel, c): v for c, v in flies.items()}, fallback=False)
                       if flies is not None else None)

    def flies(self, condition: str) -> Optional['D.Population']:
        if self._flies is None:
            self._flies = FlyCache()
        return self._flies.get(('kc_soma', self.panel, condition))

    def __call__(self, counts: pd.DataFrame, condition: str = 'on',
                 curve_counts: Optional[pd.DataFrame] = None) -> P.Population:
        r = self.reading
        X = _abbrev(counts.astype(float))
        name = f'{self.panel}/{r.name}'
        if not r.needs_flies:
            return P.Population.from_criterion(X, MODEL_SPIKES, at_or_above=True, name=name, graded=True)
        odors = list(X.columns)
        fl_on = self.flies('on')
        if fl_on is None:
            raise ValueError(f'kc_soma {self.panel}: no APL-on flies')
        fc = self.flies(condition) or fl_on
        if r.curve == 'none':
            curve = StageCurve('none')
        elif r.curve_set_on == 'on' and condition != 'on':
            if curve_counts is None:
                raise ValueError(f"reading {r.name!r}: APL off needs the APL-on counts (curve_counts) for its curve")
            curve = set_curve(pooled_values([_abbrev(curve_counts.astype(float))[odors]]), fl_on, odors, r,
                              self.seed_index)
        else:
            curve = set_curve(pooled_values([X]), fc if r.curve_set_on == 'own' else fl_on, odors, r, self.seed_index)
        return finish_population(X, 'kc_soma', condition, r, curve, fc, odors, self.seed_index, name)[0]


def soma_readings() -> Dict[str, Callable[..., P.Population]]:
    """{name: `SomaReading(name)`} for every KC-soma reading, the readings hal's archive scoring may name."""
    return {r.name: SomaReading(r.name) for r in readings_of('kc_soma')}


def reading_provenance(run: SavedRun, stage: str, panel: str, condition: str, reading: Reading,
                       curve: Optional[StageCurve], noise: Optional[NoiseModel],
                       flies: Optional['D.Population'], odors: Sequence[str], n_units_in: int, n_units_out: int,
                       seed_index: int, criterion: float, lock: Optional[dict] = None,
                       glomeruli: Optional[Sequence[str]] = None) -> dict:
    """What a population read records: the run (source, tag, model, wiring, input panel, files with sha256, whether
    settled), the stage, panel, odors and condition, the reading (`asdict`), the units before and after sampling or
    pooling and the glomeruli of a named stage, the seeds, the curve (`StageCurve.as_dict`), the noise
    (`NoiseModel.as_dict`), the flies (cohort, flies, criterion, files), the criterion used, the lock's sha256 (its
    'sha256' key, if any)."""
    return dict(
        run=dict(source=run.source, tag=run.tag, model=run.model, wiring=run.wiring, input_panel=run.input_panel,
                 settled=run.settled, files=list(run.files)),
        stage=stage, panel=panel, condition=condition, odors=list(odors), reading=asdict(reading),
        n_units_in=int(n_units_in), n_units_out=int(n_units_out),
        glomeruli=None if glomeruli is None else [str(g) for g in glomeruli],
        pool_seed=POOL_SEED + int(seed_index), noise_seed=NOISE_SEED + int(seed_index),
        curve=None if curve is None else curve.as_dict(), noise=None if noise is None else noise.as_dict(),
        noise_note=(NO_TRIALS_RULE if reading.noise != 'none' and noise is None else None),
        implementation=implementation_constants(),
        flies=None if flies is None else dict(cohort=flies.cohort, flies=list(flies.flies),
                                              criterion=float(flies.criterion), files=list(flies.files)),
        criterion=float(criterion), lock_sha256=(lock or {}).get('sha256') if isinstance(lock, Mapping) else None)
