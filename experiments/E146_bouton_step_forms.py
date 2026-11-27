#!/usr/bin/env python3
"""E146, what could make the flies' step from the PN dendrites to the PN boutons: one weight per glomerulus, a term
that follows the whole odor, or what two sets of recordings differ by (no simulation).

The entry is item E146 of experiment_priorities.md (proposed by viki 2026-11-26 20:05 from E142's outcome (e) and next
steps, E143's consistency at the PN stages, E134's megamat somata against the control flies and the Summary's pending
list; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-26 19:30; moved in under the idle-host
rule; claimed by atlas). The step, the strength view, the no-step yardstick and the halves are those of the step from the
PN dendrites to the PN boutons by itself (E142: experiments/E142_dendrite_bouton_step.py, its step 0 in
experiments/E142_step0.py, its lock experiments/E142_lock.json): strength by odor has two parts, each odor's response
rate (part 'rate') and its mean response (part 'mean'), each a multiple of the mean over the 17 megamat odors; a part's
error is the root mean square over odors of a prediction's entry minus the flies' mean entry in the flies' one-fly unit,
and the strength view's error the root mean square of its two parts' errors (E136's adopted measure). The flies'
populations, curves, noise and criteria are E134's (its lock experiments/E134_lock.json, kcq/population_fit.py,
kcq/stage_readings.py, data_prep/stage_population_data.py).

Nothing is run, fitted to a model or searched in a model: the forms below are fitted to the recorded flies only, as a
bound over classes of forms, and enter no reading and no model.

The guard: no statistic is computed before this experiment's lock exists, is committed with no uncommitted change and
equals its copy in the shared results (`require_lock`); E142's lock (its committed copy equal to its shared copy, its
data files' sha256 with the paths read on this host, `remap_path`) and E134's lock (`stage_readings.require_lock`) are
required too. Every part but 'lock' calls `require_all` first.

Parts (CLI `--part`, in the entry's order; a later part reads the earlier parts' files):
- lock: lock.json and its committed copy experiments/E146_lock.json (no statistic is read).
- identities (design item 1): E142's no-step yardstick (all flies / fit half / held-out half), its hemibrain
  count-weighted yardstick and E134's two values between Remy's megamat somata and the control flies are reproduced;
  any that fails stops the run (identities.json says which).
- power (item 2, before any fitted weight's prediction or error is read): plants of known weights on the dendrite
  flies (the hemibrain's bouton counts; deals of each bouton fly's tuning-group sizes, E128's; draws as uneven as the
  fitted weights, whose unevenness is the only quantity of a fit to the real boutons read here), each with the noise of
  a resampled set of the bouton flies; the left-out reading against the planted no-step at every penalty of the grid;
  the penalty the lock's rule picks; "the bound has power" (power.csv, power_summary.csv, power.json).
- bound (item 3): the weight bound on left-out odors at the chosen penalty, against the no-step and the permuted
  weights; beside, the weights fitted on all odors to the fit half, read against the held-out half, and the weights
  fitted to all odors in sample (bound.csv, bound_null.csv, bound_weights.csv, bound.json).
- features (item 4): the four population features of an odor, each on left-out odors, against the no-step, the
  dendrites' entry rescaled alone and the permuted best (features.csv, features_odors.csv, features_null.csv,
  features.json).
- sets (item 5): two sets at one stage: Remy's megamat somata and the control flies on the 12 shared odors, each set's
  mean view in the other's one-fly units, and the no-step yardstick, with intervals over resamples of both sets; the
  12 odors' differences between the sets beside the step's entries (sets.csv, sets_draws.csv, sets_odors.csv,
  sets.json).
- beside (item 6, classing nothing): the method on the model's own stages (the family held and A1); the flies' step
  from the boutons to the claws under items 3 and 4; each bouton fly's own step for the seven odors; the fitted weights'
  rank correlation with the connectomes' bouton counts (beside_*.csv, beside.json).
- classes: the classes and the outcome letter (classes.json, tables.md).

Outputs (<shared>/results/E146_bouton_step_forms/): lock.json and the parts' files above, provenance_<part>.json each
with `H.code_commits()`; large intermediates (the plants' fitted weights) under outputs/E146_bouton_step_forms/.

Holds (CLAUDE.md, the entry's Limits): no bouton region is given an identity (the weights are over the dendrite
recordings' 29 named glomeruli, the same for every region, every odor and every dendrite fly); nothing is fitted to a
held-out panel (the control flies' somata are read as a yardstick only); only Pratyush's hemibrain and male CNS for the
counts (E138's `all_counts`); Remy's processed data without the bouton fly 2022-10-21 fly02 and without her older
megamat flies; Yang's processed products only; no model run, no olfsysm change. Concentrations: megamat at 10^-3 on both
sides; the 12 odors shared by megamat and the APL16 panel are at the same concentrations in both sets (paa, euc and
-aPine are not among them, so no substitution enters).

Usage (from the repository root; light, no simulation):
    ./venv/bin/python -u experiments/E146_bouton_step_forms.py --part lock     # then commit experiments/E146_lock.json
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E146_bouton_step_forms.py \\
        --part identities power bound features sets beside classes
"""
from __future__ import annotations

import argparse
import json
import os
import socket
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E146_bouton_step_forms'
SCRIPT = 'experiments/E146_bouton_step_forms.py'
LOCK_NAME = 'lock.json'
LOCK_COPY = H.HERE / 'experiments' / 'E146_lock.json'
E142_EXPERIMENT = 'E142_dendrite_bouton_step'
E134_EXPERIMENT = 'E134_population_baseline'
E128_EXPERIMENT = 'E128_bouton_groups'
E119_EXPERIMENT = 'E119_stages_per_roi'
E138_EXPERIMENT = 'E138_strength_by_odor_entries'

PANEL = 'megamat'
N_ODORS = 17
STRENGTH_PARTS: Tuple[str, ...] = ('strength/rate', 'strength/mean')
PART_KEY = {'strength/rate': 'rate', 'strength/mean': 'mean'}
OTHER_PART = {'strength/rate': 'strength/mean', 'strength/mean': 'strength/rate'}
HALVES: Tuple[str, ...] = ('fit', 'held')
FLY_SETS: Tuple[str, ...] = ('all',) + HALVES

# the classes' margins (E142's): below the no-step yardstick by MARGIN against each half of the bouton flies, and
# below the NULL_PERCENTILE-th percentile of the permuted null against all five
MARGIN = 0.10
NULL_PERCENTILE = 5.0
N_PERM = 1000
PERM_SEED = 0

# the weight per glomerulus: one non-negative weight per glomerulus of the dendrite recordings (E142's lock's 29),
# fitted by L-BFGS-B from equal weights with the analytic gradient; the loss the strength view's squared error over the
# fit odors in the bouton flies' one-fly units plus PENALTY x the weights' squared coefficient of variation
PENALTIES: Tuple[float, ...] = (0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0)
W_MIN = 1e-9
FIT_OPTIONS = dict(maxiter=3000, ftol=1e-13, gtol=1e-9)

# the power check
N_HEMIBRAIN_PLANTS = 20          # the hemibrain's bouton counts, each with its own noise draw
N_DEALS_PER_FLY = 20             # E128's deals of each bouton fly's tuning-group sizes
N_UNEVEN = 20                    # draws as uneven as the fitted weights
PLANT_KINDS: Tuple[str, ...] = ('hemibrain', 'deals', 'uneven')
POWER_ELIGIBLE = 0.9             # a plant counts when its planted no-step is 0.9 or more
POWER_SHARE = 0.9                # "9 of 10"
POWER_MIN_ELIGIBLE = 10          # eligible plants over the kinds, at least
POWER_SEED = 0
GROUP_CRITERION = 2.0            # E128's / E119's A4 tuning groups
GROUP_CUT = 0.3
HEMIBRAIN = 'hemibrain-first'
CALYCES: Tuple[str, ...] = ('hemibrain-first', 'malecns-left', 'malecns-right')

# the population features of an odor
FEATURES: Tuple[str, ...] = ('dendrite_top3_share', 'orn_summed', 'orn_top3_share', 'other_part')
DENDRITE_FEATURES: Tuple[str, ...] = ('dendrite_top3_share', 'other_part')
TOP_K = 3

# two sets at one stage
N_BOOT_SETS = 1000
SETS_SEED = 0
SETS_Q: Tuple[float, float] = (5.0, 95.0)
MEGAMAT12_CRITERION = 'lock'     # E134's megamat somata at the lock's criterion (megamat's own level)

# the identities (design item 1): quoted value, digits it is quoted to, and where the reproduction is checked
IDENTITIES = dict(
    no_step_all=dict(quoted=1.052, digits=3, source='E142 no_step_entries.csv / yardsticks.csv (all flies)'),
    no_step_fit=dict(quoted=1.184, digits=3, source='E142 no_step_entries.csv (fit half)'),
    no_step_held=dict(quoted=1.084, digits=3, source='E142 no_step_entries.csv (held-out half)'),
    hemibrain_counts=dict(quoted=0.997, digits=3, source='E142 yardsticks.csv (pn_bouton, dendrites, hemibrain-first, '
                                                         'all flies)'),
    megamat_in_control_units=dict(quoted=0.85, digits=2, source="E134 flies_megamat_vs_control.csv (lock criterion, "
                                                                 "reference control, view_strength, mean_views)"),
    control_in_megamat_units=dict(quoted=0.92, digits=2, source="E134 flies_megamat_vs_control.csv (lock criterion, "
                                                                 "reference megamat, view_strength, mean_views)"),
    control_floor=dict(quoted=0.26, digits=2, source='E134 flies_megamat_vs_control.csv ref_floor (reference control)'),
    megamat_floor=dict(quoted=0.33, digits=2, source='E134 flies_megamat_vs_control.csv ref_floor (reference megamat)'),
)
IDENTITY_TOL = 1e-9              # against the saved table's value

# the seven odors whose step is outside zero on both halves in a part (E142's flies_steps_summary.csv)
SEVEN_ODORS: Tuple[str, ...] = ('IaA', 'pa', 'ep', 't2h', '1-6ol', 'benz', 'va')
# the model's own stages (item 6): E142's groups
BESIDE_GROUPS: Tuple[str, ...] = ('family_held', 'A1')

CLASSES = {
    'power': dict(
        name='the bound has power',
        rule='on at least POWER_SHARE of the plants of each kind whose planted no-step is POWER_ELIGIBLE or more, the '
             'left-out reading (against all five synthetic bouton flies) is below the planted no-step by MARGIN; a kind '
             'with no such plant is vacuous and named; the class holds when every kind that is not vacuous passes and '
             'POWER_MIN_ELIGIBLE or more plants are eligible over the kinds; read at the chosen penalty; otherwise the '
             'weight class is "not readable"',
        share=POWER_SHARE, eligible=POWER_ELIGIBLE, margin=MARGIN, min_eligible=POWER_MIN_ELIGIBLE,
        kinds=list(PLANT_KINDS)),
    'weight_makes_step': dict(
        name='a weight per glomerulus makes the step',
        rule='at the chosen penalty, the left-out reading (17 left-out predictions, weights fitted to all five bouton '
             'flies on the other 16 odors) is below the no-step by MARGIN against each half of the bouton flies '
             "(E136's halves, each half's mean in all five flies' one-fly units), and below the NULL_PERCENTILE-th "
             'percentile of the permuted weights against all five',
        margin=MARGIN, null_percentile=NULL_PERCENTILE, n_perm=N_PERM, perm_seed=PERM_SEED),
    'feature_makes_step': dict(
        name='a population feature makes the step',
        rule='for one of the four features, the left-out reading is below the no-step by MARGIN against each half and '
             'below the NULL_PERCENTILE-th percentile of the best of the four with the features permuted among the '
             'odors, against all five',
        margin=MARGIN, null_percentile=NULL_PERCENTILE, n_perm=N_PERM, perm_seed=PERM_SEED, features=list(FEATURES)),
    'not_told_apart': dict(
        name='not told apart from two sets of recordings',
        rule="the no-step's interval (SETS_Q over N_BOOT_SETS resamples of the dendrite and bouton flies) and either "
             'interval between the sets (megamat somata in the control flies\' units, the control flies in the megamat '
             "somata's) overlap; a statement of size, naming no cause",
        q=list(SETS_Q), n_boot=N_BOOT_SETS, seed=SETS_SEED),
}
OUTCOMES = {
    'a': 'a weight per glomerulus makes the step (with power): the step lies in which glomeruli the boutons\' regions '
         'weigh (by counts, by the regions a recording holds or by a gain; the three cannot be told apart without '
         'labels); whether a weighting fitted on some odors and tested on others may enter the boutons\' reading goes '
         'to Tom and Yang',
    'b': 'a population feature makes it: the step is a term at the boutons that follows the whole odor; a planner '
         'writes it as a form, and the ORN->PN step with terms per glomerulus is worth reading again',
    'c': 'neither, with power, and the step is not told apart from two sets of recordings: the lead closes; Tom and '
         'Yang decide the boutons\' weight in the megamat stage error and the yardstick for a stage recorded by another '
         'hand',
    'd': 'neither, with power, and the step is larger than two sets differ by: the lead stays, named as lying outside '
         'the recorded glomeruli (the input glomeruli the dendrite recordings do not hold, or the regions)',
    'e': 'not readable: the bound has no power on the plants (and no feature makes the step); the decision of (c) goes '
         'to Tom and Yang without it',
}
OUTCOME_RULE = ("(a) when the power holds and the weight class holds; (b) when a feature makes the step (whatever the "
                "power; both letters when both hold); otherwise, with power: (c) when not told apart, else (d); without "
                "power and no feature: (e)")

CHOICES: Tuple[str, ...] = (
    "paths: E142's lock was written on viki, so its data files' paths (under /home/yang/src/kc_queue/"
    "model_experiment_queue and /local/tensor/Tom/claude_kc_remapping) are read at this repository's root and this "
    "host's shared directory (`remap_path`), and their sha256 checked there; E146's lock records each input both "
    "absolute (this host) and relative ('repo:' or 'shared:'), and is checked by the relative path",
    "the left-out fit: weights fitted on 16 odors with both sides relative to their mean over those 16 odors (each "
    "dendrite fly's weighted parts and each bouton fly's parts re-normalized over the 16, then averaged over the "
    "flies), so the left-out odor's bouton value enters nothing in the fit; the prediction of the 17th is the weighted "
    "dendrites' entry relative over all 17 odors (as E142's count-weighted yardstick reads it, and equal to the "
    "16-odor prediction p16 converted exactly by p17 = 17 p16 / (16 + p16))",
    "the fit's loss: the mean over the two strength parts of the mean squared difference over the fit odors in the "
    "bouton flies' one-fly units (all five flies' units, also for a half's fit), plus the penalty times the weights' "
    "squared coefficient of variation over the 29 glomeruli (population variance over the squared mean; 0 at equal "
    "weights; scale-free, as the parts are); weights bounded below by W_MIN, from equal weights, L-BFGS-B with the "
    "analytic gradient, then scaled to mean 1",
    "the penalty grid PENALTIES; the penalty is chosen by the power check alone: the one whose smallest passing share "
    "over the kinds that are not vacuous is largest, ties to the larger penalty; the class of the power is read at it",
    "the uneven draws: their unevenness is, at each penalty, the coefficient of variation over the 29 glomeruli of the "
    "weights fitted to all 17 odors of all five bouton flies at that penalty, the only quantity of a fit to the real "
    "boutons read before the power check (no prediction or error of a fit is read); each draw lognormal, iid over "
    "glomeruli, its log-SD solved so that the draw's coefficient of variation equals that one exactly",
    "a plant: the truth the planted weights' weighted strength parts of the seven dendrite flies (their mean over the "
    "flies, relative over 17 odors); its noise five bouton flies drawn with replacement, each synthetic fly the truth "
    "plus that bouton fly's deviation from the bouton flies' mean (per part and odor, uncorrected); the synthetic "
    "reference in the real bouton flies' one-fly units; the planted no-step the seven dendrite flies unweighted against "
    "the synthetic flies; the left-out reading against all five synthetic flies. The hemibrain kind: the hemibrain's "
    "bouton counts on the 29 glomeruli (E138's all_counts, 'hemibrain-first'), 20 noise draws. The deals: per bouton "
    "fly 20 deals of its tuning groups (E128's fly_groups on its trial means over E119's 12 odors, criterion 2, cut "
    "0.3; the groups only, the remainder not dealt), each group to one of the 29 glomeruli at random with replacement, "
    "a glomerulus's weight the summed sizes dealt to it (E128's deal_population's weights), each deal its own noise "
    "draw: 100 plants. The uneven kind: 20 draws, each its own noise draw. Seeds: np.random.default_rng([POWER_SEED, "
    "kind index, plant index]) for the weights and the noise; the uneven draws' standard normals are shared by the "
    "penalties (only their spread follows the penalty's unevenness)",
    "beside the power check (classing nothing): at the chosen penalty, the same plants with the truth from a resampled "
    "set of the seven dendrite flies (drawn with replacement) in place of the measured set, so that the dendrites' own "
    "noise enters too",
    "'9 of 10' read as a share of POWER_SHARE or more of each kind's eligible plants; a kind with no eligible plant is "
    "vacuous; at least POWER_MIN_ELIGIBLE eligible plants over the kinds",
    "the weight bound's permuted control: the 17 left-out folds' fitted weights each permuted among the 29 glomeruli by "
    "the same permutation (N_PERM, PERM_SEED; E142_step0.permutations), the 17 predictions scored; the class reads the "
    "5th percentile against all five bouton flies, the halves' beside",
    "the left-out predictions are fitted to all five bouton flies' entries and scored against all five and against "
    "each half (E136's halves, ctx.reference: the half's mean in all five flies' one-fly units)",
    "the features: on the positive part of the responses (values clipped at 0); 'dendrite_top3_share' the share of an "
    "odor's summed dendrite response in its three strongest glomeruli, from the per-glomerulus mean over a half of the "
    "dendrite flies (the glomeruli that half recorded); 'orn_summed' the odor's summed response at the ORN terminals, "
    "the sum over glomeruli of the per-glomerulus mean over the nine megamat ORN flies (the glomeruli any of them "
    "recorded); 'orn_top3_share' the same share at the ORN terminals; 'other_part' the odor's entry in the other "
    "strength part (for the rate part the mean part's entry, and the reverse) of the upstream flies",
    "a feature's prediction, per part: the boutons' entry (each bouton fly relative over the 16 fit odors, averaged) as "
    "a + b x the upstream entry + c x the feature, least squares on the 16 odors, predicted on the 17th and converted to "
    "the 17-odor scale (p17 = 17 p16 / (16 + p16)); the ORN features with the upstream entry of all seven dendrite "
    "flies; the dendrite features both ways: the upstream entry from one of E136's dendrite halves (3 / 4 flies) and the "
    "feature from the other, the two ways' predictions averaged; 'the dendrites' entry rescaled alone' a + b x the "
    "entry of all seven flies (and, beside, both ways on the halves)",
    "the features' permuted control: N_PERM permutations of the 17 odors (PERM_SEED), each permuting every feature's "
    "values among the odors (the same permutation for the four, and for both ways), the lowest of the four left-out "
    "readings against all five per permutation; the class reads its 5th percentile",
    "two sets at one stage: each set's flies resampled with replacement, independently (N_BOOT_SETS, SETS_SEED), each "
    "resample's scored set's mean strength parts against the reference set's resampled mean in the measured reference "
    "set's one-fly units (as hal's floor keeps them); the megamat somata at E134's lock criterion (megamat's own "
    "level), the control flies at theirs, the 12 shared odors (E134's MEGAMAT12, entries relative over the 12); the "
    "no-step the same over the seven dendrite and five bouton flies, in all five bouton flies' units; the intervals the "
    "SETS_Q percentiles of the resampled errors; 'overlap' the no-step's interval meeting either between-sets interval",
    "the 12 odors' differences between the sets: per odor and part, Remy's megamat somata's mean entry minus the control "
    "flies' (relative over the 12), in the control flies' one-fly units and in the megamat somata's, beside the flies' "
    "step by odor from the dendrites to the boutons (E142's flies_steps_summary.csv, all flies) in the bouton flies' "
    "units; their root mean squares beside",
    "beside, the model's own stages: E142's groups 'family_held' (the family held at the megamat flies' 0.0915, eight "
    "wirings) and 'A1' (E136's A1 on its four T wirings), read as E142's step 0 reads them against all flies (E134's "
    "default readings, APL on); the wirings as the dendrite flies (their dendrite populations, named glomeruli) and the "
    "wirings' bouton populations as the bouton flies (their mean entries the target, the real bouton flies' one-fly "
    "units); the model's no-step and the left-out weight reading at the chosen penalty over the glomeruli the model's "
    "dendrite reading names",
    "beside, the flies' step from the boutons to the claws: item 3 with the dendrite glomeruli weighted to the four claw "
    "flies' entries (left-out, the chosen penalty, the claw flies' units; E142's claw halves 2 / 2), against the no-step "
    "from the dendrites and from the boutons; item 4 with the boutons as the upstream stage (all five bouton flies' "
    "entry; 'other_part' both ways on E136's bouton halves 2 / 3; the dendrite top-3 share from all seven dendrite "
    "flies; the ORN features as for the step to the boutons); no class",
    "beside, each bouton fly's own step for the seven odors: the fly's entry minus the seven dendrite flies' mean entry, "
    "both parts; the seven are those outside zero on both halves in a part (E142's flies_steps_summary.csv), checked "
    "against SEVEN_ODORS",
    "beside, the fitted weights' Spearman rank correlation over the 29 glomeruli with each calyx's bouton counts (the "
    "weights fitted to all 17 odors of all five bouton flies at the chosen penalty; the median over the 17 left-out "
    "folds beside)",
    "beside, classing nothing: the weight bound at every penalty of the grid (bound_by_penalty.csv), computed after the "
    "class is read at the chosen one",
)


# ---------------------------------------------------------------------------
# small helpers

def out_dir() -> Path:
    """<shared>/results/E146_bouton_step_forms/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir() -> Path:
    """outputs/E146_bouton_step_forms/ (local, created)."""
    return H.outputs_dir(EXPERIMENT)


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _e140():
    from experiments import E140_aa_claw_soma_step as E140
    return E140


def _write_csv(df: pd.DataFrame, name: str, out: Optional[Path] = None) -> Path:
    p = (out_dir() if out is None else Path(out)) / name
    tmp = p.with_name(f'.{p.name}.{_host()}.{os.getpid()}.tmp')
    df.to_csv(tmp, index=False)
    os.replace(tmp, p)
    return p


def _write_json(name: str, obj, out: Optional[Path] = None) -> Path:
    p = (out_dir() if out is None else Path(out)) / name
    _e140()._write_json(p, obj)
    return p


def remap_path(p: str) -> Path:
    """A path recorded on another host, read on this one: '<...>/model_experiment_queue/<rest>' at this repository's
    root, '<...>/Tom/claude_kc_remapping/<rest>' in this host's shared directory, 'repo:<rest>' and 'shared:<rest>'
    likewise; any other path unchanged."""
    s = str(p)
    if s.startswith('repo:'):
        return H.HERE / s[len('repo:'):]
    if s.startswith('shared:'):
        return H.shared_dir() / s[len('shared:'):]
    for key, root in (('/model_experiment_queue/', lambda: H.HERE), ('/Tom/claude_kc_remapping/', H.shared_dir)):
        i = s.find(key)
        if i >= 0:
            return root() / s[i + len(key):]
    return Path(s)


def file_entry(rel: str) -> dict:
    """{'rel', 'path', 'bytes', 'mtime', 'sha256'} of an input given as 'repo:<path>' or 'shared:<path>'."""
    p = remap_path(rel)
    return dict(rel=rel, **_e140().file_entry(p))


def check_entry(entry: Mapping) -> Path:
    """The lock's input `entry`, read by its relative path, after checking its sha256 (ValueError if it changed)."""
    p = remap_path(entry.get('rel') or entry['path'])
    got = _e140()._sha256(p)
    if got != entry['sha256']:
        raise ValueError(f"{p}: sha256 {got[:12]} differs from the lock's {entry['sha256'][:12]}")
    return p


# ---------------------------------------------------------------------------
# the lock

def lock_inputs() -> dict:
    """The inputs the lock records with sha256: 'data', checked when read (a change stops the run): the locks this
    experiment rests on and the files it reads values from; 'results', recorded (the saved tables the identities and
    the beside items compare with; a later change is reported, not refused)."""
    data = {k: file_entry(v) for k, v in dict(
        e142_lock='repo:experiments/E142_lock.json',
        e142_lock_shared=f'shared:results/{E142_EXPERIMENT}/lock.json',
        e134_lock='repo:experiments/E134_lock.json',
        e138_classes_lock='repo:experiments/E138_classes_lock.json',
        e119_a4_groups=f'shared:results/{E119_EXPERIMENT}/a4_groups.csv',
        e119_a4_summary=f'shared:results/{E119_EXPERIMENT}/a4_summary.csv').items()}
    res = {}
    for k, v in dict(
            e142_no_step_entries=f'shared:results/{E142_EXPERIMENT}/no_step_entries.csv',
            e142_yardsticks=f'shared:results/{E142_EXPERIMENT}/yardsticks.csv',
            e142_counts=f'shared:results/{E142_EXPERIMENT}/counts.csv',
            e142_flies_steps_summary=f'shared:results/{E142_EXPERIMENT}/flies_steps_summary.csv',
            e142_models_step=f'shared:results/{E142_EXPERIMENT}/models_step.csv',
            e142_models_entries=f'shared:results/{E142_EXPERIMENT}/models_entries.csv',
            e142_classes_step0=f'shared:results/{E142_EXPERIMENT}/classes_step0.json',
            e134_flies_megamat_vs_control=f'shared:results/{E134_EXPERIMENT}/flies_megamat_vs_control.csv',
            e134_flies_no_step=f'shared:results/{E134_EXPERIMENT}/flies_no_step.csv',
            e128_part0=f'shared:results/{E128_EXPERIMENT}/part0.json',
            e138_feeders=f'shared:results/{E138_EXPERIMENT}/feeders.csv').items():
        p = remap_path(v)
        res[k] = file_entry(v) if p.exists() else dict(rel=v, path=str(p), missing=True)
    return dict(data=data, results=res)


def recordings() -> dict:
    """The recordings read, names only (no value is read): the dendrite, bouton and claw flies with their halves (E142's
    lock), the megamat ORN flies (the features), Remy's megamat KC somata and the control flies (two sets at one
    stage); the exclusions."""
    from data_prep import stage_population_data as D
    l142 = json.loads((H.HERE / 'experiments' / 'E142_lock.json').read_text())
    out = {st: dict(flies=list(l142['recordings'][st]['flies']), halves=l142['recordings'][st]['halves'],
                    criterion=l142['recordings'][st]['criterion'], units=l142['recordings'][st]['units'])
           for st in ('pn_dendrite', 'pn_bouton', 'kc_claw')}
    for name, key in (('orn', ('orn', PANEL, 'wild_type')), ('kc_soma_megamat', ('kc_soma', PANEL, 'wild_type')),
                      ('kc_soma_control', ('kc_soma', 'apl16', 'control'))):
        p = D.load_population(*key)
        out[name] = dict(key=list(key), flies=list(p.flies), criterion=float(p.criterion), units=p.units)
    if any(f.startswith('2022-10-21__fly02') for f in out['pn_bouton']['flies']):
        raise RuntimeError('the excluded bouton fly is in the bouton population')
    out['excluded'] = dict(
        flies=[dict(stage=x['stage'], fly=x['fly'], reason=x['reason']) for x in D.EXCLUDED_FLIES],
        older_megamat="Remy's older megamat flies (pyolfsysm's data/from_remy/2024-11-12): never read",
        validation2='not read', kiwi_control_mixtures='not read', tnt='not read (APL output blocked flies)')
    return out


def lock_contents() -> dict:
    """The lock: the inputs (sha256), the recordings, halves and glomeruli, the odors, the definitions, seeds,
    permutation counts and penalty grid, the classes and outcomes, the identities, the choices and the holds. No
    statistic is read."""
    from experiments import E134_flies as F
    l142 = json.loads((H.HERE / 'experiments' / 'E142_lock.json').read_text())
    definitions = dict(
        weight_per_glomerulus="E142's count-weighted yardstick (E142_step0.weighted_arrays) with one non-negative weight "
                              'per glomerulus of the dendrite recordings in place of the counts, the same for every odor '
                              'and every dendrite fly: rate the weighted fraction of glomeruli responding, mean the '
                              'weighted mean response, each relative to its mean over the odors, the dendrite flies\' '
                              'mean parts the prediction',
        left_out='weights fitted on 16 odors predict the 17th odor\'s two entries; the 17 predictions make one view, '
                 'scored against all five bouton flies and against each half',
        permuted_weights='the fitted weights permuted among the glomeruli, N_PERM permutations',
        power_check='known weights planted on the dendrite flies (the hemibrain\'s bouton counts, 20 deals of each '
                    "bouton fly's tuning-group sizes, 20 draws as uneven as the fitted weights), each with the noise of "
                    'a resampled set of the bouton flies',
        features=dict(names=list(FEATURES), top_k=TOP_K, dendrite_features=list(DENDRITE_FEATURES),
                      rule="the boutons' entry as a linear function of the dendrites' entry and one feature, with the "
                           'odor left out; a dendrite feature from the half of the dendrite flies the upstream entry is '
                           'not, both ways'),
        two_sets="Remy's megamat soma flies and the control flies on the 12 shared odors (E134): each set's mean view "
                 "in the other's one-fly units, with a 5-95 % interval over N_BOOT_SETS resamples of both sets; the "
                 'no-step yardstick the same interval over the dendrite and bouton flies',
        megamat12=list(F.MEGAMAT12), megamat12_criterion=MEGAMAT12_CRITERIA_NOTE)
    holds = ['no simulation, no model run, no olfsysm change; nothing fitted to a model or searched in one',
             'no bouton region given an identity: the weights are over the dendrite recordings\' named glomeruli, the '
             'same for every region; they enter no reading and no model without Tom\'s or Yang\'s word',
             'nothing fitted to a held-out panel: the control flies\' somata are read as a yardstick only',
             "only Pratyush's hemibrain and male CNS for the counts (E138's all_counts)",
             "Remy's processed data without the bouton fly 2022-10-21 fly02 and without her older megamat flies; "
             "Yang's processed products only", 'validation2 and the mixtures not read']
    return dict(experiment=EXPERIMENT, script=SCRIPT, part='lock', written=_now(), host=_host(),
                commits=H.code_commits(), panel=PANEL, odors=list(l142['odors']),
                e142_lock_sha256=_e140()._sha256(H.HERE / 'experiments' / 'E142_lock.json'),
                inputs=lock_inputs(), recordings=recordings(),
                glomeruli=dict(dendrites=list(l142['glomeruli']['dendrites'])),
                definitions=definitions,
                constants=dict(margin=MARGIN, null_percentile=NULL_PERCENTILE, n_perm=N_PERM, perm_seed=PERM_SEED,
                               penalties=list(PENALTIES), w_min=W_MIN, fit_options=FIT_OPTIONS,
                               n_hemibrain_plants=N_HEMIBRAIN_PLANTS, n_deals_per_fly=N_DEALS_PER_FLY,
                               n_uneven=N_UNEVEN, plant_kinds=list(PLANT_KINDS), power_eligible=POWER_ELIGIBLE,
                               power_share=POWER_SHARE, power_min_eligible=POWER_MIN_ELIGIBLE, power_seed=POWER_SEED,
                               group_criterion=GROUP_CRITERION, group_cut=GROUP_CUT, hemibrain=HEMIBRAIN,
                               calyces=list(CALYCES), features=list(FEATURES), top_k=TOP_K,
                               n_boot_sets=N_BOOT_SETS, sets_seed=SETS_SEED, sets_q=list(SETS_Q),
                               megamat12_criterion=MEGAMAT12_CRITERION, seven_odors=list(SEVEN_ODORS),
                               beside_groups=list(BESIDE_GROUPS)),
                identities=IDENTITIES, identity_tol=IDENTITY_TOL, classes=CLASSES, outcomes=OUTCOMES,
                outcome_rule=OUTCOME_RULE, choices=list(CHOICES), holds=holds,
                note='written before any statistic of E146 was computed; every part but lock refuses until this file is '
                     'committed (require_lock)')


MEGAMAT12_CRITERIA_NOTE = ("the megamat somata at E134's lock criterion (MEGAMAT_CRITERIA 'lock', megamat's own level, "
                           "the reading E134's 0.85 and 0.92 are of); the control flies at theirs")


def part_lock() -> dict:
    """results/E146_bouton_step_forms/lock.json and its committed copy experiments/E146_lock.json (commit it before any
    other part); refuses to overwrite a committed lock."""
    E140 = _e140()
    if LOCK_COPY.exists() and E140._committed(LOCK_COPY):
        raise PermissionError(f'{LOCK_COPY.name} is committed already; the lock is fixed')
    lock = lock_contents()
    text = json.dumps(lock, indent=1, default=E140._json_default) + '\n'
    (out_dir() / LOCK_NAME).write_text(text)
    LOCK_COPY.write_text(text)
    print(f'wrote {out_dir() / LOCK_NAME} and {LOCK_COPY} (commit the copy before any other part)', flush=True)
    return lock


def require_lock(copy: Optional[Path] = None, shared: Optional[Path] = None) -> dict:
    """The committed lock: experiments/E146_lock.json exists, is tracked with no uncommitted change, and has the same
    bytes as results/.../lock.json. PermissionError otherwise (no statistic is read before)."""
    E140 = _e140()
    copy = LOCK_COPY if copy is None else Path(copy)
    if not copy.exists():
        raise PermissionError(f'{copy.name}: no lock yet (--part lock, then commit it)')
    if not E140._committed(copy):
        raise PermissionError(f'{copy.name} is not committed (or has uncommitted changes)')
    shared = (out_dir() / LOCK_NAME) if shared is None else Path(shared)
    if not shared.exists() or shared.read_bytes() != copy.read_bytes():
        raise PermissionError(f'{shared} differs from the committed {copy.name}')
    lock = json.loads(copy.read_text())
    lock['sha256'] = E140._sha256(copy)
    return lock


def require_all() -> Tuple[dict, dict, dict]:
    """(this experiment's lock, E142's lock, E134's lock), all required before any statistic: this lock's data files
    checked; E142's lock committed and equal to its shared copy (E142.require_lock), its data files checked at their
    paths on this host (`remap_path`); E134's lock (stage_readings.require_lock)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the context

class Ctx:
    """What every part reads: the locks, E142's step-0 Context (the flies of every set on the strength parts, all flies
    and each half; E136's Flies), the dendrite flies' matrices on the 29 glomeruli, the bouton flies' references, the
    no-step's views, the ORN flies, the bouton flies' tuning-group sizes and the connectomes' counts."""

    def __init__(self, lock: Mapping, l142: Mapping, e134: Mapping):
        raise NotImplementedError


# ---------------------------------------------------------------------------
# pure functions: the weighted entries, the fit, the left-out predictions

def weighted_entries(w: np.ndarray, mats: Sequence[Tuple[np.ndarray, np.ndarray]], norm: Sequence[int],
                     jacobian: bool = False) -> Dict[str, Tuple[np.ndarray, Optional[np.ndarray]]]:
    """The dendrite flies' mean weighted strength parts for one weighting, relative to their mean over the odors `norm`,
    with the Jacobian.

    Args:
        w: glomeruli (non-negative; at least one positive in each fly's recorded glomeruli).
        mats: per fly (values, responds), glomeruli x odors (NaN where the fly lacks a glomerulus), E142_step0's
            `fly_matrices`.
        norm: the positions of the odors each fly's parts are made relative over.
        jacobian: also return d entries / d w.

    Returns:
        {part: (entries over all odors (odors,), Jacobian (glomeruli x odors) or None)}; with `norm` all odors equal to
        E142_step0.weighted_arrays averaged over the flies.
    """
    raise NotImplementedError


def relative_parts(parts: Mapping[str, np.ndarray], norm: Sequence[int]) -> Dict[str, np.ndarray]:
    """Each fly's parts (flies x odors) relative to its mean over the odors `norm`, averaged over the flies (NaN-aware):
    the bouton flies' target for a fit on those odors."""
    raise NotImplementedError


def fit_loss(w: np.ndarray, mats, targets: Mapping[str, np.ndarray], units: Mapping[str, float], fit: Sequence[int],
             penalty: float) -> Tuple[float, np.ndarray]:
    """The fit's loss and gradient: the mean over the two parts of the mean squared difference over `fit` between the
    weighted dendrites' entries (relative over `fit`) and `targets` (relative over `fit`) in `units`, plus `penalty`
    times the weights' squared coefficient of variation."""
    raise NotImplementedError


def fit_weights(mats, targets: Mapping[str, np.ndarray], units: Mapping[str, float], fit: Sequence[int],
                penalty: float, w0: Optional[np.ndarray] = None) -> Tuple[np.ndarray, dict]:
    """The weights minimizing `fit_loss` (L-BFGS-B from equal weights, bounded below by W_MIN), scaled to mean 1, and
    the optimizer's record (loss, iterations, success)."""
    raise NotImplementedError


def left_out(mats, fly_parts: Mapping[str, np.ndarray], units: Mapping[str, float], penalty: float
             ) -> Tuple[Dict[str, np.ndarray], np.ndarray, List[dict]]:
    """The left-out predictions: for each odor, weights fitted on the other odors to the flies' parts (`fly_parts`,
    flies x odors, each fly re-normalized over the fit odors), the odor's two entries predicted (the weighted dendrites
    relative over all odors). Returns ({part: predictions (odors,)}, weights (odors x glomeruli), the fits' records)."""
    raise NotImplementedError


def predictions_of(weights: np.ndarray, mats) -> Dict[str, np.ndarray]:
    """The left-out predictions of given per-fold weights (odors x glomeruli; fold o predicts odor o), many sets at once
    when `weights` is sets x odors x glomeruli ({part: sets x odors})."""
    raise NotImplementedError


def score_entries(pred: Mapping[str, np.ndarray], ref) -> dict:
    """The strength view's error and its two parts of predictions (part -> odors, or sets x odors) against a reference
    (hal's `score_arrays`): {'view_strength', 'part_rate', 'part_mean'} (arrays for many sets)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the power check

def plant_truth(w: np.ndarray, mats) -> Dict[str, np.ndarray]:
    """The planted truth: the weighted strength parts of the dendrite flies (their mean, relative over all odors)."""
    raise NotImplementedError


def noise_fly_parts(truth: Mapping[str, np.ndarray], bouton_parts: Mapping[str, np.ndarray], rng) -> Dict[str, np.ndarray]:
    """Synthetic bouton flies: as many bouton flies drawn with replacement, each the truth plus that fly's deviation from
    the bouton flies' mean ({part: flies x odors})."""
    raise NotImplementedError


def synthetic_reference(parts: Mapping[str, np.ndarray], units: Mapping[str, float], odors: Sequence[str]):
    """A hal FlyReference of synthetic flies' parts in given one-fly units."""
    raise NotImplementedError


def deal_weights(sizes: Sequence[int], n: int, rng) -> np.ndarray:
    """One deal (E128's): each group of `sizes` to one of `n` glomeruli at random with replacement; a glomerulus's
    weight the summed sizes dealt to it."""
    raise NotImplementedError


def uneven_weights(cv: float, z: np.ndarray) -> np.ndarray:
    """A lognormal draw exp(s z) with s solved so that its coefficient of variation (population SD over mean) is `cv`
    exactly; equal weights for cv 0."""
    raise NotImplementedError


def cv_of(w: np.ndarray) -> float:
    """The coefficient of variation of weights (population SD over the mean)."""
    raise NotImplementedError


def power_rule(rows: pd.DataFrame) -> dict:
    """"The bound has power" on one penalty's plants (columns kind, no_step, left_out): per kind the eligible plants
    (no_step >= POWER_ELIGIBLE), the passing ones (left_out <= no_step - MARGIN), the share and whether the kind passes
    or is vacuous; the class."""
    raise NotImplementedError


def choose_penalty(by_penalty: Mapping[float, dict]) -> float:
    """The lock's rule: the penalty whose smallest passing share over the kinds that are not vacuous is largest, ties to
    the larger penalty."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the features

def top_share(values: pd.DataFrame, k: int = TOP_K) -> pd.Series:
    """Per odor (column), the share of its summed positive response (values clipped at 0) in its k strongest units."""
    raise NotImplementedError


def flies_mean_table(pops: Sequence, odors: Sequence[str]) -> pd.DataFrame:
    """Per unit (glomerulus) and odor, the mean over the flies that recorded it."""
    raise NotImplementedError


def linear_left_out(x: np.ndarray, f: Optional[np.ndarray], targets16: Mapping[int, np.ndarray]) -> np.ndarray:
    """Left-out predictions of a + b x (+ c f): for each odor o, least squares on the other odors against
    `targets16[o]` (the target relative over the other odors, all odors' positions), predicted at o and converted to
    the all-odor scale (p17 = n p / (n - 1 + p))."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# two sets at one stage

def resampled_errors(scored: Mapping[str, np.ndarray], reference: Mapping[str, np.ndarray],
                     units: Mapping[str, float], n_boot: int = N_BOOT_SETS, seed: int = SETS_SEED) -> np.ndarray:
    """The strength view's error of the scored set's mean parts against the reference set's mean, with both sets'
    flies resampled with replacement independently, in fixed units (n_boot,)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the classes

def weight_class(bound: Mapping, power_holds: bool) -> dict:
    """"A weight per glomerulus makes the step" from bound.json's numbers; "not readable" without power."""
    raise NotImplementedError


def feature_class(features: Mapping) -> dict:
    """"A population feature makes the step" from features.json's numbers."""
    raise NotImplementedError


def overlap(a: Sequence[float], b: Sequence[float]) -> bool:
    """Whether two intervals meet."""
    raise NotImplementedError


def outcome(power_holds: bool, weight_holds: bool, feature_holds: bool, not_told_apart: bool) -> List[str]:
    """The outcome letters by OUTCOME_RULE."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the parts

def part_identities() -> dict:
    """Design item 1 (identities.json; stops on a failure)."""
    raise NotImplementedError


def part_power() -> dict:
    """Design item 2 (power.csv, power_summary.csv, power_beside.csv, power.json)."""
    raise NotImplementedError


def part_bound() -> dict:
    """Design item 3 (bound.csv, bound_null.csv, bound_weights.csv, bound_by_penalty.csv, bound.json)."""
    raise NotImplementedError


def part_features() -> dict:
    """Design item 4 (features.csv, features_odors.csv, features_null.csv, features.json)."""
    raise NotImplementedError


def part_sets() -> dict:
    """Design item 5 (sets.csv, sets_draws.csv, sets_odors.csv, sets.json)."""
    raise NotImplementedError


def part_beside() -> dict:
    """Design item 6 (beside_models.csv, beside_claws.csv, beside_flies_seven.csv, beside_counts_rho.csv,
    beside.json)."""
    raise NotImplementedError


def part_classes() -> dict:
    """The classes and the outcome letter (classes.json, tables.md)."""
    raise NotImplementedError


PARTS: Tuple[str, ...] = ('identities', 'power', 'bound', 'features', 'sets', 'beside', 'classes')


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', nargs='+', required=True, choices=('lock',) + PARTS)
    a = ap.parse_args(argv)
    if 'lock' in a.part:
        part_lock()
        return
    fns = dict(identities=part_identities, power=part_power, bound=part_bound, features=part_features, sets=part_sets,
               beside=part_beside, classes=part_classes)
    for p in PARTS:
        if p in a.part:
            fns[p]()


if __name__ == '__main__':
    main()
