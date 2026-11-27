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
- lock: lock.json and its committed copy experiments/E146_lock.json (no statistic is read); amend: the lock's
  amendments (AMENDMENTS) written into it, refused once any statistic exists.
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

# the lock's amendments (each before any statistic of E146 was computed; `part_amend`)
AMENDMENTS: Tuple[str, ...] = (
    "2026-11-26 20:28, atlas, before any statistic of E146 was computed (only the tests on planted data had run; the "
    "results folder held lock.json alone): the left-out fit's dendrite side is the dendrite flies' mean weighted "
    "entries (each fly relative over all 17 odors, as the prediction reads them) re-normalized over the 16 fit odors, "
    "in place of each dendrite fly re-normalized over the 16 and then averaged (the lock's second choice as first "
    "written); the bouton side is unchanged (each bouton fly re-normalized over the 16, then averaged). With the form "
    "as first written a plant's own weights did not reproduce its noiseless truth on the fit odors (the mean of "
    "re-normalized flies is not the re-normalized mean), so the bound would have been charged an error that no "
    "weighting removes; with the amended form the planted weights reproduce the truth exactly. The prediction of the "
    "left-out odor is unchanged (the weighted dendrites' entry relative over all 17). Code: fit_entries, fit_loss",
)

CHOICES: Tuple[str, ...] = (
    "paths: E142's lock was written on viki, so its data files' paths (under /home/yang/src/kc_queue/"
    "model_experiment_queue and /local/tensor/Tom/claude_kc_remapping) are read at this repository's root and this "
    "host's shared directory (`remap_path`), and their sha256 checked there; E146's lock records each input both "
    "absolute (this host) and relative ('repo:' or 'shared:'), and is checked by the relative path",
    "the left-out fit (as amended 2026-11-26 20:28, AMENDMENTS): weights fitted on 16 odors with both sides relative to "
    "their mean over those 16 odors: the dendrite flies' mean weighted entries (each fly relative over all 17 odors, as "
    "the prediction reads them) re-normalized over the 16, and each bouton fly's parts re-normalized over the 16, then "
    "averaged over the flies, so the left-out odor's bouton value enters nothing in the fit; the prediction of the 17th "
    "is the weighted dendrites' entry relative over all 17 odors (as E142's count-weighted yardstick reads it, and equal "
    "to the 16-odor prediction p16 converted exactly by p17 = 17 p16 / (16 + p16))",
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
                outcome_rule=OUTCOME_RULE, choices=list(CHOICES), holds=holds, amendments=list(AMENDMENTS),
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


def part_amend() -> dict:
    """The lock's amendments (AMENDMENTS) written into the committed lock: its 'choices' replaced by CHOICES and its
    'amendments' by AMENDMENTS, every other field kept as written; refused once any part but the lock has written to
    the results folder (an amendment comes before any statistic). Commit the copy after."""
    E140 = _e140()
    if not (LOCK_COPY.exists() and E140._committed(LOCK_COPY)):
        raise PermissionError('no committed lock to amend')
    others = [p.name for p in out_dir().iterdir() if p.name != LOCK_NAME]
    if others:
        raise PermissionError(f'statistics exist already ({others[:5]}); the lock is not amended after them')
    lock = json.loads(LOCK_COPY.read_text())
    lock['choices'] = list(CHOICES)
    lock['amendments'] = list(AMENDMENTS)
    lock['amended'] = dict(time=_now(), host=_host())
    text = json.dumps(lock, indent=1, default=E140._json_default) + '\n'
    (out_dir() / LOCK_NAME).write_text(text)
    LOCK_COPY.write_text(text)
    print(f'amended {out_dir() / LOCK_NAME} and {LOCK_COPY} (commit the copy)', flush=True)
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
    from experiments import E142_dendrite_bouton_step as E142
    from kcq import stage_readings as R
    lock = require_lock()
    for e in lock['inputs']['data'].values():
        check_entry(e)
    l142 = E142.require_lock()
    if l142['sha256'] != lock['e142_lock_sha256']:
        raise PermissionError("E142's lock differs from the one E146's lock records")
    for e in l142['inputs']['data'].values():
        p = remap_path(e['path'])
        if _e140()._sha256(p) != e['sha256']:
            raise ValueError(f"{p}: differs from E142's lock")
    return lock, l142, R.require_lock()


def provenance(part: str, ctx: Optional['Ctx'] = None, lock: Optional[Mapping] = None, **extra) -> None:
    """provenance_<part>.json: the script, `H.code_commits()`, host, time, the locks' sha256 and `extra`."""
    lk = ctx.lock if ctx is not None else (lock or {})
    _write_json(f'provenance_{part}.json',
                dict(part=part, script=SCRIPT, **H.code_commits(), host=_host(), time=_now(),
                     lock_sha256=lk.get('sha256'), e142_lock_sha256=(ctx.l142.get('sha256') if ctx else None),
                     e134_lock_sha256=(ctx.e134.get('sha256') if ctx else None), **extra))


def _read_json(name: str) -> dict:
    p = out_dir() / name
    if not p.exists():
        raise FileNotFoundError(f'{name}: run its part first')
    return json.loads(p.read_text())


# ---------------------------------------------------------------------------
# the context

class Ctx:
    """What every part reads: the locks, E142's step-0 Context (the flies of every set on the strength parts, all flies
    and each half; E136's Flies), the dendrite flies' matrices on the 29 glomeruli, the bouton flies' references, the
    no-step's views, the ORN flies, the bouton flies' tuning-group sizes and the connectomes' counts."""

    def __init__(self, lock: Mapping, l142: Mapping, e134: Mapping):
        from experiments import E142_step0 as S0
        self.lock, self.l142, self.e134 = dict(lock), dict(l142), dict(e134)
        self.s0 = S0.Context(l142, e134)
        self.odors = list(lock['odors'])
        if self.s0.odors != self.odors:
            raise RuntimeError("E142's odors differ from E146's lock's")
        self.glomeruli = list(lock['glomeruli']['dendrites'])
        self.pops = self.s0.pops
        self.names = self.s0.names
        for st in ('pn_dendrite', 'pn_bouton', 'kc_claw'):
            if [p.name for p in self.pops[st]] != list(lock['recordings'][st]['flies']):
                raise RuntimeError(f"{st}: flies differ from E146's lock")
            if self.names[st] != lock['recordings'][st]['halves']:
                raise RuntimeError(f"{st}: halves differ from E146's lock")
        self.mats = S0.fly_matrices(self.pops['pn_dendrite'], self.glomeruli, self.odors)
        self.prep = prepare(self.mats)
        self.refs = {s: self.s0.reference('pn_bouton', s) for s in FLY_SETS}
        self.units = {p: float(self.refs['all'].one_fly[p]) for p in STRENGTH_PARTS}
        for s, r in self.refs.items():
            for p in STRENGTH_PARTS:
                if list(r.mean[p].index) != self.odors or float(r.one_fly[p]) != self.units[p]:
                    raise RuntimeError(f'the bouton reference {s} is not on the odors or units of all flies')
        self.bouton_parts = {p: self.refs['all'].parts[p].reindex(columns=self.odors).to_numpy(dtype=float)
                             for p in STRENGTH_PARTS}
        self.n = len(self.odors)
        self.all_idx = np.arange(self.n)
        self._cache: Dict[str, object] = {}

    def claw_refs(self) -> Dict[str, object]:
        """The claw flies' references (E136's halves; all four flies' units)."""
        return {s: self.s0.reference('kc_claw', s) for s in FLY_SETS}

    def dendrite_set(self, s: str) -> List:
        """The dendrite flies of a fly set."""
        return self.s0.set_pops('pn_dendrite', s)

    def orn_pops(self) -> List:
        """The megamat wild-type ORN flies (hal's populations), checked against the lock's names."""
        if 'orn' not in self._cache:
            from data_prep import stage_population_data as D
            from kcq import population_fit as P
            orn = P.fly_populations(D.load_population('orn', PANEL, 'wild_type'))
            if [o.name for o in orn] != list(self.lock['recordings']['orn']['flies']):
                raise RuntimeError("the ORN flies differ from E146's lock")
            self._cache['orn'] = orn
        return self._cache['orn']

    def group_sizes(self) -> Dict[str, np.ndarray]:
        """Per bouton fly, its tuning groups' sizes (E128's fly_groups on its trial means over E119's 12 odors,
        criterion 2, cut 0.3; the groups only), checked against E119's a4_summary.csv (responsive ROIs, groups)."""
        if 'sizes' not in self._cache:
            from data_prep import stage_population_data as D
            from experiments import E119_stages_per_roi as E119
            from experiments import E128_bouton_groups as E128
            bp = D.load_population('pn_bouton', PANEL, 'wild_type', trials=True)
            if list(bp.flies) != list(self.lock['recordings']['pn_bouton']['flies']):
                raise RuntimeError("the bouton flies differ from E146's lock")
            s4 = pd.read_csv(check_entry(self.lock['inputs']['data']['e119_a4_summary'])).set_index('acq')
            out = {}
            for f in bp.flies:
                M = bp.tables[f].reindex(columns=list(E119.ODORS)).to_numpy(dtype=float)
                lab = E128.fly_groups(M, GROUP_CRITERION, GROUP_CUT)
                sizes = np.unique(lab[lab > 0], return_counts=True)[1]
                if f in s4.index and (len(sizes) != int(s4.loc[f, 'n_tuning_groups'])
                                      or int((lab > 0).sum()) != int(s4.loc[f, 'n_responsive'])):
                    raise RuntimeError(f"{f}: the tuning groups differ from E119's A4")
                out[f] = sizes
            self._cache['sizes'] = out
        return self._cache['sizes']

    def counts(self) -> pd.DataFrame:
        """Each calyx's bouton counts on the 29 glomeruli (E138's all_counts; glomeruli x calyces, 0 where absent),
        checked against E142's counts.csv."""
        if 'counts' not in self._cache:
            from experiments import E138_strength_by_odor_entries as E138
            c = E138.all_counts()
            out = pd.DataFrame({cal: c[cal].n_boutons.reindex(self.glomeruli).fillna(0).astype(float)
                                for cal in CALYCES})
            e142 = pd.read_csv(remap_path(self.lock['inputs']['results']['e142_counts']['rel']))
            for cal in CALYCES:
                x = e142[e142.calyx == cal].set_index('glomerulus').n_boutons.reindex(self.glomeruli).fillna(0)
                if not np.allclose(x.to_numpy(dtype=float), out[cal].to_numpy()):
                    raise RuntimeError(f"{cal}: the counts differ from E142's counts.csv")
            self._cache['counts'] = out
        return self._cache['counts']

    def no_step_views(self) -> Dict[str, np.ndarray]:
        """The no-step's prediction: the seven dendrite flies' mean strength parts (relative over all odors)."""
        return {p: v for p, (v, _) in weighted_entries(np.ones(len(self.glomeruli)), self.prep, self.all_idx).items()}


def context() -> Ctx:
    """`require_all`, then the context."""
    return Ctx(*require_all())


# ---------------------------------------------------------------------------
# pure functions: the weighted entries, the fit, the left-out predictions

def prepare(mats: Sequence[Tuple[np.ndarray, np.ndarray]]) -> dict:
    """The dendrite flies' matrices stacked for `weighted_entries`: M (flies x glomeruli x odors, 1 where recorded) and
    per part X (the responds or the values, 0 where not recorded)."""
    if isinstance(mats, dict):
        return mats
    M = np.stack([np.isfinite(np.asarray(V, dtype=float)).astype(float) for V, _ in mats])
    Xr = np.stack([np.nan_to_num(np.asarray(Rr, dtype=float) * np.isfinite(np.asarray(V, dtype=float)))
                   for V, Rr in mats])
    Xm = np.stack([np.nan_to_num(np.asarray(V, dtype=float)) for V, _ in mats])
    return dict(M=M, X={'strength/rate': Xr, 'strength/mean': Xm})


def resample_prep(prep: Mapping, idx: Sequence[int]) -> dict:
    """The prepared matrices of a resampled set of the flies (positions `idx`)."""
    idx = np.asarray(idx, dtype=int)
    return dict(M=prep['M'][idx], X={p: X[idx] for p, X in prep['X'].items()})


def weighted_entries(w: np.ndarray, mats, norm: Sequence[int],
                     jacobian: bool = False) -> Dict[str, Tuple[np.ndarray, Optional[np.ndarray]]]:
    """The dendrite flies' mean weighted strength parts for one weighting, relative to their mean over the odors `norm`,
    with the Jacobian.

    Args:
        w: glomeruli (non-negative; at least one positive in each fly's recorded glomeruli).
        mats: per fly (values, responds), glomeruli x odors (NaN where the fly lacks a glomerulus), E142_step0's
            `fly_matrices`, or their `prepare`d stack.
        norm: the positions of the odors each fly's parts are made relative over.
        jacobian: also return d entries / d w.

    Returns:
        {part: (entries over all odors (odors,), Jacobian (glomeruli x odors) or None)}; with `norm` all odors equal to
        E142_step0.weighted_arrays averaged over the flies.
    """
    prep = prepare(mats)
    w = np.asarray(w, dtype=float)
    norm = np.asarray(norm, dtype=int)
    M = prep['M']
    D = np.einsum('g,fgo->fo', w, M)
    out = {}
    for p in STRENGTH_PARTS:
        X = prep['X'][p]
        a = np.einsum('g,fgo->fo', w, X) / D
        abar = a[:, norm].mean(axis=1)
        e = a / abar[:, None]
        E = np.nanmean(e, axis=0) if not jacobian else e.mean(axis=0)
        J = None
        if jacobian:
            Ja = (X - M * a[:, None, :]) / D[:, None, :]
            Jb = Ja[:, :, norm].mean(axis=2)
            J = ((Ja - e[:, None, :] * Jb[:, :, None]) / abar[:, None, None]).mean(axis=0)
        out[p] = (E, J)
    return out


def relative_parts(parts: Mapping[str, np.ndarray], norm: Sequence[int]) -> Dict[str, np.ndarray]:
    """Each fly's parts (flies x odors) relative to its mean over the odors `norm`, averaged over the flies (NaN-aware):
    the bouton flies' target for a fit on those odors."""
    norm = np.asarray(norm, dtype=int)
    out = {}
    for p, A in parts.items():
        A = np.asarray(A, dtype=float)
        m = np.nanmean(A[:, norm], axis=1, keepdims=True)
        out[p] = np.nanmean(A / m, axis=0)
    return out


def fit_entries(w: np.ndarray, mats, fit: Sequence[int], jacobian: bool = False
                ) -> Dict[str, Tuple[np.ndarray, Optional[np.ndarray]]]:
    """The fit's side of the dendrites (the lock's amendment of 2026-11-26 20:28): the dendrite flies' mean weighted
    entries relative over all odors (`weighted_entries`, as the prediction reads them), re-normalized over the odors
    `fit`; with the Jacobian. Returns {part: (entries over all odors, Jacobian (glomeruli x odors) or None)}."""
    prep = prepare(mats)
    fit = np.asarray(fit, dtype=int)
    ent = weighted_entries(w, prep, np.arange(prep['M'].shape[2]), jacobian=jacobian)
    out = {}
    for p, (E, J) in ent.items():
        m = float(E[fit].mean())
        e = E / m
        Jf = None
        if jacobian:
            Jf = (J - e[None, :] * J[:, fit].mean(axis=1)[:, None]) / m
        out[p] = (e, Jf)
    return out


def fit_loss(w: np.ndarray, mats, targets: Mapping[str, np.ndarray], units: Mapping[str, float], fit: Sequence[int],
             penalty: float) -> Tuple[float, np.ndarray]:
    """The fit's loss and gradient: the mean over the two parts of the mean squared difference over `fit` between the
    weighted dendrites' entries (`fit_entries`: relative over all odors, re-normalized over `fit`) and `targets`
    (relative over `fit`) in `units`, plus `penalty` times the weights' squared coefficient of variation."""
    w = np.asarray(w, dtype=float)
    fit = np.asarray(fit, dtype=int)
    ent = fit_entries(w, mats, fit, jacobian=True)
    loss, grad = 0.0, np.zeros_like(w)
    k = len(STRENGTH_PARTS)
    for p in STRENGTH_PARTS:
        E, J = ent[p]
        r = (E[fit] - np.asarray(targets[p], dtype=float)[fit]) / units[p]
        loss += float(np.mean(r ** 2)) / k
        grad += (2.0 / len(fit)) * (J[:, fit] @ r) / units[p] / k
    if penalty:
        G = len(w)
        m, q = float(w.mean()), float(np.mean(w ** 2))
        loss += penalty * (q / m ** 2 - 1.0)
        grad += penalty * (2.0 / (G * m ** 2)) * (w - q / m)
    return loss, grad


def fit_weights(mats, targets: Mapping[str, np.ndarray], units: Mapping[str, float], fit: Sequence[int],
                penalty: float, w0: Optional[np.ndarray] = None) -> Tuple[np.ndarray, dict]:
    """The weights minimizing `fit_loss` (L-BFGS-B from equal weights, bounded below by W_MIN), scaled to mean 1, and
    the optimizer's record (loss, iterations, success)."""
    from scipy.optimize import minimize
    prep = prepare(mats)
    G = prep['M'].shape[1]
    x0 = np.ones(G) if w0 is None else np.asarray(w0, dtype=float)
    res = minimize(fit_loss, x0, args=(prep, targets, units, fit, float(penalty)), jac=True, method='L-BFGS-B',
                   bounds=[(W_MIN, None)] * G, options=dict(FIT_OPTIONS))
    w = res.x / res.x.mean()
    return w, dict(loss=float(res.fun), nit=int(res.nit), success=bool(res.success), message=str(res.message))


def left_out(mats, fly_parts: Mapping[str, np.ndarray], units: Mapping[str, float], penalty: float
             ) -> Tuple[Dict[str, np.ndarray], np.ndarray, List[dict]]:
    """The left-out predictions: for each odor, weights fitted on the other odors to the flies' parts (`fly_parts`,
    flies x odors, each fly re-normalized over the fit odors), the odor's two entries predicted (the weighted dendrites
    relative over all odors). Returns ({part: predictions (odors,)}, weights (odors x glomeruli), the fits' records)."""
    prep = prepare(mats)
    n = prep['M'].shape[2]
    G = prep['M'].shape[1]
    pred = {p: np.full(n, np.nan) for p in STRENGTH_PARTS}
    W = np.zeros((n, G))
    recs = []
    for o in range(n):
        fit = np.array([i for i in range(n) if i != o])
        w, rec = fit_weights(prep, relative_parts(fly_parts, fit), units, fit, penalty)
        ent = weighted_entries(w, prep, np.arange(n))
        for p in STRENGTH_PARTS:
            pred[p][o] = ent[p][0][o]
        W[o] = w
        recs.append(dict(odor_index=o, **rec))
    return pred, W, recs


def predictions_of(weights: np.ndarray, mats) -> Dict[str, np.ndarray]:
    """The left-out predictions of given per-fold weights (odors x glomeruli; fold o predicts odor o), many sets at once
    when `weights` is sets x odors x glomeruli ({part: sets x odors})."""
    prep = prepare(mats)
    W = np.asarray(weights, dtype=float)
    one = W.ndim == 2
    W = W[None] if one else W
    M = prep['M']
    D = np.einsum('sog,fgq->sofq', W, M)
    n = M.shape[2]
    diag = np.arange(n)
    out = {}
    for p in STRENGTH_PARTS:
        a = np.einsum('sog,fgq->sofq', W, prep['X'][p]) / D
        e = a / a.mean(axis=3, keepdims=True)
        x = e[:, diag, :, diag]                     # odors x sets x flies
        out[p] = np.moveaxis(x.mean(axis=2), 0, 1)  # sets x odors
    return {p: v[0] for p, v in out.items()} if one else out


def score_entries(pred: Mapping[str, np.ndarray], ref) -> dict:
    """The strength view's error and its two parts of predictions (part -> odors, or sets x odors) against a reference
    (hal's `score_arrays`): {'view_strength', 'part_rate', 'part_mean'} (arrays for many sets)."""
    from kcq import population_fit as P
    one = np.asarray(pred[STRENGTH_PARTS[0]]).ndim == 1
    arrays = {p: np.atleast_2d(np.asarray(pred[p], dtype=float)) for p in STRENGTH_PARTS}
    df = P.score_arrays(arrays, ref)
    out = dict(view_strength=df.view_strength.to_numpy(), part_rate=df.part_strength_rate.to_numpy(),
               part_mean=df.part_strength_mean.to_numpy())
    return {k: float(v[0]) for k, v in out.items()} if one else out


# ---------------------------------------------------------------------------
# the power check

def plant_truth(w: np.ndarray, mats) -> Dict[str, np.ndarray]:
    """The planted truth: the weighted strength parts of the dendrite flies (their mean, relative over all odors)."""
    prep = prepare(mats)
    return {p: v for p, (v, _) in weighted_entries(w, prep, np.arange(prep['M'].shape[2])).items()}


def noise_fly_parts(truth: Mapping[str, np.ndarray], bouton_parts: Mapping[str, np.ndarray], rng) -> Dict[str, np.ndarray]:
    """Synthetic bouton flies: as many bouton flies drawn with replacement, each the truth plus that fly's deviation from
    the bouton flies' mean ({part: flies x odors})."""
    n = len(next(iter(bouton_parts.values())))
    idx = rng.integers(0, n, size=n)
    out = {}
    for p, A in bouton_parts.items():
        A = np.asarray(A, dtype=float)
        dev = A - np.nanmean(A, axis=0, keepdims=True)
        out[p] = np.asarray(truth[p], dtype=float)[None, :] + dev[idx]
    return out


def synthetic_reference(parts: Mapping[str, np.ndarray], units: Mapping[str, float], odors: Sequence[str],
                        edges: Sequence = ()):
    """A hal FlyReference of synthetic flies' parts in given one-fly units."""
    from kcq import population_fit as P
    n = len(next(iter(parts.values())))
    names = [f'synthetic{i}' for i in range(n)]
    tab = {p: pd.DataFrame(np.asarray(parts[p], dtype=float), index=names, columns=list(odors)) for p in STRENGTH_PARTS}
    return P._reference_from_parts(tab, names, list(odors), tuple(edges), one_fly={p: float(units[p]) for p in tab})


def deal_weights(sizes: Sequence[int], n: int, rng) -> np.ndarray:
    """One deal (E128's): each group of `sizes` to one of `n` glomeruli at random with replacement; a glomerulus's
    weight the summed sizes dealt to it."""
    to = rng.integers(0, n, size=len(sizes))
    return np.bincount(to, weights=np.asarray(sizes, dtype=float), minlength=n).astype(float)


def cv_of(w: np.ndarray) -> float:
    """The coefficient of variation of weights (population SD over the mean)."""
    w = np.asarray(w, dtype=float)
    m = float(w.mean())
    return float(w.std() / m) if m > 0 else np.nan


def uneven_weights(cv: float, z: np.ndarray) -> np.ndarray:
    """A lognormal draw exp(s z) with s solved so that its coefficient of variation (population SD over mean) is `cv`
    exactly; equal weights for cv 0 (the largest reachable, as s grows, where `cv` is beyond it)."""
    from scipy.optimize import brentq
    z = np.asarray(z, dtype=float)
    if not np.isfinite(cv) or cv <= 0:
        return np.ones_like(z)

    def f(s):
        return cv_of(np.exp(s * (z - z.max()))) - cv
    hi = 1.0
    while f(hi) < 0 and hi < 1e3:
        hi *= 2.0
    if f(hi) < 0:
        return np.exp(hi * (z - z.max()))
    s = brentq(f, 0.0, hi, xtol=1e-14, rtol=1e-14)
    return np.exp(s * (z - z.max()))


def power_rule(rows: pd.DataFrame) -> dict:
    """"The bound has power" on one penalty's plants (columns kind, no_step, left_out): per kind the eligible plants
    (no_step >= POWER_ELIGIBLE), the passing ones (no_step - left_out >= MARGIN), the share and whether the kind passes
    or is vacuous; the class."""
    kinds = {}
    for k in PLANT_KINDS:
        x = rows[rows.kind == k]
        el = x[x.no_step >= POWER_ELIGIBLE]
        ok = (el.no_step - el.left_out) >= MARGIN
        ne, npass = int(len(el)), int(ok.sum())
        share = npass / ne if ne else np.nan
        kinds[k] = dict(n=int(len(x)), n_eligible=ne, n_pass=npass, share=share, vacuous=bool(ne == 0),
                        passes=bool(ne > 0 and share >= POWER_SHARE),
                        no_step_median=float(x.no_step.median()) if len(x) else np.nan,
                        left_out_median=float(x.left_out.median()) if len(x) else np.nan)
    live = [k for k, v in kinds.items() if not v['vacuous']]
    n_el = int(sum(kinds[k]['n_eligible'] for k in live))
    min_share = float(min(kinds[k]['share'] for k in live)) if live else np.nan
    holds = bool(live and all(kinds[k]['passes'] for k in live) and n_el >= POWER_MIN_ELIGIBLE)
    return dict(kinds=kinds, vacuous=[k for k in PLANT_KINDS if kinds[k]['vacuous']], n_eligible=n_el,
                min_share=min_share, holds=holds)


def choose_penalty(by_penalty: Mapping[float, dict]) -> float:
    """The lock's rule: the penalty whose smallest passing share over the kinds that are not vacuous is largest, ties to
    the larger penalty."""
    def key(lam):
        r = by_penalty[lam]
        s = r['min_share']
        ok = r['n_eligible'] >= POWER_MIN_ELIGIBLE and s is not None and np.isfinite(s)
        return (float(s) if ok else -np.inf, float(lam))
    return float(max(by_penalty, key=key))


# ---------------------------------------------------------------------------
# the features

def top_share(values: pd.DataFrame, k: int = TOP_K) -> pd.Series:
    """Per odor (column), the share of its summed positive response (values clipped at 0) in its k strongest units."""
    X = values.clip(lower=0)
    out = {}
    for o in X.columns:
        x = X[o].dropna().to_numpy(dtype=float)
        s = float(x.sum())
        out[o] = float(np.sort(x)[::-1][:k].sum() / s) if s > 0 else np.nan
    return pd.Series(out)


def flies_mean_table(pops: Sequence, odors: Sequence[str]) -> pd.DataFrame:
    """Per unit (glomerulus) and odor, the mean over the flies that recorded it."""
    tabs = []
    for p in pops:
        V = p.values.reindex(columns=list(odors)).astype(float).copy()
        V.index = V.index.astype(str)
        tabs.append(V)
    return pd.concat(tabs).groupby(level=0, sort=True).mean()


def summed_positive(values: pd.DataFrame) -> pd.Series:
    """Per odor, the summed positive response (values clipped at 0) over the units."""
    return values.clip(lower=0).sum(axis=0, min_count=1)


def linear_left_out(x: np.ndarray, f: Optional[np.ndarray], fly_parts: np.ndarray) -> np.ndarray:
    """Left-out predictions of a + b x (+ c f): for each odor o, least squares on the other odors against the flies'
    entries relative over those odors (`fly_parts`, flies x odors, each fly re-normalized, then averaged), predicted at
    o and converted to the all-odor scale (p_all = n p / (n - 1 + p))."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    A = np.column_stack([np.ones(n), x] + ([np.asarray(f, dtype=float)] if f is not None else []))
    pred = np.full(n, np.nan)
    for o in range(n):
        fit = np.array([i for i in range(n) if i != o])
        t = relative_parts({'x': fly_parts}, fit)['x']
        coef, *_ = np.linalg.lstsq(A[fit], t[fit], rcond=None)
        p16 = float(A[o] @ coef)
        pred[o] = n * p16 / (n - 1 + p16)
    return pred


class FeatureInputs:
    """The upstream entries and the features of one step.

    Attributes:
        x: {set: {part: entries (odors,)}} of the upstream flies ('all', 'fit', 'held').
        single: {feature: {part: values (odors,)}}, taken whole (with the upstream entry of all flies).
        halves: {feature: {half: {part: values}}}, taken from the half the upstream entry is not (both ways).
    """

    def __init__(self, x, single, halves):
        self.x, self.single, self.halves = x, single, halves

    def features(self) -> List[str]:
        return [f for f in FEATURES if f in self.single or f in self.halves]


def feature_prediction(I: FeatureInputs, name: Optional[str], fly_parts: Mapping[str, np.ndarray],
                       perm: Optional[np.ndarray] = None, ways: bool = False) -> Dict[str, np.ndarray]:
    """A feature's left-out predictions per part (`linear_left_out`); `name` None for the upstream entry rescaled alone
    (`ways` True: both ways on the halves, the two predictions averaged); `perm` the feature's values permuted among the
    odors (value of odor perm[i] at odor i)."""
    out = {}
    for p in STRENGTH_PARTS:
        if name is None and not ways:
            out[p] = linear_left_out(I.x['all'][p], None, fly_parts[p])
        elif name is None:
            out[p] = np.mean([linear_left_out(I.x[h][p], None, fly_parts[p]) for h in HALVES], axis=0)
        elif name in I.single:
            f = np.asarray(I.single[name][p], dtype=float)
            out[p] = linear_left_out(I.x['all'][p], f if perm is None else f[perm], fly_parts[p])
        else:
            preds = []
            for a, b in (('fit', 'held'), ('held', 'fit')):
                f = np.asarray(I.halves[name][b][p], dtype=float)
                preds.append(linear_left_out(I.x[a][p], f if perm is None else f[perm], fly_parts[p]))
            out[p] = np.mean(preds, axis=0)
    return out


def mean_parts(pops: Sequence, odors: Sequence[str]) -> Dict[str, np.ndarray]:
    """The flies' mean strength parts (hal's population_views, relative over `odors`)."""
    from kcq import population_fit as P
    vs = [P.population_views(p, list(odors), parts=list(STRENGTH_PARTS)) for p in pops]
    return {q: np.nanmean([v[q].reindex(list(odors)).to_numpy(dtype=float) for v in vs], axis=0)
            for q in STRENGTH_PARTS}


def fly_parts_of(pops: Sequence, odors: Sequence[str]) -> Dict[str, np.ndarray]:
    """Each fly's strength parts (flies x odors)."""
    from kcq import population_fit as P
    vs = [P.population_views(p, list(odors), parts=list(STRENGTH_PARTS)) for p in pops]
    return {q: np.array([v[q].reindex(list(odors)).to_numpy(dtype=float) for v in vs]) for q in STRENGTH_PARTS}


def orn_features(ctx: Ctx) -> Dict[str, Dict[str, np.ndarray]]:
    """The two ORN features ('orn_summed', 'orn_top3_share'), the same for both parts."""
    tab = flies_mean_table(ctx.orn_pops(), ctx.odors)
    s = summed_positive(tab).reindex(ctx.odors).to_numpy(dtype=float)
    t = top_share(tab).reindex(ctx.odors).to_numpy(dtype=float)
    return {'orn_summed': {p: s for p in STRENGTH_PARTS}, 'orn_top3_share': {p: t for p in STRENGTH_PARTS}}


def dendrite_step_inputs(ctx: Ctx) -> FeatureInputs:
    """The step from the dendrites to the boutons: the dendrite flies' entries (all and each half), the ORN features
    taken whole, the dendrite top-3 share and the other part's entry from the other half."""
    x = {s: mean_parts(ctx.dendrite_set(s), ctx.odors) for s in FLY_SETS}
    top = {h: top_share(flies_mean_table(ctx.dendrite_set(h), ctx.odors)).reindex(ctx.odors).to_numpy(dtype=float)
           for h in HALVES}
    halves = {'dendrite_top3_share': {h: {p: top[h] for p in STRENGTH_PARTS} for h in HALVES},
              'other_part': {h: {p: x[h][OTHER_PART[p]] for p in STRENGTH_PARTS} for h in HALVES}}
    return FeatureInputs(x, orn_features(ctx), halves)


def bouton_step_inputs(ctx: Ctx) -> FeatureInputs:
    """The step from the boutons to the claws (beside): the bouton flies' entries (all and each half), the ORN features
    and the dendrite top-3 share (all seven dendrite flies) taken whole, the other part's entry from the other bouton
    half."""
    ref = ctx.refs['all']
    x = {'all': {p: ref.mean[p].reindex(ctx.odors).to_numpy(dtype=float) for p in STRENGTH_PARTS}}
    for h in HALVES:
        keep = list(ctx.names['pn_bouton'][h])
        x[h] = {p: ref.parts[p].loc[keep].reindex(columns=ctx.odors).mean(axis=0).to_numpy(dtype=float)
                for p in STRENGTH_PARTS}
    single = dict(orn_features(ctx))
    top = top_share(flies_mean_table(ctx.dendrite_set('all'), ctx.odors)).reindex(ctx.odors).to_numpy(dtype=float)
    single['dendrite_top3_share'] = {p: top for p in STRENGTH_PARTS}
    halves = {'other_part': {h: {p: x[h][OTHER_PART[p]] for p in STRENGTH_PARTS} for h in HALVES}}
    return FeatureInputs(x, single, halves)


def feature_readings(I: FeatureInputs, fly_parts: Mapping[str, np.ndarray], refs: Mapping[str, object],
                     no_step: Mapping[str, float], n_perm: int = N_PERM, seed: int = PERM_SEED
                     ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Each feature's left-out reading against each fly set, the upstream entry rescaled alone (all flies, and both ways
    on the halves) beside; the permuted best's null (`n_perm` permutations of the odors, the same for every feature);
    the predictions by odor. Returns (readings, null, predictions)."""
    from experiments import E142_step0 as S0
    rows, preds = [], []
    forms = [(None, False, 'rescaled_alone'), (None, True, 'rescaled_alone_halves')] + [
        (f, False, f) for f in I.features()]
    for name, ways, label in forms:
        pr = feature_prediction(I, name, fly_parts, ways=ways)
        for s, ref in refs.items():
            sc = score_entries(pr, ref)
            rows.append(dict(form=label, feature=label if name else '', fly_set=s, left_out=sc['view_strength'],
                             left_out_rate=sc['part_rate'], left_out_mean=sc['part_mean'], no_step=no_step[s],
                             below_no_step_by=no_step[s] - sc['view_strength']))
        for p in STRENGTH_PARTS:
            preds += [dict(form=label, part=PART_KEY[p], odor=o, prediction=float(v)) for o, v in
                      zip(refs['all'].mean[p].index, pr[p])]
    n = len(next(iter(I.x['all'].values())))
    perms = S0.permutations(n, n_perm, seed)
    nl = []
    feats = I.features()
    for i, pm in enumerate(perms):
        best, which = np.inf, ''
        rec = dict(perm=i)
        for f in feats:
            v = score_entries(feature_prediction(I, f, fly_parts, perm=pm), refs['all'])['view_strength']
            rec[f] = v
            if v < best:
                best, which = v, f
        rec.update(best=best, best_feature=which)
        nl.append(rec)
    return pd.DataFrame(rows), pd.DataFrame(nl), pd.DataFrame(preds)


# ---------------------------------------------------------------------------
# two sets at one stage

def resampled_errors(scored: Mapping[str, np.ndarray], reference: Mapping[str, np.ndarray],
                     units: Mapping[str, float], n_boot: int = N_BOOT_SETS, seed: int = SETS_SEED) -> np.ndarray:
    """The strength view's error of the scored set's mean parts against the reference set's mean, with both sets'
    flies resampled with replacement independently, in fixed units (n_boot,)."""
    rng = np.random.default_rng(seed)
    ns = len(next(iter(scored.values())))
    nr = len(next(iter(reference.values())))
    i = rng.integers(0, ns, size=(n_boot, ns))
    j = rng.integers(0, nr, size=(n_boot, nr))
    sq = []
    for p in STRENGTH_PARTS:
        S = np.nanmean(np.asarray(scored[p], dtype=float)[i], axis=1)
        R = np.nanmean(np.asarray(reference[p], dtype=float)[j], axis=1)
        sq.append(np.nanmean((S - R) ** 2, axis=1) / float(units[p]) ** 2)
    return np.sqrt(np.mean(sq, axis=0))


def measured_error(scored: Mapping[str, np.ndarray], reference: Mapping[str, np.ndarray],
                   units: Mapping[str, float]) -> dict:
    """The strength view's error (and parts) of the scored set's mean parts against the reference set's mean."""
    out = {}
    for p in STRENGTH_PARTS:
        d = np.nanmean(np.asarray(scored[p], dtype=float), axis=0) - np.nanmean(np.asarray(reference[p], dtype=float),
                                                                              axis=0)
        out[f'part_{PART_KEY[p]}'] = float(np.sqrt(np.nanmean(d ** 2)) / float(units[p]))
    out['view_strength'] = float(np.sqrt(np.mean([out[f'part_{PART_KEY[p]}'] ** 2 for p in STRENGTH_PARTS])))
    return out


def interval(x: np.ndarray, q: Sequence[float] = SETS_Q) -> List[float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return [float(np.percentile(x, q[0])), float(np.percentile(x, q[1]))]


def soma_sets() -> Dict[str, Dict[str, np.ndarray]]:
    """Remy's megamat somata (E134's lock criterion) and the control flies on the 12 shared odors: each fly's strength
    parts relative over the 12 ({'megamat': {part: flies x 12}, 'control': ...}), with the populations."""
    import dataclasses
    from experiments import E134_flies as F
    from kcq import population_fit as P
    odors = list(F.MEGAMAT12)
    mm = F.recorded(('kc_soma', PANEL, 'wild_type'))
    mmc = dataclasses.replace(mm, criterion=F.megamat_criterion(MEGAMAT12_CRITERION, mm))
    pops = {'megamat': P.fly_populations(mmc), 'control': F.flies_of(('kc_soma', 'apl16', 'control'))}
    return {k: dict(parts=fly_parts_of(v, odors), pops=v, odors=odors) for k, v in pops.items()}


# ---------------------------------------------------------------------------
# the classes

def weight_class(bound: Mapping, power_holds: bool) -> dict:
    """"A weight per glomerulus makes the step" from bound.json's numbers; "not readable" without power."""
    lo, ns = bound['left_out'], bound['no_step']
    margin = {h: bool(ns[h] - lo[h] >= MARGIN) for h in HALVES}
    below = bool(lo['all'] < bound['null']['all']['null_p05'])
    holds = bool(all(margin.values()) and below)
    cls = ('not readable' if not power_holds else ('makes the step' if holds else 'does not make the step'))
    return dict(cls=cls, holds=bool(holds and power_holds), holds_on_numbers=holds, margin_met=margin,
                below_no_step_by={h: float(ns[h] - lo[h]) for h in FLY_SETS}, below_null_p05=below,
                null_p05=bound['null']['all']['null_p05'], penalty=bound['penalty'], power_holds=bool(power_holds),
                n_forms=1)


def feature_class(features: Mapping) -> dict:
    """"A population feature makes the step" from features.json's numbers."""
    p05 = features['null_best_p05']
    per = {}
    for f in FEATURES:
        r = features['forms'][f]
        margin = {h: bool(r['no_step'][h] - r['left_out'][h] >= MARGIN) for h in HALVES}
        below = bool(r['left_out']['all'] < p05)
        per[f] = dict(margin_met=margin, below_null_p05=below, holds=bool(all(margin.values()) and below),
                      left_out=r['left_out'], below_no_step_by={h: float(r['no_step'][h] - r['left_out'][h])
                                                                for h in FLY_SETS})
    passing = [f for f in FEATURES if per[f]['holds']]
    return dict(cls='makes the step' if passing else 'does not make the step', holds=bool(passing), passing=passing,
                features=per, null_best_p05=p05, n_forms=len(FEATURES))


def overlap(a: Sequence[float], b: Sequence[float]) -> bool:
    """Whether two intervals meet."""
    return bool(max(a[0], b[0]) <= min(a[1], b[1]))


def outcome(power_holds: bool, weight_holds: bool, feature_holds: bool, not_told_apart: bool) -> List[str]:
    """The outcome letters by OUTCOME_RULE."""
    out = []
    if power_holds and weight_holds:
        out.append('a')
    if feature_holds:
        out.append('b')
    if out:
        return out
    if power_holds:
        return ['c'] if not_told_apart else ['d']
    return ['e']


# ---------------------------------------------------------------------------
# the parts

def _identity_row(name: str, got: float, saved: Optional[float]) -> dict:
    spec = IDENTITIES[name]
    q, dg = float(spec['quoted']), int(spec['digits'])
    ok_q = bool(np.isfinite(got) and round(float(got), dg) == round(q, dg))
    ok_s = bool(saved is None or (np.isfinite(got) and abs(float(got) - float(saved)) <= IDENTITY_TOL))
    return dict(name=name, quoted=q, digits=dg, reproduced=float(got), saved=saved, matches_quoted=ok_q,
                matches_saved=ok_s, reproduces=bool(ok_q and ok_s), source=spec['source'])


def part_identities() -> dict:
    """Design item 1 (identities.json; stops on a failure)."""
    from experiments import E142_step0 as S0
    from kcq import population_fit as P
    t0 = time.time()
    ctx = context()
    res = ctx.lock['inputs']['results']
    rows, checks = [], {}
    ns = pd.read_csv(remap_path(res['e142_no_step_entries']['rel']))
    ns = ns[(ns.upstream == 'pn_dendrite') & (ns.downstream == 'pn_bouton')]
    mine = ctx.no_step_views()
    for s in FLY_SETS:
        hal = float(P.no_step_score(ctx.pops['pn_dendrite'], ctx.refs[s]).views['strength'])
        own = score_entries(mine, ctx.refs[s])['view_strength']
        checks[f'no_step_{s}_own_code'] = dict(hal=hal, own=own, agree=bool(abs(hal - own) <= 1e-10))
        rows.append(_identity_row(f'no_step_{s}', hal, float(ns[ns.fly_set == s].view_error.iloc[0])))
    ys = pd.read_csv(remap_path(res['e142_yardsticks']['rel']))
    ys = ys[(ys.stage == 'pn_bouton') & (ys.glomeruli == 'dendrites') & (ys.calyx == HEMIBRAIN) & (ys.fly_set == 'all')]
    wh = ctx.counts()[HEMIBRAIN].to_numpy(dtype=float)
    hb = float(S0.yardstick_scores(wh[None, :], ctx.mats, ctx.refs['all']).iloc[0].view_strength)
    own = score_entries(plant_truth(wh, ctx.prep), ctx.refs['all'])['view_strength']
    checks['hemibrain_own_code'] = dict(e142=hb, own=own, agree=bool(abs(hb - own) <= 1e-10))
    rows.append(_identity_row('hemibrain_counts', hb, float(ys.counts.iloc[0])))
    sets = soma_sets()
    mc = pd.read_csv(remap_path(res['e134_flies_megamat_vs_control']['rel']))
    mc = mc[(mc.criterion == MEGAMAT12_CRITERION) & (mc.name == 'view_strength')].set_index('reference')
    for name, (scored, refside) in (('megamat_in_control_units', ('megamat', 'control')),
                                    ('control_in_megamat_units', ('control', 'megamat'))):
        ref = P.fly_reference(sets[refside]['pops'], sets[refside]['odors'],
                              _breadth(len(sets[refside]['odors'])), parts=list(STRENGTH_PARTS))
        units = {p: float(ref.one_fly[p]) for p in STRENGTH_PARTS}
        got = measured_error(sets[scored]['parts'], sets[refside]['parts'], units)['view_strength']
        rows.append(_identity_row(name, got, float(mc.loc[refside, 'mean_views'])))
        fl = float(P.floor(ref, n_boot=1000, seed=0).view_strength.mean())
        rows.append(_identity_row(f'{refside}_floor', fl, float(mc.loc[refside, 'ref_floor'])))
    df = pd.DataFrame(rows)
    ok = bool(df.reproduces.all() and all(c['agree'] for c in checks.values()))
    _write_csv(df, 'identities.csv')
    _write_json('identities.json', dict(all_reproduce=ok, rows=df.to_dict('records'), own_code_checks=checks,
                                        time=_now(), host=_host()))
    provenance('identities', ctx, seconds=round(time.time() - t0, 1), all_reproduce=ok)
    print(df[['name', 'quoted', 'reproduced', 'saved', 'reproduces']].to_string(index=False), flush=True)
    if not ok:
        raise RuntimeError('an identity does not reproduce (identities.json); stopping')
    return dict(rows=df, checks=checks)


def _breadth(n: int):
    from experiments import E134_population_baseline as B
    return B.breadth_edges(n)


def _require_identities() -> None:
    if not _read_json('identities.json').get('all_reproduce'):
        raise PermissionError('the identities do not all reproduce; no later part runs')


def plants(ctx: Ctx, cvs: Mapping[float, float]) -> List[dict]:
    """The power check's plants: per plant its kind, index, the truth's weights (per penalty for the uneven kind, all
    from the same standard normals) and its noise draw (the bouton flies' positions)."""
    G = len(ctx.glomeruli)
    n_b = len(ctx.pops['pn_bouton'])
    out = []
    wh = ctx.counts()[HEMIBRAIN].to_numpy(dtype=float)
    for i in range(N_HEMIBRAIN_PLANTS):
        rng = np.random.default_rng([POWER_SEED, 0, i])
        out.append(dict(kind='hemibrain', plant=i, fly='', weights={lam: wh for lam in PENALTIES},
                        noise=rng.integers(0, n_b, size=n_b)))
    sizes = ctx.group_sizes()
    k = 0
    for f in ctx.lock['recordings']['pn_bouton']['flies']:
        for d in range(N_DEALS_PER_FLY):
            rng = np.random.default_rng([POWER_SEED, 1, k])
            w = deal_weights(sizes[f], G, rng)
            out.append(dict(kind='deals', plant=k, fly=f, weights={lam: w for lam in PENALTIES},
                            noise=rng.integers(0, n_b, size=n_b)))
            k += 1
    for i in range(N_UNEVEN):
        rng = np.random.default_rng([POWER_SEED, 2, i])
        z = rng.standard_normal(G)
        out.append(dict(kind='uneven', plant=i, fly='', weights={lam: uneven_weights(cvs[lam], z) for lam in PENALTIES},
                        noise=rng.integers(0, n_b, size=n_b)))
    return out


def plant_reading(ctx: Ctx, w_true: np.ndarray, noise: np.ndarray, penalties: Sequence[float],
                  truth_prep: Optional[Mapping] = None) -> Tuple[float, Dict[float, dict], Dict[float, np.ndarray]]:
    """One plant: the synthetic bouton flies (the truth, from `truth_prep` or the measured dendrite flies, plus the
    bouton flies' deviations at positions `noise`), the planted no-step, and per penalty the left-out reading against
    all five synthetic flies. Returns (no-step, {penalty: reading}, {penalty: fold weights})."""
    truth = plant_truth(w_true, truth_prep if truth_prep is not None else ctx.prep)
    parts = {}
    for p, A in ctx.bouton_parts.items():
        dev = A - np.nanmean(A, axis=0, keepdims=True)
        parts[p] = truth[p][None, :] + dev[np.asarray(noise, dtype=int)]
    sref = synthetic_reference(parts, ctx.units, ctx.odors)
    ns = score_entries(ctx.no_step_views(), sref)['view_strength']
    readings, weights = {}, {}
    for lam in penalties:
        pred, W, recs = left_out(ctx.prep, parts, ctx.units, lam)
        sc = score_entries(pred, sref)
        readings[lam] = dict(left_out=sc['view_strength'], left_out_rate=sc['part_rate'], left_out_mean=sc['part_mean'],
                             n_not_converged=int(sum(not r['success'] for r in recs)),
                             max_nit=int(max(r['nit'] for r in recs)))
        weights[lam] = W
    return ns, readings, weights


def _plant_job(args):
    """One plant in a worker process (`part_power`'s pool)."""
    ctx, pl, penalties, beside_penalty = _JOB['ctx'], args, _JOB['penalties'], _JOB['beside']
    rows, wsave = [], {}
    if beside_penalty is None:
        # each penalty has its own truth for the uneven kind
        for lam in penalties:
            ns, rd, W = plant_reading(ctx, pl['weights'][lam], pl['noise'], [lam])
            rows.append(dict(penalty=lam, kind=pl['kind'], plant=pl['plant'], fly=pl['fly'],
                             weights_cv=cv_of(pl['weights'][lam]), no_step=ns, **rd[lam]))
            wsave[lam] = W
    else:
        lam = beside_penalty
        rng = np.random.default_rng([POWER_SEED, PLANT_KINDS.index(pl['kind']), pl['plant'], 1])
        nd = ctx.prep['M'].shape[0]
        idx = rng.integers(0, nd, size=nd)
        ns, rd, W = plant_reading(ctx, pl['weights'][lam], pl['noise'], [lam], truth_prep=resample_prep(ctx.prep, idx))
        rows.append(dict(penalty=lam, kind=pl['kind'], plant=pl['plant'], fly=pl['fly'],
                         weights_cv=cv_of(pl['weights'][lam]), no_step=ns, dendrite_draw=' '.join(map(str, idx)),
                         **rd[lam]))
        wsave[lam] = W
    return rows, (pl['kind'], pl['plant'], wsave)


_JOB: Dict[str, object] = {}


def _run_plants(ctx: Ctx, pls: Sequence[dict], penalties: Sequence[float], beside: Optional[float],
                workers: int) -> Tuple[List[dict], list]:
    """The plants' readings, in `workers` processes (fork; the context shared)."""
    _JOB.update(ctx=ctx, penalties=list(penalties), beside=beside)
    rows, saves = [], []
    if workers > 1:
        import multiprocessing as mp
        with mp.get_context('fork').Pool(workers) as pool:
            for i, (r, s) in enumerate(pool.imap(_plant_job, pls, chunksize=1)):
                rows += r
                saves.append(s)
                if (i + 1) % 20 == 0:
                    print(f'  plants: {i + 1} of {len(pls)}', flush=True)
    else:
        for i, pl in enumerate(pls):
            r, s = _plant_job(pl)
            rows += r
            saves.append(s)
            if (i + 1) % 20 == 0:
                print(f'  plants: {i + 1} of {len(pls)}', flush=True)
    return rows, saves


def part_power(workers: int = 2) -> dict:
    """Design item 2 (power_unevenness.csv, power.csv, power_summary.csv, power_beside.csv, power.json): the uneven
    draws' unevenness from the weights fitted to all 17 odors at each penalty (no prediction or error of a fit is
    read), the plants at every penalty, the penalty the lock's rule picks, the class at it, and beside the plants with
    the dendrite flies resampled at the chosen penalty."""
    _require_identities()
    t0 = time.time()
    ctx = context()
    cvs = {}
    for lam in PENALTIES:
        w, rec = fit_weights(ctx.prep, relative_parts(ctx.bouton_parts, ctx.all_idx), ctx.units, ctx.all_idx, lam)
        cvs[lam] = cv_of(w)
    _write_csv(pd.DataFrame([dict(penalty=k, fitted_weights_cv=v) for k, v in cvs.items()]), 'power_unevenness.csv')
    pls = plants(ctx, cvs)
    print(f'power: {len(pls)} plants x {len(PENALTIES)} penalties', flush=True)
    rows, saves = _run_plants(ctx, pls, PENALTIES, None, workers)
    df = pd.DataFrame(rows)
    df['below_no_step_by'] = df.no_step - df.left_out
    df['eligible'] = df.no_step >= POWER_ELIGIBLE
    df['passes'] = df.eligible & (df.below_no_step_by >= MARGIN)
    by = {float(lam): power_rule(df[df.penalty == lam]) for lam in PENALTIES}
    lam_star = choose_penalty(by)
    print(f'power: chosen penalty {lam_star} ({time.time() - t0:.0f} s); beside: the dendrites resampled', flush=True)
    brows, bsaves = _run_plants(ctx, pls, PENALTIES, lam_star, workers)
    bdf = pd.DataFrame(brows)
    bdf['below_no_step_by'] = bdf.no_step - bdf.left_out
    bdf['eligible'] = bdf.no_step >= POWER_ELIGIBLE
    bdf['passes'] = bdf.eligible & (bdf.below_no_step_by >= MARGIN)
    beside_rule = power_rule(bdf)
    summ = []
    for lam, r in by.items():
        for k, v in r['kinds'].items():
            summ.append(dict(reading='the plants', penalty=lam, kind=k, **v, class_holds=r['holds'],
                             min_share=r['min_share'], n_eligible_all=r['n_eligible'], chosen=bool(lam == lam_star)))
    for k, v in beside_rule['kinds'].items():
        summ.append(dict(reading='beside: the dendrites resampled', penalty=lam_star, kind=k, **v,
                         class_holds=beside_rule['holds'], min_share=beside_rule['min_share'],
                         n_eligible_all=beside_rule['n_eligible'], chosen=True))
    sm = pd.DataFrame(summ)
    _write_csv(df, 'power.csv')
    _write_csv(bdf, 'power_beside.csv')
    _write_csv(sm, 'power_summary.csv')
    np.savez_compressed(local_dir() / 'power_weights.npz',
                        truth=np.array([[pl['weights'][lam] for lam in PENALTIES] for pl in pls]),
                        kinds=np.array([pl['kind'] for pl in pls]), plants=np.array([pl['plant'] for pl in pls]),
                        penalties=np.array(PENALTIES),
                        folds=np.array([[s[2][lam] for lam in PENALTIES] for s in saves]))
    out = dict(penalty=lam_star, rule=CLASSES['power']['rule'], by_penalty={str(k): v for k, v in by.items()},
               holds=by[lam_star]['holds'], at_chosen=by[lam_star], beside_dendrites_resampled=beside_rule,
               unevenness={str(k): v for k, v in cvs.items()}, n_plants=len(pls), time=_now(), host=_host())
    _write_json('power.json', out)
    provenance('power', ctx, seconds=round(time.time() - t0, 1), workers=workers, penalty=lam_star)
    print(json.dumps(dict(penalty=lam_star, holds=out['holds'],
                          kinds={k: (v['n_eligible'], v['n_pass'], v['share']) for k, v in
                                 by[lam_star]['kinds'].items()}), default=str), flush=True)
    return out


def _chosen_penalty() -> float:
    return float(_read_json('power.json')['penalty'])


def bound_reading(ctx: Ctx, prep: Mapping, fly_parts: Mapping[str, np.ndarray], units: Mapping[str, float],
                  refs: Mapping[str, object], penalty: float, no_step: Mapping[str, float], null: bool = True
                  ) -> Tuple[dict, pd.DataFrame, np.ndarray, Dict[str, np.ndarray]]:
    """The weight bound of one step: the left-out predictions at `penalty` scored against every fly set, beside the
    no-step; the permuted weights' null against every set. Returns (record, null draws, fold weights, predictions)."""
    from experiments import E142_step0 as S0
    pred, W, recs = left_out(prep, fly_parts, units, penalty)
    lo = {s: score_entries(pred, r) for s, r in refs.items()}
    rec = dict(penalty=penalty, left_out={s: v['view_strength'] for s, v in lo.items()},
               left_out_parts={s: dict(rate=v['part_rate'], mean=v['part_mean']) for s, v in lo.items()},
               no_step=dict(no_step), n_not_converged=int(sum(not r['success'] for r in recs)),
               max_nit=int(max(r['nit'] for r in recs)))
    nl = pd.DataFrame()
    if null:
        perms = S0.permutations(W.shape[1], N_PERM, PERM_SEED)
        pw = predictions_of(np.transpose(W[:, perms], (1, 0, 2)), prep)   # perms x folds x glomeruli
        cols = {}
        rec['null'] = {}
        for s, r in refs.items():
            v = score_entries(pw, r)['view_strength']
            cols[s] = v
            rec['null'][s] = S0.null_summary(rec['left_out'][s], v, NULL_PERCENTILE)
        nl = pd.DataFrame(dict(perm=np.arange(len(perms)), **cols))
    return rec, nl, W, pred


def part_bound() -> dict:
    """Design item 3 (bound.csv, bound_null.csv, bound_weights.csv, bound_odors.csv, bound_by_penalty.csv,
    bound.json)."""
    from kcq import population_fit as P
    _require_identities()
    t0 = time.time()
    ctx = context()
    lam = _chosen_penalty()
    ns = {s: float(P.no_step_score(ctx.pops['pn_dendrite'], ctx.refs[s]).views['strength']) for s in FLY_SETS}
    rec, nl, W, pred = bound_reading(ctx, ctx.prep, ctx.bouton_parts, ctx.units, ctx.refs, lam, ns)
    # beside: fitted on all odors to the fit half, read against the held-out half; and in sample to all five
    fit_half = ctx.refs['fit']
    fp = {p: fit_half.parts[p].reindex(columns=ctx.odors).to_numpy(dtype=float) for p in STRENGTH_PARTS}
    w_half, r_half = fit_weights(ctx.prep, relative_parts(fp, ctx.all_idx), ctx.units, ctx.all_idx, lam)
    half_pred = plant_truth(w_half, ctx.prep)
    rec['fit_half_to_held'] = dict(held=score_entries(half_pred, ctx.refs['held'])['view_strength'],
                                   fit_in_sample=score_entries(half_pred, ctx.refs['fit'])['view_strength'],
                                   no_step_held=ns['held'], fit=r_half)
    w_all, r_all = fit_weights(ctx.prep, relative_parts(ctx.bouton_parts, ctx.all_idx), ctx.units, ctx.all_idx, lam)
    all_pred = plant_truth(w_all, ctx.prep)
    rec['in_sample_all'] = {s: score_entries(all_pred, ctx.refs[s])['view_strength'] for s in FLY_SETS}
    rec['in_sample_fit'] = r_all
    rec['weights_cv'] = dict(in_sample_all=cv_of(w_all), fit_half=cv_of(w_half),
                             left_out_median=float(np.median([cv_of(w) for w in W])))
    rows = [dict(reading='left_out', fly_set=s, error=rec['left_out'][s], no_step=ns[s],
                 below_no_step_by=ns[s] - rec['left_out'][s], null_p05=rec['null'][s]['null_p05'],
                 null_median=rec['null'][s]['null_median'],
                 null_percentile_of_actual=rec['null'][s]['null_percentile_of_actual'],
                 part_rate=rec['left_out_parts'][s]['rate'], part_mean=rec['left_out_parts'][s]['mean'])
            for s in FLY_SETS]
    rows += [dict(reading='in_sample_all_odors', fly_set=s, error=rec['in_sample_all'][s], no_step=ns[s],
                  below_no_step_by=ns[s] - rec['in_sample_all'][s]) for s in FLY_SETS]
    rows.append(dict(reading='fit_half_all_odors_read_on_held', fly_set='held', error=rec['fit_half_to_held']['held'],
                     no_step=ns['held'], below_no_step_by=ns['held'] - rec['fit_half_to_held']['held']))
    wt = pd.DataFrame(dict(glomerulus=ctx.glomeruli, in_sample_all=w_all, fit_half=w_half,
                           **{f'left_out_{o}': W[i] for i, o in enumerate(ctx.odors)}))
    ns_views = ctx.no_step_views()
    od = []
    for p in STRENGTH_PARTS:
        for i, o in enumerate(ctx.odors):
            od.append(dict(part=PART_KEY[p], odor=o, flies=float(ctx.refs['all'].mean[p][o]),
                           no_step=float(ns_views[p][i]), left_out=float(pred[p][i]), in_sample=float(all_pred[p][i])))
    # beside, classing nothing: every penalty of the grid
    byp = []
    for l2 in PENALTIES:
        pr, _, recs = left_out(ctx.prep, ctx.bouton_parts, ctx.units, l2)
        for s in FLY_SETS:
            byp.append(dict(penalty=l2, chosen=bool(l2 == lam), fly_set=s,
                            left_out=score_entries(pr, ctx.refs[s])['view_strength'], no_step=ns[s]))
    _write_csv(pd.DataFrame(rows), 'bound.csv')
    _write_csv(nl, 'bound_null.csv')
    _write_csv(wt, 'bound_weights.csv')
    _write_csv(pd.DataFrame(od), 'bound_odors.csv')
    _write_csv(pd.DataFrame(byp), 'bound_by_penalty.csv')
    rec.update(time=_now(), host=_host())
    _write_json('bound.json', rec)
    provenance('bound', ctx, seconds=round(time.time() - t0, 1), penalty=lam)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)
    return rec


def part_features() -> dict:
    """Design item 4 (features.csv, features_odors.csv, features_null.csv, features_values.csv, features.json)."""
    from kcq import population_fit as P
    _require_identities()
    t0 = time.time()
    ctx = context()
    ns = {s: float(P.no_step_score(ctx.pops['pn_dendrite'], ctx.refs[s]).views['strength']) for s in FLY_SETS}
    I = dendrite_step_inputs(ctx)
    rd, nl, pr = feature_readings(I, ctx.bouton_parts, ctx.refs, ns)
    p05 = float(np.percentile(nl.best.to_numpy(dtype=float), NULL_PERCENTILE))
    forms = {}
    for f, x in rd.groupby('form', sort=False):
        o = x.set_index('fly_set')
        forms[f] = dict(left_out={s: float(o.loc[s, 'left_out']) for s in FLY_SETS},
                        no_step={s: float(o.loc[s, 'no_step']) for s in FLY_SETS})
    vals = []
    for f in FEATURES:
        if f in I.single:
            vals += [dict(feature=f, source='whole', part=PART_KEY[p], odor=o, value=float(v))
                     for p in STRENGTH_PARTS for o, v in zip(ctx.odors, I.single[f][p])]
        else:
            vals += [dict(feature=f, source=h, part=PART_KEY[p], odor=o, value=float(v))
                     for h in HALVES for p in STRENGTH_PARTS for o, v in zip(ctx.odors, I.halves[f][h][p])]
    out = dict(forms=forms, null_best_p05=p05, null_best_median=float(nl.best.median()),
               null_best_feature_counts=nl.best_feature.value_counts().to_dict(), n_perm=N_PERM, time=_now(),
               host=_host())
    _write_csv(rd, 'features.csv')
    _write_csv(pr, 'features_odors.csv')
    _write_csv(nl, 'features_null.csv')
    _write_csv(pd.DataFrame(vals), 'features_values.csv')
    _write_json('features.json', out)
    provenance('features', ctx, seconds=round(time.time() - t0, 1))
    print(rd.to_string(index=False), f'\npermuted best p05 {p05:.3f}', flush=True)
    return out


def part_sets() -> dict:
    """Design item 5 (sets.csv, sets_draws.csv, sets_odors.csv, sets.json)."""
    from kcq import population_fit as P
    _require_identities()
    t0 = time.time()
    ctx = context()
    sets = soma_sets()
    refs = {k: P.fly_reference(v['pops'], v['odors'], _breadth(len(v['odors'])), parts=list(STRENGTH_PARTS))
            for k, v in sets.items()}
    units = {k: {p: float(r.one_fly[p]) for p in STRENGTH_PARTS} for k, r in refs.items()}
    comps = {'megamat_in_control_units': ('megamat', 'control'), 'control_in_megamat_units': ('control', 'megamat')}
    rows, draws = [], {}
    for name, (sc, rf) in comps.items():
        m = measured_error(sets[sc]['parts'], sets[rf]['parts'], units[rf])
        d = resampled_errors(sets[sc]['parts'], sets[rf]['parts'], units[rf])
        iv = interval(d)
        rows.append(dict(comparison=name, scored=sc, reference=rf, n_scored=len(sets[sc]['pops']),
                         n_reference=len(sets[rf]['pops']), n_odors=len(sets[rf]['odors']), measured=m['view_strength'],
                         part_rate=m['part_rate'], part_mean=m['part_mean'], lo=iv[0], hi=iv[1],
                         draws_median=float(np.median(d))))
        draws[name] = d
    den = fly_parts_of(ctx.pops['pn_dendrite'], ctx.odors)
    m = measured_error(den, ctx.bouton_parts, ctx.units)
    d = resampled_errors(den, ctx.bouton_parts, ctx.units)
    iv = interval(d)
    rows.append(dict(comparison='no_step', scored='pn_dendrite', reference='pn_bouton', n_scored=len(den['strength/rate']),
                     n_reference=len(ctx.bouton_parts['strength/rate']), n_odors=ctx.n, measured=m['view_strength'],
                     part_rate=m['part_rate'], part_mean=m['part_mean'], lo=iv[0], hi=iv[1],
                     draws_median=float(np.median(d))))
    draws['no_step'] = d
    tab = pd.DataFrame(rows).set_index('comparison')
    if abs(tab.loc['no_step', 'measured'] - float(P.no_step_score(ctx.pops['pn_dendrite'],
                                                                 ctx.refs['all']).views['strength'])) > 1e-10:
        raise RuntimeError("the sets' no-step is not E142's")
    ov = {k: overlap([tab.loc['no_step', 'lo'], tab.loc['no_step', 'hi']], [tab.loc[k, 'lo'], tab.loc[k, 'hi']])
          for k in comps}
    # the 12 odors' differences between the sets, beside the step's entries
    o12 = sets['megamat']['odors']
    st = pd.read_csv(remap_path(ctx.lock['inputs']['results']['e142_flies_steps_summary']['rel']))
    st = st[st.step == 'pn_dendrite->pn_bouton'].set_index(['part', 'odor'])
    od = []
    for p in STRENGTH_PARTS:
        k = PART_KEY[p]
        mm = np.nanmean(sets['megamat']['parts'][p], axis=0)
        cc = np.nanmean(sets['control']['parts'][p], axis=0)
        for i, o in enumerate(o12):
            ch = float(st.loc[(k, o), 'change_all'])
            od.append(dict(part=k, odor=o, megamat=float(mm[i]), control=float(cc[i]), sets_diff=float(mm[i] - cc[i]),
                           sets_diff_control_units=float((mm[i] - cc[i]) / units['control'][p]),
                           sets_diff_megamat_units=float((mm[i] - cc[i]) / units['megamat'][p]),
                           step_change=ch, step_change_bouton_units=ch / ctx.units[p],
                           step_both_halves=int(st.loc[(k, o), 'outside_both_halves'])))
    odf = pd.DataFrame(od)
    rms = []
    for k, x in odf.groupby('part', sort=False):
        all17 = st.xs(k, level='part').change_all.to_numpy(dtype=float)
        rms.append(dict(part=k, sets_diff_rms=float(np.sqrt(np.mean(x.sets_diff ** 2))),
                        step_rms_12=float(np.sqrt(np.mean(x.step_change ** 2))),
                        step_rms_17=float(np.sqrt(np.mean(all17 ** 2))),
                        r_sets_diff_vs_step_12=float(np.corrcoef(x.sets_diff, x.step_change)[0, 1])))
    out = dict(intervals={k: [float(tab.loc[k, 'lo']), float(tab.loc[k, 'hi'])] for k in tab.index},
               measured={k: float(tab.loc[k, 'measured']) for k in tab.index}, overlap=ov,
               not_told_apart=bool(any(ov.values())), odor_rms=rms, units=units, q=list(SETS_Q),
               n_boot=N_BOOT_SETS, seed=SETS_SEED, time=_now(), host=_host())
    _write_csv(tab.reset_index(), 'sets.csv')
    _write_csv(pd.DataFrame({k: v for k, v in draws.items()}).rename_axis('draw').reset_index(), 'sets_draws.csv')
    _write_csv(odf, 'sets_odors.csv')
    _write_csv(pd.DataFrame(rms), 'sets_odors_rms.csv')
    _write_json('sets.json', out)
    provenance('sets', ctx, seconds=round(time.time() - t0, 1))
    print(tab.to_string(), '\n', ov, flush=True)
    return out


def model_stage_reading(ctx: Ctx, group: str, penalty: float) -> List[dict]:
    """Beside: the method on a saved group's own stages (E142's group `group`): the wirings' dendrite populations as
    the dendrite flies, their bouton populations' mean entries the target in the bouton flies' units; the model's
    no-step, the left-out weight reading and the in-sample one."""
    from experiments import E142_dendrite_bouton_step as E142
    from experiments import E142_step0 as S0
    g = next(x for x in E142.MODEL_GROUPS if x.name == group)
    runs, missing = S0.group_runs(g)
    if not runs:
        return [dict(group=group, missing=';'.join(missing))]
    reads = S0.read_group(runs, ctx.s0, 'all')
    dpops = [reads[w].populations[('pn_dendrite', PANEL, 'on')] for w in runs]
    bpops = [reads[w].populations[('pn_bouton', PANEL, 'on')] for w in runs]
    del reads
    gl = sorted(set().union(*[set(map(str, p.values.index)) for p in dpops]))
    prep = prepare(S0.fly_matrices(dpops, gl, ctx.odors))
    bparts = fly_parts_of(bpops, ctx.odors)
    ref = synthetic_reference(bparts, ctx.units, ctx.odors)
    ns = score_entries({p: v for p, (v, _) in weighted_entries(np.ones(len(gl)), prep, ctx.all_idx).items()},
                       ref)['view_strength']
    pred, W, recs = left_out(prep, bparts, ctx.units, penalty)
    lo = score_entries(pred, ref)['view_strength']
    w_all, _ = fit_weights(prep, relative_parts(bparts, ctx.all_idx), ctx.units, ctx.all_idx, penalty)
    ins = score_entries(plant_truth(w_all, prep), ref)['view_strength']
    return [dict(group=group, n_wirings=len(runs), wirings=';'.join(runs), missing=';'.join(missing),
                 n_glomeruli=len(gl), penalty=penalty, no_step=ns, left_out=lo, below_no_step_by=ns - lo,
                 in_sample=ins, weights_cv=cv_of(w_all), n_not_converged=int(sum(not r['success'] for r in recs)))]


def part_beside(groups: Sequence[str] = BESIDE_GROUPS) -> dict:
    """Design item 6 (beside_models.csv, beside_claws.csv, beside_claws_features.csv, beside_claws_null.csv,
    beside_flies_seven.csv, beside_counts_rho.csv, beside.json): the method on the model's own stages; the flies' step
    from the boutons to the claws under items 3 and 4; each bouton fly's own step for the seven odors; the fitted
    weights' rank correlation with the connectomes' counts. Classes nothing."""
    from scipy import stats
    from kcq import population_fit as P
    _require_identities()
    t0 = time.time()
    ctx = context()
    lam = _chosen_penalty()
    out = dict(penalty=lam)
    # the seven odors
    st = pd.read_csv(remap_path(ctx.lock['inputs']['results']['e142_flies_steps_summary']['rel']))
    st = st[st.step == 'pn_dendrite->pn_bouton']
    seven = [o for o in ctx.odors if (st[st.odor == o].outside_both_halves != 0).any()]
    out['seven_odors'] = seven
    out['seven_as_locked'] = bool(seven == list(SEVEN_ODORS))
    nsv = ctx.no_step_views()
    rows = []
    for i, f in enumerate(ctx.refs['all'].flies):
        for p in STRENGTH_PARTS:
            for o in SEVEN_ODORS:
                j = ctx.odors.index(o)
                flies_step = float(ctx.refs['all'].mean[p][o] - nsv[p][j])
                own = float(ctx.refs['all'].parts[p].loc[f, o] - nsv[p][j])
                rows.append(dict(fly=f, half=('fit' if f in ctx.names['pn_bouton']['fit'] else 'held'),
                                 part=PART_KEY[p], odor=o, own_step=own, own_step_bouton_units=own / ctx.units[p],
                                 flies_step=flies_step, flies_step_bouton_units=flies_step / ctx.units[p],
                                 same_sign=bool(np.sign(own) == np.sign(flies_step))))
    seven_df = pd.DataFrame(rows)
    out['seven_same_sign'] = {f'{k[0]}|{k[1]}': int(x.same_sign.sum())
                              for k, x in seven_df.groupby(['part', 'odor'], sort=False)}
    # the fitted weights' rank correlation with the counts
    bw = pd.read_csv(out_dir() / 'bound_weights.csv').set_index('glomerulus').reindex(ctx.glomeruli)
    cnt = ctx.counts()
    rho = []
    folds = [c for c in bw.columns if c.startswith('left_out_')]
    for cal in CALYCES:
        r_all = stats.spearmanr(bw.in_sample_all, cnt[cal]).statistic
        r_folds = [stats.spearmanr(bw[c], cnt[cal]).statistic for c in folds]
        rho.append(dict(calyx=cal, rho_in_sample_all=float(r_all), rho_left_out_median=float(np.median(r_folds)),
                        rho_left_out_min=float(np.min(r_folds)), rho_left_out_max=float(np.max(r_folds)),
                        rho_fit_half=float(stats.spearmanr(bw.fit_half, cnt[cal]).statistic)))
    # the step from the boutons to the claws
    crefs = ctx.claw_refs()
    cunits = {p: float(crefs['all'].one_fly[p]) for p in STRENGTH_PARTS}
    cparts = {p: crefs['all'].parts[p].reindex(columns=ctx.odors).to_numpy(dtype=float) for p in STRENGTH_PARTS}
    ns_d = {s: float(P.no_step_score(ctx.pops['pn_dendrite'], crefs[s]).views['strength']) for s in FLY_SETS}
    ns_b = {s: float(P.no_step_score(ctx.pops['pn_bouton'], crefs[s]).views['strength']) for s in FLY_SETS}
    crec, cnull, _, _ = bound_reading(ctx, ctx.prep, cparts, cunits, crefs, lam, ns_d)
    claws = [dict(reading='weights_left_out', fly_set=s, error=crec['left_out'][s], no_step_from_dendrites=ns_d[s],
                  no_step_from_boutons=ns_b[s], null_p05=crec['null'][s]['null_p05']) for s in FLY_SETS]
    rd, fnl, _ = feature_readings(bouton_step_inputs(ctx), cparts, crefs, ns_b)
    out['claws'] = dict(weights=crec, features_null_best_p05=float(np.percentile(fnl.best, NULL_PERCENTILE)),
                        no_step_from_dendrites=ns_d, no_step_from_boutons=ns_b)
    # the model's own stages
    models = []
    for g in groups:
        print(f'beside: the model group {g}', flush=True)
        models += model_stage_reading(ctx, g, lam)
    _write_csv(pd.DataFrame(models), 'beside_models.csv')
    _write_csv(pd.DataFrame(claws), 'beside_claws.csv')
    _write_csv(rd, 'beside_claws_features.csv')
    _write_csv(fnl, 'beside_claws_null.csv')
    _write_csv(cnull, 'beside_claws_weights_null.csv')
    _write_csv(seven_df, 'beside_flies_seven.csv')
    _write_csv(pd.DataFrame(rho), 'beside_counts_rho.csv')
    out.update(models=models, counts_rho=rho, time=_now(), host=_host())
    _write_json('beside.json', out)
    provenance('beside', ctx, seconds=round(time.time() - t0, 1), penalty=lam)
    print(pd.DataFrame(models).to_string(index=False), '\n', pd.DataFrame(rho).to_string(index=False), flush=True)
    return out


def _f(x, spec: str = '.3f') -> str:
    try:
        return format(float(x), spec) if np.isfinite(float(x)) else '–'
    except (TypeError, ValueError):
        return str(x)


def part_classes() -> dict:
    """The classes and the outcome letter (classes.json, tables.md)."""
    lock = require_lock()
    pw, bd, ft, sts = (_read_json(n) for n in ('power.json', 'bound.json', 'features.json', 'sets.json'))
    wc = weight_class(bd, bool(pw['holds']))
    fc = feature_class(ft)
    letters = outcome(bool(pw['holds']), wc['holds'], fc['holds'], bool(sts['not_told_apart']))
    out = dict(power=dict(holds=bool(pw['holds']), penalty=pw['penalty'], at_chosen=pw['at_chosen']),
               weight_makes_step=wc, feature_makes_step=fc,
               not_told_apart=dict(holds=bool(sts['not_told_apart']), intervals=sts['intervals'],
                                   overlap=sts['overlap']),
               outcome=letters, outcome_text={k: OUTCOMES[k] for k in letters}, outcome_rule=OUTCOME_RULE,
               n_forms_tried=dict(weights=1, features=len(FEATURES), penalties_in_power_check=len(PENALTIES)),
               lock_sha256=lock['sha256'], time=_now(), host=_host())
    _write_json('classes.json', out)
    write_tables(out)
    provenance('classes', lock=lock)
    print(json.dumps(dict(outcome=letters, power=pw['holds'], weight=wc['cls'], feature=fc['cls'],
                          not_told_apart=sts['not_told_apart'])), flush=True)
    return out


def write_tables(cls: Mapping) -> None:
    """tables.md: the identities, the power check, the weight bound, the features, the two sets, the beside items and
    the classes."""
    d = out_dir()
    L = ["# E146: what could make the flies' step from the PN dendrites to the PN boutons (no simulation)", '',
         "Errors are the strength view's (the rate and mean parts, root mean square) in the bouton flies' one-fly units "
         "unless a row says otherwise; 'fit' and 'held' are E136's halves of the bouton flies (2 / 3), each half's mean "
         "in all five flies' units. The no-step is the seven dendrite flies' mean entries as the prediction.", '']
    idf = pd.read_csv(d / 'identities.csv')
    L += ['## The identities', '', '| value | quoted | reproduced | saved table | reproduces |', '|---|---:|---:|---:|:---:|']
    L += [f'| {r.name} | {r.quoted:g} | {_f(r.reproduced, ".4f")} | {_f(r.saved, ".4f")} | {"yes" if r.reproduces else "no"} |'
          for r in idf.itertuples()]
    sm = pd.read_csv(d / 'power_summary.csv')
    pw = json.loads((d / 'power.json').read_text())
    L += ['', f"## The power check (chosen penalty {pw['penalty']:g}; the bound has power: {pw['holds']})", '',
          'Per penalty and kind: plants, eligible (planted no-step >= 0.9), passing (left-out reading below the planted '
          'no-step by 0.10), share; medians of the planted no-step and the left-out reading.', '',
          '| reading | penalty | kind | plants | eligible | passing | share | no-step median | left-out median | '
          'class holds |', '|---|---:|---|---:|---:|---:|---:|---:|---:|:---:|']
    for r in sm.itertuples():
        L.append(f'| {r.reading} | {r.penalty:g} | {r.kind} | {r.n} | {r.n_eligible} | {r.n_pass} | {_f(r.share, ".2f")} '
                 f'| {_f(r.no_step_median)} | {_f(r.left_out_median)} | {r.class_holds} |')
    un = pd.read_csv(d / 'power_unevenness.csv')
    L += ['', "The uneven draws' coefficient of variation per penalty (the weights fitted to all 17 odors): "
          + ', '.join(f'{r.penalty:g}: {r.fitted_weights_cv:.2f}' for r in un.itertuples())]
    bd = pd.read_csv(d / 'bound.csv')
    L += ['', '## The weight bound', '', '| reading | fly set | error | no-step | below the no-step by | permuted 5th pct | '
          'percentile of the actual |', '|---|---|---:|---:|---:|---:|---:|']
    for r in bd.itertuples():
        L.append(f'| {r.reading} | {r.fly_set} | {_f(r.error)} | {_f(r.no_step)} | {_f(r.below_no_step_by, "+.3f")} | '
                 f'{_f(getattr(r, "null_p05", np.nan))} | {_f(getattr(r, "null_percentile_of_actual", np.nan), ".1f")} |')
    bp = pd.read_csv(d / 'bound_by_penalty.csv')
    L += ['', 'Beside, classing nothing, the left-out reading at every penalty (all / fit / held): ' + '; '.join(
        f'{lam:g}: ' + ' / '.join(_f(x[x.fly_set == s].left_out.iloc[0]) for s in FLY_SETS)
        for lam, x in bp.groupby('penalty', sort=True))]
    fd = pd.read_csv(d / 'features.csv')
    ftj = json.loads((d / 'features.json').read_text())
    L += ['', f"## The population features (permuted best's 5th percentile {_f(ftj['null_best_p05'])})", '',
          '| form | all | fit | held | no-step all / fit / held |', '|---|---:|---:|---:|---|']
    for f, x in fd.groupby('form', sort=False):
        o = x.set_index('fly_set')
        L.append(f'| {f} | {_f(o.loc["all", "left_out"])} | {_f(o.loc["fit", "left_out"])} | '
                 f'{_f(o.loc["held", "left_out"])} | ' + ' / '.join(_f(o.loc[s, 'no_step']) for s in FLY_SETS) + ' |')
    ss = pd.read_csv(d / 'sets.csv')
    L += ['', '## Two sets at one stage (5-95 % intervals over 1,000 resamples of both sets)', '',
          '| comparison | measured | interval | resamples median |', '|---|---:|---|---:|']
    L += [f'| {r.comparison} | {_f(r.measured)} | {_f(r.lo)} to {_f(r.hi)} | {_f(r.draws_median)} |'
          for r in ss.itertuples()]
    rr = pd.read_csv(d / 'sets_odors_rms.csv')
    L += ['', "The 12 odors' differences between the sets beside the step (root mean square over odors, entries): "
          + '; '.join(f'{r.part}: sets {_f(r.sets_diff_rms)}, step on the 12 {_f(r.step_rms_12)} (on 17 '
                      f'{_f(r.step_rms_17)}), r {_f(r.r_sets_diff_vs_step_12, ".2f")}' for r in rr.itertuples())]
    if (d / 'beside_models.csv').exists():
        bm = pd.read_csv(d / 'beside_models.csv')
        L += ['', "## Beside: the method on the model's own stages", '',
              '| group | wirings | glomeruli | no-step | left-out | in sample |', '|---|---:|---:|---:|---:|---:|']
        L += [f'| {r.group} | {getattr(r, "n_wirings", "")} | {getattr(r, "n_glomeruli", "")} | '
              f'{_f(getattr(r, "no_step", np.nan))} | {_f(getattr(r, "left_out", np.nan))} | '
              f'{_f(getattr(r, "in_sample", np.nan))} |' for r in bm.itertuples()]
        bc = pd.read_csv(d / 'beside_claws.csv')
        bcf = pd.read_csv(d / 'beside_claws_features.csv')
        L += ['', "## Beside: the step from the boutons to the claws (claw flies' units)", '',
              '| reading | all | fit | held |', '|---|---:|---:|---:|']
        L.append('| weights, left out (dendrite glomeruli) | ' + ' | '.join(
            _f(bc[bc.fly_set == s].error.iloc[0]) for s in FLY_SETS) + ' |')
        L.append('| no-step from the dendrites | ' + ' | '.join(
            _f(bc[bc.fly_set == s].no_step_from_dendrites.iloc[0]) for s in FLY_SETS) + ' |')
        L.append('| no-step from the boutons | ' + ' | '.join(
            _f(bc[bc.fly_set == s].no_step_from_boutons.iloc[0]) for s in FLY_SETS) + ' |')
        for f, x in bcf.groupby('form', sort=False):
            o = x.set_index('fly_set')
            L.append(f'| feature {f} (boutons upstream) | ' + ' | '.join(_f(o.loc[s, 'left_out']) for s in FLY_SETS)
                     + ' |')
        cr = pd.read_csv(d / 'beside_counts_rho.csv')
        L += ['', "## Beside: the fitted weights' Spearman ρ with the bouton counts", '',
              '| calyx | in sample, all odors | left-out folds: median (min to max) | fit half |', '|---|---:|---:|---:|']
        L += [f'| {r.calyx} | {_f(r.rho_in_sample_all, ".2f")} | {_f(r.rho_left_out_median, ".2f")} '
              f'({_f(r.rho_left_out_min, ".2f")} to {_f(r.rho_left_out_max, ".2f")}) | {_f(r.rho_fit_half, ".2f")} |'
              for r in cr.itertuples()]
        sv = pd.read_csv(d / 'beside_flies_seven.csv')
        L += ['', "## Beside: each bouton fly's own step for the seven odors (bouton flies' units)", '',
              '| part | odor | flies\' step | ' + ' | '.join(sv.fly.unique()) + ' | same sign |',
              '|---|---|---:|' + '---:|' * sv.fly.nunique() + '---:|']
        for (p, o), x in sv.groupby(['part', 'odor'], sort=False):
            L.append(f'| {p} | {o} | {_f(x.flies_step_bouton_units.iloc[0], "+.2f")} | '
                     + ' | '.join(_f(v, '+.2f') for v in x.own_step_bouton_units) + f' | {int(x.same_sign.sum())} of {len(x)} |')
    L += ['', '## Classes', '', '```', json.dumps(cls, indent=1, default=_e140()._json_default), '```', '']
    (d / 'tables.md').write_text('\n'.join(L) + '\n')


PARTS: Tuple[str, ...] = ('identities', 'power', 'bound', 'features', 'sets', 'beside', 'classes')


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', nargs='+', required=True, choices=('lock', 'amend') + PARTS)
    ap.add_argument('--workers', type=int, default=2, help="the power check's processes (default 2)")
    a = ap.parse_args(argv)
    if 'lock' in a.part:
        part_lock()
        return
    if 'amend' in a.part:
        part_amend()
        return
    fns = dict(identities=part_identities, power=lambda: part_power(a.workers), bound=part_bound,
               features=part_features, sets=part_sets, beside=part_beside, classes=part_classes)
    for p in PARTS:
        if p in a.part:
            fns[p]()


if __name__ == '__main__':
    main()
