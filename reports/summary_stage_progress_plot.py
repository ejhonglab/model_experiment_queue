#!/usr/bin/env python3
"""The Summary's stage-target progress plot: the best model so far on each statistic of the stage target, against
experiment completion time (CLAUDE.md, "Summary reports", "The stage-target progress plot waits for the first fit to a
stage"; its reading as ruled on 2026-11-25 by a Fable subagent, claude-fable-5-2, spawned from hal 11:35).

Series, in one-fly units (lower is better; 1: as far from the flies' mean as one fly is from the others): the megamat
stage error on strength by odor and each megamat stage's strength error (PN dendrites, PN boutons, KC claws, KC
somata); the KC somata's strength error on validation2, on the control mixtures as uploaded and on the APL16 panel with
APL on and with APL off. As descriptions, in a second panel: the APL-effect loss and its pattern part under the
calibrated readout (the APL16 panel, a held-out test).

One reading per series. A series' level is the unrefitted best family's value, the mean over its eight wirings against
all the flies (`levels`). A variant run on all eight wirings and fitted to no flies is drawn at its own value; a fitted
variant at the level plus its change from the family where it was not fitted (the T wirings and, at a fitted stage, the
held-out half of the flies, with the family read there too: the change rule 3 gates on, from its experiment's
classes). No point for a variant outside the realism bands, one read on fewer than four wirings, one with a response
rate not held at its target on any wiring, a control fit or a fit to permuted flies. Hollow points (shown; no
staircase steps on them; ruled 2026-11-26 by a Fable subagent, claude-fable-5-2, spawned from viki 19:30): every
variant but the family's own readings and those that their experiment's score classes as holding (rule 3: megamat
only, a fall at another stage's or view's cost, does not improve, no class, and an experiment scored without these
classes), a value that holds a stage the variant was fitted to or set from, on flies not held out (the PN dendrites of
the antennal-lobe fits and of the imposed recordings, and a megamat stage error that includes them), and a variant on
another brain's wiring (the male CNS's calyces), drawn from its calyx's own family's value against all the flies plus
its change from that family. An experiment's point is its lowest filled variant, else its lowest hollow one; each
series' staircase is the running best of its filled points.

The experiments read (`ADAPTERS`; one function each, from the experiment's own committed tables):
- E134, the models as populations at every recorded stage: the best family as fitted and held at the megamat flies'
  response rate (score_stage_table.csv, each stage's default reading; writeup_panels.csv; writeup_apl16_models.csv);
  its beside models are read on one wiring each (no point);
- E135, the antennal lobe held to the recorded PN dendrites: its arms on eight wirings, but the permuted controls
  (score_scores.csv, score_apl_effect.csv, score_realism.csv);
- E136, the best family fitted to all megamat stages at once: its arms but the permuted control (classes.json's changes
  on T against the held-out half, score_apl16.csv, in its measure's folder);
- E139, the weaker antennal-lobe normalization without the calyx refit: F at the four rungs, fitted to no flies, on
  all eight wirings against all the flies (score_megamat.csv's 'all' half, score_panels.csv, score_apl16.csv); B, K
  and step 3's refit A0 at 314 as fitted (classes.json, T against the held-out half); A1 beside is E136's;
- E141, the rise of the inhibition among KCs: its arms, all on E139's B at 1000, as fitted (classes.json); S outside
  a realism band and P and T, whose reads stopped after their holds failed, without a point; B, K and A1 beside are
  E139's and E136's;
- E142, the step from the PN dendrites to the PN boutons: the amended main arm, the arm with APL→bouton ×0.155 alone
  and the lock's main arm (no class: its holds fail), as fitted, from the step-1 score's tables on T against the
  held-out half (score_step1_megamat.csv, score_step1_panels.csv, score_step1_apl16.csv; step1_classes.json's classes);
  their dendrites and stage errors hollow (imposed); no point for the permuted controls, the arms stopped by claim and
  E135's arms re-scored beside;
- E144, A1 and its parts on the male CNS's two calyces: the settings read on a calyx's four wirings, each against its
  calyx's own family on the same wirings (classes.json, score_apl16.csv), drawn from that family's value against all
  the flies (score_megamat.csv, score_panels.csv) plus that change, all hollow (a variant on another brain's wiring);
  the settings read on one wiring (all of H's, F, Dabs) and those without a class (L's D5.8) without a point;
- E145, the stage fit made again in three other cells of wirings by halves of the flies: X, F and W, each read on its
  four read wirings against its read half (classes.json's cells; the APL effect from score_apl16.csv), all hollow by
  their classes (X none, at KC–KC's bound; F one stage at another's cost; W strength at another view's cost); its
  cell A1 beside is E136's arm. E145's four-factor control is a yardstick, with no point.
A variant at fixed scales that carries a fitted arm's parameters (E139's B and K, E141's arms, E142's arms, E144's
settings: A1's APL→bouton ×0.155 or KC–KC rise, fitted on S against the fit half) is read as a fitted one, where those
were not fitted: on wirings outside S and, at the megamat stages, against the held-out half; each experiment's own
classes are read there. E144's change is from each male calyx's own family (E130's fit of the family on that calyx,
as each of the hemibrain's wirings has its own fit of the family), added to that family's own value against all the
flies, and its points are hollow (ruled 2026-11-26 by a Fable subagent, claude-fable-5-2, spawned from viki 19:31).
Not read (nothing of their own to draw): E137, the held-out panels read like for like (new readings of the flies,
with the family and yardsticks re-read); E138, the odors that carry the strength error (no simulation, saved runs'
entries); E140, aa at the step from the claws to the somata (no simulation; its step 1 did not run); E143, the megamat
stage error through the flies' sampling (no simulation; resamples of the flies under saved runs).
An experiment counts once it is completed in the experiment index (`--include E<n>` dates a running one now, for a
preview). Writes, under <shared>/results/summary_progress/:

    stage_progress_points.csv   one row per experiment, variant and series (drawn or not, and why)
    stage_progress_best.csv     the running best per series
    stage_progress_plot.html    the plot (inline SVG + script), inserted by md_to_artifact.py's `%%% html <path>` line

Usage (reports/build_summary.py runs it; `--skip-stage-progress-plot` there skips it):
    ./venv/bin/python reports/summary_stage_progress_plot.py [--include E136]
"""
from __future__ import annotations

import argparse
import html
import json
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from kcq import harness as H

OUT = H.shared_dir() / 'results' / 'summary_progress'
RESULTS = H.shared_dir() / 'results'
STAGES = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')
T_WIRINGS = ('P3', 'P4', 'M2', 'M3')
MIN_WIRINGS = 4               # a variant read on fewer wirings has no point
READING = 'curve_noise'       # every stage's default reading (kcq.stage_readings.DEFAULT_READINGS)
IMPROVE_SE = 2.0              # E136's bar: a fall of this many standard errors over wirings (its IMPROVE_SE)

# series: key -> (label, panel); panel 'stage' (one-fly units) or 'apl' (the APL effect, a description)
SERIES: Dict[str, tuple] = {
    'stage_error': ('the megamat stage error, strength by odor', 'stage'),
    'pn_dendrite': ('PN dendrites, strength by odor', 'stage'),
    'pn_bouton': ('PN boutons, strength by odor', 'stage'),
    'kc_claw': ('KC claws, strength by odor', 'stage'),
    'kc_soma': ('KC somata, strength by odor', 'stage'),
    'validation2': ('KC somata on validation2', 'stage'),
    'cmix': ('KC somata on the control mixtures, as uploaded', 'stage'),
    'apl16_on': ('KC somata on the APL16 panel, APL on', 'stage'),
    'apl16_off': ('KC somata on the APL16 panel, APL off', 'stage'),
    'apl_loss': ('the APL-effect loss, calibrated readout', 'apl'),
    'apl_pattern': ('its pattern part', 'apl'),
}
SERIES_COLORS = {
    'stage_error': ('#1f2328', '#e6edf3'),
    'pn_dendrite': ('#2a78d6', '#3987e5'), 'pn_bouton': ('#83b1ed', '#8cbaf7'),
    'kc_claw': ('#1baf7a', '#199e70'), 'kc_soma': ('#008300', '#2a9d2a'),
    'validation2': ('#eda100', '#c98500'), 'cmix': ('#e87ba4', '#d55181'),
    'apl16_on': ('#bc4c00', '#e0823d'), 'apl16_off': ('#f0a875', '#f4b88f'),
    'apl_loss': ('#8250df', '#a475f9'), 'apl_pattern': ('#c2a5f0', '#c9b0fb'),
}
Y_MAX = {'stage': 2.0, 'apl': 0.12}


def _mean(x) -> float:
    x = pd.Series(list(x), dtype=float).dropna()
    return float(x.mean()) if len(x) else np.nan


def _per_wiring(rows: pd.DataFrame, col: str = 'view_strength') -> pd.DataFrame:
    """From one variant's rows (stage, wiring, `col`) at the megamat stages: per wiring each stage's value and the
    megamat stage error (the root mean square over the four stages); empty without all four."""
    p = rows.pivot_table(index='wiring', columns='stage', values=col, aggfunc='first')
    if not set(STAGES) <= set(p.columns):
        return pd.DataFrame()
    p = p[list(STAGES)].dropna()
    p['stage_error'] = np.sqrt((p[list(STAGES)] ** 2).mean(axis=1))
    return p


def _point(exp, variant, series, value, reading, drawn=True, filled=True, why='', **popup) -> dict:
    """One variant's value on one series: `value` as drawn, `drawn` False for no point, `filled` False for a hollow one
    (`why` says which and why); `popup` the fitted variant's details (arm and family readings, change, se, bar, class,
    the control's fall)."""
    return dict(exp=exp, variant=variant, series=series, value=value, reading=reading, drawn=bool(drawn),
                filled=bool(filled and drawn), why=why, popup=json.dumps(popup) if popup else '')


HOLDS = 'holds'
RULE3 = ('not a hold (rule 3: a fall that is megamat only, or bought at another stage or view, is no improvement; '
         'ruled 2026-11-26 by a Fable subagent, claude-fable-5-2, spawned from viki 19:30)')
NOT_HELD = 'a response rate not held at its target (not at the operating point)'


def _class_hollow(cls: Optional[str], reasons: Sequence[str] = ()) -> Optional[str]:
    """Why a variant's points are hollow by its class: None for a variant its experiment classes as holding (rule 3),
    else its other class or its lack of one (a staircase steps only on the family's readings and on a hold)."""
    if cls == HOLDS:
        return None
    if not cls or str(cls).startswith('no class'):
        return 'no class' + (': ' + '; '.join(reasons) if reasons else '')
    return f'{cls}: {RULE3}'


def _hold_fails(reasons: Sequence[str]) -> List[str]:
    """The reasons that say a hold failed: a variant with any has no point (NOT_HELD)."""
    return [str(r) for r in reasons or [] if str(r).startswith('a hold fails')]


def _e134_dir() -> Path:
    return RESULTS / 'E134_population_baseline'


def _e136_dir() -> Path:
    return RESULTS / 'E136_stage_population_fit' / 'measure_strength'


def levels() -> Dict[str, float]:
    """Each series' level: the unrefitted best family held at the megamat flies' response rate (E134's held runs), the
    mean over its eight wirings against all the flies; the APL16 panel from its own fits (E134's writeup_apl16_models.csv,
    curve and noise), the APL effect from E136's score of the same runs (score_apl16.csv) when present."""
    d = _e134_dir()
    t = pd.read_csv(d / 'score_stage_table.csv')
    t = t[(t.group == 'megamat_held') & (t.reading == READING) & (t.condition == 'on') & (t.panel == 'megamat')
          & t.wiring.isin(WIRINGS)]
    p = _per_wiring(t)
    out = {k: float(p[k].mean()) for k in p.columns}
    pan = pd.read_csv(d / 'writeup_panels.csv').set_index('panel')
    for k in ('validation2', 'cmix'):
        if k in pan.index:
            out[k] = float(pan.loc[k, 'held_strength'])
    a = pd.read_csv(d / 'writeup_apl16_models.csv')
    for c in ('on', 'off'):
        x = a[(a.model == 'the best family, eight wirings') & (a.condition == c)]
        if len(x):
            out[f'apl16_{c}'] = float(x.curve_noise_strength.iloc[0])
    pa = _e136_dir() / 'score_apl16.csv'
    if pa.exists():
        f = pd.read_csv(pa)
        f = f[f.arm == 'family']
        if len(f):
            out.update(apl_loss=_mean(f.cal_loss), apl_pattern=_mean(f.cal_pattern))
    return out


def e134_points() -> List[dict]:
    """E134's variants: the best family held at the megamat flies' response rate (the level itself) and as fitted, each
    on its eight wirings and fitted to no flies, so drawn at its own value; the beside models, read on one wiring each,
    have no point."""
    d = _e134_dir()
    t = pd.read_csv(d / 'score_stage_table.csv')
    t = t[(t.reading == READING) & (t.condition == 'on') & (t.panel == 'megamat') & t.stage.isin(STAGES)
          & ~t.wiring.isin(['mean', 'common'])]
    pan = pd.read_csv(d / 'writeup_panels.csv').set_index('panel')
    apl = pd.read_csv(d / 'writeup_apl16_models.csv')
    lv = levels()
    out = []
    groups = {'megamat_held': ("the best family, held at the megamat flies' response rate", 'held_strength'),
              'megamat_fitted': ('the best family, as fitted', 'fitted_strength')}
    for g, (label, pcol) in groups.items():
        p = _per_wiring(t[t.group == g])
        reading = f'{len(p)} wirings against all flies'
        for k in p.columns:
            out.append(_point('E134', label, k, float(p[k].mean()), reading))
        for k in ('validation2', 'cmix'):
            if k in pan.index:
                out.append(_point('E134', label, k, float(pan.loc[k, pcol]), 'eight wirings against all flies'))
        for c in ('on', 'off'):
            x = apl[(apl.model == 'the best family, eight wirings') & (apl.condition == c)]
            if len(x):
                out.append(_point('E134', label, f'apl16_{c}', float(x.curve_noise_strength.iloc[0]),
                                  "the APL16 panel: the family's own fits on eight wirings, against all recordings"))
        for k in ('apl_loss', 'apl_pattern'):
            if k in lv:
                out.append(_point('E134', label, k, lv[k], "the family's own APL16 fits on eight wirings (read in "
                                  "E136's score)"))
    beside = sorted(set(t.group) - set(groups))
    for g in beside:
        out.append(_point('E134', g.replace('_megamat', '').replace('_', ' '), 'stage_error', np.nan, 'one wiring',
                          drawn=False, why='read on fewer than four wirings'))
    return out


E135_ARMS = {'imposed': 'the recorded dendrites imposed', 'imposed_dff': 'the recorded dendrites imposed, table in dF/F',
             'R1': "R1 (the PN stage's global parameters fitted to the dendrites)",
             'rung': 'R3 (the highest kept rung)', 'imposed_route2_repeat': 'route 2 (the solved input), real wiring'}
# E135's dendrites are fitted to or set from all the PN flies, so its dendrite values and the stage errors that hold
# them are hollow
E135_HOLLOW = ('pn_dendrite', 'stage_error')
# and E135 was scored without rule 3's classes, so every point of it is hollow (ruled 2026-11-26)
E135_NO_CLASS = 'no class: scored without the classes of rule 3'


def e135_points() -> List[dict]:
    """E135's arms but the permuted controls: the megamat stages and the APL16 panel (score_scores.csv, the default
    reading), the calibrated APL effect (score_apl_effect.csv), on eight wirings and fitted to no KC or bouton flies, so
    drawn at their own values; every point hollow (scored without rule 3's classes; the dendrites and the stage errors
    that hold them also set from the PN flies, E135_HOLLOW); no point for an arm on fewer than four wirings or outside
    the realism bands (score_realism.csv)."""
    d = RESULTS / 'E135_antennal_lobe_dendrites'
    s = pd.read_csv(d / 'score_scores.csv')
    s = s[(s.reading == READING) & ~s.wiring.isin(['mean', 'common'])]
    real = pd.read_csv(d / 'score_realism.csv')
    ae = pd.read_csv(d / 'score_apl_effect.csv')
    out = []
    for arm, label in E135_ARMS.items():
        bad = [f for f in real[real.arm == arm].realism_flag if isinstance(f, str) and f.strip()]
        p = _per_wiring(s[(s.group == f'megamat_{arm}') & (s.condition == 'on') & (s.panel == 'megamat')
                          & s.stage.isin(STAGES)])
        drawn, why = True, ''
        if bad:
            drawn, why = False, 'outside the realism bands: ' + bad[0]
        elif len(p) < MIN_WIRINGS:
            drawn, why = False, f'read on {len(p)} wiring{"s" if len(p) != 1 else ""}'
        reading = f'{len(p)} wiring{"s" if len(p) != 1 else ""} against all flies'
        for k in p.columns:
            hw = ('the PN dendrites are fitted to or set from all the PN flies; ' if k in E135_HOLLOW else '') + \
                E135_NO_CLASS
            out.append(_point('E135', label, k, float(p[k].mean()), reading, drawn, False, why or hw))
        a = s[(s.group == f'apl16_{arm}') & (s.stage == 'kc_soma') & (s.panel == 'apl16')]
        for c in ('on', 'off'):
            x = a[a.condition == c]
            if len(x):
                out.append(_point('E135', label, f'apl16_{c}', _mean(x.view_strength),
                                  f'the APL16 panel, {x.wiring.nunique()} wirings, against all recordings',
                                  drawn and x.wiring.nunique() >= MIN_WIRINGS, False, why or E135_NO_CLASS))
        e = ae[ae.arm == arm]
        if len(e):
            for k, col in (('apl_loss', 'cal_apl_effect_rmse'), ('apl_pattern', 'cal_apl_effect_rmse_centered')):
                out.append(_point('E135', label, k, _mean(e[col]), f'{e.wiring.nunique()} wirings',
                                  drawn and e.wiring.nunique() >= MIN_WIRINGS, False, why or E135_NO_CLASS))
    return out


E136_ARMS = {'A0': 'A0, the family refitted', 'A1': 'A1, with the PN normalization searched',
             'A2': "A2, APL's inhibition of the claws as a gain", 'A3': 'A3, the per-claw threshold free',
             'A4': 'A4, without the inhibition among KCs'}
E136_PANELS = {'validation2': 'validation2/validation2/on', 'cmix': 'mixtures/cmix/on', 'apl16_on': 'apl16/apl16/on',
               'apl16_off': 'apl16/apl16/off'}


def e136_points() -> List[dict]:
    """E136's arms but the permuted control (fitted): each drawn at the level plus its change from the unrefitted family
    on T (classes.json: against the held-out half at the megamat stages, against all the panel's flies on a held-out
    panel; the APL effect from score_apl16.csv, the mean over T of the arm minus the family); hollow unless classed
    as holding (rule 3; `_class_hollow`); no point outside the realism bands (any wiring's flag in score_apl16.csv)."""
    d = _e136_dir()
    cl = json.loads((d / 'classes.json').read_text()) if (d / 'classes.json').exists() else {}
    if not cl.get('classes'):
        return []
    a = pd.read_csv(d / 'score_apl16.csv') if (d / 'score_apl16.csv').exists() else pd.DataFrame()
    lv = levels()
    cf = cl.get('control_fall')
    out = []
    for arm, label in E136_ARMS.items():
        r = cl['classes'].get(arm)
        if not r or r.get('cls') == 'reads missing':
            continue
        cls = r.get('cls')
        flags = [f for f in (a[a.arm == arm].realism_flag if len(a) else []) if isinstance(f, str) and f.strip()]
        drawn, why = (False, 'outside the realism bands: ' + flags[0]) if flags else (True, '')
        hw = _class_hollow(cls, r.get('reasons') or [])
        filled = hw is None
        if drawn and not filled:
            why = hw
        comps = {'stage_error': (r.get('megamat') or {}, 'the held-out wirings (T) against the held-out half'),
                 **{s: ((r.get('stages') or {}).get(s) or {}, 'the held-out wirings (T) against the held-out half')
                    for s in STAGES},
                 **{k: ((r.get('panels') or {}).get(key) or {}, "the held-out wirings (T) against all the panel's flies")
                    for k, key in E136_PANELS.items()}}
        for k, (c, where) in comps.items():
            if not c.get('complete') or k not in lv:
                continue
            se = c.get('se')
            out.append(_point('E136', label, k, lv[k] + float(c['change']), where, drawn, filled, why,
                              arm=c.get('arm'), family=c.get('ref'), change=c.get('change'), se=se,
                              bar=-IMPROVE_SE * se if se is not None else None, cls=cls, control_fall=cf,
                              where=where))
        if len(a):
            fam = a[(a.arm == 'family') & a.wiring.isin(T_WIRINGS)].set_index('wiring')
            x = a[(a.arm == arm) & a.wiring.isin(T_WIRINGS)].set_index('wiring')
            for k, col in (('apl_loss', 'cal_loss'), ('apl_pattern', 'cal_pattern')):
                ws = [w for w in T_WIRINGS if w in x.index and w in fam.index]
                if len(ws) < MIN_WIRINGS or k not in lv:
                    continue
                ch = _mean(x.loc[ws, col] - fam.loc[ws, col])
                out.append(_point('E136', label, k, lv[k] + ch, 'the held-out wirings (T)', drawn, filled, why,
                                  arm=_mean(x.loc[ws, col]), family=_mean(fam.loc[ws, col]), change=ch, cls=cls,
                                  control_fall=cf, where='the held-out wirings (T)'))
    return out


# ---------------------------------------------------------------------------
# the readings shared by the adapters from E139 on: a variant's comparisons with the family (E136's `compare` dicts),
# from its experiment's classes or from its score tables, and its points at the level plus each change

T_HELD = 'the held-out wirings (T) against the held-out half'
T_PANEL = "the held-out wirings (T) against all the panel's flies"
ALL8 = 'all eight wirings against all the flies (fitted to none)'
ALL8_PANEL = "all eight wirings against all the panel's flies (fitted to none)"
# the series of the held-out panels and the APL effect: (panel, condition) of the score_panels tables, and the column
# of the score_apl16 tables
PANEL_SERIES = {'validation2': ('validation2', 'on'), 'cmix': ('cmix', 'on'), 'apl16_on': ('apl16', 'on'),
                'apl16_off': ('apl16', 'off')}
APL_SERIES = {'apl_loss': 'cal_loss', 'apl_pattern': 'cal_pattern'}
# each series' spread in the classes' spreads (the family's SD over its eight wirings); the APL effect has none
SPREAD_KEYS = {'stage_error': 'stage', **{s: s for s in STAGES}, **E136_PANELS}


def _compare(arm_v: Dict[str, float], base_v: Dict[str, float], spread: float, wirings: Sequence[str]) -> dict:
    """E136's `compare` (experiments/E136_stage_population_fit.py) on table values (per wiring): the change over
    `wirings` (the variant minus the family on the same wiring), its standard error over the wirings (the larger of
    half the spread and the changes' SD over sqrt(n)) and the two means; incomplete where a wiring is missing."""
    ch = [arm_v[w] - base_v[w] for w in wirings if w in arm_v and w in base_v
          and np.isfinite(arm_v[w]) and np.isfinite(base_v[w])]
    if len(ch) < len(wirings) or not len(ch):
        return dict(n=len(ch), complete=False)
    sd = float(np.std(ch, ddof=1)) if len(ch) > 1 else np.nan
    se = float(max(spread / 2.0 if np.isfinite(spread) else 0.0, sd / np.sqrt(len(ch))))
    return dict(n=len(ch), complete=True, change=float(np.mean(ch)), se=se, spread=spread,
                arm=float(np.mean([arm_v[w] for w in wirings])), ref=float(np.mean([base_v[w] for w in wirings])))


def _by_wiring(df: pd.DataFrame, key: str, name: str, value: str, **sel) -> Dict[str, float]:
    """{wiring: value} of the rows of `df` whose `key` column is `name` (and each `sel` column its value)."""
    x = df[df[key].astype(str) == name]
    for c, v in sel.items():
        x = x[x[c].astype(str) == str(v)]
    return {str(w): float(v) for w, v in zip(x.wiring, x[value]) if str(w) not in ('mean', 'common', 'mean_S',
                                                                                    'mean_T')}


def _class_comps(r: dict, where: str = T_HELD, panel_where: str = T_PANEL) -> Dict[str, tuple]:
    """A classes record's comparisons (E136's form: 'megamat', 'stages', 'panels'), per series with where it was read."""
    return {'stage_error': (r.get('megamat') or {}, where),
            **{s: ((r.get('stages') or {}).get(s) or {}, where) for s in STAGES},
            **{k: ((r.get('panels') or {}).get(key) or {}, panel_where) for k, key in E136_PANELS.items()}}


def _table_comps(mega: pd.DataFrame, panels: pd.DataFrame, apl: pd.DataFrame, name: str, family: str,
                 wirings: Sequence[str], half: str, spreads: Dict[str, float], where: str, panel_where: str,
                 apl_where: str, key: str = 'variant', panel_key: str = 'arm', apl_key: str = 'arm') -> Dict[str, tuple]:
    """A variant's comparisons with `family` from its experiment's score tables, on `wirings`: the megamat stages and
    their stage error against `half` of the flies (score_megamat: 'stage' and each stage's column), the held-out panels'
    KC somata against all their flies (score_panels' 'error', strength by odor), the calibrated APL effect
    (score_apl16)."""
    out = {}
    if len(mega):
        for k, col in (('stage_error', 'stage'), *((s, s) for s in STAGES)):
            out[k] = (_compare(_by_wiring(mega, key, name, col, half=half),
                               _by_wiring(mega, key, family, col, half=half),
                               float(spreads.get(SPREAD_KEYS[k], np.nan)), wirings), where)
    if len(panels):
        p = panels[panels.stage.astype(str) == 'kc_soma'] if 'stage' in panels else panels
        for k, (pan, cond) in PANEL_SERIES.items():
            out[k] = (_compare(_by_wiring(p, panel_key, name, 'error', panel=pan, condition=cond),
                               _by_wiring(p, panel_key, family, 'error', panel=pan, condition=cond),
                               float(spreads.get(SPREAD_KEYS[k], np.nan)), wirings), panel_where)
    out.update(_apl_comps(apl, name, family, wirings, apl_where, apl_key))
    return out


def _apl_comps(apl: pd.DataFrame, name: str, family: str, wirings: Sequence[str], where: str,
               key: str = 'arm') -> Dict[str, tuple]:
    """The calibrated APL effect's loss and pattern part (score_apl16's cal_loss, cal_pattern): the variant minus the
    family on `wirings`, with no standard error (as E136's adapter reads it)."""
    out = {}
    if not len(apl):
        return out
    for k, col in APL_SERIES.items():
        c = _compare(_by_wiring(apl, key, name, col), _by_wiring(apl, key, family, col), np.nan,
                     wirings)
        out[k] = (dict(c, se=None), where)
    return out


def _flags(apl: pd.DataFrame, name: str, key: str = 'arm') -> List[str]:
    """The realism flags of a variant's APL16 reads (score_apl16's realism_flag), as '<wiring>: <flag>'."""
    if not len(apl) or 'realism_flag' not in apl:
        return []
    x = apl[(apl[key].astype(str) == name) & ~apl.wiring.astype(str).isin(['mean', 'common'])]
    return [f'{w}: {f}' for w, f in zip(x.wiring, x.realism_flag) if isinstance(f, str) and f.strip()]


def _emit(exp: str, label: str, comps: Dict[str, tuple], lv: Dict[str, float], drawn: bool = True, why: str = '',
          hollow: Optional[Dict[str, str]] = None, cls: Optional[str] = None, control_fall: Optional[float] = None,
          fitted: bool = True, level_note: Optional[str] = None) -> List[dict]:
    """A variant's points: each complete comparison drawn at the series' level plus its change (`lv`, `levels`), with
    the popup's readings; `hollow` the series drawn open, by reason ('*' every series); `fitted` False for a variant
    fitted to no flies, read on all eight wirings (no bar: the class is read elsewhere); `level_note` names a level
    that is not the unrefitted family's (E144's calyx families), and the popup then gives it. None complete: one row on
    the stage error, not drawn, saying why."""
    hollow = hollow or {}
    out = []
    for k, (c, where) in comps.items():
        if not c.get('complete') or k not in lv:
            continue
        se = c.get('se')
        hw = hollow.get(k, hollow.get('*'))
        out.append(_point(exp, label, k, lv[k] + float(c['change']), where, drawn, hw is None,
                          why if (why or not drawn) else (hw or ''), arm=c.get('arm'), family=c.get('ref'),
                          change=c.get('change'), se=se, bar=-IMPROVE_SE * se if fitted and se is not None else None,
                          cls=cls, control_fall=control_fall, where=where,
                          **(dict(level=lv[k], level_note=level_note) if level_note else {})))
    if not out:
        out.append(_point(exp, label, 'stage_error', np.nan, '-', drawn=False,
                          why=why or 'no comparison with the family on four wirings or more'))
    return out


def _undrawn(exp: str, label: str, why: str) -> dict:
    """A variant with no point (a control, a reference, a beside arm drawn under its own experiment), on the stage
    error, saying why."""
    return _point(exp, label, 'stage_error', np.nan, '-', drawn=False, why=why)


E139_LABELS = {'F_inhadd100': "the family's calyx, pn_inhadd 100", 'F_inhadd314': "the family's calyx, pn_inhadd 314",
               'F_inhadd1000': "the family's calyx, pn_inhadd 1000",
               'F_inhadd3141': "the family's calyx, pn_inhadd 3141",
               'B_inhadd314': 'APL→bouton ×0.155, pn_inhadd 314', 'B_inhadd1000': 'APL→bouton ×0.155, pn_inhadd 1000',
               'K_inhadd1000': 'KC–KC ×3.24, pn_inhadd 1000', 'A0_inhadd314': 'A0 refitted at pn_inhadd 314'}


def e139_points() -> List[dict]:
    """E139's settings: F at the four rungs (the family's calyx, a ladder of the normalization fitted to no flies) on
    all eight wirings against all the flies (score_megamat.csv's 'all' half, score_panels.csv, score_apl16.csv),
    the class beside as read on T; B and K (A1's APL→bouton ×0.155 or KC–KC ×3.24, fitted on S against the fit half)
    and step 3's refit A0 at 314 on T against the held-out half (classes.json; the APL effect on T from
    score_apl16.csv); hollow unless classed as holding (rule 3: every one of them), no point outside the realism bands
    or with a hold that fails; A1 beside is E136's arm."""
    d = RESULTS / 'E139_normalization_ladder'
    cl = json.loads((d / 'classes.json').read_text())
    mega, pan, apl = (pd.read_csv(d / f) for f in ('score_megamat.csv', 'score_panels.csv', 'score_apl16.csv'))
    lv, spreads = levels(), cl.get('spreads') or {}
    ref = cl.get('reference', 'family')
    out = []
    for v in [*cl.get('variants', []), *cl.get('fitted', [])]:
        r = cl['classes'].get(v) or {}
        label, cls = E139_LABELS.get(v, v), r.get('cls')
        flags = _flags(apl, v)
        if _hold_fails(r.get('reasons')):
            out.append(_undrawn('E139', label, NOT_HELD + ': ' + '; '.join(_hold_fails(r.get('reasons')))))
            continue
        drawn, why = (False, 'outside the realism bands: ' + flags[0]) if flags else (True, '')
        hw = _class_hollow(cls, r.get('reasons') or [])
        hollow = {'*': hw} if hw else None
        if r.get('setting') == 'F':
            comps = _table_comps(mega, pan, apl, v, ref, WIRINGS, 'all', spreads, ALL8, ALL8_PANEL,
                                 'all eight wirings (fitted to none)')
            out += _emit('E139', label, comps, lv, drawn, why, hollow,
                         f'{cls} (its class, read on T against the held-out half)', fitted=False)
        else:
            comps = {**_class_comps(r), **_apl_comps(apl, v, ref, T_WIRINGS, 'the held-out wirings (T)')}
            out += _emit('E139', label, comps, lv, drawn, why, hollow, cls, r.get('control_fall'))
    if cl.get('beside_arm'):
        out.append(_undrawn('E139', f"{cl['beside_arm']} (beside)", "E136's arm, drawn under E136"))
    return out


E141_LABELS = {'D1.8_inhadd1000': "D1.8: KC–KC ×1.8, on E139's B at 1000",
               'D5.8_inhadd1000': "D5.8: KC–KC ×5.83, on E139's B at 1000",
               'N_inhadd1000': "N: KC–KC ×3.24 with APL on, ×1 with APL off, on E139's B at 1000",
               'G_inhadd1000': "G: APL's inhibition of the claws as a gain, KC–KC ×1, on E139's B at 1000",
               'P_inhadd1000': "P: a uniform neuron in the matrix's place, on E139's B at 1000",
               'T_inhadd1000': "T: one neuron per KC type, on E139's B at 1000",
               'S_inhadd1000': "S: the matrix shuffled within type pairs, on E139's B at 1000"}


def e141_points() -> List[dict]:
    """E141's arms, all at E139's B at 1000 (A1's APL→bouton ×0.155, fitted on S against the fit half), so read as
    fitted: on T against the held-out half (classes.json; the APL effect on T from score_apl16.csv); hollow unless
    classed as holding (rule 3: D1.8 and D5.8 hold); no point outside the realism bands (classes.json's flagged
    wirings or score_apl16.csv's flags: S, P), with a hold that fails (T) or with no comparison on T; B and K beside are
    E139's, A1 E136's."""
    d = RESULTS / 'E141_kckc_rise_forms'
    cl = json.loads((d / 'classes.json').read_text())
    apl = pd.read_csv(d / 'score_apl16.csv')
    lv, ref = levels(), cl.get('reference', 'family')
    out = []
    for v in cl.get('variants', []):
        r = cl['classes'].get(v) or {}
        label, cls = E141_LABELS.get(v, v), r.get('cls')
        flags = [f'{w}: {f}' for w, f in (r.get('flagged_wirings') or {}).items()] or _flags(apl, v)
        comps = {**_class_comps(r), **_apl_comps(apl, v, ref, T_WIRINGS, 'the held-out wirings (T)')}
        if flags:
            drawn, why = False, 'outside the realism bands: ' + flags[0]
        elif _hold_fails(r.get('reasons')):
            drawn, why = False, NOT_HELD + ': ' + '; '.join(_hold_fails(r.get('reasons')))
        elif not (r.get('megamat') or {}).get('complete'):
            drawn, why = False, 'no reading on T: ' + '; '.join(r.get('reasons') or [str(cls)])
        else:
            drawn, why = True, ''
        hw = _class_hollow(cls, r.get('reasons') or [])
        out += _emit('E141', label, comps, lv, drawn, why, {'*': hw} if hw else None, cls, r.get('control_fall'))
    for v in cl.get('beside_e139') or []:
        out.append(_undrawn('E141', f'{v} (beside)', "E139's setting, drawn under E139"))
    if cl.get('beside_arm'):
        out.append(_undrawn('E141', f"{cl['beside_arm']} (beside)", "E136's arm, drawn under E136"))
    return out


E142_LABELS = {'E142_a1stage_imposed': "the dendrites imposed on A1's own PN stage (the amended main arm)",
               'E142_imposed_aplpn': "the dendrites imposed on the family's PN stage, APL→bouton ×0.155 alone",
               'E142_imposed_a1': "the dendrites imposed on the family's PN stage, A1's two calyx factors (the lock's "
                                  "main arm)"}
# E142's dendrites are set from all the PN flies (imposed), so its dendrite values and the stage errors that hold them
# are hollow
E142_HOLLOW = ('pn_dendrite', 'stage_error')
E142_NO_CLASS = "no class of rule 3: the step's own classes, not a hold"


def e142_points() -> List[dict]:
    """E142's step-1 arms with A1's calyx factors (fitted on S against the fit half), read as fitted: on T against the
    held-out half (score_step1_megamat.csv), the APL16 somata on T (score_step1_panels.csv) and the APL effect on T
    (score_step1_apl16.csv), each against the unrefitted family there; every point hollow (scored with the step's own
    classes, not rule 3's; the dendrites and the stage error also set from all the PN flies); no point for an arm
    whose hold fails (the lock's main arm, NOT_HELD); the class as amended (step1_classes.json) in the popup. No point
    for the permuted controls, for the arms read on fewer than four wirings (stopped by claim) and for E135's imposed
    arms re-scored beside (drawn under E135)."""
    d = RESULTS / 'E142_dendrite_bouton_step'
    sc = json.loads((d / 'step1_classes.json').read_text())
    mega, pan, apl = (pd.read_csv(d / f) for f in ('score_step1_megamat.csv', 'score_step1_panels.csv',
                                                   'score_step1_apl16.csv'))
    e136 = _e136_dir() / 'classes.json'
    spreads = {**(json.loads(e136.read_text()).get('spreads', {}) if e136.exists() else {}), **sc.get('spreads', {})}
    lv = levels()
    carry, hides = sc['classes'].get('dendrites_carry') or {}, sc['classes'].get('apl_hides_step') or {}
    n_read: Dict[str, int] = {}
    for row in sc.get('holds_table') or []:
        if str(row.get('panel', 'megamat')) == 'megamat':
            n_read[row['variant']] = int(row['n_read'])
    out = []
    for v in sc.get('variants', []):
        if v in ('family', 'A1'):
            continue
        label = E142_LABELS.get(v, v)
        if v.startswith('E135_'):
            out.append(_undrawn('E142', label, "E135's arm re-scored beside, drawn under E135"))
            continue
        if v.endswith('_perm'):
            out.append(_undrawn('E142', label, 'a permuted control'))
            continue
        if n_read.get(v, 0) < MIN_WIRINGS:
            out.append(_undrawn('E142', label, f'read on {n_read.get(v, 0)} of 8 wirings (the rest stopped by claim '
                                                'after a hold failed)'))
            continue
        bad = sc.get('holds_failed', {}).get(v) or []
        flags = _flags(apl, v, 'variant')
        if bad:
            out.append(_undrawn('E142', label, f'{NOT_HELD}: a hold fails on {len(bad)} of 8 wirings'))
            continue
        if v == carry.get('arm'):
            cls = f"the dendrites carry under the calyx that holds: {carry.get('cls')}"
        elif v in (hides.get('kept') or []):
            cls = f"APL at the boutons hides the step: {hides.get('cls')}"
        else:
            cls = 'no class'
        # scored without rule 3's classes, so hollow at every series (ruled 2026-11-26)
        hollow = {'*': E142_NO_CLASS,
                  **{k: 'the PN dendrites are set from all the PN flies (imposed); ' + E142_NO_CLASS
                     for k in E142_HOLLOW}}
        comps = _table_comps(mega, pan, apl, v, 'family', T_WIRINGS, 'held', spreads, T_HELD, T_PANEL,
                             'the held-out wirings (T)', apl_key='variant')
        drawn, why = (False, 'outside the realism bands: ' + flags[0]) if flags else (True, '')
        out += _emit('E142', label, comps, lv, drawn, why, hollow, cls)
    return out


E144_CALYCES = {'H': "the hemibrain with 'argmax' claws", 'L': "the male CNS's left calyx",
                'R': "the male CNS's right calyx"}
E144_SETTINGS = {'F': 'F (pn_inhadd 1000)', 'B': 'B (pn_inhadd 1000, APL→bouton ×0.155)',
                 'D1.8': "D1.8 (B with the calyx's KC–KC ×1.8)", 'A1': "A1 (B with the calyx's KC–KC ×3.24)",
                 'D5.8': "D5.8 (B with the calyx's KC–KC ×5.83)", 'Dabs': "Dabs (B with L's KC–KC at 85.6)"}


E144_HOLLOW = ("a variant on another brain's wiring: drawn from its calyx's own family against all the flies plus its "
               "change from that family; no staircase steps on it")


def _e144_levels(fam: str, wirings: Sequence[str], mega: pd.DataFrame, panels: pd.DataFrame,
                 apl: pd.DataFrame) -> Dict[str, float]:
    """A male calyx's levels: its own family (E130's fit of the family there) read against all the flies, the mean over
    `wirings` (the megamat stages and their stage error from score_megamat.csv's half 'all', the held-out panels' KC
    somata from score_panels.csv, the calibrated APL effect from score_apl16.csv), as `levels` reads the hemibrain's
    family over its eight wirings; a series missing on any of `wirings` has no level."""
    out = {}

    def put(k, v):
        if all(w in v and np.isfinite(v[w]) for w in wirings):
            out[k] = float(np.mean([v[w] for w in wirings]))
    if len(mega):
        for k, col in (('stage_error', 'stage'), *((s, s) for s in STAGES)):
            put(k, _by_wiring(mega, 'variant', fam, col, half='all'))
    if len(panels):
        p = panels[panels.stage.astype(str) == 'kc_soma'] if 'stage' in panels else panels
        for k, (pan, cond) in PANEL_SERIES.items():
            put(k, _by_wiring(p, 'variant', fam, 'error', panel=pan, condition=cond))
    if len(apl):
        for k, col in APL_SERIES.items():
            put(k, _by_wiring(apl, 'variant', fam, col))
    return out


def e144_points() -> List[dict]:
    """E144's settings on the male CNS's two calyces (and on the hemibrain with 'argmax' claws), each against its
    calyx's own family (E130's fit of the family there) on the same wirings. Ruled 2026-11-26 (a Fable subagent,
    claude-fable-5-2, spawned from viki 19:31): a setting on another brain's wiring is drawn from its calyx's own
    family's value against all the flies (`_e144_levels`) plus its change from that family (the megamat stages against
    the held-out half, as its settings carry A1's factors fitted on S against the fit half; the held-out panels against
    all their flies, classes.json; the APL effect on the same wirings, score_apl16.csv), and every such point is hollow
    (E144_HOLLOW). No point for a setting read on fewer than four wirings (every one on H's real wiring, F and Dabs),
    outside the realism bands, or without a class (L's D5.8, whose real wiring's hold failed); the calyces' families
    are the references."""
    d = RESULTS / 'E144_a1_malecns_calyces'
    cl = json.loads((d / 'classes.json').read_text())
    apl = pd.read_csv(d / 'score_apl16.csv')
    mega = pd.read_csv(d / 'score_megamat.csv') if (d / 'score_megamat.csv').exists() else pd.DataFrame()
    panels = pd.read_csv(d / 'score_panels.csv') if (d / 'score_panels.csv').exists() else pd.DataFrame()
    out = []
    for v, r in cl['classes'].items():
        cx, st = r.get('calyx', v.split('_')[0]), r.get('setting', v.split('_', 1)[-1])
        label = f'{E144_CALYCES.get(cx, cx)}: {E144_SETTINGS.get(st, st)}'
        ws = list(r.get('wirings') or [])
        if len(ws) < MIN_WIRINGS:
            out.append(_undrawn('E144', label, f'read on {len(ws)} wiring{"s" if len(ws) != 1 else ""}'))
            continue
        cls, ref = r.get('cls'), r.get('reference', f'{cx}_family')
        if _hold_fails(r.get('reasons')):
            out.append(_undrawn('E144', label, NOT_HELD + ': ' + '; '.join(_hold_fails(r.get('reasons')))))
            continue
        if cls == 'no class':
            out.append(_undrawn('E144', label, 'no class: ' + '; '.join(r.get('reasons') or [])))
            continue
        flags = [f'{w}: {f}' for w, f in (r.get('flagged_wirings') or {}).items()] or _flags(apl, v, 'variant')
        drawn, why = (False, 'outside the realism bands: ' + flags[0]) if flags else (True, '')
        name = E144_CALYCES.get(cx, cx)
        where = f"{name}'s {len(ws)} wirings against the held-out half, beside the family's own fit there"
        pwhere = f"{name}'s {len(ws)} wirings against all the panel's flies, beside the family's own fit there"
        comps = {**_class_comps(r, where, pwhere),
                 **_apl_comps(apl, v, ref, ws, f"{name}'s {len(ws)} wirings, beside the family's own fit there",
                              'variant')}
        lv = _e144_levels(ref, ws, mega, panels, apl)
        out += _emit('E144', label, comps, lv, drawn, why, {'*': E144_HOLLOW}, cls, r.get('control_fall'),
                     level_note=f"{name}'s own family's ({len(ws)} wirings, all flies)")
    for cx in sorted({str(r.get('calyx')) for r in cl['classes'].values()}):
        out.append(_undrawn('E144', f'{E144_CALYCES.get(cx, cx)}: the family',
                            "the reference of its calyx's settings (E130's fit of the family there)"))
    return out


S_WIRINGS = tuple(w for w in WIRINGS if w not in T_WIRINGS)
E145_CELLS = ('X', 'F', 'W')


def e145_points() -> List[dict]:
    """E145's three new cells, each A1's form fitted on one set of four wirings against one half of the flies and read
    on the other four against the other half (classes.json's `cells`, E136's form): drawn at the level plus its change
    from the unrefitted family where it is read (the megamat stages against its read half, the held-out panels against
    all their flies, the APL effect on its read wirings from score_apl16.csv), with the control fit's fall beside. Every
    cell is hollow by its class (rule 3: X has none, at KC–KC's bound; F one stage at another's cost; W strength at
    another view's cost); no point outside the realism bands or with a hold that fails. Cell A1 is E136's arm, drawn
    under E136."""
    d = RESULTS / 'E145_stage_fit_replication'
    cl = json.loads((d / 'classes.json').read_text())
    apl = pd.read_csv(d / 'score_apl16.csv') if (d / 'score_apl16.csv').exists() else pd.DataFrame()
    lv, cf = levels(), cl.get('control_fall')
    out = []
    for c in E145_CELLS:
        cell = cl.get('cells', {}).get(c)
        if not cell:
            continue
        r, label = cell.get('cls') or {}, f"{c}: {cell.get('label', c)}"
        ws = list(cell.get('read_wirings') or [])
        if len(ws) < MIN_WIRINGS:
            out.append(_undrawn('E145', label, f'read on {len(ws)} wiring{"s" if len(ws) != 1 else ""}'))
            continue
        cls, reasons = r.get('cls'), r.get('reasons') or []
        if _hold_fails(reasons):
            out.append(_undrawn('E145', label, NOT_HELD + ': ' + '; '.join(_hold_fails(reasons))))
            continue
        flags = _flags(apl, cell.get('variant', f'E145_{c}'))
        drawn, why = (False, 'outside the realism bands: ' + flags[0]) if flags else (True, '')
        name = 'T' if set(ws) == set(T_WIRINGS) else 'S' if set(ws) == set(S_WIRINGS) else f'{len(ws)} wirings'
        where = f"the read wirings ({name}) against the held-out half ({cell.get('read_half')})"
        comps = {**_class_comps(r, where, f"the read wirings ({name}) against all the panel's flies"),
                 **_apl_comps(apl, cell.get('variant', f'E145_{c}'), 'family', ws, f'the read wirings ({name})')}
        hw = _class_hollow(cls, reasons)
        out += _emit('E145', label, comps, lv, drawn, why, {'*': hw} if hw else None, cls,
                     r.get('control_fall') if r.get('control_fall') is not None else cf)
    if (cl.get('cells') or {}).get('A1'):
        out.append(_undrawn('E145', 'A1 (beside)', "E136's arm, drawn under E136"))
    return out


ADAPTERS: Dict[str, Callable[[], List[dict]]] = {'E134': e134_points, 'E135': e135_points, 'E136': e136_points,
                                                 'E139': e139_points, 'E141': e141_points, 'E142': e142_points,
                                                 'E144': e144_points, 'E145': e145_points}
# the canonical names (CLAUDE.md), for an experiment not yet in the index as completed
CANONICAL = {'E134': 'the models as populations at every recorded stage',
             'E135': 'the antennal lobe held to the recorded PN dendrites',
             'E136': 'the best family fitted to all megamat stages at once',
             'E139': 'the weaker antennal-lobe normalization without the calyx refit, along its strength',
             'E141': 'the rise of the inhibition among KCs beside the weaker antennal-lobe normalization, along its '
                     'dose and in other forms',
             'E142': 'the step from the PN dendrites to the PN boutons by itself, with the recorded dendrites under '
                     "the stage fit's calyx",
             'E144': "the weaker antennal-lobe normalization and the rise of the inhibition among KCs on the male "
                     "CNS's two calyces",
             'E145': 'the fit of all megamat stages made again on the other wirings and the other halves of the flies'}


def floors() -> Dict[str, float]:
    """The flies' floors the plot draws: the megamat stages' (E134's flies_views.csv, strength by odor, the megamat
    stage error's as the root mean square of the four) and the APL-effect loss's and pattern part's (the registry)."""
    v = pd.read_csv(_e134_dir() / 'flies_views.csv')
    v = v[(v.panel == 'megamat') & (v.level == 'view') & (v.name == 'strength') & v.stage.isin(STAGES)]
    f = {str(s): float(g.floor.iloc[0]) for s, g in v.groupby('stage')}
    reg = json.loads(H.data_targets_path().read_text())
    out = dict(stage_error=float(np.sqrt(np.mean([f[s] ** 2 for s in STAGES]))) if set(STAGES) <= set(f) else np.nan,
               stage_min=min(f.values()) if f else np.nan, stage_max=max(f.values()) if f else np.nan)
    for k, key in (('apl_loss', 'apl_effect_rmse_floor'), ('apl_pattern', 'apl_effect_pattern_floor')):
        x = reg.get(key)
        out[k] = float(x['value']) if isinstance(x, dict) else np.nan
    return out


def collect(include: Sequence[str] = ()) -> pd.DataFrame:
    """Every adapter's points, dated by the experiment's completion in the index (or now, for `include`); per
    experiment and series its point (the lowest filled drawn variant, else the lowest hollow one: 'best_of_exp'), and
    each series' running best over the filled points ('running_best', 'new_best')."""
    from reports.summary_progress_plot import experiment_index
    idx = experiment_index().set_index('exp')
    now = datetime.now().replace(second=0, microsecond=0)
    rows = []
    for exp, fn in ADAPTERS.items():
        if exp in idx.index:
            # the index's name up to its first parenthesis (the rest is the design in brief)
            when, name = idx.loc[exp, 'completed'], str(idx.loc[exp, 'name']).split(' (')[0]
        elif exp in include:
            when, name = now, CANONICAL.get(exp, exp) + ' (running)'
        else:
            continue
        for r in fn():
            rows.append(dict(r, completed=when, name=name))
    pts = pd.DataFrame(rows)
    if not len(pts):
        return pts
    pts['best_of_exp'] = False
    ok = pts.drawn & np.isfinite(pts.value.astype(float))
    for _, g in pts[ok].groupby(['exp', 'series']):
        f = g[g.filled]
        pts.loc[(f if len(f) else g).value.astype(float).idxmin(), 'best_of_exp'] = True
    pts = pts.sort_values(['completed', 'exp', 'series', 'value']).reset_index(drop=True)
    pts['running_best'] = np.nan
    pts['new_best'] = False
    for _, g in pts[pts.best_of_exp].groupby('series'):
        best = np.inf
        for i in g.sort_values('completed').index:
            if pts.loc[i, 'filled'] and float(pts.loc[i, 'value']) < best - 1e-12:
                best = float(pts.loc[i, 'value'])
                pts.loc[i, 'new_best'] = True
            pts.loc[i, 'running_best'] = best if np.isfinite(best) else np.nan
    return pts


def plot_html(pts: pd.DataFrame, fl: Dict[str, float], lv: Optional[Dict[str, float]] = None) -> str:
    """The figure: a panel per kind of series (one-fly units; the APL effect), each series' experiments' points (hollow
    ones open) and its running best as a staircase, the reference lines fixed in the legend; `lv` the series' levels
    (`levels`), named in a fitted point's popup."""
    lv = lv or {}
    best = pts[pts.best_of_exp]
    panels = []
    for pk, title, refs in (('stage', 'one-fly units (lower is better)',
                             [dict(y=1.0, label='one fly: as far from the flies\' mean as one fly is from the others',
                                   cls='ref-one'),
                              dict(y=fl.get('stage_error', np.nan), cls='ref-floor',
                                   label=f"the flies' floor (the megamat stages {fl.get('stage_min', np.nan):.2f} to "
                                         f"{fl.get('stage_max', np.nan):.2f}; the stage error's "
                                         f"{fl.get('stage_error', np.nan):.2f})")]),
                            ('apl', 'the APL effect on the APL16 panel, a description (lower is better)',
                             [dict(y=fl.get('apl_loss', np.nan), cls='ref-floor',
                                   label=f"the flies' floor of the loss ({fl.get('apl_loss', np.nan):.3f}) and of the "
                                         f"pattern part ({fl.get('apl_pattern', np.nan):.3f})"),
                              dict(y=fl.get('apl_pattern', np.nan), cls='ref-floor', label='')])):
        series = []
        for k, (label, panel) in SERIES.items():
            if panel != pk:
                continue
            sub = best[best.series == k].sort_values('completed')
            if not len(sub):
                continue
            series.append(dict(key=k, label=label, light=SERIES_COLORS[k][0], dark=SERIES_COLORS[k][1], points=[
                dict(t=r.completed.isoformat(timespec='minutes'), exp=r.exp, name=r.name, variant=r.variant,
                     reading=r.reading, value=round(float(r.value), 4), filled=bool(r.filled), why=r.why or '',
                     best=None if not np.isfinite(r.running_best) else round(float(r.running_best), 4),
                     new=bool(r.new_best), level=round(float(lv[k]), 4) if k in lv else None,
                     fit=json.loads(r.popup) if r.popup else None) for r in sub.itertuples()]))
        if series:
            panels.append(dict(key=pk, title=title, ymax=Y_MAX[pk],
                               refs=[r for r in refs if np.isfinite(r['y'])], series=series))
    data_json = json.dumps(panels).replace('&', '&amp;').replace('</', '<\\/')
    tbl = ['<table><tr><th>statistic</th><th class="num">best so far</th><th>set by</th></tr>']
    for p in panels:
        for s in p['series']:
            setters = [q for q in s['points'] if q['new']]
            if not setters:
                continue
            setter = setters[-1]
            tbl.append(f"<tr><td>{html.escape(s['label'])}</td><td class=\"num\">{setter['value']:.3f}</td>"
                       f"<td>{html.escape(setter['exp'] + ': ' + setter['name'] + ', ' + setter['variant'])} "
                       f"({html.escape(setter['reading'])})</td></tr>")
    tbl.append('</table>')
    css = """
<style>
.sprog { position: relative; --p-ink: var(--ink, #1f2328); --p-muted: var(--muted, #57606a); --p-line: var(--line, #d8dee4); --p-bg: var(--bg, #fff); margin: 14px 0 6px; }
.sprog svg { width: 100%; height: auto; display: block; font: 12px -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; }
.sprog .grid line { stroke: var(--p-line); stroke-width: 1; }
.sprog .axis text { fill: var(--p-muted); } .sprog .axis line { stroke: var(--p-line); }
.sprog .ptitle { fill: var(--p-ink); font-weight: 600; }
.sprog .step { fill: none; stroke-width: 2; stroke-linejoin: round; }
.sprog .ref-one { stroke: var(--p-ink); stroke-width: 1.5; stroke-dasharray: 2,4; }
.sprog .ref-floor { stroke: var(--p-muted); stroke-width: 1.2; stroke-dasharray: 7,4; }
.sprog .pt { stroke: var(--p-bg); stroke-width: 1.5; cursor: pointer; } .sprog .pt.hit { fill: transparent; stroke: none; }
.sprog .legend { display: flex; flex-wrap: wrap; gap: 6px 18px; margin: 8px 0 0; font-size: 0.86em; color: var(--p-ink); }
.sprog .legend .item { display: flex; align-items: center; gap: 7px; cursor: pointer; user-select: none; }
.sprog .legend .item.fixed { cursor: default; }
.sprog .legend .sw { width: 14px; height: 14px; border-radius: 3px; border: 1px solid rgba(0,0,0,0.15); }
.sprog .legend .sw.ref-one { background: none; border: 0; height: 0; width: 22px; border-top: 2px dotted var(--p-ink); border-radius: 0; }
.sprog .legend .sw.ref-floor { background: none; border: 0; height: 0; width: 22px; border-top: 2px dashed var(--p-muted); border-radius: 0; }
.sprog .legend .item.off { opacity: 0.4; } .sprog .legend .item.off .sw { background: var(--p-line) !important; }
.sprog .tip { position: absolute; pointer-events: none; background: var(--p-bg); color: var(--p-ink); border: 1px solid var(--p-line); border-radius: 6px; padding: 7px 9px; font-size: 0.82em; line-height: 1.35; box-shadow: 0 2px 8px rgba(0,0,0,0.12); max-width: 380px; display: none; z-index: 5; }
.sprog .tip .muted { color: var(--p-muted); }
.sprog .note { color: var(--p-muted); font-size: 0.86em; margin-top: 6px; }
.sprog details { margin-top: 8px; font-size: 0.9em; } .sprog details summary { cursor: pointer; color: var(--p-muted); }
</style>"""
    script = """
<script>
(function () {
  var root = document.getElementById('stage-progress-plot');
  var panels = JSON.parse(root.getAttribute('data-panels'));
  var dark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
  var NS = 'http://www.w3.org/2000/svg';
  function el(tag, attrs, parent) { var e = document.createElementNS(NS, tag); for (var k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
  var times = [];
  panels.forEach(function (P) { P.series.forEach(function (s) { s.points.forEach(function (p) { p.time = new Date(p.t).getTime(); times.push(p.time); }); }); });
  var t0 = Math.min.apply(null, times), t1 = Math.max.apply(null, times);
  var pad = (t1 - t0) * 0.06 || 6 * 3600000; t0 -= pad; t1 += pad;
  var tip = root.querySelector('.tip');
  function f3(v, sign) { return v === null || v === undefined ? '–' : (sign && v > 0 ? '+' : '') + Number(v).toFixed(3); }
  function showTip(evt, s, p, seg) {
    var h = (seg ? '<b>best so far on: ' + s.label + '</b><br>' : '<b>' + p.exp + ': ' + p.name + '</b><br>completed ' + p.t.replace('T', ' ') + '<br>') +
      s.label + ': <b>' + p.value.toFixed(3) + '</b>' + (seg ? ' (' + p.exp + ')' : '') + '<br>variant: ' + p.variant;
    if (p.fit) {
      h += '<br><span class="muted">drawn at ' + (p.fit.level_note ? 'the ' + p.fit.level_note + ' ' + f3(p.fit.level) : 'the unrefitted family’s ' + f3(p.level) + ' (eight wirings, all flies)') + ' plus its change where it was not fitted, on ' + p.fit.where +
        ': the variant ' + f3(p.fit.arm) + ', the family ' + f3(p.fit.family) + ', change ' + f3(p.fit.change, true) +
        (p.fit.se !== undefined && p.fit.se !== null ? ' (SE ' + f3(p.fit.se) + ', bar ' + f3(p.fit.bar, true) + ')' : '') +
        '; class: ' + p.fit.cls + (p.fit.control_fall !== null && p.fit.control_fall !== undefined ? '; the control fit’s fall ' + f3(p.fit.control_fall) : '') + '</span>';
    } else {
      h += '<br><span class="muted">read on ' + p.reading + '</span>';
    }
    if (!p.filled && p.why) h += '<br><span class="muted">hollow: ' + p.why + '</span>';
    if (!seg && p.new) h += '<br>a new best at the time';
    tip.innerHTML = h;
    tip.style.display = 'block';
    var r = root.getBoundingClientRect();
    var left = evt.clientX - r.left + 14, top = evt.clientY - r.top - 10;
    if (left + 400 > r.width) left = Math.max(0, evt.clientX - r.left - 410);
    tip.style.left = left + 'px'; tip.style.top = top + 'px';
  }
  function hideTip() { tip.style.display = 'none'; }
  var groups = {};
  panels.forEach(function (P, pi) {
    var W = 920, Hh = pi === 0 ? 400 : 230, m = {l: 58, r: 18, t: 26, b: 40};
    var iw = W - m.l - m.r, ih = Hh - m.t - m.b;
    function x(t) { return m.l + (t - t0) / (t1 - t0) * iw; }
    function y(v) { v = Math.max(0, Math.min(P.ymax, v)); return m.t + (P.ymax - v) / P.ymax * ih; }
    var div = document.createElement('div'); root.querySelector('.chart').appendChild(div);
    var svg = el('svg', {viewBox: '0 0 ' + W + ' ' + Hh, role: 'img', 'aria-label': P.title + ', the best model so far against completion time'}, div);
    var grid = el('g', {'class': 'grid'}, svg), axis = el('g', {'class': 'axis'}, svg);
    var tt = el('text', {x: m.l, y: 15, 'class': 'ptitle'}, svg); tt.textContent = P.title;
    var nt = 5;
    for (var i = 0; i <= nt; i++) {
      var v = P.ymax * i / nt;
      el('line', {x1: m.l, x2: W - m.r, y1: y(v), y2: y(v), 'stroke-dasharray': '2,4'}, grid);
      var tx = el('text', {x: m.l - 8, y: y(v) + 4, 'text-anchor': 'end'}, axis); tx.textContent = P.ymax < 1 ? v.toFixed(2) : v.toFixed(1);
    }
    P.refs.forEach(function (r) { el('line', {x1: m.l, x2: W - m.r, y1: y(r.y), y2: y(r.y), 'class': r.cls}, svg); });
    var d0 = new Date(t0); d0.setMinutes(0, 0, 0); d0.setHours(Math.ceil(d0.getHours() / 6) * 6);
    var first = true;
    for (var d = new Date(d0); d.getTime() <= t1; d.setHours(d.getHours() + 6)) {
      var t = d.getTime(), hh = d.getHours();
      el('line', {x1: x(t), x2: x(t), y1: m.t + ih, y2: m.t + ih + 5}, axis);
      var tl = el('text', {x: x(t), y: m.t + ih + 18, 'text-anchor': 'middle'}, axis);
      tl.textContent = (hh === 0 || first) ? (d.getMonth() + 1) + '/' + d.getDate() + ' ' + ('0' + hh).slice(-2) + ':00' : ('0' + hh).slice(-2) + ':00';
      first = false;
    }
    el('line', {x1: m.l, x2: W - m.r, y1: m.t + ih, y2: m.t + ih}, axis);
    P.series.forEach(function (s) {
      var color = dark ? s.dark : s.light;
      var g = el('g', {'data-key': s.key}, svg); groups[s.key] = g;
      var setters = s.points.filter(function (p) { return p.new; });
      setters.forEach(function (p, i) {
        var next = setters[i + 1];
        var xs = x(p.time), xe = next ? x(next.time) : W - m.r, yy = y(p.best);
        var dd = 'M' + xs + ' ' + (i === 0 ? yy : y(setters[i - 1].best)) + ' V' + yy + ' H' + xe;
        el('path', {d: dd, 'class': 'step', stroke: color}, g);
        var hit = el('path', {d: dd, 'class': 'step', stroke: 'transparent', 'stroke-width': 12}, g);
        hit.addEventListener('mousemove', function (evt) { showTip(evt, s, p, true); });
        hit.addEventListener('mouseleave', hideTip);
      });
      s.points.forEach(function (p) {
        if (p.filled) el('circle', {cx: x(p.time), cy: y(p.value), r: p.new ? 5 : 4, fill: color, 'fill-opacity': p.new ? 1 : 0.55, 'class': 'pt'}, g);
        else el('circle', {cx: x(p.time), cy: y(p.value), r: 4, fill: 'none', stroke: color, 'stroke-width': 1.6}, g);
        var hit = el('circle', {cx: x(p.time), cy: y(p.value), r: 11, 'class': 'pt hit'}, g);
        hit.addEventListener('mousemove', function (evt) { showTip(evt, s, p, false); });
        hit.addEventListener('mouseleave', hideTip);
      });
    });
    var legend = document.createElement('div'); legend.className = 'legend'; div.appendChild(legend);
    P.refs.forEach(function (r) {
      if (!r.label) return;
      var it = document.createElement('span'); it.className = 'item fixed';
      var sw = document.createElement('span'); sw.className = 'sw ' + r.cls;
      var tx = document.createElement('span'); tx.textContent = r.label;
      it.appendChild(sw); it.appendChild(tx); legend.appendChild(it);
    });
    P.series.forEach(function (s) {
      var it = document.createElement('span'); it.className = 'item';
      var sw = document.createElement('span'); sw.className = 'sw'; sw.style.background = dark ? s.dark : s.light;
      var tx = document.createElement('span'); tx.textContent = s.label;
      it.appendChild(sw); it.appendChild(tx); legend.appendChild(it);
      it.addEventListener('click', function () { var off = it.classList.toggle('off'); groups[s.key].style.display = off ? 'none' : ''; });
    });
  });
})();
</script>"""
    return (css + '\n<figure class="sprog" id="stage-progress-plot" data-panels=\'' + data_json.replace("'", '&#39;')
            + "'>\n<div class=\"chart\" style=\"position:relative\"></div><div class=\"tip\"></div>\n"
            '<div class="note">A fitted model is drawn at the unrefitted family\'s value on its eight wirings against '
            'all flies plus the model\'s change from the family where it was not fitted (the held-out wirings and, at a '
            'fitted stage, the held-out half of the flies). Each point is one experiment\'s best model on that statistic, '
            'among the models inside the realism bands and read on four wirings or more (filled and larger when it was '
            'the best so far; open when the model has no class or the value holds a stage it was fitted to or set from; '
            'points above the top of a panel are drawn at its top); the staircase is the best so far over the filled '
            'points. Hover a point or a step for the experiment, the model and its readings; click a legend entry to '
            'hide or show its line.</div>\n'
            '<details><summary>Table: the best so far on each statistic</summary>' + ''.join(tbl) + '</details>\n'
            '</figure>' + script + '\n')


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--include', nargs='*', default=[], help='running experiments to date now (a preview)')
    ap.add_argument('--out', type=Path, default=OUT, help='where to write (default results/summary_progress/)')
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    pts = collect(a.include)
    if not len(pts):
        print('no experiment with stage-target tables is completed yet; nothing written')
        return
    pts.to_csv(a.out / 'stage_progress_points.csv', index=False)
    pts[pts.new_best][['series', 'exp', 'name', 'completed', 'variant', 'value', 'reading', 'popup']].to_csv(
        a.out / 'stage_progress_best.csv', index=False)
    (a.out / 'stage_progress_plot.html').write_text(plot_html(pts, floors(), levels()))
    pd.set_option('display.width', 200)
    print(pts[pts.best_of_exp].pivot_table(index='series', columns='exp', values='value').round(3).to_string())
    print(f"wrote {a.out / 'stage_progress_plot.html'}, {int(pts.best_of_exp.sum())} points over {pts.exp.nunique()} "
          f"experiments")


if __name__ == '__main__':
    main()
