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
classes). No point for a variant outside the realism bands, one read on fewer than four wirings, a control fit or a fit
to permuted flies. Hollow points (shown; no staircase steps on them): a fitted variant its experiment's score leaves
without a class, and a value that holds a stage the variant was fitted to or set from, on flies not held out (the PN
dendrites of the antennal-lobe fits and of the imposed recordings, and a megamat stage error that includes them). An
experiment's point is its lowest filled variant, else its lowest hollow one; each series' staircase is the running best
of its filled points.

The experiments read (`ADAPTERS`; one function each, from the experiment's own committed tables):
- E134, the models as populations at every recorded stage: the best family as fitted and held at the megamat flies'
  response rate (score_stage_table.csv, each stage's default reading; writeup_panels.csv; writeup_apl16_models.csv);
  its beside models are read on one wiring each (no point);
- E135, the antennal lobe held to the recorded PN dendrites: its arms on eight wirings, but the permuted controls
  (score_scores.csv, score_apl_effect.csv, score_realism.csv);
- E136, the best family fitted to all megamat stages at once: its arms but the permuted control (classes.json's changes
  on T against the held-out half, score_apl16.csv, in its measure's folder).
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


def e135_points() -> List[dict]:
    """E135's arms but the permuted controls: the megamat stages and the APL16 panel (score_scores.csv, the default
    reading), the calibrated APL effect (score_apl_effect.csv), on eight wirings and fitted to no KC or bouton flies, so
    drawn at their own values; the dendrites and the stage errors that hold them hollow (E135_HOLLOW); no point for an
    arm on fewer than four wirings or outside the realism bands (score_realism.csv)."""
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
            hollow = k in E135_HOLLOW
            out.append(_point('E135', label, k, float(p[k].mean()), reading, drawn, not hollow,
                              why or ('the PN dendrites are fitted to or set from all the PN flies' if hollow else '')))
        a = s[(s.group == f'apl16_{arm}') & (s.stage == 'kc_soma') & (s.panel == 'apl16')]
        for c in ('on', 'off'):
            x = a[a.condition == c]
            if len(x):
                out.append(_point('E135', label, f'apl16_{c}', _mean(x.view_strength),
                                  f'the APL16 panel, {x.wiring.nunique()} wirings, against all recordings',
                                  drawn and x.wiring.nunique() >= MIN_WIRINGS, True, why))
        e = ae[ae.arm == arm]
        if len(e):
            for k, col in (('apl_loss', 'cal_apl_effect_rmse'), ('apl_pattern', 'cal_apl_effect_rmse_centered')):
                out.append(_point('E135', label, k, _mean(e[col]), f'{e.wiring.nunique()} wirings',
                                  drawn and e.wiring.nunique() >= MIN_WIRINGS, True, why))
    return out


E136_ARMS = {'A0': 'A0, the family refitted', 'A1': 'A1, with the PN normalization searched',
             'A2': "A2, APL's inhibition of the claws as a gain", 'A3': 'A3, the per-claw threshold free',
             'A4': 'A4, without the inhibition among KCs'}
E136_PANELS = {'validation2': 'validation2/validation2/on', 'cmix': 'mixtures/cmix/on', 'apl16_on': 'apl16/apl16/on',
               'apl16_off': 'apl16/apl16/off'}


def e136_points() -> List[dict]:
    """E136's arms but the permuted control (fitted): each drawn at the level plus its change from the unrefitted family
    on T (classes.json: against the held-out half at the megamat stages, against all the panel's flies on a held-out
    panel; the APL effect from score_apl16.csv, the mean over T of the arm minus the family); hollow without a class;
    no point outside the realism bands (any wiring's flag in score_apl16.csv)."""
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
        filled = cls not in ('no class',)
        if drawn and not filled:
            why = 'no class: ' + '; '.join(r.get('reasons') or [])
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


ADAPTERS: Dict[str, Callable[[], List[dict]]] = {'E134': e134_points, 'E135': e135_points, 'E136': e136_points}
# the canonical names (CLAUDE.md), for an experiment not yet in the index as completed
CANONICAL = {'E134': 'the models as populations at every recorded stage',
             'E135': 'the antennal lobe held to the recorded PN dendrites',
             'E136': 'the best family fitted to all megamat stages at once'}


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
      h += '<br><span class="muted">drawn at the unrefitted family’s ' + f3(p.level) + ' (eight wirings, all flies) plus its change where it was not fitted, on ' + p.fit.where +
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
