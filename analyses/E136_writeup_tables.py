"""The tables of E136's write-up (the best family fitted to all megamat stages at once), from the files its parts wrote
to results/E136_stage_population_fit/measure_<measure>/ (the score's score_*.csv, classes.json and
consistency_<measure>.csv; the fits' fit_<arm>.json and fit_history_<arm>.csv), from E134's flies tables (the floors
and the no-step yardstick) and, for the APL effect's amplitude and shape parts, which the score does not write, from the
APL16 reads in the archive (`--no-apl` skips them). No simulation and no fit.

Writes writeup_tables.md and writeup_values.csv (the same numbers, long: arm, stage, panel, condition, half, statistic,
value, for the summary's placeholders) beside the score's files:
- the stage table: strength by odor at each megamat stage on the held-out wirings (T) against the held-out half of the
  flies, every arm beside the unrefitted family, with the flies' floor and one fly against the others (E134, all flies),
  the flat and half-to-half yardsticks and E134's no-step yardstick; the megamat stage error against the held-out half,
  all flies and the fit half on T, and on the fitted wirings (S) against the fit half as fitted;
- the classes: each arm's change of the megamat stage error on T with its standard error and bar, each stage's change
  against its spread, the guard on the other views, the held-out panels, the control's fall, the reasons and the class;
- the fits: evaluations, restarts, the best and the classed point, the factors, the single-factor returns, bounds;
- the held-out panels (KC somata, strength and all views), the strength split (spread and order of each part), the
  two-view error and the measure's consistency (beside);
- realism: the four bands on the APL16 input per arm (flagged wirings named), the megamat operating point, the largest
  spike counts, the fitted scales beside the family's, and the APL effect (held out) with its amplitude and shape parts.

Usage (from the repository root): ./venv/bin/python -m analyses.E136_writeup_tables [--measure strength]
    [--score-dir <dir>] [--no-apl]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from experiments import E136_stage_population_fit as F

STAGES = F.LOSS_STAGES
STAGE_NAMES = {'pn_dendrite': 'PN dendrites', 'pn_bouton': 'PN boutons', 'kc_claw': 'KC claws',
               'kc_soma': 'KC somata'}
CONTROL = F.CONTROL
ARM_ORDER = (F.FAMILY, 'A0', 'A1', 'A2', 'A3', 'A4', CONTROL)
ARM_NAMES = {
    F.FAMILY: 'the family (unrefitted)',
    'A0': 'A0, the family refitted',
    'A1': "A1, with E135's antennal lobe",
    'A1_search': 'A1, the PN normalization searched (E135 handed on none)',
    'A2': "A2, APL's inhibition of the claws as a gain",
    'A3': 'A3, the per-claw threshold searched',
    'A3_after': 'A3, claws read after the threshold',
    'A4': 'A4, without the inhibition among KCs',
    CONTROL: 'the control (A0 on odor-permuted flies)',
}
PANEL_NAMES = {'validation2/validation2/on': 'validation2', 'mixtures/kmix/on': 'kiwi mixtures',
               'mixtures/cmix/on': 'control mixtures', 'apl16/apl16/on': 'APL16, APL on',
               'apl16/apl16/off': 'APL16, APL off'}
# the panels scored beside the class (the kiwi mixtures: Tom, 2026-11-25 10:42, in hal's terminal, relayed by hal)
BESIDE_KEYS = tuple('/'.join(k) for k in F.BESIDE_PANELS)
BESIDE_NOTE = ("the kiwi mixtures are scored beside the class, not in it (Tom, 2026-11-25 10:42, relayed by hal: \"go "
               "with your read on the mixture panels\"); the control mixtures stay in it, read as uploaded")
READING_NAMES = {'components_first_block': "the components' first block", 'components_second_block':
                 "the components' second block"}
FACTOR_NAMES = {'aplpn': 'APL→bouton', 'kckc': 'KC–KC', 'claw_thr': 'per-claw threshold (Hz)',
                'inhadd': 'PN normalization, additive term', 'inhsc': 'PN normalization, scale'}
SCALE_COLS = (('wAPLPN_scale', 'APL→bouton'), ('wAPLKC_scale', 'APL→KC'), ('wKCKC_scale', 'KC–KC'),
              ('fixed_thr', 'somatic threshold'), ('claw_thr', 'per-claw threshold (Hz)'),
              ('pn_inhadd', 'PN normalization, additive'), ('pn_inhsc', 'PN normalization, scale'))


# ---------------------------------------------------------------------------
# helpers

def _read_json(p: Path) -> dict:
    return json.loads(p.read_text()) if p.exists() else {}


def _mr(x, fmt: str = '.3f') -> str:
    """The mean of `x` with its range, or the one value; '–' if empty."""
    x = pd.Series(x, dtype=float).dropna()
    if not len(x):
        return '–'
    if len(x) == 1 or np.isclose(x.min(), x.max()):
        return f'{x.iloc[0]:{fmt}}'
    return f'{x.mean():{fmt}} ({x.min():{fmt}} to {x.max():{fmt}})'


def _rng(x, fmt: str = '.3f') -> str:
    x = pd.Series(x, dtype=float).dropna()
    if not len(x):
        return '–'
    return f'{x.min():{fmt}}' if np.isclose(x.min(), x.max()) else f'{x.min():{fmt}} to {x.max():{fmt}}'


def _f(v, fmt: str = '.3f') -> str:
    try:
        v = float(v)
    except (TypeError, ValueError):
        return '–'
    return '–' if not np.isfinite(v) else f'{v:{fmt}}'


def _yn(v) -> str:
    return 'yes' if bool(v) else 'no'


def _rms(x) -> float:
    x = np.asarray(list(x), dtype=float)
    return float(np.sqrt(np.mean(x ** 2))) if len(x) else np.nan


def arm_label(arm: str, fits: Mapping[str, dict]) -> str:
    """An arm's name in the tables; A1 by the mode E135's verdict gave it."""
    if arm == 'A1' and ((fits.get('A1') or {}).get('e135') or {}).get('mode') == 'search':
        return ARM_NAMES['A1_search']
    return ARM_NAMES.get(arm, arm)


class Values:
    """The long table of numbers (writeup_values.csv)."""

    def __init__(self):
        self.rows: List[dict] = []

    def add(self, arm, stage, statistic, value, panel='megamat', condition='on', half='held'):
        try:
            v = float(value)
        except (TypeError, ValueError):
            return
        if np.isfinite(v):
            self.rows.append(dict(arm=arm, stage=stage, panel=panel, condition=condition, half=half,
                                  statistic=statistic, value=v))

    def frame(self) -> pd.DataFrame:
        """One row per key (arm, stage, panel, condition, half, statistic): a number two tables both give is written
        once, so that a placeholder selects one row (reports/report_values.py refuses two). Raises ValueError if two
        tables give one key different values."""
        keys = ['arm', 'stage', 'panel', 'condition', 'half', 'statistic']
        df = pd.DataFrame(self.rows, columns=keys + ['value'])
        n = df.groupby(keys, sort=False).value.agg(lambda s: s.round(12).nunique())
        if (n > 1).any():
            raise ValueError(f'values that differ under one key: {n[n > 1].index.tolist()}')
        return df.drop_duplicates(keys).reset_index(drop=True)


# ---------------------------------------------------------------------------
# inputs

def load(score_dir: Path, fit_dir: Path) -> dict:
    """Everything the tables read: the score's files from `score_dir`, the fits' from `fit_dir`, E134's flies tables."""
    d = dict(classes=_read_json(score_dir / 'classes.json'),
             mega=pd.read_csv(score_dir / 'score_megamat.csv'),
             panels=pd.read_csv(score_dir / 'score_panels.csv'),
             apl=pd.read_csv(score_dir / 'score_apl16.csv') if (score_dir / 'score_apl16.csv').exists()
             else pd.DataFrame(),
             prov=_read_json(score_dir / 'provenance_score.json'))
    p = score_dir / f'consistency_{F.MEASURE}.csv'
    d['consistency'] = pd.read_csv(p) if p.exists() else pd.DataFrame()
    fits, hist, dets = {}, {}, {}
    for a in ARM_ORDER[1:]:
        f = _read_json(fit_dir / f'fit_{a}.json')
        if f:
            fits[a] = f
        h = fit_dir / f'fit_history_{a}.csv'
        if h.exists():
            hist[a] = pd.read_csv(h, float_precision='round_trip')
        r = fit_dir / f'determined_{a}.csv'
        if r.exists():
            dets[a] = pd.read_csv(r, float_precision='round_trip')
    d['fits'], d['hist'], d['dets'] = fits, hist, dets
    e = F.e134_dir()
    d['flies_views'] = pd.read_csv(e / 'flies_views.csv') if (e / 'flies_views.csv').exists() else pd.DataFrame()
    d['no_step'] = pd.read_csv(e / 'flies_no_step.csv') if (e / 'flies_no_step.csv').exists() else pd.DataFrame()
    return d


def arms_present(d: dict) -> List[str]:
    have = set(d['mega'].arm)
    return [a for a in ARM_ORDER if a in have]


def flies_row(d: dict) -> Dict[str, dict]:
    """Per stage, the megamat wild-type flies' floor and one fly against the others on strength (E134, all flies), and
    E134's no-step yardstick of the step into the stage."""
    out = {}
    v, n = d['flies_views'], d['no_step']
    for st in STAGES:
        rec = {}
        if len(v):
            x = v[(v.stage == st) & (v.panel == 'megamat') & (v.level == 'view') & (v.name == 'strength')]
            if len(x):
                r = x.iloc[0]
                rec.update(floor=r.floor, loo_mean=r.loo_mean, loo_max=r.loo_max, n_flies=int(r.n_flies))
        if len(n):
            x = n[(n.name == 'view_strength') & n.step.str.endswith(f'-> {st}__megamat__wild_type')]
            if len(x):
                rec['no_step'] = float(x.iloc[0].no_step)
        out[st] = rec
    return out


# ---------------------------------------------------------------------------
# the tables

def stage_table(d: dict, V: Values) -> List[str]:
    """Strength by odor per arm and megamat stage on T against the held-out half (mean over T, range), the megamat
    stage error (mean over T and the error of T's mean views), with the flies' rows and the yardsticks."""
    m, fits = d['mega'], d['fits']
    yard = ((d['classes'].get('beside') or {}).get('yardsticks') or {})
    fl = flies_row(d)
    L = ['Strength by odor in one-fly units (lower is better), megamat panel, against the held-out half of the flies, '
         'on the held-out wirings (T: P3, P4, M2, M3): the mean over the four wirings (their range). "Common" is the '
         'error of the four wirings\' mean views (the error common to the wirings).', '',
         '| model or yardstick | ' + ' | '.join(STAGE_NAMES[s] for s in STAGES) +
         ' | megamat stage error | common |', '|---|' + '---:|' * (len(STAGES) + 2)]
    if any(fl.values()):
        L.append("| the flies' floor (E134, all flies) | " +
                 ' | '.join(_f(fl[s].get('floor'), '.2f') for s in STAGES) + ' | – | – |')
        L.append('| one fly against the others, mean (largest) | ' +
                 ' | '.join(f"{_f(fl[s].get('loo_mean'), '.2f')} ({_f(fl[s].get('loo_max'), '.2f')})" for s in STAGES)
                 + ' | – | – |')
        for s in STAGES:
            V.add('flies', s, 'floor', fl[s].get('floor'), half='all')
            V.add('flies', s, 'loo_mean', fl[s].get('loo_mean'), half='all')
            V.add('flies', s, 'loo_max', fl[s].get('loo_max'), half='all')
    if yard:
        for key, name in (('flat', 'flat profile (every odor at the panel\'s mean)'),
                          ('half_to_half', 'half to half (the fit half\'s mean view)')):
            vals = [yard.get(s, {}).get(key) for s in STAGES]
            L.append(f'| {name} | ' + ' | '.join(_f(v) for v in vals) +
                     f" | {_f(_rms([v for v in vals if v is not None]))} | – |")
            for s, v in zip(STAGES, vals):
                V.add(f'yardstick_{key}', s, 'strength', v)
            V.add(f'yardstick_{key}', 'four_stages', 'stage_error', _rms([v for v in vals if v is not None]))
    if any('no_step' in fl[s] for s in STAGES):
        L.append("| no step (the upstream recordings as the prediction; E134, all flies) | " +
                 ' | '.join(_f(fl[s].get('no_step')) for s in STAGES) + ' | – | – |')
        for s in STAGES:
            V.add('yardstick_no_step', s, 'strength', fl[s].get('no_step'), half='all')
    for a in arms_present(d):
        if a == CONTROL:
            continue
        x = m[(m.arm == a) & (m.half == 'held') & m.wiring.isin(F.T_WIRINGS)]
        c = m[(m.arm == a) & (m.half == 'held') & (m.wiring == 'common')]
        if not len(x):
            continue
        cc = c.iloc[0] if len(c) and a != F.FAMILY else None
        L.append(f'| {arm_label(a, fits)} | ' + ' | '.join(_mr(x[s]) for s in STAGES) +
                 f" | {_mr(x.stage)} | {_f(cc.stage) if cc is not None else '–'} |")
        for s in STAGES:
            V.add(a, s, 'strength_mean', x[s].mean())
            V.add(a, s, 'strength_min', x[s].min())
            V.add(a, s, 'strength_max', x[s].max())
        V.add(a, 'four_stages', 'stage_error_mean', x.stage.mean())
        V.add(a, 'four_stages', 'stage_error_min', x.stage.min())
        V.add(a, 'four_stages', 'stage_error_max', x.stage.max())
        if cc is not None:
            V.add(a, 'four_stages', 'stage_error_common', cc.stage)
            for s in STAGES:
                V.add(a, s, 'strength_common', cc[s])
    fam8 = m[(m.arm == F.FAMILY) & (m.half == 'held') & m.wiring.isin(F.WIRINGS)]
    if len(fam8):
        L.append('| the family on all eight wirings (the spreads\' basis) | ' +
                 ' | '.join(_mr(fam8[s]) for s in STAGES) + f' | {_mr(fam8.stage)} | – |')
        V.add(F.FAMILY, 'four_stages', 'stage_error_mean_8', fam8.stage.mean())
    fc = m[(m.arm == F.FAMILY) & (m.half == 'held') & (m.wiring == 'common')]
    if len(fc):
        V.add(F.FAMILY, 'four_stages', 'stage_error_common_8', fc.iloc[0].stage)
    return L + ['']


def halves_table(d: dict, V: Values) -> List[str]:
    """The megamat stage error on strength per arm against the held-out half, all flies and the fit half on T (means
    over T), and on S against the fit half as fitted (the fit's classed point; the family: A0's first evaluation)."""
    m, fits, hist = d['mega'], d['fits'], d['hist']
    L = ['The megamat stage error on strength (one-fly units, lower is better): on T against each set of flies, mean '
         'over the four wirings, and on the fitted wirings (S) against the fit half, as fitted.', '',
         '| model | T, held-out half | T, all flies | T, fit half | S, fit half (as fitted) |',
         '|---|---:|---:|---:|---:|']
    fam_s = None
    if 'A0' in hist and len(hist['A0']):
        h0 = hist['A0']
        e0 = h0[h0['eval'] == 0]
        fam_s = float(e0.iloc[0].loss_stages) if len(e0) else None
    for a in arms_present(d):
        if a == CONTROL:
            continue
        cells = []
        for half in ('held', 'all', 'fit'):
            x = m[(m.arm == a) & (m.half == half) & m.wiring.isin(F.T_WIRINGS)]
            cells.append(_f(x.stage.mean()) if len(x) else '–')
            if len(x):
                V.add(a, 'four_stages', 'stage_error_mean', x.stage.mean(), half=half)
        if a == F.FAMILY:
            s = fam_s
        else:
            f = fits.get(a) or {}
            s = f.get('classed_loss') if f.get('classed_loss') is not None else (f.get('best') or {}).get('loss_sim')
        cells.append(_f(s))
        V.add(a, 'four_stages', 'stage_error_S_fit', s, half='fit_S')
        L.append(f'| {arm_label(a, fits)} | ' + ' | '.join(cells) + ' |')
    return L + ['']


def classes_table(d: dict, V: Values) -> List[str]:
    """Per arm, the class and what it rests on (classes.json)."""
    c = d['classes']
    cls, fits = c.get('classes') or {}, d['fits']
    cf = c.get('control_fall')
    L = [f"The classes, on T against the held-out half, against the unrefitted family there (strength by odor; the "
         f"bar is {F.IMPROVE_SE:g} standard errors over wirings; the control fit's fall against the permuted held-out "
         f"half is {_f(cf)}); {BESIDE_NOTE}.", '',
         '| arm | change of the megamat stage error (SE) | bar | improves | stages worse than their spread | '
         'guard: stages worse on all views | held-out panels worse | fall vs the control\'s | adds on A0 | reasons '
         '| class |', '|---|---:|---:|---|---|---|---|---|---|---|---|']
    V.add(CONTROL, 'four_stages', 'control_fall', cf)
    for a in ARM_ORDER:
        r = cls.get(a)
        if not r:
            continue
        mg = r.get('megamat') or {}
        if not mg.get('complete'):
            L.append(f"| {arm_label(a, fits)} | – | – | – | – | – | – | – | – | {', '.join(r.get('reasons') or [])} | "
                     f"{r.get('cls')} |")
            continue
        ad = r.get('adds_vs_A0')
        adds = '–' if ad is None else f"{_f(ad.get('change'), '+.3f')} ({_yn(r.get('adds'))})"
        L.append(f"| {arm_label(a, fits)} | {_f(mg['change'], '+.3f')} ({_f(mg['se'])}) | {_f(-F.IMPROVE_SE * mg['se'], '+.3f')}"
                 f" | {_yn(mg.get('improves'))} | {', '.join(STAGE_NAMES[s] for s in r.get('stage_worse') or []) or 'none'}"
                 f" | {', '.join(STAGE_NAMES.get(s, s) for s in r.get('guard_worse') or []) or 'none'}"
                 f" | {', '.join(PANEL_NAMES.get(k, k) for k in r.get('panel_worse') or []) or 'none'}"
                 f" | {_f(r.get('fall'))} vs {_f(cf)} | {adds} | {'; '.join(r.get('reasons') or []) or '–'}"
                 f" | {r.get('cls')} |")
        V.add(a, 'four_stages', 'change', mg['change'])
        V.add(a, 'four_stages', 'se', mg['se'])
        V.add(a, 'four_stages', 'bar', F.IMPROVE_SE * mg['se'])
        V.add(a, 'four_stages', 'fall', r.get('fall'))
        if ad is not None:
            V.add(a, 'four_stages', 'change_vs_A0', ad.get('change'))
        alt = r.get('claws_after_threshold')
        if alt:
            am = alt.get('megamat') or {}
            L.append(f"| {ARM_NAMES['A3_after']} | {_f(am.get('change'), '+.3f')} ({_f(am.get('se'))}) | – | "
                     f"{_yn(am.get('improves'))} | – | – | – | – | – | – | {alt.get('cls')} |")
            V.add('A3_after', 'four_stages', 'change', am.get('change'))
    return L + ['']


def stage_changes_table(d: dict, V: Values) -> List[str]:
    """Per arm and stage, the change of the strength error on T (SE) against the stage's spread, and the guard's change
    of the error over all views against its spread."""
    c = d['classes']
    cls, fits = c.get('classes') or {}, d['fits']
    sp = c.get('spreads') or {}
    spa = ((c.get('beside') or {}).get('spreads_all_views') or {})
    L = ['Each stage on T against the held-out half: the change of the strength error from the unrefitted family '
         '(standard error), then the change of the population error over all views (the guard of rule 3); the spreads '
         'are the unrefitted family\'s SD over its eight wirings.', '',
         '| arm | ' + ' | '.join(f'{STAGE_NAMES[s]}, strength | all views' for s in STAGES) + ' |',
         '|---|' + '---:|' * (2 * len(STAGES)),
         '| spread | ' + ' | '.join(f"{_f(sp.get(s))} | {_f(spa.get(s))}" for s in STAGES) + ' |']
    for s in STAGES:
        V.add(F.FAMILY, s, 'spread', sp.get(s))
        V.add(F.FAMILY, s, 'spread_all_views', spa.get(s))
    V.add(F.FAMILY, 'four_stages', 'spread', sp.get('stage'))
    V.add(F.FAMILY, 'four_stages', 'spread_all_views', spa.get('stage'))
    for a in ARM_ORDER:
        r = cls.get(a)
        if not r or not (r.get('megamat') or {}).get('complete'):
            continue
        cells = []
        for s in STAGES:
            x = (r.get('stages') or {}).get(s) or {}
            g = (r.get('guard') or {}).get(s) or {}
            mark = ' ↑' if x.get('worse_than_spread') else (' ↓' if x.get('improves') else '')
            gmark = ' ↑' if g.get('worse_than_spread') else ''
            cells.append(f"{_f(x.get('change'), '+.3f')} ({_f(x.get('se'))}){mark} | {_f(g.get('change'), '+.3f')}{gmark}")
            V.add(a, s, 'change', x.get('change'))
            V.add(a, s, 'change_all_views', g.get('change'))
        L.append(f'| {arm_label(a, fits)} | ' + ' | '.join(cells) + ' |')
    return L + ['↑ higher by more than the spread (worse); ↓ lower by more than the bar.', '']


def _panel(r: Mapping, k: str) -> dict:
    """An arm's comparison on one held-out panel: in the class ('panels') or beside it ('panels_beside'; a score from
    before 2026-11-25 10:42 has the kiwi mixtures in 'panels')."""
    return (r.get('panels') or {}).get(k) or (r.get('panels_beside') or {}).get(k) or {}


def panels_table(d: dict, V: Values) -> List[str]:
    """The KC somata's strength error on each held-out panel on T (mean over T), the unrefitted family's and each arm's
    with its change, the spread over the family's eight wirings, and the guard's change on all views; the panels scored
    beside the class (F.BESIDE_PANELS) marked."""
    c = d['classes']
    cls, fits, sp = c.get('classes') or {}, d['fits'], c.get('spreads') or {}
    spa = ((c.get('beside') or {}).get('spreads_all_views') or {})
    keys = ['/'.join(k) for k in F.SCORED_PANELS]
    L = ['The KC somata on the held-out panels, on T, strength by odor (one-fly units, the mean over T); each arm\'s '
         'cell is its value and its change from the unrefitted family, then the change on all views; the column marked '
         '"(beside)" is scored beside the class and gates nothing (Tom, 2026-11-25 10:42).', '',
         '| model | ' + ' | '.join(PANEL_NAMES[k] + (' (beside)' if k in BESIDE_KEYS else '') for k in keys) + ' |',
         '|---|' + '---:|' * len(keys)]
    ref = None
    for a in ARM_ORDER:
        r = cls.get(a)
        if r and (r.get('panels') or {}):
            ref = r
            break
    if ref:
        L.append(f'| {ARM_NAMES[F.FAMILY]} | ' + ' | '.join(_f(_panel(ref, k).get('ref')) for k in keys) + ' |')
        for k in keys:
            V.add(F.FAMILY, 'kc_soma', 'strength_mean', _panel(ref, k).get('ref'), panel=k)
    L.append('| spread (strength / all views) | ' + ' | '.join(f'{_f(sp.get(k))} / {_f(spa.get(k))}' for k in keys)
             + ' |')
    for k in keys:
        V.add(F.FAMILY, 'kc_soma', 'spread', sp.get(k), panel=k)
        V.add(F.FAMILY, 'kc_soma', 'spread_all_views', spa.get(k), panel=k)
    for a in ARM_ORDER:
        r = cls.get(a)
        if not r or not (r.get('panels') or {}):
            continue
        cells = []
        for k in keys:
            x = _panel(r, k)
            g = (r.get('guard') or {}).get(k) or {}
            if not x.get('complete'):
                cells.append('–')
                continue
            mark = ' ↑' if x.get('worse_than_spread') else ''
            gmark = ' ↑' if g.get('worse_than_spread') else ''
            cells.append(f"{_f(x.get('arm'))} ({_f(x.get('change'), '+.3f')}{mark}; {_f(g.get('change'), '+.3f')}{gmark})")
            V.add(a, 'kc_soma', 'strength_mean', x.get('arm'), panel=k)
            V.add(a, 'kc_soma', 'change', x.get('change'), panel=k)
            V.add(a, 'kc_soma', 'change_all_views', g.get('change'), panel=k)
        L.append(f'| {arm_label(a, fits)} | ' + ' | '.join(cells) + ' |')
    return L + ['↑ higher than the unrefitted family by more than the spread (worse).', '']


def e137_table(d: dict, V: Values) -> List[str]:
    """The mixture panels' KC somata under E137's readings of the flies, beside the class (classes.json's
    beside.e137_readings: `compare` on T against the unrefitted family): per model the strength error (mean over T) and
    its change, then the change on all views; the kiwi mixtures' two readings are the bounds Tom asked for."""
    c = d['classes']
    cls, fits = c.get('classes') or {}, d['fits']
    e = (c.get('beside') or {}).get('e137_readings') or {}
    if not e:
        return []
    keys = [k for k in (f'{p}/{r}' for p, rs in F.E137_BESIDE.items() for r in rs)]
    meas = '+'.join(F.MEASURES[F.MEASURE] or ()) or 'all'
    first = next(iter(e.values()))
    L = ["The mixture panels' KC somata under E137's readings of the flies (its lock: baseline-noise units, the flies' "
         "trial noise added to the model after the curve), on T, beside the class: the strength error (one-fly units "
         "of the reading's own flies, the mean over T), its change from the unrefitted family, then the change on all "
         "views. The kiwi mixtures' two readings bound the components' fall across the session.", '',
         '| model | ' + ' | '.join(f"{PANEL_NAMES['mixtures/' + k.split('/')[0] + '/on']}, "
                                   f"{READING_NAMES.get(k.split('/')[1], k.split('/')[1])}" for k in keys) + ' |',
         '|---|' + '---:|' * len(keys),
         f'| {ARM_NAMES[F.FAMILY]} | ' + ' | '.join(_f((first.get(f'{k}/{meas}') or {}).get('ref')) for k in keys) + ' |',
         '| spread (strength / all views) | ' + ' | '.join(
             f"{_f((first.get(f'{k}/{meas}') or {}).get('spread'))} / {_f((first.get(f'{k}/all') or {}).get('spread'))}"
             for k in keys) + ' |']
    for k in keys:
        pan, rd = k.split('/')
        V.add(F.FAMILY, 'kc_soma', 'strength_mean', (first.get(f'{k}/{meas}') or {}).get('ref'),
              panel=f'mixtures/{pan}/on', half=rd)
        V.add(F.FAMILY, 'kc_soma', 'spread', (first.get(f'{k}/{meas}') or {}).get('spread'), panel=f'mixtures/{pan}/on',
              half=rd)
    for a in ARM_ORDER:
        if a not in e or a not in cls:
            continue
        cells = []
        for k in keys:
            pan, rd = k.split('/')
            x, g = e[a].get(f'{k}/{meas}') or {}, e[a].get(f'{k}/all') or {}
            if not x.get('complete'):
                cells.append('–')
                continue
            mark = ' ↑' if x.get('worse_than_spread') else ''
            gmark = ' ↑' if g.get('worse_than_spread') else ''
            cells.append(f"{_f(x.get('arm'))} ({_f(x.get('change'), '+.3f')}{mark}; {_f(g.get('change'), '+.3f')}{gmark})")
            V.add(a, 'kc_soma', 'strength_mean', x.get('arm'), panel=f'mixtures/{pan}/on', half=rd)
            V.add(a, 'kc_soma', 'change', x.get('change'), panel=f'mixtures/{pan}/on', half=rd)
            V.add(a, 'kc_soma', 'change_all_views', g.get('change'), panel=f'mixtures/{pan}/on', half=rd)
        L.append(f'| {arm_label(a, fits)} | ' + ' | '.join(cells) + ' |')
    return L + ['↑ higher than the unrefitted family by more than the spread (worse).', '']


def fits_table(d: dict, V: Values) -> List[str]:
    """Per fit: searched coordinates, evaluations against the budget, restarts and how they ended, the best evaluation,
    the single-factor returns (the loss's rise when a factor is returned to the family's value; where the return leaves
    a response rate that can no longer be held, the rise of the stage error alone beside it, the rest being the hold's
    penalty, from determined_<arm>.csv), the classed point and its factors, and the bounds."""
    fits = d['fits']
    L = ['The fits, on S against the fit half (megamat stage error on strength; the family starts at its first '
         'evaluation).', '',
         '| arm | searched | evaluations / budget | ended | best (eval) | returns: rise of the loss | classed point, '
         'factors | classed loss | at a bound |', '|---|---|---:|---|---:|---|---|---:|---|']
    for a in ARM_ORDER[1:]:
        f = fits.get(a)
        if not f:
            continue
        best = (f.get('best') or {}).get('loss_sim')
        det = f.get('determined') or {}
        split = return_split(f, d.get('dets', {}).get(a))
        rets = '; '.join(f"{FACTOR_NAMES.get(k, k)} {_f(v.get('rise'), '+.4f')}"
                         + (f" (stages {_f(split[k][0], '+.4f')}, the hold's penalty {_f(split[k][1], '.4f')})"
                            if k in split and split[k][1] > 0 else '')
                         + f"{'' if v.get('determined') else ' (not determined)'}" for k, v in det.items()) or '–'
        cz = f.get('classed_z') or f.get('best_z') or {}
        fac = ', '.join(f'{FACTOR_NAMES.get(k, k)} ×{F.coord_value(k, v):.3g}' if k != 'claw_thr'
                        else f'{FACTOR_NAMES[k]} {F.coord_value(k, v):.1f}' for k, v in cz.items())
        cl = f.get('classed_loss') if f.get('classed_loss') is not None else best
        bounds = f.get('at_bound_classed') if f.get('at_bound_classed') is not None else f.get('at_bound')
        ended = '; '.join(str(m).replace('Optimization terminated successfully.', 'converged')
                          .replace('Maximum number of function evaluations has been exceeded.', 'budget')
                          for m in f.get('messages') or [])
        L.append(f"| {arm_label(a, fits)} | {', '.join(FACTOR_NAMES.get(k, k) for k in f.get('searched') or [])} | "
                 f"{f.get('n_evals')} / {f.get('budget')} | {ended} | {_f(best, '.4f')} ({f.get('best_eval')}) | {rets} | "
                 f"{fac} | {_f(cl, '.4f')} | {', '.join(bounds or []) or 'none'} |")
        V.add(a, 'four_stages', 'best_loss', best, half='fit_S')
        V.add(a, 'four_stages', 'classed_loss', cl, half='fit_S')
        V.add(a, 'four_stages', 'n_evals', f.get('n_evals'), half='fit_S')
        for k, v in cz.items():
            V.add(a, 'factor', k, F.coord_value(k, v), half='fit_S')
        for k, v in det.items():
            V.add(a, 'factor', f'{k}_return_rise', v.get('rise'), half='fit_S')
            if k in split:
                V.add(a, 'factor', f'{k}_return_rise_stages', split[k][0], half='fit_S')
                V.add(a, 'factor', f'{k}_return_penalty', split[k][1], half='fit_S')
    if any(s[1] > 0 for a in fits for s in return_split(fits[a], d.get('dets', {}).get(a)).values()):
        L += ['', 'Where a return leaves a response rate that can no longer be held, the rise of the loss is the rise of '
              'the stage error plus the hold\'s penalty, given beside it.']
    return L + ['']


def return_split(fit: dict, det: Optional[pd.DataFrame]) -> Dict[str, Tuple[float, float]]:
    """Per searched factor returned alone (determined_<arm>.csv's '<factor>0' rows): the rise of the stage error from
    the fit's best and the hold's penalty at the returned point, whose sum is the rise of the loss."""
    if det is None or det.empty or 'loss_stages' not in det.columns:
        return {}
    best = float((fit.get('best') or {}).get('loss_stages', np.nan))
    out = {}
    for k in fit.get('determined') or {}:
        r = det[det['label'].astype(str) == f'{k}0']
        if len(r):
            out[k] = (float(r['loss_stages'].iloc[0]) - best, float(r['hold_penalty'].iloc[0]))
    return out


def split_table(d: dict, V: Values) -> List[str]:
    """The strength split of the wirings' common error on T against the held-out half: per stage, the spread and order
    parts of the response-rate part and of the mean-response part, per arm."""
    m, fits = d['mega'], d['fits']
    L = ['The strength split on T, the common error against the held-out half: spread / order of the response-rate part '
         'and of the mean-response part, per stage (level is 0 by construction: each odor is a multiple of the panel\'s '
         'mean).', '',
         '| model | ' + ' | '.join(f'{STAGE_NAMES[s]}, rates | means' for s in STAGES) + ' |',
         '|---|' + '---:|' * (2 * len(STAGES))]
    for a in arms_present(d):
        if a == CONTROL:
            continue
        c = m[(m.arm == a) & (m.half == 'held') & (m.wiring == 'common')]
        if not len(c):
            continue
        r = c.iloc[0]
        cells = []
        for s in STAGES:
            for part in ('rate', 'mean'):
                sp, od = r.get(f'{s}_split_strength_{part}_spread'), r.get(f'{s}_split_strength_{part}_order')
                cells.append(f'{_f(sp, ".2f")} / {_f(od, ".2f")}')
                V.add(a, s, f'split_{part}_spread', sp)
                V.add(a, s, f'split_{part}_order', od)
                V.add(a, s, f'split_{part}_sd_ratio', r.get(f'{s}_split_strength_{part}_sd_ratio'))
        L.append(f'| {arm_label(a, fits)} | ' + ' | '.join(cells) + ' |')
    return L + ['']


def beside_table(d: dict, V: Values) -> List[str]:
    """The two-view error's change on T (beside) and the measure's consistency over the saved evaluations."""
    c = d['classes']
    bs, fits = c.get('beside') or {}, d['fits']
    L = []
    tv = bs.get('two_view') or {}
    if tv:
        L += ['The two-view error (strength and sharing; beside, no class), on T against the held-out half, mean over '
              'T.', '', '| arm | the arm | the unrefitted family | change (SE) |', '|---|---:|---:|---:|']
        for a in ARM_ORDER:
            x = tv.get(a)
            if not x or not x.get('complete'):
                continue
            L.append(f"| {arm_label(a, fits)} | {_f(x.get('arm'))} | {_f(x.get('ref'))} | "
                     f"{_f(x.get('change'), '+.3f')} ({_f(x.get('se'))}) |")
            V.add(a, 'four_stages', 'two_view_mean', x.get('arm'))
            V.add(a, 'four_stages', 'two_view_change', x.get('change'))
            V.add(F.FAMILY, 'four_stages', 'two_view_mean', x.get('ref'))
        L.append('')
    k = d['consistency']
    if len(k):
        L += ["The measure's consistency over the saved evaluations of every arm (design item 7, a description): ρ "
              'over the evaluations between the error against the fit half and against the held-out half.', '',
              '| statistic | evaluations | ρ | yardstick |', '|---|---:|---:|---:|']
        for r in k.itertuples(index=False):
            L.append(f'| {r.statistic} | {int(r.n_evals)} | {_f(r.rho)} | {_f(r.yardstick, ".2f")} |')
            V.add('all_arms', str(r.statistic), 'consistency_rho', r.rho, half='fit_vs_held')
        L.append('')
    return L


def data_targets() -> Dict[str, float]:
    """The data-target registry's values (results/data_targets/data_targets.json), by key; empty if not built."""
    p = F.out_dir().parent.parent / 'data_targets' / 'data_targets.json'
    if not p.exists():
        p = F.out_dir().parent / 'data_targets' / 'data_targets.json'
    if not p.exists():
        return {}
    d = json.loads(p.read_text())
    d = d.get('values', d)
    return {k: v['value'] for k, v in d.items() if isinstance(v, dict) and isinstance(v.get('value'), (int, float))}


def realism_table(d: dict, V: Values) -> List[str]:
    """The four bands on the APL16 input per arm (range over the eight wirings; flagged wirings named), and the largest
    spike counts."""
    a, fits = d['apl'], d['fits']
    if not len(a):
        return []
    L = ['The APL16 input (held out), range over the eight wirings. Bands: response rate with APL on 0.045 to 0.085; '
         'APL-off / APL-on ratio 1.6 to 2.9; response rate ρ with APL on at least 0.37 (the reference model\'s 0.47 '
         'less 0.10); odors per KC within 0.155 (APL on) and 0.223 (APL off) of the flies.', '',
         '| model | response rate, APL on | APL off | ratio | response rate ρ, APL on | APL off | odors per KC, APL on '
         '| APL off | largest spike count, APL on / off | flagged |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    t = data_targets()
    if t:
        g = lambda k, f='.3f': _f(t.get(k), f)
        L.append(f"| the KC data (control flies / APL output blocked (TNT), noise criterion) | {g('noise_ctrl_frac')} | "
                 f"{g('noise_tnt_frac')} | {g('noise_rate_ratio_fraction', '.2f')} | – | – | band "
                 f"{g('noise_ctrl_odors_per_kc_tv_max')} | band {g('noise_tnt_odors_per_kc_tv_max')} | – | – |")
    for arm in ARM_ORDER:
        x = a[a.arm == arm]
        if not len(x):
            continue
        flag = x.realism_flag.fillna('').astype(str)
        named = '; '.join(f'{w}: {fl}' for w, fl in zip(x.wiring, flag) if fl) or 'none'
        L.append(f'| {arm_label(arm, fits)} | {_rng(x.rate_on)} | {_rng(x.rate_off)} | {_rng(x.rate_ratio, ".2f")} | '
                 f'{_rng(x.rate_rho_on, ".2f")} | {_rng(x.rate_rho_off, ".2f")} | {_rng(x.opk_on_tv)} | '
                 f'{_rng(x.opk_off_tv)} | {_rng(x.max_count_on, ".0f")} / {_rng(x.max_count_off, ".0f")} | {named} |')
        for col in ('rate_on', 'rate_off', 'rate_ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv',
                    'max_count_on', 'max_count_off'):
            if col in x:
                V.add(arm, 'kc_soma', f'{col}_min', x[col].min(), panel='apl16', condition='realism', half='all')
                V.add(arm, 'kc_soma', f'{col}_max', x[col].max(), panel='apl16', condition='realism', half='all')
        V.add(arm, 'kc_soma', 'n_flagged', int((flag != '').sum()), panel='apl16', condition='realism', half='all')
    return L + ['']


def operating_table(d: dict, V: Values) -> List[str]:
    """The megamat operating point and the fitted scales on S at each fit's best evaluation (absolute values, range over
    the four wirings), the family's (A0's first evaluation) first."""
    from experiments import E134_population_baseline as B
    from kcq import stage_readings as R
    hist, fits = d['hist'], d['fits']
    try:
        target = B.held_rate(F.PANEL, R.require_lock())
    except Exception:   # a description only
        target = np.nan
    rows = []
    if 'A0' in hist:
        h = hist['A0']
        e0 = h[h['eval'] == 0]
        if len(e0):
            rows.append((F.FAMILY, e0.iloc[0]))
    for a in ARM_ORDER[1:]:
        f, h = fits.get(a), hist.get(a)
        if f and h is not None and len(h):
            b = h[h['eval'] == f.get('best_eval')]
            if len(b):
                rows.append((a, b.iloc[0]))
    if not rows:
        return []
    cols = [(c, n) for c, n in SCALE_COLS if any(f'{c}_{w}' in r.index for _, r in rows for w in F.S_WIRINGS)]
    L = [f"The megamat operating point and the fitted scales on S at each fit's best evaluation (range over the four "
         f"wirings; the flies' response rate with APL on {_f(target, '.4f')}, held to within 0.5 %).", '',
         '| model | response rate, APL on | APL off | ' + ' | '.join(n for _, n in cols) + ' |',
         '|---|---:|---:|' + '---:|' * len(cols)]
    for a, r in rows:
        on = [r.get(f'rate_on_{w}') for w in F.S_WIRINGS]
        off = [r.get(f'rate_off_{w}') for w in F.S_WIRINGS]
        cells = [_rng(on, '.4f'), _rng(off, '.4f')]
        for c, _ in cols:
            vals = [r.get(f'{c}_{w}') for w in F.S_WIRINGS]
            cells.append(_rng(vals, '.4g'))
            V.add(a, 'S', f'{c}_min', pd.Series(vals, dtype=float).min(), half='fit_S')
            V.add(a, 'S', f'{c}_max', pd.Series(vals, dtype=float).max(), half='fit_S')
        L.append(f'| {arm_label(a, fits)} | ' + ' | '.join(cells) + ' |')
    return L + ['']


def apl_effect_parts(arms: Sequence[str]) -> pd.DataFrame:
    """Per arm and wiring, the APL effect's amplitude and shape parts, the spread ratio and Pearson r on the per-pair
    APL effect through the calibrated readout (E67's `full_row` calibration), from the APL16 reads in the archive."""
    from experiments.E67_apl_effect_fits import refs
    from kcq import harness as H
    from kcq import stage_readings as R
    ref = refs()
    rows = []
    for arm in arms:
        for w in F.WIRINGS:
            r = F.load_read(F.ReadPoint(arm, 'apl16', w, None))
            if r is None:
                continue
            s_on, s_off = (r.conditions[c]['spike_counts'] for c in R.CONDITIONS)
            K = H.calibrate_readout(s_on, 'saturating', ref['cal_target'], stat='participation_ratio',
                                    responders_min=2.0)
            on = H.readout_pair_distances(s_on, 'saturating', K)
            off = H.readout_pair_distances(s_off, 'saturating', K)
            sh = H.apl_effect_spread_shape(on, off, ref['graded'])
            st = H.apl_effect_stats(on, off, ref['graded'])
            rows.append(dict(arm=arm, wiring=w, **{k: float(v) for k, v in sh.items()},
                             apl_effect_pearson=float(st.get('apl_effect_pearson', np.nan))))
    return pd.DataFrame(rows)


DIVISOR_HIGH = 10.0   # E92's: a shunt above about tenfold at the responders' claws


def wiring_W(wiring: str) -> pd.DataFrame:
    """A wiring's claw x bouton matrix: E94's local cache (outputs/E94_bouton_share/W/), else copied there from E94's
    model-output archive."""
    import shutil
    from experiments import E94_bouton_share as E94
    from kcq import harness as H
    got = E94.load_wiring(wiring)
    if got is None:
        src = H.shared_dir() / 'data' / 'model_outputs' / 'E94_bouton_share' / 'W'
        for k in ('W_long', 'W_index', 'W_columns', 'wAPLPN'):
            p = src / f'{wiring}_{k}.parquet'
            if p.exists():
                shutil.copyfile(p, E94._w_dir() / p.name)
        got = E94.load_wiring(wiring)
    if got is None:
        raise FileNotFoundError(f"wiring {wiring}'s claw x bouton matrix is neither in E94's cache nor in its archive")
    return got[0]


def sites_row(W: pd.DataFrame, on: Mapping, off: Mapping, divisive: bool) -> dict:
    """Where APL's inhibition is made in one run, from its window means (conditions' 'bouton_win', 'claw_win' and, under
    a divisive claw form, 'responses'): E94's `claw_loss_shares` ('drive_removed', the share of the claws' summed drive
    APL removes; 'bouton_share', the share of that loss made at the boutons); with `divisive`, the divisor read on
    window means, x / state at the claws of each odor's responders (x = W @ the boutons' rates with APL on; E92's divisor
    1 + i is per time step, so this is its window-mean reading): 'divisor_p50', 'divisor_p90' and 'divisor_above10',
    the share above DIVISOR_HIGH."""
    from experiments import E94_bouton_share as E94
    from kcq import harness as H
    from kcq import stage_readings as R
    odors = [str(c) for c in on['bouton_win'].columns]
    E94.check_stage_shapes(W, on)
    b_on = on['bouton_win'][odors].to_numpy(dtype=float)
    c_on = on['claw_win'][odors].to_numpy(dtype=float)
    sh = E94.claw_loss_shares(W, off['bouton_win'][odors].to_numpy(dtype=float), b_on, c_on)
    row = dict(bouton_share=sh['bouton_share'], drive_removed=sh['drive_removed'])
    if not divisive:
        return row
    x = np.asarray(W, dtype=float) @ b_on
    kc = E94.kc_of_claw(W)
    resp = R._abbrev(on['responses'])[odors]
    ids = np.asarray(resp.index.get_level_values(H.KC_ID_LEVEL) if isinstance(resp.index, pd.MultiIndex)
                     else resp.index)
    rr = resp.to_numpy() > 0
    pool = []
    for j in range(len(odors)):
        m = np.isin(kc, ids[rr[:, j]]) & (c_on[:, j] > 0) & (x[:, j] > 0)
        pool.append(x[m, j] / c_on[m, j])
    v = np.concatenate(pool) if pool else np.zeros(0)
    if v.size:
        row.update(divisor_p50=float(np.quantile(v, 0.5)), divisor_p90=float(np.quantile(v, 0.9)),
                   divisor_above10=float((v > DIVISOR_HIGH).mean()))
    return row


def apl_sites(arms: Sequence[str]) -> pd.DataFrame:
    """Per arm (the unrefitted family first), panel and wiring, `sites_row` of its reads: megamat on T against the held-out
    half, APL16 on the eight wirings."""
    from kcq import stage_readings as R
    rows = []
    for arm in [F.FAMILY] + [a for a in arms if a != F.FAMILY]:
        div = bool(F.arm_kws(arm).get('apl_inh_divisive')) if arm != F.FAMILY else False
        for panel, ws, half in ((F.PANEL, F.T_WIRINGS, 'held'), ('apl16', F.WIRINGS, None)):
            for w in ws:
                r = F.load_read(F.ReadPoint(arm, panel, w, half))
                if r is None:
                    continue
                on, off = (r.conditions[c] for c in R.CONDITIONS)
                rows.append(dict(arm=arm, panel=panel, wiring=w, **sites_row(wiring_W(w), on, off, div)))
    return pd.DataFrame(rows)


def apl_sites_table(d: dict, sites: Optional[pd.DataFrame], V: Values) -> List[str]:
    """Where APL's inhibition is made (named in advance in the entry's Realism): the bouton share against the
    connectome's, and the divisor at the responders' claws of the divisive arm, mean (range) over the wirings."""
    from experiments import E94_bouton_share as E94
    if sites is None or not len(sites):
        return []
    fits = d['fits']
    div = 'divisor_above10' in sites and sites['divisor_above10'].notna().any()
    L = [f"Where APL's inhibition is made, window means, mean over the wirings (range): megamat on T, APL16 on the eight. "
         f"\"Removed\" is the share of the claws' summed drive APL removes; \"at the boutons\" the share of that loss "
         f"made at the PN boutons (E94's), against the connectome's {E94.ANATOMICAL_SHARE:g} (APL's calyx synapses onto "
         f"boutons over those onto boutons and KCs)."
         + (f" Under the divisive claw form, the divisor at the claws of each odor's responders read on window means "
            f"(E92's is per time step): its median and the share above {DIVISOR_HIGH:g}." if div else ''), '',
         '| model | megamat: removed | at the boutons | APL16: removed | at the boutons |'
         + (' megamat: divisor, median | above tenfold | APL16: divisor, median | above tenfold |' if div else ''),
         '|---|---:|---:|---:|---:|' + ('---:|---:|---:|---:|' if div else '')]
    for arm in ARM_ORDER:
        x = sites[sites.arm == arm]
        if not len(x):
            continue
        cells = []
        for panel in (F.PANEL, 'apl16'):
            p = x[x.panel == panel]
            for col in ('drive_removed', 'bouton_share'):
                cells.append(_mr(p[col]) if len(p) else '–')
                if len(p):
                    V.add(arm, 'calyx', f'{col}_mean', p[col].mean(), panel=panel)
                    V.add(arm, 'calyx', f'{col}_min', p[col].min(), panel=panel)
                    V.add(arm, 'calyx', f'{col}_max', p[col].max(), panel=panel)
        if div:
            for panel in (F.PANEL, 'apl16'):
                p = x[x.panel == panel]
                has = len(p) and 'divisor_p50' in p and p['divisor_p50'].notna().any()
                cells += [_mr(p['divisor_p50'], '.2f'), _mr(p['divisor_above10'])] if has else ['–', '–']
                if has:
                    for col in ('divisor_p50', 'divisor_p90', 'divisor_above10'):
                        V.add(arm, 'calyx', f'{col}_mean', p[col].mean(), panel=panel)
                        V.add(arm, 'calyx', f'{col}_max', p[col].max(), panel=panel)
        L.append(f'| {arm_label(arm, fits)} | ' + ' | '.join(cells) + ' |')
    V.add('connectome', 'calyx', 'bouton_share', E94.ANATOMICAL_SHARE)
    return L + ['']


def apl_effect_table(d: dict, parts: Optional[pd.DataFrame], V: Values) -> List[str]:
    """The APL effect on the APL16 input (a held-out test, reported, gating nothing), calibrated readout, mean over the
    eight wirings: the loss, the mean APL effect, the pattern part with its amplitude and shape parts and the spread
    ratio, ρ and Pearson r on the per-pair APL effect."""
    a, fits = d['apl'], d['fits']
    if not len(a):
        return []
    L = ['The APL effect on odor-odor distances, APL16 input (held out), calibrated readout, mean over the eight '
         'wirings (the flies\' mean APL effect +0.18; the pattern part\'s floor about 0.029).', '',
         '| model | the APL-effect loss | the mean APL effect | the pattern part | amplitude part | shape part | '
         'spread ratio | ρ on the per-pair APL effect | Pearson r |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for arm in ARM_ORDER:
        x = a[a.arm == arm]
        if not len(x):
            continue
        p = parts[parts.arm == arm] if parts is not None and len(parts) else pd.DataFrame()
        g = lambda col: p[col].mean() if len(p) and col in p else np.nan
        L.append(f"| {arm_label(arm, fits)} | {_f(x.cal_loss.mean())} | {_f(x.cal_mean.mean(), '+.3f')} | "
                 f"{_f(x.cal_pattern.mean())} | {_f(g('apl_effect_amplitude'))} | {_f(g('apl_effect_shape'))} | "
                 f"{_f(g('apl_effect_sd_ratio'), '.2f')} | {_f(x.cal_rho.mean(), '.2f')} | "
                 f"{_f(g('apl_effect_pearson'), '.2f')} |")
        for col, name in (('cal_loss', 'apl_effect_loss'), ('cal_mean', 'apl_effect_mean'),
                          ('cal_pattern', 'apl_effect_pattern'), ('cal_rho', 'apl_effect_rho')):
            V.add(arm, 'kc_soma', f'{name}_mean', x[col].mean(), panel='apl16', condition='apl_effect', half='all')
        for col in ('apl_effect_amplitude', 'apl_effect_shape', 'apl_effect_sd_ratio', 'apl_effect_pearson'):
            V.add(arm, 'kc_soma', f'{col}_mean', g(col), panel='apl16', condition='apl_effect', half='all')
    return L + ['']


def main(score_dir: Optional[Path] = None, fit_dir: Optional[Path] = None, with_apl: bool = True) -> Path:
    """Writes writeup_tables.md and writeup_values.csv into `score_dir` (default the measure's out_dir); `with_apl`
    False skips what reads the saved runs (the APL effect's parts and where APL's inhibition is made)."""
    fit_dir = F.out_dir() if fit_dir is None else Path(fit_dir)
    score_dir = fit_dir if score_dir is None else Path(score_dir)
    d = load(score_dir, fit_dir)
    V = Values()
    c, pv = d['classes'], d['prov']
    arms = [a for a in arms_present(d)]
    parts = apl_effect_parts(arms) if with_apl and len(d['apl']) else None
    sites = apl_sites(arms) if with_apl else None
    L = ['# E136 write-up tables', '',
         f"From the score of {pv.get('time')} ({pv.get('host')}, queue {str(pv.get('queue_commit'))[:7]}): measure "
         f"{c.get('measure')!r}, stages {', '.join(c.get('stages') or [])}; arms scored: "
         f"{', '.join(pv.get('arms') or [])}.", '',
         '## The stage table', ''] + stage_table(d, V) + halves_table(d, V)
    L += ['## The classes', ''] + classes_table(d, V) + stage_changes_table(d, V)
    L += ['## The held-out panels', ''] + panels_table(d, V) + e137_table(d, V)
    L += ['## The fits', ''] + fits_table(d, V)
    L += ['## The strength split', ''] + split_table(d, V)
    L += ['## Beside', ''] + beside_table(d, V)
    L += ['## Realism', ''] + realism_table(d, V) + operating_table(d, V) + apl_sites_table(d, sites, V)
    L += ['## The APL effect (held out)', ''] + apl_effect_table(d, parts, V)
    V.frame().to_csv(score_dir / 'writeup_values.csv', index=False)
    out = score_dir / 'writeup_tables.md'
    out.write_text('\n'.join(L) + '\n')
    print(out)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--measure', choices=tuple(F.MEASURES), default='strength')
    ap.add_argument('--score-dir', type=Path, help="where the score's files are (default the measure's out_dir)")
    ap.add_argument('--no-apl', action='store_true', help="skip what reads the saved runs: the APL effect's amplitude "
                    "and shape parts, and where APL's inhibition is made")
    a = ap.parse_args()
    F.set_measure(a.measure)
    main(a.score_dir, with_apl=not a.no_apl)
