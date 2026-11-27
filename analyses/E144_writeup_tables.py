"""The tables of E144's write-up (the weaker antennal-lobe normalization and the rise of the inhibition among KCs on the
male CNS's two calyces), from the files its parts wrote to results/E144_a1_malecns_calyces/ (the score's classes.json
and score_*.csv; E143's resampled reading, e143_*.csv and e143_intervals.json; step 0's step0_*.csv and
prediction.json; the reads' reads_<host>_<lane>.csv), from E134's flies tables (the floors), E141's score (the
hemibrain family's eight wirings), E130's APL synapse table and the data-target registry (the KC data's values beside
the realism bands). No simulation and no fit.

Writes writeup_tables.md and writeup_values.csv (the same numbers, long: variant, calyx, setting, stage, panel,
condition, half, statistic, value) beside the score's files.

Usage (from the repository root): ./venv/bin/python -m analyses.E144_writeup_tables [--score-dir <dir>]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from analyses import E136_writeup_tables as W
from experiments import E136_stage_population_fit as F
from experiments import E144_a1_malecns_calyces as X
from experiments import E144_score as S

STAGES = W.STAGES
STAGE_NAMES = W.STAGE_NAMES
CALYX_NAMES = {'H': "H: the hemibrain, 'argmax' claws", 'L': 'L: the male CNS, left calyx',
               'R': 'R: the male CNS, right calyx'}
SETTING_NAMES = {'family': "the calyx's family (E130's fit, pn_inhadd 31.4)",
                 'F': 'F: pn_inhadd 1000',
                 'B': 'B: F with APL→bouton ×0.155',
                 'D1.8': 'D1.8: B with KC–KC ×1.8',
                 'A1': 'A1: B with KC–KC ×3.24',
                 'D5.8': 'D5.8: B with KC–KC ×5.83',
                 'Dabs': "Dabs: B with KC–KC at 85.6 (A1's on the hemibrain)"}
PANEL_KEYS = tuple('/'.join(k) for k in F.HOLDS_PANELS) + tuple('/'.join(k) for k in F.BESIDE_PANELS)
PANEL_NAMES = {'validation2/validation2/on': 'validation2', 'mixtures/cmix/on': 'control mixtures',
               'apl16/apl16/on': 'APL16, APL on', 'apl16/apl16/off': 'APL16, APL off',
               'mixtures/kmix/on': 'kiwi mixtures (beside)'}
_f, _mr = W._f, W._mr


def _rng(x, fmt: str = '.3f') -> str:
    """A range over wirings as 'lo to hi', or one value where both print alike."""
    if x is None or not len(x):
        return '–'
    lo, hi = _f(float(np.nanmin(x)), fmt), _f(float(np.nanmax(x)), fmt)
    return lo if lo == hi else f'{lo} to {hi}'


def _names(xs: Sequence[str]) -> str:
    """Stages and panels by name, '; '-separated (a panel's name may hold a comma); 'none' for none."""
    return '; '.join(STAGE_NAMES.get(x, PANEL_NAMES.get(x, x)) for x in xs) or 'none'


class Values:
    """The numbers of the tables, long (writeup_values.csv)."""

    def __init__(self):
        self.rows: List[dict] = []

    def add(self, variant, stage, statistic, value, panel='megamat', condition='on', half='held'):
        if value is None or (isinstance(value, float) and not np.isfinite(value)):
            return
        c, s = (str(variant).split('_', 1) + [''])[:2] if '_' in str(variant) else ('', '')
        self.rows.append(dict(variant=variant, calyx=c, setting=s, stage=stage, panel=panel, condition=condition,
                              half=half, statistic=statistic, value=float(value)))

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)


def label(v: str) -> str:
    c, s = v.split('_', 1)
    return f'{c}, {SETTING_NAMES[s.split("@")[0]]}'


def variants() -> List[str]:
    return [X.variant(c, s) for c in X.CALYCES for s in X.SETTINGS if c in X.SETTINGS[s].calyces]


def load(score_dir: Path) -> dict:
    """Everything the tables read."""
    from experiments import E141_kckc_rise_forms as E141
    rd = lambda n: pd.read_csv(score_dir / n) if (score_dir / n).exists() else pd.DataFrame()
    reads = pd.concat([pd.read_csv(p) for p in sorted(score_dir.glob('reads_*.csv'))], ignore_index=True)
    e = F.e134_dir()
    return dict(classes=W._read_json(score_dir / 'classes.json'), mega=rd('score_megamat.csv'),
                panels=rd('score_panels.csv'), apl=rd('score_apl16.csv'), sites=rd('score_sites.csv'),
                without=rd(f'score_without_{S.WITHOUT_ODOR}.csv'), e137=rd('score_panels_e137.csv'),
                shares=rd('e143_class_shares.csv'), changes=rd('e143_intervals_changes.csv'),
                shortcut=rd('e143_shortcut_check.csv'), e143=W._read_json(score_dir / 'e143_intervals.json'),
                prediction=W._read_json(score_dir / 'prediction.json'), reads=reads,
                step0_m=rd('step0_matrices.csv'), step0_s=rd('step0_scales.csv'),
                prov=W._read_json(score_dir / 'provenance_score.json'),
                prov_r=W._read_json(score_dir / 'provenance_resampled.json'),
                e141_mega=pd.read_csv(E141.out_dir() / 'score_megamat.csv'),
                flies_views=pd.read_csv(e / 'flies_views.csv') if (e / 'flies_views.csv').exists() else pd.DataFrame(),
                no_step=pd.read_csv(e / 'flies_no_step.csv') if (e / 'flies_no_step.csv').exists() else pd.DataFrame())


# ---------------------------------------------------------------------------
# the stage table and the megamat reads

def stage_table(d: dict, V: Values) -> List[str]:
    """Strength by odor per megamat stage against the held-out half: the flies' rows, the hemibrain family on its eight
    wirings (E141's score), and every calyx's family and A1 (mean over the calyx's wirings, range; the real wiring)."""
    m = d['mega']
    fl = W.flies_row(d)
    L = ['Strength by odor in one-fly units (lower is better), megamat panel, against the held-out half of the flies. '
         'A male calyx is read on its four wirings (its real wiring and three within-type claw permutations): the '
         'mean over them (their range), and the real wiring alone. H is its real wiring only. "Common" is the error '
         'of the wirings\' mean views.', '',
         '| model or yardstick | wirings | ' + ' | '.join(STAGE_NAMES[s] for s in STAGES) +
         ' | megamat stage error | common |', '|---|---|' + '---:|' * (len(STAGES) + 2)]
    L.append("| the flies' floor (E134, all flies) | – | " + ' | '.join(_f(fl[s].get('floor'), '.2f') for s in STAGES)
             + ' | – | – |')
    L.append('| one fly against the others, mean (largest) | – | ' +
             ' | '.join(f"{_f(fl[s].get('loo_mean'), '.2f')} ({_f(fl[s].get('loo_max'), '.2f')})" for s in STAGES)
             + ' | – | – |')
    for s in STAGES:
        V.add('flies', s, 'floor', fl[s].get('floor'), half='all')
        V.add('flies', s, 'loo_mean', fl[s].get('loo_mean'), half='all')
    e = d['e141_mega']
    fam8 = e[(e.variant == 'family') & (e.half == 'held') & e.wiring.isin(F.WIRINGS)]
    L.append("| the hemibrain family on its eight wirings ('first' claws; E141's score) | 8 | " +
             ' | '.join(_mr(fam8[s]) for s in STAGES) + f' | {_mr(fam8.stage)} | – |')
    for s in ('stage',) + tuple(STAGES):
        V.add('hemibrain_family8', s, 'min', fam8[s].min())
        V.add('hemibrain_family8', s, 'max', fam8[s].max())
        V.add('hemibrain_family8', s, 'mean', fam8[s].mean())
    for c in X.CALYCES:
        for st in ('family', 'A1'):
            v = X.variant(c, st)
            x = m[(m.variant == v) & (m.half == 'held') & m.wiring.isin(X.WIRINGS_OF[c])]
            r = m[(m.variant == v) & (m.half == 'held') & (m.wiring == 'real')]
            cc = m[(m.variant == v) & (m.half == 'held') & (m.wiring == 'common')]
            if not len(x):
                continue
            if c != 'H':
                L.append(f'| {label(v)} | {len(x)} | ' + ' | '.join(_mr(x[s]) for s in STAGES) +
                         f" | {_mr(x.stage)} | {_f(cc.iloc[0].stage) if len(cc) else '–'} |")
            L.append(f'| {label(v)}, real wiring | real | ' + ' | '.join(_f(r.iloc[0][s]) for s in STAGES) +
                     f' | {_f(r.iloc[0].stage)} | – |')
            for s in ('stage',) + tuple(STAGES):
                V.add(v, s, 'mean', x[s].mean())
                V.add(v, s, 'min', x[s].min())
                V.add(v, s, 'max', x[s].max())
                V.add(v, s, 'real', r.iloc[0][s])
    return L + ['']


def inside_table(d: dict, V: Values) -> List[str]:
    """"The calyx's family is inside the hemibrain's wirings": its strength error on its real wiring against the held-
    out half beside the hemibrain family's range over its eight wirings, per stage."""
    cal = d['classes'].get('calyces') or {}
    L = ["The calyx's family on its real wiring against the hemibrain family's range over its eight wirings (strength "
         'by odor, held-out half). "Outside" is the entry\'s class, read at the boutons, the claws and the somata.', '',
         '| stage | the hemibrain family, eight wirings (range) | ' + ' | '.join(f'{c} family, real' for c in X.CALYCES)
         + ' |', '|---|---:|' + '---:|' * len(X.CALYCES)]
    for st in ('stage',) + tuple(STAGES):
        lo, hi = (cal['H']['family_inside'][st]['hemibrain_min'], cal['H']['family_inside'][st]['hemibrain_max'])
        cells = []
        for c in X.CALYCES:
            r = cal[c]['family_inside'][st]
            cells.append(f"{_f(r['value'])}{'' if r['inside'] else ' (outside)'}")
            V.add(X.variant(c, 'family'), st, 'inside', float(bool(r['inside'])))
        name = 'megamat stage error' if st == 'stage' else STAGE_NAMES[st]
        L.append(f'| {name} | {_f(lo)} to {_f(hi)} | ' + ' | '.join(cells) + ' |')
    return L + ['']


def all_views_table(d: dict, V: Values) -> List[str]:
    """Each calyx's family over all three views (strength, sharing, breadth; rule 3's guard), mean over the calyx's
    wirings, from the guard's reference of the calyx's A1 class (the family on A1's wirings), with the hemibrain
    family's spreads on all views beside."""
    cl = d['classes'].get('classes') or {}
    gsp = d['classes']['spreads_all_views']
    keys = tuple(STAGES) + PANEL_KEYS
    L = ["Each calyx's family over all three views (strength by odor, sharing between odors, breadth among units; "
         "one-fly units, held-out half for megamat, all flies for the held-out panels), mean over the calyx's wirings. "
         "Beside, as yardstick, the hemibrain family's SD over its eight wirings on the same error.", '',
         '| model | wirings | ' + ' | '.join(STAGE_NAMES.get(k, PANEL_NAMES.get(k, k)) for k in keys) + ' |',
         '|---|---|' + '---:|' * len(keys),
         "| the hemibrain family's spread (SD over eight wirings) | – | " + ' | '.join(_f(gsp.get(k)) for k in keys)
         + ' |']
    for c in X.CALYCES:
        r = cl.get(X.variant(c, 'A1')) or {}
        g = r.get('guard') or {}
        cells = [_f((g.get(k) or {}).get('ref')) if (g.get(k) or {}).get('complete') else '–' for k in keys]
        L.append(f"| {label(X.variant(c, 'family'))} | {', '.join(r.get('wirings') or [])} | " + ' | '.join(cells) + ' |')
        for k in keys:
            if (g.get(k) or {}).get('complete'):
                st, pan, cond = (k, 'megamat', 'on') if '/' not in k else ('kc_soma', *k.split('/')[1:])
                V.add(X.variant(c, 'family'), st, 'all_views_mean', g[k]['ref'], panel=pan, condition=cond)
    return L + ['']


def halves_table(d: dict, V: Values) -> List[str]:
    """The megamat stage error per variant against the held-out half, all flies and the fit half (mean over the
    variant's read wirings)."""
    m = d['mega']
    L = ['The megamat stage error on strength (one-fly units, lower is better) per setting, mean over the wirings it '
         'was read on, against each set of flies. The class reads the held-out half.', '',
         '| calyx and setting | wirings | held-out half | all flies | fit half |', '|---|---|---:|---:|---:|']
    for v in variants():
        c = v.split('_')[0]
        cells, n = [], 0
        for h in ('held', 'all', 'fit'):
            x = m[(m.variant == v) & (m.half == h) & m.wiring.isin(X.WIRINGS_OF[c])]
            n = max(n, len(x))
            cells.append(_f(x.stage.mean()) if len(x) else '–')
            if len(x):
                V.add(v, 'stage', 'mean', x.stage.mean(), half=h)
        L.append(f'| {label(v)} | {n} | ' + ' | '.join(cells) + ' |')
    return L + ['']


# ---------------------------------------------------------------------------
# the classes

def _ch(x: Optional[dict], fmt='+.3f') -> str:
    return _f(x.get('change'), fmt) if x and x.get('complete') else '–'


def classes_table(d: dict, V: Values) -> List[str]:
    """Each setting's class against its calyx's family on the wirings it was read on, with the megamat change and its
    bar, what failed it, E139's price and E143's share of resamples giving the scored class; the descriptions on the
    wirings whose holds held beside ('@held')."""
    cl = d['classes'].get('classes') or {}
    sh = d['shares'].set_index('arm') if len(d['shares']) else pd.DataFrame()
    L = ["Each setting against its calyx's family on the same wirings, megamat against the held-out half, strength by "
         'odor in one-fly units (lower is better). The bar is twice the standard error over the wirings (the spread '
         "where one wiring is read); the spreads are the hemibrain family's over its eight wirings. "
         '"Other views" is rule 3\'s guard: the error over all views. The price is E139\'s: the calibrated pattern part '
         'of the APL effect above the family\'s by more than 0.003 against all recordings or a half. The last column is '
         "E143's: the share of 1,000 resamples of the flies that give the scored class.", '',
         '| calyx and setting | wirings | class | megamat change (bar) | stages worse than their spread | other views '
         'worse | held-out panels worse | price on the APL effect | resamples giving the class |',
         '|---|---|---|---:|---|---|---|---:|---:|']
    for v in variants():
        if v.endswith('_family'):
            continue
        r = cl.get(v)
        if not r:
            continue
        for key, rec in ((v, r), (f'{v}@held', r.get('on_wirings_held'))):
            if not rec:
                continue
            m = rec.get('megamat') or {}
            pr = r.get('price_on_apl_effect') or {} if key == v else {}
            share = sh.loc[key, 'share_scored'] if key in getattr(sh, 'index', []) else np.nan
            names = _names
            cls = rec.get('cls')
            if cls == 'no class' and rec.get('reasons'):
                cls = f"no class ({'; '.join(rec['reasons'])})"
            lab = label(v) + (', on the wirings whose holds held (a description)' if key != v else '')
            L.append(f"| {lab} | {', '.join(rec.get('wirings') or [])} | {cls} | "
                     f"{_ch(m)} ({_f(rec.get('bar'), '.3f')}) | {names(rec.get('stage_worse') or [])} | "
                     f"{names(rec.get('guard_worse') or [])} | {names(rec.get('panel_worse') or [])} | "
                     f"{('yes' if pr.get('price') else 'no') + ' (' + _f(pr.get('price_value'), '+.4f') + ')' if pr.get('wirings') else '–'} | "
                     f'{_f(share, ".3f")} |')
            V.add(key, 'stage', 'change', m.get('change'))
            V.add(key, 'stage', 'bar', rec.get('bar'))
            V.add(key, 'stage', 'resampled_share_scored', share)
            if pr.get('price_value') is not None:
                V.add(key, 'apl_effect', 'price_value', pr.get('price_value'), panel='apl16', condition='on_off')
    return L + ['']


def stage_changes_table(d: dict, V: Values) -> List[str]:
    """Each setting's change of each megamat stage (mean over its wirings), and on the guard (all views) at the somata,
    against the spreads."""
    cl = d['classes'].get('classes') or {}
    sp, gsp = d['classes']['spreads'], d['classes']['spreads_all_views']
    L = ["Each setting's change of each megamat stage's strength error from its calyx's family on the same wirings "
         '(mean over the wirings; held-out half; lower is better), and the somata\'s error over all views (rule 3\'s '
         'guard). The spreads are the hemibrain family\'s SD over its eight wirings.', '',
         '| calyx and setting | wirings | ' + ' | '.join(STAGE_NAMES[s] for s in STAGES) +
         ' | KC somata, all views |', '|---|---|' + '---:|' * (len(STAGES) + 1),
         '| the spread | – | ' + ' | '.join(_f(sp[s]) for s in STAGES) + f" | {_f(gsp['kc_soma'])} |"]
    for v in variants():
        r = cl.get(v)
        if not r or 'stages' not in r:
            continue
        L.append(f"| {label(v)} | {', '.join(r['wirings'])} | " +
                 ' | '.join(_ch(r['stages'].get(s)) for s in STAGES) +
                 f" | {_ch((r.get('guard') or {}).get('kc_soma'))} |")
        for s in STAGES:
            V.add(v, s, 'change', (r['stages'].get(s) or {}).get('change'))
        V.add(v, 'kc_soma', 'change_all_views', ((r.get('guard') or {}).get('kc_soma') or {}).get('change'))
    return L + ['']


def panels_table(d: dict, V: Values) -> List[str]:
    """Each setting's change of each held-out panel's KC population error from its calyx's family (mean over the
    wirings), the spreads above."""
    cl = d['classes'].get('classes') or {}
    sp = d['classes']['spreads']
    L = ["The held-out panels' KC somata: each setting's change of the strength error from its calyx's family on the "
         "same wirings (mean over them; lower is better). A change above the panel's spread fails the class; the kiwi "
         'mixtures are beside and gate nothing.', '',
         '| calyx and setting | wirings | ' + ' | '.join(PANEL_NAMES[k] for k in PANEL_KEYS) + ' |',
         '|---|---|' + '---:|' * len(PANEL_KEYS),
         '| the spread | – | ' + ' | '.join(_f(sp[k]) for k in PANEL_KEYS) + ' |']
    for v in variants():
        r = cl.get(v)
        if not r or 'panels' not in r:
            continue
        allp = dict(r['panels'], **(r.get('panels_beside') or {}))
        L.append(f"| {label(v)} | {', '.join(r['wirings'])} | " + ' | '.join(_ch(allp.get(k)) for k in PANEL_KEYS) +
                 ' |')
        for k in PANEL_KEYS:
            inp, pan, cond = k.split('/')
            V.add(v, 'kc_soma', 'change', (allp.get(k) or {}).get('change'), panel=pan, condition=cond)
    return L + ['']


def ladder_table(d: dict, V: Values) -> List[str]:
    """The dose ladder (a description): validation2's change from the calyx's family along B, D1.8, A1 and D5.8 on the
    real wiring, beside the hemibrain's (E141's, on its 'first' real wiring), and on the male calyces' four wirings."""
    cal = d['classes'].get('calyces') or {}
    cl = d['classes'].get('classes') or {}
    L = ["Validation2's KC somata along the dose: the change of the strength error from the calyx's family (lower is "
         'better; the spread 0.038). The hemibrain\'s row is the entry\'s, from E141: the mean over the '
         'hemibrain\'s four held-out wirings (P3, P4, M2, M3).', '',
         '| calyx | wirings | ' + ' | '.join(S.DOSES) + ' |', '|---|---|' + '---:|' * len(S.DOSES),
         "| the hemibrain ('first' claws; E141) | four held out | " + ' | '.join(_f(S.HEMIBRAIN_LADDER[s], '+.3f')
                                                                     for s in S.DOSES) + ' |']
    for c in X.CALYCES:
        lad = cal[c]['dose_ladder']
        L.append(f'| {CALYX_NAMES[c]} | real | ' + ' | '.join(_f(lad[s]['change'], '+.3f') for s in S.DOSES) + ' |')
        for s in S.DOSES:
            V.add(X.variant(c, s), 'kc_soma', 'ladder_real', lad[s]['change'], panel='validation2')
        if c in X.MALE:
            cells = []
            for s in S.DOSES:
                r = cl.get(X.variant(c, s)) or {}
                p = (r.get('panels') or {}).get('validation2/validation2/on')
                ok = r.get('wirings') and len(r['wirings']) == 4
                cells.append(_ch(p) if ok else '–')
                if ok and p:
                    V.add(X.variant(c, s), 'kc_soma', 'ladder_four', p.get('change'), panel='validation2')
            L.append(f'| {CALYX_NAMES[c]} | four | ' + ' | '.join(cells) + ' |')
    return L + ['']


def step2_table(d: dict, step1: Optional[dict], V: Values) -> List[str]:
    """Step 2: D1.8 and D5.8 on the male calyces, their class on the real wiring (step 1's score) beside their class on
    the four wirings (the final score), with the megamat change, the held-out panels that fail and validation2's
    change."""
    cl = d['classes'].get('classes') or {}
    c1 = (step1 or {}).get('classes') or {}
    L = ["The doses on a male calyx's four wirings (step 2 adds the three permutations), beside the real wiring alone "
         "(step 1's score), and A1 on the same four wirings. Changes from the calyx's family, mean over the wirings.",
         '',
         '| calyx and setting | wirings | class | megamat change (bar) | validation2 | APL16, APL on | other views '
         'worse | held-out panels worse |', '|---|---|---|---:|---:|---:|---|---|']
    for c in X.MALE:
        for s in ('D1.8', 'A1', 'D5.8'):
            v = X.variant(c, s)
            recs = []
            if s != 'A1' and v in c1:
                recs.append(('real (step 1)', c1[v]))
            if v in cl:
                recs.append(('', cl[v]))
                if cl[v].get('on_wirings_held'):
                    recs.append(('held', cl[v]['on_wirings_held']))
            for tag, r in recs:
                m = r.get('megamat') or {}
                names = _names
                cls = r.get('cls')
                if cls == 'no class' and r.get('reasons'):
                    cls = f"no class ({'; '.join(r['reasons'])})"
                lab = label(v) + {'real (step 1)': " (step 1's score)", '': '',
                                  'held': ', on the wirings whose holds held (a description)'}[tag]
                L.append(f"| {lab} | {', '.join(r.get('wirings') or [])} | {cls} | {_ch(m)} "
                         f"({_f(r.get('bar'))}) | {_ch((r.get('panels') or {}).get('validation2/validation2/on'))} | "
                         f"{_ch((r.get('panels') or {}).get('apl16/apl16/on'))} | {names(r.get('guard_worse') or [])} | "
                         f"{names(r.get('panel_worse') or [])} |")
                if tag == 'real (step 1)':
                    V.add(v, 'stage', 'change_step1', m.get('change'))
    return L + ['']


def shares_table(d: dict, V: Values) -> List[str]:
    """E143's resampled reading: per class the share of resamples giving the scored class and the classes' shares;
    the derived classes and the outcomes (a) to (d)."""
    sh, e = d['shares'], d['e143']
    if not len(sh):
        return []
    L = ["E143's 1,000 resamples of the flies of every megamat stage and held-out panel (the model read once against "
         'all the flies; a description, gating nothing). The share of resamples giving each class.', '',
         '| calyx and setting | wirings | class as scored | against all the flies | share giving the scored class | '
         'the classes\' shares |', '|---|---|---|---|---:|---|']
    for r in sh.itertuples():
        shares = json.loads(r.shares) if isinstance(r.shares, str) else {}
        txt = '; '.join(f'{k} {v:.3f}' for k, v in sorted(shares.items(), key=lambda kv: -kv[1]))
        lab = label(r.arm) + (', on the wirings whose holds held' if '@held' in r.arm else '')
        L.append(f'| {lab} | {r.wirings} | {r.cls_scored} | {r.cls_all_flies} | {_f(r.share_scored)} | {txt} |')
    L += ['', 'The derived classes and the outcomes through the same resamples: the value as scored, against all the '
          'flies, and the share of resamples on which it is true.', '',
          '| calyx | statistic | as scored | against all the flies | share true |', '|---|---|---|---|---:|']
    names = dict(a1_holds='A1 holds', a1_megamat_falls='megamat falls with A1', b_holds='B holds',
                 rise_needed='the rise is needed')
    for c, x in (e.get('derived') or {}).items():
        for k, y in x.items():
            L.append(f"| {c} | {names.get(k, k)} | {y.get('scored')} | {y.get('all_flies')} | "
                     f"{_f(y.get('share_true'))} |")
            V.add(f'{c}_derived', k, 'share_true', y.get('share_true'), half='resampled')
    for k, y in (e.get('outcomes') or {}).items():
        L.append(f"| – | outcome ({k}) | {y.get('scored')} | {y.get('all_flies')} | {_f(y.get('share_true'))} |")
        V.add('outcome', k, 'share_true', y.get('share_true'), half='resampled')
    sc = e.get('shortcut_check') or {}
    L += ['', f"E143's read-once check on these reads (the model read once against all the flies, a half changing only "
          f"the reference, against E136's own reading of the held-out half): largest difference of a megamat stage "
          f"error {_f(sc.get('megamat_stage_max_abs_diff'))}, of any value {_f(sc.get('any_max_abs_diff'))}, against "
          f"E143's tolerance {_f(sc.get('tol'), '.2f')} (within it on every megamat stage error: "
          f"{sc.get('megamat_stage_within_tol')}).", '']
    V.add('shortcut', 'stage', 'max_abs_diff', sc.get('megamat_stage_max_abs_diff'), half='resampled')
    V.add('shortcut', 'any', 'max_abs_diff', sc.get('any_max_abs_diff'), half='resampled')
    return L


def shortcut_table(d: dict, V: Values) -> List[str]:
    """E143's read-once check by calyx: the largest difference of the megamat stage error and of any statistic."""
    s = d['shortcut']
    if not len(s):
        return []
    s = s.assign(calyx=s.variant.str.split('_').str[0])
    L = ['The read-once check by calyx: the largest difference between the shortcut and the score\'s reading, on the '
         "megamat stage error and on any statistic (E143's tolerance 0.02).", '',
         '| calyx | megamat stage error | any statistic | statistics beyond 0.02 |', '|---|---:|---:|---:|']
    for c, x in s.groupby('calyx', sort=False):
        st = x[x.statistic == 'stage']
        L.append(f"| {c} | {_f(st['diff'].abs().max())} | {_f(x['diff'].abs().max())} | "
                 f"{int((~x.within_tol.astype(bool)).sum())} of {len(x)} |")
        V.add(f'{c}_shortcut', 'stage', 'max_abs_diff', st['diff'].abs().max(), half='resampled')
    return L + ['']


def claw_default_table(d: dict, V: Values) -> List[str]:
    """H against the hemibrain's 'first' real wiring: each setting's megamat stages and held-out panels, the difference
    and the spread; the class on one wiring under each."""
    cd = d['classes'].get('claw_default') or {}
    keys = ('stage',) + tuple(STAGES) + PANEL_KEYS
    L = ["The claw default: H ('argmax' claws, the hemibrain's real wiring) minus the same setting on the 'first' real "
         "wiring (E141's and E139's saved reads), strength by odor against the held-out half. A difference above the "
         'spread changes it (the entry\'s class); the last columns give the class on one wiring under each.', '',
         '| setting | ' + ' | '.join('megamat stage error' if k == 'stage' else STAGE_NAMES.get(k, PANEL_NAMES.get(k, k))
                                    for k in keys) + ' | class on H | class on first |',
         '|---|' + '---:|' * len(keys) + '---|---|',
         '| the spread | ' + ' | '.join(_f(d['classes']['spreads'][k]) for k in keys) + ' | – | – |']
    for s in X.SETTINGS:
        r = cd.get(s)
        if not r or not r.get('complete'):
            continue
        cells = []
        for k in keys:
            x = r['diffs'].get(k)
            cells.append('–' if x is None else _f(x['diff'], '+.3f') + (' (beyond)' if x['differs'] else ''))
            if x is not None:
                st, pan, cond = (k, 'megamat', 'on') if '/' not in k else ('kc_soma', *k.split('/')[1:])
                V.add(X.variant('H', s), st, 'diff_first', x['diff'], panel=pan, condition=cond)
        L.append(f"| {SETTING_NAMES[s]} | " + ' | '.join(cells) +
                 f" | {r.get('class_on_H', '–')} | {r.get('class_on_first', '–')} |")
    return L + ['']


# ---------------------------------------------------------------------------
# the holds, realism, the scales, where APL's inhibition is made, the APL effect

def holds_table(d: dict, V: Values) -> List[str]:
    """The reads at a bound and the doses' megamat reads on L (APL's factor against its x10 bound, the threshold's,
    the response rates reached against their targets and the absolute scales)."""
    r = d['reads']
    x = r[(r.calyx == 'L') & r.setting.isin(['A1', 'D5.8', 'Dabs']) & (r.panel == 'megamat')].copy()
    x['order'] = x.setting.map({'A1': 0, 'D5.8': 1, 'Dabs': 2})
    x = x.sort_values(['order', 'wiring'])
    L = ["The left calyx's megamat holds at the highest doses (E136's bracket: APL's scale from ×0.1 to ×10 of the "
         'start, the start being the calyx\'s fit with the setting\'s factors). Both response rates are held to within '
         '0.5 % of their targets (APL on 0.0915, APL off the calyx\'s fit\'s ratio).', '',
         '| setting | wiring | APL factor | at the bound | threshold factor | response rate, APL on (target) | APL off '
         '(target) | APL→KC | KC–KC | somatic threshold (Hz) |', '|---|---|---:|---|---:|---:|---:|---:|---:|---:|']
    for t in x.itertuples():
        L.append(f'| {SETTING_NAMES[t.setting]} | {t.wiring} | {_f(t.apl_factor, ".2f")} | '
                 f"{'yes' if str(t.at_bound_apl) == 'True' else 'no'} | {_f(t.thr_factor, '.3f')} | "
                 f'{_f(t.rate_on, ".4f")} ({_f(t.frac_on_target, ".4f")}) | {_f(t.rate_off, ".4f")} '
                 f'({_f(t.frac_off_target, ".4f")}) | {_f(t.wAPLKC_scale, ".3f")} | {_f(t.wKCKC_scale, ".1f")} | '
                 f'{_f(t.fixed_thr, ".1f")} |')
        for col in ('apl_factor', 'thr_factor', 'rate_on', 'rate_off', 'wAPLKC_scale', 'wKCKC_scale', 'fixed_thr'):
            V.add(f'L_{t.setting}', 'reads', col, getattr(t, col), condition=t.wiring)
    return L + ['']


def realism_table(d: dict, V: Values) -> List[str]:
    """The four bands on the APL16 input per calyx and setting (range over its wirings), with the KC data's values."""
    a = d['apl']
    t = W.data_targets()
    L = ['The APL16 input (held out), range over each setting\'s wirings. Bands: response rate with APL on 0.045 to '
         '0.085; APL-off / APL-on ratio 1.6 to 2.9; response rate ρ with APL on at least 0.37 (the reference model\'s '
         '0.47 less 0.10); odors per KC within 0.155 (APL on) and 0.223 (APL off) of the flies (total variation '
         'distance).', '',
         '| model | wirings | response rate, APL on | APL off | ratio | response rate ρ, APL on | APL off | odors per KC, '
         'APL on | APL off | largest spike count, APL on / off | flagged |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    if t:
        g = lambda k, f='.3f': _f(t.get(k), f)
        L.append(f"| the KC data (control flies / APL output blocked (TNT), noise criterion) | – | "
                 f"{g('noise_ctrl_frac')} | {g('noise_tnt_frac')} | {g('noise_rate_ratio_fraction', '.2f')} | – | – | "
                 f"band {g('noise_ctrl_odors_per_kc_tv_max')} | band {g('noise_tnt_odors_per_kc_tv_max')} | – | – |")
    for v in variants():
        x = a[a.variant == v]
        if not len(x):
            continue
        flag = x.realism_flag.fillna('').astype(str)
        named = '; '.join(f'{w}: {fl}' for w, fl in zip(x.wiring, flag) if fl) or 'none'
        L.append(f'| {label(v)} | {len(x)} | {_rng(x.rate_on)} | {_rng(x.rate_off)} | {_rng(x.rate_ratio, ".2f")} | '
                 f'{_rng(x.rate_rho_on, ".2f")} | {_rng(x.rate_rho_off, ".2f")} | {_rng(x.opk_on_tv)} | '
                 f'{_rng(x.opk_off_tv)} | {_rng(x.max_count_on, ".0f")} / {_rng(x.max_count_off, ".0f")} | {named} |')
        for col in ('rate_on', 'rate_off', 'rate_ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv',
                    'max_count_on', 'max_count_off'):
            V.add(v, 'kc_soma', f'{col}_min', x[col].min(), panel='apl16', condition='realism', half='all')
            V.add(v, 'kc_soma', f'{col}_max', x[col].max(), panel='apl16', condition='realism', half='all')
        V.add(v, 'kc_soma', 'n_flagged', int((flag != '').sum()), panel='apl16', condition='realism', half='all')
    return L + ['']


def scales_table(d: dict, V: Values) -> List[str]:
    """The megamat operating point and the absolute scales the holds reached (range over each setting's wirings), with
    the largest megamat spike counts."""
    r = d['reads']
    L = ["The megamat reads' operating point and the absolute scales the holds reached, range over each setting's "
         "wirings (the flies' response rate with APL on 0.0915).", '',
         '| model | wirings | response rate, APL on | APL off | APL→bouton | APL→KC | KC–KC | somatic threshold (Hz) | '
         'PN normalization, additive | largest spike count, APL on / off |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for v in variants():
        x = r[(r.variant == v) & (r.panel == 'megamat')]
        if not len(x):
            continue
        add = x['pn_inhadd'] if 'pn_inhadd' in x and x['pn_inhadd'].notna().any() else None
        L.append(f'| {label(v)} | {len(x)} | {_rng(x.rate_on, ".4f")} | {_rng(x.rate_off, ".4f")} | '
                 f'{_rng(x.wAPLPN_scale, ".4f")} | {_rng(x.wAPLKC_scale, ".3f")} | {_rng(x.wKCKC_scale, ".1f")} | '
                 f'{_rng(x.fixed_thr, ".1f")} | {_rng(add, ".0f") if add is not None else "31.4 (default)"} | '
                 f'{_rng(x.max_count_on, ".0f")} / {_rng(x.max_count_off, ".0f")} |')
        for col in ('rate_on', 'rate_off', 'wAPLPN_scale', 'wAPLKC_scale', 'wKCKC_scale', 'fixed_thr', 'max_count_on',
                    'max_count_off'):
            V.add(v, 'reads', f'{col}_min', x[col].min(), condition='scales', half='all')
            V.add(v, 'reads', f'{col}_max', x[col].max(), condition='scales', half='all')
    return L + ['']


def sites_table(d: dict, V: Values) -> List[str]:
    """Where APL's inhibition is made (score_sites.csv): the share of the claws' summed drive APL removes and the share
    of that loss made at the PN boutons, mean over the wirings (range), beside the calyx's share of APL's calyx synapses
    at the boutons (E130's)."""
    s = d['sites']
    L = ['Where APL\'s inhibition is made, window means, mean over each setting\'s wirings (range). "Removed" is the '
         "share of the claws' summed drive APL removes; \"at the boutons\" the share of that loss made at the PN "
         "boutons, beside the calyx's share of APL's calyx synapses that are onto boutons (E130's step 0).", '',
         '| model | wirings | connectome\'s share at the boutons | megamat: removed | at the boutons | APL16: removed | '
         'at the boutons |', '|---|---:|---:|---:|---:|---:|---:|']
    for v in variants():
        x = s[s.variant == v]
        if not len(x):
            continue
        mm, aa = x[x.panel == 'megamat'], x[x.panel == 'apl16']
        L.append(f'| {label(v)} | {len(mm)} | {_f(x.anatomical_bouton_share.iloc[0])} | {_mr(mm.drive_removed)} | '
                 f'{_mr(mm.bouton_share)} | {_mr(aa.drive_removed)} | {_mr(aa.bouton_share)} |')
        for pan, y in (('megamat', mm), ('apl16', aa)):
            if len(y):
                V.add(v, 'sites', 'bouton_share_mean', y.bouton_share.mean(), panel=pan, half='all')
                V.add(v, 'sites', 'drive_removed_mean', y.drive_removed.mean(), panel=pan, half='all')
        V.add(v, 'sites', 'anatomical_bouton_share', x.anatomical_bouton_share.iloc[0], half='all')
    return L + ['']


def apl_effect_table(d: dict, V: Values) -> List[str]:
    """The APL effect, held out (score_apl16.csv's price columns: the calibrated readout against all recordings, z8),
    mean over each setting's wirings: the loss, the mean APL effect, the pattern part with its amplitude and shape parts
    and the spread ratio, ρ and Pearson r on the per-pair APL effect."""
    a = d['apl']
    L = ['The APL effect on odor-odor distances (a held-out test; gates nothing), under the calibrated readout against '
         "all recordings, mean over each setting's wirings. The flies' mean APL effect is +0.178; the floor of the "
         'pattern part about 0.029.', '',
         '| model | wirings | loss | mean APL effect | pattern part | amplitude part | shape part | spread ratio | ρ on '
         'the per-pair APL effect | Pearson r |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    cols = ('price_z8_loss', 'price_mean_effect', 'price_z8_pattern', 'price_z8_amplitude', 'price_z8_shape',
            'price_z8_sd_ratio', 'price_z8_rho', 'price_z8_pearson')
    for v in variants():
        x = a[a.variant == v]
        if not len(x):
            continue
        L.append(f'| {label(v)} | {len(x)} | ' + ' | '.join(_f(x[c].mean(), '+.3f' if c == 'price_mean_effect' else
                                                                  '.2f' if c in ('price_z8_sd_ratio', 'price_z8_rho',
                                                                                 'price_z8_pearson') else '.4f')
                                                         for c in cols) + ' |')
        for c in cols:
            V.add(v, 'apl_effect', c, x[c].mean(), panel='apl16', condition='on_off', half='all')
    return L + ['']


def without_aa_table(d: dict, V: Values) -> List[str]:
    """The megamat somata's strength error with and without aa's entry (mean over each setting's wirings)."""
    w = d['without']
    if not len(w):
        return []
    L = ["The megamat KC somata's strength error (held-out half) with and without aa's entry, mean over each setting's "
         "wirings, and aa's share of the error.", '',
         '| model | wirings | with aa | without aa | aa\'s share |', '|---|---:|---:|---:|---:|']
    for v in variants():
        x = w[w.variant == v]
        x = x[x['row'].isin(X.WIRINGS_OF[v.split('_')[0]])]
        if not len(x):
            continue
        L.append(f'| {label(v)} | {len(x)} | {_f(x.view_error.mean())} | {_f(x.view_error_without.mean())} | '
                 f'{_f(x.view_share.mean())} |')
        V.add(v, 'kc_soma', 'with_aa', x.view_error.mean())
        V.add(v, 'kc_soma', 'without_aa', x.view_error_without.mean())
    return L + ['']


def e137_table(d: dict, V: Values) -> List[str]:
    """The mixtures under E137's readings of the mixture flies (strength by odor; mean over each setting's wirings)."""
    e = d['e137']
    if not len(e):
        return []
    e = e[e['views'] != 'all']
    combos = sorted({(p, r) for p, r in zip(e.panel, e.reading)})
    nm = {'kmix': 'kiwi mixtures', 'cmix': 'control mixtures'}
    rn = W.READING_NAMES
    L = ["The mixtures' KC somata under E137's readings of the mixture flies (each reading's own units and trial noise), "
         "strength by odor, mean over each setting's wirings. Beside; gates nothing.", '',
         '| model | ' + ' | '.join(f'{nm.get(p, p)}, {rn.get(r, r)}' for p, r in combos) + ' |',
         '|---|' + '---:|' * len(combos)]
    for v in variants():
        x = e[e.variant == v]
        if not len(x):
            continue
        cells = []
        for p, r in combos:
            y = x[(x.panel == p) & (x.reading == r)]
            cells.append(_f(y.error.mean()) if len(y) else '–')
            if len(y):
                V.add(v, 'kc_soma', f'e137_{r}', y.error.mean(), panel=p)
        L.append(f'| {label(v)} | ' + ' | '.join(cells) + ' |')
    return L + ['']


# ---------------------------------------------------------------------------
# step 0, the prediction, the reads

def matrix_table(d: dict, V: Values) -> List[str]:
    """Step 0's KC-KC matrices (step0_matrices.csv): the statistics the within-type-pair shuffle changes."""
    m = d['step0_m'].set_index('matrix')
    rows = (('row_sum_cv', 'unevenness: CV of the KCs\' summed input', '.3f'),
            ('reciprocal_share', 'reciprocity: nonzero entries whose transpose is nonzero', '.3f'),
            ('shared_glomerulus_ratio', 'weight between KCs that share an input glomerulus, over their share of pairs '
                                        '(1: no preference)', '.3f'),
            ('within_type_weight_share', 'weight within KC type (share of the total)', '.3f'),
            ('inputs_per_kc', 'KC–KC inputs per KC (nonzero entries)', '.2f'))
    cols = (('hemibrain', "the hemibrain ('first')"), ('hemibrain_SH', 'the hemibrain shuffled within type pairs'),
            ('L', 'L'), ('R', 'R'))
    L = ['The KC–KC matrices (step 0, each normalized to a mean row sum of 1).', '',
         '| statistic | ' + ' | '.join(n for _, n in cols) + ' |', '|---|' + '---:|' * len(cols)]
    for key, name, fmt in rows:
        L.append(f'| {name} | ' + ' | '.join(_f(m.loc[c, key], fmt) for c, _ in cols) + ' |')
        for c, _ in cols:
            V.add(f'matrix_{c}', 'step0', key, m.loc[c, key], half='all')
    return L + ['']


def prediction_table(d: dict, V: Values) -> List[str]:
    """The planner's prediction per calyx (prediction.json, written before any read) beside the outcome."""
    p = (d['prediction'].get('prediction') or {})
    cal = d['classes'].get('calyces') or {}
    cl = d['classes'].get('classes') or {}
    L = ["The planner's prediction per calyx (written 11:10, before any read) and the score.", '',
         '| calyx | megamat falls with A1, predicted | scored | A1 holds, predicted | A1\'s class | what was predicted to '
         'fail it | what failed it |', '|---|---|---|---|---|---|---|']
    names = lambda xs: _names(xs) if xs else ''
    for c in X.CALYCES:
        a = cl.get(X.variant(c, 'A1')) or {}
        failed = '; '.join(x for x in (names(a.get('stage_worse') or []),
                                        ('other views at ' + names(a.get('guard_worse'))) if a.get('guard_worse') else '',
                                        names(a.get('panel_worse') or [])) if x) or 'nothing'
        L.append(f"| {c} | {p.get(c, {}).get('megamat_falls')} | {'yes' if cal.get(c, {}).get('a1_megamat_falls') else 'no'} "
                 f"| {p.get(c, {}).get('a1_holds')} | {a.get('cls')} | {p.get(c, {}).get('would_fail_it')} | {failed} |")
    return L + ['']


def outcomes_table(d: dict, V: Values) -> List[str]:
    o = d['classes'].get('outcomes') or {}
    L = ["The entry's outcomes as scored.", '', '| outcome | met |', '|---|---|']
    for k in 'abcdef':
        L.append(f'| ({k}) | {o.get(k)} |')
    return L + ['']


def lanes_table(d: dict, V: Values) -> List[str]:
    """The reads by host and lane: count, step, commits, dirty."""
    r = d['reads']
    step = lambda t: '2' if (t.split('_')[0] in X.MALE and t.split('_')[1] in ('D1.8', 'D5.8')
                             and not t.endswith('_real')) else '1'
    r = r.assign(step=r.tag.map(step))
    L = ['The reads by lane (reads_<host>_<lane>.csv).', '',
         '| host | lane | reads | steps | queue | pyolfsysm | hong2p | any checkout dirty |',
         '|---|---|---:|---|---|---|---|---|']
    for (h, ln), x in r.groupby(['host', 'lane'], sort=True):
        dirty = any(str(v) == 'True' for col in ('queue_dirty', 'pyolfsysm_dirty', 'al_analysis_dirty', 'hong2p_dirty')
                    for v in x[col])
        L.append(f"| {h} | {ln} | {len(x)} | {', '.join(sorted(set(x.step)))} | "
                 f"{', '.join(sorted({str(c)[:7] for c in x.queue_commit}))} | "
                 f"{', '.join(sorted({str(c)[:7] for c in x.pyolfsysm_commit}))} | "
                 f"{', '.join(sorted({str(c)[:7] for c in x.hong2p_commit}))} | {'yes' if dirty else 'no'} |")
    L.append(f'| all | – | {len(r)} ({r.tag.nunique()} tags) | 1, 2 | – | – | – | – |')
    V.add('reads', 'all', 'n', len(r), half='all')
    return L + ['']


def main(score_dir: Optional[Path] = None) -> Path:
    score_dir = Path(score_dir) if score_dir is not None else X.out_dir()
    d = load(score_dir)
    step1 = W._read_json(score_dir / 'step1_score' / 'classes.json')
    V = Values()
    prov = d['prov']
    head = [f'# E144 write-up tables', '',
            f"From results/E144_a1_malecns_calyces/ (the score at queue {str(prov.get('queue_commit'))[:7]}, "
            f"{prov.get('time')}; E143's resampled reading at queue {str(d['prov_r'].get('queue_commit'))[:7]}, "
            f"{d['prov_r'].get('time')}); written by analyses/E144_writeup_tables.py.", '']
    parts = [('The stage table', stage_table), ("The calyx's family against the hemibrain's wirings", inside_table),
             ('The families over all views', all_views_table),
             ('The megamat stage error against each set of flies', halves_table), ('The classes', classes_table),
             ('The megamat stages', stage_changes_table), ('The held-out panels', panels_table),
             ('The dose ladder at validation2', ladder_table),
             ('Step 2: the doses on the permutations', lambda d_, V_: step2_table(d_, step1, V_)),
             ('The resampled reading', shares_table), ('The read-once check by calyx', shortcut_table),
             ('The claw default', claw_default_table), ('The holds at the bound', holds_table),
             ('Realism: the bands', realism_table), ('Realism: the operating point and the scales', scales_table),
             ("Where APL's inhibition is made", sites_table), ('The APL effect, held out', apl_effect_table),
             ('The megamat somata without aa', without_aa_table), ("The mixtures under E137's readings", e137_table),
             ('Step 0: the KC-KC matrices', matrix_table), ('The prediction and the outcome', prediction_table),
             ('The outcomes', outcomes_table), ('The reads', lanes_table)]
    lines = head
    for title, fn in parts:
        lines += [f'## {title}', ''] + fn(d, V)
    out = score_dir / 'writeup_tables.md'
    out.write_text('\n'.join(lines) + '\n')
    V.frame().to_csv(score_dir / 'writeup_values.csv', index=False)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--score-dir', type=Path, default=None)
    a = ap.parse_args()
    print(main(a.score_dir))
