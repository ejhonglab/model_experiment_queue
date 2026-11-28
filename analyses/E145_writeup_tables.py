"""The tables of E145's write-up (the fit of all megamat stages made again on the other wirings and the other halves of
the flies), from the files its parts wrote to results/E145_stage_fit_replication/: the score's classes.json and
score_*.csv; E143's resampled reading (e143_*.csv, e143_intervals.json); step 0's step0_cells.csv and prediction.json;
each cell's fit record, history, returns, rescore, consistency and reads in its folder (<cell>/); and E136's own A1 (its
fit, history, rescore, consistency and reads on the adopted measure), with E134's flies tables (the floors), E94's
connectome share and the data-target registry (the KC data's values beside the realism bands). No simulation and no fit.

Writes writeup_tables.md and writeup_values.csv (the same numbers, long: arm, stage, panel, condition, half, statistic,
value; the arm is the cell's name, 'family@<cell>' for the unrefitted family where the cell is read, 'family' for it on
the eight wirings) beside the score's files.

Usage (from the repository root): ./venv/bin/python -m analyses.E145_writeup_tables [--score-dir <dir>]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from analyses import E136_writeup_tables as W
from experiments import E136_stage_population_fit as F0
from experiments import E145_stage_fit_replication as E

STAGES = W.STAGES
STAGE_NAMES = W.STAGE_NAMES
CELLS = ('A1',) + tuple(E.FIT_CELLS)
CELL_NAMES = {'A1': 'A1 (E136): S, half a', 'X': 'X, the exchange: T, half b', 'F': 'F, the other flies: S, half b',
              'W': 'W, the other wirings: T, half a'}
FACTORS = ('aplpn', 'kckc', 'inhadd', 'inhsc')
FACTOR_NAMES = {'aplpn': 'APL→bouton', 'kckc': 'KC–KC', 'inhadd': 'PN normalization, additive term (pn_inhadd)',
                'inhsc': 'gain on pn_inhsc'}
PANEL_KEYS = tuple('/'.join(k) for k in F0.HOLDS_PANELS)
BESIDE_KEYS = tuple('/'.join(k) for k in F0.BESIDE_PANELS)
PANEL_NAMES = W.PANEL_NAMES
OUTCOMES = {'a': 'the fall, the hold and all three factors repeat',
            'b': 'the fall and the two factors before the calyx repeat, KC–KC does not (two or more cells megamat only)',
            'c': 'the hold repeats at other factors', 'd': 'the cells disagree',
            'e': 'the fall does not repeat in a cell', 'f': 'a fit ends at a bound'}
_f, _mr = W._f, W._mr


def _rng(x, fmt: str = '.3f') -> str:
    """A range over wirings as 'lo to hi', or one value where both print alike."""
    x = pd.Series(x, dtype=float).dropna()
    if not len(x):
        return '–'
    lo, hi = _f(x.min(), fmt), _f(x.max(), fmt)
    return lo if lo == hi else f'{lo} to {hi}'


def _names(xs: Sequence[str]) -> str:
    """Stages and panels by name, '; '-separated; 'none' for none."""
    return '; '.join(STAGE_NAMES.get(x, PANEL_NAMES.get(x, x)) for x in xs) or 'none'


def _ch(x: Optional[Mapping], fmt: str = '+.3f') -> str:
    return _f(x.get('change'), fmt) if x and x.get('complete') else '–'


def variant(cell: str) -> str:
    """A cell's fitted point as the score names it."""
    return f'E145_{cell}'


def factor_value(f: str, z: float) -> float:
    """A searched factor's value at coordinate `z`: the additive term as pn_inhadd itself, the others as factors."""
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    return float(PN_INHADD_DEFAULT * np.exp(z)) if f == 'inhadd' else float(np.exp(z))


def _fmt_factor(f: str, v) -> str:
    return _f(v, '.0f') if f == 'inhadd' else ('×' + _f(v, '.3f') if v is not None else '–')


# ---------------------------------------------------------------------------
# inputs

def _rd(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, float_precision='round_trip') if p.exists() else pd.DataFrame()


def cell_dir(cell: str, score_dir: Path) -> Path:
    """Where a cell's fit and reads are: E136's results folder on the adopted measure for A1, <score_dir>/<cell>/
    otherwise."""
    if cell == 'A1':
        with E.e136() as F:
            return F.out_dir()
    return score_dir / cell


def cell_files(cell: str, score_dir: Path) -> dict:
    """A cell's fit record, history, returns, consistency, rescore and reads (A1's reads only from E136's)."""
    p = cell_dir(cell, score_dir)
    reads = [pd.read_csv(x) for x in sorted(p.glob('reads_*.csv'))]
    r = pd.concat(reads, ignore_index=True) if reads else pd.DataFrame()
    if len(r) and 'arm' in r:
        r = r[r.arm.astype(str) == E.ARM]
    return dict(dir=p, fit=W._read_json(p / f'fit_{E.ARM}.json'), hist=_rd(p / f'fit_history_{E.ARM}.csv'),
                det=_rd(p / f'determined_{E.ARM}.csv'), cons=_rd(p / 'consistency_strength.csv'),
                rescore=_rd(p / f'rescore_{E.ARM}.csv'), reads=r)


def load(score_dir: Path) -> dict:
    """Everything the tables read."""
    d = dict(classes=W._read_json(score_dir / 'classes.json'), mega=_rd(score_dir / 'score_megamat.csv'),
             panels=_rd(score_dir / 'score_panels.csv'), apl=_rd(score_dir / 'score_apl16.csv'),
             sites=_rd(score_dir / 'score_sites.csv'), without=_rd(score_dir / 'score_without_aa.csv'),
             e137=_rd(score_dir / 'score_e137.csv'), own=_rd(score_dir / 'score_own_wirings.csv'),
             shares=_rd(score_dir / 'e143_class_shares.csv'), changes=_rd(score_dir / 'e143_intervals_changes.csv'),
             shortcut=_rd(score_dir / 'e143_shortcut_check.csv'), e143=W._read_json(score_dir / 'e143_intervals.json'),
             prediction=W._read_json(score_dir / 'prediction.json'), step0=_rd(score_dir / 'step0_cells.csv'),
             prov=W._read_json(score_dir / 'provenance_score.json'),
             prov_r=W._read_json(score_dir / 'provenance_resampled.json'))
    d['cells'] = {c: cell_files(c, score_dir) for c in CELLS}
    e = F0.e134_dir()
    d['flies_views'] = _rd(e / 'flies_views.csv')
    d['no_step'] = _rd(e / 'flies_no_step.csv')
    return d


def recs(d: dict) -> Dict[str, dict]:
    return (d['classes'] or {}).get('cells') or {}


# ---------------------------------------------------------------------------
# the cells and the fits

def cells_table(d: dict, V: W.Values) -> List[str]:
    """The four cells (step 0's table): what each searches and is fitted against, where it is read, the unrefitted
    family and A1's point there, and the cell's fitted point there (the score)."""
    s0 = d['step0'].set_index('cell') if len(d['step0']) else pd.DataFrame()
    rc = recs(d)
    L = ['The cells of wirings by halves of the flies. S is real, P1, P2, M1; T is P3, P4, M2, M3. Megamat stage error '
         'on strength by odor in one-fly units (lower is better), mean over the four read wirings, against the read '
         "half. A1's point is E136's (step 0, no simulation); the cell's point is its own fit's classed point (the "
         'score).', '',
         '| cell | searched | fitted against | read on | the family where read | A1\'s point there | '
         "the cell's point there |", '|---|---|---|---|---:|---:|---:|']
    for c in CELLS:
        cell = E.get_cell(c)
        r = s0.loc[c] if c in getattr(s0, 'index', []) else {}
        m = ((rc.get(c) or {}).get('cls') or {}).get('megamat') or {}
        fam, a1 = (r.get('family_where_read'), r.get('a1_where_read')) if len(r) else (None, None)
        L.append(f'| {CELL_NAMES[c]} | {cell.searched} | half {cell.fit_half} | {cell.read_on}, half {cell.read_half} | '
                 f"{_f(fam)} | {_f(a1)} | {_f(m.get('arm')) if m.get('complete') else '–'} |")
        V.add(f'family@{c}', 'stage', 'step0_where_read', fam)
        V.add(f'A1@{c}', 'stage', 'step0_where_read', a1)
        V.add(c, 'stage', 'where_read', m.get('arm') if m.get('complete') else None)
    return L + ['']


def _start_loss(h: pd.DataFrame) -> Optional[float]:
    if not len(h) or 'loss_sim' not in h:
        return None
    x = h[h['eval'] == 0]
    return float(x.loss_sim.iloc[0]) if len(x) else None


def own_consistency(rescore: pd.DataFrame, col: str = 'stage_error', measure: str = 'strength') -> Optional[float]:
    """A fit's own consistency between the halves: the Spearman ρ over its evaluations between `col` against the fit
    half and against the held-out half (E136's `consistency`, on this fit's rescore alone). E136's consistency file
    pools every arm's evaluations (483 for its five arms), so A1's own value is read from A1's rescore, not from it."""
    from scipy import stats
    r = rescore[rescore.measure == measure] if len(rescore) and 'measure' in rescore else rescore
    if not len(r) or col not in r:
        return None
    a = r[r.half == 'fit'].set_index('eval')[col]
    b = r[r.half == 'held'].set_index('eval')[col]
    idx = a.index.intersection(b.index)
    return float(stats.spearmanr(a.loc[idx], b.loc[idx])[0]) if len(idx) > 2 else None


def fits_table(d: dict, V: W.Values) -> List[str]:
    """Per cell, the fit: evaluations, start and best loss (on the searched wirings against the fit half), the best
    point's factors, any bound, and the consistency of the strength stage error between the halves over the
    evaluations (rho; E136's yardstick 0.8)."""
    L = ['Each cell\'s fit (E136\'s A1 in its search mode: 140 evaluations with one restart; the loss is the megamat stage '
         'error on strength on the searched wirings against the fit half, lower is better). "Consistency" is ρ over the '
         "fit's own evaluations between their megamat stage errors against the fit half and against the read half "
         "(E136's yardstick 0.8; E136's own 0.823 pools its five arms' 483 evaluations).", '',
         '| cell | evaluations | start | best (evaluation) | APL→bouton | KC–KC | pn_inhadd | gain on pn_inhsc | '
         'at a bound | consistency ρ |', '|---|---:|---:|---:|---:|---:|---:|---:|---|---:|']
    for c in CELLS:
        x = d['cells'][c]
        f = x['fit']
        if not f:
            continue
        bz = f.get('best_z') or {}
        vals = {k: factor_value(k, bz[k]) for k in FACTORS if k in bz}
        rho = own_consistency(x['rescore'])
        best = None
        h = x['hist']
        if len(h) and f.get('best_eval') is not None:
            y = h[h['eval'] == f['best_eval']]
            best = float(y.loss_sim.iloc[0]) if len(y) else None
        start = _start_loss(h)
        L.append(f"| {CELL_NAMES[c]} | {f.get('n_evals')} | {_f(start)} | {_f(best)} ({f.get('best_eval')}) | "
                 + ' | '.join(_fmt_factor(k, vals.get(k)) for k in FACTORS)
                 + f" | {', '.join(f.get('at_bound') or []) or 'none'} | {_f(rho)} |")
        V.add(c, 'fit', 'n_evals', f.get('n_evals'), half=E.get_cell(c).fit_half)
        V.add(c, 'fit', 'start_loss', start, half=E.get_cell(c).fit_half)
        V.add(c, 'fit', 'best_loss', best, half=E.get_cell(c).fit_half)
        V.add(c, 'fit', 'consistency_rho', rho, half='both')
        for k, v in vals.items():
            V.add(c, 'fit', f'best_{k}', v, half=E.get_cell(c).fit_half)
    return L + ['']


def returns_table(d: dict, V: W.Values) -> List[str]:
    """Per cell, the single-factor returns (each searched factor returned alone to the start: the loss's rise, and
    whether it passes the 0.010 rule), the classed point's factors and the entry's repeat rules on them (the score's
    repeats)."""
    rc = recs(d)
    L = ["The single-factor returns and the classed point. A factor is determined when returning it alone to the start "
         "raises the loss by 0.010 or more; an undetermined factor returns to the start in the classed point. \"Repeats\" "
         "is the entry's rule: pn_inhadd between 562 and 1,778; APL→bouton determined, ×0.078 to ×0.31; KC–KC "
         "determined, ×1.8 or more; the gain on pn_inhsc is not read.", '',
         '| cell | factor | rise when returned | determined | classed value | repeats |', '|---|---|---:|---|---:|---|']
    for c in CELLS:
        f = d['cells'][c]['fit']
        det = f.get('determined') or {}
        rep = (rc.get(c) or {}).get('repeats') or {}
        cz = f.get('classed_z') or f.get('best_z') or {}
        for k in FACTORS:
            if k not in cz:
                continue
            dk = det.get(k) or {}
            rk = rep.get('inhsc_gain' if k == 'inhsc' else k) or {}
            val = factor_value(k, cz[k])
            r = rk.get('repeats')
            L.append(f"| {CELL_NAMES[c]} | {FACTOR_NAMES[k]} | {_f(dk.get('rise'), '+.4f')} | "
                     f"{W._yn(dk.get('determined')) if dk else '–'} | {_fmt_factor(k, val)} | "
                     f"{'not read' if k == 'inhsc' else W._yn(r) if r is not None else '–'} |")
            V.add(c, 'fit', f'rise_{k}', dk.get('rise'), half=E.get_cell(c).fit_half)
            V.add(c, 'fit', f'classed_{k}', val, half=E.get_cell(c).fit_half)
            if k != 'inhsc' and r is not None:
                V.add(c, 'fit', f'repeats_{k}', float(bool(r)), half='all')
        V.add(c, 'fit', 'classed_loss', f.get('classed_loss'), half=E.get_cell(c).fit_half)
    return L + ['']


def rescore_table(d: dict, V: W.Values) -> List[str]:
    """Per cell, the rescore: every evaluation on the searched wirings against the read half, the start, the best
    evaluation and the read half's lowest."""
    L = ["The rescore (every evaluation on the searched wirings, against the read half; a description). The fit's best "
         "evaluation is the fit half's lowest.", '',
         '| cell | the start, fit half / read half | the best, fit half / read half | the read half\'s lowest (evaluation) |',
         '|---|---:|---:|---:|']
    for c in CELLS:
        x = d['cells'][c]
        r, f = x['rescore'], x['fit']
        if not len(r) or not f:
            continue
        r = r[r.get('measure', 'strength') == 'strength'] if 'measure' in r else r
        fh, rh = (('fit', 'held') if set(r.half) <= {'fit', 'held'} else (E.get_cell(c).fit_half, E.get_cell(c).read_half))
        p = r.pivot_table(index='eval', columns='half', values='stage_error')
        if fh not in p or rh not in p:
            continue
        b = f.get('best_eval')
        lo = p[rh].idxmin()
        L.append(f'| {CELL_NAMES[c]} | {_f(p.loc[0, fh])} / {_f(p.loc[0, rh])} | '
                 f'{_f(p.loc[b, fh]) if b in p.index else "–"} / {_f(p.loc[b, rh]) if b in p.index else "–"} | '
                 f'{_f(p[rh].min())} ({lo}) |')
        V.add(c, 'stage', 'rescore_start', p.loc[0, rh], half='read')
        if b in p.index:
            V.add(c, 'stage', 'rescore_best', p.loc[b, rh], half='read')
        V.add(c, 'stage', 'rescore_lowest', p[rh].min(), half='read')
    return L + ['']


KCKC_BINS = ((0.0, 1.8, 'below ×1.8'), (1.8, 3.24, '×1.8 to ×3.24'), (3.24, 5.83, '×3.24 to ×5.83'),
             (5.83, 7.3, '×5.83 to ×7.3'), (7.3, np.inf, 'at ×7.3 (the bound)'))
SLOPE_WINDOWS = (0.005, 0.01, 0.02)


def kckc_bin(k: float, top: float) -> str:
    """The KCKC_BINS label of a KC–KC factor `k`, the last for one at the bound `top`."""
    if np.isclose(k, top, rtol=1e-9, atol=0):
        return KCKC_BINS[-1][2]
    return next(lab for a, b, lab in KCKC_BINS[:-1] if a <= k < b or (b == KCKC_BINS[-2][1] and k >= a))


def class_without_bound(rec: Mapping) -> Optional[str]:
    """A cell's class with the bound rule left out, from its scored record (E136's `classify`'s order on the record's
    fields: other reasons, then holds, megamat only, improves with panels unread, one stage at another's cost, strength
    at another view's cost, does not improve); None if the record has no megamat comparison."""
    r = rec.get('cls') or {}
    m = r.get('megamat') or {}
    if not m:
        return None
    if [x for x in (r.get('reasons') or []) if not x.startswith('at a bound')]:
        return 'no class'
    imp, pw, pm = r.get('improves'), r.get('panel_worse') or [], r.get('panel_missing') or []
    if imp and not pw and not pm:
        return 'holds'
    if imp and pw:
        return 'megamat only'
    if imp:
        return 'improves (held-out panels not all read)'
    if m.get('improves') and r.get('stage_worse'):
        return "one stage at another's cost"
    if m.get('improves') and r.get('guard_worse'):
        return "strength at another view's cost"
    return 'does not improve'


def local_slope(h: pd.DataFrame, window: float, coord: str = 'z_kckc') -> dict:
    """The fit half's loss regressed on the four searched coordinates (least squares, with an intercept) over the
    evaluations whose loss is within `window` of the best: the number of evaluations, the coefficient on `coord` (the
    loss's change per e-fold of the factor) and its standard error."""
    cols = [f'z_{k}' for k in FACTORS]
    x = h[h.loss_sim <= h.loss_sim.min() + window]
    A = np.column_stack([np.ones(len(x)), x[cols].to_numpy(float)])
    y = x.loss_sim.to_numpy(float)
    if len(y) <= A.shape[1]:
        return dict(n=len(y), slope=np.nan, se=np.nan)
    beta = np.linalg.lstsq(A, y, rcond=None)[0]
    s2 = float(((y - A @ beta) ** 2).sum() / (len(y) - A.shape[1]))
    se = float(np.sqrt(s2 * np.linalg.pinv(A.T @ A)[1 + cols.index(coord), 1 + cols.index(coord)]))
    return dict(n=int(len(y)), slope=float(beta[1 + cols.index(coord)]), se=se)


def kckc_table(d: dict, V: W.Values, cell: str = 'X') -> List[str]:
    """A cell along KC–KC (for X, whose fit ends at the KC–KC bound; outcome (f)): the evaluations at the bound; the
    fit half's local slope on the KC–KC coordinate near the best, over three windows, with the gain a doubling of the
    factor would buy beside the 0.010 rule; the evaluations binned by the KC–KC factor, with the best fit-half loss and
    the claws' error from the start against each half (the rescore, on the searched wirings); and the cell's class
    with the bound rule left out."""
    x = d['cells'][cell]
    h, r = x['hist'], x['rescore']
    if not len(h):
        return []
    lo, hi = F0.BOUNDS['kckc']
    at = h[np.isclose(h.z_kckc, hi, atol=1e-9, rtol=0)]
    f = x['fit']
    L = [f"Cell {cell} along KC–KC (a description; the bound rule gives {cell} no class). The fit's loss is on its "
         'searched wirings against its fit half; the rescore reads the same evaluations against its read half.', '',
         f"- Evaluations at the KC–KC bound (×{np.exp(hi):.1f}): {len(at)} of {len(h)} "
         f"({', '.join(str(int(e)) for e in at['eval'])}); the best evaluation ({f.get('best_eval')}) is "
         f"{'' if f.get('best_eval') in set(at['eval']) else 'not '}one of them.", '']
    V.add(cell, 'fit', 'n_at_kckc_bound', len(at), half=E.get_cell(cell).fit_half)
    L += ['The fit half\'s loss near the best, regressed on the four coordinates: its change per e-fold of KC–KC, and '
          'what a doubling of the factor would buy (the slope × ln 2; a factor is determined at 0.010).', '',
          '| evaluations within … of the best | evaluations | change per e-fold of KC–KC | a doubling buys |',
          '|---:|---:|---:|---:|']
    for w in SLOPE_WINDOWS:
        s = local_slope(h, w)
        L.append(f"| {w:.3f} | {s['n']} | {_f(s['slope'], '+.4f')} ± {_f(s['se'], '.4f')} | "
                 f"{_f(-s['slope'] * np.log(2), '.4f')} |")
        V.add(cell, 'fit', f'kckc_slope_w{w:g}', s['slope'], half=E.get_cell(cell).fit_half)
        V.add(cell, 'fit', f'kckc_slope_se_w{w:g}', s['se'], half=E.get_cell(cell).fit_half)
    rs = r[r.get('measure', 'strength') == 'strength'] if len(r) and 'measure' in r else r
    if len(rs):
        p = rs.pivot_table(index='eval', columns='half', values='kc_claw_error')
        start = p.loc[0] if 0 in p.index else None
        k = np.exp(h.set_index('eval').z_kckc)
        L += ['', "The evaluations by the KC–KC factor: the best fit-half loss, and the claws' strength error from the "
              "start's (evaluation 0) against each half, mean (and the lowest against the read half).", '',
              '| KC–KC factor | evaluations | best loss, fit half | claws, fit half: mean change | claws, read half: mean '
              'change (lowest) |', '|---|---:|---:|---:|---:|']
        bins = k.map(lambda v: kckc_bin(v, float(np.exp(hi))))
        for _, _, lab in KCKC_BINS:
            ev = [e for e in bins.index[bins == lab] if e in p.index]
            if not ev or start is None:
                continue
            cf, ch = p.loc[ev, 'fit'] - start['fit'], p.loc[ev, 'held'] - start['held']
            L.append(f"| {lab} | {len(ev)} | {_f(h.set_index('eval').loc[ev, 'loss_sim'].min())} | {_f(cf.mean(), '+.3f')} | "
                     f"{_f(ch.mean(), '+.3f')} ({_f(ch.min(), '+.3f')}) |")
            V.add(cell, 'kc_claw', f'change_mean_{lab.split(" (")[0]}', ch.mean(), half='read')
    wb = class_without_bound((recs(d).get(cell) or {}))
    L += ['', f"- {cell}'s class with the bound rule left out: {wb or '–'}."]
    return L + ['']


# ---------------------------------------------------------------------------
# the score

def stage_table(d: dict, V: W.Values) -> List[str]:
    """Strength by odor per megamat stage where each cell is read (its four read wirings, its read half): the flies'
    rows, then per cell the unrefitted family and the cell's point (mean over the wirings, with their range)."""
    m = d['mega']
    fl = W.flies_row(d)
    L = ['Strength by odor in one-fly units (lower is better), megamat panel, where each cell is read: on its four read '
         'wirings against its read half, the mean over the wirings (range). The flies\' rows are E134\'s on all flies.',
         '', '| model or yardstick | ' + ' | '.join(STAGE_NAMES[s] for s in STAGES) + ' | megamat stage error |',
         '|---|' + '---:|' * (len(STAGES) + 1)]
    L.append("| the flies' floor | " + ' | '.join(_f(fl[s].get('floor'), '.2f') for s in STAGES) + ' | – |')
    L.append('| one fly against the others, mean (largest) | ' +
             ' | '.join(f"{_f(fl[s].get('loo_mean'), '.2f')} ({_f(fl[s].get('loo_max'), '.2f')})" for s in STAGES)
             + ' | – |')
    for s in STAGES:
        V.add('flies', s, 'floor', fl[s].get('floor'), half='all')
    for c in CELLS:
        cell = E.get_cell(c)
        for name, lab, key in (('family', 'the family', f'family@{c}'), (variant(c), "the cell's point", c)):
            x = m[(m.cell == c) & (m.name == name) & (m.half == cell.read_half) & m.wiring.isin(cell.read_wirings())]
            if not len(x):
                continue
            L.append(f'| {CELL_NAMES[c]}: {lab} | ' + ' | '.join(_mr(x[s]) for s in STAGES) + f' | {_mr(x.stage)} |')
            for s in ('stage',) + tuple(STAGES):
                V.add(key, s, 'mean', x[s].mean(), half='read')
    return L + ['']


def changes_table(d: dict, V: W.Values) -> List[str]:
    """Each cell's change of each megamat stage's strength error and of the megamat stage error from the family where
    it is read (the score), with the spreads (the family's SD over its eight wirings), and the claws' and somata's
    error over all views (rule 3's guard)."""
    rc = recs(d)
    any_rec = next((r['cls'] for r in rc.values() if (r.get('cls') or {}).get('stages')), None)
    if any_rec is None:
        return []
    L = ["Each cell's change from the family where it is read (mean over its four read wirings; read half; lower is "
         "better), per megamat stage on strength by odor, and over all three views at the claws and the somata (rule "
         "3's guard). The spreads are the family's SD over its eight wirings.", '',
         '| cell | megamat stage error | ' + ' | '.join(STAGE_NAMES[s] for s in STAGES) +
         ' | KC claws, all views | KC somata, all views |', '|---|' + '---:|' * (len(STAGES) + 3)]
    sp = {s: (any_rec['stages'].get(s) or {}).get('spread') for s in STAGES}
    gsp = {s: ((any_rec.get('guard') or {}).get(s) or {}).get('spread') for s in ('kc_claw', 'kc_soma')}
    L.append(f"| the spread | {_f((any_rec.get('megamat') or {}).get('spread'))} | " +
             ' | '.join(_f(sp[s]) for s in STAGES) + f" | {_f(gsp['kc_claw'])} | {_f(gsp['kc_soma'])} |")
    for c in CELLS:
        r = (rc.get(c) or {}).get('cls') or {}
        if not r.get('stages'):
            continue
        g = r.get('guard') or {}
        L.append(f"| {CELL_NAMES[c]} | {_ch(r.get('megamat'))} | " + ' | '.join(_ch(r['stages'].get(s)) for s in STAGES)
                 + f" | {_ch(g.get('kc_claw'))} | {_ch(g.get('kc_soma'))} |")
        V.add(c, 'stage', 'change', (r.get('megamat') or {}).get('change'), half='read')
        for s in STAGES:
            V.add(c, s, 'change', (r['stages'].get(s) or {}).get('change'), half='read')
        for s in ('kc_claw', 'kc_soma'):
            V.add(c, s, 'change_all_views', (g.get(s) or {}).get('change'), half='read')
    return L + ['']


def classes_table(d: dict, V: W.Values) -> List[str]:
    """Each cell's class (E136's classes with E139's price), with what it rests on and E143's share of resamples that
    give it."""
    rc = recs(d)
    sh = d['shares'].set_index('cell') if len(d['shares']) else pd.DataFrame()
    cf = (d['classes'] or {}).get('control_fall')
    L = ["Each cell's class against the family where it is read. The bar is twice the standard error of the megamat "
         f"change over the four read wirings; a fall must also exceed the control fit's {_f(cf)}. \"Other views\" is rule "
         "3's guard (the error over all views). The price is E139's: the calibrated pattern part of the APL effect "
         "above the family's by more than 0.003. The last column is E143's: the share of 1,000 resamples of the flies "
         "that give the scored class.", '',
         '| cell | class | megamat change (bar) | stages worse than their spread | other views worse | held-out panels '
         'worse | price on the APL effect | resamples giving the class |', '|---|---|---:|---|---|---|---:|---:|']
    for c in CELLS:
        rec = rc.get(c)
        if not rec:
            continue
        r = rec.get('cls') or {}
        m = r.get('megamat') or {}
        bar = 2.0 * float(m['se']) if m.get('se') is not None else None
        cls = r.get('cls')
        if cls == 'no class' and r.get('reasons'):
            cls = f"no class ({'; '.join(r['reasons'])})"
        pr = rec.get('price') or {}
        share = sh.loc[c, 'share_scored'] if c in getattr(sh, 'index', []) else None
        price = (f"{'yes' if pr.get('price') else 'no'} ({_f(pr.get('price_value'), '+.4f')})"
                 if pr.get('complete') else '–')
        L.append(f'| {CELL_NAMES[c]} | {cls} | {_ch(m)} ({_f(bar)}) | {_names(r.get("stage_worse") or [])} | '
                 f'{_names(r.get("guard_worse") or [])} | {_names(r.get("panel_worse") or [])} | {price} | '
                 f'{_f(share)} |')
        V.add(c, 'stage', 'bar', bar, half='read')
        V.add(c, 'stage', 'resampled_share_scored', share, half='read')
        V.add(c, 'apl_effect', 'price_value', pr.get('price_value'), panel='apl16', condition='on_off', half='all')
    return L + ['']


def panels_table(d: dict, V: W.Values) -> List[str]:
    """Each cell's change of each held-out panel's KC somata strength error from the family (mean over the cell's four
    read wirings, all flies, as the class compares them), with the spreads; the kiwi mixtures beside."""
    rc = recs(d)
    keys = PANEL_KEYS + BESIDE_KEYS
    any_rec = next((r['cls'] for r in rc.values() if (r.get('cls') or {}).get('panels')), None)
    if any_rec is None:
        return []
    get = lambda r, k: (r.get('panels') or {}).get(k) or (r.get('panels_beside') or {}).get(k)
    L = ["Each cell's change of the held-out panels' KC somata strength error from the family (mean over the cell's "
         'four read wirings, all flies, as the class compares them; lower is better), with the family\'s spread over '
         'its eight wirings. The kiwi mixtures are beside and gate nothing.', '',
         '| cell | ' + ' | '.join(PANEL_NAMES.get(k, k) + (' (beside)' if k in BESIDE_KEYS else '') for k in keys) + ' |',
         '|---|' + '---:|' * len(keys),
         '| the spread | ' + ' | '.join(_f((get(any_rec, k) or {}).get('spread')) for k in keys) + ' |']
    for c in CELLS:
        r = (rc.get(c) or {}).get('cls') or {}
        if not r.get('panels'):
            continue
        cells = []
        for k in keys:
            x = get(r, k) or {}
            cells.append(_ch(x) + (' (worse)' if x.get('complete') and x.get('worse_than_spread') else ''))
            inp, pan, cond = k.split('/')
            V.add(c, 'kc_soma', 'change', x.get('change') if x.get('complete') else None, panel=pan, condition=cond,
                  half='all')
        L.append(f'| {CELL_NAMES[c]} | ' + ' | '.join(cells) + ' |')
    return L + ['']


def experiment_table(d: dict, V: W.Values) -> List[str]:
    """The entry's classes over the three cells and the outcome letters, with E143's share of resamples that agree
    with each."""
    ex = (d['classes'] or {}).get('experiment') or {}
    der = (d['e143'] or {}).get('derived') or {}
    sh = lambda k: _f((der.get(k) or {}).get('share_scored') if isinstance(der.get(k), dict) else der.get(k))
    fx = ex.get('megamat_fixes') or {}
    rows = [('the fall repeats (all three cells beyond the bar and the control fit\'s fall)', ex.get('fall_repeats'),
             'fall_repeats'),
            ('megamat fixes pn_inhadd', fx.get('inhadd'), None), ('megamat fixes APL→bouton', fx.get('aplpn'), None),
            ('megamat fixes KC–KC', fx.get('kckc'), None),
            ('the hold repeats (X holds, and F or W)', ex.get('hold_repeats'), 'hold_repeats'),
            ("the flies' (F differs from A1, W does not)", ex.get('the_flies'), 'the_flies'),
            ("the wirings' (W differs from A1, F does not)", ex.get('the_wirings'), 'the_wirings')]
    L = ["The entry's classes over the three cells, as scored, and E143's share of the 1,000 resamples that agree.", '',
         '| class | scored | resamples agreeing |', '|---|---|---:|']
    for name, v, key in rows:
        L.append(f"| {name} | {W._yn(v) if v is not None else '–'} | {sh(key) if key else '–'} |")
    L.append(f"| cells with no class | {', '.join(ex.get('no_class_cells') or []) or 'none'} | – |")
    L.append(f"| cells that hold | {ex.get('n_cells_hold')} | – |")
    L += ['', '| outcome | met | resamples agreeing |', '|---|---|---:|']
    met = set(ex.get('outcome') or [])
    for k, name in OUTCOMES.items():
        L.append(f"| ({k}) {name} | {W._yn(k in met) if ex.get('complete') else '–'} | "
                 f"{sh(f'outcome_{k}') if f'outcome_{k}' in der else '–'} |")
        V.add('experiment', 'outcome', k, float(k in met) if ex.get('complete') else None, half='read')
    for name, v, key in rows:
        if key:
            V.add('experiment', 'class', key, float(bool(v)) if v is not None else None, half='read')
    return L + ['']


def own_wirings_table(d: dict, V: W.Values) -> List[str]:
    """Each fitted point on its own searched wirings against the fit half and the read half (score_own_wirings.csv;
    the megamat stage error of the wirings' mean views, and the mean over the wirings)."""
    o = d['own']
    if not len(o):
        return []
    L = ['Each fitted point on its own searched wirings (a description): the megamat stage error on strength against '
         "the half it was fitted to and against the other, per stage the mean over the wirings.", '',
         '| cell | half | megamat stage error | mean over the wirings | ' + ' | '.join(STAGE_NAMES[s] for s in STAGES)
         + ' |', '|---|---|---:|---:|' + '---:|' * len(STAGES)]
    for c in CELLS:
        cell = E.get_cell(c)
        for h, lab in ((cell.fit_half, 'fitted'), (cell.read_half, 'read')):
            x = o[(o.cell == c) & (o.half == h)]
            if not len(x):
                continue
            r = x.iloc[0]
            L.append(f'| {CELL_NAMES[c]} | {h} ({lab}) | {_f(r.stage)} | {_f(r.mean_over_wirings)} | '
                     + ' | '.join(_f(r[s]) for s in STAGES) + ' |')
            V.add(c, 'stage', 'own_wirings', r.stage, half=lab)
    return L + ['']


def shortcut_table(d: dict, V: W.Values) -> List[str]:
    """E143's read-once check per cell: the largest difference between the score's reading and the shortcut's (the
    model read once against the resampled flies), on the megamat stage error and on any statistic, and how many of the
    checked values are beyond E143's tolerance. Beyond it, the cell's resampled shares are approximate (as E144 read
    them); the resampled reading gates nothing."""
    s = d['shortcut']
    if not len(s):
        return []
    tol = ((d['e143'] or {}).get('shortcut_check') or {}).get('tol')
    L = [f"E143's read-once check on each cell's read half: the largest difference between the score's reading and the "
         f"shortcut's, on the megamat stage error and on any statistic, and the values beyond E143's tolerance of "
         f"{_f(tol, '.2f')}. Where they are, the cell's resampled shares are approximate.", '',
         '| cell | megamat stage error | any statistic | values beyond the tolerance |', '|---|---:|---:|---:|']
    for c in CELLS:
        x = s[s.cell == c]
        if not len(x):
            continue
        st = x[x.statistic == 'stage']
        n_out = int((~x.within_tol.astype(bool)).sum())
        L.append(f"| {CELL_NAMES[c]} | {_f(st['diff'].abs().max(), '.4f')} | {_f(x['diff'].abs().max(), '.4f')} | "
                 f"{n_out} of {len(x)} |")
        V.add(c, 'stage', 'shortcut_max_diff', st['diff'].abs().max(), half='read')
        V.add(c, 'all', 'shortcut_n_beyond_tol', n_out, half='read')
    return L + ['']


# ---------------------------------------------------------------------------
# realism, the sites and the APL effect

def realism_table(d: dict, V: W.Values) -> List[str]:
    """The four bands on the APL16 input per cell's point and for the family (range over the eight wirings), with the
    KC data's values and the flagged wirings."""
    a = d['apl']
    if not len(a):
        return []
    t = W.data_targets()
    L = ["The APL16 input (held out), range over the eight wirings. Bands: response rate with APL on 0.045 to 0.085; "
         "APL-off / APL-on ratio 1.6 to 2.9; response rate ρ with APL on at least 0.37; odors per KC within 0.155 (APL "
         "on) and 0.223 (APL off) of the flies (total variation distance).", '',
         '| model | response rate, APL on | APL off | ratio | response rate ρ, APL on | APL off | odors per KC, APL on | '
         'APL off | largest spike count, APL on / off | flagged |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    if t:
        g = lambda k, fmt='.3f': _f(t.get(k), fmt)
        L.append(f"| the KC data (control flies / APL output blocked (TNT), noise criterion) | {g('noise_ctrl_frac')} | "
                 f"{g('noise_tnt_frac')} | {g('noise_rate_ratio_fraction', '.2f')} | – | – | "
                 f"band {g('noise_ctrl_odors_per_kc_tv_max')} | band {g('noise_tnt_odors_per_kc_tv_max')} | – | – |")
    fam = a[(a.arm == 'family') & (a.cell == a.cell.iloc[0])]
    for key, lab, x in [('family', 'the family (unrefitted)', fam)] + \
            [(c, CELL_NAMES[c], a[a.arm == variant(c)]) for c in CELLS]:
        if not len(x):
            continue
        flag = x.realism_flag.fillna('').astype(str) if 'realism_flag' in x else pd.Series([''] * len(x))
        named = '; '.join(f'{w}: {fl}' for w, fl in zip(x.wiring, flag) if fl) or 'none'
        L.append(f'| {lab} | {_rng(x.rate_on)} | {_rng(x.rate_off)} | {_rng(x.rate_ratio, ".2f")} | '
                 f'{_rng(x.rate_rho_on, ".2f")} | {_rng(x.rate_rho_off, ".2f")} | {_rng(x.opk_on_tv)} | '
                 f'{_rng(x.opk_off_tv)} | {_rng(x.max_count_on, ".0f")} / {_rng(x.max_count_off, ".0f")} | {named} |')
        for col in ('rate_on', 'rate_off', 'rate_ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv',
                    'max_count_on', 'max_count_off'):
            V.add(key, 'kc_soma', f'{col}_min', x[col].min(), panel='apl16', condition='realism', half='all')
            V.add(key, 'kc_soma', f'{col}_max', x[col].max(), panel='apl16', condition='realism', half='all')
        V.add(key, 'kc_soma', 'n_flagged', int((flag != '').sum()), panel='apl16', condition='realism', half='all')
    return L + ['']


def scales_table(d: dict, V: W.Values, panel: str = 'megamat') -> List[str]:
    """A panel's operating point and the absolute scales each cell's holds reached (range over its reads), with the
    largest spike counts."""
    name = dict(apl16='APL16').get(panel, panel)
    L = [f"The {name} reads' operating point and the absolute scales the holds reached, range over each cell's reads "
         + ("(megamat on its four read wirings; the flies' response rate with APL on 0.0915)." if panel == 'megamat'
            else '(the eight wirings).'), '',
         '| cell | reads | response rate, APL on | APL off | APL→bouton | APL→KC | KC–KC | somatic threshold | '
         'pn_inhadd | pn_inhsc | largest spike count, APL on / off |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for c in CELLS:
        r = d['cells'][c]['reads']
        if not len(r):
            continue
        x = r[r.panel == panel]
        if not len(x):
            continue
        col = lambda k, fmt: _rng(x[k], fmt) if k in x else '–'
        L.append(f'| {CELL_NAMES[c]} | {len(x)} | {col("rate_on", ".4f")} | {col("rate_off", ".4f")} | '
                 f'{col("wAPLPN_scale", ".4f")} | {col("wAPLKC_scale", ".3f")} | {col("wKCKC_scale", ".1f")} | '
                 f'{col("fixed_thr", ".1f")} | {col("pn_inhadd", ".0f")} | {col("pn_inhsc", ".1f")} | '
                 f'{col("max_count_on", ".0f")} / {col("max_count_off", ".0f")} |')
        for k in ('rate_on', 'rate_off', 'wAPLPN_scale', 'wAPLKC_scale', 'wKCKC_scale', 'fixed_thr', 'pn_inhadd',
                  'max_count_on', 'max_count_off'):
            if k in x:
                V.add(c, 'reads', f'{k}_min', x[k].min(), panel=panel, condition='scales', half='all')
                V.add(c, 'reads', f'{k}_max', x[k].max(), panel=panel, condition='scales', half='all')
    return L + ['']


def sites_table(d: dict, V: W.Values) -> List[str]:
    """Where APL's inhibition is made (score_sites.csv): the share of the claws' summed drive APL removes and the share
    of that loss made at the PN boutons, mean over the wirings (range), beside the connectome's share (E94)."""
    from experiments import E94_bouton_share as E94
    s = d['sites']
    if not len(s):
        return []
    L = ["Where APL's inhibition is made, window means, mean over the wirings (range): megamat on the cell's four read "
         "wirings, APL16 on the eight. \"Removed\" is the share of the claws' summed drive APL removes; \"at the "
         f"boutons\" the share of that loss made at the PN boutons, against the connectome's {E94.ANATOMICAL_SHARE:g} "
         "(APL's calyx synapses onto boutons over those onto boutons and KCs).", '',
         '| model | megamat: removed | at the boutons | APL16: removed | at the boutons |', '|---|---:|---:|---:|---:|']
    for c in CELLS:
        x = s[s.variant == variant(c)]
        if not len(x):
            continue
        mm, aa = x[x.panel == 'megamat'], x[x.panel == 'apl16']
        L.append(f'| {CELL_NAMES[c]} | {_mr(mm.drive_removed)} | {_mr(mm.bouton_share)} | {_mr(aa.drive_removed)} | '
                 f'{_mr(aa.bouton_share)} |')
        for pan, y in (('megamat', mm), ('apl16', aa)):
            if len(y):
                V.add(c, 'calyx', 'bouton_share_mean', y.bouton_share.mean(), panel=pan, half='all')
                V.add(c, 'calyx', 'drive_removed_mean', y.drive_removed.mean(), panel=pan, half='all')
    V.add('connectome', 'calyx', 'bouton_share', E94.ANATOMICAL_SHARE, half='all')
    return L + ['']


def apl_effect_table(d: dict, V: W.Values) -> List[str]:
    """The APL effect, held out (score_apl16.csv, the calibrated readout against all recordings), mean over the eight
    wirings: the loss, the mean APL effect, the pattern part and ρ on the per-pair APL effect, for the family and each
    cell's point, with the KC data's floors."""
    a = d['apl']
    if not len(a):
        return []
    t = W.data_targets()
    cols = (('cal_loss', 'loss', '.4f'), ('cal_mean', 'mean APL effect', '+.3f'), ('cal_pattern', 'pattern part', '.4f'),
            ('cal_rho', 'ρ on the per-pair APL effect', '.2f'))
    L = ['The APL effect on odor-odor distances (a held-out test; gates nothing), under the calibrated readout against '
         'all recordings, mean over the eight wirings. The KC data\'s row gives the flies\' mean APL effect, the floors '
         'of the loss and of the pattern part (a resampled set of recordings against the measured one) and the ceiling '
         'of ρ.', '', '| model | ' + ' | '.join(n for _, n, _ in cols) + ' |', '|---|' + '---:|' * len(cols)]
    if t:
        L.append(f"| the KC data (floors; ceiling) | {_f(t.get('apl_effect_rmse_floor'), '.4f')} | "
                 f"{_f(t.get('apl_effect_mean'), '+.3f')} | {_f(t.get('apl_effect_pattern_floor'), '.4f')} | "
                 f"{_f(t.get('apl_effect_rho_ceiling'), '.2f')} |")
    fam = a[(a.arm == 'family') & (a.cell == a.cell.iloc[0])]
    for key, lab, x in [('family', 'the family (unrefitted)', fam)] + \
            [(c, CELL_NAMES[c], a[a.arm == variant(c)]) for c in CELLS]:
        if not len(x):
            continue
        L.append(f'| {lab} | ' + ' | '.join(_f(x[k].mean(), fmt) for k, _, fmt in cols) + ' |')
        for k, _, _ in cols:
            V.add(key, 'apl_effect', k, x[k].mean(), panel='apl16', condition='on_off', half='all')
    return L + ['']


def without_aa_table(d: dict, V: W.Values) -> List[str]:
    """The megamat somata's strength error with and without aa's entry where each cell is read (mean over its read
    wirings), for the family and the cell's point."""
    w = d['without']
    if not len(w):
        return []
    L = ["The megamat KC somata's strength error where each cell is read (its read wirings, its read half), with and "
         "without aa's entry, mean over the wirings, and aa's share of the error.", '',
         '| cell | model | with aa | without aa | aa\'s share |', '|---|---|---:|---:|---:|']
    for c in CELLS:
        cell = E.get_cell(c)
        for name, lab, key in (('family', 'the family', f'family@{c}'), (variant(c), "the cell's point", c)):
            x = w[(w.cell == c) & (w.arm == name) & w['row'].isin(cell.read_wirings())]
            if not len(x):
                continue
            L.append(f'| {CELL_NAMES[c]} | {lab} | {_f(x.view_error.mean())} | {_f(x.view_error_without.mean())} | '
                     f'{_f(x.view_share.mean())} |')
            V.add(key, 'kc_soma', 'with_aa', x.view_error.mean(), half='read')
            V.add(key, 'kc_soma', 'without_aa', x.view_error_without.mean(), half='read')
    return L + ['']


def e137_table(d: dict, V: W.Values) -> List[str]:
    """The mixtures under E137's readings of the mixture flies (strength by odor; mean over the eight wirings), for the
    family and each cell's point."""
    e = d['e137']
    if not len(e):
        return []
    e = e[e['views'].astype(str) != 'all'] if 'views' in e else e
    combos = sorted({(p, r) for p, r in zip(e.panel, e.reading)})
    nm = {'kmix': 'kiwi mixtures', 'cmix': 'control mixtures'}
    L = ["The mixtures' KC somata under E137's readings of the mixture flies (each reading's own units and trial noise), "
         'strength by odor, mean over the eight wirings. Beside; gates nothing.', '',
         '| model | ' + ' | '.join(f'{nm.get(p, p)}, {W.READING_NAMES.get(r, r)}' for p, r in combos) + ' |',
         '|---|' + '---:|' * len(combos)]
    fam = e[(e.arm == 'family') & (e.cell == e.cell.iloc[0])]
    for key, lab, x in [('family', 'the family (unrefitted)', fam)] + \
            [(c, CELL_NAMES[c], e[e.arm == variant(c)]) for c in CELLS]:
        if not len(x):
            continue
        cells = []
        for p, r in combos:
            y = x[(x.panel == p) & (x.reading == r)]
            cells.append(_f(y.error.mean()) if len(y) else '–')
            if len(y):
                V.add(key, 'kc_soma', f'e137_{r}', y.error.mean(), panel=p, half='all')
        L.append(f'| {lab} | ' + ' | '.join(cells) + ' |')
    return L + ['']


# ---------------------------------------------------------------------------
# the prediction and the reads

def prediction_table(d: dict, V: W.Values) -> List[str]:
    """The planner's chances (prediction.json, written before any fit) beside what was scored."""
    p = d['prediction'] or {}
    ex = (d['classes'] or {}).get('experiment') or {}
    cls = ex.get('cls') or {}
    rc = recs(d)
    kk = lambda c: (((rc.get(c) or {}).get('repeats') or {}).get('kckc') or {}).get('repeats')
    scored = {'fall_repeats_all_three': W._yn(ex.get('fall_repeats')) if 'fall_repeats' in ex else '–',
              'holds_in_X': W._yn(cls.get('X') == 'holds') if 'X' in cls else '–',
              'holds_in_all_three': W._yn(all(cls.get(c) == 'holds' for c in E.FIT_CELLS)) if cls else '–',
              'kckc_determined_at_1.8_or_more_in_all_three': W._yn(all(kk(c) for c in E.FIT_CELLS)) if rc else '–',
              'kckc_determined_at_1.8_or_more_in_none': W._yn(not any(kk(c) for c in E.FIT_CELLS)) if rc else '–'}
    names = {'fall_repeats_all_three': 'the megamat fall repeats in all three cells',
             'holds_in_X': '"holds" in cell X', 'holds_in_all_three': '"holds" in all three cells',
             'kckc_determined_at_1.8_or_more_in_all_three': 'KC–KC determined at ×1.8 or more in all three',
             'kckc_determined_at_1.8_or_more_in_none': 'KC–KC determined at ×1.8 or more in none'}
    L = ["The planner's chances, written before any fit (prediction.json), and what was scored.", '',
         '| event | the planner\'s chance | scored |', '|---|---|---|']
    for k, n in names.items():
        L.append(f'| {n} | {p.get(k, "–")} | {scored[k]} |')
    return L + ['']


def lanes_table(d: dict, V: W.Values) -> List[str]:
    """The reads of the three fitted cells by host and lane: count, commits, dirty."""
    rows = []
    for c in E.FIT_CELLS:
        r = d['cells'][c]['reads']
        if len(r):
            rows.append(r.assign(cell=c))
    if not rows:
        return []
    r = pd.concat(rows, ignore_index=True)
    L = ["The fitted cells' reads by lane (<cell>/reads_<host>_<lane>.csv).", '',
         '| cell | host | lane | reads | queue | pyolfsysm | any checkout dirty | all held |',
         '|---|---|---|---:|---|---|---|---|']
    for (c, h, ln), x in r.groupby(['cell', 'host', 'lane'], sort=True):
        dirty = any(str(v) == 'True' for col in ('queue_dirty', 'pyolfsysm_dirty', 'al_analysis_dirty', 'hong2p_dirty')
                    if col in x for v in x[col])
        held = None
        if {'frac_on_target', 'rate_on'} <= set(x):
            held = bool((np.abs(x.rate_on / x.frac_on_target - 1) <= F0.HOLD_TOL).all())
        L.append(f"| {c} | {h} | {ln} | {len(x)} | {', '.join(sorted({str(v)[:7] for v in x.queue_commit}))} | "
                 f"{', '.join(sorted({str(v)[:7] for v in x.pyolfsysm_commit}))} | {W._yn(dirty)} | "
                 f"{W._yn(held) if held is not None else '–'} |")
        V.add(c, 'reads', 'n', len(x), condition=f'{h}/{ln}', half='all')
    return L + ['']


def main(score_dir: Optional[Path] = None) -> Path:
    score_dir = Path(score_dir) if score_dir is not None else E.out_dir()
    d = load(score_dir)
    V = W.Values()
    prov, prov_r = d['prov'], d['prov_r']
    head = ['# E145 write-up tables', '',
            f"From results/E145_stage_fit_replication/ (the score at queue {str(prov.get('queue_commit'))[:7]}, "
            f"{prov.get('time')}; E143's resampled reading at queue {str(prov_r.get('queue_commit'))[:7]}, "
            f"{prov_r.get('time')}); written by analyses/E145_writeup_tables.py.", '']
    parts = [('The cells', cells_table), ('The fits', fits_table), ('The returns and the classed points', returns_table),
             ('The rescore', rescore_table), ('X along KC–KC', kckc_table), ('The stage table', stage_table), ('The megamat stages', changes_table),
             ('The classes', classes_table), ('The held-out panels', panels_table),
             ("The experiment's classes and outcomes", experiment_table),
             ('Each point on its own searched wirings', own_wirings_table),
             ("E143's read-once check", shortcut_table), ('Realism: the bands', realism_table),
             ('Realism: the operating point and the scales', scales_table),
             ('Realism: the scales on the APL16 input', lambda d_, V_: scales_table(d_, V_, 'apl16')),
             ("Where APL's inhibition is made", sites_table), ('The APL effect, held out', apl_effect_table),
             ('The megamat somata without aa', without_aa_table), ("The mixtures under E137's readings", e137_table),
             ('The prediction', prediction_table), ('The reads', lanes_table)]
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
