"""E145's score: each cell's fitted point classed by E136's classes with E139's price against the unrefitted family
where the cell is read, and the experiment's classes and outcome (experiment_priorities.md, E145; the lock's CLASSES,
OUTCOMES and REPEAT).

Each cell is scored inside E136's `cell` (E145's `e136_cell`), so that E136's own functions read it: its fitted point is
registered with E136 as a variant ('E145_<cell>', `register_cell`), whose runs are the fit's best evaluation on the
cell's searched wirings (as E136 reads a fitted arm on S) and E145's reads everywhere else; E136's `score_megamat`,
`score_panels`, `apl16_rows` and `classify` then read the cell's read wirings (E136's T_WIRINGS inside the cell) and
its read half (E136's 'held' inside the cell). The spreads are the family's over its eight wirings as E136 scored them
(E136's classes.json: 'spreads', and 'beside' 'spreads_all_views' for rule 3's guard); the control's fall is E136's
(the lock's yardstick). E139's price reads the APL16 reads' calibrated pattern part against the family on the cell's
read wirings (deciding) and on all eight (beside). E136's own cell (A1) is scored the same way from E136's runs, and
must reproduce E136's class and megamat change (`check_e136`), which checks the machinery.

Beside, classing nothing: each fitted point on its own searched wirings against both halves; the megamat somata with
and without aa's entry; the mixtures under E137's readings and the kiwi mixtures; where APL's inhibition is made. E143's
resampled reading per cell is `part_resampled` (written once the cells' reads are in)."""
from __future__ import annotations

import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E145_stage_fit_replication as E

# the tolerance of the check against E136's own score
CHECK_TOL = 1e-9


def variant(cell: str) -> str:
    """A cell's fitted point as an E136 variant name, 'E145_<cell>'."""
    E.get_cell(cell)
    return f'E145_{cell}'


def _e136_dirs() -> Dict[str, List[Path]]:
    """E136's fit and read folders on the adopted measure (its A1 cell's runs), read outside any cell."""
    with E.e136() as F:
        return dict(fit=F.fit_dirs(E.ARM), read=F.read_dirs())


def cell_fit(cell: str) -> dict:
    """A cell's fit record: E136's fit_A1.json for E136's own cell, E145's for the others ({} before its fit ends)."""
    if cell == 'A1':
        with E.e136() as F:
            return E._read_json(F.out_dir() / f'fit_{E.ARM}.json')
    return E.fit_record(cell)


def cell_loader(cell: str):
    """The loader of a cell's variant for E136's `register_variant`: on the megamat input on the cell's searched
    wirings the fit's best evaluation, every other point the cell's read (E145's tags; E136's own for its cell);
    None if absent. The folders and the fit record are read when the loader is made, which must be outside any E136
    cell (inside one, E136's folders are E145's)."""
    from kcq import stage_readings as R
    from experiments import E136_stage_population_fit as F0
    if F0.PATHS:
        raise RuntimeError('cell_loader is made outside an E136 cell')
    c = E.get_cell(cell)
    fit = cell_fit(cell)
    if cell == 'A1':
        d = _e136_dirs()
        fit_dirs, read_dirs = d['fit'], d['read']
        tag = lambda panel, w: f'{E.ARM}_{panel}_{w}'
    else:
        fit_dirs = [E.local_dir(cell, f'fit_{E.ARM}'), E.archive_dir(cell) / f'fit_{E.ARM}']
        read_dirs = E.read_dirs(cell)
        tag = lambda panel, w: E.point_tag(E.Point(cell, panel, w))
    searched = set(c.searched_wirings())

    def load(pt, claw_after: bool = False):
        try:
            if pt.panel == 'megamat' and pt.wiring in searched:
                if not fit:
                    return None
                return R.load_step3_run(f"{E.ARM}_e{fit['best_eval']}_{pt.wiring}", dirs=fit_dirs)
            return R.load_step3_run(tag(pt.panel, pt.wiring), dirs=read_dirs)
        except FileNotFoundError:
            return None
    return load


def register_cell(F, cell: str, loader=None) -> str:
    """Registers a cell's variant with E136 (`register_variant`; `loader`, default `cell_loader`); returns its name."""
    v = variant(cell)
    F.register_variant(v, loader if loader is not None else cell_loader(cell))
    return v


def e136_classes() -> dict:
    """E136's classes.json on the adopted measure: the spreads, the guard's spreads, the control's fall and A1's
    class."""
    with E.e136() as F:
        return E._read_json(F.out_dir() / 'classes.json')


def repeats(cell: str, fit: Optional[Mapping] = None) -> dict:
    """Whether each factor of REPEAT repeats in a cell, at its classed point (fit record 'classed_z', else 'best_z';
    pn_inhadd by E136's `pn_norm_kws`); the gain on pn_inhsc reported beside, reading nothing. `fit` (default
    `cell_fit`, read outside any E136 cell) is the cell's fit record."""
    from experiments import E136_stage_population_fit as F
    fit = cell_fit(cell) if fit is None else fit
    if not fit:
        return dict(complete=False)
    z = {k: float(v) for k, v in (fit.get('classed_z') or fit['best_z']).items()}
    det = fit.get('determined') or {}
    pn = F.pn_norm_kws(z)
    out = dict(complete='determined' in fit, z=z)
    inhadd = float(pn.get('pn_inhadd', np.nan))
    r = E.REPEAT['inhadd']
    out['inhadd'] = dict(value=inhadd, repeats=bool(r['lo'] <= inhadd <= r['hi']),
                         determined=bool((det.get('inhadd') or {}).get('determined')))
    for k in ('aplpn', 'kckc'):
        r = E.REPEAT[k]
        x = float(np.exp(z[k]))
        d = bool((det.get(k) or {}).get('determined'))
        inside = x >= r['lo'] and (r['hi'] is None or x <= r['hi'])
        out[k] = dict(value=x, determined=d, repeats=bool(d and inside))
    out['inhsc_gain'] = dict(value=float(np.exp(z.get('inhsc', 0.0))),
                             determined=bool((det.get('inhsc') or {}).get('determined')), read=False)
    out['pn_inhsc'] = float(pn.get('pn_inhsc', np.nan))
    return out


def apl_parts(F, names: Sequence[str]) -> Dict[str, Dict[str, dict]]:
    """E139's `apl_effect_parts` of each name's APL16 reads, per wiring of the eight."""
    from experiments import E139_normalization_ladder as E139
    out = {}
    for n in names:
        out[n] = {}
        for w in F.WIRINGS:
            r = F.load_read(F.ReadPoint(n, 'apl16', w, None))
            if r is not None:
                out[n][w] = E139.apl_effect_parts(E139.run_outputs(r))
    return out


def sites(F, v: str) -> pd.DataFrame:
    """Where APL's inhibition is made in a cell's megamat reads (its read wirings) and APL16 reads (all eight wirings):
    analyses/E136_writeup_tables.py's `sites_row` on the wiring's claw x bouton matrix."""
    from analyses import E136_writeup_tables as T
    rows = []
    for panel, ws in (('megamat', F.T_WIRINGS), ('apl16', F.WIRINGS)):
        for w in ws:
            r = F.load_read(F.ReadPoint(v, panel, w, None))
            if r is None:
                continue
            on, off = (r.conditions[c] for c in ('on', 'off'))
            rows.append(dict(variant=v, panel=panel, wiring=w, **T.sites_row(T.wiring_W(w), on, off, False)))
    return pd.DataFrame(rows)


def score_cell(cell: str, beside: bool = True) -> dict:
    """One cell scored inside its E136 cell: the family and the cell's variant on the megamat stages against the read
    half (and against all flies and the fit half beside), on the read wirings; the held-out panels; the APL16 rows and
    their realism flags; rule 3's guard on every view; E136's `classify` with the spreads and the control's fall as E136
    scored them; E139's price; the repeats; with `beside`, the point on its own searched wirings against both halves,
    the somata without aa's entry, E137's readings of the mixtures and the sites. Returns dict(rec, tables)."""
    from experiments import E139_normalization_ladder as E139
    from kcq import stage_readings as R
    lock = E.require_lock()
    c = E.get_cell(cell)
    ec = e136_classes()
    spreads = {k: float(v) for k, v in ec['spreads'].items()}
    g_spreads = {k: float(v) for k, v in (ec.get('beside') or {}).get('spreads_all_views', {}).items()}
    control_fall = float(lock['yardsticks']['control_fall'])
    fit = cell_fit(cell)
    loader = cell_loader(cell)
    with E.e136_cell(cell) as F:
        v = register_cell(F, cell, loader)
        stages = F.readable_stages(F.require_gate()['classes'])
        names = [F.FAMILY, v]
        T, S = F.T_WIRINGS, F.S_WIRINGS
        mrows, mega, g_mega = [], {}, {}
        for n in names:
            for half in ('held', None, 'fit'):
                sc = F.score_megamat(n, stages, half, T)
                if sc is None:
                    continue
                vals = F.megamat_values(sc, stages, T)
                if half == 'held':
                    mega[n] = vals
                hname = {'held': c.read_half, 'fit': c.fit_half, None: 'all'}[half]
                for w, x in vals.items():
                    mrows.append(dict(cell=cell, name=n, half=hname, wiring=w, **x))
                mrows.append(dict(cell=cell, name=n, half=hname, wiring='common', stage=sc['stage_error'],
                                  downstream=sc['stage_error_downstream'], **{s: sc[f'{s}_error'] for s in stages}))
            sc = F.score_megamat(n, stages, 'held', T, views=None)
            if sc is not None:
                g_mega[n] = F.megamat_values(sc, stages, T)
        p_rows = [F.score_panels(n) for n in names]
        panels = pd.concat([p for p in p_rows if len(p)], ignore_index=True) if any(len(p) for p in p_rows) \
            else pd.DataFrame()
        g_rows = [F.score_panels(n, views=None) for n in names]
        g_panels = pd.concat([p for p in g_rows if len(p)], ignore_index=True) if any(len(p) for p in g_rows) \
            else pd.DataFrame()
        pv = {n: F.panel_values(panels, n) for n in names}
        g_pv = {n: F.panel_values(g_panels, n) for n in names}
        a_rows = [F.apl16_rows(n) for n in names]
        apl = pd.concat([a for a in a_rows if len(a)], ignore_index=True) if any(len(a) for a in a_rows) \
            else pd.DataFrame()
        realism = {n: F.realism_of(apl[apl.arm == n]) if len(apl) else {} for n in names}
        guard = dict(mega=g_mega, panels=g_pv, spreads=g_spreads)
        if v in mega and F.FAMILY in mega:
            cls = F.classify(v, stages, mega, pv, spreads, realism, fit, control_fall, guard)
        else:
            cls = dict(cls='reads missing')
        parts = apl_parts(F, names)
        price = E139.price_of(parts, v, F.FAMILY)
        rec = dict(cell=cell, label=c.label, variant=v, searched=list(S), read_wirings=list(T), fit_half=c.fit_half,
                   read_half=c.read_half, stages=list(stages), cls=cls, price=price, repeats=repeats(cell, fit),
                   fit={k: fit.get(k) for k in ('best_eval', 'n_evals', 'factors', 'at_bound', 'at_bound_classed',
                                                 'classed_z', 'classed_loss', 'pn_kws')})
        tables = dict(megamat=pd.DataFrame(mrows), panels=panels, panels_all_views=g_panels, apl16=apl)
        if beside:
            own = []
            for half in ('fit', 'held'):
                sc = F.score_megamat(v, stages, half, S)
                if sc is not None:
                    hname = c.fit_half if half == 'fit' else c.read_half
                    own.append(dict(cell=cell, name=v, wirings='searched', half=hname, stage=sc['stage_error'],
                                    mean_over_wirings=float(np.mean([sc[f'stage_error_{w}'] for w in S])),
                                    **{s: sc[f'{s}_error'] for s in stages}))
            tables['own_wirings'] = pd.DataFrame(own)
            wo = [F.stage_without_odor(n, 'kc_soma', 'aa', 'held', T) for n in names]
            tables['without_aa'] = pd.concat([x.assign(cell=cell) for x in wo if x is not None], ignore_index=True) \
                if any(x is not None for x in wo) else pd.DataFrame()
            e137 = [F.score_e137_beside(n) for n in names]
            e137 = pd.concat([x for x in e137 if len(x)], ignore_index=True) if any(len(x) for x in e137) \
                else pd.DataFrame()
            tables['e137'] = e137.assign(cell=cell) if len(e137) else e137
            rec['e137_beside'] = F.e137_beside(e137, [v]) if len(e137) else {}
            try:
                tables['sites'] = sites(F, v).assign(cell=cell)
            except Exception as e:   # a description only
                rec['sites_error'] = str(e)
    for k in ('megamat', 'panels', 'panels_all_views', 'apl16'):
        if len(tables[k]):
            tables[k].insert(0, 'cell', cell) if 'cell' not in tables[k] else None
    return dict(rec=rec, tables=tables)


def check_e136(rec: Mapping) -> dict:
    """E136's own cell scored here against E136's classes.json: the class and the megamat change, value and family
    value on T against half b, within CHECK_TOL. RuntimeError otherwise."""
    a = (e136_classes().get('classes') or {}).get(E.ARM) or {}
    mine = rec['cls']
    out = dict(cls=(mine.get('cls'), a.get('cls')))
    diffs = {k: abs(float(mine['megamat'][k]) - float(a['megamat'][k])) for k in ('change', 'arm', 'ref', 'se')}
    out['diffs'] = diffs
    out['ok'] = bool(mine.get('cls') == a.get('cls') and max(diffs.values()) <= CHECK_TOL)
    if not out['ok']:
        raise RuntimeError(f"E145's score does not reproduce E136's A1: {out}")
    return out


def fall_repeats_in(rec: Mapping, control_fall: float) -> Optional[bool]:
    """Whether a cell's megamat stage error falls from the family's by more than the bar (twice its standard error)
    and by more than the control's fall; None where its megamat reads are incomplete."""
    m = (rec['cls'] or {}).get('megamat') or {}
    if not m.get('complete'):
        return None
    fall = -float(m['change'])
    return bool(fall > 2.0 * float(m['se']) and fall > control_fall)


def experiment_classes(recs: Mapping[str, Mapping], control_fall: float) -> dict:
    """The entry's classes over the three cells (CLASSES) and the outcome letters whose conditions hold (OUTCOMES), from
    the cells' records (`score_cell`); every class None while a cell is unscored."""
    cells = list(E.FIT_CELLS)
    have = [k for k in cells if k in recs and (recs[k]['cls'] or {}).get('cls') not in (None, 'reads missing')]
    cls = {k: recs[k]['cls'].get('cls') for k in recs}
    fall = {k: fall_repeats_in(recs[k], control_fall) for k in recs}
    rep = {k: recs[k]['repeats'] for k in recs}
    out = dict(cells_scored=have, cls=cls, fall=fall)
    if len(have) < len(cells):
        out.update(complete=False)
        return out
    kk = lambda k: bool(rep[k].get('kckc', {}).get('repeats'))
    out['fall_repeats'] = all(fall[k] for k in cells)
    out['megamat_fixes'] = {f: all(rep[k].get(f, {}).get('repeats') for k in cells) for f in ('inhadd', 'aplpn', 'kckc')}
    out['hold_repeats'] = cls['X'] == 'holds' and (cls['F'] == 'holds' or cls['W'] == 'holds')
    a1 = recs.get('A1')
    if a1 is not None:
        diff = {k: cls[k] != cls['A1'] or kk(k) != kk('A1') for k in ('F', 'W')}
        out['the_flies'] = bool(diff['F'] and not diff['W'])
        out['the_wirings'] = bool(diff['W'] and not diff['F'])
    out['no_class_cells'] = [k for k in cells if cls[k] == 'no class']
    out['n_cells_hold'] = sum(cls[k] == 'holds' for k in cells)
    at_bound = [k for k in cells if (recs[k]['fit'].get('at_bound_classed') if 'at_bound_classed' in recs[k]['fit']
                                     else recs[k]['fit'].get('at_bound'))]
    letters = []
    fx = out['megamat_fixes']
    if out['fall_repeats'] and out['hold_repeats'] and all(fx.values()):
        letters.append('a')
    if out['fall_repeats'] and fx['inhadd'] and fx['aplpn'] and \
            sum((not kk(k)) and cls[k] == 'megamat only' for k in cells) >= 2:
        letters.append('b')
    if out['hold_repeats'] and not all(fx.values()):
        letters.append('c')
    if out.get('the_flies') or out.get('the_wirings') or (not letters and len({cls[k] for k in cells}) > 1):
        letters.append('d')
    if not out['fall_repeats']:
        letters.append('e')
    if at_bound:
        letters.append('f')
    out.update(complete=True, outcome=letters, at_bound_cells=at_bound)
    return out


def part_score(cells: Sequence[str] = ('A1',) + E.FIT_CELLS) -> dict:
    """The score: every cell of `cells` scored (`score_cell`), E136's own checked against E136's (`check_e136`), the
    experiment's classes (`experiment_classes`); writes score_megamat.csv, score_panels.csv, score_panels_all_views.csv,
    score_apl16.csv, score_own_wirings.csv, score_without_aa.csv, score_panels_e137.csv, score_sites.csv, classes.json
    and provenance_score.json in the results folder. A cell whose fit has not ended is left out."""
    t0 = time.time()
    lock = E.require_lock()
    control_fall = float(lock['yardsticks']['control_fall'])
    recs, tabs = {}, {}
    for k in cells:
        if k != 'A1' and not cell_fit(k):
            continue
        s = score_cell(k)
        recs[k] = s['rec']
        for name, t in s['tables'].items():
            tabs.setdefault(name, []).append(t)
    check = check_e136(recs['A1']) if 'A1' in recs else None
    out = E.out_dir()
    for name, ts in tabs.items():
        ts = [t for t in ts if len(t)]
        if ts:
            pd.concat(ts, ignore_index=True).to_csv(out / f'score_{name}.csv', index=False)
    classes = dict(cells=recs, experiment=experiment_classes(recs, control_fall), check_e136=check,
                   spreads_source="E136's classes.json (measure_strength)", control_fall=control_fall)
    E._write_json(out / 'classes.json', classes)
    E._write_json(out / 'provenance_score.json', dict(part='score', script='experiments/E145_score.py',
                  cells=list(recs), lock_written=lock.get('written'), seconds=round(time.time() - t0, 1),
                  host=E._host(), time=E._now(),
                  **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))}))
    return dict(cells={k: r['cls'].get('cls') for k, r in recs.items()}, experiment=classes['experiment'],
                check_e136=check)
