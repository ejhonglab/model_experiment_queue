"""E144's score (design item 3, with item 4's readings beside): the reads of experiments/E144_a1_malecns_calyces.py
scored per calyx, each setting against the calyx's own family on the same wirings, by E136's classes and E139's price
(E141's score, on the calyx's wirings); and beside it E143's resampled reading of its class tables (design item 4's
"E143's intervals over flies on every change", the lock's first amendment; `part_resampled`). Kept out of the read
module, which the read lanes load, so that this file can change while lanes run (E108's `_still_clean` guard checks
only the files a lane has loaded).

Usage (from the repository root): ./venv/bin/python -u -m experiments.E144_a1_malecns_calyces --part score
[--part resampled, after the score]
"""
from __future__ import annotations

import contextlib
import json
import os
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments.E144_a1_malecns_calyces import (CALYCES, DOSES, EXPERIMENT, MALE, PANELS, SCRIPT, SETTINGS,
                                                 WIRING_SEEDS, WIRINGS_OF, Point, _host, _now, _read_json, _write_json,
                                                 local_dir, lock_path, out_dir, point_tag, read_dirs, reads_table,
                                                 require_lock, seed_index, tables, variant)


# the 'first' real wiring's counterpart of each setting on H (E139's and E141's scored reads of the hemibrain's real
# wiring: A1 is E136's fitted A1, which E139 and E141 take as B with the rise)
FIRST_COUNTERPART: Dict[str, str] = {'family': 'family', 'F': 'F_inhadd1000', 'B': 'B_inhadd1000',
                                     'D1.8': 'D1.8_inhadd1000', 'A1': 'A1', 'D5.8': 'D5.8_inhadd1000'}
FIRST_SOURCES: Tuple[str, ...] = ('E141_kckc_rise_forms', 'E139_normalization_ladder')
# the dose ladder's hemibrain values (the entry's: validation2's change along B, D1.8, A1 and D5.8 on the real wiring)
HEMIBRAIN_LADDER: Dict[str, float] = {'B': 0.080, 'D1.8': 0.021, 'A1': -0.016, 'D5.8': -0.036}
FAMILY_INSIDE_STAGES: Tuple[str, ...] = ('pn_bouton', 'kc_claw', 'kc_soma')
RISES: Tuple[str, ...] = ('D1.8', 'A1', 'D5.8')
WITHOUT_ODOR = 'aa'


@contextlib.contextmanager
def e136_for(calyx: str):
    """E136 on the adopted measure (E139's `e136`) with its `seed_index` giving this experiment's index of the calyx's
    wirings (the lock's: real 0, P1 1, P2 2, P3 3), so that E136's StageLoss reads each run with its own seeds; restored
    after."""
    from experiments import E139_normalization_ladder as E139
    with E139.e136() as F:
        old = F.seed_index
        F.seed_index = lambda w: seed_index(calyx, w)
        try:
            yield F
        finally:
            F.seed_index = old


def load_run(pt: Point):
    """A read's saved run (`read_dirs`), or None."""
    from kcq import stage_readings as R
    try:
        return R.load_step3_run(point_tag(pt), dirs=read_dirs())
    except FileNotFoundError:
        return None


def has_run(pt: Point) -> bool:
    """Whether a read's run is saved (its <tag>_run.json in `read_dirs`), without loading it."""
    return any((Path(d) / f'{point_tag(pt)}_run.json').exists() for d in read_dirs())


def expected_points(calyx: str, setting: str) -> List[Point]:
    """The reads a setting's class on a calyx waits for: its step-1 points (`planned_points`), and its step-2 points
    (D1.8 and D5.8 on a male calyx's permutations) once any of them is saved."""
    from experiments.E144_a1_malecns_calyces import planned_points
    pts = [p for p in planned_points(['1']) if p.calyx == calyx and p.setting == setting]
    if calyx in MALE and setting in ('D1.8', 'D5.8'):
        p2 = [p for p in planned_points(['2'], calyces=[calyx]) if p.setting == setting]
        if any(has_run(p) for p in p2):
            pts += p2
    return pts


def missing_points(calyx: str, setting: str) -> List[str]:
    """The tags of a setting's expected reads on a calyx (`expected_points`) that are not saved yet."""
    return [point_tag(p) for p in expected_points(calyx, setting) if not has_run(p)]


def read_wirings(calyx: str, setting: str, panel: str) -> List[str]:
    """The wirings of the calyx on which the setting's read of `panel` is saved, in WIRINGS_OF's order."""
    return [w for w in WIRINGS_OF[calyx] if has_run(Point(calyx, setting, w, panel))]


def hemibrain_spreads() -> dict:
    """The spreads of every class (the entry's): the hemibrain family's SD over its eight wirings on the measure
    ('spreads') and on every view ('spreads_all_views', rule 3's guard), as E141's score computed them (its classes.json;
    E136's on the measure beside, 'max_abs_diff_from_e136')."""
    from experiments import E136_stage_population_fit as F
    from experiments import E141_kckc_rise_forms as E141
    c = _read_json(E141.out_dir() / 'classes.json')
    if not c.get('spreads') or not (c.get('beside') or {}).get('spreads_all_views'):
        raise RuntimeError("E141's classes.json has no spreads")
    with_m = _read_json(H.results_dir(F.EXPERIMENT) / 'measure_strength' / 'classes.json').get('spreads') or {}
    diff = [abs(float(v) - float(with_m[k])) for k, v in c['spreads'].items() if k in with_m]
    return dict(spreads={k: float(v) for k, v in c['spreads'].items()},
                spreads_all_views={k: float(v) for k, v in c['beside']['spreads_all_views'].items()},
                source=str(E141.out_dir() / 'classes.json'), max_abs_diff_from_e136=max(diff) if diff else None)


def megamat_of(F, calyx: str, setting: str, ws: Sequence[str], stages: Sequence[str], half: Optional[str],
               views='measure') -> Optional[dict]:
    """A setting's megamat stages on the calyx's wirings `ws` against a half of the flies (all for None): E136's
    `StageLoss` (split) on the saved reads, through `e136_for`'s seeds; None while a read is missing."""
    from kcq import stage_readings as R
    runs = {w: load_run(Point(calyx, setting, w, 'megamat')) for w in ws}
    if not runs or any(r is None for r in runs.values()):
        return None
    return F.StageLoss(stages, F.Flies(half, stages, views=views), R.require_lock(), split=True)(runs)


def panels_of(F, calyx: str, setting: str, ws: Sequence[str], views='measure', flies=None) -> pd.DataFrame:
    """E136's `score_panels` on the calyx's wirings `ws`: the KC somata of the validation2, mixture and APL16 reads, each
    panel and condition with flies, every wiring's population error and the error of their mean views ('common'), with
    the views and the strength split."""
    from experiments import E134_population_baseline as B
    from kcq import population_fit as P
    from kcq import stage_readings as R
    flies = flies if flies is not None else R.FlyCache()
    lock = R.require_lock()
    views = F.MEASURES[F.MEASURE] if views == 'measure' else views
    rows = []
    for panel in ('validation2', 'mixtures', 'apl16'):
        runs = {w: load_run(Point(calyx, setting, w, panel)) for w in ws}
        if not runs or any(r is None for r in runs.values()):
            continue
        reads = {w: R.read_run(r, stages=['kc_soma'], lock=lock, flies=flies, seed_index=seed_index(calyx, w))
                 for w, r in runs.items()}
        for key in sorted({k for rd in reads.values() for k in rd.populations}):
            st, pan, cond = key
            fl = flies.get(key)
            if fl is None:
                continue
            pops = [reads[w].populations[key] for w in ws]
            fp, have = B.fly_side_read(fl)
            od = [o for o in pops[0].values.columns if o in set(have)]
            ref = P.fly_reference(fp, od, B.breadth_edges(len(od)), parts=F.measure_parts(views))
            tab = P.wirings_as_flies(pops, ref)
            names = {p_.name: w for p_, w in zip(pops, ws)}
            sp = F.strength_split(pops, ref).add_prefix('split_')
            for name, row in tab.iterrows():
                rows.append(dict(variant=variant(calyx, setting), calyx=calyx, setting=setting, input=panel, stage=st,
                                 panel=pan, condition=cond, wiring=names.get(name, str(name)),
                                 views='all' if views is None else '+'.join(views), n_odors=len(od), n_flies=len(fp),
                                 **{k: float(v) for k, v in row.items()},
                                 **({k: float(v) for k, v in sp.loc[name].items()} if name in sp.index else {})))
    return pd.DataFrame(rows)


def panel_values_of(rows: pd.DataFrame, v: str, ws: Sequence[str]) -> Dict[Tuple[str, str, str], Dict[str, float]]:
    """From `panels_of`' rows, per (input, panel, condition) of E136's SCORED_PANELS each wiring's KC population
    error."""
    from experiments import E136_stage_population_fit as F
    out = {}
    for inp, pan, cond in F.SCORED_PANELS:
        x = rows[(rows.variant == v) & (rows.input == inp) & (rows.panel == pan) & (rows.condition == cond)
                 & rows.wiring.isin(list(ws))] if len(rows) else rows
        out[(inp, pan, cond)] = {str(r.wiring): float(r.error) for r in x.itertuples()} if len(x) else {}
    return out


def apl16_of(F, calyx: str, setting: str, ws: Sequence[str]) -> Tuple[pd.DataFrame, dict, dict]:
    """A setting's APL16 reads on the calyx's wirings: per wiring E67's `full_row`, E136's `apl16_realism` (the bands)
    and the largest counts (as E136's `apl16_rows`), and E139's `apl_effect_parts` (the price's statistics, under the
    calibrated readout). Returns (rows, parts per wiring, realism flag per wiring)."""
    from experiments import E139_normalization_ladder as E139
    from experiments.E67_apl_effect_fits import full_row
    rows, parts, flags = [], {}, {}
    for w in ws:
        r = load_run(Point(calyx, setting, w, 'apl16'))
        if r is None:
            continue
        o = E139.run_outputs(r)
        row = dict(variant=variant(calyx, setting), calyx=calyx, setting=setting, wiring=w, **full_row(o))
        row.update(F.apl16_realism(o['on'][0], o['off'][0]))
        for c in ('on', 'off'):
            row[f'max_count_{c}'] = float(o[c][1].to_numpy(dtype=float).max())
        parts[w] = E139.apl_effect_parts(o)
        row.update({f'price_{k}': x for k, x in parts[w].items()})
        flags[w] = row['realism_flag'] if isinstance(row['realism_flag'], str) else ''
        rows.append(row)
    return pd.DataFrame(rows), parts, flags


def price_on(parts: Mapping[str, Mapping[str, Mapping[str, float]]], v: str, ref: str, ws: Sequence[str]) -> dict:
    """E139's price on the APL effect of a setting against the calyx's family on the wirings `ws` both have: per target
    of E139's PRICE_TARGETS the mean over them of the per-wiring changes of the calibrated pattern part; 'price' where
    one is above E139's PRICE."""
    from experiments import E139_normalization_ladder as E139
    a, f = parts.get(v, {}), parts.get(ref, {})
    ww = [w for w in ws if w in a and w in f]
    rec = dict(wirings=ww, complete=len(ww) == len(ws) and bool(ww), threshold=E139.PRICE, by_target={})
    if not ww:
        rec.update(price=None, price_value=None)
        return rec
    for t in E139.PRICE_TARGETS:
        ch = {w: a[w][f'{t}_pattern'] - f[w][f'{t}_pattern'] for w in ww}
        rec['by_target'][t] = dict(changes=ch, change=float(np.mean(list(ch.values()))),
                                   value=float(np.mean([a[w][f'{t}_pattern'] for w in ww])),
                                   family=float(np.mean([f[w][f'{t}_pattern'] for w in ww])))
    vals = [x['change'] for x in rec['by_target'].values()]
    rec.update(price=bool(max(vals) > E139.PRICE), price_value=float(max(vals)))
    return rec


def classify_on(F, v: str, ref: str, ws: Sequence[str], stages: Sequence[str], mega: Mapping, pv: Mapping,
                spreads: Mapping[str, float], flags: Mapping[str, Mapping[str, str]], guard: Optional[Mapping]) -> dict:
    """E136's `classify` of one setting against the calyx's family (`ref`) on the calyx's wirings `ws` it was read on
    (E136's `compare` over `ws`: the standard error the larger of half the spread and the SD of the changes over sqrt(n),
    so that a single wiring's bar is the spread; the spreads the hemibrain family's), on the held-out half: 'holds',
    'megamat only', 'one stage at another's cost', 'strength at another view's cost', 'does not improve', or 'no class'
    where a read is missing or a wiring of the setting leaves a realism band (any band; the family's flags beside).
    Rule 3's guard on every view as E136's. Nothing is fitted, so there is no bound and no permuted control."""
    a, f = mega.get(v, {}), mega.get(ref, {})
    m = F.compare({w: x['stage'] for w, x in a.items()}, {w: x['stage'] for w, x in f.items()}, spreads['stage'], ws)
    st = {s: F.compare({w: x[s] for w, x in a.items()}, {w: x[s] for w, x in f.items()}, spreads[s], ws) for s in stages}
    pn, pb = {}, {}
    for keys, out in ((F.HOLDS_PANELS, pn), (F.BESIDE_PANELS, pb)):
        for key in keys:
            k = '/'.join(key)
            out[k] = F.compare(pv.get(v, {}).get(key, {}), pv.get(ref, {}).get(key, {}), spreads.get(k, np.nan), ws)
    rec = dict(variant=v, reference=ref, wirings=list(ws), megamat=m, stages=st, panels=pn, panels_beside=pb)
    reasons = []
    if not m.get('complete') or not all(x.get('complete') for x in st.values()):
        reasons.append('reads missing')
    left = [w for w in ws if flags.get(v, {}).get(w)]
    if left:
        reasons.append('leaves a realism band on ' + ', '.join(left))
    rec['family_flagged'] = [w for w in ws if flags.get(ref, {}).get(w)]
    stage_worse = [s for s, x in st.items() if x.get('worse_than_spread')]
    panel_worse = [k for k, x in pn.items() if x.get('complete') and x.get('worse_than_spread')]
    panel_missing = [k for k, x in pn.items() if not x.get('complete')]
    guard_worse, guard_rec = [], {}
    if guard is not None:
        ga, gf = guard['mega'].get(v, {}), guard['mega'].get(ref, {})
        for s in stages:
            c = F.compare({w: x[s] for w, x in ga.items()}, {w: x[s] for w, x in gf.items()},
                          guard['spreads'].get(s, np.nan), ws)
            guard_rec[s] = c
            if c.get('worse_than_spread'):
                guard_worse.append(s)
        for key in F.SCORED_PANELS:
            k = '/'.join(key)
            c = F.compare(guard['panels'].get(v, {}).get(key, {}), guard['panels'].get(ref, {}).get(key, {}),
                          guard['spreads'].get(k, np.nan), ws)
            guard_rec[k] = c
            if key in F.HOLDS_PANELS and c.get('complete') and c.get('worse_than_spread') and k not in panel_worse:
                panel_worse.append(k)
    improves = bool(m.get('improves')) and not stage_worse and not guard_worse
    rec.update(improves=improves, stage_worse=stage_worse, guard_worse=guard_worse, guard=guard_rec,
               panel_worse=panel_worse, panel_missing=panel_missing, bar=2.0 * m['se'] if m.get('complete') else None)
    if reasons:
        cls = 'no class'
    elif improves and not panel_worse and not panel_missing:
        cls = 'holds'
    elif improves and panel_worse:
        cls = 'megamat only'
    elif improves:
        cls = 'improves (held-out panels not all read)'
    elif m.get('improves') and stage_worse:
        cls = "one stage at another's cost"
    elif m.get('improves') and guard_worse:
        cls = "strength at another view's cost"
    else:
        cls = 'does not improve'
    rec.update(cls=cls, reasons=reasons)
    return rec


def sites_of(calyx: str, setting: str, ws: Sequence[str]) -> pd.DataFrame:
    """Where APL's inhibition is made in a setting's megamat and APL16 reads (analyses/E136_writeup_tables.py's
    `sites_row` on the calyx's own claw x bouton wiring: the share of the claws' summed drive APL removes and the share of
    that made at the boutons), per wiring, beside the calyx's share of APL's calyx synapses at the boutons (E130's
    step0_apl.csv)."""
    from analyses import E136_writeup_tables as T
    from experiments import E130_malecns_family as E130
    ap = pd.read_csv(E130.out_dir() / 'step0_apl.csv').set_index('calyx')['bouton_share_apl_output']
    rows = []
    for panel in ('megamat', 'apl16'):
        for w in ws:
            r = load_run(Point(calyx, setting, w, panel))
            if r is None:
                continue
            on, off = (r.conditions[c] for c in ('on', 'off'))
            rows.append(dict(variant=variant(calyx, setting), calyx=calyx, setting=setting, panel=panel, wiring=w,
                             anatomical_bouton_share=float(ap.get(CALYCES[calyx], np.nan)),
                             **T.sites_row(tables(calyx, w)['W'], on, off, False)))
    return pd.DataFrame(rows)


def without_odor_of(F, calyx: str, setting: str, ws: Sequence[str], odor: str = WITHOUT_ODOR) -> pd.DataFrame:
    """The megamat somata's strength error with and without one odor's entry (E136's `stage_without_odor` on the
    calyx's wirings: E138's `odor_entries`, against the held-out half in all flies' units)."""
    from experiments import E138_strength_by_odor_entries as E138
    from kcq import stage_readings as R
    runs = {w: load_run(Point(calyx, setting, w, 'megamat')) for w in ws}
    if not runs or any(r is None for r in runs.values()):
        return pd.DataFrame()
    flies = F.Flies('held', ['kc_soma'])
    lock = R.require_lock()
    pops = [R.read_run(runs[w], stages=['kc_soma'], lock=lock, flies=flies.cache, seed_index=seed_index(calyx, w)
                       ).populations[('kc_soma', F.PANEL, 'on')] for w in ws]
    ref = flies.reference('kc_soma', list(pops[0].values.columns))
    ent = E138.odor_entries(pops, ref, names=list(ws))
    ent = ent[ent['odor'].astype(str) == odor]
    out = ent[['row', 'view_error', 'view_error_without', 'view_share']].copy()
    for k, x in (('odor', odor), ('setting', setting), ('calyx', calyx), ('variant', variant(calyx, setting))):
        out.insert(0, k, x)
    return out.reset_index(drop=True)


def e137_of(F, calyx: str, setting: str, ws: Sequence[str]) -> pd.DataFrame:
    """The mixture reads scored under E137's readings of the mixture flies (E136's `score_e137_beside` on the calyx's
    wirings, seeds by their order, which is E137's)."""
    from data_prep import stage_population_data as D
    from experiments import E137_heldout_readings as E137
    from kcq import population_fit as P
    runs = [load_run(Point(calyx, setting, w, 'mixtures')) for w in ws]
    if not runs or any(r is None for r in runs):
        return pd.DataFrame()
    lock = E137.load_lock()
    rows = []
    for panel, readings in F.E137_BESIDE.items():
        od = list(D.panel_odors(panel))
        for reading in readings:
            units = E137.PRIMARY_UNITS[reading]
            pop, noise = E137.reading_population(panel, reading, units, lock)
            fps = P.fly_populations(pop)
            pops = E137.read_model_panel(runs, panel, pop, noise)
            names = {p_.name: w for p_, w in zip(pops, ws)}
            for views in dict.fromkeys((F.MEASURES[F.MEASURE], None)):
                ref = P.fly_reference(fps, od, D.breadth_bins(len(od)), parts=F.measure_parts(views))
                for name, r in P.wirings_as_flies(pops, ref).iterrows():
                    rows.append(dict(variant=variant(calyx, setting), calyx=calyx, setting=setting, input='mixtures',
                                     stage='kc_soma', panel=panel, condition='on', reading=reading, units=units,
                                     views='all' if views is None else '+'.join(views),
                                     wiring=names.get(name, str(name)), n_odors=len(od), n_flies=len(fps),
                                     **{k: float(x) for k, x in r.items()}))
    return pd.DataFrame(rows)


def first_values(stages: Sequence[str]) -> Tuple[Dict[str, dict], Dict[str, dict], Dict[str, str]]:
    """The 'first' real wiring's values of each setting's counterpart (FIRST_COUNTERPART), from E141's and E139's score
    tables (the first that has it): its megamat stages against the held-out half ('stage' and each stage), its held-out
    panels' KC population errors (per E136 SCORED_PANELS key 'input/panel/condition'); and where each came from."""
    from experiments import E136_stage_population_fit as F
    mega, pan, src = {}, {}, {}
    for s, cp in FIRST_COUNTERPART.items():
        for exp in FIRST_SOURCES:
            d = H.results_dir(exp)
            if not (d / 'score_megamat.csv').exists():
                continue
            m = pd.read_csv(d / 'score_megamat.csv')
            x = m[(m.variant == cp) & (m.half == 'held') & (m.wiring == 'real')]
            if not len(x):
                continue
            r = x.iloc[-1]
            mega[s] = dict(stage=float(r['stage']), **{st: float(r[st]) for st in stages})
            p = pd.read_csv(d / 'score_panels.csv')
            pan[s] = {}
            for inp, pn_, cond in F.SCORED_PANELS:
                y = p[(p.arm == cp) & (p.input == inp) & (p.panel == pn_) & (p.condition == cond) & (p.wiring == 'real')]
                if len(y):
                    pan[s]['/'.join((inp, pn_, cond))] = float(y.iloc[-1]['error'])
            src[s] = f'{exp}:{cp}'
            break
    return mega, pan, src


def check_against_first(stages: Sequence[str]) -> pd.DataFrame:
    """The score's megamat reading checked where it overlaps E141's and E139's tables: each setting's counterpart on
    the 'first' real wiring (FIRST_COUNTERPART: E134's held run for the family, E139's and E141's reads, E136's A1 its
    fit's best evaluation on the S wiring 'real') read through this score's path (E136's `StageLoss` with the split,
    against the held-out half and all flies, inside `e136_for('H')`: the real wiring's seed index 0, as E136's), beside
    the value of the table `first_values` takes it from (score_megamat.csv, wiring 'real'); each difference should be
    0, a wiring's errors not depending on the wirings read with it. One row per counterpart, half and statistic."""
    from experiments import E139_normalization_ladder as E139
    from experiments import E141_kckc_rise_forms as E141
    from kcq import stage_readings as R
    _, _, src = first_values(stages)
    rows = []
    with e136_for('H') as F:
        for s, cp in FIRST_COUNTERPART.items():
            exp = (src.get(s) or ':').split(':')[0]
            if not exp:
                continue
            if cp not in ('family', 'A1'):
                try:
                    E141.parse_variant(cp)
                    E141.register_variants(F, [cp])
                except ValueError:
                    E139.register_variants(F, [cp])
            runs = F.megamat_runs(cp, ['real'], 'held')
            tab = pd.read_csv(H.results_dir(exp) / 'score_megamat.csv')
            for half in ('held', None):
                if runs is None:
                    rows.append(dict(setting=s, counterpart=cp, table=exp, half=half or 'all', statistic='stage',
                                     status='run not on this host'))
                    continue
                sc = F.StageLoss(stages, F.Flies(half, stages), R.require_lock(), split=True)(runs)
                got = F.megamat_values(sc, stages, ['real'])['real']
                x = tab[(tab.variant == cp) & (tab.half == (half or 'all')) & (tab.wiring == 'real')]
                for k in ('stage', *stages):
                    want = float(x.iloc[-1][k]) if len(x) else np.nan
                    rows.append(dict(setting=s, counterpart=cp, table=exp, half=half or 'all', statistic=k,
                                     score_path=float(got[k]), table_value=want, diff=float(got[k]) - want,
                                     status='read'))
    return pd.DataFrame(rows)


def hemibrain_family_range(stages: Sequence[str]) -> Dict[str, Tuple[float, float]]:
    """The hemibrain family's range over its eight wirings of each megamat stage's error against the held-out half (E141's
    score_megamat.csv, variant 'family'), for "the calyx's family is inside the hemibrain's wirings"."""
    from experiments import E136_stage_population_fit as F
    from experiments import E141_kckc_rise_forms as E141
    m = pd.read_csv(E141.out_dir() / 'score_megamat.csv')
    x = m[(m.variant == 'family') & (m.half == 'held') & m.wiring.isin(F.WIRINGS)]
    return {st: (float(x[st].min()), float(x[st].max())) for st in ('stage', *stages)}


def stages_read(F) -> Tuple[str, ...]:
    """The stages the score reads (E136's `stage_classes_on_measure` on the measure in force, written to this
    experiment's stage_classes.json); RuntimeError if a stage is not readable."""
    cl = F.stage_classes_on_measure(F.MEASURES[F.MEASURE])
    _write_json(out_dir() / 'stage_classes.json', dict(measure=F.MEASURE, classes=cl, time=_now(), host=_host(),
                                                       **F._commits()))
    bad = [s for s in F.LOSS_STAGES if cl[s]['cls'] in ('not readable', 'not scored')]
    if bad:
        raise RuntimeError(f'on the measure {F.MEASURE!r} the stages {bad} are not readable: {cl}')
    return tuple(F.LOSS_STAGES)


def _known(cls: Optional[str], value: bool) -> Optional[bool]:
    """A derived class from a setting's class: None while the setting's class is not read ('reads missing', absent)."""
    return None if cls in (None, 'reads missing') else bool(value)


def rise_needed(cls: Mapping[str, Optional[str]]) -> Optional[bool]:
    """"The rise is needed there": B does not hold and D1.8, A1 or D5.8 does; None while B's class is not read, or no
    rise holds and one of them is not read."""
    if cls.get('B') in (None, 'reads missing'):
        return None
    rises = [cls.get(s) for s in RISES]
    if any(x == 'holds' for x in rises):
        return cls.get('B') != 'holds'
    return None if any(x in (None, 'reads missing') for x in rises) else False


def outcome_letters(cal: Mapping[str, dict], h_default: Mapping[str, dict]) -> dict:
    """The entry's outcomes read off the calyces' classes: (a) A1 holds on both male calyces; (b) megamat falls with A1
    on both and A1 holds on neither; (c) A1 holds on one; (d) B holds on a male calyx; (e) a male calyx's family is
    outside the hemibrain's wirings at a megamat stage (beside any of the others); (f) the claw default changes a class
    on H. None where a calyx's classes are not all in."""
    def b(x):
        return None if x is None else bool(x)
    holds = {c: cal.get(c, {}).get('a1_holds') for c in MALE}
    falls = {c: cal.get(c, {}).get('a1_megamat_falls') for c in MALE}
    known = all(x is not None for x in holds.values())
    out = dict(a=b(known and all(holds.values())) if known else None,
               b=b(known and not any(holds.values()) and all(bool(falls[c]) for c in MALE)) if known else None,
               c=b(known and sum(bool(x) for x in holds.values()) == 1) if known else None,
               d=True if any(cal.get(c, {}).get('b_holds') is True for c in MALE) else
               None if any(cal.get(c, {}).get('b_holds') is None for c in MALE) else False,
               e=any(not all(cal.get(c, {}).get('family_inside', {}).get(st, {}).get('inside', True)
                             for st in FAMILY_INSIDE_STAGES) for c in MALE),
               f=any(bool(x.get('class_differs')) for x in h_default.values()))
    return out


def part_score() -> dict:
    """The score (design item 3, with item 4's readings beside), per calyx, each setting against the calyx's own family
    on the same wirings, on E136's adopted measure (strength by odor; E139's `e136`), its stages E136's readable ones:
    - the megamat stages on each setting's read wirings against the held-out half (deciding), all flies and the fit half
      beside (score_megamat.csv; per wiring, 'common' the wirings' mean views, 'mean' the mean of each wiring's errors);
    - the held-out panels' KC population errors on the measure and on every view (score_panels.csv,
      score_panels_all_views.csv), the kiwi mixtures beside; E137's readings of the mixture flies (score_panels_e137.csv);
    - the APL16 reads with E67's statistics, the realism bands and the price's statistics (score_apl16.csv);
    - E136's classes on the calyx's wirings (`classify_on`; the spreads the hemibrain family's, `hemibrain_spreads`; the
      calyx's own SD over its wirings beside), no class where a hold fails (E139's `holds_failed`) or a wiring leaves a
      band; E139's price against the calyx's family (`price_on`);
    - per calyx: "A1 replicates" (A1 holds), "the rise is needed" (B does not hold and D1.8, A1 or D5.8 does), "the
      family is inside the hemibrain's wirings" (its strength error on the real wiring, against the held-out half, inside
      the hemibrain family's range over its eight wirings at the boutons, the claws and the somata), the dose ladder
      (validation2's change along B, D1.8, A1 and D5.8 on the real wiring, beside the hemibrain's), the number of
      settings read and holding;
    - on H, "the claw default changes it": each setting's megamat stages and held-out panels against the 'first' real
      wiring's (`first_values`), differing where by more than the spread, and whether its class on the real wiring
      differs from the 'first' counterpart's against the 'first' family on the same values;
    - beside: where APL's inhibition is made (score_sites.csv), the megamat somata with and without aa's entry
      (score_without_aa.csv); the check of this path against E141's and E139's tables on the 'first' real wiring
      (`check_against_first`, score_check_first.csv); E143's intervals over flies are `part_resampled`'s.
    classes.json (with the outcomes' letters) and provenance_score.json. A setting whose expected reads
    (`expected_points`) are not all in, or whose calyx's family's are not, gets no class ('reads missing', with the
    missing tags); its values are scored on what is in."""
    from experiments import E139_normalization_ladder as E139
    from kcq import stage_readings as R
    t0 = time.time()
    require_lock()
    reads = reads_table()
    sp = hemibrain_spreads()
    spreads, g_spreads = sp['spreads'], sp['spreads_all_views']
    mrows, prow, grow, arows, srows, wrows, erows = [], [], [], [], [], [], []
    classes, calyces = {}, {}
    with E139.e136() as F0:
        F0.require_gate()
        stages = stages_read(F0)
    fam_range = hemibrain_family_range(stages)
    first_mega, first_pan, first_src = first_values(stages)
    check = check_against_first(stages)
    check.to_csv(out_dir() / 'score_check_first.csv', index=False)
    h_default = {}
    for c in CALYCES:
        with e136_for(c) as F:
            names = [s for s in SETTINGS if c in SETTINGS[s].calyces]
            ws_of = {s: read_wirings(c, s, 'megamat') for s in names}
            missing = {s: missing_points(c, s) for s in names}
            mega, mega_all, pv, g_mega, g_pv, parts, flags = {}, {}, {}, {}, {}, {}, {}
            for s in names:
                v, ws = variant(c, s), ws_of[s]
                if not ws:
                    continue
                for half in ('held', None, 'fit'):
                    sc = megamat_of(F, c, s, ws, stages, half)
                    if sc is None:
                        continue
                    vals = F.megamat_values(sc, stages, ws)
                    if half == 'held':
                        mega[v] = vals
                    elif half is None:
                        mega_all[v] = vals
                    base = dict(variant=v, calyx=c, setting=s, half=half or 'all')
                    mrows += [dict(base, wiring=w, **x) for w, x in vals.items()]
                    mrows.append(dict(base, wiring='mean', **pd.DataFrame([vals[w] for w in ws]).mean().to_dict()))
                    mrows.append(dict(base, wiring='common', stage=sc['stage_error'],
                                      downstream=sc['stage_error_downstream'], **{st: sc[f'{st}_error'] for st in stages}))
                sc = megamat_of(F, c, s, ws, stages, 'held', views=None)
                if sc is not None:
                    g_mega[v] = F.megamat_values(sc, stages, ws)
                pws = [w for w in ws if all(has_run(Point(c, s, w, p)) for p in ('validation2', 'mixtures', 'apl16'))]
                pr = panels_of(F, c, s, pws)
                gr = panels_of(F, c, s, pws, views=None)
                prow.append(pr)
                grow.append(gr)
                pv[v] = panel_values_of(pr, v, pws) if len(pr) else {}
                g_pv[v] = panel_values_of(gr, v, pws) if len(gr) else {}
                ar, parts[v], flags[v] = apl16_of(F, c, s, ws)
                arows.append(ar)
                srows.append(sites_of(c, s, ws))
                wrows.append(without_odor_of(F, c, s, ws))
                erows.append(e137_of(F, c, s, pws))
            ref = variant(c, 'family')
            guard = dict(mega=g_mega, panels=g_pv, spreads=g_spreads)
            own_sd = {k: F._sd(x['stage'] for x in mega.get(ref, {}).values()) for k in ('stage',)}
            own_sd.update({st: F._sd(x[st] for x in mega.get(ref, {}).values()) for st in stages})
            for s in names:
                v = variant(c, s)
                if s == 'family':
                    continue
                ws = [w for w in ws_of[s] if w in ws_of.get('family', [])]
                if missing[s] or missing['family'] or v not in mega or ref not in mega or not ws:
                    # a class reads only complete settings (every expected read of the setting and of the family)
                    classes[v] = dict(cls='reads missing', calyx=c, setting=s, reference=ref,
                                      n_reads=int((reads.get('variant', pd.Series(dtype=object)) == v).sum()),
                                      n_expected=len(expected_points(c, s)), missing=missing[s],
                                      family_missing=missing['family'])
                    continue
                rec = classify_on(F, v, ref, ws, stages, mega, pv, spreads, flags, guard)
                failed = E139.holds_failed(v, reads, {}, F.HOLD_TOL)
                if failed:
                    rec['reasons'].append('a hold fails: ' + ', '.join(failed))
                    rec['cls'] = 'no class'
                rec.update(holds_failed=failed, calyx=c, setting=s, flagged_wirings={w: f for w, f in flags.get(v, {})
                                                                                     .items() if f},
                           price_on_apl_effect=price_on(parts, v, ref, ws))
                classes[v] = rec
            cls = {s: classes.get(variant(c, s), {}).get('cls') for s in names if s != 'family'}
            ladder = {}
            for s in DOSES:
                k = 'validation2/validation2/on'
                key = ('validation2', 'validation2', 'on')
                a_, f_ = pv.get(variant(c, s), {}).get(key, {}), pv.get(ref, {}).get(key, {})
                ladder[s] = dict(change=float(a_['real'] - f_['real']) if 'real' in a_ and 'real' in f_ else None,
                                 hemibrain=HEMIBRAIN_LADDER[s], panel=k)
            inside = {}
            fr = mega.get(ref, {}).get('real')
            for st in ('stage', *stages):
                lo, hi = fam_range[st]
                val = None if fr is None else float(fr[st])
                inside[st] = dict(value=val, hemibrain_min=lo, hemibrain_max=hi,
                                  inside=None if val is None else bool(lo <= val <= hi))
            a1 = classes.get(variant(c, 'A1'), {})
            calyces[c] = dict(classes=cls, a1_holds=None if a1.get('cls') in (None, 'reads missing') else
                              bool(a1.get('cls') == 'holds'),
                              a1_megamat_falls=None if a1.get('cls') in (None, 'reads missing') else
                              (a1.get('megamat') or {}).get('improves'),
                              b_holds=_known(cls.get('B'), cls.get('B') == 'holds'),
                              rise_needed=rise_needed(cls),
                              family_inside=inside,
                              family_inside_all=all(bool(inside[st]['inside']) for st in FAMILY_INSIDE_STAGES
                                                    if inside[st]['inside'] is not None),
                              dose_ladder=ladder, own_sd_family=own_sd, settings_read=sum(1 for s in names if ws_of[s]),
                              settings_complete=[s for s in names if not missing[s]],
                              settings_holding=[s for s, x in cls.items() if x == 'holds'],
                              missing={s: missing[s] for s in names if missing[s]},
                              wirings={s: ws_of[s] for s in names})
            if c == 'H':
                fm = first_mega.get('family')
                for s in names:
                    v = variant(c, s)
                    hv = mega.get(v, {}).get('real')
                    fv = first_mega.get(s)
                    if hv is None or fv is None or missing[s] or missing['family']:
                        h_default[s] = dict(complete=False, first=first_src.get(s), missing=missing[s],
                                            family_missing=missing['family'])
                        continue
                    diffs = {k: dict(H=float(hv[k]), first=float(fv[k]), diff=float(hv[k] - fv[k]),
                                     spread=float(spreads[k]), differs=bool(abs(hv[k] - fv[k]) > spreads[k]))
                             for k in ('stage', *stages)}
                    hp = {'/'.join(key): x.get('real') for key, x in pv.get(v, {}).items()}
                    for k, x in first_pan.get(s, {}).items():
                        if hp.get(k) is not None:
                            diffs[k] = dict(H=float(hp[k]), first=float(x), diff=float(hp[k] - x),
                                            spread=float(spreads.get(k, np.nan)),
                                            differs=bool(abs(hp[k] - x) > spreads.get(k, np.inf)))
                    differs = [k for k, x in diffs.items() if x['differs']]
                    # the entry's class (a stage or panel differs by more than its spread) beside outcome (f)'s reading
                    # (a class differs, `class_differs`)
                    rec = dict(first=first_src.get(s), diffs=diffs, differs=differs, changes_it=bool(differs),
                               complete=True)
                    if s != 'family' and fm is not None:
                        # the class on one wiring with H's values and with the 'first' values (the same rule)
                        def one(mv, fmv, pvv, fpv):
                            m = F.compare({'real': mv['stage']}, {'real': fmv['stage']}, spreads['stage'], ['real'])
                            stw = [st for st in stages if F.compare({'real': mv[st]}, {'real': fmv[st]}, spreads[st],
                                                                    ['real']).get('worse_than_spread')]
                            pw = [k for k in ('/'.join(x) for x in F.HOLDS_PANELS) if k in pvv and k in fpv and
                                  F.compare({'real': pvv[k]}, {'real': fpv[k]}, spreads[k], ['real'])
                                  .get('worse_than_spread')]
                            imp = bool(m.get('improves')) and not stw
                            return 'holds' if imp and not pw else 'megamat only' if imp else 'does not improve'
                        href = {'/'.join(key): x.get('real') for key, x in pv.get(ref, {}).items()}
                        hc = one(hv, mega.get(ref, {}).get('real', fm), hp, href) if ref in mega else None
                        fc = one(fv, fm, first_pan.get(s, {}), first_pan.get('family', {}))
                        rec.update(class_on_H=hc, class_on_first=fc, class_differs=bool(hc is not None and hc != fc))
                    h_default[s] = rec
    cat = lambda xs: pd.concat([x for x in xs if len(x)], ignore_index=True, sort=False) if any(len(x) for x in xs) \
        else pd.DataFrame()
    pd.DataFrame(mrows).to_csv(out_dir() / 'score_megamat.csv', index=False)
    for name, xs in (('score_panels', prow), ('score_panels_all_views', grow), ('score_apl16', arows),
                     ('score_sites', srows), (f'score_without_{WITHOUT_ODOR}', wrows), ('score_panels_e137', erows)):
        cat(xs).to_csv(out_dir() / f'{name}.csv', index=False)
    out = dict(stages=list(stages), measure='strength', reference="each calyx's family (E130's fit)",
               spreads=spreads, spreads_all_views=g_spreads, spreads_source=sp['source'],
               spreads_max_abs_diff_from_e136=sp['max_abs_diff_from_e136'],
               price_rule=dict(threshold=E139.PRICE, targets=list(E139.PRICE_TARGETS),
                               decided_on="the calyx's wirings the setting was read on"),
               hemibrain_family_range=fam_range, first_counterparts=first_src, claw_default=h_default,
               check_first=dict(file='score_check_first.csv', n_read=int((check.get('status') == 'read').sum()),
                                max_abs_diff=float(check['diff'].abs().max()) if 'diff' in check and len(check)
                                else None),
               calyces=calyces, outcomes=outcome_letters(calyces, h_default), classes=classes,
               n_settings_read=sum(x['settings_read'] for x in calyces.values()),
               n_calyces_read=sum(1 for x in calyces.values() if x['settings_read']))
    _write_json(out_dir() / 'classes.json', out)
    _write_json(out_dir() / 'provenance_score.json', dict(
        part='score', script=SCRIPT, lock=_read_json(lock_path()).get('written'), e134_lock=R.require_lock().get('sha256'),
        seconds=round(time.time() - t0, 1), host=_host(), time=_now(), **H.code_commits()))
    return {v: x.get('cls') for v, x in classes.items()}


# ---------------------------------------------------------------------------
# E143's resampled reading of the class tables (design item 4; the lock's amendment 'e143_resampled')

ALL_WIRINGS: Tuple[str, ...] = tuple(WIRING_SEEDS)
AMENDMENT_ID = 'e143_resampled'
# the derived classes recomputed on each resample (classes.json's 'calyces' -> <calyx> -> key)
DERIVED: Tuple[str, ...] = ('a1_holds', 'a1_megamat_falls', 'b_holds', 'rise_needed')


def views_dir() -> Path:
    """outputs/E144_a1_malecns_calyces/views/: the reads read once against all the flies (E143's layout)."""
    return local_dir('views')


def resampled_names() -> List[str]:
    """Every variant of the design (calyx and setting), in CALYCES' and SETTINGS' order."""
    return [variant(c, s) for c in CALYCES for s in SETTINGS if c in SETTINGS[s].calyces]


def read_views(names: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Each variant's saved reads on each input read once against all the flies, as E143 reads saved runs (E143's
    `read_runs` inside its `views_at` on `views_dir`: E134's default readings with the curve, noise and criterion of
    all the stage's or panel's flies, this experiment's seed index through `e136_for`), one file per variant and input
    (E143's `views_path`; the flies' parts beside, E143's `flies_path`), read again only where its wirings changed.
    Returns a summary row per variant and input."""
    from experiments import E143_megamat_fly_sampling as E143
    want = set(names) if names is not None else None
    out = []
    with E143.views_at(views_dir()):
        for c in CALYCES:
            with e136_for(c):
                for s in SETTINGS:
                    if c not in SETTINGS[s].calyces or (want is not None and variant(c, s) not in want):
                        continue
                    v = variant(c, s)
                    for panel in PANELS:
                        ws = [w for w in WIRINGS_OF[c] if has_run(Point(c, s, w, panel))]
                        if not ws:
                            continue
                        path = E143.views_path(v, panel)
                        if path.exists() and sorted(set(pd.read_parquet(path, columns=['wiring'])['wiring'])) \
                                == sorted(ws):
                            out.append(dict(variant=v, input=panel, wirings=','.join(ws), status='saved before'))
                            continue
                        t0 = time.time()
                        runs = {w: load_run(Point(c, s, w, panel)) for w in ws}
                        df = E143.read_runs(v, panel, runs)
                        df['n_wirings'] = len(runs)
                        tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
                        df.to_parquet(tmp, index=False)
                        os.replace(tmp, path)
                        out.append(dict(variant=v, input=panel, wirings=','.join(ws), status='read',
                                        seconds=round(time.time() - t0, 1)))
                        print(f'E144 views {v} {panel}: {len(ws)} wirings, {time.time() - t0:.0f} s', flush=True)
    return pd.DataFrame(out)


def resampled_sides(names: Sequence[str], lock: Mapping) -> Dict[str, object]:
    """E143's sides (per megamat stage 'stage/megamat/on' and held-out panel 'kc_soma/<panel>/<condition>') of the
    variants `names` from `views_dir` (E143's `load_side`), each side's flies checked against E143's lock."""
    from experiments import E143_megamat_fly_sampling as E143
    keys = {E143.megamat_key(s): E143.megamat_key(s) for s in E143.LOSS_STAGES}
    keys.update({E143.side_key(k): k for k in E143.HOLDS_PANELS + E143.BESIDE_PANELS})
    out = {}
    with E143.views_at(views_dir()):
        for key, lock_key in keys.items():
            if not E143.flies_path(key).exists():
                raise RuntimeError(f'{E143.flies_path(key)}: no flies read for {key} (read_views first)')
            out[key] = E143.load_side(key, list(names))
            E143.check_side_flies(out[key], E143.lock_flies(lock, lock_key))
    return out


def restrict(rv, v: str, ref: str, ws: Sequence[str]):
    """E143's `Resampled` (on ALL_WIRINGS) kept to one comparison, as E143's `reclass` reads it: the variant `v` under
    its own name and the calyx's family `ref` as 'family' (reclass's reference), both on the wirings `ws`; a value
    either lacks is left out (its comparison is then not given)."""
    from experiments import E143_megamat_fly_sampling as E143
    idx = [ALL_WIRINGS.index(w) for w in ws]

    def pick(d: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
        return {v: d[v][:, idx], 'family': d[ref][:, idx]} if v in d and ref in d else {}
    return E143.Resampled(mega={k: pick(d) for k, d in rv.mega.items()}, guard={k: pick(d) for k, d in rv.guard.items()},
                          panels={k: pick(d) for k, d in rv.panels.items()},
                          guard_panels={k: pick(d) for k, d in rv.guard_panels.items()})


def resampled_classes(classes: Mapping[str, dict], rv) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, np.ndarray],
                                                                  Dict[str, dict]]:
    """Every class of classes.json's 'classes' through E143's resamples (`rv`, on ALL_WIRINGS, the last reference all
    the flies): E143's `class_interval_rows` on the comparison (`restrict`: the setting against its calyx's family on
    the wirings it was scored on; E144's records carry E136's comparisons, so E143's reclass reads them as they are),
    with no joint splits (megamat's 'stands' None). Returns (the class rows, the change rows, each variant's class per
    reference and its comparisons, for the derived classes)."""
    from experiments import E143_megamat_fly_sampling as E143
    krows, crows, cls_of, comps_of = [], [], {}, {}
    for v, rec in classes.items():
        c = rec.get('calyx') or str(v).split('_', 1)[0]
        ref = rec.get('reference') or variant(c, 'family')
        sub = restrict(rv, v, ref, rec.get('wirings') or []) if 'megamat' in rec else rv
        krow, rows = E143.class_interval_rows(EXPERIMENT, v, v, rec, sub, None)
        cls_of[v], comps_of[v] = E143.reclass(rec, v, sub, EXPERIMENT)
        label = dict(calyx=c, setting=rec.get('setting') or str(v).split('_', 1)[1], reference=ref,
                     wirings=','.join(rec.get('wirings') or []))
        krows.append(dict(label, **krow))
        crows += [dict(label, **r) for r in rows]
    kdf = pd.DataFrame(krows)
    if len(kdf):
        kdf['shares'] = kdf['shares'].map(json.dumps)
    return kdf, pd.DataFrame(crows), cls_of, comps_of


def derived_per_reference(cls_of: Mapping[str, np.ndarray], comps_of: Mapping[str, dict], K: int
                          ) -> Dict[str, Dict[str, Optional[np.ndarray]]]:
    """Per calyx the derived classes on every reference, as the score derives them from the classes (`part_score`):
    'a1_holds' (A1 holds), 'a1_megamat_falls' (A1's megamat stage error improves), 'b_holds' (B holds), 'rise_needed'
    (B does not hold and D1.8, A1 or D5.8 does); None where A1's (or B's) class is not read."""
    def holds(x: Optional[np.ndarray]) -> Optional[np.ndarray]:
        return None if x is None or x[0] == 'reads missing' else np.asarray(x == 'holds', bool)
    none = np.zeros(K, bool)
    out = {}
    for c in CALYCES:
        g = lambda s: holds(cls_of.get(variant(c, s))) if c in SETTINGS[s].calyces else None
        mm = (comps_of.get(variant(c, 'A1')) or {}).get('megamat')
        # as the score: a setting not read does not hold
        bh = g('B') if g('B') is not None else none
        rises = [x if x is not None else none for x in (g(s) for s in RISES)]
        out[c] = dict(a1_holds=g('A1'), a1_megamat_falls=None if mm is None else np.asarray(mm['improves'], bool),
                      b_holds=bh, rise_needed=~bh & np.any(rises, axis=0))
    return out


def outcomes_per_reference(der: Mapping[str, Mapping[str, Optional[np.ndarray]]], K: int
                           ) -> Dict[str, Optional[np.ndarray]]:
    """The outcomes (a) to (d) on every reference from the male calyces' derived classes (`outcome_letters`' rules);
    None where A1's class is not read on a male calyx. (e) and (f) are not class changes and are not resampled."""
    hold = [der.get(c, {}).get('a1_holds') for c in MALE]
    falls = [der.get(c, {}).get('a1_megamat_falls') for c in MALE]
    known = all(x is not None for x in hold)
    n = np.sum(hold, axis=0) if known else None
    d = np.any([der.get(c, {}).get('b_holds') if der.get(c, {}).get('b_holds') is not None else np.zeros(K, bool)
                for c in MALE], axis=0)
    return dict(a=(n == len(MALE)) if known else None,
                b=((n == 0) & np.all([f if f is not None else np.zeros(K, bool) for f in falls], axis=0))
                if known else None,
                c=(n == 1) if known else None, d=d)


def _share_of(x: Optional[np.ndarray], scored) -> dict:
    """A derived class through the resamples: its value against all the flies (the last reference), the share of
    resamples on which it is True and the share giving the scored value."""
    if x is None:
        return dict(scored=scored, all_flies=None, share_true=None, share_scored=None)
    x = np.asarray(x, bool)
    return dict(scored=scored, all_flies=bool(x[-1]), share_true=float(np.mean(x[:-1])),
                share_scored=None if scored is None else float(np.mean(x[:-1] == bool(scored))))


def shortcut_check(classes: Mapping[str, dict], e143_lock: Mapping, sides: Mapping[str, object]) -> pd.DataFrame:
    """E143's read-once check (its design item 1) on this experiment's reads: the shortcut (the model read once against
    all the flies, a half changing only the reference) against E136's own held-out half of each megamat stage (E143's
    `e136_split`; the held-out panels against all their flies, as scored), beside each class's scored values (its
    comparison's 'arm' and 'ref', the means over the wirings scored); within E143's CHECK_TOL it stands for the score's
    reading. A description of the method: it gates nothing."""
    from experiments import E143_megamat_fly_sampling as E143
    W = {}
    for st in E143.LOSS_STAGES:
        n = len(E143.lock_flies(e143_lock, E143.megamat_key(st)))
        W[E143.megamat_key(st)] = E143.half_weights(n, [E143.e136_split(n)[1]])
    for k in E143.HOLDS_PANELS + E143.BESIDE_PANELS:
        n = len(E143.lock_flies(e143_lock, k))
        W[k] = E143.half_weights(n, [tuple(range(n))])
    rv = E143.weighted_values(e143_lock, W, side_of=sides.__getitem__, wirings=ALL_WIRINGS)
    rows = []
    for v, rec in classes.items():
        if 'megamat' not in rec or not (rec.get('megamat') or {}).get('complete'):
            continue
        ref, ws = rec['reference'], list(rec['wirings'])
        idx = [ALL_WIRINGS.index(w) for w in ws]
        comps = [('stage', rec['megamat']), *((st, x) for st, x in rec['stages'].items()),
                 *((k, x) for k, x in rec['panels'].items())]
        for what, sc in comps:
            if not sc.get('complete'):
                continue
            src = rv.mega.get(what) if what in rv.mega else rv.panels.get(what, {})
            for side_name, nm in (('arm', v), ('ref', ref)):
                got = float(np.mean(src[nm][0, idx])) if nm in src else np.nan
                rows.append(dict(variant=v, statistic=what, side=side_name, name=nm, scored=float(sc[side_name]),
                                 shortcut=got, diff=got - float(sc[side_name])))
    df = pd.DataFrame(rows).drop_duplicates(['statistic', 'name', 'side', 'variant'])
    if len(df):
        df['within_tol'] = df['diff'].abs() <= E143.CHECK_TOL
    return df


def amended_lock() -> dict:
    """E144's lock as amended for the resampled reading: `require_lock` (results, its committed copy and the code
    agree), the amendment AMENDMENT_ID in it, and the committed copy tracked with no uncommitted change (E140's
    `_committed`), so that no resampled number is read before the amendment is committed; RuntimeError otherwise."""
    from experiments import E140_aa_claw_soma_step as E140
    from experiments.E144_a1_malecns_calyces import committed_lock_path
    lk = require_lock()
    if not any(isinstance(a, dict) and a.get('id') == AMENDMENT_ID for a in lk.get('amendments') or []):
        raise RuntimeError(f"the lock has no amendment {AMENDMENT_ID!r}: write it (--part lock) and commit it first")
    if not E140._committed(committed_lock_path()):
        raise RuntimeError(f'{committed_lock_path()} is not committed (or has uncommitted changes)')
    return lk


def part_resampled() -> dict:
    """E143's intervals over flies on every change of the score's class tables (design item 4, the lock's amendment
    'e143_resampled'; it gates nothing): the reads read once against all the flies (`read_views`), E143's resamples of
    its lock's flies (`resamples`: N_BOOT, seed BOOT_SEED, megamat's four stages and the held-out panels from one
    generator) as references with all the flies last (`resampled_values` on this experiment's sides and wirings), and
    every class of classes.json recomputed on each (`resampled_classes`), with the derived classes and outcomes (a) to
    (d) (`derived_per_reference`, `outcomes_per_reference`); beside, E143's read-once check on these reads
    (`shortcut_check`). Writes e143_intervals_changes.csv, e143_class_shares.csv, e143_shortcut_check.csv,
    e143_intervals.json and provenance_resampled.json beside the score's tables."""
    from experiments import E143_megamat_fly_sampling as E143
    t0 = time.time()
    lk = amended_lock()
    cl = _read_json(out_dir() / 'classes.json')
    if not cl.get('classes'):
        raise RuntimeError('no score yet (classes.json): run --part score first')
    e143_lock = E143.require_lock()
    names = resampled_names()
    summary = read_views(names)
    sides = resampled_sides(names, e143_lock)
    draws = E143.resamples(e143_lock)
    rv = E143.resampled_values(e143_lock, draws, side_of=sides.__getitem__, wirings=ALL_WIRINGS)
    K = int(next(iter(rv.mega['stage'].values())).shape[0])
    kdf, cdf, cls_of, comps_of = resampled_classes(cl['classes'], rv)
    chk = shortcut_check(cl['classes'], e143_lock, sides)
    chk.to_csv(out_dir() / 'e143_shortcut_check.csv', index=False)
    der = derived_per_reference(cls_of, comps_of, K)
    oc = outcomes_per_reference(der, K)
    scored_c = cl.get('calyces') or {}
    derived = {c: {k: _share_of(der[c][k], (scored_c.get(c) or {}).get(k)) for k in DERIVED} for c in CALYCES}
    outcomes = {k: _share_of(x, (cl.get('outcomes') or {}).get(k)) for k, x in oc.items()}
    cdf.to_csv(out_dir() / 'e143_intervals_changes.csv', index=False)
    kdf.to_csv(out_dir() / 'e143_class_shares.csv', index=False)
    out = dict(amendment=AMENDMENT_ID, gates='nothing', n_boot=int(K - 1), boot_seed=E143.BOOT_SEED, ci=list(E143.CI),
               class_stands_share=E143.CLASS_SHARE, wirings=list(ALL_WIRINGS),
               reference="E143's resamples of all the stage's or panel's flies (the last reference all the flies, "
                         "unresampled); the model read once against all the flies",
               scored_against='the held-out half (classes.json), as scored', n_changes=int(len(cdf)),
               n_classes=int(len(kdf)), derived=derived, outcomes=outcomes,
               shortcut_check=dict(file='e143_shortcut_check.csv', tol=E143.CHECK_TOL,
                                   megamat_stage_max_abs_diff=float(chk[chk.statistic == 'stage']['diff'].abs().max())
                                   if len(chk) else None,
                                   any_max_abs_diff=float(chk['diff'].abs().max()) if len(chk) else None,
                                   megamat_stage_within_tol=bool(chk[chk.statistic == 'stage']['within_tol'].all())
                                   if len(chk) else None),
               e143_lock_sha256=e143_lock.get('sha256'), views=str(views_dir()),
               views_summary=json.loads(summary.to_json(orient='records')) if len(summary) else [])
    _write_json(out_dir() / 'e143_intervals.json', out)
    _write_json(out_dir() / 'provenance_resampled.json', dict(
        part='resampled', script='experiments/E144_score.py', lock=lk.get('written'), lock_amended=lk.get('amended'),
        e143_lock_sha256=e143_lock.get('sha256'), seconds=round(time.time() - t0, 1), host=_host(), time=_now(),
        **H.code_commits()))
    return dict(n_changes=len(cdf), n_classes=len(kdf), outcomes={k: v['share_scored'] for k, v in outcomes.items()})
