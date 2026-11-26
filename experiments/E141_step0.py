#!/usr/bin/env python3
"""E141's step 0, no simulation (experiments/E141_kckc_rise_forms.py `--part step0`, the entry's design item 1): the
saved runs E139's score reads, read again for the unrefitted family, E139's F, B and K at 1000 and E136's A1.

- `soma_splits`: each panel's KC somata, read as E136's score reads them, with the strength error split into level,
  spread and order by E136's `strength_split`, on all flies and on each half of the flies, per wiring and over the
  wirings, with each odor's entry (E138's `odor_entries`). `spread_table` pivots it into the entry's lead table (the
  model's spread of the odors' response rates as a multiple of the flies', mean over T).
- `holds`: what the holds reach per variant, panel and wiring: the somatic threshold and the scales with APL on and
  APL off, the held response rate against its target, the response rates by KC type, the largest spike counts and the
  number of odors each KC responds to in the realism bins.
- `a1_by_kckc`: A1's saved evaluations binned by the KC-KC factor (`kckc_bins`), built as E139's `a1_by_rung` bins them
  by the additive term: per bin and half, the best evaluation on the fit half's loss and the bin's medians, each stage's
  error beside. A lead: the other factors differ between evaluations.

The runs are those E139's score takes (E139's `register_variants` with E136's `load_read` and `megamat_runs`): the
family's megamat and mixture runs are E134's held runs, its APL16 runs E90's part B and its validation2 runs E136's
reads; F, B and K are E139's reads, F's APL16 runs E102's saved ladder; A1 is E136's reads, its S wirings on megamat its
fit's best evaluation. The measure is E136's adopted one (strength by odor), set by E139's `e136`.

Judgments the entry leaves open (the narrowest reading):
- the panels are exactly those E136's score reads at the KC somata: megamat from the megamat input (E136's
  `score_megamat`'s soma stage) and every (input, panel, condition) of E136's `score_panels` with flies to read it
  against: validation2, the kiwi ('kmix') and control ('cmix') mixtures as uploaded (the default reading), the APL16
  panel with APL on against the control flies and APL off against APL output blocked (TNT), and the APL16 input's 12
  megamat odors against Remy's megamat flies;
- "all flies" is E136's reading (`score_panels`: the model read against all the panel's flies, the reference those of
  E134's `fly_side_read`, so validation2's six flies on its 15 odors). A half is E136's split (`split_names`: hal's
  `split_halves(n, 1, 0)`, 'fit' the first part and 'held' the rest) of those flies, read as E136 reads megamat against
  a half (`Flies`): the model read against the half alone (curve, noise and criterion), the reference the half's mean
  views in all the flies' one-fly units (`subset_reference(..., keep_units=True)`). On megamat that is E136's own split;
  on the held-out panels nothing was fitted to either half, and 'fit' only names the first part. E138 read its halves
  with each half's own one-fly units; the spread as a multiple of the flies' and r do not depend on the unit, the
  level, spread and order do. The APL16 panel's halves split recordings (hemispheres), as E136's populations list
  them, so a fly can have a hemisphere in each half; the control half is paired with the TNT half of the same name;
- the rows over wirings: 'mean' (all eight), 'mean_S' and 'mean_T' (the means of each wiring's value) and 'common'
  (the eight wirings' mean views, E136's); each odor's entry per wiring, 'common' and 'mean_T';
- the odors per KC are counted on each panel's own odors in the fixed realism bins (0, 1, 2-3, 4-7, 8 or more;
  kcq.apl_effect's `odors_per_kc`) whatever the panel's size, so they are comparable with the bands on the APL16 panel
  only (`opk_tv`, against the flies in baseline-noise units, beside the band);
- the KC-KC bins: N_KCKC_BINS equal widths in ln of the factor from x1 to the largest factor evaluated (the last edge
  inclusive), and one bin below x1 for the evaluations there (`kckc_bins`).
"""
from __future__ import annotations

import contextlib
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

# the family, E139's three settings at 1000 and E136's A1
FAMILY = 'family'
A1 = 'A1'
VARIANTS: Tuple[str, ...] = (FAMILY, 'F_inhadd1000', 'B_inhadd1000', 'K_inhadd1000', A1)
HALVES: Tuple[str, ...] = ('all', 'fit', 'held')
# the summary rows over wirings (beside each wiring's): of the panel statistics, of each odor's entry
SUMMARY_ROWS: Tuple[str, ...] = ('mean', 'mean_S', 'mean_T', 'common')
ODOR_SUMMARY_ROWS: Tuple[str, ...] = ('mean_T', 'common')
# each odor's entry (E138's `odor_entries` columns): the model's entry and the flies' mean entry in each strength part,
# the residual in one-fly units, and the odor's share of the view's squared error
ODOR_STATS: Tuple[str, ...] = ('rate_model', 'rate_flies', 'rate_e', 'mean_model', 'mean_flies', 'mean_e', 'view_share')
# the KC types of the responses' index level and their column suffixes in `holds`
KC_TYPES: Dict[str, str] = {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}
N_KCKC_BINS = 6
STAGE_COLS: Tuple[str, ...] = ('stage_error', 'stage_error_downstream', 'pn_dendrite_error', 'pn_bouton_error',
                               'kc_claw_error', 'kc_soma_error')


def _is_e139(name: str) -> bool:
    from experiments import E139_normalization_ladder as L
    try:
        L.parse_variant(name)
    except ValueError:
        return False
    return True


@contextlib.contextmanager
def e136(variants: Sequence[str] = VARIANTS) -> Iterator[object]:
    """E136's module on E139's measure (E139's `e136`), with E139's variants among `variants` registered with its score
    functions (E139's `register_variants`)."""
    from experiments import E139_normalization_ladder as L
    with L.e136() as F:
        L.register_variants(F, [v for v in variants if _is_e139(v)])
        yield F


def variant_runs(F, name: str, panel: str, wirings: Optional[Sequence[str]] = None) -> Optional[Dict[str, object]]:
    """A variant's saved runs on one input, per wiring, as E136's score takes them: megamat by `megamat_runs` (a fitted
    arm's S wirings its fit's best evaluation), every other input by `load_read`; None while one is missing."""
    wirings = tuple(wirings or F.WIRINGS)
    if panel == F.PANEL:
        return F.megamat_runs(name, wirings, 'held')
    runs = {w: F.load_read(F.ReadPoint(name, panel, w, None)) for w in wirings}
    return None if any(r is None for r in runs.values()) else runs


# ---------------------------------------------------------------------------
# the flies and their halves

def panel_half_names(fl) -> Dict[str, List[str]]:
    """E136's split (`split_names`) of the flies a recorded population is read on (E134's `fly_side_read`: all its
    flies, validation2's six), {'fit': [...], 'held': [...]}."""
    from experiments import E134_population_baseline as B
    from experiments import E136_stage_population_fit as F
    from kcq import stage_readings as R
    fp, _ = B.fly_side_read(fl)
    names = [p.name for p in fp]
    if names == [str(f) for f in fl.flies]:
        return F.split_names(fl)
    return F.split_names(R.subset_flies(fl, names))


class HalfFlies:
    """A FlyCache-like mapping (stage, panel, condition) -> one half of `base`'s flies (`panel_half_names`), each
    subset made once so that its noise model is computed once (kcq.stage_readings caches it on the population)."""

    def __init__(self, base, half: str):
        if half not in ('fit', 'held'):
            raise ValueError(f'{half=}')
        self.base, self.half = base, half
        self._got: Dict[tuple, object] = {}

    def get(self, key, default=None):
        from kcq import stage_readings as R
        if key not in self._got:
            fl = self.base.get(key)
            self._got[key] = None if fl is None else R.subset_flies(fl, panel_half_names(fl)[self.half])
        return self._got[key] if self._got[key] is not None else default


class References:
    """The flies' references of the KC somata as E136's `score_panels` builds them (all the flies of E134's
    `fly_side_read` on the model's odors they have, the measure's parts), and of a half: the half's mean views in all
    the flies' one-fly units (`subset_reference(..., keep_units=True)`, as E136's `Flies` reads megamat). Cached."""

    def __init__(self, base, views: Optional[Sequence[str]]):
        self.base, self.views = base, views
        self._full: Dict[tuple, object] = {}

    def __call__(self, key, odors: Sequence[str], half: str = 'all'):
        from experiments import E134_population_baseline as B
        from experiments import E136_stage_population_fit as F
        from kcq import population_fit as P
        fl = self.base.get(key)
        if fl is None:
            return None
        k = (key, tuple(odors))
        if k not in self._full:
            fp, have = B.fly_side_read(fl)
            od = [o for o in odors if o in set(have)]
            self._full[k] = P.fly_reference(fp, od, B.breadth_edges(len(od)), parts=F.measure_parts(self.views))
        ref = self._full[k]
        if half == 'all':
            return ref
        return P.subset_reference(ref, panel_half_names(fl)[half], keep_units=True)


# ---------------------------------------------------------------------------
# the somata's strength split

def _renamed(pops: Sequence, names: Sequence[str]) -> list:
    from kcq import population_fit as P
    return [P.Population(p.values, p.responds, name=str(n), named_units=p.named_units, graded=p.graded)
            for p, n in zip(pops, names)]


def _with_means(df: pd.DataFrame, wirings: Sequence[str], s_wirings: Sequence[str], t_wirings: Sequence[str],
                rows: Sequence[str]) -> pd.DataFrame:
    """`df` (one row per wiring, plus E136's 'mean' and 'common' where it has them) with the summary rows `rows`:
    'mean', 'mean_S', 'mean_T' the means of each wiring's row, 'common' kept as given."""
    out = [df.loc[list(wirings)]]
    for r in rows:
        if r == 'common':
            if 'common' in df.index:
                out.append(df.loc[['common']])
            continue
        ws = {'mean': wirings, 'mean_S': [w for w in wirings if w in s_wirings],
              'mean_T': [w for w in wirings if w in t_wirings]}[r]
        if ws:
            out.append(df.loc[list(ws)].mean().to_frame(r).T)
    return pd.concat(out)


def split_rows(pops: Sequence, ref, wirings: Sequence[str], s_wirings: Sequence[str], t_wirings: Sequence[str],
               odor_entries: bool = True) -> pd.DataFrame:
    """One panel and condition's rows of `soma_splits` (long: wiring, statistic, odor, value): per wiring and
    SUMMARY_ROWS, E136's population error and its strength view and parts (hal's `wirings_as_flies`: 'error',
    'view_strength', 'part_strength_rate', 'part_strength_mean') and E136's `strength_split` ('strength_<part>_<term>',
    the terms level, spread, order, sd_ratio and r; 'strength_rate_sd_ratio' is the spread of the odors' response rates
    as a multiple of the flies'); with `odor_entries`, per wiring and ODOR_SUMMARY_ROWS each odor's ODOR_STATS (E138's
    `odor_entries`).

    Args:
        pops: hal's `Population`s, one per wiring, in the order of `wirings`.
        ref: the flies' reference they are scored against.
        wirings, s_wirings, t_wirings: the wirings' names, and those of S and T among them.
        odor_entries: add each odor's entry.
    """
    from experiments import E136_stage_population_fit as F
    from experiments import E138_strength_by_odor_entries as E138
    from kcq import population_fit as P
    wirings = [str(w) for w in wirings]
    pops = _renamed(pops, wirings)
    tab = P.wirings_as_flies(pops, ref)
    tab = tab[[c for c in tab.columns if c == 'error' or c.startswith(('view_', 'part_'))]]
    sp = F.strength_split(pops, ref)
    wide = _with_means(pd.concat([tab, sp], axis=1), wirings, s_wirings, t_wirings, SUMMARY_ROWS)
    long = wide.rename_axis('wiring').reset_index().melt(id_vars='wiring', var_name='statistic', value_name='value')
    long.insert(2, 'odor', None)
    parts = [long]
    if odor_entries:
        ent = E138.odor_entries(pops, ref, names=wirings)
        per = []
        for odor, g in ent.groupby('odor', sort=False):
            g = g.set_index('row')[list(ODOR_STATS)]
            g = _with_means(g, wirings, s_wirings, t_wirings, ODOR_SUMMARY_ROWS)
            g = g.rename_axis('wiring').reset_index().melt(id_vars='wiring', var_name='statistic', value_name='value')
            g.insert(2, 'odor', str(odor))
            per.append(g)
        parts += per
    return pd.concat(parts, ignore_index=True)


def soma_splits(variants: Sequence[str] = VARIANTS, halves: Sequence[str] = HALVES,
                inputs: Optional[Sequence[str]] = None, odor_entries: bool = True) -> pd.DataFrame:
    """Step 0's first table: the KC somata of every panel E136's score reads (the module docstring), each variant's
    saved runs on the eight wirings read on all flies and on each half (`HalfFlies`, `References`), scored by
    `split_rows`.

    Args:
        variants: the variants (E139's names, 'family', 'A1').
        halves: of HALVES.
        inputs: the inputs read (default E136's READ_PANELS).
        odor_entries: add each odor's entry.

    Returns:
        Long: variant, input, panel, condition, half, wiring (a wiring, 'mean', 'mean_S', 'mean_T' or 'common'),
        statistic, odor (None for the panel's statistics), value, n_odors, n_flies (the half's). A variant missing a run
        on an input has no rows for it (printed).
    """
    from kcq import stage_readings as R
    out = []
    with e136(variants) as F:
        lock = R.require_lock()
        base = R.FlyCache()
        caches = {'all': base, **{h: HalfFlies(base, h) for h in halves if h != 'all'}}
        refs = References(base, F.MEASURES[F.MEASURE])
        for v in variants:
            for inp in (inputs or F.READ_PANELS):
                runs = variant_runs(F, v, inp)
                if runs is None:
                    print(f'E141 step 0: {v} has no saved run on {inp} for every wiring; skipped', flush=True)
                    continue
                for half in halves:
                    reads = {w: R.read_run(r, stages=['kc_soma'], lock=lock, flies=caches[half],
                                           seed_index=F.seed_index(w)) for w, r in runs.items()}
                    keys = sorted({k for rd in reads.values() for k in rd.populations})
                    for key in keys:
                        if base.get(key) is None or not all(key in rd.populations for rd in reads.values()):
                            continue
                        pops = [reads[w].populations[key] for w in runs]
                        ref = refs(key, list(pops[0].values.columns), half)
                        df = split_rows(pops, ref, list(runs), F.S_WIRINGS, F.T_WIRINGS, odor_entries)
                        _, panel, cond = key
                        for c, x in (('n_flies', len(ref.flies)), ('n_odors', len(ref.odors)), ('half', half),
                                     ('condition', cond), ('panel', panel), ('input', inp), ('variant', v)):
                            if c in ('n_flies', 'n_odors'):
                                df[c] = x
                            else:
                                df.insert(0, c, x)
                        out.append(df)
    cols = ['variant', 'input', 'panel', 'condition', 'half', 'wiring', 'statistic', 'odor', 'value', 'n_odors',
            'n_flies']
    return pd.concat(out, ignore_index=True)[cols] if out else pd.DataFrame(columns=cols)


def spread_table(splits: pd.DataFrame, half: str = 'all', wiring: str = 'mean_T',
                 statistic: str = 'strength_rate_sd_ratio') -> pd.DataFrame:
    """`soma_splits`'s panel statistic as the entry's lead table: (input, panel, condition) x variant, by default the
    model's spread of the odors' response rates as a multiple of the flies', the mean over T, on all flies."""
    s = splits[(splits.half == half) & (splits.wiring == wiring) & (splits.statistic == statistic)
               & splits.odor.isna()]
    t = s.pivot_table(index=['input', 'panel', 'condition'], columns='variant', values='value', aggfunc='first')
    return t[[v for v in VARIANTS if v in t.columns] + [v for v in t.columns if v not in VARIANTS]]


# ---------------------------------------------------------------------------
# what the holds reach

def _params(F, run) -> Dict[str, Dict[str, float]]:
    """{'on', 'off'}: a run's absolute parameters (F.VARIANT_PARAMS): its info's 'params' / 'params_off' (E134's,
    E136's, E139's and E102's runs); for E90's part B (the family on the APL16 input, which records none) its wiring's
    fit's (E83's `saved_fit`), APL off with E136's OFF_ZERO scales at 0."""
    info = getattr(run, 'info', None) or {}
    if info.get('params'):
        on = dict(info['params'])
        off = dict(info.get('params_off') or {k: (0.0 if k in F.OFF_ZERO else x) for k, x in on.items()})
    elif str(getattr(run, 'source', '')).startswith('E90:'):
        from experiments import E83_cross_wiring_fit as E83
        on = dict(E83.saved_fit(run.wiring)['p'])
        off = {k: (0.0 if k in F.OFF_ZERO else x) for k, x in on.items()}
    else:
        on = off = {}
    return {c: {k: float(p[k]) if k in p else np.nan for k in F.VARIANT_PARAMS} for c, p in (('on', on), ('off', off))}


def _bounds(run, reads: Mapping[str, Mapping]) -> Dict[str, float]:
    """Whether the hold ended at a bracket's bound (the threshold's with APL off, the APL->KC scale's with APL on): from
    E139's or E136's reads tables by the run's tag, else its info (E102's points, E136's fit evaluations); NaN where
    nothing records it (E134's held runs, E90's)."""
    tag = str(getattr(run, 'tag', ''))
    info = getattr(run, 'info', None) or {}
    for tab in reads.values():
        if tag in tab:
            r = tab[tag]
            return dict(at_bound_thr=float(str(r.get('at_bound_thr')).strip().lower() in ('true', '1', '1.0')),
                        at_bound_apl=float(str(r.get('at_bound_apl')).strip().lower() in ('true', '1', '1.0')))
    if 'at_bound_thr' in info or 'at_bound_apl' in info:
        return dict(at_bound_thr=float(bool(info.get('at_bound_thr'))),
                    at_bound_apl=float(bool(info.get('at_bound_apl'))))
    h = info.get('hold')
    if isinstance(h, Mapping) and isinstance(h.get('off'), Mapping) and isinstance(h.get('on'), Mapping):
        return dict(at_bound_thr=float(bool(h['off'].get('at_end'))), at_bound_apl=float(bool(h['on'].get('at_end'))))
    return dict(at_bound_thr=np.nan, at_bound_apl=np.nan)


def panel_holds(responses: pd.DataFrame, counts: pd.DataFrame, odors: Sequence[str]) -> Dict[str, float]:
    """One condition's somata on one panel's odors: the response rate overall and by KC type ('rate', 'rate_<type>' with
    KC_TYPES' suffixes, the responses' KC-type index level), the largest spike count ('max_count'), the fractions of KCs
    responding to 0, 1, 2-3, 4-7 and 8 or more of the odors ('opk_<bin>', kcq.apl_effect's `odors_per_kc`), and the
    numbers of KCs and odors."""
    from kcq import apl_effect as A
    from kcq import stage_readings as R
    r = R._abbrev(responses)
    c = R._abbrev(counts)
    od = [o for o in odors if o in set(r.columns)]
    r, c = r[od], c[od]
    resp = r.to_numpy(dtype=float) > 0.5
    out = dict(n_kcs=int(r.shape[0]), n_odors=len(od), rate=float(resp.mean()),
               max_count=float(np.nanmax(c.to_numpy(dtype=float))))
    types = r.index.get_level_values(H.KC_TYPE_LEVEL).astype(str) if H.KC_TYPE_LEVEL in (r.index.names or []) else None
    for t, suf in KC_TYPES.items():
        m = np.asarray(types == t) if types is not None else np.zeros(r.shape[0], dtype=bool)
        out[f'rate_{suf}'] = float(resp[m].mean()) if m.any() else np.nan
        out[f'n_kcs_{suf}'] = int(m.sum())
    out.update({f'opk_{b}': float(x) for b, x in zip(A.OPK_BINS, A.odors_per_kc(r))})
    return out


def run_holds(F, run, panel: str, targets: Optional[Tuple[float, float]] = None,
              reads: Optional[Mapping[str, Mapping]] = None) -> List[dict]:
    """A run's rows of `holds`: per panel it is read on at the somata (kcq.stage_readings' `run_panels`) and condition,
    the parameters (`_params`), the response rate over all the input's odors (what the hold holds, 'rate_held') beside
    its target ('frac_target', from `targets` = (APL off, APL on) fractions) and their relative miss, the hold's bounds
    (`_bounds`), `panel_holds`, and on the APL16 panel the odors-per-KC distance from the flies (baseline-noise units,
    kcq.apl_effect's) beside its band."""
    from kcq import apl_effect as A
    from kcq import stage_readings as R
    prm = _params(F, run)
    bnd = _bounds(run, reads or {})
    rows = []
    for pan, odors in R.run_panels(run, 'kc_soma'):
        for c in ('on', 'off'):
            if c not in run.conditions:
                continue
            resp, cnt = run.conditions[c]['responses'], run.conditions[c]['spike_counts']
            held = float((resp.to_numpy(dtype=float) > 0.5).mean())
            tgt = np.nan if targets is None else float(targets[0] if c == 'off' else targets[1])
            row = dict(input=panel, panel=pan, condition=c, source=str(getattr(run, 'source', '')),
                       tag=str(getattr(run, 'tag', '')), **prm[c], rate_held=held, frac_target=tgt,
                       hold_miss=abs(held - tgt) / tgt if tgt and np.isfinite(tgt) else np.nan, **bnd,
                       **panel_holds(resp, cnt, odors))
            if pan == 'apl16':
                d = np.array([row[f'opk_{b}'] for b in A.OPK_BINS])
                row.update(opk_tv=float(0.5 * np.abs(d - A.opk_data(c, 'noise')).sum()),
                           opk_tv_band=float(A.opk_band(c, 'noise')))
            rows.append(row)
    return rows


def _reads_by_tag() -> Dict[str, Dict[str, Mapping]]:
    """E139's and E136's reads tables (on the adopted measure), each as {tag: row}."""
    from experiments import E139_normalization_ladder as L
    from experiments import E136_stage_population_fit as F
    out = {}
    for name, tab in (('E139', L.reads_table()), ('E136', F.reads_table())):
        out[name] = {str(r['tag']): r for r in tab.to_dict('records')} if len(tab) else {}
    return out


def holds(variants: Sequence[str] = VARIANTS, inputs: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Step 0's second table: what the holds reach, per variant, input, panel, condition and wiring (`run_holds` on the
    runs E136's score takes, `variant_runs`), the targets E136's `read_targets` (both rates at each panel's operating
    point), with the rows 'mean_S', 'mean_T', 'mean', 'min' and 'max' over the wirings.

    Returns:
        Wide: variant, input, panel, condition, wiring, source, tag, the parameters of that condition (F.VARIANT_PARAMS:
        fixed_thr, wAPLKC_scale, wKCAPL_scale, wAPLPN_scale, wPNAPL_scale, wKCKC_scale), rate_held, frac_target,
        hold_miss, at_bound_thr, at_bound_apl, n_kcs, n_odors, rate, max_count, rate_<type> and n_kcs_<type> (KC_TYPES'
        suffixes: g, apbp for a'b', ab), opk_0, opk_1, opk_2_3, opk_4_7, opk_8plus, and on the APL16 panel opk_tv and
        opk_tv_band.
    """
    rows = []
    with e136(variants) as F:
        reads = _reads_by_tag()
        for v in variants:
            for inp in (inputs or F.READ_PANELS):
                runs = variant_runs(F, v, inp)
                if runs is None:
                    print(f'E141 step 0: {v} has no saved run on {inp} for every wiring; skipped', flush=True)
                    continue
                for w, run in runs.items():
                    n_od = run.conditions['on']['spike_counts'].shape[1]
                    t_off, t_on, f_off, f_on = F.read_targets(w, inp, n_od)
                    for r in run_holds(F, run, inp, (f_off, f_on), reads):
                        rows.append(dict(variant=v, wiring=w, **r))
    df = pd.DataFrame(rows)
    if not len(df):
        return df
    keys = ['variant', 'input', 'panel', 'condition']
    num = [c for c in df.columns
           if c not in keys + ['wiring', 'source', 'tag'] and pd.api.types.is_numeric_dtype(df[c])]
    summ = []
    for k, g in df.groupby(keys, sort=False):
        for lab, ws, fn in (('mean_S', F.S_WIRINGS, 'mean'), ('mean_T', F.T_WIRINGS, 'mean'),
                            ('mean', F.WIRINGS, 'mean'), ('min', F.WIRINGS, 'min'), ('max', F.WIRINGS, 'max')):
            sub = g[g.wiring.isin(ws)]
            if len(sub):
                summ.append(dict(zip(keys, k), wiring=lab, **getattr(sub[num], fn)().to_dict()))
    first = keys + ['wiring', 'source', 'tag']
    out = pd.concat([df, pd.DataFrame(summ)], ignore_index=True)
    return out[first + [c for c in out.columns if c not in first]]


# ---------------------------------------------------------------------------
# A1's evaluations along the KC-KC factor

def kckc_bins(factor: np.ndarray, n_bins: int = N_KCKC_BINS) -> Tuple[np.ndarray, np.ndarray]:
    """The KC-KC factors' bins: `n_bins` equal widths in ln of the factor from 1 to the largest factor (the last bin
    closed above), and one bin from the smallest factor to 1 where any is below 1. Returns (each factor's bin index, 0
    the bin below 1 where there is one; the edges, n_bins + 1 of them, or n_bins + 2 with the bin below 1)."""
    f = np.asarray(factor, dtype=float)
    hi = float(np.nanmax(f))
    if not hi > 1.0:
        raise ValueError(f'the largest KC-KC factor is {hi}, not above 1')
    edges = np.exp(np.linspace(0.0, np.log(hi), n_bins + 1))
    lo = float(np.nanmin(f))
    if lo < 1.0:
        edges = np.r_[lo, edges]
    idx = np.searchsorted(edges, f, side='right') - 1
    return np.clip(idx, 0, len(edges) - 2), edges


def _e136_results() -> Path:
    """E136's results on the adopted measure (read only)."""
    from experiments import E139_normalization_ladder as L
    with L.e136() as F:
        return F.out_dir()


def a1_by_kckc(results: Optional[Path] = None, n_bins: int = N_KCKC_BINS) -> pd.DataFrame:
    """Step 0's third table: A1's saved evaluations (fit_history_A1.csv) and their rescore against both halves
    (rescore_A1.csv, strength by odor; the fit half's equals the history's), binned by the KC-KC factor (`kckc_bins`),
    as E139's `a1_by_rung` bins them by the additive term: per bin and half, 'best' (the bin's evaluation with the
    lowest megamat stage error against the fit half, the S wirings' mean views, read against that half) and 'median'
    (the median over the bin's evaluations of each column), with the bin's edges and number of evaluations, the
    coordinates as factors (pn_inhadd absolute), the megamat stage error, the error over the boutons, claws and somata
    (stage_error_downstream) and each stage's error. Single evaluations whose other factors differ: a lead."""
    from experiments import E139_normalization_ladder as L
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    d = Path(results) if results is not None else _e136_results()
    hist = pd.read_csv(d / f'fit_history_{A1}.csv', float_precision='round_trip')
    resc = pd.read_csv(d / f'rescore_{A1}.csv', float_precision='round_trip')
    resc = resc[resc.measure == L.MEASURE]
    z = hist.set_index('eval')[[f'z_{k}' for k in ('aplpn', 'kckc', 'inhadd', 'inhsc')]]
    by_half = {h: resc[resc.half == h].set_index('eval')[list(STAGE_COLS)] for h in ('fit', 'held')}
    kckc = np.exp(z['z_kckc'])
    idx, edges = kckc_bins(kckc.to_numpy(), n_bins)
    bins = pd.Series(idx, index=z.index)
    inhadd = PN_INHADD_DEFAULT * np.exp(z['z_inhadd'])
    rows = []
    for b in sorted(set(bins)):
        ev = [i for i in bins.index[bins == b] if i in by_half['fit'].index]
        if not ev:
            continue
        best = int(by_half['fit'].loc[ev, 'stage_error'].idxmin())
        for half, tab in by_half.items():
            have = [i for i in ev if i in tab.index]
            base = dict(bin=int(b), bin_lo=float(edges[b]), bin_hi=float(edges[b + 1]),
                        bin_center=float(np.sqrt(edges[b] * edges[b + 1])), half=half, n_evals=len(have))
            for stat in ('best', 'median'):
                if stat == 'best':
                    if best not in tab.index:
                        continue
                    vals, zz = tab.loc[best], z.loc[best]
                    row = dict(base, statistic='best', eval=best, pn_inhadd=float(inhadd.loc[best]))
                else:
                    vals, zz = tab.loc[have].median(), z.loc[have].median()
                    row = dict(base, statistic='median', eval=np.nan, pn_inhadd=float(inhadd.loc[have].median()))
                row.update(kckc_factor=float(np.exp(zz['z_kckc'])), aplpn_factor=float(np.exp(zz['z_aplpn'])),
                           inhsc_gain=float(np.exp(zz['z_inhsc'])))
                row.update({(c.replace('_error', '') if c.endswith('_error') and c != 'stage_error' else c):
                            float(vals[c]) for c in STAGE_COLS})
                rows.append(row)
    return pd.DataFrame(rows)
