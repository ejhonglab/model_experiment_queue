#!/usr/bin/env python3
"""E127 step 0: an expansive output at the PN boutons, read on the family's saved stages (no simulation, no olfsysm
change).

The entry is item E127 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from viki
2026-11-20 21:59; moved in under the idle-host rule 2026-11-22 12:28). "The family" is the breadth-2 APL-effect fit on
E72 part B's eight wirings (S: real, P1, P2, M1; T: P3, P4, M2, M3), each at its own fit.

The change. Each bouton's output becomes r_ref (r / r_ref)^q, r its PN's rate (clipped at 0), q from 1 to 3 (QS),
r_ref the wiring's mean PN rate over the panel (`r_ref_of`: the window mean of E90's saved APL-on PN stage over all
glomeruli and the 16 panel odors), with half and twice that rate beside (REF_FACTORS). The expansion sits before APL's
bouton site (the entry), so:
- APL off (exact for the bouton stage): a bouton is its PN's rate, so its output is f(r) in the window and before the
  odor, the change their difference.
- APL on (first order): the saved APL inhibition at each bouton, I = r - bouton (window mean and pre-odor mean
  separately, E90's saved stages at the settled pre-odor level, E119's `model_stage`), is held and subtracted from f(r),
  clipped at 0. APL's feedback is not solved again, and f acts on window means, not on each step (by Jensen's inequality
  the per-step reading would expand the strong odors a little more).
r_ref is a common factor of f, so with APL off every statistic is independent of it; with APL on it sets the size of f(r)
against the held inhibition and so enters every APL-on number (REF_FACTORS beside).

Readings (E119's, as E120 reads them; `E120.read_units`): the per-odor totals' CV over the 12 megamat odors, the
responders' CV and the top-1 % share of the boutons of the 28 shared glomeruli, unpooled, at E119's criteria (1, 2, 3
baseline-noise SDs), the threshold matched to the flies' mean overall response rate, in E119's three readings ('rate',
'session_z', 'noisy'), against the flies' range (Remy's 5 bouton flies; fly 2022-10-21 fly02 excluded by E119's part 0).
The claws (E119's claw stage, the drive alone, groups of 14 of one glomerulus) beside: each claw's drive changed by the
wiring's PN->claw matrix (the rebuild archive's wPNKC) times its boutons' change, the APL term held, clipped at 0;
deciding nothing.

E98's APL ratios (`apl_ratios`), predicted. On E101's A3 runs (E98's extended input, the eight wirings at their fits,
APL on, stages saved; data/model_outputs/E101_apl_vs_kc_claws/A3/): the boutons expanded as above, each claw's drive
changed by W times its boutons' change (APL's term at the claw held), the summed drive change of the claws of the 36
input glomeruli (E101's 'claw_input_00_change'), and APL's evoked response through E101's line for that wiring (log2 APL
on log2 claws over the 16 panel odors, reduced major axis, results/E101_apl_vs_kc_claws/a3_lines.csv): a ratio is
(x_num / x_den)^slope, the intercept cancelling. Two readings: 'line' (the prediction itself) and 'anchored' (E98's own
linear ratio of that wiring, results/E98_apl_response_strength/c_wirings.csv, times line(q) / line(1)); the anchored
reading is the one classified, the line at q 1 beside as the relation's error. Every carrying comparison of E98 is
computed; V1 (paa over ms near 10^-3, inside if inside under any estimate of ms's input, E98's `v1_verdict`) and V2 (paa
at 10^-5 over 10^-4, validation2/diag) are E98's two that decide; E98's `classify` and 6 of 8 wirings.

The gate (the entry's step 0): open if some q up to 3 brings the totals' CV inside the flies' range on at least 6 of the
8 wirings without taking an amplitude statistic (the responders' CV, the top-1 % share) beyond the flies'. Fixed here,
where the entry gives no rule (flagged for a planner in verdict.json):
- "beyond the flies'" is E119's per-ROI rule: at least 6 of the 8 wirings above the flies' range on that statistic
  (`E119.scalar_differs`); the stricter reading (any wiring above) is beside;
- the gate is read with APL on (the flies' boutons have APL), r_ref at the wiring's mean, E119's verdict criterion
  (2 SDs), in the 'rate' reading (E120's level rule); the other readings, APL off (exact) and r_ref at half and twice
  are beside and decide nothing; whether any of them would open it is recorded;
- the q the boutons ask for: if the gate opens, the smallest q on the grid that passes; otherwise none, and per
  statistic the range of q at which it is inside on 6 of 8 wirings is reported.

Outputs (<shared>/results/E127_bouton_expansion/): step0_bouton_rows.csv (per q, r_ref factor, condition, reading,
criterion, wiring), step0_bouton_summary.csv, step0_claw_summary.csv, step0_apl_ratios.csv (per comparison, wiring, q,
factor), step0_apl_summary.csv, step0_rates.csv (r_ref per wiring, the expanded boutons' largest rates and their share
above 300 Hz), step0_identity.json, verdict.json, comparison_table.md, step0_figure.png, provenance.json.

Usage (from the repository root; light, no simulation):
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E127_bouton_expansion.py --part step0
"""
from __future__ import annotations

import argparse
import json
import time
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E127_bouton_expansion'
QS = tuple(float(x) for x in np.round(np.arange(1.0, 3.0 + 1e-9, 0.125), 3))
REF_FACTORS = (0.5, 1.0, 2.0)
MIN_IN = 6                                  # of the 8 wirings
GATE_READING = 'rate'
HIGH_RATE_HZ = 300.0
AMPLITUDE = ('resp_cv', 'top1')
THREE = ('totals_cv',) + AMPLITUDE
TAG = 'rate0.068_aplboth_k2_cap3'
# wiring label (E90's / E101's) -> (archive experiment folder, rebuild tag), as E121
REBUILD = {'real': ('E67_apl_effect_fits', TAG)}
REBUILD.update({f'P{i}': ('E72_wiring_population', f'wiring-claw_permutation_within_type-seed{i}_{TAG}')
                for i in range(1, 5)})
REBUILD.update({f'M{i}': ('E72_wiring_population', f'wiring-curveball_within_type_matched-seed{i}_{TAG}')
                for i in range(1, 4)})
DECIDING = {'V1': 'V1, ', 'V2': 'V2, '}     # E98's comparisons that decide (prefix of `comparison`, decides True)


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def archive(experiment: str) -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / experiment


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


# ---------------------------------------------------------------------------
# the transform

def expand(r, r_ref: float, q: float):
    """The expansive output r_ref (max(r, 0) / r_ref)^q (q 1: r clipped at 0)."""
    x = np.clip(np.asarray(r, dtype=float), 0, None)
    return r_ref * (x / r_ref) ** q


def expanded_boutons(pn_win: np.ndarray, pn_pre: np.ndarray, bouton_win: np.ndarray, bouton_pre: np.ndarray,
                     r_ref: float, q: float, cond: str) -> Tuple[np.ndarray, np.ndarray]:
    """Boutons x stimuli window and pre-odor outputs after the expansion. `pn_*` are each bouton's PN's rates (boutons x
    stimuli), `bouton_*` the saved bouton stage of the same run. APL off ('off'): f(r) (the saved boutons are their PNs'
    rates; not used). APL on ('on', first order): max(f(r) - I, 0) with I = r - bouton, the saved APL inhibition at the
    bouton, held."""
    out = []
    for r, b in ((pn_win, bouton_win), (pn_pre, bouton_pre)):
        f = expand(r, r_ref, q)
        if cond == 'on':
            f = np.clip(f - (np.asarray(r, dtype=float) - np.asarray(b, dtype=float)), 0, None)
        elif cond != 'off':
            raise ValueError(cond)
        out.append(f)
    return out[0], out[1]


def claws_changed(claw: np.ndarray, W: np.ndarray, d_bouton: np.ndarray) -> np.ndarray:
    """Claws x stimuli drive with its boutons changed by `d_bouton` (boutons x stimuli), APL's term held:
    max(claw + W d_bouton, 0)."""
    return np.clip(np.asarray(claw, dtype=float) + W @ d_bouton, 0, None)


# ---------------------------------------------------------------------------
# the saved stages

@lru_cache(maxsize=None)
def stage(wiring: str, cc: str) -> dict:
    """E90's saved stage run (the panel's 16 odors) at the settled pre-odor level (E119's `model_stage`)."""
    from experiments import E119_stages_per_roi as E119
    return E119.model_stage(wiring, cc)


def r_ref_of(wiring: str) -> float:
    """The wiring's mean PN rate over the panel (E90's saved APL-on PN stage, window mean, all glomeruli, 16 odors)."""
    return float(stage(wiring, 'on')['pn_win'][list(H.APL16_ODORS)].to_numpy(dtype=float).mean())


@lru_cache(maxsize=None)
def wiring_matrix(wiring: str) -> pd.DataFrame:
    """The wiring's claw x bouton matrix from its rebuild archive (wPNKC.parquet)."""
    exp, tag = REBUILD[wiring]
    return pd.read_parquet(archive(exp) / 'rebuild' / tag / 'wPNKC.parquet')


def _aligned_W(wiring: str, claws: pd.Index, boutons: pd.Index) -> np.ndarray:
    W = wiring_matrix(wiring)
    if not (W.index.equals(claws) and W.columns.equals(boutons)):
        W = W.reindex(index=claws, columns=boutons)
        if W.isna().any().any():
            raise RuntimeError(f'{wiring}: the rebuild wiring does not hold the saved claws and boutons')
    return W.to_numpy(dtype=float)


def _pn_of_boutons(pn: pd.DataFrame, gl: np.ndarray) -> np.ndarray:
    missing = sorted(set(gl) - set(pn.index))
    if missing:
        raise RuntimeError(f'boutons of glomeruli without a PN row: {missing[:5]}')
    return pn.loc[gl].to_numpy(dtype=float)


def identity(wirings: Sequence[str]) -> dict:
    """Checks on the saved stages: the PN stage equal with APL on and off; with APL off the boutons equal their PNs'
    rates and the claws equal W times the boutons (the rebuild wiring is the run's); per wiring the largest differences."""
    from experiments import E119_stages_per_roi as E119
    out = {}
    for w in wirings:
        on, off = stage(w, 'on'), stage(w, 'acute_off')
        gl = np.asarray(E119.unit_meta(w)['bouton_glomerulus'])
        Wm = _aligned_W(w, off['claw_win'].index, off['bouton_win'].index)
        out[w] = dict(pn_on_minus_off=float(np.abs(on['pn_win'].to_numpy() - off['pn_win'].to_numpy()).max()),
                      bouton_off_minus_pn=float(np.abs(off['bouton_win'].to_numpy() - _pn_of_boutons(off['pn_win'], gl)).max()),
                      claw_off_minus_Wb=float(np.abs(off['claw_win'].to_numpy() - Wm @ off['bouton_win'].to_numpy()).max()),
                      max_claw=float(off['claw_win'].to_numpy().max()))
        out[w]['exact'] = bool(out[w]['pn_on_minus_off'] == 0 and out[w]['bouton_off_minus_pn'] < 1e-9
                               and out[w]['claw_off_minus_Wb'] < 1e-6 * max(1.0, out[w]['max_claw']))
    return out


# ---------------------------------------------------------------------------
# the bouton (and claw) statistics

def panel_units(wiring: str, q: float, factor: float, cond: str, with_claws: bool = False) -> dict:
    """One wiring's expanded boutons (and claws) on the panel as units x odors changes (odor columns abbreviated), with
    the expanded boutons' window rates (for the realism counts)."""
    from experiments import E119_stages_per_roi as E119
    st = stage(wiring, 'on' if cond == 'on' else 'acute_off')
    gl = np.asarray(E119.unit_meta(wiring)['bouton_glomerulus'])
    cols = list(H.APL16_ODORS)
    r_ref = factor * r_ref_of(wiring)
    bw, bp = st['bouton_win'][cols].to_numpy(dtype=float), st['bouton_pre'][cols].to_numpy(dtype=float)
    ew, ep = expanded_boutons(_pn_of_boutons(st['pn_win'][cols], gl), _pn_of_boutons(st['pn_pre'][cols], gl), bw, bp,
                              r_ref, q, cond)
    ab = [H.odor_abbrev(str(c)) for c in cols]
    out = dict(bouton=pd.DataFrame(ew - ep, columns=ab), bouton_win=ew, r_ref=r_ref)
    if with_claws:
        Wm = _aligned_W(wiring, st['claw_win'].index, st['bouton_win'].index)
        cw = claws_changed(st['claw_win'][cols].to_numpy(dtype=float), Wm, ew - bw)
        cp = claws_changed(st['claw_pre'][cols].to_numpy(dtype=float), Wm, ep - bp)
        out['claw'] = pd.DataFrame(cw - cp, columns=ab)
    return out


def bouton_stats(qs: Sequence[float] = QS, wirings: Optional[Sequence[str]] = None, factors: Sequence[float] = REF_FACTORS,
                 readings: Optional[Sequence[str]] = None, claws: bool = True) -> Dict[str, pd.DataFrame]:
    """E119's per-ROI statistics of the expanded boutons (and the claws beside) per q, r_ref factor, condition (APL on
    at every factor, APL off at factor 1, its numbers not depending on it), reading, criterion and wiring; and the
    expanded boutons' rates. Returns {'bouton', 'claw', 'rates'}."""
    from experiments import E120_bouton_stage_measured as E120
    wirings = tuple(wirings or E120.WIRINGS)
    readings = tuple(readings or E120.READINGS)
    rows_b, rows_c, rates = [], [], []
    for q in qs:
        for cond, fs in (('on', factors), ('off', (1.0,))):
            for fac in fs:
                with_c = claws and cond == 'on' and fac == 1.0
                pu = {w: panel_units(w, q, fac, cond, with_c) for w in wirings}
                keys = dict(q=float(q), factor=float(fac), cond=cond)
                ub = {w: E120.stage_subset(pu[w]['bouton'], 'bouton', w) for w in wirings}
                rows_b.append(E120.read_units(ub, 'bouton', readings, **keys))
                if with_c:
                    uc = {w: E120.stage_subset(pu[w]['claw'], 'claw', w) for w in wirings}
                    rows_c.append(E120.read_units(uc, 'claw', ('rate',), **keys))
                for w in wirings:
                    x = pu[w]['bouton_win']
                    rates.append(dict(**keys, wiring=w, r_ref=pu[w]['r_ref'], max_rate=float(x.max()),
                                      share_above_300=float((x > HIGH_RATE_HZ).mean())))
        print(f'q {q:g} done', flush=True)
    return dict(bouton=pd.concat(rows_b, ignore_index=True),
                claw=pd.concat(rows_c, ignore_index=True) if rows_c else pd.DataFrame(), rates=pd.DataFrame(rates))


def summarize(rows: pd.DataFrame, stats: Sequence[str] = THREE) -> pd.DataFrame:
    """Per q, factor, condition, reading and criterion: each statistic's mean, min and max over the wirings and the
    numbers of wirings below / inside / above the flies' range (`E120.summarize`)."""
    from experiments import E120_bouton_stage_measured as E120
    return E120.summarize(rows, ('q', 'factor', 'cond', 'reading', 'criterion'), stats)


# ---------------------------------------------------------------------------
# E98's APL ratios through E101's line

@lru_cache(maxsize=None)
def a3_run(wiring: str) -> dict:
    from experiments import E101_apl_vs_kc_claws as E101
    run = E101.load_a3_run(wiring, archive('E101_apl_vs_kc_claws') / 'A3')
    if run is None:
        raise RuntimeError(f"E101's A3 run of {wiring} is not in the archive")
    return run


def claw_input_change(wiring: str, q: float, factor: float) -> pd.Series:
    """Per stimulus of E101's A3 run (APL on), the summed drive change of the claws of the 36 input glomeruli with the
    boutons expanded (APL's terms at the bouton and the claw held)."""
    run = a3_run(wiring)
    on, meta = run['on'], run['meta']
    stims = list(on['claw_win'].columns)
    gl = np.asarray(meta['bouton_glomerulus']).astype(str)
    pn_w, pn_p = on['pn_win'].copy(), on['pn_pre'].copy()
    pn_w.index, pn_p.index = pn_w.index.astype(str), pn_p.index.astype(str)
    bw, bp = on['bouton_win'][stims].to_numpy(dtype=float), on['bouton_pre'][stims].to_numpy(dtype=float)
    ew, ep = expanded_boutons(_pn_of_boutons(pn_w[stims], gl), _pn_of_boutons(pn_p[stims], gl), bw, bp,
                              factor * r_ref_of(wiring), q, 'on')
    Wm = _aligned_W(wiring, on['claw_win'].index, on['bouton_win'].index)
    m = np.isin(np.asarray(meta['claw_glomerulus']).astype(str), [str(g) for g in H.apl16_orn_deltas().index])
    cw = claws_changed(on['claw_win'][stims].to_numpy(dtype=float)[m], Wm[m], ew - bw)
    cp = claws_changed(on['claw_pre'][stims].to_numpy(dtype=float)[m], Wm[m], ep - bp)
    return pd.Series((cw - cp).sum(axis=0), index=stims)


def e98_tables() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """E98's carrying comparisons (c_comparisons.csv) and their per-wiring linear ratios (c_wirings.csv)."""
    from experiments import E98_apl_response_strength as E98
    d = E98.out_dir() if hasattr(E98, 'out_dir') else H.results_dir('E98_apl_response_strength')
    c = pd.read_csv(Path(d) / 'c_comparisons.csv')
    c = c[c.reading_set == 'carrying'].reset_index(drop=True)
    w = pd.read_csv(Path(d) / 'c_wirings.csv')
    w = w[w.reading_set == 'carrying']
    return c, w


def e101_lines() -> pd.DataFrame:
    from experiments import E101_apl_vs_kc_claws as E101
    l = pd.read_csv(E101.out_dir() / 'a3_lines.csv')
    return l.groupby('wiring', sort=False)[['slope_rma', 'intercept_rma', 'r', 'n_odors']].first()


def apl_ratios(qs: Sequence[float] = QS, wirings: Optional[Sequence[str]] = None,
               factors: Sequence[float] = REF_FACTORS) -> pd.DataFrame:
    """E98's carrying ratios predicted per comparison, wiring, q and r_ref factor: 'line' (x_num / x_den)^slope, the
    claws' summed drive change through E101's line; 'anchored' E98's linear ratio times line(q) / line(1); the class of
    the anchored ratio against the flies' interval (E98's `classify`); E101's identity at q 1 (the claw total against
    a3_family.csv's claw_input_00_change)."""
    from experiments import E98_apl_response_strength as E98
    from experiments import E101_apl_vs_kc_claws as E101
    from experiments import E120_bouton_stage_measured as E120
    wirings = tuple(wirings or E120.WIRINGS)
    comps, cw = e98_tables()
    lines = e101_lines()
    fam = pd.read_csv(E101.out_dir() / 'a3_family.csv')
    fam = fam[fam.cond == 'on'].set_index(['wiring', 'stimulus'])['claw_input_00_change']
    rows = []
    for w in wirings:
        slope = float(lines.loc[w, 'slope_rma'])
        x1 = {fac: claw_input_change(w, 1.0, fac) for fac in factors}
        ident = float(np.max(np.abs(x1[1.0].to_numpy() - fam.loc[w].reindex(x1[1.0].index).to_numpy())))
        for q in qs:
            for fac in factors:
                x = x1[fac] if q == 1.0 else claw_input_change(w, q, fac)
                for _, c in comps.iterrows():
                    num, den = str(c.numerator), str(c.denominator)
                    if num not in x.index or den not in x.index:
                        raise RuntimeError(f'{w}: {num} or {den} not in A3 run')
                    line = (x[num] / x[den]) ** slope if x[num] > 0 and x[den] > 0 else np.nan
                    line1 = (x1[fac][num] / x1[fac][den]) ** slope
                    e = cw[(cw.comparison == c.comparison) & (cw.wiring == w)]
                    e98 = float(e.ratio_linear.iloc[0]) if len(e) else np.nan
                    anch = e98 * line / line1
                    cls = (E98.classify([anch], c.flies_lo90, c.flies_hi90) if np.isfinite(anch) and anch > 0
                           else 'undefined')
                    rows.append(dict(comparison=c.comparison, decides=bool(c.decides), numerator=num, denominator=den,
                                     wiring=w, q=float(q), factor=float(fac), slope=slope, x_num=float(x[num]),
                                     x_den=float(x[den]), line=float(line), e98_linear=e98, anchored=float(anch),
                                     flies_lo90=float(c.flies_lo90), flies_hi90=float(c.flies_hi90), cls=cls,
                                     claw_identity_q1=ident))
        print(f'APL ratios {w} done', flush=True)
    return pd.DataFrame(rows)


def apl_summary(ar: pd.DataFrame) -> pd.DataFrame:
    """Per comparison, q and factor: the anchored and line ratios' range over the wirings, the class counts and E98's
    comparison verdict (6 of 8)."""
    from experiments import E98_apl_response_strength as E98
    out = []
    for (comp, q, fac), g in ar.groupby(['comparison', 'q', 'factor'], sort=False):
        cl = list(g.cls)
        out.append(dict(comparison=comp, decides=bool(g.decides.iloc[0]), q=q, factor=fac,
                        flies_lo90=g.flies_lo90.iloc[0], flies_hi90=g.flies_hi90.iloc[0],
                        anchored_min=g.anchored.min(), anchored_max=g.anchored.max(), line_min=g.line.min(),
                        line_max=g.line.max(), n_inside=cl.count('inside'), n_flatter=cl.count('flatter'),
                        n_steeper=cl.count('steeper'), n_mixed=cl.count('mixed'),
                        category=E98.comparison_verdict(cl, MIN_IN)))
    return pd.DataFrame(out)


def apl_verdicts(asumm: pd.DataFrame) -> pd.DataFrame:
    """Per q and factor, E98's V1 (`v1_verdict` over the estimates of ms's input) and V2 (validation2/diag)."""
    from experiments import E98_apl_response_strength as E98
    d = asumm[asumm.decides]
    out = []
    for (q, fac), g in d.groupby(['q', 'factor'], sort=False):
        v1 = {c: k for c, k in zip(g.comparison, g.category) if c.startswith(DECIDING['V1'])}
        v2 = [k for c, k in zip(g.comparison, g.category) if c.startswith(DECIDING['V2'])]
        out.append(dict(q=q, factor=fac, V1=E98.v1_verdict(v1), V2=v2[0] if v2 else None))
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# the gate

def gate_rows(summ: pd.DataFrame, criterion: float, reading: str, cond: str, factor: float,
              min_in: int = MIN_IN) -> pd.DataFrame:
    """Per q under one reading: the totals' CV inside on `min_in` wirings, each amplitude statistic beyond the flies'
    (E119's rule: `min_in` wirings above; and the strict reading, any wiring above), and whether q passes."""
    s = summ[(summ.criterion == criterion) & (summ.reading == reading) & (summ.cond == cond)
             & np.isclose(summ.factor, factor)].sort_values('q')
    out = s[['q', 'totals_cv', 'totals_cv_in', 'resp_cv', 'resp_cv_in', 'resp_cv_above', 'top1', 'top1_in',
             'top1_above']].copy()
    out['totals_in'] = out.totals_cv_in >= min_in
    out['amp_beyond'] = np.any([out[f'{a}_above'] >= min_in for a in AMPLITUDE], axis=0)
    out['amp_beyond_strict'] = np.any([out[f'{a}_above'] >= 1 for a in AMPLITUDE], axis=0)
    out['passes'] = out.totals_in & ~out.amp_beyond
    out['passes_strict'] = out.totals_in & ~out.amp_beyond_strict
    return out.reset_index(drop=True)


def q_range(g: pd.DataFrame, col: str, min_in: int = MIN_IN) -> Optional[List[float]]:
    """The smallest and largest q at which `col` (a count of wirings inside) reaches `min_in`; None if none."""
    q = g.q[g[col] >= min_in]
    return [float(q.min()), float(q.max())] if len(q) else None


def gate(summ: pd.DataFrame, criterion: Optional[float] = None) -> Dict:
    """The step-0 gate (module docstring): the primary reading's verdict, the q the boutons ask for, and every other
    reading's (decides nothing)."""
    from experiments import E120_bouton_stage_measured as E120
    criterion = E120.VERDICT_CRITERION if criterion is None else criterion
    prim = gate_rows(summ, criterion, GATE_READING, 'on', 1.0)
    passing = prim.q[prim.passes]
    beside = {}
    for (rd, cond, fac), _ in summ[summ.criterion == criterion].groupby(['reading', 'cond', 'factor'], sort=False):
        g = gate_rows(summ, criterion, rd, cond, fac)
        beside[f'{rd}|{cond}|r_ref x{fac:g}'] = dict(
            opens=bool(g.passes.any()), opens_strict=bool(g.passes_strict.any()),
            q_pass=[float(x) for x in g.q[g.passes]],
            q_in=dict(totals_cv=q_range(g, 'totals_cv_in'), resp_cv=q_range(g, 'resp_cv_in'),
                      top1=q_range(g, 'top1_in')),
            totals_cv_at_q_max=float(g.totals_cv.iloc[-1]), max_totals_cv_in=int(g.totals_cv_in.max()))
    return dict(criterion=criterion, reading=GATE_READING, cond='on', factor=1.0, opens=bool(len(passing)),
                q=float(passing.min()) if len(passing) else None, opens_strict=bool(prim.passes_strict.any()),
                q_in=dict(totals_cv=q_range(prim, 'totals_cv_in'), resp_cv=q_range(prim, 'resp_cv_in'),
                          top1=q_range(prim, 'top1_in')),
                any_other_reading_opens=any(v['opens'] for v in beside.values()),
                per_q=prim.to_dict('records'), beside=beside)


# ---------------------------------------------------------------------------
# the part

def figure(summ: pd.DataFrame, asumm: pd.DataFrame, path: Path, criterion: float) -> None:
    """The three bouton statistics (APL on, r_ref at the mean, 'rate'; APL off beside) and the anchored V2 / V1 ratios
    against q, each against the flies' range."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from experiments import E120_bouton_stage_measured as E120
    sp = E120.fly_spreads('bouton')[criterion]
    fig, axs = plt.subplots(1, 5, figsize=(19, 3.8))
    names = dict(totals_cv="per-odor totals' CV", resp_cv="responders' CV", top1='top-1 % share')
    for ax, s in zip(axs, THREE):
        for cond, ls in (('on', '-'), ('off', '--')):
            g = summ[(summ.criterion == criterion) & (summ.reading == GATE_READING) & (summ.cond == cond)
                     & np.isclose(summ.factor, 1.0)].sort_values('q')
            ax.plot(g.q, g[s], ls, color='C0', label=f'model, APL {cond} (mean of 8)')
            ax.fill_between(g.q, g[f'{s}_lo_w'], g[f'{s}_hi_w'], color='C0', alpha=0.15 if cond == 'on' else 0.07)
        ax.axhspan(sp[s]['lo'], sp[s]['hi'], color='k', alpha=0.15, label="flies' range")
        ax.set_title(f'boutons: {names[s]}')
        ax.set_xlabel('exponent q')
    axs[0].legend(fontsize=7)
    for ax, key in zip(axs[3:], ('V2, ', 'V1, ms\'s input from megamat/megamat')):
        g = asumm[asumm.comparison.str.startswith(key) & np.isclose(asumm.factor, 1.0)].sort_values('q')
        if not len(g):
            continue
        ax.fill_between(g.q, g.anchored_min, g.anchored_max, color='C1', alpha=0.4, label='model (8 wirings)')
        for fac, ls in ((0.5, ':'), (2.0, '--')):
            h = asumm[(asumm.comparison == g.comparison.iloc[0]) & np.isclose(asumm.factor, fac)].sort_values('q')
            ax.plot(h.q, h.anchored_min, ls, color='C1', lw=0.8, label=f'r_ref x{fac:g} (min, max)')
            ax.plot(h.q, h.anchored_max, ls, color='C1', lw=0.8)
        ax.axhspan(g.flies_lo90.iloc[0], g.flies_hi90.iloc[0], color='k', alpha=0.15, label="flies' 90 % interval")
        ax.set_title(("APL: paa at 10^-5 over 10^-4" if key.startswith('V2') else 'APL: paa over ms near 10^-3')
                     + '\n(predicted, anchored)')
        ax.set_xlabel('exponent q')
    axs[3].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def comparison_table(verdict: dict, summ: pd.DataFrame, avs: pd.DataFrame, asumm: pd.DataFrame,
                     rates: pd.DataFrame, criterion: float) -> str:
    from experiments import E120_bouton_stage_measured as E120
    sp = E120.fly_spreads('bouton')[criterion]
    show = [1.0, 1.5, 2.0, 2.5, 3.0]
    L = ['# E127 step 0: an expansive output at the PN boutons, read on the family\'s saved stages', '',
         f"Gate: **{'open' if verdict['opens'] else 'closed'}** (APL on, first order; r_ref the wiring's mean PN rate; "
         f"'{GATE_READING}' reading, criterion {criterion:g} SDs; 6 of 8 wirings). q the boutons ask for: "
         f"{verdict['q']}. Another reading opens it: {verdict['any_other_reading_opens']}.", '',
         '## The bouton statistics (mean over the 8 wirings; wirings inside / above the flies\' range)', '',
         "| statistic | flies' range | " + ' | '.join(f'q {q:g}' for q in show) + ' |',
         '|---|---|' + '---|' * len(show)]
    g = summ[(summ.criterion == criterion) & (summ.reading == GATE_READING) & (summ.cond == 'on')
             & np.isclose(summ.factor, 1.0)].set_index('q')
    for s, nm in (('totals_cv', "per-odor totals' CV"), ('resp_cv', "responders' CV"), ('top1', 'top-1 % share')):
        L.append(f"| {nm}, APL on | {sp[s]['lo']:.3f} to {sp[s]['hi']:.3f} | " + ' | '.join(
            f"{g.loc[q, s]:.3f} ({int(g.loc[q, f'{s}_in'])} / {int(g.loc[q, f'{s}_above'])})" for q in show) + ' |')
    g = summ[(summ.criterion == criterion) & (summ.reading == GATE_READING) & (summ.cond == 'off')].set_index('q')
    L.append(f"| per-odor totals' CV, APL off (exact) | {sp['totals_cv']['lo']:.3f} to {sp['totals_cv']['hi']:.3f} | "
             + ' | '.join(f"{g.loc[q, 'totals_cv']:.3f} ({int(g.loc[q, 'totals_cv_in'])} / "
                          f"{int(g.loc[q, 'totals_cv_above'])})" for q in show) + ' |')
    L += ['', "## E98's APL ratios, predicted (anchored; range over the 8 wirings; r_ref at the mean)", '',
          "| comparison | flies' 90 % interval | " + ' | '.join(f'q {q:g}' for q in show) + ' |',
          '|---|---|' + '---|' * len(show)]
    a = asumm[asumm.decides & np.isclose(asumm.factor, 1.0)]
    for comp, h in a.groupby('comparison', sort=False):
        h = h.set_index('q')
        L.append(f"| {comp} | {h.flies_lo90.iloc[0]:.2f} to {h.flies_hi90.iloc[0]:.2f} | " + ' | '.join(
            f"{h.loc[q, 'anchored_min']:.2f} to {h.loc[q, 'anchored_max']:.2f} ({h.loc[q, 'category']})" for q in show)
            + ' |')
    L += ['', '| q | r_ref factor | V1 | V2 |', '|---|---|---|---|']
    for _, r in avs[avs.q.isin(show)].iterrows():
        L.append(f'| {r.q:g} | {r.factor:g} | {r.V1} | {r.V2} |')
    rr = rates[(rates.cond == 'on') & np.isclose(rates.factor, 1.0)].groupby('q')
    L += ['', '## Realism: the expanded boutons\' window rates (APL on, r_ref at the mean; panel odors)', '',
          '| q | largest rate (Hz, over wirings) | share of bouton x odor above 300 Hz (max over wirings) |', '|---|---|---|']
    for q, h in rr:
        if q in show:
            L.append(f'| {q:g} | {h.max_rate.max():.0f} | {h.share_above_300.max():.4f} |')
    L += ['', 'r_ref per wiring (Hz): ' + ', '.join(f'{w} {v:.2f}' for w, v in
                                                     rates[np.isclose(rates.factor, 1.0)].drop_duplicates('wiring').set_index('wiring').r_ref.items())]
    return '\n'.join(L) + '\n'


def part_step0(qs: Sequence[float] = QS, wirings: Optional[Sequence[str]] = None) -> Dict:
    """Step 0 whole: the identities, the bouton and claw statistics, the APL ratios, the gate; every output file."""
    from experiments import E120_bouton_stage_measured as E120
    t0 = time.time()
    wirings = tuple(wirings or E120.WIRINGS)
    od = out_dir()
    missing = []
    for w in wirings:
        try:
            stage(w, 'on'), stage(w, 'acute_off'), a3_run(w), wiring_matrix(w)
        except Exception as e:                                       # a saved stage is missing: named, gate on the rest
            missing.append(dict(wiring=w, error=repr(e)))
    use = tuple(w for w in wirings if w not in {m['wiring'] for m in missing})
    ident = identity(use)
    _write_json(od / 'step0_identity.json', dict(stages=ident, missing=missing))
    if not all(v['exact'] for v in ident.values()):
        raise RuntimeError(f'identity failed: {ident}')
    bs = bouton_stats(qs, use)
    bs['bouton'].to_csv(od / 'step0_bouton_rows.csv', index=False)
    summ = summarize(bs['bouton'])
    summ.to_csv(od / 'step0_bouton_summary.csv', index=False)
    if len(bs['claw']):
        summarize(bs['claw']).to_csv(od / 'step0_claw_summary.csv', index=False)
    bs['rates'].to_csv(od / 'step0_rates.csv', index=False)
    ar = apl_ratios(qs, use)
    ar.to_csv(od / 'step0_apl_ratios.csv', index=False)
    asumm = apl_summary(ar)
    asumm.to_csv(od / 'step0_apl_summary.csv', index=False)
    avs = apl_verdicts(asumm)
    avs.to_csv(od / 'step0_apl_verdicts.csv', index=False)
    v = gate(summ)
    v1 = avs.set_index(['q', 'factor'])
    v.update(n_wirings=len(use), missing=missing,
             apl_at_q=None if v['q'] is None else dict(v1.loc[(v['q'], 1.0)]),
             apl_q_inside=dict(V2=[float(q) for q in avs.q[(avs.factor == 1.0) & (avs.V2 == 'inside')]],
                               V1=[float(q) for q in avs.q[(avs.factor == 1.0) & (avs.V1 == 'inside')]]),
             rules_chosen=[
                 "beyond the flies' = at least 6 of 8 wirings above the flies' range (E119's per-ROI rule); "
                 "the strict reading (any wiring above) beside",
                 "gate read with APL on (first order), r_ref at the wiring's mean, 'rate' reading, criterion 2 SDs; "
                 "other readings, APL off and r_ref at half / twice beside",
                 'the q the boutons ask for: the smallest q on the grid that passes',
                 "APL ratios: classified on the anchored reading (E98's linear ratio x line(q) / line(1))"],
             flag_for_planner='the entry gives no rule for "beyond the flies\'" nor for choosing q among passing values; '
                              'the rules above were fixed by the code subagent',
             outcome='(d) step 0 fails: no simulation is spent' if not v['opens'] else 'step 0 passes: steps 1-2 at q')
    _write_json(od / 'verdict.json', v)
    (od / 'comparison_table.md').write_text(comparison_table(v, summ, avs, asumm, bs['rates'], v['criterion']))
    figure(summ, asumm, od / 'step0_figure.png', v['criterion'])
    _write_json(od / 'provenance.json', dict(part='step0', script=H.script_path(__file__), **H.code_commits(),
                generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'), seconds=round(time.time() - t0, 1),
                settings=dict(QS=list(qs), REF_FACTORS=REF_FACTORS, MIN_IN=MIN_IN, GATE_READING=GATE_READING,
                              wirings=list(use), REBUILD=REBUILD,
                              sources=['data/model_outputs/E90_model_stages/B', 'results/E90_model_stages/settled_preodor',
                                       'data/model_outputs/E101_apl_vs_kc_claws/A3', 'results/E101_apl_vs_kc_claws/a3_lines.csv',
                                       'results/E101_apl_vs_kc_claws/a3_family.csv',
                                       'results/E98_apl_response_strength/c_comparisons.csv, c_wirings.csv',
                                       'data/model_outputs/<exp>/rebuild/<tag>/wPNKC.parquet',
                                       "E119's per-ROI flies (part 0)"])))
    print(f"gate {'open' if v['opens'] else 'closed'}; q {v['q']}; {time.time() - t0:.0f} s", flush=True)
    return v


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', choices=['step0'], required=True)
    ap.add_argument('--wirings', default=None, help='comma-separated (default the eight)')
    a = ap.parse_args()
    part_step0(wirings=a.wirings.split(',') if a.wirings else None)


if __name__ == '__main__':
    main()
