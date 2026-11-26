"""The tables of E139's write-up (the weaker antennal-lobe normalization without the calyx refit, along its strength),
from the files its parts wrote to results/E139_normalization_ladder/ (step 0's step0_*.csv and .json; the score's
score_*.csv, classes.json and provenance_score.json; the reads' reads_<host>_<lane>.csv), from E134's flies tables (the
floors and the no-step yardstick) and, for the normalization's share of the PN gain's denominator, from E135's port of
the family's PN stage (seconds; `--no-port` skips it). No simulation of the mushroom body and no fit.

The megamat, held-out-panel, without-aa and realism tables are E136's write-up tables (analyses/E136_writeup_tables.py)
on E139's files, with E139's variants in place of E136's arms (`as_e136`); the step-0 tables, the classes with the price
on the APL effect, the APL effect by target, the beside step and the reads by lane are E139's own.

Writes writeup_tables.md and writeup_values.csv (the same numbers, long: arm, stage, panel, condition, half, statistic,
value, 'arm' being the variant, for the summary's placeholders) beside the score's files.

Usage (from the repository root): ./venv/bin/python -m analyses.E139_writeup_tables [--score-dir <dir>] [--no-port]
"""
from __future__ import annotations

import argparse
import contextlib
import json
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from analyses import E136_writeup_tables as W
from experiments import E136_stage_population_fit as F
from experiments import E139_normalization_ladder as L

STAGES = W.STAGES
STAGE_NAMES = W.STAGE_NAMES
SETTING_NAMES = {'F': "the family's calyx", 'B': 'APL→bouton ×0.155', 'K': 'KC–KC ×3.24'}
TARGET_NAMES = {'z8': 'all recordings', 'halfA': 'half A', 'halfB': 'half B'}
_f, _mr, _rng, _yn = W._f, W._mr, W._rng, W._yn


def variant_order(names: Sequence[str]) -> List[str]:
    """The family, then F, B and K by rung, step 3's arms, then A1."""
    kinds = ('F', 'B', 'K', L.FIT_PREFIX)

    def key(v):
        s = L._setting_of(v)
        return kinds.index(s) if s in kinds else len(kinds), L._rung_of(v)
    rest = sorted((v for v in names if v not in (F.FAMILY, L.A1)), key=key)
    return [F.FAMILY] + rest + ([L.A1] if L.A1 in names else [])


def variant_label(v: str) -> str:
    """A variant's name in the tables."""
    if v == F.FAMILY:
        return 'the family (unrefitted, pn_inhadd 31.4)'
    if v == L.A1:
        return 'A1 (E136): pn_inhadd 1017, APL→bouton ×0.155, KC–KC ×3.24'
    s, r = L._setting_of(v), L._rung_of(v)
    if s == L.FIT_PREFIX:
        return f'A0 refitted at pn_inhadd {r:g}'
    return f'{SETTING_NAMES.get(s, s)}, pn_inhadd {r:g}'


@contextlib.contextmanager
def as_e136(order: Sequence[str]) -> Iterator[None]:
    """For the block, E136's write-up tables iterate `order` (E139's variants) and name them by `variant_label` (their
    module globals ARM_ORDER and ARM_NAMES, restored after)."""
    old = W.ARM_ORDER, W.ARM_NAMES
    W.ARM_ORDER = tuple(order)
    W.ARM_NAMES = dict(old[1], **{v: variant_label(v) for v in order})
    try:
        yield
    finally:
        W.ARM_ORDER, W.ARM_NAMES = old


# ---------------------------------------------------------------------------
# inputs

def _csv(p: Path) -> pd.DataFrame:
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def load(score_dir: Path, results_dir: Optional[Path] = None) -> dict:
    """Everything the tables read, in the form E136's tables take ('arm' the variant), with E139's own files beside:
    the score's from `score_dir`, step 0's and the reads' from `results_dir` (default E139's out_dir)."""
    results_dir = L.out_dir() if results_dir is None else Path(results_dir)
    def arm_col(df):
        if len(df) and 'variant' in df.columns:
            df = df.drop(columns=[c for c in ('arm',) if c in df.columns]).rename(columns={'variant': 'arm'})
        return df
    d = dict(classes=W._read_json(score_dir / 'classes.json'),
             mega=arm_col(_csv(score_dir / 'score_megamat.csv')),
             panels=_csv(score_dir / 'score_panels.csv'),
             apl=_csv(score_dir / 'score_apl16.csv'),
             without=_csv(score_dir / f'score_without_{L.WITHOUT_ODOR}.csv'),
             extended=_csv(score_dir / 'score_beside_extended.csv'),
             extended_w=_csv(score_dir / 'score_beside_extended_wirings.csv'),
             prov=W._read_json(score_dir / 'provenance_score.json'), fits={}, hist={}, dets={})
    for name in ('a1_by_rung', 'e102_ladder', 'pn_stage_ladder', 'recorded_step_signature'):
        d[f'step0_{name}'] = _csv(results_dir / f'step0_{name}.csv')
    for name in ('pn_stage_ladder', 'recorded_step_signature'):
        d[f'step0_{name}_json'] = W._read_json(results_dir / f'step0_{name}.json')
    reads = [pd.read_csv(p) for p in sorted(results_dir.glob('reads_*.csv')) if 'smoke' not in p.name]
    d['reads'] = pd.concat(reads, ignore_index=True) if reads else pd.DataFrame()
    e = F.e134_dir()
    d['flies_views'] = _csv(e / 'flies_views.csv')
    d['no_step'] = _csv(e / 'flies_no_step.csv')
    return d


def present(d: dict) -> List[str]:
    return variant_order(list(dict.fromkeys(d['mega'].arm))) if len(d['mega']) else []


# ---------------------------------------------------------------------------
# step 0

def a1_by_rung_table(d: dict, V: W.Values) -> List[str]:
    """A1's saved evaluations binned by rung (step0_a1_by_rung.csv): per bin and half, the best evaluation's and the
    median evaluation's stage errors."""
    a = d['step0_a1_by_rung']
    if not len(a):
        return []
    L_ = ["A1's saved evaluations on the S wirings, binned by rung of the additive term (the evaluation's pn_inhadd "
          "nearest the rung on a log scale; single evaluations whose other factors differ, so a lead): the best "
          "evaluation of each bin on the fit half's loss and the bin's medians, against the fit half and the held-out "
          "half (strength by odor, one-fly units).", '',
          '| rung (bin) | half | evaluations | statistic | pn_inhadd | APL→bouton | KC–KC | megamat stage error | '
          + ' | '.join(STAGE_NAMES[s] for s in STAGES) + ' |', '|---:|---|---:|---|---:|---:|---:|---:|'
          + '---:|' * len(STAGES)]
    for r in a.itertuples():
        L_.append(f'| {r.bin_rung:g} | {r.half} | {int(r.n_evals)} | {r.statistic} | {_f(r.pn_inhadd, ".0f")} | '
                  f'{_f(r.aplpn_factor)} | {_f(r.kckc_factor, ".2f")} | {_f(r.stage_error)} | '
                  + ' | '.join(_f(getattr(r, s)) for s in STAGES) + ' |')
        for col in ('stage_error', 'pn_inhadd', *STAGES):
            V.add('A1_evals', f'rung{r.bin_rung:g}', f'{r.statistic}_{col}', getattr(r, col), half=f'{r.half}_S')
    return L_ + ['']


def pn_ladder_table(d: dict, V: W.Values) -> List[str]:
    """The PN stage alone along the ladder (step0_pn_stage_ladder.csv; E135's port): against the held-out half of the
    PN-dendrite flies, the mean over the eight wirings' input sampling, the strength error, every view and the spread
    of the odors' totals beside the flies'."""
    p = d['step0_pn_stage_ladder']
    if not len(p):
        return []
    x = p[(p.half == 'held') & (p.row == 'mean') & (p.views == 'all')]
    L_ = ["The PN stage alone along the ladder (E135's port of the family's PN stage, the calyx not run), against the "
          "held-out half of the PN-dendrite flies, mean over the eight wirings: each view's error in one-fly units "
          "(lower is better) and the spread of the odors' totals (the coefficient of variation over the 12 megamat odors "
          "of each odor's mean evoked rate; the flies' "
          f"{_f(L.FLIES_TOTALS_CV)}).", '',
          '| point | pn_inhadd | pn_inhsc | strength by odor | sharing | breadth | the matrix itself | all three views '
          '| with the matrix | spread of the odors\' totals |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in x.itertuples():
        L_.append(f'| {r.point} | {_f(r.pn_inhadd, ".1f")} | {_f(r.pn_inhsc, ".1f")} | {_f(r.view_strength)} | '
                  f'{_f(r.view_sharing)} | {_f(r.view_breadth)} | {_f(r.view_matrix)} | {_f(r.error)} | '
                  f'{_f(r.error_with_matrix)} | {_f(r.totals_cv)} |')
        for col in ('view_strength', 'view_sharing', 'view_breadth', 'view_matrix', 'error', 'error_with_matrix',
                    'totals_cv', 'pn_inhadd', 'pn_inhsc'):
            V.add(f'pn_port_{r.point}', 'pn_dendrite', col, getattr(r, col))
    V.add('flies', 'pn_dendrite', 'totals_cv', L.FLIES_TOTALS_CV, half='all')
    return L_ + ['']


def signature_table(d: dict, V: W.Values) -> List[str]:
    """The recorded step's signature (step0_recorded_step_signature.*): each odor's residual off a proportional step
    against its summed ORN input, the flies' (with their interval) beside each point's, and the summary json's
    correlations and its check that the family's and A1's stages can be told apart."""
    s, j = d['step0_recorded_step_signature'], d['step0_recorded_step_signature_json']
    if not len(s):
        return []
    pts = [c for c in s.columns if c not in ('odor', 'total_orn', 'flies', 'flies_lo', 'flies_hi')]
    L_ = ["The step from the ORNs to the PN dendrites, odor by odor: each odor's residual off a proportional step (the "
          "log of its summed PN response over its summed ORN input, the panel's mean removed), the flies' with their "
          f"{L.SIGNATURE_CI:.0%} bootstrap interval over flies, and each point's on E135's port.", '',
          '| odor | summed ORN input | the flies (interval) | ' + ' | '.join(pts) + ' |',
          '|---|---:|---:|' + '---:|' * len(pts)]
    for r in s.itertuples():
        L_.append(f'| {r.odor} | {_f(r.total_orn, ".1f")} | {_f(r.flies, "+.3f")} ({_f(r.flies_lo, "+.3f")} to '
                  f'{_f(r.flies_hi, "+.3f")}) | ' + ' | '.join(_f(getattr(r, p), '+.3f') for p in pts) + ' |')
    rho = {p: float(pd.Series(s[p]).corr(s.total_orn, method='spearman')) for p in ['flies'] + pts}
    L_ += ['', 'ρ over the odors between the residual and the summed ORN input (the normalization\'s signature: '
           'negative where an odor with more input is turned down more): ' +
           ', '.join(f'{p} {rho[p]:+.2f}' for p in ['flies'] + pts) + '.']
    for p, v in rho.items():
        V.add(f'signature_{p}', 'orn_to_pn_dendrite', 'rho_residual_vs_total', v, half='all')
    sep = {k: v for k, v in (j or {}).items() if isinstance(v, (int, float, str, bool)) and k not in ('identity_against',)}
    if sep:
        L_ += ['', 'From step0_recorded_step_signature.json: ' + '; '.join(f'{k} {v}' for k, v in sep.items()) + '.']
    return L_ + ['']


def e102_ladder_table(d: dict, V: W.Values) -> List[str]:
    """E102's saved ladder scored as E136 scores a read (step0_e102_ladder.csv; the family's calyx with E102's breadth
    recomputed at each rung, the APL16 input): per step the mean over the eight wirings of the APL16 somata's strength
    and all-view errors with APL on and off, the calibrated pattern part against each target and its change from
    E102's identity point, and the KCs at E102's low-concentration stimulus."""
    e = d['step0_e102_ladder']
    if not len(e):
        return []
    w = e[e.wiring.isin(F.WIRINGS)]
    g = w.groupby(['variant', 'step'], sort=False)
    L_ = ["E102's saved ladder (the family's calyx on the APL16 input, its breadth recomputed at each rung) scored as "
          "E136 scores a read, mean over the eight wirings: the APL16 somata's strength error and error over all views "
          "with APL on and off (one-fly units), the calibrated pattern part against all recordings and its change from "
          "the identity point against each target, and the KCs' response rates at E102's low-concentration stimulus "
          "beside the flies'.", '',
          '| step | strength, APL on | APL off | all views, APL on | APL off | pattern part (all recordings) | change: '
          'all recordings | half A | half B | KCs at low concentration, APL on / off (the flies) |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    for (v, st), x in g:
        lo = (f"{_f(x.low_rate_on.median(), '.4f')} / {_f(x.low_rate_off.median(), '.4f')} "
              f"({_f(x.low_flies_on.iloc[0], '.4f')} / {_f(x.low_flies_off.iloc[0], '.4f')})"
              if 'low_rate_on' in x and x.low_rate_on.notna().any() else '–')
        L_.append(f'| {st} | {_f(x.strength_error_on.mean())} | {_f(x.strength_error_off.mean())} | '
                  f'{_f(x.all_error_on.mean())} | {_f(x.all_error_off.mean())} | {_f(x.price_z8_pattern.mean())} | '
                  f'{_f(x.d_price_z8_pattern.mean(), "+.3f")} | {_f(x.d_price_halfA_pattern.mean(), "+.3f")} | '
                  f'{_f(x.d_price_halfB_pattern.mean(), "+.3f")} | {lo} |')
        for col in ('strength_error_on', 'strength_error_off', 'all_error_on', 'all_error_off', 'price_z8_pattern',
                    'd_price_z8_pattern', 'd_price_halfA_pattern', 'd_price_halfB_pattern'):
            V.add(f'e102_{st}', 'kc_soma', f'{col}_mean8', x[col].mean(), panel='apl16', half='all')
    return L_ + ['']


# ---------------------------------------------------------------------------
# the score

def halves_table(d: dict, V: W.Values) -> List[str]:
    """The megamat stage error per variant on T against the held-out half, all flies and the fit half, and on S
    against the fit half (where a rung is named), means over the four wirings, with each setting's named rung."""
    m, c = d['mega'], d['classes']
    named = c.get('named_rungs') or {}
    L_ = ['The megamat stage error on strength (one-fly units, lower is better), mean over four wirings: on the held-out '
          'wirings T against each set of flies, and on the searched wirings S against the fit half (where each '
          "setting's rung is named; nothing is fitted, so S and T differ only in the wirings).", '',
          '| model | T, held-out half | T, all flies | T, fit half | S, fit half |', '|---|---:|---:|---:|---:|']
    for v in present(d):
        cells = []
        for half, ws in (('held', F.T_WIRINGS), ('all', F.T_WIRINGS), ('fit', F.T_WIRINGS), ('fit', F.S_WIRINGS)):
            x = m[(m.arm == v) & (m.half == half) & m.wiring.isin(ws)]
            cells.append(_f(x.stage.mean()) if len(x) else '–')
            if len(x):
                V.add(v, 'four_stages', 'stage_error_mean', x.stage.mean(),
                      half=half if ws is F.T_WIRINGS else f'{half}_S')
        L_.append(f'| {variant_label(v)} | ' + ' | '.join(cells) + ' |')
    if named:
        L_ += ['', 'The rung each setting names on S against the fit half (the lowest mean stage error over S): ' +
               '; '.join(f"{SETTING_NAMES.get(s, s)}: {rec.get('named'):g} (" +
                         ', '.join(f'{k} {val:.3f}' for k, val in (rec.get('values') or {}).items()) + ')'
                         for s, rec in named.items()) + '.']
        for s, rec in named.items():
            V.add(f'setting_{s}', 'four_stages', 'named_rung', rec.get('named'), half='fit_S')
    return L_ + ['']


def classes_table(d: dict, V: W.Values) -> List[str]:
    """Per variant, the class and what it rests on (classes.json): the change of the megamat stage error on T against
    the held-out half with its standard error and bar, the stages worse than their spread, the guard, the held-out
    panels worse, the holds, the price on the APL effect and "serves both"."""
    c = d['classes']
    cls = c.get('classes') or {}
    L_ = [f"The classes, on T against the held-out half, against the unrefitted family there (strength by odor; the bar "
          f"is {F.IMPROVE_SE:g} standard errors over wirings). The price is the largest of the three targets' mean "
          f"change over T of the calibrated pattern part; above {L.PRICE:g} it is a price. "
          f"{W.BESIDE_NOTE[0].upper() + W.BESIDE_NOTE[1:]}.", '',
          '| model | change of the megamat stage error (SE) | bar | improves | stages worse than their spread | guard: '
          'stages worse on all views | held-out panels worse | holds fail | price on the APL effect | serves both | '
          'reasons | class |', '|---|---:|---:|---|---|---|---|---|---:|---|---|---|']
    for v in present(d):
        r = cls.get(v)
        if not r:
            continue
        mg = r.get('megamat') or {}
        pr = r.get('price_on_apl_effect') or {}
        pv = pr.get('price_value')
        price = '–' if pv is None else f"{_f(pv, '+.4f')} ({_yn(pr.get('price'))})"
        if not mg.get('complete'):
            L_.append(f"| {variant_label(v)} | – | – | – | – | – | – | – | {price} | – | "
                      f"{'; '.join(r.get('reasons') or [])} | {r.get('cls')} |")
            continue
        L_.append(f"| {variant_label(v)} | {_f(mg['change'], '+.3f')} ({_f(mg['se'])}) | "
                  f"{_f(-F.IMPROVE_SE * mg['se'], '+.3f')} | {_yn(mg.get('improves'))} | "
                  f"{', '.join(STAGE_NAMES[s] for s in r.get('stage_worse') or []) or 'none'} | "
                  f"{', '.join(STAGE_NAMES.get(s, s) for s in r.get('guard_worse') or []) or 'none'} | "
                  f"{', '.join(W.PANEL_NAMES.get(k, k) for k in r.get('panel_worse') or []) or 'none'} | "
                  f"{len(r.get('holds_failed') or []) or 'none'} | {price} | {_yn(r.get('serves_both'))} | "
                  f"{'; '.join(r.get('reasons') or []) or '–'} | {r.get('cls')} |")
        V.add(v, 'four_stages', 'change', mg['change'])
        V.add(v, 'four_stages', 'se', mg['se'])
        V.add(v, 'four_stages', 'bar', F.IMPROVE_SE * mg['se'])
        V.add(v, 'four_stages', 'n_holds_failed', len(r.get('holds_failed') or []))
        if pv is not None:
            V.add(v, 'kc_soma', 'price_value', pv, panel='apl16', condition='apl_effect', half='T')
    cands = c.get('step3_candidates') or []
    L_ += ['', "Step 3's candidates (a rung below "
           f"{L.CANDIDATE_BELOW:g} at F or B that improves, with a price under {L.CANDIDATE_PRICE_SHARE:g} of A1's): "
           + ('; '.join(f"{x['variant']} (price {x['price_value']:+.4f} against A1's {x['a1_price_value']:+.4f})"
                        for x in cands) if cands else 'none') + '.']
    V.add('step3', 'four_stages', 'n_candidates', len(cands))
    return L_ + ['']


def price_table(d: dict, V: W.Values) -> List[str]:
    """The APL effect on the APL16 input by target (classes.json's price records and score_apl16.csv's price_*
    columns; E67's calibrated readout, K set per run): per variant the calibrated pattern part against all recordings,
    half A and half B, mean over T and over the eight wirings, the changes from the family over T, and beside them the
    amplitude and shape parts, the spread ratio, ρ and Pearson r on the per-pair APL effect against all recordings and
    the mean APL effect (means over the eight wirings)."""
    c, a = d['classes'], d['apl']
    cls = c.get('classes') or {}
    if not len(a):
        return []
    L_ = ["The APL effect on odor-odor distances on the APL16 input (held out; reported, and a price is a class beside "
          "the megamat class), calibrated readout. The pattern part is the loss with both means removed (the flies' "
          "floor about 0.029); \"T\" is the mean over the four held-out wirings, \"8\" over all eight; the changes are "
          "from the unrefitted family on the same wirings.", '',
          '| model | pattern part, all recordings: T / 8 | change on T: all recordings | half A | half B | amplitude '
          'part | shape part | spread ratio | ρ on the per-pair APL effect | Pearson r | the mean APL effect |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for v in present(d):
        x = a[a.arm == v]
        if not len(x):
            continue
        xt = x[x.wiring.isin(F.T_WIRINGS)]
        bt = ((cls.get(v) or {}).get('price_on_apl_effect') or {}).get('by_target') or {}
        g = lambda col, y=x: y[col].mean() if col in y and len(y) else np.nan
        ch = {t: (bt.get(t) or {}).get('change_T') for t in L.PRICE_TARGETS}
        L_.append(f"| {variant_label(v)} | {_f(g('price_z8_pattern', xt))} / {_f(g('price_z8_pattern'))} | "
                  + ' | '.join(_f(ch[t], '+.4f') if v != F.FAMILY else '–' for t in L.PRICE_TARGETS) + ' | '
                  f"{_f(g('price_z8_amplitude'))} | {_f(g('price_z8_shape'))} | {_f(g('price_z8_sd_ratio'), '.2f')} | "
                  f"{_f(g('price_z8_rho'), '.2f')} | {_f(g('price_z8_pearson'), '.2f')} | "
                  f"{_f(g('price_mean_effect'), '+.3f')} |")
        for t in L.PRICE_TARGETS:
            V.add(v, 'kc_soma', f'pattern_{t}_mean_T', g(f'price_{t}_pattern', xt), panel='apl16',
                  condition='apl_effect', half='T')
            V.add(v, 'kc_soma', f'pattern_{t}_mean_8', g(f'price_{t}_pattern'), panel='apl16', condition='apl_effect',
                  half='all')
            if v != F.FAMILY:
                V.add(v, 'kc_soma', f'pattern_{t}_change_T', ch[t], panel='apl16', condition='apl_effect', half='T')
                V.add(v, 'kc_soma', f'pattern_{t}_change_8', (bt.get(t) or {}).get('change_8'), panel='apl16',
                      condition='apl_effect', half='all')
        for col in ('amplitude', 'shape', 'sd_ratio', 'rho', 'pearson'):
            V.add(v, 'kc_soma', f'z8_{col}_mean_8', g(f'price_z8_{col}'), panel='apl16', condition='apl_effect',
                  half='all')
        V.add(v, 'kc_soma', 'mean_effect_mean_8', g('price_mean_effect'), panel='apl16', condition='apl_effect',
              half='all')
    return L_ + ['']


def extended_table(d: dict, V: W.Values) -> List[str]:
    """The beside step (score_beside_extended.csv; E102's readings of E98's extended input): per variant E98's V2 (held
    and its linear reading), V1, E100's W1 and the KCs' response rates at E102's carrying low-concentration stimulus
    beside the flies', with their classes against C2's; gating nothing."""
    e = d['extended']
    if not len(e):
        return []
    cols = [c for c in ('v2_held', 'v2_linear', 'v1', 'w1') if c in e.columns]
    L_ = ["Beside the class: the flies' APL along concentration and the KCs at low concentration on E98's extended "
          "input, read by E102's rules (medians over the wirings read; E102's own rows of its identity point and of the "
          "family's calyx at 1000, with E102's breadth, beside).", '',
          '| model | source | wirings | breadth | ' + ' | '.join(cols) + ' | KCs at low concentration, APL on / off '
          '(the flies) | class, APL on / off |', '|---|---|---:|---|' + '---|' * len(cols) + '---|---|']
    for r in e.to_dict('records'):
        lo = (f"{_f(r.get('low_on'), '.4f')} / {_f(r.get('low_off'), '.4f')} ({_f(r.get('low_flies_on'), '.4f')} / "
              f"{_f(r.get('low_flies_off'), '.4f')})")
        v = str(r['variant'])
        # E139's own variants by their labels; E102's rows ('family', 'F_inhadd1000 (E102)') as the score names them
        name = variant_label(v) if r.get('source') == 'E139:beside' else v
        L_.append(f"| {name} | {r.get('source')} | {int(r.get('n_wirings') or 0)} | {r.get('breadth', '–')} | "
                  + ' | '.join(str(r.get(c, '–')) if not isinstance(r.get(c), float) else _f(r.get(c), '.3f')
                               for c in cols)
                  + f" | {lo} | {r.get('low_class_on', '–')} / {r.get('low_class_off', '–')} |")
        for c in ('low_on', 'low_off', 'low_flies_on', 'low_flies_off', 'v2_linear_min', 'v2_linear_max'):
            V.add(str(r['variant']), 'kc_soma', f'extended_{c}', r.get(c), panel='extended', half='all')
    return L_ + ['']


def read_step(row: Mapping) -> str:
    """The step a read belongs to: 'beside' for the extended input, step 3's 'fit' for its arms, else the numbered step
    whose settings hold the read's setting (STEPS)."""
    if str(row.get('panel')) == L.EXTENDED or str(row.get('step')) == 'beside':
        return 'beside'
    s = str(row.get('setting'))
    if s == L.FIT_PREFIX:
        return 'fit'
    for k, (pairs, _) in L.STEPS.items():
        if any(s == ps for ps, _ in pairs):
            return k
    return '?'


def lanes_table(d: dict, V: W.Values) -> List[str]:
    """Which lane read what (reads_<host>_<lane>.csv): per host and lane, the steps, the number of reads, the queue,
    pyolfsysm and hong2p commits and whether any checkout was dirty."""
    r = d['reads']
    if not len(r):
        return []
    r = r.assign(step=[read_step(x) for x in r.to_dict('records')])
    L_ = ['The reads by lane (reads_<host>_<lane>.csv; the smoke lane left out).', '',
          '| host | lane | reads | steps | queue | pyolfsysm | hong2p | any checkout dirty |', '|---|---|---:|---|---|---|---|---|']
    for (h, ln), x in r.groupby(['host', 'lane'], sort=True):
        dirty = any(L._truthy(v) for col in ('queue_dirty', 'pyolfsysm_dirty', 'al_analysis_dirty', 'hong2p_dirty')
                    if col in x for v in x[col])
        steps = sorted(set(x.step))
        L_.append(f"| {h} | {ln} | {len(x)} | {', '.join(steps)} | "
                  f"{', '.join(sorted({str(c)[:7] for c in x.queue_commit}))} | "
                  f"{', '.join(sorted({str(c)[:7] for c in x.pyolfsysm_commit}))} | "
                  f"{', '.join(sorted({str(c)[:7] for c in x.hong2p_commit}))} | {_yn(dirty)} |")
    V.add('reads', 'all', 'n_reads', len(r), half='all')
    return L_ + ['']


def normalization_share(points: Mapping[str, Mapping[str, float]]) -> pd.DataFrame:
    """The additive term's share of the PN gain's denominator (olfsysm's sim_PN_layer: the gain is pn_inhsc / (pn_inhadd
    + 0.25 inhA + 0.75 inhB), inhA and inhB the LN's two filtered responses), per point on E135's port of the family's
    PN stage on the megamat input (the LN does not depend on the PN's normalization): its mean over the 17 odors in the
    stage reading's window, with the range over odors."""
    from experiments import E135_antennal_lobe_dendrites as AL
    st = L.family_pn_stage()
    inhA, inhB = st.ln(st.sp)
    win, _ = AL.windows(st.sp, 'stage')
    rows = []
    for name, p in points.items():
        add = float(p['pn_inhadd'])
        share = (add / (add + 0.25 * inhA[:, win] + 0.75 * inhB[:, win])).mean(axis=1)
        rows.append(dict(point=name, pn_inhadd=add, pn_inhsc=float(p['pn_inhsc']), share_mean=float(share.mean()),
                         share_min=float(share.min()), share_max=float(share.max())))
    return pd.DataFrame(rows)


def share_table(share: Optional[pd.DataFrame], V: W.Values) -> List[str]:
    if share is None or not len(share):
        return []
    L_ = ["The additive term's share of the PN gain's denominator (pn_inhadd over pn_inhadd plus the LN's filtered "
          "response), on the megamat input in the stage reading's window: the mean over the 17 odors (range over odors). "
          "olfsysm's default is 31.4.", '', '| point | pn_inhadd | pn_inhsc | share of the denominator |',
          '|---|---:|---:|---:|']
    for r in share.itertuples():
        L_.append(f'| {r.point} | {_f(r.pn_inhadd, ".1f")} | {_f(r.pn_inhsc, ".1f")} | {_f(r.share_mean)} '
                  f'({_f(r.share_min)} to {_f(r.share_max)}) |')
        V.add(f'pn_port_{r.point}', 'pn_dendrite', 'inhadd_share_mean', r.share_mean, half='all')
        V.add(f'pn_port_{r.point}', 'pn_dendrite', 'inhadd_share_min', r.share_min, half='all')
        V.add(f'pn_port_{r.point}', 'pn_dendrite', 'inhadd_share_max', r.share_max, half='all')
    return L_ + ['']


def operating_table(d: dict, V: W.Values) -> List[str]:
    """The megamat operating point of E139's own reads (reads_*.csv, panel megamat): per variant the response rates
    with APL on and off (range over wirings and halves) and the absolute scales the holds reached."""
    r = d['reads']
    if not len(r):
        return []
    x = r[r.panel == F.PANEL]
    if not len(x):
        return []
    cols = [(c, n) for c, n in W.SCALE_COLS if c in x.columns and x[c].notna().any()]
    L_ = ["The megamat operating point of E139's reads (range over the eight wirings and the halves; the flies' "
          "response rate with APL on 0.0915, held to within E136's tolerance) and the absolute scales the holds reached.",
          '', '| model | response rate, APL on | APL off | ' + ' | '.join(n for _, n in cols) + ' |',
          '|---|---:|---:|' + '---:|' * len(cols)]
    for v in [v for v in variant_order(list(dict.fromkeys(x.variant))) if v in set(x.variant)]:
        y = x[x.variant == v]
        L_.append(f"| {variant_label(v)} | {_rng(y.rate_on, '.4f')} | {_rng(y.rate_off, '.4f')} | "
                  + ' | '.join(_rng(y[c], '.4g') for c, _ in cols) + ' |')
        for c, _ in cols:
            V.add(v, 'megamat_reads', f'{c}_min', y[c].min(), half='all')
            V.add(v, 'megamat_reads', f'{c}_max', y[c].max(), half='all')
        V.add(v, 'megamat_reads', 'rate_on_min', y.rate_on.min(), half='all')
        V.add(v, 'megamat_reads', 'rate_on_max', y.rate_on.max(), half='all')
    return L_ + ['']


def main(score_dir: Optional[Path] = None, with_port: bool = True) -> Path:
    """Writes writeup_tables.md and writeup_values.csv into `score_dir` (default E139's out_dir)."""
    score_dir = L.out_dir() if score_dir is None else Path(score_dir)
    d = load(score_dir)
    V = W.Values()
    c, pv = d['classes'], d['prov']
    order = present(d)
    with L.e136() as _:
        share = normalization_share(L.pn_stage_points()) if with_port else None
        with as_e136(order):
            Lines = ['# E139 write-up tables', '',
                     f"From the score of {pv.get('time')} ({pv.get('host')}, queue {str(pv.get('queue_commit'))[:7]}): "
                     f"measure {c.get('measure')!r}, stages {', '.join(c.get('stages') or [])}; variants scored: "
                     f"{', '.join(pv.get('variants') or [])}; fitted: {', '.join(c.get('fitted') or []) or 'none'}.", '',
                     '## Step 0', ''] + a1_by_rung_table(d, V) + pn_ladder_table(d, V) + signature_table(d, V) \
                + e102_ladder_table(d, V)
            Lines += ['## The stage table', ''] + W.stage_table(d, V) + halves_table(d, V)
            Lines += ['## The classes', ''] + classes_table(d, V) + W.stage_changes_table(d, V)
            Lines += W.without_odor_table(d, d['without'] if len(d['without']) else None, V)
            Lines += ['## The held-out panels', ''] + W.panels_table(d, V) + W.e137_table(d, V)
            Lines += ['## The APL effect (held out) and the price', ''] + price_table(d, V)
            Lines += ['## Beside: the flies\' APL along concentration and the KCs at low concentration', ''] \
                + extended_table(d, V)
            Lines += ['## Realism', ''] + W.realism_table(d, V) + operating_table(d, V) + share_table(share, V)
            Lines += ['## The reads', ''] + lanes_table(d, V)
    V.frame().to_csv(score_dir / 'writeup_values.csv', index=False)
    out = score_dir / 'writeup_tables.md'
    out.write_text('\n'.join(Lines) + '\n')
    print(out)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--score-dir', type=Path, help="where the score's files are (default E139's out_dir)")
    ap.add_argument('--no-port', action='store_true', help="skip the normalization's share of the PN gain's "
                    "denominator (E135's port)")
    a = ap.parse_args()
    main(a.score_dir, with_port=not a.no_port)
