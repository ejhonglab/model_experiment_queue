"""The tables of E141's write-up (the rise of the inhibition among KCs beside the weaker antennal-lobe normalization,
along its dose and in other forms), from the files its parts wrote to results/E141_kckc_rise_forms/ (step 0's
step0_*.csv; the score's score_*.csv, classes.json and provenance_score.json; the reads' reads_<host>_<lane>.csv), from
E134's flies tables (the floors and the no-step yardstick), from E139's reads of B and K at 1000 (the operating point
beside), and for where APL's inhibition is made from the saved runs of the reads (`--no-sites` skips it). No simulation
of the mushroom body and no fit.

The stage table, the stage changes, the held-out panels, the mixture readings, the megamat somata without aa and the
realism bands are E136's write-up tables (analyses/E136_writeup_tables.py) on E141's files, and the halves, the APL
effect by target and the operating point are E139's (analyses/E139_writeup_tables.py), each with E141's variants in
place of their arms (`as_e136`, `as_e139`); step 0's tables, the classes with the entry's added classes and the
comparisons with A1, the dose ladder, the spread at the somata, the beside step and the reads by lane are E141's own.

Writes writeup_tables.md and writeup_values.csv (the same numbers, long: arm, stage, panel, condition, half, statistic,
value, 'arm' being the variant, for the summary's placeholders) beside the score's files.

Usage (from the repository root): ./venv/bin/python -m analyses.E141_writeup_tables [--score-dir <dir>] [--no-sites]
"""
from __future__ import annotations

import argparse
import contextlib
from pathlib import Path
from typing import Iterator, List, Optional, Sequence

import numpy as np
import pandas as pd

from analyses import E136_writeup_tables as W
from analyses import E139_writeup_tables as W139
from experiments import E136_stage_population_fit as F
from experiments import E139_normalization_ladder as E139
from experiments import E141_kckc_rise_forms as L

STAGES = W.STAGES
STAGE_NAMES = W.STAGE_NAMES
_f, _mr, _rng, _yn = W._f, W._mr, W._rng, W._yn

# the dose ladder (KC-KC x1 is E139's B at 1000, x3.24 is A1), then the forms, then E139's K at 1000 beside
DOSE_LADDER = (('B_inhadd1000', 1.0), (L.variant('D1.8'), L.ARMS['D1.8'].kckc), (L.A1, L.RISE),
               (L.variant('D5.8'), L.ARMS['D5.8'].kckc))
FORM_ORDER = ('N', 'G', 'S', 'T', 'P', 'M')
VARIANT_ORDER = ((F.FAMILY, 'F_inhadd1000') + tuple(v for v, _ in DOSE_LADDER) + tuple(L.variant(a) for a in FORM_ORDER)
                 + ('K_inhadd1000',))
# the somata's panels of the spread tables: (input, panel, condition, name); 'megamat' from the APL16 input is the
# megamat flies on the APL16 panel's 12 megamat odors
SOMA_PANELS = (('validation2', 'validation2', 'on', 'validation2'), ('mixtures', 'cmix', 'on', 'control mixtures'),
               ('mixtures', 'kmix', 'on', 'kiwi mixtures'), ('apl16', 'apl16', 'on', 'APL16, APL on'),
               ('apl16', 'apl16', 'off', 'APL16, APL off'),
               ('apl16', 'megamat', 'on', "megamat flies, the APL16 input's 12 odors"))


def variant_label(v: str) -> str:
    """A variant's name in the tables."""
    if v == F.FAMILY:
        return 'the family (unrefitted, pn_inhadd 31.4)'
    if v == L.A1:
        return 'A1 (E136): pn_inhadd 1017, APL→bouton ×0.155, KC–KC ×3.24'
    if v == 'B_inhadd1000':
        return "KC–KC ×1 (E139's B at 1000)"
    if v == 'K_inhadd1000':
        return "KC–KC ×3.24 at the family's APL→bouton (E139's K at 1000)"
    if v == 'F_inhadd1000':
        return "the family's calyx at 1000 (E139's F)"
    try:
        a = L.parse_variant(v)
    except ValueError:
        return v
    names = {'D1.8': f"D1.8: KC–KC ×{L.ARMS['D1.8'].kckc:.3g}", 'D5.8': f"D5.8: KC–KC ×{L.ARMS['D5.8'].kckc:.3g}",
             'N': f'N: KC–KC ×{L.RISE:g} with APL on, ×1 with APL off',
             'G': "G: APL's inhibition of the claws as a gain, KC–KC ×1",
             'P': f'P: one uniform neuron at the summed weight of ×{L.RISE:g}',
             'T': f'T: one neuron per KC type at the summed weight of ×{L.RISE:g}',
             'S': f'S: the matrix shuffled within type pairs at the summed weight of ×{L.RISE:g}',
             'M': "M: MB-C1's two cells beside KC–KC ×1"}
    return names.get(a, v)


def variant_order(names: Sequence[str]) -> List[str]:
    """`names` in VARIANT_ORDER's order, others after it."""
    have = list(dict.fromkeys(names))
    return [v for v in VARIANT_ORDER if v in have] + [v for v in have if v not in VARIANT_ORDER]


def present(d: dict) -> List[str]:
    return variant_order(list(d['mega'].arm)) if len(d['mega']) else []


@contextlib.contextmanager
def as_e136(order: Sequence[str]) -> Iterator[None]:
    """For the block, E136's write-up tables iterate `order` and name the variants by `variant_label`."""
    old = W.ARM_ORDER, W.ARM_NAMES
    W.ARM_ORDER = tuple(order)
    W.ARM_NAMES = dict(old[1], **{v: variant_label(v) for v in order})
    try:
        yield
    finally:
        W.ARM_ORDER, W.ARM_NAMES = old


@contextlib.contextmanager
def as_e139() -> Iterator[None]:
    """For the block, E139's write-up tables order and name the variants as E141's (`variant_order`, `present`,
    `variant_label`)."""
    old = W139.variant_order, W139.present, W139.variant_label
    W139.variant_order, W139.present, W139.variant_label = variant_order, present, variant_label
    try:
        yield
    finally:
        W139.variant_order, W139.present, W139.variant_label = old


# ---------------------------------------------------------------------------
# inputs

def load(score_dir: Path, results_dir: Optional[Path] = None) -> dict:
    """Everything the tables read, in the form E136's and E139's tables take ('arm' the variant): the score's files from
    `score_dir`, step 0's and the reads' from `results_dir` (default E141's out_dir), with E139's reads of its variants
    scored beside (`L.E139_BESIDE`) in 'reads_e139'."""
    results_dir = L.out_dir() if results_dir is None else Path(results_dir)
    d = W139.load(score_dir, results_dir)
    for name in ('soma_splits', 'holds', 'a1_by_kckc'):
        d[f'step0_{name}'] = W139._csv(results_dir / f'step0_{name}.csv')
    r = E139.reads_table()
    d['reads_e139'] = r[r.variant.isin(L.E139_BESIDE)] if len(r) and 'variant' in r else pd.DataFrame()
    return d


# ---------------------------------------------------------------------------
# step 0

def soma_spread_table(d: dict, V: W.Values) -> List[str]:
    """Step 0: E139's saved reads at each panel's somata (step0_soma_splits.csv; all flies, the mean over T): the model's
    spread of the odors' response rates as a multiple of the flies' (sd_ratio) and r over the odors."""
    s = d['step0_soma_splits']
    if not len(s):
        return []
    s = s[(s.half == 'all') & (s.wiring == 'mean_T') & s.odor.isna()]
    L_ = ["Step 0, the spread at the somata: the model's SD over odors of the odors' response rates as a multiple of the "
          "flies' (and r over the odors), E139's saved reads, all flies, the mean over T. It describes; it gates nothing.",
          '', '| model | ' + ' | '.join(n for *_, n in SOMA_PANELS) + ' |', '|---|' + '---:|' * len(SOMA_PANELS)]
    for v in variant_order(list(s.variant)):
        cells = []
        for inp, pan, cond, _ in SOMA_PANELS:
            x = s[(s.variant == v) & (s.input == inp) & (s.panel == pan) & (s.condition == cond)].set_index(
                'statistic').value
            sd, r = x.get('strength_rate_sd_ratio'), x.get('strength_rate_r')
            cells.append(f"{_f(sd, '.2f')} ({_f(r, '.2f')})" if sd is not None else '–')
            V.add(v, 'kc_soma', 'step0_sd_ratio_mean_T', sd, panel=f'{inp}/{pan}', condition=cond, half='all')
            V.add(v, 'kc_soma', 'step0_r_mean_T', r, panel=f'{inp}/{pan}', condition=cond, half='all')
        L_.append(f'| {variant_label(v)} | ' + ' | '.join(cells) + ' |')
    return L_ + ['']


def soma_split_table(d: dict, V: W.Values) -> List[str]:
    """Step 0: the somata's strength error on response rates split into level, spread and order (E136's
    `strength_split`; E139's saved reads, all flies, the mean over T)."""
    s = d['step0_soma_splits']
    if not len(s):
        return []
    s = s[(s.half == 'all') & (s.wiring == 'mean_T') & s.odor.isna()]
    L_ = ["Step 0, the somata's strength error on response rates split into its level, spread and order parts (E136's "
          "`strength_split`, squared parts adding to the squared error), all flies, the mean over T: level / spread / "
          "order.", '', '| model | ' + ' | '.join(n for *_, n in SOMA_PANELS) + ' |', '|---|' + '---:|' * len(SOMA_PANELS)]
    for v in variant_order(list(s.variant)):
        cells = []
        for inp, pan, cond, _ in SOMA_PANELS:
            x = s[(s.variant == v) & (s.input == inp) & (s.panel == pan) & (s.condition == cond)].set_index(
                'statistic').value
            parts = [x.get(f'strength_rate_{k}') for k in ('level', 'spread', 'order')]
            cells.append(' / '.join(_f(p, '.2f') for p in parts) if parts[0] is not None else '–')
            for k, p in zip(('level', 'spread', 'order'), parts):
                V.add(v, 'kc_soma', f'step0_rate_{k}_mean_T', p, panel=f'{inp}/{pan}', condition=cond, half='all')
        L_.append(f'| {variant_label(v)} | ' + ' | '.join(cells) + ' |')
    return L_ + ['']


def holds_table(d: dict, V: W.Values) -> List[str]:
    """Step 0: what E139's holds reach (step0_holds.csv), range over the wirings: the somatic threshold, the KC-KC and
    APL->KC scales, the response rates of the three KC types and the odors per KC's distance from the flies', on the
    megamat input with APL on and on the APL16 input with APL on and off."""
    h = d['step0_holds']
    if not len(h):
        return []
    L_ = ["Step 0, what E139's holds reach (absolute scales; range over the wirings): megamat with APL on (the held-out "
          "half's reads), APL16 with APL on and off.", '',
          '| model | input, condition | somatic threshold | KC–KC | APL→KC | response rate: γ / α′β′ / αβ | odors per KC, '
          'distance from the flies |', '|---|---|---:|---:|---:|---|---:|']
    for v in variant_order(list(h.variant)):
        for inp, cond in (('megamat', 'on'), ('apl16', 'on'), ('apl16', 'off')):
            x = h[(h.variant == v) & (h.input == inp) & (h.condition == cond)]
            if not len(x):
                continue
            L_.append(f"| {variant_label(v)} | {inp}, APL {cond} | {_rng(x.fixed_thr, '.1f')} | {_rng(x.wKCKC_scale, '.1f')} | "
                      f"{_rng(x.wAPLKC_scale, '.3f')} | {_rng(x.rate_g, '.3f')} / {_rng(x.rate_apbp, '.3f')} / "
                      f"{_rng(x.rate_ab, '.3f')} | {_rng(x.opk_tv, '.3f')} |")
            for col in ('fixed_thr', 'wKCKC_scale', 'wAPLKC_scale', 'rate_g', 'rate_apbp', 'rate_ab', 'opk_tv'):
                if col in x and x[col].notna().any():
                    V.add(v, 'step0_holds', f'{col}_min', x[col].min(), panel=inp, condition=cond, half='all')
                    V.add(v, 'step0_holds', f'{col}_max', x[col].max(), panel=inp, condition=cond, half='all')
    return L_ + ['']


def a1_by_kckc_table(d: dict, V: W.Values) -> List[str]:
    """Step 0: A1's saved evaluations binned by the KC-KC factor (step0_a1_by_kckc.csv): per bin the number of
    evaluations and the best and median megamat stage error on S against the fit half and on T against the held-out
    half, with the KC somata's (a lead: the other factors differ along the bins)."""
    a = d['step0_a1_by_kckc']
    if not len(a):
        return []
    L_ = ["Step 0, A1's saved evaluations binned by the KC–KC factor (a lead: the other factors differ along the bins): "
          "the megamat stage error, best / median in the bin (the KC somata's in brackets).", '',
          '| KC–KC factor | evaluations | S against the fit half | T against the held-out half |', '|---|---:|---:|---:|']
    for b, x in a.groupby('bin', sort=True):
        def cell(half):
            y = x[x.half == half].set_index('statistic')
            if not len(y):
                return '–'
            return (f"{_f(y.stage_error.get('best'))} / {_f(y.stage_error.get('median'))} "
                    f"({_f(y.kc_soma.get('best'))} / {_f(y.kc_soma.get('median'))})")
        r = x.iloc[0]
        L_.append(f'| {r.bin_lo:.2f} to {r.bin_hi:.2f} | {int(r.n_evals)} | {cell("fit")} | {cell("held")} |')
        for half in ('fit', 'held'):
            y = x[x.half == half].set_index('statistic')
            for st in ('best', 'median'):
                if st in y.index:
                    V.add(f'A1_kckc_bin{int(b)}', 'four_stages', f'stage_error_{st}', y.stage_error[st], half=half)
                    V.add(f'A1_kckc_bin{int(b)}', 'kc_soma', f'strength_{st}', y.kc_soma[st], half=half)
        V.add(f'A1_kckc_bin{int(b)}', 'four_stages', 'n_evals', r.n_evals, half='fit')
    return L_ + ['']


# ---------------------------------------------------------------------------
# the classes

def _keys() -> List[str]:
    return ['/'.join(k) for k in F.SCORED_PANELS if '/'.join(k) not in W.BESIDE_KEYS]


def classes_table(d: dict, V: W.Values) -> List[str]:
    """Per variant, E136's class and what it rests on (classes.json), the held-out panel nearest its spread (the spread
    less the change; a near-tie when small), the entry's added classes and the number of arms read."""
    c = d['classes']
    cls = c.get('classes') or {}
    keys = _keys()
    L_ = [f"The classes, on T against the held-out half, against the unrefitted family there (strength by odor; the bar is "
          f"{F.IMPROVE_SE:g} standard errors over wirings). \"Nearest\" is the held-out panel in the class whose change "
          f"comes nearest its spread, with the spread less the change (a near-tie when small). "
          f"{W.BESIDE_NOTE[0].upper() + W.BESIDE_NOTE[1:]}. Arms read: {c.get('arms_read')}; holding: "
          f"{', '.join(c.get('arms_holding') or []) or 'none'}.", '',
          '| model | change of the megamat stage error (SE) | improves | guard: stages worse on all views | held-out panels '
          'worse | nearest | holds fail | stands in | does more than A1 | a lower price | reasons | class |',
          '|---|---:|---|---|---|---|---|---|---|---|---|---|']
    for v in variant_order(list(cls)):
        r = cls[v]
        mg = r.get('megamat') or {}
        extra = ' | '.join(_yn(r.get(k)) if r.get(k) is not None else '–'
                           for k in ('stands_in', 'does_more_than_a1', 'lower_price'))
        if not mg.get('complete'):
            L_.append(f"| {variant_label(v)} | – | – | – | – | – | {len(r.get('holds_failed') or []) or 'none'} | {extra} | "
                      f"{'; '.join(r.get('reasons') or []) or '–'} | {r.get('cls')} |")
            V.add(v, 'four_stages', 'n_holds_failed', len(r.get('holds_failed') or []))
            V.add(v, 'four_stages', 'n_reads', r.get('n_reads'))
            continue
        margins = {k: (r['panels'][k]['spread'] - r['panels'][k]['change']) for k in keys
                   if (r.get('panels') or {}).get(k, {}).get('complete')}
        near = min(margins, key=margins.get) if margins else None
        L_.append(f"| {variant_label(v)} | {_f(mg['change'], '+.3f')} ({_f(mg['se'])}) | {_yn(mg.get('improves'))} | "
                  f"{', '.join(STAGE_NAMES.get(s, W.PANEL_NAMES.get(s, s)) for s in r.get('guard_worse') or []) or 'none'} | "
                  f"{', '.join(W.PANEL_NAMES.get(k, k) for k in r.get('panel_worse') or []) or 'none'} | "
                  + (f"{W.PANEL_NAMES.get(near, near)}: {_f(margins[near], '+.5f')}" if near else '–') + ' | '
                  f"{len(r.get('holds_failed') or []) or 'none'} | {extra} | {'; '.join(r.get('reasons') or []) or '–'} | "
                  f"{r.get('cls')} |")
        V.add(v, 'four_stages', 'change', mg['change'])
        V.add(v, 'four_stages', 'se', mg['se'])
        V.add(v, 'four_stages', 'bar', F.IMPROVE_SE * mg['se'])
        for k, m in margins.items():
            V.add(v, 'kc_soma', 'margin_to_spread', m, panel=k)
    V.add('arms', 'all', 'n_arms_read', c.get('arms_read'), half='all')
    V.add('arms', 'all', 'n_arms_holding', len(c.get('arms_holding') or []), half='all')
    return L_ + ['']


def against_a1_table(d: dict, V: W.Values) -> List[str]:
    """Per arm, the comparisons with A1 (classes.json's 'against_a1'): the change of the megamat stage error from A1
    on T against the held-out half (SE), each held-out panel's change from A1 (better or worse than its spread), and the
    price's change from A1 by target with "a lower price"."""
    cls = d['classes'].get('classes') or {}
    keys = _keys()
    rows = [(v, cls[v]['against_a1']) for v in variant_order(list(cls)) if (cls[v].get('against_a1') or {})]
    if not rows:
        return []
    L_ = ["Each arm against A1 on T against the held-out half: the change of the megamat stage error from A1 (standard "
          "error), each held-out panel's change from A1 (↓ better and ↑ worse than the panel's spread), and the calibrated "
          f"pattern part's change from A1 against all recordings / half A / half B (\"lower\" below −{L.LOWER_PRICE:g} in "
          "all three).", '',
          '| arm | megamat stage error | ' + ' | '.join(W.PANEL_NAMES[k] for k in keys) + ' | the price, from A1 | lower |',
          '|---|---:|' + '---:|' * len(keys) + '---:|---|']
    for v, a in rows:
        mg = a.get('megamat') or {}
        cells = []
        for k in keys:
            x = (a.get('panels') or {}).get(k) or {}
            mark = ' ↓' if x.get('better') else (' ↑' if x.get('worse') else '')
            cells.append(f"{_f(x.get('change'), '+.3f')}{mark}" if x else '–')
            V.add(v, 'kc_soma', 'change_from_a1', x.get('change'), panel=k)
        pr = a.get('price') or {}
        bt = pr.get('by_target') or {}
        L_.append(f"| {variant_label(v)} | {_f(mg.get('change'), '+.4f')} ({_f(mg.get('se'))}) | " + ' | '.join(cells)
                  + f" | {' / '.join(_f(bt.get(t), '+.4f') for t in E139.PRICE_TARGETS)} | {_yn(pr.get('lower'))} |")
        V.add(v, 'four_stages', 'change_from_a1', mg.get('change'))
        for t in E139.PRICE_TARGETS:
            V.add(v, 'kc_soma', f'price_{t}_change_from_a1', bt.get(t), panel='apl16', condition='apl_effect', half='T')
    return L_ + ['']


def dose_table(d: dict, V: W.Values) -> List[str]:
    """The dose ladder (KC-KC x1, x1.8, x3.24, x5.83 at the base point; x1 is E139's B at 1000 and x3.24 A1): the
    megamat stage error and each held-out panel's strength error on T (mean over T; change from the family), the
    KC-KC scale the APL16 holds reached and the class."""
    c = d['classes']
    cls, sp = c.get('classes') or {}, c.get('spreads') or {}
    keys = _keys()
    r = pd.concat([d['reads'], d['reads_e139']], ignore_index=True, sort=False) if len(d['reads']) else d['reads_e139']
    L_ = ["The dose ladder at the base point (pn_inhadd 1000, APL→bouton ×0.155): strength by odor on T, the mean over T "
          "with the change from the unrefitted family (the spreads in the second row; ↑ worse than the spread). A1 is "
          "its own fit's point (pn_inhadd 1017). The KC–KC scale is the absolute scale the APL16 holds reached (range "
          "over the eight wirings; A1's from step 0's table).", '',
          '| KC–KC | model | megamat stage error | ' + ' | '.join(W.PANEL_NAMES[k] for k in keys)
          + ' | KC–KC scale, APL16 | class |', '|---:|---|---:|' + '---:|' * len(keys) + '---:|---|',
          '| | spread | ' + f"{_f(sp.get('stage'))} | " + ' | '.join(_f(sp.get(k)) for k in keys) + ' | | |']
    h = d['step0_holds']
    for v, k_ in DOSE_LADDER:
        rec = cls.get(v) or {}
        mg = rec.get('megamat') or {}
        cells = []
        for k in keys:
            x = (rec.get('panels') or {}).get(k) or {}
            mark = ' ↑' if x.get('worse_than_spread') else ''
            cells.append(f"{_f(x.get('arm'))} ({_f(x.get('change'), '+.3f')}{mark})" if x.get('complete') else '–')
        if v == L.A1 and len(h):
            y = h[(h.variant == v) & (h.input == 'apl16') & (h.condition == 'on')].wKCKC_scale
        else:
            y = r[(r.variant == v) & (r.panel == 'apl16')].wKCKC_scale if len(r) else pd.Series(dtype=float)
        L_.append(f"| ×{k_:.3g} | {variant_label(v)} | {_f(mg.get('arm'))} ({_f(mg.get('change'), '+.3f')}) | "
                  + ' | '.join(cells) + f" | {_rng(y, '.1f')} | {rec.get('cls', '–')} |")
        V.add(v, 'four_stages', 'kckc_factor', k_)
        V.add(v, 'four_stages', 'stage_error_mean', mg.get('arm'))
        if len(y):
            V.add(v, 'kc_soma', 'wKCKC_scale_min', y.min(), panel='apl16', half='all')
            V.add(v, 'kc_soma', 'wKCKC_scale_max', y.max(), panel='apl16', half='all')
    return L_ + ['']


def spread_table(d: dict, V: W.Values) -> List[str]:
    """The spread at the somata of every scored model (score_panels.csv's strength split; the mean over T): the
    model's SD over odors of the odors' response rates as a multiple of the flies', with r over the odors."""
    p = d['panels']
    if not len(p) or 'split_strength_rate_sd_ratio' not in p:
        return []
    p = p[p.stage == 'kc_soma'] if 'stage' in p else p
    p = p[p.wiring.isin(F.T_WIRINGS)]
    L_ = ["The spread at the somata (the entry's definition; it describes and gates nothing): the model's SD over odors "
          "of the odors' response rates as a multiple of the flies' (r over the odors), the mean over T.", '',
          '| model | ' + ' | '.join(n for *_, n in SOMA_PANELS) + ' |', '|---|' + '---:|' * len(SOMA_PANELS)]
    for v in variant_order(list(p.arm)):
        cells = []
        for inp, pan, cond, _ in SOMA_PANELS:
            x = p[(p.arm == v) & (p.input == inp) & (p.panel == pan) & (p.condition == cond)]
            sd, r = (x.split_strength_rate_sd_ratio.mean(), x.split_strength_rate_r.mean()) if len(x) else (None, None)
            cells.append(f"{_f(sd, '.2f')} ({_f(r, '.2f')})" if len(x) else '–')
            if len(x):
                V.add(v, 'kc_soma', 'sd_ratio_mean_T', sd, panel=f'{inp}/{pan}', condition=cond, half='all')
                V.add(v, 'kc_soma', 'r_mean_T', r, panel=f'{inp}/{pan}', condition=cond, half='all')
        L_.append(f'| {variant_label(v)} | ' + ' | '.join(cells) + ' |')
    return L_ + ['']


# ---------------------------------------------------------------------------
# beside, realism and the reads

def extended_table(d: dict, V: W.Values) -> List[str]:
    """The beside step (score_beside_extended.csv; E102's readings of E98's extended input), as E139's table but with
    E141's names: per variant E98's V2 (held and its linear reading), V1, E100's W1 and the KCs' response rates at
    E102's carrying low-concentration stimulus beside the flies'; gating nothing."""
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
        v = str(r['variant'])
        lo = (f"{_f(r.get('low_on'), '.4f')} / {_f(r.get('low_off'), '.4f')} ({_f(r.get('low_flies_on'), '.4f')} / "
              f"{_f(r.get('low_flies_off'), '.4f')})")
        name = variant_label(v) if str(r.get('source')) in ('E141:beside', 'E139:beside') else v
        L_.append(f"| {name} | {r.get('source')} | {int(r.get('n_wirings') or 0)} | {r.get('breadth', '–')} | "
                  + ' | '.join(str(r.get(c, '–')) if not isinstance(r.get(c), float) else _f(r.get(c), '.3f')
                               for c in cols)
                  + f" | {lo} | {r.get('low_class_on', '–')} / {r.get('low_class_off', '–')} |")
        for c in ('low_on', 'low_off', 'low_flies_on', 'low_flies_off', 'v2_linear_min', 'v2_linear_max', 'n_wirings'):
            V.add(v, 'kc_soma', f'extended_{c}', r.get(c), panel='extended', half='all')
    return L_ + ['']


def operating_table(d: dict, V: W.Values) -> List[str]:
    """E139's operating table (the megamat operating point and the absolute scales the holds reached) on E141's reads
    joined by E139's reads of the variants scored beside."""
    r = pd.concat([d['reads'], d['reads_e139']], ignore_index=True, sort=False)
    out = W139.operating_table(dict(d, reads=r), V)
    if out:
        out[0] = out[0].replace("of E139's reads", "of the reads (E141's, with E139's of B and K at 1000)")
    return out


def apl_sites(names: Sequence[str]) -> pd.DataFrame:
    """Where APL's inhibition is made (E136's write-up's `sites_row`, E94's shares) in each model's reads: megamat on T
    against the held-out half, APL16 on the eight wirings; E141's and E139's variants registered as the score registers
    them, G (the gain at the claws) with the divisor at the responders' claws. Runs without stage windows are left out."""
    from kcq import stage_readings as R
    L.register_variants(F, [v for v in names if _is_e141(v)])
    E139.register_variants(F, [v for v in names if v in L.E139_BESIDE])
    rows = []
    for v in names:
        div = _is_e141(v) and L.ARMS[L.parse_variant(v)].form == 'gain'
        for panel, ws, half in ((F.PANEL, F.T_WIRINGS, 'held'), ('apl16', F.WIRINGS, None)):
            for w in ws:
                r = F.load_read(F.ReadPoint(v, panel, w, half))
                if r is None:
                    continue
                on, off = (r.conditions[c] for c in R.CONDITIONS)
                if not all(k in on and k in off for k in ('bouton_win', 'claw_win')):
                    continue
                rows.append(dict(arm=v, panel=panel, wiring=w, **W.sites_row(W.wiring_W(w), on, off, div)))
    return pd.DataFrame(rows)


def _is_e141(v: str) -> bool:
    try:
        L.parse_variant(v)
        return True
    except ValueError:
        return False


def read_step(row) -> str:
    """The step a read belongs to: 'beside' for the extended input, else step 1 or 2 by its arm."""
    tag = str(row.get('tag'))
    if str(row.get('panel')) == L.EXTENDED or f'_{L.EXTENDED}_' in tag:
        return 'beside'
    a = tag.split('_')[0]
    return '1' if a in L.STEP1 else ('2' if a in L.STEP2_IF_P + L.STEP2_ELSE else '?')


def lanes_table(d: dict, V: W.Values) -> List[str]:
    """Which lane read what (reads_<host>_<lane>.csv): per host and lane, the reads, the steps and arms, the queue,
    pyolfsysm and hong2p commits and whether any checkout was dirty."""
    r = d['reads']
    if not len(r):
        return []
    r = r.assign(step=[read_step(x) for x in r.to_dict('records')], a=r.tag.str.split('_').str[0])
    L_ = ['The reads by lane (reads_<host>_<lane>.csv).', '',
          '| host | lane | reads | steps | arms | queue | pyolfsysm | hong2p | any checkout dirty |',
          '|---|---|---:|---|---|---|---|---|---|']
    for (h, ln), x in r.groupby(['host', 'lane'], sort=True):
        dirty = any(E139._truthy(v) for col in ('queue_dirty', 'pyolfsysm_dirty', 'al_analysis_dirty', 'hong2p_dirty')
                    if col in x for v in x[col])
        L_.append(f"| {h} | {ln} | {len(x)} | {', '.join(sorted(set(x.step)))} | {', '.join(sorted(set(x.a)))} | "
                  f"{', '.join(sorted({str(c)[:7] for c in x.queue_commit}))} | "
                  f"{', '.join(sorted({str(c)[:7] for c in x.pyolfsysm_commit}))} | "
                  f"{', '.join(sorted({str(c)[:7] for c in x.hong2p_commit}))} | {_yn(dirty)} |")
    V.add('reads', 'all', 'n_reads', len(r), half='all')
    for s, x in r.groupby('step'):
        V.add('reads', 'all', f'n_reads_step_{s}', len(x), half='all')
    return L_ + ['']


def main(score_dir: Optional[Path] = None, with_sites: bool = True) -> Path:
    """Writes writeup_tables.md and writeup_values.csv into `score_dir` (default E141's out_dir)."""
    score_dir = L.out_dir() if score_dir is None else Path(score_dir)
    d = load(score_dir)
    V = W.Values()
    c, pv = d['classes'], d['prov']
    order = present(d)
    with E139.e136() as _:
        sites = apl_sites(order) if with_sites else None
        with as_e136(order), as_e139():
            Lines = ['# E141 write-up tables', '',
                     f"From the score of {pv.get('time')} ({pv.get('host')}, queue {str(pv.get('queue_commit'))[:7]}): "
                     f"measure {c.get('measure')!r}, stages {', '.join(c.get('stages') or [])}; variants scored: "
                     f"{', '.join(pv.get('variants') or [])}; E139's beside: {', '.join(c.get('beside_e139') or [])}.",
                     '', '## Step 0', ''] + soma_spread_table(d, V) + soma_split_table(d, V) + holds_table(d, V) \
                + a1_by_kckc_table(d, V)
            Lines += ['## The stage table', ''] + W.stage_table(d, V) + W139.halves_table(d, V)
            Lines += ['## The classes', ''] + classes_table(d, V) + against_a1_table(d, V) + W.stage_changes_table(d, V)
            Lines += ['## The dose ladder and the spread at the somata', ''] + dose_table(d, V) + spread_table(d, V)
            Lines += W.without_odor_table(d, d['without'] if len(d['without']) else None, V)
            Lines += ['## The held-out panels', ''] + W.panels_table(d, V) + W.e137_table(d, V)
            Lines += ['## The APL effect (held out) and the price', ''] + W139.price_table(d, V)
            Lines += ["## Beside: the flies' APL along concentration and the KCs at low concentration", ''] \
                + extended_table(d, V)
            Lines += ['## Realism', ''] + W.realism_table(d, V) + operating_table(d, V)
            Lines += W.apl_sites_table(d, sites, V)
            Lines += ['## The reads', ''] + lanes_table(d, V)
    V.frame().to_csv(score_dir / 'writeup_values.csv', index=False)
    out = score_dir / 'writeup_tables.md'
    out.write_text('\n'.join(Lines) + '\n')
    print(out)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--score-dir', type=Path, help="where the score's files are (default E141's out_dir)")
    ap.add_argument('--no-sites', action='store_true', help="skip where APL's inhibition is made (reads the saved runs)")
    a = ap.parse_args()
    main(a.score_dir, with_sites=not a.no_sites)
