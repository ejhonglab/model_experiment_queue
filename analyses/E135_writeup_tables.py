"""The tables of E135's write-up (the antennal lobe held to the recorded PN dendrites), from the files its parts wrote to
results/E135_antennal_lobe_dendrites/ and from E134's tables; no simulation and no scoring.

Writes results/E135_antennal_lobe_dendrites/writeup_tables.md, and writeup_values.csv (the same numbers, long, for the
summary's placeholders):
- the stage table: strength by odor at each megamat stage (one-fly units, the default reading, the mean over the
  eight wirings with their range) for E134's held runs and every step-2 arm, and the megamat stage error;
- the views beside (sharing, breadth, the three-view population error) for the same arms;
- the APL16 somata, APL on and off, against the family's own runs (E134's step 2, E90:B);
- step 1's rungs on the held-out flies and odors, both directions, and the controls;
- ruling 5's scale factors, and R1's and the kept rung's as-ruled tables against their fits on all 17 odors;
- the imposed checks (ruling 8), the routes (ruling 9), the realism of the APL16 runs with the family's own beside, and
  their APL effect (held out).

Usage (from the repository root): ./venv/bin/python -m analyses.E135_writeup_tables
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from experiments import E135_antennal_lobe_dendrites as E

STAGES = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
STAGE_NAMES = {'pn_dendrite': 'PN dendrites', 'pn_bouton': 'PN boutons', 'kc_claw': 'KC claws',
               'kc_soma': 'KC somata'}
ARM_NAMES = {
    'family': 'the family (its own fits)',
    'imposed': 'the recorded dendrites imposed',
    'imposed_glom_perm': 'the same, glomeruli permuted',
    'imposed_odor_perm': 'the same, odors permuted',
    'R1': 'R1 (global parameters fitted)',
    'rung': 'R3 (the highest kept rung)',
    'rung_glom_perm': 'R3, glomeruli permuted',
    'imposed_route2_repeat': 'route 2 (solved input), real wiring only',
    'imposed_dff': 'the recorded dendrites imposed, table in dF/F',
    'imposed_glom_perm_dff': 'the same in dF/F, glomeruli permuted',
    'imposed_odor_perm_dff': 'the same in dF/F, odors permuted',
}
ARM_ORDER = ('imposed', 'imposed_glom_perm', 'imposed_odor_perm', 'imposed_dff', 'imposed_glom_perm_dff',
             'imposed_odor_perm_dff', 'R1', 'rung', 'rung_glom_perm', 'imposed_route2_repeat')


def _mr(x: pd.Series, fmt: str = '.2f') -> str:
    """The mean of `x` with its range, or the one value."""
    x = pd.Series(x).dropna()
    if not len(x):
        return '–'
    if len(x) == 1 or x.min() == x.max():
        return f'{x.iloc[0]:{fmt}}'
    return f'{x.mean():{fmt}} ({x.min():{fmt}} to {x.max():{fmt}})'


def _arms(rows: pd.DataFrame) -> List[str]:
    have = set(rows.arm)
    return [a for a in ARM_ORDER if a in have]


def stage_table(rows: pd.DataFrame, stage_error: pd.DataFrame) -> List[str]:
    """Strength by odor per arm (rows) and megamat stage (columns), mean and range over the wirings (excluded wirings
    left out), E134's held runs first, with the megamat stage error on strength (four stages) and on three views."""
    arms = _arms(rows)
    k = rows[~rows.excluded.astype(bool)]
    se = stage_error[stage_error.input_panel == 'megamat']
    L = ['Strength by odor in one-fly units (lower is better), the default reading, mean over the wirings (range).',
         '', '| arm | ' + ' | '.join(STAGE_NAMES[st] for st in STAGES) +
         ' | the megamat stage error, strength (four stages) | three views (beside) |',
         '|---|' + '---:|' * (len(STAGES) + 2)]
    ref = 'imposed' if (k.arm == 'imposed').any() else k.arm.iloc[0]
    b, sb = k[k.arm == ref], se[se.arm == ref]
    L.append("| E134's held runs (the family) | " + ' | '.join(_mr(b[b.stage == st].strength_base) for st in STAGES) +
             f" | {_mr(sb.stage_error_strength_base, '.3f')} | {_mr(sb.stage_error_three_view_base, '.3f')} |")
    for a in arms:
        x, y = k[k.arm == a], se[se.arm == a]
        L.append(f'| {ARM_NAMES[a]} | ' + ' | '.join(_mr(x[x.stage == st].strength) for st in STAGES) +
                 f" | {_mr(y.stage_error_strength, '.3f')} | {_mr(y.stage_error_three_view, '.3f')} |")
    return L + ['']


def views_table(rows: pd.DataFrame, arms: Sequence[str] = ('imposed', 'imposed_dff', 'R1', 'rung')) -> List[str]:
    """Sharing, breadth and the three-view population error per stage, E134's held runs beside each arm (means over
    the wirings)."""
    k = rows[~rows.excluded.astype(bool)]
    arms = [a for a in arms if (k.arm == a).any()]
    L = ['Sharing / breadth / the three-view population error, one-fly units, mean over the wirings.', '',
         '| stage | E134\'s held runs | ' + ' | '.join(ARM_NAMES[a] for a in arms) + ' |',
         '|---|' + '---:|' * (len(arms) + 1)]
    for st in STAGES:
        x = k[k.stage == st]
        b = x[x.arm == arms[0]]
        cells = [f'{b.view_sharing_base.mean():.2f} / {b.view_breadth_base.mean():.2f} / {b.error_base.mean():.2f}']
        for a in arms:
            y = x[x.arm == a]
            cells.append(f'{y.view_sharing.mean():.2f} / {y.view_breadth.mean():.2f} / {y.error.mean():.2f}')
        L.append(f'| {STAGE_NAMES[st]} | ' + ' | '.join(cells) + ' |')
    return L + ['']


def apl16_table(scores: pd.DataFrame, e134_step2: Optional[pd.DataFrame]) -> List[str]:
    """The APL16 somata (the input of the APL16 panel), strength by odor and the three-view population error, APL on
    against the control flies and APL off against APL output blocked (TNT), per arm, beside the family's own runs."""
    s = scores[(scores.input_panel == 'apl16') & (scores.stage == 'kc_soma') & (scores.panel == 'apl16')
               & (scores.reading == 'curve_noise') & scores.wiring.isin(E.WIRINGS)]
    L = ['APL16 KC somata, strength by odor / the three-view population error, one-fly units, mean over the wirings.',
         '', '| arm | APL on | APL off |', '|---|---:|---:|']
    if e134_step2 is not None:
        b = e134_step2[(e134_step2.source == 'E90:B') & (e134_step2.stage == 'kc_soma') &
                       (e134_step2.panel == 'apl16') & (e134_step2.reading == 'curve_noise') &
                       e134_step2.wiring.isin(E.WIRINGS)]
        c = {cond: b[b.condition == cond] for cond in ('on', 'off')}
        L.append('| the family (E134\'s step 2, E90:B) | ' + ' | '.join(
            f'{c[x].view_strength.mean():.2f} / {c[x].error.mean():.2f}' for x in ('on', 'off')) + ' |')
    for a in ('R1', 'rung', 'imposed'):
        y = s[s.arm == a]
        if not len(y):
            continue
        L.append(f'| {ARM_NAMES[a]} | ' + ' | '.join(
            f'{y[y.condition == x].view_strength.mean():.2f} / {y[y.condition == x].error.mean():.2f}'
            for x in ('on', 'off')) + ' |')
    return L + ['']


def step1_table(rungs: pd.DataFrame, controls: Optional[pd.DataFrame]) -> List[str]:
    """Step 1's rungs on the held-out flies and odors: the matrix view (what keeps a rung), strength by odor and the
    other views, in both directions, with whether each was kept and in the lock's affine_noise reading."""
    L = ['Step 1, held-out PN dendrites, one-fly units (half a -> half b / half b -> half a).', '',
         '| rung | the matrix view (keeps a rung) | strength by odor | sharing | breadth | kept | kept, affine_noise |',
         '|---|---:|---:|---:|---:|---|---|']
    d = ('half_a->half_b', 'half_b->half_a')
    for _, g in rungs.iterrows():
        f = lambda v: ' / '.join(f"{g[f'{x}_view_{v}']:.2f}" for x in d)
        ka = g.get('kept_affine_noise')
        L.append(f"| {g['rung']} | {f('matrix')} | {f('strength')} | {f('sharing')} | {f('breadth')} | "
                 f"{'yes' if g['kept'] else 'no'} | {'–' if pd.isna(ka) else ('yes' if ka else 'no')} |")
    L.append('')
    if controls is not None and len(controls):
        L += ['The controls (odor labels permuted before the fit), the point read like the rungs:', '',
              '| control | R1 | R2m | R2f |', '|---|---:|---:|---:|']
        for r in controls.itertuples(index=False):
            L.append(f'| {r.perm} | {r.R1:.2f} | {r.R2m:.2f} | {r.R2f:.2f} |')
        L.append('')
    return L


def checks_table(checks: pd.DataFrame, verdict: dict) -> List[str]:
    """Ruling 8 (the imposed tables against the PN rates and the PN matrix error) and ruling 9 (the routes)."""
    L = []
    if len(checks):
        L += ['Ruling 8, per imposed run: the largest difference between the table and the PN rates in E134\'s '
              'reading, and the PN-dendrite matrix error (one-fly units; the flies\' floor about 0.31).', '',
              '| panel, arm | runs | largest difference (Hz) | PN matrix error | entries clipped |',
              '|---|---:|---:|---:|---:|']
        for (p, a), x in checks.groupby(['panel', 'arm'], sort=False):
            L.append(f'| {p}, {a} | {len(x)} | {x.imposed_check_max_hz.max():.3f} | {_mr(x.pn_matrix_error)} | '
                     f'{int(x.clip_n_clipped.sum()) if "clip_n_clipped" in x else 0} |')
        L.append('')
    ro = verdict.get('routes') or {}
    if ro.get('available'):
        L += ['Ruling 9, the real wiring, strength by odor (one-fly units) by route:', '',
              '| stage | route 1 | route 2 | difference | threshold | agree |', '|---|---:|---:|---:|---:|---|']
        for st, x in ro['stages'].items():
            L.append(f"| {STAGE_NAMES[st]} | {x['route1']:.3f} | {x['route2']:.3f} | {x['diff']:.3f} | "
                     f"{x['threshold']:.3f} | {'yes' if x['agree'] else 'no'} |")
        sv = ro.get('solve') or {}
        L += ['', f"Route 2's solve: converged {sv.get('solve_converged')}, largest miss "
                  f"{sv.get('solve_max_abs_diff_hz', np.nan):.1f} Hz, the glomeruli without input drift up to "
                  f"{sv.get('solve_no_input_max_drift_hz', np.nan):.1f} Hz.", '']
    return L


def scale_table(fits_rec: dict) -> List[str]:
    """Ruling 5's scale factors per fit (the fold that held out the named odors, and the fit on all odors) and rung,
    those outside 0.5 to 2 named."""
    sc = fits_rec.get('scale') or {}
    fac = sc.get('factors') or {}
    if not fac:
        return []
    af = fits_rec.get('all_flies') or {}
    keys = list(next(iter(fac.values())))
    L = ['Ruling 5, the rate-scale factor of each fit on all seven PN flies (outside 0.5 to 2 named).', '',
         '| fit (odors held out) | ' + ' | '.join(fac) + ' |', '|---|' + '---:|' * len(fac)]
    for k in keys:
        held = ', '.join(af.get(k, {}).get('held_out', [])) or 'none (all 17 odors)'
        L.append(f'| {k} ({held}) | ' + ' | '.join(f'{fac[r][k]:.3g}' for r in fac) + ' |')
    named = sc.get('named') or {}
    return L + ['', f"Named: {', '.join(f'{r}: folds {list(v)}' for r, v in named.items())}.", '']


def as_ruled_table(fits_rec: dict, wiring: str = 'real') -> List[str]:
    """The as-ruled tables of R1 and the kept rung on the megamat input (each megamat odor from the fold fit that left
    it out, at its scale factor; what step 2 ran) against the same rung's fit on all 17 odors: per odor, the as-ruled
    total over glomeruli as a multiple of the all-odor fit's, and Pearson r over all entries. Builds the family's PN
    stage on `wiring` (a model build, no simulation of the KCs)."""
    st, _ = E.family_stage('megamat', wiring)
    fits = fits_rec['all_flies']
    held = {o: k for k, f in fits.items() if k != 'all' for o in f['held_out']}
    L = ['As ruled (each odor from the fold that left it out, rescaled) against the fit on all 17 odors, megamat '
         f'input, {wiring} wiring, before clipping.', '',
         '| rung | r over entries | odor totals, as ruled / all-odor fit (range) | the held-out odors of folds 0, 2, 3 |',
         '|---|---:|---:|---|']
    for rung in ('R1', E.step2_rung_of(fits_rec) or 'R3'):
        r = {'R2m': 'R2', 'R2f': 'R2'}.get(rung, rung)
        a = E.rung_table(st, fits, r)
        b = E.rung_table(st, fits, r, fold_odors=())
        a, b = a.loc[:, b.columns], b
        rr = float(np.corrcoef(a.to_numpy().ravel(), b.to_numpy().ravel())[0, 1])
        tot = a.sum() / b.sum()
        sel = [o for o in tot.index if held.get(o) in ('0', '2', '3')]
        L.append(f"| {rung} | {rr:.3f} | {tot.min():.2f} to {tot.max():.2f} | "
                 f"{', '.join(f'{o} {tot[o]:.2f}' for o in sel)} |")
    return L + ['']


def data_targets() -> Dict[str, float]:
    """The data-target registry's values (results/data_targets/data_targets.json, `kcq.data_targets`), by key; empty if
    it is not built."""
    p = E.out_dir().parent / 'data_targets' / 'data_targets.json'
    if not p.exists():
        return {}
    d = json.loads(p.read_text())
    d = d.get('values', d)
    return {k: v['value'] for k, v in d.items() if isinstance(v, dict) and isinstance(v.get('value'), (int, float))}


def realism_table(realism: pd.DataFrame) -> List[str]:
    """The APL16 runs' realism per arm (range over the wirings): response rates with APL on and off, their ratio,
    the response rate ρ with APL on and off, odors per KC (total variation distance), and flags."""
    a = realism[realism.input_panel == 'apl16']
    L = ['APL16 runs, range over the eight wirings (bands: response rate with APL on 0.045 to 0.085; ratio 1.6 to 2.9; '
         'response rate ρ with APL on at least 0.37; odors per KC 0.155 with APL on, 0.223 with APL off).', '',
         '| arm | response rate, APL on | APL off | ratio | response rate ρ, APL on | APL off | odors per KC, APL on '
         '| APL off | flagged |', '|---|---:|---:|---:|---:|---:|---:|---:|---|']
    rng = lambda x, f='.3f': f'{x.min():{f}} to {x.max():{f}}'
    t = data_targets()
    if t:
        g = lambda k, f='.3f': f'{t[k]:{f}}' if k in t else '–'
        L.append(f"| the KC data (control flies / APL output blocked (TNT), noise criterion) | {g('noise_ctrl_frac')} | "
                 f"{g('noise_tnt_frac')} | {g('noise_rate_ratio_fraction', '.2f')} | – | – | band "
                 f"{g('noise_ctrl_odors_per_kc_tv_max')} | band {g('noise_tnt_odors_per_kc_tv_max')} | – |")
    for arm, x in a.groupby('arm', sort=False):
        flag = x.realism_flag.fillna('').astype(str)
        L.append(f'| {ARM_NAMES.get(arm, arm)} | {rng(x.rate_on)} | {rng(x.rate_off)} | {rng(x.rate_ratio, ".2f")} | '
                 f'{rng(x.rate_rho_on, ".2f")} | {rng(x.rate_rho_off, ".2f")} | {rng(x.opk_on_tv)} | '
                 f'{rng(x.opk_off_tv)} | {int((flag != "").sum())} of {len(x)} |')
    return L + ['']


def family_realism(wirings: Sequence[str] = E.WIRINGS) -> pd.DataFrame:
    """The family's own APL16 runs (its saved fit on each wiring, `family_saved`) through the same bands
    (`run_bands`), as rows of realism_table's columns (arm 'family')."""
    rows = []
    for w in wirings:
        b = E.run_bands(E.family_saved(w))
        rows.append(dict(input_panel='apl16', arm='family', wiring=w, rate_on=b.get('rate_on'),
                         rate_off=b.get('rate_off'), rate_ratio=b.get('ratio', b.get('rate_ratio')),
                         rate_rho_on=b.get('rate_rho_on'), rate_rho_off=b.get('rate_rho_off'),
                         opk_on_tv=b.get('opk_on_tv'), opk_off_tv=b.get('opk_off_tv'),
                         realism_flag=b.get('realism_flag') or ''))
    return pd.DataFrame(rows)


def apl_effect_table(effect: pd.DataFrame) -> List[str]:
    """The APL effect on odor-odor distances of the APL16 runs (a held-out test, reported, not gated), through the
    calibrated readout: the loss, the mean APL effect (the flies' +0.18), the pattern part with its amplitude and shape
    parts, ρ and Pearson r on the per-pair APL effect; means over the wirings, the family's own fits beside."""
    cols = [('cal_apl_effect_rmse', 'the APL-effect loss'), ('cal_apl_effect_mean', 'the mean APL effect'),
            ('cal_apl_effect_rmse_centered', 'the pattern part'), ('cal_apl_effect_amplitude', 'amplitude part'),
            ('cal_apl_effect_shape', 'shape part'), ('cal_apl_effect_spearman', 'ρ on the per-pair APL effect'),
            ('cal_apl_effect_pearson', 'Pearson r on the per-pair APL effect')]
    cols = [(c, n) for c, n in cols if c in effect.columns]
    L = ['The APL effect on the APL16 panel, calibrated readout, mean over the eight wirings.', '',
         '| arm | ' + ' | '.join(n for _, n in cols) + ' |', '|---|' + '---:|' * len(cols)]
    for arm in ('family', 'R1', 'rung', 'imposed'):
        x = effect[effect.arm == arm]
        if len(x):
            L.append(f'| {ARM_NAMES.get(arm, "the family (its own fits)")} | ' +
                     ' | '.join(f'{x[c].mean():+.3f}' if 'mean' in c else f'{x[c].mean():.3f}' for c, _ in cols) + ' |')
    return L + ['']


def values_table(rows: pd.DataFrame, stage_error: pd.DataFrame, scores: pd.DataFrame,
                 e134_step2: Optional[pd.DataFrame], realism: pd.DataFrame, effect: Optional[pd.DataFrame],
                 rungs: pd.DataFrame) -> pd.DataFrame:
    """The numbers of the tables above as one long table, for the summary's placeholders
    (`E135_antennal_lobe_dendrites/writeup_values.csv:arm=...;stage=...;statistic=...:value`): columns arm, stage,
    panel, condition, statistic, value. Arms are the step-2 arms, 'family_held' (E134's held runs), 'family' (the
    family's own APL16 runs) and step 1's rungs (stage pn_dendrite, statistics per direction)."""
    out = []

    def add(arm, stage, statistic, value, panel='megamat', condition='on'):
        if value is not None and not (isinstance(value, float) and np.isnan(value)):
            out.append(dict(arm=arm, stage=stage, panel=panel, condition=condition, statistic=statistic,
                            value=float(value)))

    k = rows[~rows.excluded.astype(bool)]
    ref = k[k.arm == 'imposed'] if (k.arm == 'imposed').any() else k
    for st in STAGES:
        b = ref[ref.stage == st]
        for col, name in (('strength_base', 'strength'), ('view_sharing_base', 'sharing'),
                          ('view_breadth_base', 'breadth'), ('error_base', 'error')):
            add('family_held', st, f'{name}_mean', b[col].mean())
        add('family_held', st, 'strength_min', b.strength_base.min())
        add('family_held', st, 'strength_max', b.strength_base.max())
    for (arm, st), x in k.groupby(['arm', 'stage']):
        add(arm, st, 'strength_mean', x.strength.mean())
        add(arm, st, 'strength_min', x.strength.min())
        add(arm, st, 'strength_max', x.strength.max())
        add(arm, st, 'strength_ratio_mean', x.strength_ratio.mean())
        add(arm, st, 'strength_ratio_min', x.strength_ratio.min())
        add(arm, st, 'strength_ratio_max', x.strength_ratio.max())
        for col in ('view_sharing', 'view_breadth', 'error'):
            add(arm, st, f'{col.replace("view_", "")}_mean', x[col].mean())
    se = stage_error[stage_error.input_panel == 'megamat']
    rs = se[se.arm == 'imposed'] if (se.arm == 'imposed').any() else se
    add('family_held', 'four_stages', 'stage_error_strength_mean', rs.stage_error_strength_base.mean())
    add('family_held', 'four_stages', 'stage_error_three_view_mean', rs.stage_error_three_view_base.mean())
    for arm, x in se.groupby('arm'):
        for col in ('stage_error_strength', 'stage_error_three_view'):
            add(arm, 'four_stages', f'{col}_mean', x[col].mean())
            add(arm, 'four_stages', f'{col}_min', x[col].min())
            add(arm, 'four_stages', f'{col}_max', x[col].max())
    s = scores[(scores.input_panel == 'apl16') & (scores.stage == 'kc_soma') & (scores.panel == 'apl16')
               & (scores.reading == 'curve_noise') & scores.wiring.isin(E.WIRINGS)]
    for (arm, cond), x in s.groupby(['arm', 'condition']):
        add(arm, 'kc_soma', 'strength_mean', x.view_strength.mean(), 'apl16', cond)
        add(arm, 'kc_soma', 'error_mean', x.error.mean(), 'apl16', cond)
    if e134_step2 is not None:
        b = e134_step2[(e134_step2.source == 'E90:B') & (e134_step2.stage == 'kc_soma') &
                       (e134_step2.panel == 'apl16') & (e134_step2.reading == 'curve_noise') &
                       e134_step2.wiring.isin(E.WIRINGS)]
        for cond, x in b.groupby('condition'):
            add('family', 'kc_soma', 'strength_mean', x.view_strength.mean(), 'apl16', cond)
            add('family', 'kc_soma', 'error_mean', x.error.mean(), 'apl16', cond)
    m = realism[realism.input_panel == 'megamat'] if 'input_panel' in realism else realism.iloc[:0]
    for arm, x in m.groupby('arm'):
        add(arm, 'pn_dendrite', 'totals_cv_mean', x.totals_cv_pn.mean())
        add(arm, 'pn_bouton', 'totals_cv_mean', x.totals_cv_bouton.mean())
    a = realism[realism.input_panel == 'apl16'] if 'input_panel' in realism else realism.iloc[:0]
    for arm, x in a.groupby('arm'):
        for col in ('rate_on', 'rate_off', 'rate_ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv'):
            if col in x:
                add(arm, 'kc_soma', f'{col}_min', x[col].min(), 'apl16', 'realism')
                add(arm, 'kc_soma', f'{col}_max', x[col].max(), 'apl16', 'realism')
    if effect is not None:
        for arm, x in effect.groupby('arm'):
            for col in ('cal_apl_effect_rmse', 'cal_apl_effect_mean', 'cal_apl_effect_rmse_centered',
                        'cal_apl_effect_amplitude', 'cal_apl_effect_shape', 'cal_apl_effect_spearman',
                        'cal_apl_effect_pearson'):
                if col in x:
                    add(arm, 'kc_soma', f'{col}_mean', x[col].mean(), 'apl16', 'apl_effect')
            if 'cal_apl_effect_sd_ratio' in x:
                add(arm, 'kc_soma', 'cal_apl_effect_sd_ratio_min', x.cal_apl_effect_sd_ratio.min(), 'apl16', 'apl_effect')
                add(arm, 'kc_soma', 'cal_apl_effect_sd_ratio_max', x.cal_apl_effect_sd_ratio.max(), 'apl16', 'apl_effect')
    for _, g in rungs.iterrows():
        for dname, dkey in (('ab', 'half_a->half_b'), ('ba', 'half_b->half_a')):
            for v in ('matrix', 'strength', 'sharing', 'breadth'):
                add(str(g['rung']), 'pn_dendrite', f'{v}_{dname}', g[f'{dkey}_view_{v}'], 'megamat', 'step1')
    return pd.DataFrame(out)


def main(root: Optional[Path] = None, with_tables: bool = True) -> Path:
    """Writes writeup_tables.md; `with_tables` False skips the as-ruled comparison (which builds the family's PN stage)."""
    d = E.out_dir() if root is None else Path(root)
    v = json.loads((d / E.SCORE_VERDICT_FILE).read_text())
    rows = pd.read_csv(d / 'score_class_rows.csv')
    se = pd.read_csv(d / 'score_stage_error.csv')
    scores = pd.read_csv(d / 'score_scores.csv')
    try:
        e134 = E.read_e134_table(E.E134_STEP2_SCORES)
    except FileNotFoundError:
        e134 = None
    checks = pd.read_csv(d / 'score_imposed_checks.csv')
    realism = pd.read_csv(d / 'score_realism.csv')
    rungs = pd.read_csv(d / 'step1_rungs.csv')
    controls = pd.read_csv(d / 'step1_controls.csv') if (d / 'step1_controls.csv').exists() else None
    oc = v.get('outcome') or {}
    L = ['# E135 write-up tables', '',
         f"From score_verdict.json of {v.get('time')} ({v.get('host')}): outcome {oc.get('outcome')!r}, class "
         f"{v.get('cls')!r}, second reading {oc.get('cls_second')!r}, dF/F repeat {oc.get('cls_dff', '–')!r}; "
         f"{v.get('n_runs')} runs ({v.get('n_planned')} planned, and the dF/F repeat's 24 when asked); dF/F repeat complete "
         f"{(v.get('dff_repeat') or {}).get('complete')}.", '',
         '## The stage table', ''] + stage_table(rows, se)
    L += ['## The views beside', ''] + views_table(rows)
    L += ['## The APL16 somata', ''] + apl16_table(scores, e134)
    L += ['## Step 1', ''] + step1_table(rungs, controls)
    fits_rec = E.read_step1_fits(d)
    L += ['## The scale factors and the as-ruled tables', ''] + scale_table(fits_rec)
    if with_tables:
        L += as_ruled_table(fits_rec)
    L += ['## Checks', ''] + checks_table(checks, v)
    if with_tables:
        realism = pd.concat([family_realism(), realism], ignore_index=True)
    L += ['## Realism', ''] + realism_table(realism)
    if (d / 'score_apl_effect.csv').exists():
        L += ['## The APL effect (held out)', ''] + apl_effect_table(pd.read_csv(d / 'score_apl_effect.csv'))
    eff = pd.read_csv(d / 'score_apl_effect.csv') if (d / 'score_apl_effect.csv').exists() else None
    values_table(rows, se, scores, e134, realism, eff, rungs).to_csv(d / 'writeup_values.csv', index=False)
    out = d / 'writeup_tables.md'
    out.write_text('\n'.join(L) + '\n')
    print(out)
    return out


if __name__ == '__main__':
    import sys
    main(with_tables='--no-tables' not in sys.argv)
