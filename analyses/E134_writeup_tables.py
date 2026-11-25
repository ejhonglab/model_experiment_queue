"""Write-up tables for E134 (the models as populations at every recorded stage), from E134's results folder, the
lock, the APL-effect table and the data-target registry; no simulation.

The numbers E134's note gives that the score's own tables (score_table.md, score_classes.json, flies_table.md,
measure_gate.json) do not give directly, each checked against the score's and the gate's own files (a difference
raises):
- `stage_table`: per stage, panel and condition, the flies (their number, floor, one fly against the others and the
  most atypical fly, for the population error and for strength by odor) beside the family's eight wirings (their mean,
  range and the error common to the wirings: the population error, with the matrix view where the units are named,
  and each view), as fitted and held, in each stage's default reading and its second; the classes read by the entry's
  rule (inside: no further from the flies' mean than the most atypical fly on 6 of 8 wirings; far: above 2 on 6 of 8;
  not readable: the two readings' classes differ) on the three-view error, checked against score_classes.json where
  the score classes the stage, and the same rule on strength by odor, the measure adopted on 2026-11-24 (with the
  flies' most atypical fly on strength as the line);
- `stage_error`: the megamat stage error per wiring, the mean over the eight wirings and the error of the wirings'
  mean views, on the three views as scored (the PN dendrites on four), on strength by odor alone (the adopted measure),
  on sharing alone, on breadth alone and on the two views strength and sharing, for the family as fitted and held and
  for the two beside models that saved every stage; checked against score_stage_error.csv, score_classes.json and the
  gate's read-only read in measure_gate.json;
- `wiring_se`: three standard errors over the eight wirings of each megamat stage's strength-by-odor and three-view
  error in the held runs, the thresholds E135 and E136 froze from this table;
- `beside`: the five beside models at every stage they saved, on the megamat, validation2 and mixture inputs, with
  their megamat stage error where all four stages exist;
- `steps`: per megamat step (as fitted and held), its class as score_classes.json gives it, re-read here from the step
  table's 'mean_distance' rows against its thresholds; per view, inside on the wirings' mean change, and the wirings
  inside under the corrected 'sets' null and under 'one_fly'; the no-step yardstick beside the model's downstream error
  on the step's odors, on the population error and on strength by odor (the PN dendrites on the 15 odors they share with
  the ORNs re-scored from the saved populations, outputs/E134_population_baseline/populations/score/, after checking
  that the same code gives the score's own numbers on all 17; the step is left without the strength column where those
  local files are missing); the flies' own recorded change beside;
- `panels`: the KC somata of each panel, the family as fitted and held, the population error and strength by odor beside
  the APL16 panel's value (step 2, the family as fitted, APL on) and the ORN recordings taken as the prediction; the
  entry's panel class on the population error (checked against score_classes.json) and the same rule read on strength;
- `consistency`: the archive's split-half ρ of every measure, version, reading, condition and variant set
  (archive_consistency_<reading>.csv), the gate's numbers checked against measure_gate.json and score_classes.json;
- `breadth`: the spread of the archived variants' view errors inside the bands at the APL16 somata beside the flies'
  floor, both readings (why breadth's order does not hold between halves of the recordings);
- `tradeoff`: design step 7, the best variant and family by the population error with APL on and off and by the
  pattern part, with each measure's ρ over the variants with the calibrated pattern part (checked against
  archive_tradeoff_best_<reading>.csv), and the same for strength by odor as a description;
- `soma_criteria`: the KC somata's response criterion per set of recordings, the footing factors, the noise bridges and
  the response rates, from the lock (checked against its stages block and the per-fly rates);
- `apl16_models`: the planner's probe of 2026-11-24 09:11 (raw counts, strength by its rate part only) reproduced from
  the committed parts, beside the committed population error (both parts of strength) on raw counts and on the default
  reading, for the family (step 2's E90 part B) and the probe's other rows (the archive's variants);
- `realism_apl16`: the four bands for the family's eight wirings and the three beside models that have an APL16 run,
  from results/apl_effect/by_variant.csv (checked against score_apl16_bands.csv), beside the KC data's values;
- `operating_points`: per step-3 group the response rates with APL on and off and their ratio, the largest spike
  counts, the threshold and the held factors beside the fit's own, and the flies' response rate (score_realism.csv
  checked against variants.csv);
- `planner_check`: the numbers of E134's readings (the planner's of 11:11, viki's of 12:05, the gate planner's of
  15:03 and measure_gate.json) and of the hand-offs that quote E134's tables (hal's E135 thresholds of 15:53, disco's
  E136 yardsticks of 15:44, hal's two-view table of 16:00), beside the script's value at the precision written, and
  whether it reproduces (reported, not raised).

Writes writeup_*.csv (one per table above), writeup_tables.md (all of them as markdown) and
provenance_writeup_tables.json to the results folder.

Usage (from the repository root; light, no simulation):
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 nice -n 19 ./venv/bin/python analyses/E134_writeup_tables.py [--out DIR]
"""
import argparse
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E134_population_baseline'
SCRIPT = 'analyses/E134_writeup_tables.py'
W8 = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')
MEGAMAT_STAGES = ('orn', 'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
ERROR_STAGES = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
# each stage's default reading and its second (the lock's; `check_readings` checks them against kcq/stage_readings)
READINGS = {'orn': ('curve_noise', 'affine_noise'), 'pn_dendrite': ('curve_noise', 'affine_noise'),
            'pn_bouton': ('curve_noise', 'curve_noise_g3'), 'kc_claw': ('curve_noise', 'curve_noise_b14'),
            'kc_soma': ('curve_noise', 'raw')}
STAGE_NAME = {'orn': 'ORNs', 'pn_dendrite': 'PN dendrites', 'pn_bouton': 'PN boutons', 'kc_claw': 'KC claws',
              'kc_soma': 'KC somata'}
PANEL_NAME = {'megamat': 'megamat', 'apl16': 'APL16 panel', 'validation2': 'validation2', 'kmix': 'kiwi mixtures',
              'cmix': 'control mixtures'}
STEPS = ('orn->pn_dendrite', 'pn_dendrite->pn_bouton', 'pn_bouton->kc_claw', 'kc_claw->kc_soma')
VIEWS = ('strength', 'sharing', 'breadth', 'matrix')
BESIDE = ('published_hemibrain', 'uniform_draw', 'bouton_baseline', 'reference', 'firing_ceiling')
BESIDE_NAME = {'published_hemibrain': "the lab's published hemibrain model", 'uniform_draw': 'a uniform draw',
               'bouton_baseline': 'the bouton baseline', 'reference': 'the current reference model',
               'firing_ceiling': 'the firing-ceiling fit'}
# the beside models' APL16 variants in the archive (the runner's replays: E3's bouton_pnapl_x1, E44's operating point
# that E69's ref40 rebuilds, E69's fit 1); the published hemibrain model and the uniform draw were tuned on megamat
BESIDE_APL16 = {'bouton_baseline': ('E3_drive_dependent_gain', 'bouton_pnapl_x1'),
                'reference': ('E44_operating_point', 'thr40_gain1.5_rate0.068_k2_cap3'),
                'firing_ceiling': ('E69_kc_firing_ceiling', 'rate0.068_aplboth_refr3_adapt0.04142t1000_k2')}
# the planner's probe rows (experiment_priorities.md, E134's entry): label, archive variants, written values (APL on,
# APL off; None where not written)
PROBE_ROWS = (
    ('the firing-ceiling fit', [('E69_kc_firing_ceiling', 'rate0.068_aplboth_refr3_adapt0.04142t1000_k2')], '1.02', '0.99'),
    ('the first joint fit', [('E10_joint_fit', 'k2')], '1.15', '1.49'),
    ('the current reference model', [('E44_operating_point', 'thr40_gain1.5_rate0.068_k2_cap3')], '1.21', '1.19'),
    ('the bouton baseline', [('E3_drive_dependent_gain', 'bouton_pnapl_x1')], '2.08', '1.76'),
    ('the per-claw threshold at 70 Hz', [('E16_claw_saturation', 'bouton_thr70')], '2.32', '1.87'),
    ('random uniform wiring, eight draws', [('E4_substrate', f'uniform_sp0.05_seed{s}') for s in range(10, 18)],
     '2.54 (2.44 to 2.62)', None),
)
INSIDE_MIN = 6
FAR = 2.0
PANEL_WITHIN = 0.2
TOL = 1e-12


def _sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _mean(x) -> float:
    x = [float(v) for v in x]
    return math.fsum(x) / len(x)


def _rms(x) -> float:
    x = np.asarray(list(x), dtype=float)
    return float(np.sqrt(np.mean(np.square(x))))


class Checks:
    """Every comparison with a saved number; `eq` raises on a difference."""

    def __init__(self):
        self.n = 0

    def eq(self, what: str, a, b, tol: float = TOL) -> None:
        self.n += 1
        if isinstance(a, str) or isinstance(b, str) or isinstance(a, (bool, np.bool_)):
            if a != b:
                raise AssertionError(f'{what}: {a!r} != {b!r}')
            return
        a, b = float(a), float(b)
        if not (np.isfinite(a) and np.isfinite(b)) or abs(a - b) > tol * max(1.0, abs(a), abs(b)):
            raise AssertionError(f'{what}: {a!r} != {b!r}')

    def rounds(self, what: str, value: float, written: float, decimals: int) -> None:
        """A number written at `decimals` in a file that is itself a reading (measure_gate.json)."""
        self.n += 1
        if round(float(value), decimals) != round(float(written), decimals):
            raise AssertionError(f'{what}: {value!r} does not round to {written!r}')


# ---------------------------------------------------------------------------
# inputs

def load(res: Path) -> dict:
    j = lambda n: json.loads((res / n).read_text())
    c = lambda n: pd.read_csv(res / n, low_memory=False)
    reg = json.loads((H.shared_dir() / 'results' / 'data_targets' / 'data_targets.json').read_text())
    by_variant = pd.read_csv(H.shared_dir() / 'results' / 'apl_effect' / 'by_variant.csv', low_memory=False)
    I = dict(sc=c('score_scores.csv'), cls=j('score_classes.json'), stp=c('score_step_table.csv'),
             se=c('score_stage_error.csv'), real=c('score_realism.csv'), bands=c('score_apl16_bands.csv'),
             s2=c('step2_scores.csv'), g2=c('step2_gate_wirings.csv'), gate=j('measure_gate.json'),
             fv=c('flies_views.csv'), floo=c('flies_loo.csv'), fns=c('flies_no_step.csv'),
             fmc=c('flies_megamat_vs_control.csv'), fss=c('flies_steps_summary.csv'), lock=j('lock.json'),
             var=c('variants.csv'), reg=reg, bv=by_variant)
    for rd in ('curve_noise', 'raw'):
        I[f'cons_{rd}'] = c(f'archive_consistency_{rd}.csv')
        I[f'arch_{rd}'] = c(f'archive_scores_{rd}.csv')
        I[f'best_{rd}'] = c(f'archive_tradeoff_best_{rd}.csv')
    return I


def check_readings(ck: Checks) -> None:
    """The stages' default and second readings here are kcq/stage_readings' (which the lock records)."""
    from kcq import stage_readings as R
    for st, (d, s) in READINGS.items():
        rs = R.readings_of(st)
        ck.eq(f'{st} default reading', rs[0].name, d)
        ck.eq(f'{st} second reading', rs[1].name, s)


def cohort(panel: str, condition: str) -> str:
    if panel == 'apl16':
        return 'control' if condition == 'on' else 'tnt'
    return 'wild_type'


def pop_name(stage: str, panel: str, condition: str) -> str:
    return f'{stage}__{panel}__{cohort(panel, condition)}'


def fly_lines(I: dict, ck: Checks, stage: str, panel: str, condition: str) -> dict:
    """The flies' floor, one fly against the others and the most atypical fly, for the population error, its matrix
    version and each view, from step 1 (flies_views.csv; the maxima checked against flies_loo.csv)."""
    name = pop_name(stage, panel, condition)
    fv = I['fv'][I['fv']['population'] == name]
    lo = I['floo'][I['floo']['population'] == name]
    out = dict(population=name, n_flies=int(fv['n_flies'].iloc[0]), n_odors=int(fv['n_odors'].iloc[0]))
    for nm, key in (('error', 'error'), ('error_with_matrix', 'error_with_matrix'), ('strength', 'view_strength'),
                    ('sharing', 'view_sharing'), ('breadth', 'view_breadth'), ('matrix', 'view_matrix')):
        r = fv[fv['name'] == nm]
        if not len(r) or not np.isfinite(float(r['loo_max'].iloc[0])):
            continue
        r = r.iloc[0]
        out[f'{nm}_floor'] = float(r['floor'])
        out[f'{nm}_loo_mean'] = float(r['loo_mean'])
        out[f'{nm}_loo_max'] = float(r['loo_max'])
        ck.eq(f'{name} {nm} loo max vs flies_loo', out[f'{nm}_loo_max'], float(lo[key].max()), 1e-9)
    return out


# ---------------------------------------------------------------------------
# the stage table

def _rows(df: pd.DataFrame, **kw) -> pd.DataFrame:
    m = np.ones(len(df), dtype=bool)
    for k, v in kw.items():
        m &= (df[k] == v).to_numpy()
    return df[m]


def family_rows(I: dict, stage: str, panel: str, condition: str, group: str, reading: str) -> pd.DataFrame:
    """The family's rows (eight wirings, 'mean', 'common') of a group in a reading: step 3's and E107's from the
    score, the APL16 panel's from step 2 (E90's part B)."""
    if group == 'apl16':
        x = _rows(I['s2'], source='E90:B', stage=stage, panel=panel, condition=condition, reading=reading)
    else:
        x = _rows(I['sc'], group=group, stage=stage, panel=panel, condition=condition, reading=reading)
    return x.set_index('wiring')


def classify(vals: Sequence[float], line: float) -> str:
    vals = np.asarray(vals, dtype=float)
    if int((vals <= line).sum()) >= INSIDE_MIN:
        return 'inside'
    if int((vals > FAR).sum()) >= INSIDE_MIN:
        return 'far'
    return 'between'


def combine(by_reading: Dict[str, str]) -> str:
    v = list(by_reading.values())
    if len(v) == 2 and len(set(v)) == 1:
        return v[0]
    return 'not readable' if len(set(v)) > 1 else f'{v[0]} (one reading)'


STAGE_GROUPS = (
    [('megamat_fitted', st, 'megamat', 'on') for st in MEGAMAT_STAGES]
    + [('megamat_held', st, 'megamat', 'on') for st in MEGAMAT_STAGES]
    + [('apl16', 'kc_soma', 'apl16', 'on'), ('apl16', 'kc_soma', 'apl16', 'off')]
    + [(g, 'kc_soma', 'validation2', 'on') for g in ('validation2_fitted', 'validation2_held')]
    + [(g, 'kc_soma', p, 'on') for g in ('mixtures_fitted', 'mixtures_held') for p in ('kmix', 'cmix')]
)


def stage_table(I: dict, ck: Checks):
    rows, cls_out = [], {}
    for group, stage, panel, cond in STAGE_GROUPS:
        fl = fly_lines(I, ck, stage, panel, cond)
        rec = dict(group=group, stage=stage, panel=panel, condition=cond, **fl)
        by3, bys, by4 = {}, {}, {}
        for k, rd in enumerate(READINGS[stage]):
            x = family_rows(I, stage, panel, cond, group, rd)
            if not len(x):
                continue
            w = x.loc[list(W8)]
            tag = 'default' if k == 0 else 'second'
            rec[f'{tag}_reading'] = rd
            for col in ('error', 'error_with_matrix', 'view_strength', 'view_sharing', 'view_breadth', 'view_matrix',
                        'part_strength_rate', 'part_strength_mean', 'part_sharing_who', 'part_sharing_how_strongly',
                        'energy_ratio'):
                v = w[col].astype(float)
                if not v.notna().any():
                    continue
                rec[f'{tag}_{col}_mean'] = float(x.loc['mean', col])
                ck.eq(f'{group} {stage} {panel} {cond} {rd} {col} mean of wirings', rec[f'{tag}_{col}_mean'],
                      _mean(v), 1e-9)
                rec[f'{tag}_{col}_min'], rec[f'{tag}_{col}_max'] = float(v.min()), float(v.max())
                rec[f'{tag}_{col}_common'] = float(x.loc['common', col])
            # the score's own lines are the step-1 flies'
            ck.eq(f'{group} {stage} {panel} {cond} {rd} the flies\' most atypical fly', float(x['flies_loo_max'].iloc[0]),
                  fl['error_loo_max'], 1e-9)
            ck.eq(f'{group} {stage} {panel} {cond} {rd} flies', int(x['n_flies'].iloc[0]), fl['n_flies'])
            n_in3 = int((w['error'].astype(float) <= fl['error_loo_max']).sum())
            n_ins = int((w['view_strength'].astype(float) <= fl['strength_loo_max']).sum())
            rec[f'{tag}_n_inside_three'], rec[f'{tag}_n_inside_strength'] = n_in3, n_ins
            by3[rd] = classify(w['error'], fl['error_loo_max'])
            bys[rd] = classify(w['view_strength'], fl['strength_loo_max'])
            rec[f'{tag}_class_three'], rec[f'{tag}_class_strength'] = by3[rd], bys[rd]
            # named glomeruli: the class on four views too (the lock's caution: the write-up says whether the PN
            # dendrites' class differs with and without the matrix view)
            if 'error_with_matrix_loo_max' in fl and w['error_with_matrix'].notna().all():
                rec[f'{tag}_n_inside_four'] = int((w['error_with_matrix'].astype(float)
                                                   <= fl['error_with_matrix_loo_max']).sum())
                by4[rd] = classify(w['error_with_matrix'], fl['error_with_matrix_loo_max'])
                rec[f'{tag}_class_four'] = by4[rd]
        rec['class_three'], rec['class_strength'] = combine(by3), combine(bys)
        if by4:
            rec['class_four'] = combine(by4)
        # the score's classes (the megamat stages)
        sc = I['cls']['stages'].get(group, {}).get(stage)
        if sc is not None:
            ck.eq(f'{group} {stage} class vs score_classes.json', rec['class_three'], sc['cls'])
            for rd, c in sc['by_reading'].items():
                ck.eq(f'{group} {stage} {rd} class vs score_classes.json', by3[rd], c)
        rows.append(rec)
        cls_out[(group, stage, panel, cond)] = (by3, bys)
    t = pd.DataFrame(rows)
    # the largest strength value of the family at the megamat stages of the stage error (E136's planner's reading)
    return t


# ---------------------------------------------------------------------------
# the megamat stage error

STAGE_ERROR_GROUPS = ('megamat_fitted', 'megamat_held', 'firing_ceiling_megamat', 'bouton_baseline_megamat')


def stage_error(I: dict, ck: Checks) -> pd.DataFrame:
    sc, se = I['sc'], I['se']
    rows = []
    for g in STAGE_ERROR_GROUPS:
        x = _rows(sc, group=g, panel='megamat', condition='on')
        ws = [w for w in dict.fromkeys(x['wiring']) if w not in ('mean', 'common')]
        for w in ws + (['common'] if len(ws) > 1 else []):
            rec = dict(group=g, wiring=w)
            per = {k: [] for k in ('three', 'three_flat', 'strength', 'sharing', 'breadth', 'two')}
            for st in ERROR_STAGES:
                r = _rows(x, stage=st, wiring=w, reading=READINGS[st][0])
                if len(r) != 1:
                    raise AssertionError(f'{g} {st} {w}: {len(r)} rows')
                r = r.iloc[0]
                e3 = float(r['error'])
                e4 = float(r['error_with_matrix']) if st == 'pn_dendrite' else e3
                s, h, b = float(r['view_strength']), float(r['view_sharing']), float(r['view_breadth'])
                per['three'].append(e4)
                per['three_flat'].append(e3)
                per['strength'].append(s)
                per['sharing'].append(h)
                per['breadth'].append(b)
                per['two'].append(_rms([s, h]))
                rec[f'{st}_strength'] = s
            for k, v in per.items():
                rec[f'stage_error_{k}'] = _rms(v)
            rows.append(rec)
        df = pd.DataFrame([r for r in rows if r['group'] == g and r['wiring'] != 'common'])
        m = {k: _mean(df[k]) for k in df.columns if k.startswith('stage_error_') or k.endswith('_strength')}
        rows.append(dict(group=g, wiring='model', **m))
        for k in ('stage_error_strength', 'stage_error_three'):
            rows[-1][f'{k}_min'], rows[-1][f'{k}_max'] = float(df[k].min()), float(df[k].max())
    out = pd.DataFrame(rows)
    # the score's own stage error (three views, the PN dendrites on four; and on three)
    for _, r in out.iterrows():
        s = _rows(se, group=r['group'], wiring=r['wiring'])
        if r['wiring'] == 'model' and r['group'] not in ('megamat_fitted', 'megamat_held'):
            s = _rows(se, group=r['group'], wiring='model')
        if not len(s):
            continue
        ck.eq(f"{r['group']} {r['wiring']} stage error vs score_stage_error.csv", r['stage_error_three'],
              float(s['stage_error'].iloc[0]), 1e-9)
        ck.eq(f"{r['group']} {r['wiring']} stage error (three) vs score_stage_error.csv", r['stage_error_three_flat'],
              float(s['stage_error_three'].iloc[0]), 1e-9)
    for g, v in I['cls']['megamat_stage_error'].items():
        if g in STAGE_ERROR_GROUPS:
            ck.eq(f'{g} stage error vs score_classes.json', float(_rows(out, group=g, wiring='model')
                                                                 ['stage_error_three'].iloc[0]), v, 1e-9)
    cut = I['gate']['family_on_the_cut_measure']
    for g in ('megamat_fitted', 'megamat_held'):
        r = _rows(out, group=g, wiring='model').iloc[0]
        ck.rounds(f'{g} strength-only stage error vs measure_gate.json', r['stage_error_strength'],
                  cut['megamat_stage_error_strength_only'][g], 3)
        ck.rounds(f'{g} stage error vs measure_gate.json', r['stage_error_three'],
                  cut['megamat_stage_error_as_scored'][g], 3)
    return out


def wiring_se(I: dict, ck: Checks) -> pd.DataFrame:
    """Three standard errors over the eight wirings (SD with n - 1, over the square root of 8) of each megamat stage's
    strength-by-odor error and three-view error (the PN dendrites also on four views), in the held runs and the default
    reading: the thresholds E135's step 2 and E136's classes froze from this table."""
    rows = []
    for g in ('megamat_held', 'megamat_fitted'):
        for st in MEGAMAT_STAGES:
            x = family_rows(I, st, 'megamat', 'on', g, READINGS[st][0]).loc[list(W8)]
            rec = dict(group=g, stage=st)
            for col in ('view_strength', 'error', 'error_with_matrix'):
                v = x[col].astype(float)
                if v.notna().all():
                    rec[f'{col}_sd'] = float(v.std(ddof=1))
                    rec[f'{col}_3se'] = 3 * float(v.std(ddof=1)) / math.sqrt(len(v))
            rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the beside models

def beside(I: dict, ck: Checks, se: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for m in BESIDE:
        for inp in ('megamat', 'validation2', 'mixtures'):
            g = f'{m}_{inp}'
            x = _rows(I['sc'], group=g)
            if not len(x):
                continue
            for (st, panel), y in x.groupby(['stage', 'panel'], sort=False):
                for k, rd in enumerate(READINGS[st]):
                    r = _rows(y, reading=rd, wiring='mean')
                    if not len(r):
                        continue
                    r = r.iloc[0]
                    c = _rows(y, reading=rd, wiring='common').iloc[0]
                    ck.eq(f'{g} {st} {panel} {rd}: one wiring, mean = common', float(r['error']), float(c['error']), 1e-12)
                    rows.append(dict(model=m, input=inp, stage=st, panel=panel, reading=rd,
                                     which='default' if k == 0 else 'second', error=float(r['error']),
                                     error_with_matrix=float(r['error_with_matrix']),
                                     strength=float(r['view_strength']), sharing=float(r['view_sharing']),
                                     breadth=float(r['view_breadth']), matrix=float(r['view_matrix']),
                                     flies_loo_max=float(r['flies_loo_max'])))
    out = pd.DataFrame(rows)
    out['stage_error_three'] = np.nan
    out['stage_error_strength'] = np.nan
    for m in BESIDE:
        s = _rows(se, group=f'{m}_megamat', wiring='model')
        if len(s):
            i = (out['model'] == m) & (out['input'] == 'megamat')
            out.loc[i, 'stage_error_three'] = float(s['stage_error_three'].iloc[0])
            out.loc[i, 'stage_error_strength'] = float(s['stage_error_strength'].iloc[0])
    return out


# ---------------------------------------------------------------------------
# the steps

def dendrites_on_orn_odors(I: dict, ck: Checks) -> Optional[pd.DataFrame]:
    """The family's PN dendrites scored on the odors they share with the ORN flies (the step's odors), from the saved
    populations of the score (outputs/E134_population_baseline/populations/score/, local to the host that scored),
    after checking that the same code gives the score's own numbers on all 17 odors and the step table's population
    row on the step's odors. None where those files are missing."""
    d = Path('outputs') / EXPERIMENT / 'populations' / 'score'
    files = {(g, w): d / f'{g}__C2__{w}__pn_dendrite__megamat__on.parquet' for g in ('megamat_fitted', 'megamat_held')
             for w in W8}
    if not all(p.exists() for p in files.values()):
        return None
    from kcq import population_fit as P
    from kcq import stage_readings as R
    from experiments import E134_population_baseline as B
    flies = R.FlyCache()
    fd, fo = flies.get(('pn_dendrite', 'megamat', 'on')), flies.get(('orn', 'megamat', 'on'))
    rows = []
    for g in ('megamat_fitted', 'megamat_held'):
        pops = [P.Population.from_criterion(pd.read_parquet(files[(g, w)]), fd.criterion, name=f'C2/{w}',
                                            named_units=True) for w in W8]
        fp, ref, od = B.fly_reference_for(fd, list(pops[0].values.columns))
        t17 = P.wirings_as_flies(pops, ref)
        x = family_rows(I, 'pn_dendrite', 'megamat', 'on', g, 'curve_noise')
        for w in W8 + ('mean', 'common'):
            k = f'C2/{w}' if w in W8 else w
            for col in ('error', 'error_with_matrix', 'view_strength', 'view_sharing', 'view_breadth'):
                ck.eq(f'{g} dendrites {w} {col}, 17 odors, recomputed', float(t17.loc[k, col]), float(x.loc[w, col]), 1e-9)
        _, have = B.fly_side_read(fo)
        od15 = [o for o in od if o in set(have)]
        ref15 = P.fly_reference(fp, od15, B.breadth_edges(len(od15)))
        t15 = P.wirings_as_flies(pops, ref15)
        pr = _rows(I['stp'], group=g, step='orn->pn_dendrite', level='population')
        ck.eq(f'{g} dendrites on the step\'s odors vs the step table', float(t15.loc['mean', 'error']),
              float(pr['model_downstream_error'].iloc[0]), 1e-9)
        for w in ('mean', 'common'):
            rows.append(dict(group=g, wiring=w, n_odors=len(od15), odors=' '.join(od15),
                             **{c: float(t15.loc[w, c]) for c in ('error', 'error_with_matrix', 'view_strength',
                                                                   'view_sharing', 'view_breadth', 'view_matrix')}))
    return pd.DataFrame(rows)


def steps(I: dict, ck: Checks, st_tab: pd.DataFrame, dend15: Optional[pd.DataFrame]) -> pd.DataFrame:
    cls, stp, fns, fss = I['cls']['steps'], I['stp'], I['fns'], I['fss']
    rows = []
    for g in ('megamat_fitted', 'megamat_held'):
        for step in STEPS:
            c = cls[g][step]
            up, down = step.split('->')
            x = _rows(stp, group=g, step=step)
            pop = _rows(x, level='population').iloc[0]
            mean = _rows(x, level='mean_distance')
            # re-read the class: a part inside when the wirings' mean change is within its threshold
            views_in = {}
            for _, r in mean.iterrows():
                thr = c['thresholds'].get(r['part'])
                ck.eq(f'{g} {step} {r["part"]} distance vs score_classes.json', float(r['distance']),
                      c['views'][r['part'].split('/')[0]]['parts'][r['part']]['distance'], 1e-12)
                v = r['part'].split('/')[0]
                views_in[v] = views_in.get(v, True) and bool(thr is not None and float(r['distance']) <= thr)
            for v, inside in views_in.items():
                ck.eq(f'{g} {step} {v} inside', inside, bool(c['views'][v]['inside']))
            ns, md = float(pop['no_step_error']), float(pop['model_downstream_error'])
            cls_re = ('no better than no step' if md >= ns else 'explained' if all(views_in.values()) else 'in part')
            ck.eq(f'{g} {step} class re-read', cls_re, c['cls'])
            up_name, down_name = f'{up}__megamat__wild_type', f'{down}__megamat__wild_type'
            nsr = fns[(fns['step'] == f'{up_name} -> {down_name}')].set_index('name')
            ck.eq(f'{g} {step} no-step vs flies_no_step.csv', ns, float(nsr.loc['error', 'no_step']), 1e-12)
            rec = dict(group=g, step=step, cls=c['cls'], read_on=c['read_on'], no_step_error=ns,
                       model_downstream_error=md, no_step_strength=float(nsr.loc['view_strength', 'no_step']),
                       no_step_sharing=float(nsr.loc['view_sharing', 'no_step']),
                       no_step_breadth=float(nsr.loc['view_breadth', 'no_step']), n_odors=int(nsr['n_odors'].iloc[0]))
            # the model's downstream views on the step's odors
            if step == 'orn->pn_dendrite':
                if dend15 is not None:
                    d = _rows(dend15, group=g, wiring='mean').iloc[0]
                    rec.update(model_strength=float(d['view_strength']), model_sharing=float(d['view_sharing']),
                               model_breadth=float(d['view_breadth']), model_on='the step\'s 15 odors, recomputed')
            else:
                s = _rows(st_tab, group=g, stage=down, panel='megamat').iloc[0]
                ck.eq(f'{g} {step} downstream error = the stage table\'s (17 odors)', md, s['default_error_mean'], 1e-9)
                rec.update(model_strength=s['default_view_strength_mean'], model_sharing=s['default_view_sharing_mean'],
                           model_breadth=s['default_view_breadth_mean'], model_on='17 odors, the stage table')
            for v in ('strength', 'sharing', 'breadth', 'matrix'):
                if v not in c['views']:
                    continue
                rec[f'{v}_inside_mean_change'] = bool(c['views'][v]['inside'])
                for part, d in c['views'][v]['parts'].items():
                    rec[f'{part}_distance'], rec[f'{part}_threshold'] = d['distance'], d['threshold']
                    f = fss[(fss['step'] == f'{up_name} -> {down_name}') & (fss['part'] == part)]
                    if len(f):
                        rec[f'{part}_flies_rms_change'] = float(f['rms_change'].iloc[0])
                    rec[f'{part}_entries_inside'] = c['entries_inside'].get(part)
                for null in ('sets', 'one_fly'):
                    n = c['wirings_inside'][null][v]
                    rec[f'{v}_{null}_n_inside'] = int(n['n_wirings_inside'])
            rec['one_fly_every_view_6of8'] = all(c['wirings_inside']['one_fly'][v]['n_wirings_inside'] >= INSIDE_MIN
                                                 for v in c['wirings_inside']['one_fly'])
            rows.append(rec)
    ck.eq('step null', I['cls']['step_null']['draw'], 'sets')
    ck.eq('step null corrected', bool(I['cls']['step_null']['correct']), True)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the panels

def panels(I: dict, ck: Checks, st_tab: pd.DataFrame) -> pd.DataFrame:
    apl = _rows(st_tab, group='apl16', condition='on').iloc[0]
    ref3, refs = apl['default_error_mean'], apl['default_view_strength_mean']
    rows = []
    for panel, gf, gh, from_ in (('megamat', 'megamat_fitted', 'megamat_held', 'kc_claw'),
                                 ('validation2', 'validation2_fitted', 'validation2_held', 'orn'),
                                 ('kmix', 'mixtures_fitted', 'mixtures_held', 'orn'),
                                 ('cmix', 'mixtures_fitted', 'mixtures_held', 'orn')):
        f = _rows(st_tab, group=gf, stage='kc_soma', panel=panel).iloc[0]
        h = _rows(st_tab, group=gh, stage='kc_soma', panel=panel).iloc[0]
        e, s = f['default_error_mean'], f['default_view_strength_mean']
        cls3 = 'holds' if abs(e - ref3) <= PANEL_WITHIN else 'panel-specific'
        clss = 'holds' if abs(s - refs) <= PANEL_WITHIN else 'panel-specific'
        sc = I['cls']['panels'][panel]
        ck.eq(f'{panel} panel error vs score_classes.json', e, sc['error'], 1e-12)
        ck.eq(f'{panel} APL16 reference vs score_classes.json', ref3, sc['apl16'], 1e-12)
        ck.eq(f'{panel} panel class vs score_classes.json', cls3, sc['cls'])
        ns = I['fns'][I['fns']['step'] == f'{from_}__{panel}__wild_type -> kc_soma__{panel}__wild_type'].set_index('name')
        rows.append(dict(panel=panel, n_flies=int(f['n_flies']), n_odors=int(f['n_odors']), floor=f['error_floor'],
                         loo_max=f['error_loo_max'], strength_loo_max=f['strength_loo_max'],
                         fitted_error=e, fitted_error_min=f['default_error_min'], fitted_error_max=f['default_error_max'],
                         held_error=h['default_error_mean'], fitted_strength=s, held_strength=h['default_view_strength_mean'],
                         fitted_sharing=f['default_view_sharing_mean'], fitted_breadth=f['default_view_breadth_mean'],
                         apl16_error=ref3, apl16_strength=refs, cls_three=cls3, cls_strength=clss,
                         no_step_from=from_, no_step_error=float(ns.loc['error', 'no_step']),
                         no_step_strength=float(ns.loc['view_strength', 'no_step']),
                         held_rate_note={'megamat': "Remy's four flies' 0.0915", 'validation2': "the control flies' 6.8 % (E107's hold)",
                                         'kmix': 'the mixture flies\' 0.0537', 'cmix': 'the mixture flies\' 0.0537'}[panel],
                         model_noise='none: the mixture somata have no single trials' if panel in ('kmix', 'cmix')
                         else 'the flies\' trial noise'))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the archive: consistency and breadth

def consistency(I: dict, ck: Checks) -> pd.DataFrame:
    t = pd.concat([I[f'cons_{rd}'].assign(reading=rd) for rd in ('curve_noise', 'raw')], ignore_index=True)
    gate, g = I['gate'], I['cls']['gate']
    inb = t[t['variant_set'] == 'inside_bands']

    def val(rd, version, measure, cond):
        x = inb[(inb['reading'] == rd) & (inb['version'] == version) & (inb['measure'] == measure)
                & (inb['condition'] == cond)]
        if len(x) != 1:
            raise AssertionError(f'{rd} {version} {measure} {cond}: {len(x)} rows')
        return float(x['split_half_rho'].iloc[0])
    for cond in ('on', 'off'):
        ck.rounds(f'gate B {cond}', val('curve_noise', 'half', 'error', cond), gate['B_consistency']['population_error'][cond], 3)
        ck.eq(f'score gate {cond}', val('curve_noise', 'half', 'error', cond), g['curve_noise']['split_half_rho'][cond], 1e-12)
        ck.eq(f'score gate raw {cond}', val('raw', 'half', 'error', cond), g['raw']['split_half_rho'][cond], 1e-12)
        for v, d in g['curve_noise']['other_versions'].items():
            ck.eq(f'score gate {v} {cond}', val('curve_noise', v, 'error', cond), d[cond], 1e-12)
        for k, d in gate['B_consistency']['beside'].items():
            rd, version = [s.strip() for s in k.split(',')]
            ck.rounds(f'gate beside {k} {cond}', val(rd, version, 'error', cond), d[cond], 3)
        for v in ('strength', 'sharing', 'breadth'):
            ck.rounds(f'gate view {v} {cond}', val('curve_noise', 'half', f'view_{v}', cond),
                      gate['views_singly'][v][cond], 3)
    ck.eq('gate B n_variants', int(inb[(inb['reading'] == 'curve_noise')]['n_variants'].iloc[0]),
          gate['B_consistency']['n_variants'])
    ck.eq('gate B n_splits', int(inb[(inb['reading'] == 'curve_noise')]['n_splits'].iloc[0]),
          gate['B_consistency']['n_splits'])
    return t


def breadth_spread(I: dict) -> pd.DataFrame:
    rows = []
    for rd in ('curve_noise', 'raw'):
        a = I[f'arch_{rd}']
        a = a[(a['inside_bands'].astype(str).str.lower() == 'true') & (a['duplicate'].astype(str).str.lower() == 'false')]
        for cond in ('on', 'off'):
            x = a[a['condition'] == cond]
            fl = I['fv'][I['fv']['population'] == pop_name('kc_soma', 'apl16', cond)].set_index('name')
            for col, nm in (('error', 'error'), ('view_strength', 'strength'), ('view_sharing', 'sharing'),
                            ('view_breadth', 'breadth')):
                v = x[col].astype(float)
                rows.append(dict(reading=rd, condition=cond, measure=nm, n_variants=int(v.notna().sum()),
                                 q05=float(v.quantile(0.05)), q25=float(v.quantile(0.25)), median=float(v.median()),
                                 q75=float(v.quantile(0.75)), q95=float(v.quantile(0.95)),
                                 iqr=float(v.quantile(0.75) - v.quantile(0.25)),
                                 flies_floor=float(fl.loc[nm, 'floor']), flies_loo_max=float(fl.loc[nm, 'loo_max'])))
    return pd.DataFrame(rows)


def tradeoff(I: dict, ck: Checks) -> pd.DataFrame:
    """Design step 7, over the archived variants inside the bands (archive_tradeoff_<reading>.csv): the best variant
    and family by the population error with APL on and with APL off and by the pattern part (calibrated and raw), as
    the archive scoring wrote them (archive_tradeoff_best_<reading>.csv), each measure's Spearman ρ over the variants
    with the calibrated pattern part re-computed here and checked; and beside them, as a description, the same for
    strength by odor at the APL16 somata (from archive_scores_<reading>.csv). No variant is chosen by any of these."""
    rows = []
    for rd in ('curve_noise', 'raw'):
        t = pd.read_csv(H.results_dir(EXPERIMENT) / f'archive_tradeoff_{rd}.csv', low_memory=False)
        best = I[f'best_{rd}'].set_index('measure')
        for m in ('error_on', 'error_off', 'cal_pattern_part', 'raw_pattern_part'):
            rho = float(t[m].rank().corr(t['cal_pattern_part'].rank()))
            ck.eq(f'tradeoff {rd} {m} rho with the calibrated pattern part', rho,
                  float(best.loc[m, 'rho_with_cal_pattern_part']), 1e-9)
            i = t[m].idxmin()
            ck.eq(f'tradeoff {rd} {m} best value', float(t.loc[i, m]), float(best.loc[m, 'best_value']), 1e-12)
            rows.append(dict(reading=rd, measure=m, best=f"{t.loc[i, 'folder']}/{t.loc[i, 'tag']}",
                             best_value=float(t.loc[i, m]), best_family=best.loc[m, 'best_family'],
                             second_family=best.loc[m, 'second_family'], rho_with_cal_pattern_part=rho,
                             n_variants=int(best.loc[m, 'n_variants'])))
        a = I[f'arch_{rd}']
        a = a[(a['inside_bands'].astype(str).str.lower() == 'true') & (a['duplicate'].astype(str).str.lower() == 'false')]
        for cond in ('on', 'off'):
            x = a[a['condition'] == cond][['folder', 'tag', 'view_strength']].merge(t[['folder', 'tag', 'cal_pattern_part']],
                                                                                    on=['folder', 'tag'])
            ck.eq(f'tradeoff {rd} strength {cond} variants', len(x), len(t))
            i = x['view_strength'].idxmin()
            rows.append(dict(reading=rd, measure=f'strength_{cond}', best=f"{x.loc[i, 'folder']}/{x.loc[i, 'tag']}",
                             best_value=float(x.loc[i, 'view_strength']), best_family=x.loc[i, 'folder'],
                             rho_with_cal_pattern_part=float(x['view_strength'].rank().corr(x['cal_pattern_part'].rank())),
                             n_variants=len(x)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the soma criteria

def soma_criteria(I: dict, ck: Checks) -> pd.DataFrame:
    sc = I['lock']['soma_criteria']
    crit = I['lock']['stages']['kc_soma_criteria']
    for k, v in sc['criteria'].items():
        ck.eq(f'criterion {k}', v, crit[k])
    rows = []
    for r in sc['table']:
        rows.append(dict(r))
    t = pd.DataFrame(rows)
    mg = sc['megamat']
    ck.eq('megamat rate mean', _mean(mg['per_fly_rate'].values()), mg['rate_mean'], 1e-12)
    ck.eq('megamat rate at 0.0915', round(mg['rate_mean'], 4), sc['target_rate'])
    v2 = sc['validation2']
    for grp, d in v2['groups'].items():
        ck.eq(f'v2 {grp} factor', d['factor'], mg['noise_bridge_median'] / d['noise_bridge_median'], 1e-12)
        ck.eq(f'v2 {grp} n', d['n_flies'], len(d['flies']))
    ck.eq('megamat noise bridge = mean of the two middle flies', mg['noise_bridge_median'],
          float(np.median(list(mg['noise_bridge_per_fly'].values()))), 1e-12)
    ck.eq('same fly factor', sc['same_fly']['factor_per_fly']['2023-04-11/3/megamat1_redo'],
          mg['noise_bridge_median'] / sc['same_fly']['noise_bridge_per_fly']['2023-04-11/3/megamat1_redo'], 1e-12)
    # extra rows: the groups' bridges, Remy's flies in the order of their bridges
    extra = []
    for grp, d in v2['groups'].items():
        extra.append(dict(cohort=f'validation2 {grp} noise bridge', n_flies=d['n_flies'],
                          bridge_median=d['noise_bridge_median'], bridge_min=d['noise_bridge_min'],
                          bridge_max=d['noise_bridge_max'], factor=d['factor']))
    for fly, b in sorted(mg['noise_bridge_per_fly'].items(), key=lambda kv: kv[1]):
        extra.append(dict(cohort=f"megamat fly {fly}", n_flies=1, bridge_median=b,
                          rate_mean=mg['per_fly_rate'][fly], rate_mean_at_apl16_criterion=mg['per_fly_rate_at_apl16_criterion'][fly]))
    sf = sc['same_fly']
    extra.append(dict(cohort='megamat, the same fly, noise bridge', n_flies=1,
                      bridge_median=list(sf['noise_bridge_per_fly'].values())[0]))
    mx = sc['mixtures']
    extra.append(dict(cohort='the mixtures\' z bridge', criterion=mx['z_criterion'],
                      criterion_before=mx['z_criterion_at_apl16_criterion']))
    return pd.concat([t, pd.DataFrame(extra)], ignore_index=True)


# ---------------------------------------------------------------------------
# the APL16 somata: the planner's probe and the committed numbers

def _probe(df: pd.DataFrame) -> pd.Series:
    """The probe's population error: strength by its rate part only, sharing the RMS of its two parts, breadth."""
    return np.sqrt((df['part_strength_rate'].astype(float) ** 2 + df['view_sharing'].astype(float) ** 2
                    + df['view_breadth'].astype(float) ** 2) / 3)


def apl16_models(I: dict, ck: Checks) -> pd.DataFrame:
    rows = []
    s2 = I['s2']
    for cond in ('on', 'off'):
        rec = dict(model='the best family, eight wirings', source='step 2 (E90 part B)', condition=cond)
        for rd in ('raw', 'curve_noise'):
            x = _rows(s2, source='E90:B', stage='kc_soma', panel='apl16', condition=cond, reading=rd).set_index('wiring')
            w = x.loc[list(W8)]
            p = _probe(w)
            rec.update({f'{rd}_probe_mean': _mean(p), f'{rd}_probe_min': float(p.min()), f'{rd}_probe_max': float(p.max()),
                        f'{rd}_error_mean': float(x.loc['mean', 'error']), f'{rd}_error_min': float(w['error'].min()),
                        f'{rd}_error_max': float(w['error'].max()),
                        f'{rd}_strength_rate': _mean(w['part_strength_rate']), f'{rd}_strength': float(x.loc['mean', 'view_strength']),
                        f'{rd}_sharing_who': _mean(w['part_sharing_who']),
                        f'{rd}_sharing_how_strongly': _mean(w['part_sharing_how_strongly']),
                        f'{rd}_sharing': float(x.loc['mean', 'view_sharing']), f'{rd}_breadth': float(x.loc['mean', 'view_breadth'])})
        rec['written'] = {'on': '1.01 (0.94 to 1.10)', 'off': '1.01 (0.94 to 1.08)'}[cond]
        rows.append(rec)
    for label, members, w_on, w_off in PROBE_ROWS:
        for cond, written in (('on', w_on), ('off', w_off)):
            rec = dict(model=label, source='the archive: ' + ', '.join(f'{f}/{t}' for f, t in members[:2])
                       + (f' and {len(members) - 2} more' if len(members) > 2 else ''), condition=cond, written=written)
            for rd in ('raw', 'curve_noise'):
                a = I[f'arch_{rd}']
                x = pd.concat([a[(a['folder'] == f) & (a['tag'] == t) & (a['condition'] == cond)] for f, t in members])
                if len(x) != len(members):
                    raise AssertionError(f'{label} {rd} {cond}: {len(x)} rows for {len(members)} variants')
                p = _probe(x)
                rec.update({f'{rd}_probe_mean': _mean(p), f'{rd}_probe_min': float(p.min()), f'{rd}_probe_max': float(p.max()),
                            f'{rd}_error_mean': _mean(x['error']), f'{rd}_error_min': float(x['error'].min()),
                            f'{rd}_error_max': float(x['error'].max()), f'{rd}_strength': _mean(x['view_strength']),
                            f'{rd}_strength_rate': _mean(x['part_strength_rate']), f'{rd}_sharing': _mean(x['view_sharing']),
                            f'{rd}_breadth': _mean(x['view_breadth'])})
                rec['inside_bands'] = ' / '.join(sorted(set(x['inside_bands'].astype(str))))
            rows.append(rec)
    # the flies: one recording against the others, the probe's way and the committed way
    for cond in ('on', 'off'):
        lo = I['floo'][I['floo']['population'] == pop_name('kc_soma', 'apl16', cond)]
        p = _probe(lo)
        fl = I['fv'][(I['fv']['population'] == pop_name('kc_soma', 'apl16', cond)) & (I['fv']['name'] == 'error')].iloc[0]
        rows.append(dict(model='one recording against the others', source='step 1 (flies_loo.csv)', condition=cond,
                         written={'on': '0.98 (0.56 to 1.29)', 'off': '0.93 (0.59 to 1.62)'}[cond],
                         raw_probe_mean=_mean(p), raw_probe_min=float(p.min()), raw_probe_max=float(p.max()),
                         raw_error_mean=_mean(lo['error']), raw_error_min=float(lo['error'].min()),
                         raw_error_max=float(lo['error'].max()), flies_floor=float(fl['floor'])))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# realism

def _reg(I: dict, k: str) -> float:
    return float(I['reg'][k]['value'])


def realism_apl16(I: dict, ck: Checks) -> pd.DataFrame:
    from kcq import apl_effect as AE
    bv = I['bv']
    rows = []
    cols = ('rate_on', 'rate_off', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv')
    for _, b in I['bands'].iterrows():
        x = bv[(bv['folder'] == b['folder']) & (bv['tag'] == b['tag'])]
        if len(x) != 1:
            raise AssertionError(f"{b['folder']}/{b['tag']}: {len(x)} rows in by_variant.csv")
        x = x.iloc[0]
        for k in ('rate_on', 'rate_off', 'opk_on_tv', 'opk_off_tv'):
            ck.eq(f"band {b['wiring']} {k} vs score_apl16_bands.csv", float(x[k]), float(b[k]), 1e-12)
        rows.append(dict(model='the best family', wiring=b['wiring'], folder=b['folder'], tag=b['tag'],
                         **{k: float(x[k]) for k in cols}, realism=str(x['realism'])))
    for m, (f, t) in BESIDE_APL16.items():
        x = bv[(bv['folder'] == f) & (bv['tag'] == t)].iloc[0]
        rows.append(dict(model=BESIDE_NAME[m], wiring='real', folder=f, tag=t, **{k: float(x[k]) for k in cols},
                         realism=str(x['realism'])))
    t = pd.DataFrame(rows)
    t['ratio'] = t['rate_off'] / t['rate_on']
    t['inside_bands'] = t['realism'].eq('inside the bands')
    data = dict(model='KC data', wiring='', folder='registry', tag='', rate_on=_reg(I, 'noise_ctrl_frac'),
                rate_off=_reg(I, 'noise_tnt_frac'), ratio=_reg(I, 'noise_rate_ratio_fraction'),
                rate_rho_on=_reg(I, 'ctrl_rate_rho_ceiling'), rate_rho_off=_reg(I, 'tnt_rate_rho_ceiling'),
                opk_on_tv=_reg(I, 'noise_ctrl_odors_per_kc_tv_max'), opk_off_tv=_reg(I, 'noise_tnt_odors_per_kc_tv_max'),
                realism=(f'bands: response rate with APL on {AE.RATE_ON_BAND[0]} to {AE.RATE_ON_BAND[1]}; ratio '
                         f'{AE.RATE_RATIO_BAND[0]} to {AE.RATE_RATIO_BAND[1]} (90 % interval '
                         + ' to '.join(f'{v:.2f}' for v in I['reg']['noise_rate_ratio_fraction']['ci'])
                         + f'); response rate ρ with APL on at least {AE.RATE_RHO_MIN:.2f}; odors per KC '
                         f"{_reg(I, 'noise_ctrl_odors_per_kc_tv_max'):.3f} / {_reg(I, 'noise_tnt_odors_per_kc_tv_max'):.3f}"))
    return pd.concat([t, pd.DataFrame([data])], ignore_index=True)


def operating_points(I: dict, ck: Checks) -> pd.DataFrame:
    real, var = I['real'], I['var'].set_index('tag')
    for _, r in real.iterrows():
        v = var.loc[r['tag']]
        for k in ('rate_on', 'rate_off', 'max_count_on', 'max_count_off', 'fixed_thr'):
            ck.eq(f"{r['tag']} {k} realism vs variants", float(r[k]), float(v[k]), 1e-12)
    real = real.merge(I['var'][['tag', 'reading', 'thr_factor', 'apl_factor', 'frac_on_target', 'frac_off_target',
                                'fit_fixed_thr', 'wAPLKC_scale']], on='tag')
    pops = I['lock']['populations']
    flies_rate = {'megamat': I['lock']['soma_criteria']['target_rate'],
                  'validation2': pops['kc_soma__validation2__wild_type']['response_rate']['mean_over_flies'],
                  'mixtures': I['var'].loc[I['var']['tag'] == 'mixtures_held_real', 'frac_on_target'].iloc[0]}
    rows = []
    for (model, panel, reading), x in real.groupby(['model', 'panel', 'reading'], sort=False):
        rec = dict(model=model, panel=panel, reading=reading, n=len(x), flies_rate_on=flies_rate[panel])
        for k in ('rate_on', 'rate_off', 'rate_ratio', 'max_count_on', 'max_count_off', 'fixed_thr', 'thr_factor',
                  'apl_factor', 'frac_on_target', 'frac_off_target'):
            v = pd.to_numeric(x[k], errors='coerce')
            rec[f'{k}_min'], rec[f'{k}_max'] = float(v.min()), float(v.max())
        if reading == 'held':
            f = real[(real['model'] == model) & (real['panel'] == panel) & (real['reading'] == 'fitted')].set_index('wiring')
            mc_on = x.set_index('wiring')['max_count_on'] / f['max_count_on']
            mc_off = x.set_index('wiring')['max_count_off'] / f['max_count_off']
            rec['max_count_vs_fitted_max'] = float(np.maximum(mc_on, mc_off).max())
            thr = x.set_index('wiring')['fixed_thr'] / f['fixed_thr']
            ck.eq(f'{model} {panel} thr factor', float((thr - x.set_index("wiring")['thr_factor']).abs().max()), 0.0, 1e-9)
        rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the readings' numbers

def _decimals(s: str) -> int:
    s = s.strip().lstrip('+−-')
    return len(s.split('.')[1]) if '.' in s else 0


def planner_check(I, st, se, wse, stp, pan, con, br, soma, a16, ops) -> pd.DataFrame:
    rows = []

    def add(reading, what, written, value):
        d = _decimals(written)
        wv = float(written.replace('−', '-'))
        rows.append(dict(reading=reading, what=what, written=written, script=float(value),
                         reproduced=bool(round(float(value), d) == round(wv, d)),
                         truncated_matches=bool(math.floor(abs(float(value)) * 10 ** d) / 10 ** d == abs(wv))))

    def cons(rd, version, measure, cond, vs='inside_bands'):
        x = con[(con['reading'] == rd) & (con['version'] == version) & (con['measure'] == measure)
                & (con['condition'] == cond) & (con['variant_set'] == vs)]
        return float(x['split_half_rho'].iloc[0]) if len(x) else np.nan

    def stv(group, stage, panel, col, cond='on'):
        return float(_rows(st, group=group, stage=stage, panel=panel, condition=cond).iloc[0][col])

    lock = I['lock']['soma_criteria']
    fmc = I['fmc']
    R1 = 'planner, open readings (11:11; chat 11:34)'
    m = fmc[(fmc['name'] == 'error') & (fmc['reference'] == 'control') & (fmc['scored'] == 'megamat')].set_index('criterion')
    add(R1, "Remy's megamat flies against the control flies, 12 odors, lock criterion", '0.98', m.loc['lock', 'mean_views'])
    add(R1, "the same at the APL16 criterion", '3.13', m.loc['apl16', 'mean_views'])
    add(R1, "Remy's megamat criterion (baseline-noise SDs)", '5.95', lock['criteria']['megamat/wild_type'])
    g = lock['validation2']['groups']
    add(R1, 'tsub03 noise bridge, lowest', '0.77', g['tsub03']['noise_bridge_min'])
    add(R1, 'tsub03 noise bridge, highest', '1.09', g['tsub03']['noise_bridge_max'])
    add(R1, 'other recordings noise bridge, lowest', '0.46', g['other']['noise_bridge_min'])
    add(R1, 'other recordings noise bridge, highest', '0.60', g['other']['noise_bridge_max'])
    add(R1, 'kmix held rate, about', '0.056', lock['mixtures']['panels']['kmix']['rate_mean'])
    add(R1, 'cmix held rate, about', '0.052', lock['mixtures']['panels']['cmix']['rate_mean'])
    b = br[(br['reading'] == 'curve_noise') & (br['condition'] == 'on') & (br['measure'] == 'breadth')].iloc[0]
    add(R1, 'median breadth error inside the bands, APL on (archive of 11:32)', '0.26', b['median'])
    add(R1, "the control flies' breadth floor", '0.23', b['flies_floor'])
    R2 = 'viki, lock rewrite and step 1 (12:05)'
    mg = lock['megamat']
    add(R2, 'megamat criterion', '5.9486', lock['criteria']['megamat/wild_type'])
    for fly, w in (('2022-10-10/1', '0.062'), ('2022-10-10/2', '0.054'), ('2022-10-11/1', '0.124'), ('2022-11-10/1', '0.126')):
        add(R2, f'megamat fly {fly} response rate', w, mg['per_fly_rate'][fly])
    add(R2, 'tsub03 factor', '1.032', g['tsub03']['factor'])
    add(R2, 'other factor', '1.767', g['other']['factor'])
    add(R2, 'same-fly factor', '1.084', list(lock['same_fly']['factor_per_fly'].values())[0])
    add(R2, 'z bridge before', '0.7057', lock['mixtures']['z_criterion_at_apl16_criterion'])
    add(R2, 'z bridge now', '1.9086', lock['mixtures']['z_criterion'])
    add(R2, 'validation2 tsub03 rate', '0.049', g['tsub03']['rate_mean'])
    add(R2, 'validation2 other rate', '0.053', g['other']['rate_mean'])
    add(R2, 'validation2 tsub03 rate before', '0.147', g['tsub03']['rate_mean_at_apl16_criterion'])
    add(R2, 'validation2 other rate before', '0.091', g['other']['rate_mean_at_apl16_criterion'])
    add(R2, 'kmix rate', '0.0558', lock['mixtures']['panels']['kmix']['rate_mean'])
    add(R2, 'cmix rate', '0.0516', lock['mixtures']['panels']['cmix']['rate_mean'])
    add(R2, 'same-fly rate', '0.069', lock['same_fly']['rate_mean'])
    add(R2, 'bridge ratio tsub03 / other', '1.71', g['tsub03']['noise_bridge_median'] / g['other']['noise_bridge_median'])
    add(R2, 'binned fly 2022-11-10/1 bridge', '1.134', mg['noise_bridge_per_fly']['2022-11-10/1'])
    add(R2, 'megamat reference bridge', '0.9107', mg['noise_bridge_median'])
    fns = I['fns'][I['fns']['name'] == 'error'].set_index('downstream')
    for k, w in (('kc_soma__megamat__wild_type', '2.28'), ('kc_soma__validation2__wild_type', '2.91'),
                 ('kc_soma__kmix__wild_type', '5.72'), ('kc_soma__cmix__wild_type', '7.04')):
        add(R2, f'no-step yardstick into {k}', w, fns.loc[k, 'no_step'])
    R3 = 'gate planner (15:03; chat 15:19, measure_gate.json)'
    add(R3, 'population error split-half ρ, APL on', '0.79', cons('curve_noise', 'half', 'error', 'on'))
    add(R3, 'population error split-half ρ, APL off', '0.52', cons('curve_noise', 'half', 'error', 'off'))
    add(R3, 'APL off, highest over versions and readings', '0.79',
        max(cons(rd, v, 'error', 'off') for rd in ('curve_noise', 'raw') for v in ('half', 'shared_curve', 'all_units',
                                                                                  'shared_curve_all_units')))
    g2 = I['g2']
    g2 = g2[g2.apply(lambda r: r['reading'] == READINGS[r['stage']][0], axis=1)].drop_duplicates(
        ['source', 'stage', 'panel', 'condition'])
    add(R3, 'wiring check, lowest mean', '0.96', g2['mean_error'].min())
    add(R3, 'wiring check, highest mean', '1.05', g2['mean_error'].max())
    for v, on, off in (('strength', '0.87', '0.91'), ('sharing', '0.87', '0.73'), ('breadth', '0.35', '0.04')):
        add(R3, f'{v} split-half ρ, APL on', on, cons('curve_noise', 'half', f'view_{v}', 'on'))
        add(R3, f'{v} split-half ρ, APL off', off, cons('curve_noise', 'half', f'view_{v}', 'off'))
    add(R3, 'sharing at the megamat PN dendrites (as fitted)', '1.23', stv('megamat_fitted', 'pn_dendrite', 'megamat', 'default_view_sharing_mean'))
    add(R3, 'sharing at the megamat somata, as fitted', '1.13', stv('megamat_fitted', 'kc_soma', 'megamat', 'default_view_sharing_mean'))
    add(R3, 'sharing at the megamat somata, held', '1.41', stv('megamat_held', 'kc_soma', 'megamat', 'default_view_sharing_mean'))
    s = lambda gr, col='stage_error_three': float(_rows(se, group=gr, wiring='model')[col].iloc[0])
    add(R3, 'megamat stage error as fitted', '0.99', s('megamat_fitted'))
    add(R3, 'megamat stage error held', '1.02', s('megamat_held'))
    add(R3, 'the firing-ceiling fit', '1.00', s('firing_ceiling_megamat'))
    add(R3, 'the bouton baseline', '1.23', s('bouton_baseline_megamat'))
    add(R3, "steps passing under 'one_fly' (as fitted)", '3', int(_rows(stp, group='megamat_fitted')['one_fly_every_view_6of8'].sum()))
    db = _rows(stp, group='megamat_fitted', step='pn_dendrite->pn_bouton').iloc[0]
    add(R3, 'PN dendrites to boutons: the model downstream', '1.05', db['model_downstream_error'])
    add(R3, 'PN dendrites to boutons: no step', '1.09', db['no_step_error'])
    for p, w in (('megamat', '0.92'), ('validation2', '1.18'), ('kmix', '2.69'), ('cmix', '2.36')):
        add(R3, f'{p} KC somata population error', w, _rows(pan, panel=p).iloc[0]['fitted_error'])
    add(R3, 'the APL16 panel\'s value', '0.87', _rows(pan, panel='megamat').iloc[0]['apl16_error'])
    add(R3, 'strength-only megamat stage error as fitted', '0.93', s('megamat_fitted', 'stage_error_strength'))
    add(R3, 'strength-only megamat stage error held', '0.92', s('megamat_held', 'stage_error_strength'))
    n_ins = int(sum(r[f'{k}_class_strength'] == 'inside' for _, r in st[st['group'].isin(['megamat_fitted', 'megamat_held'])
                                                                          & st['stage'].isin(ERROR_STAGES)].iterrows()
                    for k in ('default', 'second')))
    add(R3, 'megamat stage-error stages inside on strength, both readings, fitted and held (of 16)', '16', n_ins)
    R4 = 'measure adopted (Tom 15:21, Yang 15:28; measure_gate.json)'
    add(R4, 'strength split-half ρ, APL on', '0.867', cons('curve_noise', 'half', 'view_strength', 'on'))
    add(R4, 'strength split-half ρ, APL off', '0.905', cons('curve_noise', 'half', 'view_strength', 'off'))
    add(R4, 'strength-only stage error as fitted (gate planner\'s read)', '0.927', s('megamat_fitted', 'stage_error_strength'))
    add(R4, 'strength-only stage error held (gate planner\'s read)', '0.918', s('megamat_held', 'stage_error_strength'))
    R5 = "hal's E135 subagent (15:53): 3 SE over the wirings, held runs"
    w3 = wse[wse['group'] == 'megamat_held'].set_index('stage')
    for stg, ws, we in (('pn_dendrite', '0.073', '0.193'), ('pn_bouton', '0.022', '0.036'), ('kc_claw', '0.015', '0.012'),
                        ('kc_soma', '0.034', '0.051')):
        add(R5, f'{stg} strength', ws, w3.loc[stg, 'view_strength_3se'])
        add(R5, f'{stg} three-view error', we, w3.loc[stg, 'error_3se'])
    add(R5, 'held stage error on strength', '0.918', s('megamat_held', 'stage_error_strength'))
    add(R5, 'held stage error on three views', '1.017', s('megamat_held'))
    R6 = "disco's E136 planner (15:44) and E136's entry"
    for stg, w in (('pn_dendrite', '1.270'), ('pn_bouton', '1.276'), ('kc_claw', '1.246'), ('kc_soma', '1.233')):
        add(R6, f'the most atypical fly on strength, {stg}', w, stv('megamat_fitted', stg, 'megamat', 'strength_loo_max'))
    fam = st[st['group'].isin(['megamat_fitted', 'megamat_held']) & st['stage'].isin(ERROR_STAGES)]
    big = max(float(r[f'{k}_view_strength_max']) for _, r in fam.iterrows() for k in ('default', 'second'))
    add(R6, 'the family\'s largest strength value at the four stages (both readings, fitted and held)', '1.17', big)
    for p, w in (('validation2', '1.04'), ('kmix', '2.52'), ('cmix', '0.75')):
        add(R6, f'{p} somata on strength, held', w, _rows(pan, panel=p).iloc[0]['held_strength'])
    for p, w in (('validation2', '1.01'), ('kmix', '2.53'), ('cmix', '0.98')):
        add(R6, f'{p} ORN recordings as the prediction, strength', w, _rows(pan, panel=p).iloc[0]['no_step_strength'])
    R7 = "hal's two-view table (16:00)"
    for rd, v, on, off in (('curve_noise', 'shared_curve', '0.81', '0.78'), ('curve_noise', 'shared_curve_all_units', '0.88', '0.83'),
                           ('raw', 'half', '0.81', '0.79'), ('raw', 'all_units', '0.87', '0.83')):
        add(R7, f'two-view error, {rd}, {v}, APL on', on, cons(rd, v, 'error_strength_sharing', 'on'))
        add(R7, f'two-view error, {rd}, {v}, APL off', off, cons(rd, v, 'error_strength_sharing', 'off'))
    R0 = "the planner's probe (09:11; the entry's table; raw counts, strength by its rate part only)"
    for _, r in a16.iterrows():
        if not isinstance(r.get('written'), str):
            continue
        w = r['written']
        main = w.split(' (')[0]
        add(R0, f"{r['model']}, APL {r['condition']}", main, r['raw_probe_mean'])
        if '(' in w:
            lo, hi = w.split('(')[1].rstrip(')').split(' to ')
            add(R0, f"{r['model']}, APL {r['condition']}, lowest", lo, r['raw_probe_min'])
            add(R0, f"{r['model']}, APL {r['condition']}, highest", hi, r['raw_probe_max'])
    fam = a16[a16['model'] == 'the best family, eight wirings'].set_index('condition')
    for cond, vals in (('on', ('1.29', '1.13', '1.15', '0.28')), ('off', ('1.12', '0.98', '1.23', '0.74'))):
        for col, w in zip(('raw_strength_rate', 'raw_sharing_who', 'raw_sharing_how_strongly', 'raw_breadth'), vals):
            add(R0, f'the family, APL {cond}, {col[4:]}', w, fam.loc[cond, col])
    fl = a16[a16['model'] == 'one recording against the others'].set_index('condition')
    add(R0, "the flies' floor, APL on", '0.25', fl.loc['on', 'flies_floor'])
    add(R0, "the flies' floor, APL off", '0.28', fl.loc['off', 'flies_floor'])
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# markdown

def _fmt(x, sp='.3f') -> str:
    if isinstance(x, str):
        return x
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return '–'
    return format(float(x), sp).replace('-', '−')


def _rng(a, b, sp='.3f') -> str:
    return _fmt(a, sp) if format(float(a), sp) == format(float(b), sp) else f'{_fmt(a, sp)} to {_fmt(b, sp)}'


ROW_LABEL = {
    ('megamat_fitted', 'megamat'): 'megamat, as fitted', ('megamat_held', 'megamat'): 'megamat, held at 0.0915',
    ('apl16', 'apl16'): 'APL16 panel, as fitted', ('validation2_fitted', 'validation2'): 'validation2, as fitted',
    ('validation2_held', 'validation2'): 'validation2, held at 6.8 %',
    ('mixtures_fitted', 'kmix'): 'kiwi mixtures, as fitted', ('mixtures_held', 'kmix'): 'kiwi mixtures, held at 0.0537',
    ('mixtures_fitted', 'cmix'): 'control mixtures, as fitted', ('mixtures_held', 'cmix'): 'control mixtures, held at 0.0537',
}


def markdown(st, se, wse, bes, stp, dend15, pan, con, br, tro, soma, a16, r16, ops, pc) -> str:
    L = ['# E134 write-up tables (generated by analyses/E134_writeup_tables.py; no simulation)', '']
    L += ['## The stage table (default readings; the family\'s eight wirings: mean, range, common)', '',
          '| row | stage | condition | flies | floor | most atypical fly: error / strength | population error | range | common | '
          'strength | sharing | breadth | matrix | error with matrix | energy guard | inside of 8 (three views / strength), '
          'default ; second | class: three views / strength |', '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    for _, r in st.iterrows():
        L.append(f"| {ROW_LABEL[(r['group'], r['panel'])]} | {STAGE_NAME[r['stage']]} | {r['condition']} | {int(r['n_flies'])} | "
                 f"{_fmt(r['error_floor'])} | {_fmt(r['error_loo_max'])} / {_fmt(r['strength_loo_max'])} | "
                 f"{_fmt(r['default_error_mean'])} | {_rng(r['default_error_min'], r['default_error_max'])} | "
                 f"{_fmt(r['default_error_common'])} | {_fmt(r['default_view_strength_mean'])} | "
                 f"{_fmt(r['default_view_sharing_mean'])} | {_fmt(r['default_view_breadth_mean'])} | "
                 f"{_fmt(r.get('default_view_matrix_mean', np.nan))} | {_fmt(r.get('default_error_with_matrix_mean', np.nan))} | "
                 f"{_fmt(r.get('default_energy_ratio_mean', np.nan), '.2f')} | "
                 f"{int(r['default_n_inside_three'])} / {int(r['default_n_inside_strength'])} ; "
                 f"{int(r['second_n_inside_three'])} / {int(r['second_n_inside_strength'])} ({r['second_reading']}) | "
                 f"{r['class_three']} / {r['class_strength']}"
                 + (f" (four views: {r['class_four']})" if isinstance(r.get('class_four'), str) else '') + ' |')
    L += ['', '### The second readings', '', '| row | stage | condition | reading | population error | range | strength | range | sharing | breadth |',
          '|---|---|---|---|---:|---:|---:|---:|---:|---:|']
    for _, r in st.iterrows():
        L.append(f"| {ROW_LABEL[(r['group'], r['panel'])]} | {STAGE_NAME[r['stage']]} | {r['condition']} | {r['second_reading']} | "
                 f"{_fmt(r['second_error_mean'])} | {_rng(r['second_error_min'], r['second_error_max'])} | "
                 f"{_fmt(r['second_view_strength_mean'])} | {_rng(r['second_view_strength_min'], r['second_view_strength_max'])} | "
                 f"{_fmt(r['second_view_sharing_mean'])} | {_fmt(r['second_view_breadth_mean'])} |")
    L += ['', '## The megamat stage error (the mean over the wirings; the range; the common error)', '',
          '| group | three views (dendrites on four) | range | on three views | strength alone | range | sharing alone | '
          'breadth alone | strength and sharing | common: three views / strength |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for g in STAGE_ERROR_GROUPS:
        r = _rows(se, group=g, wiring='model').iloc[0]
        c = _rows(se, group=g, wiring='common')
        com = f"{_fmt(c['stage_error_three'].iloc[0])} / {_fmt(c['stage_error_strength'].iloc[0])}" if len(c) else '–'
        L.append(f"| {g} | {_fmt(r['stage_error_three'])} | {_rng(r['stage_error_three_min'], r['stage_error_three_max'])} | "
                 f"{_fmt(r['stage_error_three_flat'])} | {_fmt(r['stage_error_strength'])} | "
                 f"{_rng(r['stage_error_strength_min'], r['stage_error_strength_max'])} | {_fmt(r['stage_error_sharing'])} | "
                 f"{_fmt(r['stage_error_breadth'])} | {_fmt(r['stage_error_two'])} | {com} |")
    L += ['', '### Three standard errors over the eight wirings (default reading)', '',
          '| group | stage | strength | three views | four views |', '|---|---|---:|---:|---:|']
    for _, r in wse.iterrows():
        L.append(f"| {r['group']} | {STAGE_NAME[r['stage']]} | {_fmt(r.get('view_strength_3se'))} | {_fmt(r.get('error_3se'))} | "
                 f"{_fmt(r.get('error_with_matrix_3se', np.nan))} |")
    L += ['', '## The beside models (one wiring each; default readings unless named)', '',
          '| model | input | stage | panel | reading | population error | strength | sharing | breadth | most atypical fly | '
          'megamat stage error: three views / strength |', '|---|---|---|---|---|---:|---:|---:|---:|---:|---:|']
    for _, r in bes.iterrows():
        L.append(f"| {r['model']} | {r['input']} | {STAGE_NAME[r['stage']]} | {r['panel']} | {r['reading']} | {_fmt(r['error'])} | "
                 f"{_fmt(r['strength'])} | {_fmt(r['sharing'])} | {_fmt(r['breadth'])} | {_fmt(r['flies_loo_max'])} | "
                 f"{_fmt(r['stage_error_three'])} / {_fmt(r['stage_error_strength'])} |")
    L += ['', '## The steps (the class on the wirings\' mean change, the corrected \'sets\' null)', '',
          '| group | step | class | no step / the model: population error | no step / the model: strength | model read on | '
          'inside on the mean change: strength / sharing / breadth / matrix | wirings inside, sets: s / sh / b / m | '
          'wirings inside, one fly: s / sh / b / m |', '|---|---|---|---:|---:|---|---|---|---|']
    yn = lambda x: '–' if x is None or (isinstance(x, float) and not np.isfinite(x)) else ('yes' if x else 'no')
    for _, r in stp.iterrows():
        vs = [v for v in VIEWS if f'{v}_sets_n_inside' in r and np.isfinite(r.get(f'{v}_sets_n_inside', np.nan))]
        L.append(f"| {r['group']} | {r['step']} | {r['cls']} | {_fmt(r['no_step_error'])} / {_fmt(r['model_downstream_error'])} | "
                 f"{_fmt(r['no_step_strength'])} / {_fmt(r.get('model_strength', np.nan))} | {r.get('model_on', '–')} | "
                 + ' / '.join(yn(r.get(f'{v}_inside_mean_change')) for v in vs) + ' | '
                 + ' / '.join(str(int(r[f'{v}_sets_n_inside'])) for v in vs) + ' | '
                 + ' / '.join(str(int(r[f'{v}_one_fly_n_inside'])) for v in vs) + ' |')
    L += ['', '### Each part: the wirings\' mean change from the flies\' (distance) against its threshold; the flies\' own change', '',
          '| group | step | part | distance | threshold | inside | flies\' RMS change | entries inside their interval |', '|---|---|---|---:|---:|---|---:|---|']
    for _, r in stp.iterrows():
        for part in ('strength/rate', 'strength/mean', 'sharing/who', 'sharing/how_strongly', 'breadth/bins', 'matrix/values'):
            if f'{part}_distance' not in r or not np.isfinite(r.get(f'{part}_distance', np.nan)):
                continue
            L.append(f"| {r['group']} | {r['step']} | {part} | {_fmt(r[f'{part}_distance'])} | {_fmt(r[f'{part}_threshold'])} | "
                     f"{yn(r[f'{part}_distance'] <= r[f'{part}_threshold'])} | {_fmt(r.get(f'{part}_flies_rms_change', np.nan))} | "
                     f"{r.get(f'{part}_entries_inside', '–')} |")
    if dend15 is not None:
        L += ['', '### The PN dendrites on the 15 odors they share with the ORN flies', '',
              '| group | wiring | population error | four views | strength | sharing | breadth | matrix |', '|---|---|---:|---:|---:|---:|---:|---:|']
        for _, r in dend15.iterrows():
            L.append(f"| {r['group']} | {r['wiring']} | {_fmt(r['error'])} | {_fmt(r['error_with_matrix'])} | {_fmt(r['view_strength'])} | "
                     f"{_fmt(r['view_sharing'])} | {_fmt(r['view_breadth'])} | {_fmt(r['view_matrix'])} |")
    L += ['', '## The panels (KC somata; the family\'s mean over eight wirings)', '',
          '| panel | flies | odors | floor | most atypical fly: error / strength | as fitted: error (range) | held: error | '
          'as fitted: strength | held: strength | APL16: error / strength | class: error / strength | ORNs as the prediction: '
          'error / strength | noise on the model |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|']
    for _, r in pan.iterrows():
        L.append(f"| {r['panel']} | {int(r['n_flies'])} | {int(r['n_odors'])} | {_fmt(r['floor'])} | {_fmt(r['loo_max'])} / "
                 f"{_fmt(r['strength_loo_max'])} | {_fmt(r['fitted_error'])} ({_rng(r['fitted_error_min'], r['fitted_error_max'])}) | "
                 f"{_fmt(r['held_error'])} | {_fmt(r['fitted_strength'])} | {_fmt(r['held_strength'])} | "
                 f"{_fmt(r['apl16_error'])} / {_fmt(r['apl16_strength'])} | {r['cls_three']} / {r['cls_strength']} | "
                 f"{_fmt(r['no_step_error'])} / {_fmt(r['no_step_strength'])} (from {STAGE_NAME[r['no_step_from']]}) | {r['model_noise']} |")
    L += ['', '## The measure\'s consistency (split-half ρ over the archived variants inside the bands; mean over splits)', '']
    inb = con[con['variant_set'] == 'inside_bands']
    measures = ('error', 'error_strength_sharing', 'view_strength', 'view_sharing', 'view_breadth', 'part_strength_rate',
                'part_strength_mean', 'part_sharing_who', 'part_sharing_how_strongly', 'dist_rho', 'rate_rho', 'energy_ratio')
    versions = ('half', 'shared_curve', 'all_units', 'shared_curve_all_units')
    for rd in ('curve_noise', 'raw'):
        x = inb[inb['reading'] == rd]
        n_var, n_spl = int(x['n_variants'].iloc[0]), int(x['n_splits'].iloc[0])
        L += [f'### {rd} ({n_var} variants, {n_spl} splits; APL on / APL off)', '',
              '| measure | ' + ' | '.join(versions) + ' |', '|---|' + '---:|' * len(versions)]
        for m in measures:
            cells = []
            for v in versions:
                y = x[(x['measure'] == m) & (x['version'] == v)].set_index('condition')['split_half_rho']
                cells.append(f"{_fmt(y.get('on', np.nan))} / {_fmt(y.get('off', np.nan))}")
            L.append(f'| {m} | ' + ' | '.join(cells) + ' |')
        L.append('')
    L += ['## The archived variants\' view errors inside the bands at the APL16 somata, beside the flies', '',
          '| reading | condition | measure | variants | 5 % | 25 % | median | 75 % | 95 % | flies\' floor | most atypical fly |',
          '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for _, r in br.iterrows():
        L.append(f"| {r['reading']} | {r['condition']} | {r['measure']} | {int(r['n_variants'])} | {_fmt(r['q05'])} | {_fmt(r['q25'])} | "
                 f"{_fmt(r['median'])} | {_fmt(r['q75'])} | {_fmt(r['q95'])} | {_fmt(r['flies_floor'])} | {_fmt(r['flies_loo_max'])} |")
    L += ['', '## The trade-off over the archived variants inside the bands (a description; no variant is chosen by it)', '',
          '| reading | measure | best variant | its value | best family | second family | ρ over the variants with the '
          'calibrated pattern part | variants |', '|---|---|---|---:|---|---|---:|---:|']
    for _, r in tro.iterrows():
        L.append(f"| {r['reading']} | {r['measure']} | {r['best']} | {_fmt(r['best_value'])} | {r['best_family']} | "
                 f"{r.get('second_family') if isinstance(r.get('second_family'), str) else '–'} | "
                 f"{_fmt(r['rho_with_cal_pattern_part'], '.2f')} | {int(r['n_variants'])} |")
    L += ['', '## The KC somata\'s criteria per set of recordings (the lock)', '',
          '| cohort | flies | criterion | values times | response rate (range) | at 1.9684 / before | noise bridge (range) |',
          '|---|---:|---:|---:|---:|---:|---:|']
    for _, r in soma.iterrows():
        rate = (f"{_fmt(r.get('rate_mean'), '.4f')} ({_rng(r['rate_min'], r['rate_max'])})"
                if np.isfinite(r.get('rate_min', np.nan)) else _fmt(r.get('rate_mean', np.nan), '.4f'))
        br_ = (f"{_fmt(r.get('bridge_median'))} ({_rng(r['bridge_min'], r['bridge_max'])})"
               if np.isfinite(r.get('bridge_min', np.nan)) else _fmt(r.get('bridge_median', np.nan)))
        before = r.get('rate_mean_at_apl16_criterion', np.nan)
        if not np.isfinite(before if isinstance(before, float) else np.nan) and np.isfinite(r.get('criterion_before', np.nan)):
            before = r['criterion_before']
        L.append(f"| {r['cohort']} | {_fmt(r.get('n_flies'), '.0f')} | {_fmt(r.get('criterion'), '.4f')} | "
                 f"{_fmt(r.get('factor'), '.3f')} | {rate} | {_fmt(before, '.4f')} | {br_} |")
    L += ['', '## The APL16 somata: the planner\'s probe reproduced, and the committed population error', '',
          '| model | condition | written | probe, raw: mean (range) | committed, raw: mean (range) | committed, default reading | '
          'raw: strength (rate part) / sharing / breadth | inside the bands |', '|---|---|---:|---:|---:|---:|---:|---|']
    for _, r in a16.iterrows():
        wr = r.get('written') if isinstance(r.get('written'), str) else '–'
        L.append(f"| {r['model']} | {r['condition']} | {wr} | {_fmt(r['raw_probe_mean'])} "
                 f"({_rng(r['raw_probe_min'], r['raw_probe_max'])}) | {_fmt(r['raw_error_mean'])} ({_rng(r['raw_error_min'], r['raw_error_max'])}) | "
                 f"{_fmt(r.get('curve_noise_error_mean', np.nan))} | {_fmt(r.get('raw_strength', np.nan))} "
                 f"({_fmt(r.get('raw_strength_rate', np.nan))}) / {_fmt(r.get('raw_sharing', np.nan))} / {_fmt(r.get('raw_breadth', np.nan))} | "
                 f"{r.get('inside_bands', '–') if isinstance(r.get('inside_bands'), str) else '–'} |")
    L += ['', '## Realism on the APL16 panel (results/apl_effect/by_variant.csv)', '',
          '| model | wiring | response rate, APL on | APL off | off / on | response rate ρ, on | off | odors per KC (TV), on / off | bands |',
          '|---|---|---:|---:|---:|---:|---:|---:|---|']
    for _, r in r16.iterrows():
        L.append(f"| {r['model']} | {r['wiring']} | {_fmt(r['rate_on'], '.4f')} | {_fmt(r['rate_off'], '.4f')} | {_fmt(r['ratio'], '.2f')} | "
                 f"{_fmt(r['rate_rho_on'], '.2f')} | {_fmt(r['rate_rho_off'], '.2f')} | {_fmt(r['opk_on_tv'])} / {_fmt(r['opk_off_tv'])} | {r['realism']} |")
    L += ['', '## The operating points of step 3\'s runs (score_realism.csv, variants.csv)', '',
          '| model | panel | input | runs | the flies\' response rate | rate, APL on | rate, APL off | off / on | largest count, on / off | '
          'threshold | threshold × fit\'s | tied APL scales × fit\'s | largest count × fitted |', '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for _, r in ops.iterrows():
        L.append(f"| {r['model']} | {r['panel']} | {r['reading']} | {int(r['n'])} | {_fmt(r['flies_rate_on'], '.4f')} | "
                 f"{_rng(r['rate_on_min'], r['rate_on_max'], '.4f')} | {_rng(r['rate_off_min'], r['rate_off_max'], '.4f')} | "
                 f"{_rng(r['rate_ratio_min'], r['rate_ratio_max'], '.2f')} | {_fmt(r['max_count_on_max'], '.0f')} / "
                 f"{_fmt(r['max_count_off_max'], '.0f')} | {_rng(r['fixed_thr_min'], r['fixed_thr_max'], '.1f')} | "
                 f"{_rng(r['thr_factor_min'], r['thr_factor_max'], '.2f') if np.isfinite(r['thr_factor_min']) else '–'} | "
                 f"{_rng(r['apl_factor_min'], r['apl_factor_max'], '.2f') if np.isfinite(r['apl_factor_min']) else '–'} | "
                 f"{_fmt(r.get('max_count_vs_fitted_max', np.nan), '.2f')} |")
    n_ok = int(pc['reproduced'].sum())
    L += ['', "## The readings' numbers", '', f'{n_ok} of {len(pc)} reproduce at the precision written ("cut off" where the '
          "written value is the script's cut off at that precision instead of rounded).", '',
          '| reading | what | written | script | reproduced |', '|---|---|---:|---:|---|']
    for _, r in pc.iterrows():
        y = 'yes' if r['reproduced'] else ('no (cut off)' if r['truncated_matches'] else 'no')
        L.append(f"| {r['reading']} | {r['what']} | {r['written']} | {_fmt(r['script'], '.6f')} | {y} |")
    return '\n'.join(L) + '\n'


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', type=Path, default=None, help='output folder (default: the results folder)')
    a = ap.parse_args(argv)
    res = H.results_dir(EXPERIMENT)
    out = a.out or res
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    I = load(res)
    ck = Checks()
    check_readings(ck)
    st = stage_table(I, ck)
    se = stage_error(I, ck)
    wse = wiring_se(I, ck)
    bes = beside(I, ck, se)
    dend15 = dendrites_on_orn_odors(I, ck)
    stp = steps(I, ck, st, dend15)
    pan = panels(I, ck, st)
    con = consistency(I, ck)
    br = breadth_spread(I)
    tro = tradeoff(I, ck)
    soma = soma_criteria(I, ck)
    a16 = apl16_models(I, ck)
    r16 = realism_apl16(I, ck)
    ops = operating_points(I, ck)
    pc = planner_check(I, st, se, wse, stp, pan, con, br, soma, a16, ops)
    tables = dict(stage_table=st, stage_error=se, wiring_se=wse, beside=bes, steps=stp, panels=pan, consistency=con,
                  breadth=br, tradeoff=tro, soma_criteria=soma, apl16_models=a16, realism_apl16=r16, operating_points=ops,
                  planner_check=pc)
    if dend15 is not None:
        tables['dendrites_15_odors'] = dend15
    for k, df in tables.items():
        df.to_csv(out / f'writeup_{k}.csv', index=False)
    (out / 'writeup_tables.md').write_text(markdown(st, se, wse, bes, stp, dend15, pan, con, br, tro, soma, a16, r16, ops,
                                                    pc))
    two_view_half = con[(con['reading'] == 'curve_noise') & (con['version'] == 'half')
                        & (con['measure'] == 'error_strength_sharing') & (con['variant_set'] == 'inside_bands')]
    inputs = ['score_scores.csv', 'score_classes.json', 'score_step_table.csv', 'score_stage_error.csv',
              'score_realism.csv', 'score_apl16_bands.csv', 'step2_scores.csv', 'step2_gate_wirings.csv',
              'measure_gate.json', 'flies_views.csv', 'flies_loo.csv', 'flies_no_step.csv',
              'flies_megamat_vs_control.csv', 'flies_steps_summary.csv', 'lock.json', 'variants.csv',
              'archive_consistency_curve_noise.csv', 'archive_consistency_raw.csv', 'archive_scores_curve_noise.csv',
              'archive_scores_raw.csv', 'archive_tradeoff_best_curve_noise.csv', 'archive_tradeoff_best_raw.csv',
              'archive_tradeoff_curve_noise.csv', 'archive_tradeoff_raw.csv']
    prov = dict(script=SCRIPT, generated=time.strftime('%Y-%m-%d %H:%M'), seconds=round(time.time() - t0, 1),
                checks_passed=ck.n, planner_numbers=int(len(pc)), planner_reproduced=int(pc['reproduced'].sum()),
                not_reproduced=pc[~pc['reproduced']][['reading', 'what', 'written', 'script',
                                                      'truncated_matches']].to_dict('records'),
                two_view_half_curve_noise=('written' if two_view_half['split_half_rho'].notna().any()
                                           else 'not yet written (NaN in archive_consistency_curve_noise.csv)'),
                dendrites_on_orn_odors='recomputed from outputs/E134_population_baseline/populations/score/'
                if dend15 is not None else 'not computed: the saved populations are not on this host',
                commits=H.code_commits(),
                inputs={**{f'results/{EXPERIMENT}/{n}': _sha256(res / n) for n in inputs},
                        'results/apl_effect/by_variant.csv': _sha256(H.shared_dir() / 'results' / 'apl_effect' / 'by_variant.csv'),
                        'results/data_targets/data_targets.json':
                            _sha256(H.shared_dir() / 'results' / 'data_targets' / 'data_targets.json')},
                note='No simulation. Every table is read from the score\'s, step 1\'s, step 2\'s, the archive scoring\'s '
                     'and the gate\'s saved files; the PN dendrites on the ORN odors are re-scored from the score\'s saved '
                     'populations by kcq.population_fit after reproducing the score\'s own numbers.')
    (out / 'provenance_writeup_tables.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    print(f'{ck.n} checks passed; readings\' numbers reproduced {int(pc["reproduced"].sum())} of {len(pc)}; wrote {out}',
          flush=True)
    for _, r in pc[~pc['reproduced']].iterrows():
        cut = ' (cut off, not rounded)' if r['truncated_matches'] else ''
        print(f"  not reproduced: {r['reading']}: {r['what']}: written {r['written']}, script {r['script']:.6f}{cut}",
              flush=True)


if __name__ == '__main__':
    main()
