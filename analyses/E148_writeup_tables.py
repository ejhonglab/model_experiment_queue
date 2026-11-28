"""The tables of E148's write-up (the megamat fits' saved evaluations read on sharing and breadth and under other
losses), from the files its parts wrote to results/E148_megamat_other_views/: step 0's step0_probe.csv, identities.csv,
step0.json and prediction.json; the room (room_summary.csv); which views the flies decide (views_decided_summary.csv,
sharing_split.csv); the losses (losses_summary.csv, losses_picks.csv); the claws (claws.csv, claws_evals.csv,
claws_odor_leads.csv); without the boutons (no_boutons.json); beside (beside_*.csv); the score (classes.json); the
views' summaries (views/<cell>_summary.json). Beside them: each cell's fit history (the picks' thresholds, scales and
response rates), E134's flies tables (the floors) and E134's lock (the flies' megamat response rate). No simulation, no
fit and no number typed.

Writes writeup_tables.md and writeup_values.csv (the same numbers, long: arm, stage, panel, condition, half, statistic,
value; the arm is a cell, a pooled set, a set of the saved settings or a saved setting's name, the condition a view or a
loss, the half 'a', 'b', 'all', a split's summary or a reference of flies) beside the parts' files.

Usage (from the repository root): ./venv/bin/python -m analyses.E148_writeup_tables [--score-dir <dir>]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from analyses import E136_writeup_tables as W
from experiments import E148_megamat_other_views as M

STAGES = W.STAGES
STAGE_NAMES = W.STAGE_NAMES
CELLS = tuple(M.CELLS)
CELL_NAMES = {'A1': 'A1 (S, half a)', 'X': 'X (T, half b)', 'F': 'F (S, half b)', 'W': 'W (T, half a)'}
SETS = ('A1', 'X', 'F', 'W', 'pool_S', 'pool_T', 'settings', 'settings_without_e135', 'inhadd1000', 'family_norm')
SET_NAMES = {'A1': 'A1', 'X': 'X', 'F': 'F', 'W': 'W', 'pool_S': 'S pooled (A1, F)', 'pool_T': 'T pooled (X, W)',
             'settings': 'the 23 settings', 'settings_without_e135': "without E135's four",
             'inhadd1000': 'the nine at pn_inhadd 1000', 'family_norm': "the five at the family's normalization"}
VIEW_NAMES = {'strength': 'strength by odor', 'sharing': 'sharing between odors', 'breadth': 'breadth among units',
              'matrix': 'the matrix', 'all': 'all views'}
ROOM_MODELS = {('family', 'T'): 'the family on T', ('family', 'S'): 'the family on S',
               ('E136:A1', 'T'): "A1's classed point on T (E136's reads)", ('E145_X', 'S'): "X's classed point on S", ('E145_F', 'T'): "F's classed point on T",
               ('E145_W', 'S'): "W's classed point on S"}
FACTORS = ('aplpn', 'kckc', 'inhadd', 'inhsc')
FACTOR_NAMES = {'aplpn': 'APL→bouton', 'kckc': 'KC–KC', 'inhadd': 'additive term (× family)',
                'inhsc': 'gain on pn_inhsc'}
SCALES = ('fixed_thr', 'wKCKC_scale', 'wAPLKC_scale', 'wAPLPN_scale')
SCALE_FMT = {'fixed_thr': '.1f', 'wKCKC_scale': '.1f', 'wAPLKC_scale': '.3f', 'wAPLPN_scale': '.4f'}
SCALE_NAMES = {'fixed_thr': 'somatic threshold', 'wKCKC_scale': 'KC–KC scale', 'wAPLKC_scale': 'APL→KC scale',
               'wAPLPN_scale': 'APL→bouton scale'}
_f, _mr, _yn = W._f, W._mr, W._yn


def _rng(x, fmt: str = '.3f') -> str:
    """A range as 'lo to hi', or one value where both print alike; '–' if empty."""
    x = pd.Series(x, dtype=float).dropna()
    if not len(x):
        return '–'
    lo, hi = _f(x.min(), fmt), _f(x.max(), fmt)
    return lo if lo == hi else f'{lo} to {hi}'


def _truth(v) -> bool:
    return str(v) == 'True'


def _views_of(stage: str) -> List[str]:
    return list(M.VIEWS) + (['matrix'] if stage == 'pn_dendrite' else [])


# ---------------------------------------------------------------------------
# inputs

def _rd(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, float_precision='round_trip') if p.exists() else pd.DataFrame()


def load(score_dir: Path, histories: bool = True) -> dict:
    """Everything the tables read: the parts' files from `score_dir`, the cells' fit histories (`histories`), E134's
    flies tables and the flies' megamat response rate."""
    names = ('step0_probe', 'identities', 'room_summary', 'views_decided_summary', 'sharing_split', 'losses_summary',
             'losses_picks', 'claws', 'claws_evals', 'claws_odor_leads', 'no_boutons_joint', 'beside_settings',
             'beside_regressions', 'beside_e147_spreads', 'beside_e147_populations')
    d = {n: _rd(score_dir / f'{n}.csv') for n in names}
    for n in ('classes', 'no_boutons', 'prediction', 'step0'):
        d[n] = W._read_json(score_dir / f'{n}.json')
    for n in ('score', 'step0', 'views_decided', 'losses', 'claws', 'no_boutons', 'beside', 'room'):
        d[f'prov_{n}'] = W._read_json(score_dir / f'provenance_{n}.json')
    d['views'] = {c: W._read_json(score_dir / 'views' / f'{c}_summary.json') for c in M.VIEW_CELLS}
    d['hist'] = {}
    d['flies_rate'] = None
    d['flies_views'] = pd.DataFrame()
    if histories:
        from experiments import E134_population_baseline as B
        for c in CELLS:
            p = M.cell_results(c) / f'fit_history_{M.ARM}.csv'
            if p.exists():
                d['hist'][c] = pd.read_csv(p, float_precision='round_trip')
        d['flies_rate'] = float(B.held_rate('megamat'))
        e = H_results(M.E134)
        d['flies_views'] = _rd(e / 'flies_views.csv')
    return d


def H_results(experiment: str) -> Path:
    from kcq import harness as H
    return H.results_dir(experiment)


def flies_rows(d: dict) -> Dict[tuple, dict]:
    """Per (stage, view): E134's floor and one fly against the others (mean, largest) on megamat's wild-type flies."""
    fv = d['flies_views']
    out = {}
    if not len(fv):
        return out
    x = fv[(fv.panel == 'megamat') & (fv.level == 'view')]
    if 'cohort' in x:
        x = x[x.cohort.astype(str) == 'wild_type']
    for r in x.itertuples():
        out[(str(r.stage), str(r.name))] = dict(floor=r.floor, loo_mean=r.loo_mean, loo_max=r.loo_max)
    return out


# ---------------------------------------------------------------------------
# step 0

def probe_table(d: dict, V: W.Values) -> List[str]:
    """The planner's probe as printed against the values reproduced from E136's and E145's rescore tables (every
    evaluation, as the planner read them), the held evaluations beside."""
    t = d['step0_probe']
    if not len(t):
        return []
    L = ["The planner's probe as the entry printed it, against the values step 0 reproduced from E136's and E145's "
         "rescore tables and fit histories (every evaluation, as the planner read them; the held evaluations beside). "
         "Pick factors: APL→bouton, KC–KC and the additive term as multiples of the family's.", '',
         '| statistic | cell | printed | reproduced | held evaluations | reproduced? | note |', '|---|---|---|---|---|---|---|']

    def fmt(s, raw=False, counts=False):
        try:
            v = json.loads(s)
        except (TypeError, ValueError):
            return str(s)
        xs = v if isinstance(v, list) else [v]
        if raw:
            return ' / '.join(str(x) for x in xs)
        return ' / '.join(str(int(x)) if counts else _f(x) for x in xs)
    for r in t.itertuples():
        note = r.note if isinstance(getattr(r, 'note', None), str) else ''
        counts = str(r.statistic).startswith(('holds_part', 'claws_worse', 'w_somata'))
        if r.statistic == 'consistency_L0' and r.cell == 'A1':
            note = (note + '; ' if note else '') + ("this row is reproduced over A1's 120 held evaluations, and the "
                                                    "column beside is over all 124")
        L.append(f'| {r.statistic} | {r.cell} | {fmt(r.printed, raw=True)} | {fmt(r.value, counts=counts)} | '
                 f'{fmt(r.value_held, counts=counts)} | '
                 f'{_yn(_truth(r.reproduced))} | {note} |')
        try:
            v = json.loads(r.value)
        except (TypeError, ValueError):
            continue
        for i, x in enumerate(v if isinstance(v, list) else [v]):
            V.add(r.cell, 'megamat', f'probe_{r.statistic}' + (f'_{i}' if isinstance(v, list) else ''), x,
                  condition='probe', half='all')
    n = int(t['reproduced'].map(_truth).sum())
    L += ['', f'{n} of {len(t)} printed values reproduced. The 0.040 bar of megamat\'s part of "holds" is the entry\'s '
          'and the probe\'s; at the control fit\'s 0.0403 A1\'s evaluation 9 (a fall of 0.04026 against half b) would '
          'not count (the lock\'s judgments).']
    return L + ['']


def identities_table(d: dict, V: W.Values) -> List[str]:
    """Every identity of step 0 by cell (identities.csv), the read-once check apart."""
    ids = d['identities']
    if not len(ids):
        return []
    x = ids[ids['kind'].astype(str).isin(['identity', 'beside'])]
    L = ['The identities (identities.csv): each cell\'s rescore table reproduced from its views by E136\'s own reading '
         '(the largest difference over every evaluation, measure, half and column), the consistencies, and E143\'s '
         'consistency table on strength from the saved settings\' views.', '',
         '| cell | identity | value | expected | passed |', '|---|---|---:|---:|---|']
    for r in x.itertuples():
        val = r.value
        fmt = '.2e' if 'reproduced' in str(r.identity) or 'table on strength' in str(r.identity) else '.4f'
        exp = _f(r.expected, '.0e' if fmt == '.2e' else '.4f')
        L.append(f'| {r.cell} | {r.identity} | {_f(val, fmt)} | {exp} | '
                 f"{_yn(_truth(r.passed)) if str(r.passed) in ('True', 'False') else 'beside'} |")
        V.add(r.cell, 'megamat', 'identity:' + str(r.identity), val, condition='step0', half='all')
    return L + ['']


def read_once_table(d: dict, V: W.Values) -> List[str]:
    """E143's read-once check by cell: on E136's split, the megamat stage error on strength at evaluation 0 and at the
    fit's best evaluation by E136's own reading and read once, with the largest difference over every evaluation; and
    the consistency on strength read each way."""
    ids = d['identities']
    if not len(ids):
        return []
    L = ["E143's read-once check (tolerance 0.02): the megamat stage error on strength in one-fly units by E136's own "
         "reading (each half's flies set the curve, the noise and the criterion) and read once (all the flies set them), "
         "on E136's split.", '',
         '| cell | half | evaluation | own reading | read once | difference | largest over every evaluation | passed |',
         '|---|---|---|---:|---:|---:|---:|---|']
    cons = {}
    for r in ids.itertuples():
        if str(r.identity).startswith('consistency') and 'from the views' in str(r.identity):
            cons[(r.cell, 'own' if 'own reading' in str(r.identity) else 'once')] = r.value
    for r in ids[ids['kind'].astype(str) == 'check'].itertuples():
        chk = json.loads(r.note)
        mx = chk.get('max_diff_every_eval') or {}
        for c in chk.get('checks') or []:
            L.append(f"| {r.cell} | {c['half']} | {c['eval']} | {_f(c['own'])} | {_f(c['once'])} | {_f(c['diff'])} | "
                     f"{_f(mx.get(c['half'] + '/stage_error_strength'))} | {_yn(c['diff'] <= M.CHECK_TOL)} |")
            V.add(r.cell, 'megamat', 'read_once_diff', c['diff'], condition=c['eval'], half=c['half'])
    L += ['', "The megamat stage error on strength's consistency between the halves of E136's split over each cell's "
          "evaluations, by E136's own reading and read once.", '',
          '| cell | own reading | read once |', '|---|---:|---:|']
    for c in CELLS:
        L.append(f"| {c} | {_f(cons.get((c, 'own')))} | {_f(cons.get((c, 'once')))} |")
        V.add(c, 'megamat', 'consistency_L0_own', cons.get((c, 'own')), condition='e136_split', half='all')
        V.add(c, 'megamat', 'consistency_L0_read_once', cons.get((c, 'once')), condition='e136_split', half='all')
    return L + ['']


def views_read_table(d: dict, V: W.Values) -> List[str]:
    """Where and how each cell's views were read (views/<cell>_summary.json)."""
    L = ["The views (design 2): every evaluation's views per stage, wiring and reading, read once per cell on the host "
         "that holds its stages.", '', '| cell | host | read | evaluations found / wanted | seconds | queue | dirty |',
         '|---|---|---|---:|---:|---|---|']
    for c in M.VIEW_CELLS:
        s = d['views'].get(c) or {}
        if not s:
            continue
        L.append(f"| {c} | {s.get('host')} | {s.get('time')} | {s.get('n_found')} / {s.get('n_wanted')} | "
                 f"{_f(s.get('seconds'), '.0f')} | {str(s.get('queue_commit'))[:7]} | {_yn(s.get('queue_dirty'))} |")
        V.add(c, 'views', 'n_found', s.get('n_found'), condition='views', half='all')
    return L + ['']


# ---------------------------------------------------------------------------
# design 3: the room

def stage_table(d: dict, V: W.Values) -> List[str]:
    """The room at the four megamat stages: per view, the flies' floor and one fly against the others (E134), the family
    and the classed points against all the flies (read once, mean over the wirings), and one half's mean view against
    the other (half to half, the median over the stage's splits)."""
    rs = d['room_summary']
    if not len(rs):
        return []
    fl = flies_rows(d)
    L = ["Each view's error in one-fly units at the four megamat stages (lower is better). The models against all the "
         "stage's flies, read once, the mean over the wirings named of each wiring's error; half to half is one half's "
         "mean view scored against the other half's, the median over the stage's splits (E143's: 35, 10, 3, 3). The "
         "flies' rows are E134's.", '']
    for v in ('strength', 'sharing', 'breadth'):
        L += [f'**{VIEW_NAMES[v]}**', '', '| model or yardstick | ' + ' | '.join(STAGE_NAMES[s] for s in STAGES) + ' |',
              '|---|' + '---:|' * len(STAGES)]
        L.append("| the flies' floor | " + ' | '.join(_f((fl.get((s, v)) or {}).get('floor'), '.2f') for s in STAGES)
                 + ' |')
        L.append('| one fly against the others, mean (largest) | ' + ' | '.join(
            f"{_f((fl.get((s, v)) or {}).get('loo_mean'), '.2f')} ({_f((fl.get((s, v)) or {}).get('loo_max'), '.2f')})"
            for s in STAGES) + ' |')
        x = rs[rs.view == v]
        hh = {s: x[x.stage == s].median_half_to_half.iloc[0] if len(x[x.stage == s]) else np.nan for s in STAGES}
        L.append('| half to half, median over splits | ' + ' | '.join(_f(hh[s]) for s in STAGES) + ' |')
        for s in STAGES:
            V.add('half_to_half', s, 'median', hh[s], condition=v, half='splits')
        for (m, wh), lab in ROOM_MODELS.items():
            y = x[(x.model == m) & (x.wirings == wh)].set_index('stage')
            if not len(y):
                continue
            L.append(f'| {lab} | ' + ' | '.join(_f(y.against_all.get(s)) for s in STAGES) + ' |')
            for s in STAGES:
                V.add(f'{m}@{wh}', s, 'against_all', y.against_all.get(s), condition=v, half='all')
        L.append('')
    return L


def room_table(d: dict, V: W.Values) -> List[str]:
    """The room per stage and view: the family and the classed points against each half of E136's split and the median
    over the stage's splits of each half, beside half to half."""
    rs = d['room_summary']
    if not len(rs):
        return []
    L = ["Against the halves (read once, one-fly units): on E136's split against half a and half b, and the median over "
         "the stage's splits; half to half beside (one half's mean view against the other's). A1's classed point on T "
         "is E136's reads there, the runs E145_A1 reads on T.", '']
    for s in STAGES:
        L += [f'**{STAGE_NAMES[s]}**', '',
              "| model | view | E136's split: half a | half b | median over splits: half a | half b | half to half, "
              "E136's split | median |", '|---|---|---:|---:|---:|---:|---:|---:|']
        x = rs[rs.stage == s]
        for v in _views_of(s) + ['all']:
            for (m, wh), lab in ROOM_MODELS.items():
                y = x[(x.view == v) & (x.model == m) & (x.wirings == wh)]
                if not len(y):
                    continue
                r = y.iloc[0]
                L.append(f'| {lab} | {VIEW_NAMES[v]} | {_f(r.e136_against_a)} | {_f(r.e136_against_b)} | '
                         f'{_f(r.median_a)} | {_f(r.median_b)} | {_f(r.e136_half_to_half)} | '
                         f'{_f(r.median_half_to_half)} |')
                for k in ('e136_against_a', 'e136_against_b', 'median_a', 'median_b'):
                    V.add(f'{m}@{wh}', s, k, getattr(r, k), condition=v, half=k.split('_')[-1])
        L.append('')
    return L


# ---------------------------------------------------------------------------
# design 4: which views the flies decide

def decided_table(d: dict, V: W.Values) -> List[str]:
    """Consistency per stage and view over every set, on every split of the stage: median / minimum, the class."""
    s = d['views_decided_summary']
    if not len(s):
        return []
    L = ['Consistency (ρ over a set\'s models between their errors against one half of a stage\'s flies and against the '
         'other) per stage and view, on every split of the stage: median / minimum. **Bold**: "is decided" (0.8 or more '
         'at the median, 0.5 or more on every split); "nr": "no reading" (the set\'s errors against all the flies range '
         'over less than twice the view\'s spread over the family\'s eight wirings).', '',
         'Every cell is read by E136\'s own reading (its read-once check failed), so its PN stages are read on E136\'s '
         'split alone and its claws and somata on their three splits; the settings are read once, on every split.', '',
         '| stage | view | ' + ' | '.join(SET_NAMES[k] for k in SETS) + ' |', '|---|---|' + '---:|' * len(SETS)]
    for st in STAGES:
        for v in _views_of(st) + ['all']:
            x = s[(s.stage == st) & (s.view == v)].set_index('set')
            cells = []
            for k in SETS:
                if k not in x.index:
                    cells.append('–')
                    continue
                r = x.loc[k]
                txt = f'{_f(r["median"], ".2f")} / {_f(r["min"], ".2f")}'
                if int(r['n_splits']) == 1:
                    txt = f'{_f(r["e136"], ".2f")} (one split)'
                if _truth(r['decided']):
                    txt = f'**{txt}**'
                if _truth(r['no_reading']):
                    txt += ' nr'
                cells.append(txt)
                for col in ('median', 'min', 'e136', 'range_all_flies', 'spread'):
                    V.add(k, st, f'consistency_{col}' if col in ('median', 'min', 'e136') else col, r[col],
                          condition=v, half='splits')
                V.add(k, st, 'decided', float(_truth(r['decided'])), condition=v, half='splits')
            L.append(f'| {STAGE_NAMES[st]} | {VIEW_NAMES[v]} | ' + ' | '.join(cells) + ' |')
    n = s[['set', 'read_once', 'n_models']].drop_duplicates('set').set_index('set')
    L += ['', 'The sets: ' + '; '.join(f'{SET_NAMES[k]} {int(n.loc[k, "n_models"])} models'
                                       for k in SETS if k in n.index) + '.']
    return L + ['']


def sharing_table(d: dict, V: W.Values) -> List[str]:
    """The sharing parts' errors split into level, spread and order against each half of E136's split: the consistency
    of each term over the set's models and its median against each half."""
    sh = d['sharing_split']
    if not len(sh):
        return []
    L = ["Sharing between odors split into level, spread and order (the strength split's terms on each sharing part's "
         "entries, in the part's one-fly unit), against each half of E136's split: the consistency of each term over the "
         "set's models, and the median of the order term against half a and half b.", '',
         '| stage | part | set | ρ level | ρ spread | ρ order | order, median against a / b |',
         '|---|---|---|---:|---:|---:|---:|']
    for st in STAGES:
        for part in ('sharing/who', 'sharing/how_strongly'):
            for k in ('A1', 'X', 'F', 'W', 'pool_S', 'pool_T', 'settings'):
                x = sh[(sh.stage == st) & (sh.part == part) & (sh.set == k)].set_index('term')
                if not len(x):
                    continue
                o = x.loc['order'] if 'order' in x.index else None
                L.append(f'| {STAGE_NAMES[st]} | {part.split("/")[1].replace("_", " ")} | {SET_NAMES[k]} | ' +
                         ' | '.join(_f(x.rho.get(t), '.2f') for t in ('level', 'spread', 'order')) +
                         f' | {_f(o.median_a if o is not None else None)} / {_f(o.median_b if o is not None else None)} |')
                for t in ('level', 'spread', 'order'):
                    V.add(k, st, f'sharing_{t}_rho', x.rho.get(t), condition=part, half='e136')
    return L + ['']


# ---------------------------------------------------------------------------
# design 5: the losses

def losses_table(d: dict, V: W.Values) -> List[str]:
    """L0 to L6: the consistency on the joint splits (median / minimum, decided) per cell and pooled set."""
    s = d['losses_summary']
    if not len(s):
        return []
    sets = [k for k in ('A1', 'X', 'F', 'W', 'pool_S', 'pool_T') if k in set(s['set'])]
    L = ["The losses' consistency between the halves on the joint splits: median / minimum, E136's split in brackets; "
         "**bold** where it \"can be decided\" (0.8 or more at the median, 0.5 or more on every joint split). Each set "
         "is read by E136's own reading, on the nine joint splits it saved (E136's split at the PN stages with each of "
         "the claws' and somata's three).", '',
         '| loss | what it reads | ' + ' | '.join(SET_NAMES[k] for k in sets) + ' | can be decided |',
         '|---|---|' + '---:|' * len(sets) + '---|']
    cl = (d['classes'] or {}).get('losses') or {}
    for Lname, (desc, _) in M.LOSSES.items():
        x = s[s.loss == Lname].set_index('set')
        if not len(x):
            continue
        cells = []
        for k in sets:
            if k not in x.index:
                cells.append('–')
                continue
            r = x.loc[k]
            t = f'{_f(r["median"], ".2f")} / {_f(r["min"], ".2f")} ({_f(r["e136"], ".2f")})'
            cells.append(f'**{t}**' if _truth(r['decided']) else t)
            for col in ('median', 'min', 'e136'):
                V.add(k, 'megamat', f'consistency_{col}', r[col], condition=Lname, half='joint')
            V.add(k, 'megamat', 'n_joint', r['n_joint'], condition=Lname, half='joint')
        L.append(f'| {Lname} | {desc} | ' + ' | '.join(cells) + f" | {_yn((cl.get(Lname) or {}).get('decided'))} |")
    return L + ['']


def picks_table(d: dict, V: W.Values) -> List[str]:
    """Each cell's pick under each loss on its fit half of E136's split: the evaluation, its four factors, its regret
    against the other half, whether it transfers and megamat's part of "holds" there; the joint splits' median regret,
    share transferring and share holding beside."""
    p, s = d['losses_picks'], d['losses_summary']
    if not len(p):
        return []
    L = ['Each cell\'s pick (its held evaluation with the lowest loss against its fit half of E136\'s split), its four '
         'factors (APL→bouton, KC–KC, the additive term and the gain on pn_inhsc as multiples of the family\'s; '
         'pn_inhadd absolute), its regret against the other half ((the pick\'s loss minus that half\'s lowest) over '
         '(evaluation 0\'s minus that lowest); it "transfers" at 0.25 or less) and megamat\'s part of "holds" at the pick '
         'against the other half; over the nine joint splits, the median regret and the shares transferring and holding.',
         '', '| cell | loss | pick | APL→bouton | KC–KC | additive term | pn_inhadd | gain on pn_inhsc | regret | '
         'transfers | holds part | median regret, joint | share transferring | share holding |',
         '|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---:|---:|']
    for c in CELLS:
        for Lname in M.LOSSES:
            x = p[(p.cell == c) & (p.loss == Lname)]
            if not len(x):
                continue
            r = x.iloc[0]
            q = s[(s['set'] == c) & (s.loss == Lname)]
            q = q.iloc[0] if len(q) else pd.Series(dtype=float)
            L.append(f'| {CELL_NAMES[c]} | {Lname} | {int(r["pick"])} | ×{_f(r.aplpn)} | ×{_f(r.kckc, ".2f")} | '
                     f'×{_f(r.inhadd, ".1f")} | {_f(r.pn_inhadd, ".0f")} | ×{_f(r.inhsc, ".2f")} | {_f(r.regret)} | '
                     f'{_yn(_truth(r.transfers))} | {_yn(_truth(r.holds_part_other_half))} | '
                     f'{_f(q.get("regret_median"))} | {_f(q.get("share_transfers"), ".2f")} | '
                     f'{_f(q.get("share_holds_part"), ".2f")} |')
            for k in ('pick', 'aplpn', 'kckc', 'inhadd', 'pn_inhadd', 'inhsc', 'regret'):
                V.add(c, 'megamat', f'pick_{k}', r[k], condition=Lname, half='e136')
            V.add(c, 'megamat', 'holds_part_at_pick', float(_truth(r.holds_part_other_half)), condition=Lname,
                  half='e136')
            for k in ('regret_median', 'share_transfers', 'share_holds_part'):
                V.add(c, 'megamat', k, q.get(k), condition=Lname, half='joint')
    return L + ['']


def _at(h: pd.DataFrame, i: int, c: str) -> Dict[str, object]:
    """An evaluation's per-wiring scales, both response rates and their ratio from a fit history."""
    ws = M.get_cell(c).wirings
    r = h[h['eval'] == int(i)].iloc[0]
    out = {k: [float(r[f'{k}_{w}']) for w in ws if f'{k}_{w}' in r.index] for k in SCALES}
    out.update(rate_on=float(r['rate_on']), rate_off=float(r['rate_off']), rate_ratio=float(r['rate_ratio']))
    return out


def calyx_table(d: dict, V: W.Values) -> List[str]:
    """Where L1's picks put the calyx beside L0's: per cell, evaluation 0 (the family), the L0 pick and the L1 pick, with
    their four factors and, from the fit history, the per-wiring somatic thresholds and absolute KC–KC, APL→KC and
    APL→bouton scales (range over the four searched wirings), and both megamat response rates (mean over the wirings)
    beside the flies'."""
    p = d['losses_picks']
    if not len(p) or not d['hist']:
        return []
    L = ["Where the picks of strength and sharing (L1) put the calyx beside those of strength alone (L0), on each cell's "
         "fit half of E136's split. Factors as multiples of the family's (evaluation 0); the absolute scales are the "
         "range over the cell's four searched wirings (from its fit history: each wiring's threshold and APL→KC scale "
         "are set by the hold); the response rates are the mean over the wirings (two or more spikes), the flies' "
         f"megamat rate with APL on {_f(d['flies_rate'], '.4f')}.", '',
         '| cell | point | eval | APL→bouton | KC–KC | additive term | ' +
         ' | '.join(SCALE_NAMES[k] for k in SCALES) + ' | rate, APL on | APL off | ratio |',
         '|---|---|---:|---:|---:|---:|' + '---:|' * len(SCALES) + '---:|---:|---:|']
    V.add('flies', 'kc_soma', 'rate_on', d['flies_rate'], half='all')
    for c in CELLS:
        h = d['hist'].get(c)
        if h is None:
            continue
        pts = [('evaluation 0 (the family)', 0)]
        for Lname in ('L0', 'L1'):
            x = p[(p.cell == c) & (p.loss == Lname)]
            if len(x):
                pts.append((f'the {Lname} pick', int(x.iloc[0]['pick'])))
        for lab, i in pts:
            a = _at(h, i, c)
            row = h[h['eval'] == i].iloc[0]
            fac = {f: float(np.exp(row[f'z_{f}'])) for f in FACTORS if f'z_{f}' in row.index}
            L.append(f'| {CELL_NAMES[c]} | {lab} | {i} | ×{_f(fac.get("aplpn"))} | ×{_f(fac.get("kckc"), ".2f")} | '
                     f'×{_f(fac.get("inhadd"), ".1f")} | ' +
                     ' | '.join(_rng(a[k], SCALE_FMT[k]) for k in SCALES) +
                     f' | {_f(a["rate_on"], ".4f")} | {_f(a["rate_off"], ".4f")} | {_f(a["rate_ratio"], ".2f")} |')
            key = f'{c}:e{i}'
            for k in SCALES:
                V.add(key, 'kc_soma', f'{k}_min', min(a[k]) if a[k] else None, half='searched')
                V.add(key, 'kc_soma', f'{k}_max', max(a[k]) if a[k] else None, half='searched')
            for k in ('rate_on', 'rate_off', 'rate_ratio'):
                V.add(key, 'kc_soma', k, a[k], half='searched')
    return L + ['']


# ---------------------------------------------------------------------------
# design 6: the claws

def claws_table(d: dict, V: W.Values) -> List[str]:
    """The classed points' claws and somata on strength against the family where each cell is read: against each half of
    E136's split (E136's own reading, as scored; read once beside) and against each claw (soma) fly alone."""
    c = d['claws']
    if not len(c):
        return []
    L = ["Each classed point's change from the family on strength by odor (mean over the cell's four read wirings of the "
         "per-wiring changes; one-fly units; lower is better), against each half of E136's split and against each fly "
         "alone; \"worse\" marks a rise above the stage's spread as scored (E136's: claws 0.017, somata 0.057).", '']
    for st in M.CALYX_STAGES:
        x = c[c.stage == st]
        refs = ['e136_a', 'e136_b'] + sorted({r for r in x.reference if str(r).startswith('fly')})
        L += [f'**{STAGE_NAMES[st]}**', '', '| cell | reading | ' + ' | '.join(
            {'e136_a': 'half a', 'e136_b': 'half b'}.get(r, r.replace('fly', 'fly ')) for r in refs) + ' |',
            '|---|---|' + '---:|' * len(refs)]
        for cell in CELLS:
            for rd in ("E136's own", 'read once'):
                y = x[(x.cell == cell) & (x.reading == rd)].set_index('reference')
                if not len(y):
                    continue
                cells = []
                for r in refs:
                    if r not in y.index:
                        cells.append('–')
                        continue
                    t = _f(y.loc[r, 'change'], '+.3f')
                    cells.append(t + (' worse' if _truth(y.loc[r, 'worse_than_spread']) else ''))
                    V.add(cell, st, 'classed_change', y.loc[r, 'change'], condition=rd, half=r)
                L.append(f'| {CELL_NAMES[cell]} | {rd} | ' + ' | '.join(cells) + ' |')
        L.append('')
    cl = (d['classes'] or {}).get('claws') or {}
    if cl:
        L += [f"The claws' class: the read half's: {_yn(cl.get('read_half'))} (the half: {cl.get('half')}); "
              f"{cl.get('kind')}; a fly alone gives the sign in all four cells: "
              + ', '.join(f"{k} {_yn(v)}" for k, v in (cl.get('fly_sign') or {}).items()) + '.', '']
    return L


def claws_evals_table(d: dict, V: W.Values) -> List[str]:
    """Each cell's held evaluations at the claws and somata: the share whose strength error rises from evaluation 0 by
    more than the spread, and the median change, against each reference."""
    e = d['claws_evals']
    if not len(e):
        return []
    refs = [r for r in dict.fromkeys(e.reference) if e.loc[e.reference == r, 'share_worse'].notna().any()]
    L = ["Each cell's held evaluations: the share whose strength error rises from evaluation 0 (the family) by more than "
         "the stage's spread, against each reference (E136's halves and each split's halves, by E136's own reading; "
         "a fly alone has no reading of its own, since no cell is read once).", '', '| stage | cell | ' + ' | '.join(refs) + ' |', '|---|---|' + '---:|' * len(refs)]
    for st in M.CALYX_STAGES:
        for c in CELLS:
            y = e[(e.stage == st) & (e.cell == c)].set_index('reference')
            if not len(y):
                continue
            L.append(f'| {STAGE_NAMES[st]} | {c} | ' + ' | '.join(_f(y.share_worse.get(r), '.2f') for r in refs) + ' |')
            for r in refs:
                V.add(c, st, 'share_evals_worse', y.share_worse.get(r), condition='strength', half=r)
    return L + ['']


def odor_leads_table(d: dict, V: W.Values) -> List[str]:
    """The odors on which the halves of the claw flies differ most (leads, never results)."""
    o = d['claws_odor_leads']
    if not len(o):
        return []
    L = ["Leads only: per split of the four claw flies, the five odors whose strength parts differ most between the "
         "halves (half a's mean minus half b's, in the part's one-fly unit).", '',
         '| split | halves (flies a / b) | part | odors (a − b) |', '|---|---|---|---|']
    for (sp, part), x in o.groupby(['split', 'part'], sort=True):
        x = x.sort_values('rank')
        e136 = ' (E136\'s)' if _truth(x.e136.iloc[0]) else ''
        L.append(f"| {int(sp)}{e136} | {x.half_a.iloc[0]} / {x.half_b.iloc[0]} | {part} | " +
                 ', '.join(f'{r.odor} {_f(r.a_minus_b, "+.2f")}' for r in x.itertuples()) + ' |')
        for r in x.itertuples():
            V.add('claw_flies', 'kc_claw', f'lead_{r.odor}', r.a_minus_b, condition=part, half=f'split{int(sp)}')
    return L + ['']


# ---------------------------------------------------------------------------
# design 7: without the boutons

def no_boutons_table(d: dict, V: W.Values) -> List[str]:
    """L4 against L0 over the saved settings on the 200 joint splits (with and without E135's four), and the classes
    recomputed with the boutons left out."""
    nb = d['no_boutons'] or {}
    sm = nb.get('summary') or {}
    if not sm:
        return []
    L = ["The megamat stage error without the boutons (L4) and with them (L0) over the saved settings (read once; each "
         "setting the mean over T of each wiring's root mean square of its stages' strength errors): consistency on "
         "E143's 200 joint splits, median / minimum, E136's split beside; \"orders\" is E143's yardstick (0.8 or more at "
         "the median, 0.5 or more on every split).", '',
         '| loss | settings | median | minimum | E136\'s split | orders |', '|---|---|---:|---:|---:|---|']
    for Lname in ('L0', 'L4'):
        r = sm.get(Lname) or {}
        for suf, lab in (('', 'the 23'), ('_without_e135', "without E135's four (19)")):
            L.append(f"| {Lname} | {lab} | {_f(r.get('median' + suf))} | {_f(r.get('min' + suf))} | "
                     f"{_f(r.get('e136' + suf))} | {_yn(r.get('orders' + suf))} |")
            for k in ('median', 'min', 'e136'):
                V.add('settings' + suf, 'megamat', f'consistency_{k}', r.get(k + suf), condition=Lname, half='joint200')
    L += ['', f"E143's own joint-split table on L0 reproduced to {_f(nb.get('e143_l0_reproduced_max_abs_diff'), '.1e')}.",
          '']
    rc = nb.get('reclasses') or {}
    wb, df = (rc.get('without_boutons') or {}).get('classes') or {}, (rc.get('default') or {}).get('classes') or {}
    if wb:
        L += ["E141's and E145's classes recomputed over the three other stages (E142's reclass: the boutons out of the "
              "stage error, the stage comparisons and the guard; the spreads from the family's eight wirings against "
              "E136's held-out half), with the default reading's re-read beside, which must give the classes as scored.",
              '', '| arm | as scored | default re-read | without the boutons |', '|---|---|---|---|']
        for k, v in wb.items():
            L.append(f"| {k} | {v.get('original_cls')} | {(df.get(k) or {}).get('cls')} | {v.get('cls')} |")
        L.append('')
    return L


# ---------------------------------------------------------------------------
# design 8: beside

def beside_tables(d: dict, V: W.Values) -> List[str]:
    """Beside, classing nothing: the settings' somata views beside their held-out panels and E139's price; the views
    regressed on the four factors; E145's classes with E147's spreads; E147's populations' sharing and breadth."""
    L = []
    b = d['beside_settings']
    if len(b):
        cols = [c for c in b.columns if c.startswith('panel_') and not c.endswith('_change')]
        L += ["Megamat's somata against all the flies (read once, the mean over T) beside each setting's held-out panels' "
              "strength errors (the class tables' comparisons, the mean over T) and E139's price on the APL-effect "
              "pattern part. A description only.", '',
              '| setting | somata: strength | sharing | breadth | ' +
              ' | '.join(W.PANEL_NAMES.get(c[6:], c[6:]) for c in cols) + " | E139's price |",
              '|---|---:|---:|---:|' + '---:|' * len(cols) + '---:|']
        for _, r in b[b.role == 'setting'].iterrows():
            L.append(f"| {r['setting']} | {_f(r['soma_strength'])} | {_f(r['soma_sharing'])} | {_f(r['soma_breadth'])} | "
                     + ' | '.join(_f(r[c]) for c in cols) + f" | {_f(r['e139_price_value'], '+.4f')} |")
            for v in ('strength', 'sharing', 'breadth'):
                V.add(r['setting'], 'kc_soma', 'against_all', r[f'soma_{v}'], condition=v, half='all')
            for c in cols:
                V.add(r['setting'], 'kc_soma', 'panel_strength_error', r[c], panel=c[6:], half='T')
        sub = b[b.role == 'setting']
        if len(sub) > 2:
            from scipy import stats
            rho = {c: stats.spearmanr(sub.soma_sharing, sub[c], nan_policy='omit')[0] for c in cols}
            L += ['', 'ρ over the settings between their somata\'s sharing error and each held-out panel\'s strength '
                  'error: ' + '; '.join(f'{W.PANEL_NAMES.get(c[6:], c[6:])} {_f(v, ".2f")}' for c, v in rho.items())
                  + '.']
            for c, v in rho.items():
                V.add('settings', 'kc_soma', 'rho_sharing_vs_panel', v, panel=c[6:], half='all')
        L.append('')
    rg = d['beside_regressions']
    if len(rg):
        L += ["Each stage's views regressed on the four log factors over each cell's held evaluations (least squares, "
              "against each half of E136's split): R² against half a / half b.", '',
              '| stage | view | ' + ' | '.join(CELLS) + ' |', '|---|---|' + '---:|' * len(CELLS)]
        for st in STAGES:
            for v in _views_of(st):
                x = rg[(rg.stage == st) & (rg.view == v)]
                cells = []
                for c in CELLS:
                    y = x[x.cell == c].set_index('half')
                    cells.append(f"{_f(y.r2.get('a'), '.2f')} / {_f(y.r2.get('b'), '.2f')}" if len(y) else '–')
                    for h in ('a', 'b'):
                        V.add(c, st, 'r2_factors', y.r2.get(h) if len(y) else None, condition=v, half=h)
                L.append(f'| {STAGE_NAMES[st]} | {VIEW_NAMES[v]} | ' + ' | '.join(cells) + ' |')
        L.append('')
    e = d['beside_e147_spreads']
    if len(e):
        L += ["E145's classes recomputed with E147's spreads over model flies (E147's classes_with_spreads).", '',
              '| cell | as scored | as scored, re-read | with E147\'s spreads |', '|---|---|---|---|']
        for r in e.itertuples():
            L.append(f'| {r.arm} | {r.cls_scored} | {r.cls_as_scored} | {r.cls_new} |')
        L.append('')
    p = d['beside_e147_populations']
    if len(p):
        x = p[p.reading.astype(str) == 'all'] if 'reading' in p else p
        L += ["E147's saved populations of model flies (each ORN fly's own input) at the megamat stages, against all "
              "the flies: sharing and breadth, beside the same at the lab's averaged input.", '',
              '| population | stage | sharing | at the lab\'s input | breadth | at the lab\'s input |',
              '|---|---|---:|---:|---:|---:|']
        for r in x.itertuples():
            if str(r.stage) not in STAGES:
                continue
            L.append(f'| {r.population} | {STAGE_NAMES[r.stage]} | {_f(r.pop_view_sharing)} | '
                     f'{_f(r.labpop_view_sharing)} | {_f(r.pop_view_breadth)} | {_f(r.labpop_view_breadth)} |')
            V.add(f'E147:{r.population}', r.stage, 'pop_view_sharing', r.pop_view_sharing, half='all')
            V.add(f'E147:{r.population}', r.stage, 'pop_view_breadth', r.pop_view_breadth, half='all')
        L.append('')
    return L


# ---------------------------------------------------------------------------
# the prediction and the outcomes

def prediction_table(d: dict, V: W.Values) -> List[str]:
    """The prediction written before any view (prediction.json) against what was scored: per stage and view (predicted
    decided over every cell from E136's split, against decided over all four cells on every split), and the planner's
    chances against the classes."""
    p = d['prediction'] or {}
    s = d['views_decided_summary']
    cl = d['classes'] or {}
    if not p:
        return []
    L = [f"The prediction (prediction.json, written {p.get('written')} on {p.get('host')}, before any view was read): a "
         "view predicted decided over every cell where its consistency on E136's split over the held evaluations was "
         "0.8 or more in all four cells; scored: decided over each of the four cells on every split.", '',
         '| stage | view | ' + ' | '.join(f'{c} on E136\'s split' for c in CELLS) + ' | predicted | scored |',
         '|---|---|' + '---:|' * len(CELLS) + '---|---|']
    for k, v in (p.get('per_stage_view') or {}).items():
        st, view = k.split('/')
        each = v.get('each') or {}
        x = s[(s.stage == st) & (s.view == view)].set_index('set') if len(s) else pd.DataFrame()
        scored = all(c in x.index and _truth(x.loc[c, 'decided']) for c in CELLS) if len(x) else None
        L.append(f'| {STAGE_NAMES[st]} | {VIEW_NAMES[view]} | ' + ' | '.join(_f(each.get(c), '.2f') for c in CELLS) +
                 f" | {_yn(v.get('predicted_decided_over_every_cell'))} | {_yn(scored) if scored is not None else '–'} |")
    rc, ls, claws, nb = cl.get('reads_calyx') or {}, cl.get('loss_with_sharing_decided') or [], cl.get('claws') or {}, \
        cl.get('no_boutons') or {}
    scored = {'reads_calyx_over_every_fit_on_all_three_splits': bool(rc.get('decided_over_the_fits')),
              'reads_calyx_over_the_nine_settings_too': bool(rc.get('holds')),
              'a_loss_with_sharing_can_be_decided': bool(ls),
              'claws_cost_is_the_read_half': bool(claws.get('read_half')),
              'claws_cost_is_one_fly': claws.get('kind') == "one fly's",
              'without_boutons_orders_the_settings': bool(nb.get('orders'))}
    names = {'reads_calyx_over_every_fit_on_all_three_splits': 'sharing or breadth decided at the claws or the somata '
             'over every fit\'s evaluations', 'reads_calyx_over_the_nine_settings_too': '… and over the nine settings '
             'at pn_inhadd 1000', 'a_loss_with_sharing_can_be_decided': 'a loss with sharing can be decided',
             'claws_cost_is_the_read_half': "the claws' cost is the read half's", 'claws_cost_is_one_fly': '… and '
             "one fly's", 'without_boutons_orders_the_settings': 'without the boutons the megamat stage error orders '
             'the settings'}
    L += ['', "The planner's chances (the entry's, in prediction.json) and what was scored (classes.json).", '',
          "| event | the planner's chance | scored |", '|---|---|---|']
    for k, n in names.items():
        L.append(f'| {n} | {p.get(k, "–")} | {_yn(scored[k])} |')
        V.add('experiment', 'megamat', f'event_{k}', float(scored[k]), condition='score', half='all')
    return L + ['']


def outcomes_table(d: dict, V: W.Values) -> List[str]:
    """The score (classes.json): each cell's identities and reading, the classes and the outcome letters."""
    cl = d['classes'] or {}
    if not cl:
        return []
    L = ['The score (classes.json' + (f", {cl.get('time')} on {cl.get('host')}" if cl.get('time') else '') + ').', '',
         '| cell | identities pass | read once | read by |', '|---|---|---|---|']
    for c, r in (cl.get('cells') or {}).items():
        L.append(f"| {c} | {_yn(r.get('identities'))} | {_yn(r.get('read_once'))} | {r.get('read_by')} |")
    rc = cl.get('reads_calyx') or {}
    L += ['', '| class | scored |', '|---|---|',
          f"| megamat reads the calyx beside strength | {_yn(rc.get('holds'))} ({', '.join(rc.get('stage_views') or []) or 'none'}) |",
          f"| decided over every fit's evaluations and both pools | {', '.join(rc.get('decided_over_the_fits') or []) or 'none'} |",
          f"| at 0.8 at the median over the fits' sets | {', '.join(rc.get('median_over_the_fits') or []) or 'none'} |",
          f"| losses with sharing that can be decided | {', '.join(cl.get('loss_with_sharing_decided') or []) or 'none'} |",
          f"| of them, picks that transfer in every cell | {', '.join(cl.get('picks_transfer') or []) or 'none'} |",
          f"| the claws' cost | the read half's: {_yn((cl.get('claws') or {}).get('read_half'))}, half "
          f"{(cl.get('claws') or {}).get('half')}, {(cl.get('claws') or {}).get('kind')} |",
          f"| without the boutons the settings are ordered | {_yn((cl.get('no_boutons') or {}).get('orders'))} |",
          f"| cells with no class | {', '.join(cl.get('no_class') or []) or 'none'} |", '',
          'Outcomes: ' + ', '.join(cl.get('outcomes') or []) + '.', '']
    for k in cl.get('outcomes') or []:
        L.append(f'- ({k}) {(cl.get("outcome_text") or {}).get(k, M.OUTCOMES.get(k, ""))}')
        V.add('experiment', 'megamat', f'outcome_{k}', 1.0, condition='score', half='all')
    return L + ['']


def main(score_dir: Optional[Path] = None, histories: bool = True) -> Path:
    """Writes writeup_tables.md and writeup_values.csv into `score_dir` (default E148's results folder); `histories`
    False leaves out the fit histories and E134's files (planted folders in tests)."""
    score_dir = Path(score_dir) if score_dir is not None else M.results_dir()
    d = load(score_dir, histories)
    V = W.Values()
    ps = d['prov_score']
    head = ['# E148 write-up tables', '',
            f"From results/E148_megamat_other_views/ (the score at queue {str(ps.get('queue_commit'))[:7]}, "
            f"{ps.get('time')} on {ps.get('host')}); written by analyses/E148_writeup_tables.py. One-fly units, lower is "
            'better, unless a table says otherwise.', '']
    parts = [('The room at the four megamat stages', stage_table), ("Step 0: the planner's probe", probe_table),
             ('Step 0: the identities', identities_table), ("Step 0: E143's read-once check", read_once_table),
             ('The views', views_read_table), ('The room: against the halves', room_table),
             ('Which views the flies decide', decided_table), ('Sharing: level, spread and order', sharing_table),
             ('Which loss the flies decide', losses_table), ('The picks', picks_table),
             ("Where L1's picks put the calyx beside L0's", calyx_table), ("The claws' cost", claws_table),
             ("The claws' cost: the evaluations", claws_evals_table), ('The claw flies: odor leads', odor_leads_table),
             ('Without the boutons', no_boutons_table), ('Beside', beside_tables),
             ('The prediction', prediction_table), ('The outcomes', outcomes_table)]
    lines = head
    for title, fn in parts:
        body = fn(d, V)
        if body:
            lines += [f'## {title}', ''] + body
    out = score_dir / 'writeup_tables.md'
    out.write_text('\n'.join(lines) + '\n')
    V.frame().to_csv(score_dir / 'writeup_values.csv', index=False)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--score-dir', type=Path, default=None)
    a = ap.parse_args()
    print(main(a.score_dir))
