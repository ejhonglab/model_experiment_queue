"""Tests of experiments/E140_aa_claw_soma_step.py (E140, aa at the step from the claws to the somata, step 0): the lock
(committed, its recordings, exclusions, halves and classes as the module fixes them; the guard refuses an uncommitted or
missing copy), the concentration and the ranks and thirds, the per-odor statistics against hal's views, the class
functions and the outcomes on planted tables, the leave-aa-out fit, and a known value: E138's stage table entries of aa
(the flies at the claws 0.92 and the somata 1.35, the family held 1.18 and 0.65) reproduced from the same inputs."""
import json

import numpy as np
import pandas as pd
import pytest

from experiments import E140_aa_claw_soma_step as E
from kcq import harness as H
from kcq import population_fit as P

ODORS = [f'o{i}' for i in range(16)] + ['aa']


def _shared_ok() -> bool:
    try:
        return (H.shared_dir() / 'data' / 'stage_populations' / 'populations').is_dir()
    except FileNotFoundError:
        return False


def planted_pop(name, seed, n_units=400, scale=None, odors=ODORS, criterion=1.5):
    rng = np.random.default_rng(seed)
    sc = np.ones(len(odors)) if scale is None else np.asarray(scale, float)
    v = pd.DataFrame(rng.gamma(1.0, 1.0, size=(n_units, len(odors))) * sc[None, :],
                     index=[str(i) for i in range(n_units)], columns=odors)
    return P.Population.from_criterion(v, criterion, name=name)


# ---------------------------------------------------------------------------
# the lock

def test_lock_is_committed_and_fixes_the_recordings():
    assert E.LOCK_COPY.exists()
    assert E._committed(E.LOCK_COPY)
    lock = json.loads(E.LOCK_COPY.read_text())
    assert lock['experiment'] == E.EXPERIMENT and lock['odor'] == 'aa' and len(lock['odors']) == 17
    claws = lock['recordings']['claw_recordings']['recordings']
    assert len(claws) == 5
    assert [c['cohort'] for c in claws].count('wild_type') == 4 and [c['cohort'] for c in claws].count('same_fly') == 1
    assert not any(f.startswith('2022-10-21__fly02') for f in lock['recordings']['pn_bouton__wild_type']['flies'])
    assert len(lock['recordings']['kc_soma__wild_type']['flies']) == 4
    # Remy's older megamat flies are never an input
    assert '2024-11-12' not in json.dumps(lock['inputs'])
    # every input that is a file carries its sha256
    for k in ('e134_lock', 'e138_classes_lock', 'frames_reduced', 'olfdb_dataset_1', 'hallem', 'megamat_input'):
        assert len(lock['inputs'][k]['sha256']) == 64
    for f in lock['inputs']['soma_frames'].values():
        assert len(f['trials']['sha256']) == 64 and len(f['respvec']['sha256']) == 64


def test_lock_classes_and_definitions_are_the_modules():
    lock = json.loads(E.LOCK_COPY.read_text())
    assert lock['classes'] == json.loads(json.dumps(E.CLASSES))
    c = lock['classes']
    assert (c['i']['flies_min'], c['i']['flies_n'], c['i']['model_min'], c['i']['model_n']) == (4, 5, 6, 8)
    assert c['ii']['top'] == 3 and c['ii']['part'] == 'sharing/who'
    assert c['iii']['min_wirings'] == 3 and c['iii']['wirings'] == list(E.T_WIRINGS)
    d = lock['definitions']
    assert d['concentration']['top_share'] == 0.1 and d['ranks']['third'] == 6
    assert d['soma_halves'] == E.soma_halves(lock['recordings']['kc_soma__wild_type']['flies'])
    counted = [r['name'] for r in d['soma_readings'] if r['counts']]
    assert counted == ['lock', 'presentation_1', 'presentation_2', 'presentation_3', 'window_early', 'window_late',
                       'criterion_apl16', 'criterion_remy']


def test_require_lock_refuses(tmp_path):
    with pytest.raises(PermissionError):
        E.require_lock(copy=tmp_path / 'E140_lock.json')
    p = tmp_path / 'E140_lock.json'
    p.write_text('{}\n')
    with pytest.raises(PermissionError):
        E.require_lock(copy=p)


def test_part_lock_refuses_to_overwrite(monkeypatch):
    def never():
        raise AssertionError('the lock would be rewritten')
    monkeypatch.setattr(E, 'lock_contents', never)
    with pytest.raises(PermissionError):
        E.part_lock()


def test_half_windows_partition_the_window():
    t = np.round(np.arange(-10, 45, 0.25), 2)
    full = (t >= E.FULL_WINDOW[0]) & (t <= E.FULL_WINDOW[1])
    early = (t >= E.HALF_WINDOWS['early'][0]) & (t <= E.HALF_WINDOWS['early'][1])
    late = (t >= E.HALF_WINDOWS['late'][0]) & (t <= E.HALF_WINDOWS['late'][1])
    assert np.array_equal(full, early | late) and not (early & late).any()
    assert (early.sum(), late.sum()) == (12, 13)


# ---------------------------------------------------------------------------
# the statistics

def test_concentration():
    assert np.isclose(E.concentration(np.ones(50)), 0.1)
    v = np.r_[10.0, np.ones(14)]                         # n 15: the top 1.5 units
    assert np.isclose(E.concentration(v), (10 + 0.5) / 24)
    assert np.isnan(E.concentration(np.ones(4)))         # fewer than MIN_RESPONDERS
    assert np.isnan(E.concentration(np.zeros(10)))
    assert np.isclose(E.concentration(np.r_[np.ones(10), np.nan]), 0.1)
    assert np.isclose(E.responders_cv(np.array([1.0, 3.0, 1.0, 3.0, 2.0])), np.std([1, 3, 1, 3, 2]) / 2)


def test_ranks_and_thirds():
    r = E.ranks(np.arange(17, dtype=float))
    assert list(r) == list(range(1, 18))
    th = [E.third(v, 17) for v in r]
    assert th[:6] == ['lowest'] * 6 and th[6:11] == ['middle'] * 5 and th[11:] == ['highest'] * 6
    r = E.ranks(np.r_[np.nan, np.arange(16, dtype=float)])
    assert np.isnan(r.iloc[0]) and r.iloc[-1] == 16
    assert E.third(np.nan, 16) == '' and E.third(11, 16) == 'highest' and E.third(10, 16) == 'middle'


def test_odor_stats_against_hals_views():
    sc = np.linspace(0.5, 2.0, len(ODORS))
    pop = planted_pop('f0', 1, scale=sc)
    s = E.odor_stats(pop, ODORS).set_index('odor')
    v = P.population_views(pop, ODORS)
    assert np.allclose(s.rate, v['strength/rate'].reindex(ODORS))
    assert np.allclose(s['mean'], v['strength/mean'].reindex(ODORS))
    resp = pop.values['aa'][pop.responds['aa'] > 0.5]
    assert np.isclose(s.loc['aa', 'concentration'], E.concentration(resp.to_numpy()))
    assert s.loc['aa', 'rate_third'] == 'highest' and s.loc['o0', 'rate_third'] == 'lowest'
    m = E.mean_rows({'a': E.odor_stats(pop, ODORS), 'b': E.odor_stats(planted_pop('f1', 2, scale=sc), ODORS)})
    ref = P.fly_reference([pop, planted_pop('f1', 2, scale=sc)], ODORS, parts=list(E.STRENGTH_PARTS))
    assert np.allclose(m.set_index('odor').rate, ref.mean['strength/rate'].reindex(ODORS))


def test_sharing_rows():
    pops = [planted_pop(f'f{i}', i) for i in range(2)]
    parts = [E.sharing_parts(p, ODORS) for p in pops]
    mp = E.mean_parts(parts)
    row = E.odor_row(mp['sharing/who'])
    assert len(row) == 16 and set(row.index) == set(ODORS) - {'aa'}
    s = E.sharing_summary(mp, ODORS)
    x = s[(s.part == 'sharing/who') & (s.odor == 'aa')]
    assert np.isclose(x.mean_distance.iloc[0], row.mean())
    assert int(s[s.part == 'sharing/who'].among_top.sum()) == 3


def test_class_functions_and_outcomes():
    def rows(n, rate, conc, key):
        return pd.DataFrame({key: [f'r{i}' for i in range(n)], 'rate_third': rate, 'concentration_third': conc})
    f = rows(5, ['lowest'] * 4 + ['middle'], ['highest'] * 5, 'recording')
    m = rows(8, ['highest'] * 6 + ['middle'] * 2, ['lowest'] * 8, 'wiring')
    assert E.class_i(f, m)['holds']
    m2 = rows(8, ['highest'] * 5 + ['middle'] * 3, ['lowest'] * 8, 'wiring')
    c = E.class_i(f, m2)
    assert not c['holds'] and c['flies_part'] and not c['model_part']
    s = pd.DataFrame(dict(part='sharing/who', odor=ODORS, mean_distance=np.arange(17.0),
                          rank=np.arange(1.0, 18.0), n_ranked=17, among_top=[False] * 14 + [True] * 3))
    assert E.class_ii({'fit': s, 'held': s})['holds']
    s2 = s.assign(mean_distance=s.mean_distance[::-1].to_numpy(), rank=s['rank'][::-1].to_numpy(),
                  among_top=s.among_top[::-1].to_numpy())
    assert not E.class_ii({'fit': s, 'held': s2})['holds']
    rd = pd.DataFrame(dict(reading=['lock', 'beside'], counts=[True, False], rate=[1.3, 0.8], rate_rank=[16.0, 9.0],
                           rate_n_ranked=[17, 17], rate_min=[1.0, 0.6], rate_max=[1.6, 1.0], mean=[1.2, 0.9],
                           mean_rank=[15.0, 8.0]))
    c4 = E.class_iv(rd)
    assert not c4['holds'] and c4['readings']['beside']['leaves_highest_third']   # a beside reading does not count
    assert E.class_iv(rd.assign(counts=True))['holds']
    assert E.outcomes_of(False, None, False)['outcomes'] == ['c']
    assert E.outcomes_of(True, None, True)['outcomes'] == ['d', 'a or b (step 1 reads (iii))']
    assert E.outcomes_of(True, None, False)['step1_called']
    assert E.outcomes_of(True, True, False)['outcomes'] == ['b']


def test_loo_fit():
    x = pd.Series(np.arange(17.0), index=ODORS)
    y = 2 * x + 1
    y['aa'] += 5.0
    f = E.loo_fit(x, y)
    assert np.isclose(f['slope'], 2) and np.isclose(f['intercept'], 1)
    assert np.isclose(f['resid_aa'], 5.0) and f['loo_rank_aa'] == 1
    assert np.isclose(f['r'], 1.0)


def test_input_features():
    inp = pd.DataFrame({'a': [10.0, 10.0, -5.0, 0.0], 'b': [40.0, 0.0, 0.0, 0.0]}, index=list('wxyz'))
    f = E.input_features(inp).set_index('odor')
    assert np.isclose(f.loc['a', 'input_sum'], 15) and np.isclose(f.loc['a', 'input_breadth'], 2.0)
    assert np.isclose(f.loc['b', 'input_breadth'], 1.0) and f.loc['b', 'input_argmax'] == 'w'


# ---------------------------------------------------------------------------
# a known value: E138's stage table

@pytest.mark.skipif(not _shared_ok(), reason='the shared directory is not mounted')
def test_known_value_e138_stage_table():
    """The stage table's aa entries in E138's note (stages.csv): the flies at the claws 0.92 (0.919957) and the somata
    1.35 (1.353566); the family held 1.18 (1.180121) and 0.65 (0.651361), the wirings' mean views."""
    from kcq import stage_readings as R
    lock, e134 = E.require_all()
    odors = list(lock['odors'])
    for stage, want in (('kc_claw', 0.919957), ('kc_soma', 1.353566)):
        pop = E.recorded(stage, 'wild_type', lock)
        m = E.mean_rows({fp.name: E.odor_stats(fp, odors) for fp in P.fly_populations(pop)}).set_index('odor')
        assert abs(m.loc['aa', 'rate'] - want) < 5e-6
    runs, missing = E.group_runs(E.MODEL_GROUPS[0])
    if missing:
        pytest.skip(f'saved runs missing: {missing}')
    rd = E.read_runs(runs, e134, R.FlyCache(), stages=['kc_claw', 'kc_soma'])
    for stage, want in (('kc_claw', 1.180121), ('kc_soma', 0.651361)):
        d = {w: E.odor_stats(m.populations[(stage, 'megamat', 'on')], odors) for w, m in rd.items()}
        assert abs(E.mean_rows(d).set_index('odor').loc['aa', 'rate'] - want) < 5e-6
