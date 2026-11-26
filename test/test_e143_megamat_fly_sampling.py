"""Tests of experiments/E143_megamat_fly_sampling.py on planted values (no saved run is read)."""
import numpy as np
import pandas as pd
import pytest

from experiments import E143_megamat_fly_sampling as E
from kcq import population_fit as P


def test_stage_splits_counts():
    assert [len(E.stage_splits(n)) for n in (7, 5, 4)] == [35, 10, 3]
    for n in (4, 5, 7):
        for a, b in E.stage_splits(n):
            assert len(a) == n // 2 and sorted(a + b) == list(range(n))
    assert all(0 in a for a, _ in E.stage_splits(4))


def test_e136_split_is_split_halves_seed0():
    for n in (4, 5, 7):
        a, b = P.split_halves(n, 1, 0)[0]
        assert E.e136_split(n) == (tuple(a), tuple(b))
        i = E.split_index(n, E.e136_split(n))
        x, y = E.stage_splits(n)[i]
        assert {x, y} == {tuple(a), tuple(b)}


def test_joint_splits_distinct_and_e136_last():
    js = E.joint_splits(dict(pn_dendrite=7, pn_bouton=5, kc_claw=4, kc_soma=4))
    assert len(js) == E.N_JOINT + 1 and js['e136'].sum() == 1 and bool(js['e136'].iloc[-1])
    drawn = js[~js['e136']][list(E.LOSS_STAGES)]
    assert not drawn.duplicated().any()
    assert drawn['pn_dendrite'].between(0, 34).all() and drawn['kc_soma'].between(0, 2).all()
    # seeded: the same draw twice
    assert js.equals(E.joint_splits(dict(pn_dendrite=7, pn_bouton=5, kc_claw=4, kc_soma=4)))


def _planted_pops(n_flies=5, n_units=40, odors=('a', 'b', 'c', 'd', 'e', 'f'), seed=0, scale=1.0):
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n_flies):
        v = pd.DataFrame(rng.gamma(2.0, scale, size=(n_units, len(odors))), columns=list(odors))
        out.append(P.Population.from_criterion(v, 2.0, name=f'f{i}'))
    return out


def _side_of(flies, models, key='kc_soma/megamat/on'):
    """A Side from planted populations, read as E143 reads saved runs (reference_views against all flies)."""
    ref = P.fly_reference(flies, None, (0, 1, 2, 4))
    frows = [pd.DataFrame(dict(part=q, fly=f, entry=list(df.columns), value=df.loc[f].to_numpy(dtype=float),
                               unit=float(ref.one_fly[q]))) for q, df in ref.parts.items() for f in df.index]
    mrows = []
    for (name, w), pop in models.items():
        for q, x in P.reference_views(pop, ref).items():
            mrows.append(pd.DataFrame(dict(setting=name, wiring=w, part=q, entry=list(x.index),
                                           value=x.to_numpy(dtype=float))))
    E._ENTRIES[(key, 'strength/rate')] = list(ref.parts['strength/rate'].columns)
    for q, df in ref.parts.items():
        E._ENTRIES[(key, q)] = list(df.columns)
    return E.side_from_tables(key, pd.concat(frows), pd.concat(mrows)), ref


def test_side_errors_match_population_fit_on_halves_and_resamples():
    flies = _planted_pops()
    models = {(f's{k}', w): _planted_pops(1, 30, seed=10 + 4 * k + j, scale=1.0 + 0.2 * k)[0]
              for k in range(3) for j, w in enumerate(E.T_WIRINGS)}
    side, ref = _side_of(flies, models)
    halves = [(0, 1), (2, 3, 4)]
    W = E.half_weights(5, halves)
    for views, col in ((('strength',), 'error'), (None, 'error')):
        err = E.side_errors(side, W, views, col)
        for i, h in enumerate(halves):
            sub = P.subset_reference(ref, list(h), keep_units=True)
            if views is not None:
                sub = P.FlyReference(parts={q: d for q, d in sub.parts.items() if q.startswith('strength/')},
                                     mean={q: d for q, d in sub.mean.items() if q.startswith('strength/')},
                                     one_fly={q: d for q, d in sub.one_fly.items() if q.startswith('strength/')},
                                     flies=sub.flies, odors=sub.odors, breadth_edges=sub.breadth_edges)
            for r, (name, w) in enumerate(side.rows):
                want = P.score(models[(name, w)], sub).error
                assert err[i, r] == pytest.approx(want, rel=1e-10)
    # a resample with replacement: the weighted mean equals the mean of the drawn rows
    draws = np.array([[0, 0, 3, 4, 4]])
    Wb = E.boot_weights(5, draws)
    err = E.side_errors(side, Wb)
    rq = ref.parts['strength/rate'].to_numpy()[draws[0]].mean(axis=0)
    mq = ref.parts['strength/mean'].to_numpy()[draws[0]].mean(axis=0)
    name, w = side.rows[0]
    v = P.reference_views(models[(name, w)], ref)
    e_r = np.sqrt(np.mean((v['strength/rate'].to_numpy() - rq) ** 2)) / ref.one_fly['strength/rate']
    e_m = np.sqrt(np.mean((v['strength/mean'].to_numpy() - mq) ** 2)) / ref.one_fly['strength/mean']
    assert err[0, 0] == pytest.approx(np.sqrt((e_r ** 2 + e_m ** 2) / 2), rel=1e-10)
    vals = E.setting_values(E.side_errors(side, W), side)
    assert list(vals.columns) == ['s0', 's1', 's2'] and vals.notna().all().all()


def test_drop_entry_equals_odor_entry_at_the_flies():
    flies = _planted_pops()
    pop = _planted_pops(1, 30, seed=99)[0]
    side, ref = _side_of(flies, {('s', w): pop for w in E.T_WIRINGS})
    W = E.half_weights(5, [(0, 1, 2, 3, 4)])
    got = E.side_errors(side, W, drop=('c',))[0, 0]
    v = P.reference_views(pop, ref)
    for q in ('strength/rate', 'strength/mean'):
        v[q] = v[q].copy()
        v[q]['c'] = ref.mean[q]['c']
    want = P.score_views({q: v[q] for q in ('strength/rate', 'strength/mean')},
                         P.FlyReference(parts={q: ref.parts[q] for q in ('strength/rate', 'strength/mean')},
                                        mean={q: ref.mean[q] for q in ('strength/rate', 'strength/mean')},
                                        one_fly={q: ref.one_fly[q] for q in ('strength/rate', 'strength/mean')},
                                        flies=ref.flies, odors=ref.odors)).error
    assert got == pytest.approx(want, rel=1e-10)
