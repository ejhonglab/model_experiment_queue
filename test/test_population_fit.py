"""Tests of kcq/population_fit.py on planted populations, whose answers are known by construction (E134)."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from kcq import population_fit as P

ODORS = [f'o{i}' for i in range(16)]
CRIT = 1.0


class Gen:
    """Planted flies: units load on a few odor factors shared by every fly, rectified (so most units are at the noise
    and the responders are positive, as in the recordings); each fly has its own odor strengths (a multiplicative
    jitter of the shared ones) and its own units and noise."""

    def __init__(self, seed=0, n_odors=16, k=4, jitter=0.15):
        rng = np.random.default_rng(seed)
        self.odors = ODORS[:n_odors]
        self.F = rng.normal(size=(k, n_odors))
        self.s = rng.uniform(0.6, 1.6, size=n_odors)
        self.jitter = jitter
        self.rng = rng

    def fly(self, name, n_units=400, named=False):
        r = self.rng
        s = self.s * np.exp(r.normal(0, self.jitter, size=len(self.s)))
        W = r.normal(size=(n_units, self.F.shape[0]))
        V = 2 * np.maximum(0, (W @ self.F) * s[None, :] * 0.7 - 0.3) + r.normal(0, 0.5, size=(n_units, len(self.s)))
        idx = [f'g{i}' for i in range(n_units)] if named else range(n_units)
        return P.Population.from_criterion(pd.DataFrame(V, index=idx, columns=self.odors), CRIT, name=name,
            named_units=named)


@pytest.fixture(scope='module')
def gen():
    return Gen(0)


@pytest.fixture(scope='module')
def flies(gen):
    return [gen.fly(f'f{i}') for i in range(11)]


@pytest.fixture(scope='module')
def ref(flies):
    return P.fly_reference(flies)


def test_breadth_bins_edges_and_pairs():
    assert P.breadth_bins((0, 1, 2, 4, 8)) == [(0, 0.0), (1, 1.0), (2, 3.0), (4, 7.0), (8, np.inf)]
    assert P.breadth_labels() == ['0', '1', '2_3', '4_7', '8plus']
    pairs = ((0, 0), (1, 1), (2, 3), (4, 7), (8, 16))
    assert P.breadth_labels(pairs) == ['0', '1', '2_3', '4_7', '8_16']
    R = pd.DataFrame(np.array([[0] * 16, [1] + [0] * 15, [1] * 3 + [0] * 13, [1] * 9 + [0] * 7]), columns=ODORS)
    c = P.breadth_counts(R)
    assert c.to_dict() == {'0': 0.25, '1': 0.25, '2_3': 0.25, '4_7': 0.0, '8plus': 0.25}
    assert P.breadth_counts(R, pairs).sum() == pytest.approx(1.0)
    with pytest.raises(ValueError):
        P.breadth_bins((1, 2))


def test_views_by_hand():
    V = pd.DataFrame([[2.0, 0.0, 3.0], [2.0, 2.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 0.0]], columns=['a', 'b', 'c'])
    pop = P.Population.from_criterion(V, 1.0)
    v = P.population_views(pop, breadth_edges=(0, 1, 2))
    rate = np.array([0.5, 0.5, 0.25])
    assert np.allclose(v['strength/rate'], rate / rate.mean())
    mean = V.mean().to_numpy()
    assert np.allclose(v['strength/mean'], mean / mean.mean())
    assert list(v['sharing/who'].index) == ['a|b', 'a|c', 'b|c']
    # who responds: a = (1,1,0,0), b = (0,1,1,0), c = (1,0,0,0)
    assert np.allclose(v['sharing/who'], [1 - 1 / 2, 1 - 1 / np.sqrt(2), 1.0])
    # how strongly: a = (2,2,0,0), b = (0,2,4,0)
    assert v['sharing/how_strongly']['a|b'] == pytest.approx(1 - 4 / (np.sqrt(8) * np.sqrt(20)))
    assert v['breadth/bins'].to_dict() == {'0': 0.25, '1': 0.25, '2plus': 0.5}
    assert 'matrix/values' not in v


def test_unit_order_and_common_scale_do_not_matter(gen, ref):
    f = gen.fly('x')
    perm = P.Population.from_criterion(f.values.sample(frac=1, random_state=1), CRIT, name='x')
    scaled = P.Population.from_criterion(f.values * 3.7, CRIT * 3.7, name='x')
    a, b, c = (P.score(p, ref).as_series() for p in (f, perm, scaled))
    pd.testing.assert_series_equal(a, b, rtol=1e-12)
    pd.testing.assert_series_equal(a, c, rtol=1e-12)


def test_scaling_values_leaves_the_mean_part(gen, ref):
    f = gen.fly('x')
    up = P.Population.from_criterion(f.values * 2.0, CRIT, name='x')
    va, vb = P.population_views(f), P.population_views(up)
    pd.testing.assert_series_equal(va['strength/mean'], vb['strength/mean'])
    # more units pass the criterion, so breadth moves up
    assert vb['breadth/bins']['0'] < va['breadth/bins']['0']


def test_one_recording_against_the_others_is_one_by_construction(flies, ref):
    loo = P.leave_one_out(flies, ref=ref)
    assert list(loo.index) == ref.flies
    for p in ref.parts:
        col = 'part_' + p.replace('/', '_')
        assert np.sqrt((loo[col] ** 2).mean()) == pytest.approx(1.0, rel=1e-12)
    # and the same from the flies directly
    pd.testing.assert_frame_equal(loo, P.leave_one_out(flies))


def test_flies_mean_scores_zero_and_pooled_flies_below_one(flies, ref):
    s = P.score_views(ref.mean, ref)
    assert s.error == 0.0 and all(v == 0.0 for v in s.parts.values())
    pooled = P.Population.from_criterion(pd.concat([f.values for f in flies], ignore_index=True), CRIT, name='pool')
    assert P.score(pooled, ref).error < 0.5


def test_a_new_fly_scores_about_one(gen, ref):
    new = [gen.fly(f'n{i}') for i in range(40)]
    sc = P.score_many(new, ref)
    rms = float(np.sqrt((sc['error'] ** 2).mean()))
    assert 0.75 < rms < 1.3
    for v in P.VIEWS:
        assert 0.6 < float(np.sqrt((sc[f'view_{v}'] ** 2).mean())) < 1.5


def test_floor_is_below_one_and_falls_with_more_flies(gen, ref):
    fl = P.floor_summary(P.floor(ref, n_boot=400, seed=1))
    assert 0.15 < fl['error'] < 0.45
    # Gaussian expectation for a mean over n flies, in one-fly units: about (n - 1) / n**1.5
    assert fl['error'] == pytest.approx(10 / 11 ** 1.5, abs=0.08)
    small = P.fly_reference([gen.fly(f's{i}') for i in range(5)])
    big = P.fly_reference([gen.fly(f'b{i}') for i in range(25)])
    assert P.floor_summary(P.floor(small, 300))['error'] > P.floor_summary(P.floor(big, 300))['error']


def test_permuted_and_random_populations_score_far(gen, ref, flies):
    f = gen.fly('x')
    own = P.score(f, ref)
    shuffled = f.values.copy()
    shuffled.columns = list(np.random.default_rng(3).permutation(f.values.columns))
    perm = P.score(P.Population.from_criterion(shuffled[ODORS], CRIT), ref)
    assert perm.error > 2 * own.error
    assert perm.views['breadth'] == pytest.approx(own.views['breadth'])   # breadth does not see odor labels
    rng = np.random.default_rng(4)
    vals = f.values.to_numpy()
    rand = pd.DataFrame(rng.permutation(vals.ravel()).reshape(vals.shape), columns=ODORS)
    r = P.score(P.Population.from_criterion(rand, CRIT), ref)
    assert r.error > 2.0
    assert r.views['sharing'] > 2.0


def test_undefined_entries_are_counted_not_scored(ref):
    V = pd.DataFrame(np.zeros((50, 16)), columns=ODORS)
    V.iloc[:10, :8] = 5.0
    s = P.score(P.Population.from_criterion(V, CRIT), ref)
    assert s.missing['strength/rate'] == 0
    assert s.missing['sharing/who'] == 120 - 28
    assert np.isfinite(s.error)
    dead = P.score(P.Population.from_criterion(pd.DataFrame(np.zeros((50, 16)), columns=ODORS), CRIT), ref)
    assert np.isnan(dead.parts['strength/rate']) and np.isnan(dead.error)
    assert dead.parts['breadth/bins'] > 1


def test_wirings_common_error_is_at_most_their_mean(gen, ref):
    wirings = [gen.fly(f'w{i}') for i in range(8)]
    w = P.wirings_as_flies(wirings, ref)
    assert list(w.index[-2:]) == ['mean', 'common']
    for c in [c for c in w.columns if c.startswith('part_')]:
        assert w.loc['common', c] <= w.loc['mean', c] + 1e-12
    gate = P.leave_one_out(wirings)
    assert gate['error'].max() < 1.6


def test_subset_reference_units(ref, flies):
    half = P.subset_reference(ref, [0, 2, 4, 6, 8])
    direct = P.fly_reference([flies[i] for i in (0, 2, 4, 6, 8)])
    for p in ref.parts:
        assert half.one_fly[p] == pytest.approx(direct.one_fly[p])
        pd.testing.assert_series_equal(half.mean[p], direct.mean[p])
    kept = P.subset_reference(ref, ['f0', 'f2', 'f4'], keep_units=True)
    assert kept.one_fly == ref.one_fly and kept.flies == ['f0', 'f2', 'f4']
    halves = P.split_halves(11, 5, seed=2)
    assert len(halves) == 5
    for a, b in halves:
        assert len(a) == 5 and len(b) == 6 and sorted(np.r_[a, b]) == list(range(11))


def test_named_units_and_missing_glomeruli(gen):
    named = [gen.fly(f'f{i}', n_units=30, named=True) for i in range(6)]
    named[1] = P.Population(named[1].values.drop(index=['g0', 'g1']), named[1].responds.drop(index=['g0', 'g1']),
        name='f1', named_units=True)
    ref = P.fly_reference(named)
    assert 'matrix/values' in ref.parts
    assert ref.parts['matrix/values'].shape == (6, 30 * 16)
    assert ref.parts['matrix/values'].loc['f1', 'g0|o0'] != ref.parts['matrix/values'].loc['f1', 'g0|o0']
    s = P.score(gen.fly('m', n_units=30, named=True), ref)
    assert np.isfinite(s.error_with_matrix) and 'matrix' in s.views
    assert s.error == pytest.approx(np.sqrt(np.mean([s.views[v] ** 2 for v in P.VIEWS])))
    assert s.error_with_matrix == pytest.approx(np.sqrt(np.mean([s.views[v] ** 2 for v in (*P.VIEWS, 'matrix')])))


def test_binary_populations_and_part_subsets(gen):
    bins = [P.Population(f.responds, f.responds, name=f.name, graded=False) for f in (gen.fly(f'f{i}') for i in range(5))]
    ref = P.fly_reference(bins)
    assert set(ref.parts) == {'strength/rate', 'sharing/who', 'breadth/bins'}
    rate_only = P.fly_reference([gen.fly(f'f{i}') for i in range(5)], parts=['strength/rate', 'sharing/who',
        'sharing/how_strongly', 'breadth/bins'])
    assert 'strength/mean' not in rate_only.parts
    s = P.score(gen.fly('x'), rate_only)
    assert s.views['strength'] == s.parts['strength/rate']
    assert s.views['sharing'] == pytest.approx(np.sqrt((s.parts['sharing/who'] ** 2 + s.parts['sharing/how_strongly'] ** 2) / 2))


def test_partial_panels(gen, flies):
    f = gen.fly('partial')
    v = f.values.copy()
    v['o3'] = np.nan
    part = P.Population.from_criterion(v, CRIT, name='partial')
    assert 'o3' not in P.common_odors([*flies, part])
    r = P.fly_reference([*flies, part])
    assert 'o3' not in r.odors and len(r.parts['sharing/who'].columns) == 105
    r16 = P.fly_reference([*flies, part], odors=ODORS)
    assert r16.parts['breadth/bins'].loc['partial'].isna().all()
    assert np.isfinite(P.score(gen.fly('y'), r16).error)


def test_steps_and_the_no_step_yardstick(gen, flies, ref):
    same = P.step_interval(flies, flies, n_boot=200)
    for p, df in same.items():
        assert np.allclose(df['change'], 0.0)
        assert ((df['lo'] <= 1e-12) & (df['hi'] >= -1e-12)).all()
    assert P.no_step_score(flies, ref).error == pytest.approx(0.0, abs=1e-12)
    # a downstream stage that sparsens every fly's responses: fewer responders, a change that leaves the interval
    down = [P.Population.from_criterion(f.values, 2.0, name=f.name) for f in flies]
    st = P.step_interval(flies, down, n_boot=200)
    b = st['breadth/bins']
    assert b.loc['0', 'change'] > 0 and b.loc['0', 'lo'] > 0
    dref = P.fly_reference(down)
    assert P.no_step_score(flies, dref).error > 1.0
    with pytest.raises(ValueError):
        P.no_step_score([P.Population.from_criterion(f.values.drop(columns='o0'), CRIT) for f in flies], dref)
    ch = P.step_change(P.mean_views([P.population_views(f) for f in flies]), dref.mean)
    assert np.allclose(ch['breadth/bins'], st['breadth/bins']['change'])


def test_energy_distance():
    rng = np.random.default_rng(5)
    u, v = rng.normal(size=300), rng.normal(0.5, 1.2, size=200)
    assert P.energy_distance(u[:, None], v[:, None]) == pytest.approx(stats.energy_distance(u, v), rel=1e-10)
    x = rng.normal(size=(100, 4))
    assert P.energy_distance(x, x) == pytest.approx(0.0, abs=1e-9)
    assert P.energy_distance(x, x + 1) > 0.5
    big = rng.normal(size=(500, 3))
    # a set's draw depends only on the set and the seed: a set against itself is 0, against another draw of it is not
    assert P.energy_distance(big, big, max_units=200, seed=1) == pytest.approx(0.0, abs=1e-9)
    other = P.energy_distance(big, big[::-1], max_units=200, seed=1)
    assert 0 < other < 0.3 * P.energy_distance(big, big + 1, max_units=200, seed=1)
    x, y = P.energy_sample(big, 200, 1), P.energy_sample(big + 1, 200, 1)
    assert P.energy_from_samples(x, y, P.mean_self_distance(x), P.mean_self_distance(y)) == pytest.approx(
        P.energy_distance(big, big + 1, max_units=200, seed=1))


def test_energy_guard(gen, flies):
    small = [P.Population.from_criterion(f.values.iloc[:150], CRIT, name=f.name) for f in flies[:5]]
    fly = P.energy_guard(small[0], small)
    assert fly['to_flies'] > 0 and 0.5 < fly['ratio'] < 1.6
    vals = small[0].values.to_numpy()
    rand = pd.DataFrame(np.random.default_rng(6).permutation(vals.ravel()).reshape(vals.shape), columns=ODORS)
    far = P.energy_guard(P.Population.from_criterion(rand, CRIT, name='r'), small)
    assert far['ratio'] > fly['ratio']
    ref = np.random.default_rng(7).exponential(size=1000)
    mapped = P.energy_guard(small[0], small, reference=ref)
    assert np.isfinite(mapped['ratio'])


def test_rank_map():
    rng = np.random.default_rng(8)
    v = pd.DataFrame(rng.normal(size=(50, 4)))
    v.iloc[0, 0] = np.nan
    ref = rng.exponential(size=199)
    m = P.rank_map(v, ref)
    assert np.isnan(m.iloc[0, 0])
    ok = v.notna().to_numpy()
    assert stats.spearmanr(v.to_numpy()[ok], m.to_numpy()[ok]).statistic == pytest.approx(1.0)
    assert m.min().min() >= ref.min() and m.max().max() <= ref.max()


def test_fly_populations_adapter(gen):
    f = gen.fly('a')

    def responses(fly):
        raise NotImplementedError

    rec = SimpleNamespace(flies=('a',), tables={'a': f.values}, criterion=CRIT, named_units=False, graded=True,
        responses=responses)
    (p,) = P.fly_populations(rec)
    assert p.name == 'a'
    pd.testing.assert_frame_equal(p.responds, f.responds)
    with pytest.raises(ValueError):
        P.as_population(f.values)
    assert P.as_population(f) is f
    assert P.as_population(f.values, 2.0, at_or_above=True).responds.sum().sum() == (f.values >= 2.0).sum().sum()
