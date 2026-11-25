"""Tests of E137 (the held-out panels read like for like): the interface and the readings' rules on made-up tables."""
import numpy as np
import pandas as pd
import pytest

from experiments import E137_heldout_readings as E


def test_interface():
    assert E.READINGS[0] == 'as_uploaded'
    assert len(E.READINGS) == 5
    assert E.POSITION_READINGS == E.READINGS[1:]
    assert E.V2_RULES == ('common_mean', 'fly_scale')
    assert set(E.PRIMARY_UNITS) == set(E.READINGS)
    for name in ('load_lock', 'fly_mixture_view', 'part_lock', 'part_position', 'part_readings', 'part_rescore',
                 'part_rise', 'part_validation2', 'part_order', 'part_classes'):
        assert callable(getattr(E, name))


def test_session_blocks():
    s = ['pfo @ 0.0'] * 3 + ['a @ -2.0'] * 3 + ['m @ 0.0'] * 3 + ['a @ -2.0'] * 3 + ['pfo @ 0.0'] * 3
    b = E.session_blocks(s)
    assert list(b.block) == [0] * 9 + [1] * 6
    assert list(b.rep) == [0, 1, 2] * 5
    a = b[b.stim == 'a @ -2.0']
    assert list(a.presentation) == [0, 1, 2, 3, 4, 5]


def test_orders_match():
    cfg = ['etoh @ -2.0', 'etoh @ -2.0', 'control mix @ -2.0', 'cmix @ 0.0']
    file = ['EtOH @ -2.0', 'EtOH @ -2.0', 'cmix @ -2.0', 'cmix @ 0.0']
    assert E.orders_match(cfg, file)
    assert not E.orders_match(cfg, file[::-1])
    assert not E.orders_match(cfg[:3], file)
    # two config names onto one file stimulus at the same concentration is not one-to-one
    assert not E.orders_match(['a @ -2.0', 'b @ -2.0'], ['x @ -2.0', 'x @ -2.0'])


def _reduced():
    """Two cells, a session of components c1, c2 (two blocks each) and one mixture, in a made-up layout."""
    stims = (['c1 @ -2.0'] * 3 + ['c2 @ -3.0'] * 3 + ['kmix @ 0.0'] * 3 + ['c1 @ -2.0'] * 3 + ['c2 @ -3.0'] * 3)
    b = E.session_blocks(stims)
    lab = {'c1 @ -2.0': 'c1', 'c2 @ -3.0': 'c2', 'kmix @ 0.0': 'kmix0'}
    rows = []
    for c in (0, 1):
        for t in b.itertuples():
            rows.append(dict(cell=c, trial=t.trial, stim=t.stim, label=lab[t.stim], in_panel=True, block=t.block,
                             rep=t.rep, presentation=t.presentation, Fc=10.0 * (t.block + 1) + t.rep + c,
                             z=float(t.block + 1), baseline_sd=2.0))
    return pd.DataFrame(rows)


def test_select_trials():
    red = _reduced()
    first = E.select_trials(red, 'first_presentations')
    u = first[first.use]
    assert set(zip(u.label, u.block, u.rep)) == {('c1', 0, 0), ('c2', 0, 0), ('kmix0', 0, 0)}
    assert first[first.noise_block].block.max() == 0
    second = E.select_trials(red, 'components_second_block')
    u = second[second.use]
    assert set(u[u.label == 'c1'].block) == {1} and set(u[u.label == 'kmix0'].block) == {0}
    both = E.select_trials(red, 'components_both_blocks')
    u = both[both.use]
    assert set(u[u.label == 'c2'].block) == {0, 1} and len(u[(u.label == 'c2') & (u.cell == 0)]) == 6
    up = E.select_trials(red, 'as_uploaded')
    assert set(up[up.use].block) == {0}
    with pytest.raises(ValueError):
        E.select_trials(red, 'nonsense')


def test_noise_values():
    red = _reduced()
    v = E.noise_values(red, factor=2.0)
    np.testing.assert_allclose(v, red.Fc.to_numpy() / 2.0 * 2.0)


def test_trial_noise_add():
    tn = E.TrialNoise(sd=2.0, n={'a': 4, 'b': 1})
    V = pd.DataFrame(np.zeros((20000, 2)), columns=['a', 'b'])
    out = tn.add(V, np.random.default_rng(0))
    sd = out.std(axis=0)
    assert abs(sd['a'] - 1.0) < 0.03 and abs(sd['b'] - 2.0) < 0.06


def test_stage_class():
    x = pd.DataFrame(dict(view_strength=[0.5] * 6 + [3.0] * 2))
    assert E.stage_class(x, 1.0) == 'inside'
    x = pd.DataFrame(dict(view_strength=[2.5] * 6 + [0.5] * 2))
    assert E.stage_class(x, 1.0) == 'far'
    x = pd.DataFrame(dict(view_strength=[1.5] * 8))
    assert E.stage_class(x, 1.0) == 'between'


def test_order_share_no_position_effect():
    rng = np.random.default_rng(1)
    rows = []
    for rep in range(3):
        order = rng.permutation(4)
        for pos, o in enumerate(order):
            t = rep * 4 + pos
            for c in range(50):
                rows.append(dict(fly='f', cell=c, odor=f'o{o}', trial=t, pos=pos, value=float(o) + 0.0 * pos))
    sh = E.order_share(pd.DataFrame(rows), crit=1.5)
    assert np.allclose(sh.slope, 0.0, atol=1e-9)
    assert np.allclose(sh.rms_shift, 0.0, atol=1e-9)


def test_order_share_with_position_effect():
    rows = []
    for rep in range(3):
        order = [0, 1, 2, 3] if rep % 2 == 0 else [3, 2, 1, 0]
        for pos, o in enumerate(order):
            t = rep * 4 + pos
            for c in range(10):
                rows.append(dict(fly='f', cell=c, odor=f'o{o}', trial=t, pos=pos, value=1.0 + o - 0.2 * pos))
    sh = E.order_share(pd.DataFrame(rows), crit=-99.0).set_index('part')
    assert sh.loc['strength/mean', 'slope'] == pytest.approx(-0.2)
    assert sh.loc['strength/mean', 'rms_shift'] > 0


def test_rise():
    t = pd.DataFrame({'c1': [1.0, 1.0], 'kmix-2': [1.0, 0.0], 'kmix-1': [1.0, 1.0], 'kmix0': [2.0, 2.0]})
    r = (t > 0.5).astype(float)
    out = E._rise(t, r, 'kmix')
    assert out['rate_0_over_m2'] == pytest.approx(2.0)
    assert out['mean_0_over_m2'] == pytest.approx(4.0)
    assert out['rate_mix0_over_components'] == pytest.approx(1.0)


def _scores(kmix_strength):
    """A made-up rescore table: the family's eight wirings, their mean and common rows, per panel and reading in its
    primary units, with the flies' most atypical fly at 1.5."""
    rows = []
    for panel in E.MIXTURE_PANELS:
        for reading in E.READINGS:
            v = kmix_strength[reading] if panel == 'kmix' else 0.5
            for g in E.FAMILY_GROUPS:
                for w in list(E.FAMILY_WIRINGS) + ['mean', 'common']:
                    rows.append(dict(panel=panel, reading=reading, units=E.PRIMARY_UNITS[reading], noise=True,
                                     primary=True, group=g, wiring=w, view_strength=v, view_sharing=1.0,
                                     view_breadth=1.0, error=1.0, loo_max_view_strength=1.5,
                                     floor_view_strength=0.3))
    return pd.DataFrame(rows)


def _lock():
    fl = {'a/1': dict(panel='kmix', found=True, reproduction=dict(passed=True)),
          'b/1': dict(panel='cmix', found=True, reproduction=dict(passed=True))}
    return dict(flies=fl)


def test_classify_and_outcome():
    far = {r: 2.5 for r in E.READINGS}
    c = E.classify(_scores(far), _lock())
    assert c['kmix/fitted']['primary']['readable'] and c['cmix/fitted']['primary']['readable']
    assert E.outcome(c)['outcomes'] == ['b']
    mixed = dict(far, components_second_block=0.7)
    c = E.classify(_scores(mixed), _lock())
    assert not c['kmix/fitted']['primary']['readable']
    assert c['kmix/fitted']['primary']['bounds'] == [0.7, 2.5]
    assert E.outcome(c)['outcomes'] == ['c']
    near = {r: 1.0 for r in E.READINGS}
    c = E.classify(_scores(near), _lock())
    assert E.outcome(c, {'fitted': {'readable': False}})['outcomes'] == ['a', 'd']
    bad = _lock()
    bad['flies']['a/1']['reproduction']['passed'] = False
    c = E.classify(_scores(near), bad)
    assert 'e' in E.outcome(c)['outcomes'] and not c['kmix/fitted']['primary']['replaces']['first_presentations']
