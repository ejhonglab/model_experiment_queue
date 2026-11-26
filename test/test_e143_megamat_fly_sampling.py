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


# ---------------------------------------------------------------------------
# the classes on planted values

def test_tells_apart():
    assert E.tells_apart([0.9, 0.85, 0.5, 0.95])
    assert not E.tells_apart([0.9, 0.85, 0.49, 0.95])      # one split under 0.5
    assert not E.tells_apart([0.79, 0.7, 0.9])              # median under 0.8
    assert not E.tells_apart([np.nan])


def test_change_stands_and_class_stands():
    assert E.change_stands(-0.1, -0.02, 0.8)
    assert not E.change_stands(-0.1, 0.01, 1.0)             # the interval holds zero
    assert not E.change_stands(-0.1, -0.02, 0.79)           # beyond the bar on too few joint splits
    assert E.change_stands(0.01, 0.2)                       # a held-out panel: the interval alone
    assert not E.change_stands(-0.01, 0.2)
    assert E.class_stands(0.8) and not E.class_stands(0.799) and not E.class_stands(np.nan)


def test_beyond_bar():
    spread = 0.02                                           # se at least 0.01, the bar 0.02
    ch = np.array([[-0.03, -0.03, -0.03, -0.03],            # beyond, in the scored direction
                   [-0.015, -0.015, -0.015, -0.015],        # inside the bar
                   [0.03, 0.03, 0.03, 0.03],                # the wrong direction
                   [-0.1, 0.0, -0.1, 0.0]])                 # mean -0.05, SD over sqrt(4) 0.029: bar 0.058
    assert E.beyond_bar(ch, -1.0, spread).tolist() == [True, False, False, False]
    assert E.beyond_bar(-ch, 1.0, spread).tolist() == [True, False, False, False]


def test_count_class():
    def rows(key, meds):
        return [dict(key=key, half=k, median=m) for k, m in meds.items()]
    df = pd.DataFrame(rows('p1', {2: 0.3, 3: 0.6, 4: 0.85}) + rows('p2', {2: 0.4, 3: 0.9}) + rows('p3', {2: 0.9, 3: 0.95}))
    c = E._count_class(df, ('p1', 'p2', 'p3'))
    assert c['count'] and not c['panel']
    assert c['per_panel']['p1']['first_half_at_high'] == 4 and c['per_panel']['p3']['panel']
    df = pd.DataFrame(rows('p1', {2: 0.85, 3: 0.9}) + rows('p2', {2: 0.82}) + rows('p3', {2: 0.3, 3: 0.5}))
    c = E._count_class(df, ('p1', 'p2', 'p3'))
    assert c['panel'] and not c['count']
    assert not c['per_panel']['p2']['count']                # halves of 2 only: never the count's


def test_outcome_letters():
    tells = dict(pn_dendrite=True, pn_bouton=True, kc_claw=False, kc_soma=False)
    kw = dict(reads_calyx=False, a1_fall=[True], a1_holds=[True], near_shares=[0.9], count=False, panel=False)
    assert E.outcome_letters(tells, **kw) == ['a']
    assert E.outcome_letters(tells, **dict(kw, reads_calyx=True)) == []
    assert E.outcome_letters({s: True for s in E.LOSS_STAGES}, **kw) == ['b']
    assert E.outcome_letters(tells, **dict(kw, a1_holds=[True, False])) == ['a', 'c']
    assert E.outcome_letters(tells, **dict(kw, a1_fall=[False])) == ['a', 'c']
    assert E.outcome_letters(tells, **dict(kw, near_shares=[0.9, 0.79])) == ['a', 'd']
    assert E.outcome_letters(tells, **dict(kw, panel=True)) == ['a', 'e']
    assert E.outcome_letters(dict(tells, kc_soma=True), **kw) == []


def _rv(arm_stage, fam_stage=1.0, stage_ch=None, panel_ch=None, guard_ch=None, K=1):
    """A planted Resampled: every value of the family 1.0 on the four T wirings; the arm's megamat stage error
    `arm_stage`, its stages' and panels' changes as given (0 elsewhere)."""
    w = len(E.T_WIRINGS)
    full = lambda v: np.full((K, w), float(v))
    stage_ch, panel_ch, guard_ch = stage_ch or {}, panel_ch or {}, guard_ch or {}
    mega = {'stage': {'family': full(fam_stage), 'arm': full(arm_stage)}}
    for st in E.LOSS_STAGES:
        mega[st] = {'family': full(1.0), 'arm': full(1.0 + stage_ch.get(st, 0.0))}
    guard = {st: {'family': full(1.0), 'arm': full(1.0 + guard_ch.get(st, 0.0))} for st in E.LOSS_STAGES}
    panels = {k: {'family': full(1.0), 'arm': full(1.0 + panel_ch.get(k, 0.0))} for k in E.HOLDS_PANELS + E.BESIDE_PANELS}
    gpanels = {k: {'family': full(1.0), 'arm': full(1.0 + guard_ch.get(k, 0.0))} for k in E.HOLDS_PANELS + E.BESIDE_PANELS}
    return E.Resampled(mega, guard, panels, gpanels)


def _rec(reasons=(), control_fall=None, guard=True):
    sc = lambda spread: dict(complete=True, se=0.01, spread=spread)
    rec = dict(megamat=sc(0.02), stages={st: sc(0.02) for st in E.LOSS_STAGES},
               panels={k: sc(0.03) for k in E.HOLDS_PANELS}, panels_beside={k: sc(0.1) for k in E.BESIDE_PANELS},
               reasons=list(reasons), control_fall=control_fall, cls='holds')
    if guard:
        rec['guard'] = {k: sc(0.02) for k in list(E.LOSS_STAGES) + list(E.HOLDS_PANELS + E.BESIDE_PANELS)}
    return rec


def test_reclass_ladder():
    cls = lambda rv, rec=None: E.reclass(rec or _rec(), 'arm', rv, E.E136)[0][0]
    assert cls(_rv(0.9)) == 'holds'                                                  # -0.1 < -2 x 0.01
    assert cls(_rv(0.99)) == 'does not improve'                                      # -0.01: inside the bar
    assert cls(_rv(0.9, panel_ch={'apl16/apl16/on': 0.05})) == 'megamat only'        # a panel worse than 0.03
    assert cls(_rv(0.9, panel_ch={'mixtures/kmix/on': 0.5})) == 'holds'              # the beside panel gates nothing
    assert cls(_rv(0.9, stage_ch={'kc_soma': 0.05})) == "one stage at another's cost"
    assert cls(_rv(0.9, guard_ch={'kc_claw': 0.05})) == "strength at another view's cost"
    assert cls(_rv(0.9, guard_ch={'validation2/validation2/on': 0.05})) == 'megamat only'  # the guard's panel
    assert cls(_rv(0.9), _rec(reasons=['at a bound: kckc'])) == 'no class'           # a fixed reason
    assert cls(_rv(0.9), _rec(reasons=['no better than the permuted control'], control_fall=0.05)) == 'holds'
    assert cls(_rv(0.97), _rec(reasons=['no better than the permuted control'], control_fall=0.04)) == 'no class'
    assert E.reclass(dict(cls='no class', reasons=['reads stopped']), 'arm', _rv(0.9), E.E141)[0][0] == 'no class'


def test_reclass_on_many_references_and_share():
    K = 4
    rv = _rv(0.9, K=K)
    rv.mega['stage']['arm'][2:] = 0.995                     # two of four resamples lose the fall
    cls, comps = E.reclass(_rec(), 'arm', rv, E.E136)
    assert cls.tolist() == ['holds', 'holds', 'does not improve', 'does not improve']
    assert np.allclose(comps['megamat']['change'], [-0.1, -0.1, -0.005, -0.005])
    assert E.class_stands(float(np.mean(cls == 'holds'))) is False


def test_e141_entry_does_more_than_a1():
    cls = np.array(['holds', 'holds', 'megamat only'], dtype=object)
    vs = {'vs_a1/megamat': dict(improves=np.array([True, False, True])),
          'vs_a1/panel/x': dict(better=np.array([False, True, True]), worse=np.array([False, False, False]))}
    out = E.e141_entry(dict(does_more_than_a1=False, stands_in=None), cls, vs)
    assert out['does_more_than_a1'].tolist() == [True, True, False] and 'stands_in' not in out
    out = E.e141_entry(dict(does_more_than_a1=None, stands_in=False), cls, vs)
    assert out['stands_in'].tolist() == [True, True, False]


def test_disjoint_halves():
    hv = E.disjoint_halves(11, 5, n_draws=50, seed=0)
    assert len(hv) == 50 and all(len(a) == len(b) == 5 and not set(a) & set(b) for a, b in hv)
    assert hv == E.disjoint_halves(11, 5, n_draws=50, seed=0)
    with pytest.raises(ValueError):
        E.disjoint_halves(4, 3)


def test_setting_of_and_settings():
    assert E.setting_of(E.E136, 'A1') == 'E136:A1' and E.setting_of(E.E139, 'A1') == 'E136:A1'
    assert E.setting_of(E.E141, 'K_inhadd1000') == 'E139:K_inhadd1000'
    assert E.setting_of(E.E141, 'D1.8_inhadd1000') == 'E141:D1.8_inhadd1000'
    assert E.setting_of(E.E139, 'A0_inhadd314') == 'E139:A0_inhadd314'
    assert len(E.CALYX_GROUPS['inhadd1000']) == 9 and len(E.CALYX_GROUPS['family_norm']) == 5
    assert len(E.setting_names()) == 23 and len(E.setting_names('yardstick')) == 5
    assert all(E.setting_of(t, a) in E.SETTING_BY_NAME for t, a in E.NEAR_TIES)
