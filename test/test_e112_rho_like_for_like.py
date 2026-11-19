"""E112 (the two ρ's on odor-odor distances read as the recordings are made), on made-up counts and tables: the
pseudo-recordings of a session z-score reading (no noise gives the cosine of the reading, the draws are seeded and
converge on the reading in expectation), the pipeline table read as R4 (the matched levels, a shared level counted for
both), the scores and E106's bound at the data itself, the verdict "the reading moves ρ", the substitution of re-read
columns into a saved readings table and the identity check of rebuilt rows."""
import numpy as np
import pandas as pd
import pytest

from kcq import harness as H
from experiments import E112_rho_like_for_like as E

ODORS = list(H.APL16_ODORS)


def pairs_series(values) -> pd.Series:
    idx = pd.MultiIndex.from_tuples(H.pair_list(ODORS), names=['odor_a', 'odor_b'])
    return pd.Series(np.asarray(values, dtype=float), index=idx)


def made_up_counts(seed=0, n=300) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    x = rng.poisson(rng.gamma(0.5, 2.0, (n, 16))) * (rng.random((n, 16)) < 0.2)
    return pd.DataFrame(x.astype(float), columns=ODORS)


def z_of(center: np.ndarray, trial_sd: np.ndarray) -> dict:
    return dict(center=center[None], trial_sd=trial_sd[None], n_trials=3, odors=None)


def test_interface():
    assert E.READINGS == ('R0', 'R1', 'R2', 'R3', 'R4')
    assert E.N_PSEUDO == {'on': 11, 'off': 8}
    for f in ('reading_distances', 'pseudo_recording_distances', 'trial_noise_distances', 'reading_scores', 'bound_scores',
              'moves_verdict', 'reread_gate', 'e109_holds', 'part_a', 'part_b_prepare', 'part_b_score'):
        assert callable(getattr(E, f))


def test_pseudo_recordings_without_noise_are_the_reading():
    c = made_up_counts(1).to_numpy() / 3.0
    z = z_of(c, np.zeros(c.shape[0]))
    d = E.pseudo_recording_distances(z, 4, seed=3)
    want = E._cs(H.pair_distances(pd.DataFrame(c, columns=ODORS), 'cosine'))
    assert np.allclose(d.reindex(want.index).to_numpy(), want.to_numpy(), atol=1e-12)


def test_pseudo_recordings_seeded_and_converge_on_the_expectation():
    c = made_up_counts(2, 400).to_numpy() / 2.0
    sd = np.full(c.shape[0], 0.4)
    z = z_of(c, sd)
    a, b = E.pseudo_recording_distances(z, 3, seed=7), E.pseudo_recording_distances(z, 3, seed=7)
    assert a.equals(b)
    assert not a.equals(E.pseudo_recording_distances(z, 3, seed=8))
    # in expectation (many pseudo-recordings) the sampled distances approach the closed form of the reading
    many = E.pseudo_recording_distances(z, 200, seed=1)
    want = E._cs(H.zcell_pair_distances(dict(z, odors=ODORS), 'graded'))
    assert np.abs(many.reindex(want.index).to_numpy() - want.to_numpy()).max() < 0.01
    # noise raises the distances on average (the pairs are farther apart through the noise)
    noiseless = E._cs(H.pair_distances(pd.DataFrame(c, columns=ODORS), 'cosine'))
    assert many.mean() > noiseless.mean()


def test_r3_conditions_draw_their_own_recordings():
    c = made_up_counts(3).to_numpy() / 3.0
    z = z_of(c, np.full(c.shape[0], 0.3))
    r = E.r3_from_z(z, z, seed=0)
    assert set(r) == {'on', 'off'}
    assert not r['on'].equals(r['off'])           # 11 against 8 pseudo-recordings, other draws
    assert E.r3_from_z(z, z, seed=0)['on'].equals(r['on'])


def pipeline_table(levels) -> pd.DataFrame:
    rows = []
    pl = H.pair_list(ODORS)
    rng = np.random.default_rng(0)
    for model, lv in levels.items():
        for level in lv:
            for cond in ('on', 'acute_off'):
                for tr in ('graded', 'thresholded'):
                    for a, b in pl:
                        rows.append(dict(model=model, level=level, condition=cond, transform=tr, odor_a=a, odor_b=b,
                                         raw=rng.random(), det=0.5,
                                         tag=f'{level}|{cond}|{tr}'))
    return pd.DataFrame(rows)


def test_trial_noise_distances_take_the_matched_levels():
    t = pipeline_table({'m1': ('match_on', 'match_off', 'lo'), 'm2': ('match_on+match_off',)})
    r = E.trial_noise_distances(t)
    assert set(r) == {'m1', 'm2'}
    for m, lv_on, lv_off in (('m1', 'match_on', 'match_off'), ('m2', 'match_on+match_off', 'match_on+match_off')):
        for c, cond, lv in (('on', 'on', lv_on), ('off', 'acute_off', lv_off)):
            x = t[(t.model == m) & (t.level == lv) & (t.condition == cond) & (t.transform == 'graded')]
            want = E._cs(pd.Series(x['raw'].to_numpy(), index=pd.MultiIndex.from_arrays([x.odor_a, x.odor_b])))
            assert np.allclose(r[m][c].reindex(want.index).to_numpy(), want.to_numpy())
    assert (E.trial_noise_distances(t, 'det')['m1']['on'] == 0.5).all()


def made_up_data(seed=0):
    rng = np.random.default_rng(seed)
    on = pairs_series(rng.uniform(0.3, 0.9, 120))
    off = pairs_series(on.to_numpy() - rng.uniform(0.05, 0.3, 120))
    return {'on': E._cs(on), 'off': E._cs(off)}


def test_scores_and_bound_at_the_data():
    D = made_up_data(1)
    eff = D['on'] - D['off']
    s = E.reading_scores(D, data=D, effect=eff)
    for c in ('on', 'off'):
        for sub in E.SUBSETS:
            assert s[f'rho_{c}_{sub}'] == pytest.approx(1.0)
        assert s[f'paa_bmyr_{c}'] == pytest.approx(E._pair_value(D[c], 'paa', 'B-myr'))
    assert s['loss'] == pytest.approx(0.0, abs=1e-12) and s['pattern'] == pytest.approx(0.0, abs=1e-12)
    assert s['rho_effect'] == pytest.approx(1.0)
    b = E.bound_scores(D, data=D)
    assert all(v == pytest.approx(1.0) for k, v in b.items() if k.startswith(('at_', 'model_', 'flies_'))
               and not k.startswith('model_pair_optimal_odor'))
    # a model off the data: the bound with paa's and B-myr's pairs at the data's values is at least as high on those
    rng = np.random.default_rng(5)
    M = {c: E._cs(D[c] + pd.Series(rng.normal(0, 0.15, 120), index=D[c].index)) for c in ('on', 'off')}
    b = E.bound_scores(M, data=D)
    s = E.reading_scores(M, data=D, effect=eff)
    assert b['at_paa_bmyr_2but_raw_on'] >= s['rho_on_all120'] - 0.05
    assert b['model_on'] == pytest.approx(s['rho_on_all120'])


def rho_rows(base, new, wirings=('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')):
    rows = []
    for r, vals in (('R0', base), ('R3', new)):
        for u in list(wirings) + ['all8']:
            on, off = vals[u]
            rows.append(dict(reading=r, unit=u, kind='population' if u == 'all8' else 'wiring', rho_on_all120=on,
                             rho_off_all120=off))
    return pd.DataFrame(rows)


def test_moves_verdict():
    ws = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')
    base = {u: (0.66, 0.74) for u in ws + ('all8',)}
    up = {u: (0.69, 0.742) for u in ws + ('all8',)}
    v = E.moves_verdict(rho_rows(base, up), 'R3')
    assert v['moves'] and v['moves_on'] and not v['moves_off'] and v['n_wirings_move_on'] == 8 and not v['provisional']
    # only five wirings up: not 6 of 8
    some = dict(up, **{w: (0.665, 0.74) for w in ('P3', 'P4', 'M2')})
    v = E.moves_verdict(rho_rows(base, some), 'R3')
    assert v['n_wirings_move_on'] == 5 and not v['moves']
    # the other condition lower by more than 0.005
    worse = {u: (0.69, 0.73) for u in ws + ('all8',)}
    assert not E.moves_verdict(rho_rows(base, worse), 'R3')['moves']
    # the sign split
    split = dict(up, **{w: (0.64, 0.74) for w in ('P3', 'P4')})
    assert E.moves_verdict(rho_rows(base, split), 'R3')['sign_split']
    assert E.wirings_needed(8) == 6 and E.wirings_needed(4) == 3


def test_substitute_and_compare_rows():
    ex = 'E109_input_search_rho'
    saved = pd.DataFrame(dict(change=['I0', 'I0', 'A1', 'A1'], unit=['T', 'T', 'T', 'T'], readout=['raw', 'cal', 'raw', 'cal'],
                              rho_on_all=[0.6, 0.61, 0.65, 0.66], all_pattern=[0.05, 0.04, 0.05, 0.04], n_flagged=[0, 0, 0, 0],
                              run_rate_on=[0.07] * 4))
    new = pd.DataFrame(dict(change=['I0', 'A1'], unit=['T', 'T'], readout=['R3', 'R3'], rho_on_all=[0.7, 0.72],
                            all_pattern=[0.03, 0.035]))
    r = E.substitute(saved, new, ex, 'rho')
    assert list(r['rho_on_all']) == [0.7, 0.61, 0.72, 0.66]
    assert list(r['all_pattern']) == [0.05, 0.04, 0.05, 0.04]
    a = E.substitute(saved, new, ex, 'all')
    assert list(a['rho_on_all']) == [0.7, 0.7, 0.72, 0.72] and list(a['all_pattern']) == [0.03, 0.03, 0.035, 0.035]
    assert list(a['run_rate_on']) == [0.07] * 4
    same = saved[saved.readout == 'raw'].drop(columns=['n_flagged', 'run_rate_on'])
    c = E.compare_rows(saved, same, ex, 'raw')
    assert c['max_abs_diff'] == 0.0 and c['n_rows'] == 2
    off = same.assign(rho_on_all=same['rho_on_all'] + 0.01)
    assert E.compare_rows(saved, off, ex, 'raw')['max_abs_diff'] == pytest.approx(0.01)
