#!/usr/bin/env python3
"""Tests for E86 (the reduced searched input as a gated candidate; E86_reduced_input_candidate.py) and for the arguments
E86 added to E77's `input_oracle` (the split's seed, a tag for its files, a checkpoint per generation with resume).

Without simulations: B1's permuted inputs are E83's `permuted_factors` draws at seeds 1 and 2 and keep each odor's
factors as a set; the inverted input's and the candidate's factors multiply to one; halves C and D are disjoint, cover
all recordings and keep the two strong recordings apart; the searches' generator seeds; `input_oracle`'s defaults name
E77's files; a search interrupted and resumed from its checkpoint continues exactly (on a stand-in simulator, so it runs
in seconds); `change_stats`' group shares; the verdict on made-up tables. The first-generation check of part C
(`input_oracle` at its defaults reproduces E77's oracle_history_A.csv to the last digit) simulates a generation of 32
candidates (minutes): it runs here only with E86_SIMULATE=1, otherwise as `--part C-check`. Needs the shared directory;
skipped without it."""

import os

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
E86 = pytest.importorskip('experiments.E86_reduced_input_candidate')
E83 = pytest.importorskip('experiments.E83_cross_wiring_fit')
E77 = pytest.importorskip('experiments.E77_input_side_common_error')


def _independent_permutation(F: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Each odor's factors permuted over its responding glomeruli, one generator, odors in the lab input's order (E83's
    part D control, written out again)."""
    orn = H.abbrev_columns(H.apl16_orn_deltas())
    rng = np.random.default_rng(seed)
    out = F.copy()
    for o in orn.columns:
        resp = np.flatnonzero(orn[o].to_numpy(dtype=float) != 0)
        v = F[o].to_numpy()[resp]
        col = F[o].to_numpy().copy()
        col[resp] = v[rng.permutation(len(v))]
        out[o] = col
    return out


@pytest.mark.parametrize('half', ['A', 'B'])
@pytest.mark.parametrize('seed', E86.B1_SEEDS)
def test_b1_inputs_are_e83_permuted_factors(half, seed):
    assert E86.B1_SEEDS == (1, 2) and E86.B3_SEEDS == (3, 4)
    F = E86.candidate_factors(half)
    P, size = E86.control_factors('perm', half, seed)
    assert size == E86.SIZE == E83.GATE_SIZE == 0.6
    pd.testing.assert_frame_equal(P, E83.permuted_factors(F, seed))
    pd.testing.assert_frame_equal(P, _independent_permutation(F, seed))
    # the planned point's input is that draw at 0.6
    pt = [p for p in E86.plan_b('B1') if p['half'] == half and p['seed'] == seed][0]
    pd.testing.assert_frame_equal(E86.point_input(None, pt), E83.input_at(P, 0.6))
    # the two seeds differ, and differ from the candidate
    other = E86.control_factors('perm', half, 3 - seed)[0]
    assert not P.equals(other) and not P.equals(F)


@pytest.mark.parametrize('half', ['A', 'B'])
def test_permuted_keeps_each_odors_factors_as_a_set(half):
    F = E86.candidate_factors(half)
    orn = H.abbrev_columns(H.apl16_orn_deltas())
    searched = [o for o in F.columns if (F[o] != 1).any()]
    assert sorted(searched) == sorted(E77.E_ODORS)
    for seed in E86.B1_SEEDS + E86.B3_SEEDS:
        P = E86.control_factors('perm', half, seed)[0]
        for o in F.columns:
            assert np.array_equal(np.sort(P[o].to_numpy()), np.sort(F[o].to_numpy())), (seed, o)
            # glomeruli without a response keep their factor
            silent = orn[o].to_numpy(dtype=float) == 0
            assert np.array_equal(P[o].to_numpy()[silent], F[o].to_numpy()[silent])
        moved = sum(int((P[o].to_numpy() != F[o].to_numpy()).sum()) for o in searched)
        assert moved > 0.6 * 36 * len(searched)


@pytest.mark.parametrize('half', ['A', 'B'])
def test_inverted_and_candidate_multiply_to_one(half):
    F, size = E86.control_factors('inverted', half)
    assert size == E86.INVERT_SIZE == -0.6
    pd.testing.assert_frame_equal(F, E86.candidate_factors(half))
    orn = H.apl16_orn_deltas().astype(float)
    cand = E83.input_at(F, E86.SIZE)
    inv = E83.input_at(F, E86.INVERT_SIZE)
    nz = orn.to_numpy() != 0
    prod = (cand.to_numpy()[nz] / orn.to_numpy()[nz]) * (inv.to_numpy()[nz] / orn.to_numpy()[nz])
    assert np.allclose(prod, 1.0, rtol=1e-12, atol=0)
    assert np.allclose(np.power(F.to_numpy(), 0.6) * np.power(F.to_numpy(), -0.6), 1.0, rtol=1e-12, atol=0)
    # and the inverted input moves every searched entry the other way
    ch = np.log(cand.to_numpy()[nz] / orn.to_numpy()[nz])
    iv = np.log(inv.to_numpy()[nz] / orn.to_numpy()[nz])
    assert np.allclose(ch, -iv, atol=1e-12)


def test_halves_c_and_d():
    hv = E86.halves()
    keys = lambda h: {(g, tuple(map(str, r))) for g in ('ctrl', 'tnt') for r in hv[h][g]}
    ctrl = H.yang_recording_rdms_values('TNT3c', 'graded')
    tnt = H.yang_recording_rdms_values('TNT', 'graded')
    every = {('ctrl', tuple(map(str, r))) for r in ctrl.index} | {('tnt', tuple(map(str, r))) for r in tnt.index}
    for a, b in (('A', 'B'), ('C', 'D')):
        assert not keys(a) & keys(b)
        assert keys(a) | keys(b) == every
        strong = [{r for r in hv[h]['tnt'] if (str(r[0]), int(r[1]), str(r[2])) in E77.STRONG_TNT} for h in (a, b)]
        assert len(strong[0]) == 1 and len(strong[1]) == 1
        # a fly's hemispheres stay together
        for g in ('ctrl', 'tnt'):
            fa = {(str(r[0]), int(r[1])) for r in hv[a][g]}
            fb = {(str(r[0]), int(r[1])) for r in hv[b][g]}
            assert not fa & fb
    r1 = E77.recording_halves(1)
    assert hv['C'] == r1['A'] and hv['D'] == r1['B']
    r0 = E77.recording_halves()
    assert hv['A'] == r0['A'] and hv['B'] == r0['B']
    # the entry's description of the second split: D holds 7 of A's 8 recordings, C 8 of B's 11
    assert len(keys('D') & keys('A')) == 7 and len(keys('A')) == 8
    assert len(keys('C') & keys('B')) == 8 and len(keys('B')) == 11
    t = E86.targets()
    for h in 'ABCD':
        assert len(t[f'half{h}']) == 120 and not t[f'half{h}'].isna().any()
    assert {E86.UNSEEN[h] for h in 'CD'} == {'C', 'D'} and E86.split_of('C') == E86.split_of('D') == 1


def test_search_seeds():
    # E77's own searches: generator seeds 0 (half A) and 1 (half B) at input_oracle's default seed
    assert E77.oracle_generator_seed(0, 'A') == 0 and E77.oracle_generator_seed(0, 'B') == 1
    calls = {n: E86.search_call(n) for n in E86.SEARCHES}
    assert {n: E77.oracle_generator_seed(c['seed'], c['half_name']) for n, c in calls.items()} == \
        {'C1_A': 10, 'C1_B': 11, 'C2_C': 20, 'C2_D': 21}
    assert [(c['half_name'], c['split_seed']) for c in calls.values()] == [('A', 0), ('B', 0), ('A', 1), ('B', 1)]
    assert all(c['tag'] == n for n, c in calls.items())


def test_oracle_paths(tmp_path):
    p = E77.oracle_paths(tmp_path, 'A')
    assert p == dict(history=tmp_path / 'oracle_history_A.csv', factors=tmp_path / 'oracle_factors_A.csv',
                     input=tmp_path / 'inputs' / 'oracle_halfA.csv')
    q = E77.oracle_paths(tmp_path, 'B', 'C2_D')
    assert q == dict(history=tmp_path / 'oracle_history_C2_D.csv', factors=tmp_path / 'oracle_factors_C2_D.csv',
                     input=tmp_path / 'inputs' / 'oracle_C2_D.csv')


# ---------------------------------------------------------------------------
# input_oracle's checkpoint and resume, on a stand-in for C2's simulator

class _Killed(Exception):
    pass


def _fake_parts():
    """A stand-in for C2's pieces: the lab's input, a deterministic column-wise model of 80 KCs (counts from a fixed
    random wiring; APL on lowers them), its saved runs at the lab's input."""
    orn = H.apl16_orn_deltas()
    rng = np.random.default_rng(5)
    W = rng.gamma(0.6, 1.0, (80, len(orn.index))) * (rng.random((80, len(orn.index))) < 0.2)
    index = pd.MultiIndex.from_arrays([np.arange(80) + 1000, ['ab'] * 80], names=['kc_id', 'kc_type'])

    def sim_for(o: pd.DataFrame):
        cols = pd.Index(list(o.columns.get_level_values('odor')), name='odor')

        def sim(p):
            drive = W @ o.to_numpy(dtype=float) / 40.0
            inh = 2.0 if p['wAPLKC_scale'] > 0 else 0.0
            counts = np.floor(np.clip(drive - p['fixed_thr'] - inh, 0, None))
            s = pd.DataFrame(counts, index=index, columns=cols)
            return (s >= 2).astype(float), s
        return sim

    p = dict(fixed_thr=1.0, wAPLKC_scale=1.0, wAPLPN_scale=1.0)
    lab = sim_for(orn)
    saved = {'on': lab(p), 'off': lab(dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0))}
    return dict(orn=orn, p=p, breadth=None, saved=saved), sim_for


@pytest.fixture
def fake_c2(monkeypatch):
    from experiments import E77_input_screen as S
    parts, sim_for = _fake_parts()
    calls = dict(n=0, kill_at=None)

    def build(extra_kws=None, *, orn=None, wPNKC=None, breadth=None, params=None):
        calls['n'] += 1
        if calls['kill_at'] is not None and calls['n'] >= calls['kill_at']:
            raise _Killed(calls['n'])
        return sim_for(orn), dict(orn=orn)

    monkeypatch.setattr(S, 'c2_parts', lambda: parts)
    monkeypatch.setattr(S, 'build_c2_simulator', build)
    monkeypatch.setattr(E77, 'c2_readout_K', lambda: 4.0)
    return calls


def _read_outputs(out, half, tag):
    p = E77.oracle_paths(out, half, tag)
    h = pd.read_csv(p['history'], float_precision='round_trip').drop(columns='seconds')
    f = pd.read_csv(p['factors'], float_precision='round_trip')
    i = pd.read_csv(p['input'], float_precision='round_trip')
    return h, f, i


ORACLE_KW = dict(odors=('eb', 'ms'), pop=4, gens=3, seed=20, split_seed=1)


def test_input_oracle_defaults_name_e77s_files(tmp_path, fake_c2):
    E77.input_oracle(tmp_path, 'A', odors=('eb',), pop=3, gens=1)
    assert sorted(p.name for p in tmp_path.iterdir()) == ['inputs', 'oracle_factors_A.csv', 'oracle_history_A.csv']
    assert [p.name for p in (tmp_path / 'inputs').iterdir()] == ['oracle_halfA.csv']


@pytest.mark.parametrize('kill_at', [2, 4, 5, 6])
def test_input_oracle_resumes_exactly(tmp_path, fake_c2, kill_at):
    # builds: 1 the factor-1 check, then one per generation (eb: 2-4, ms: 5-7); a kill at build k leaves the
    # checkpoint of the generation before it (none at k = 2)
    ref = tmp_path / 'ref'
    ref.mkdir()
    E77.input_oracle(ref, 'B', tag='t', **ORACLE_KW)
    ref_out = _read_outputs(ref, 'B', 't')
    assert len(ref_out[0]) == 6
    # with a checkpoint and no interruption: the same
    ck_run = tmp_path / 'ck'
    ck_run.mkdir()
    fake_c2['n'] = 0
    E77.input_oracle(ck_run, 'B', tag='t', checkpoint=tmp_path / 'ck.pkl', **ORACLE_KW)
    for a, b in zip(_read_outputs(ck_run, 'B', 't'), ref_out):
        pd.testing.assert_frame_equal(a, b)
    # interrupted, then resumed
    run = tmp_path / 'run'
    run.mkdir()
    ck = tmp_path / 'run.pkl'
    fake_c2.update(n=0, kill_at=kill_at)
    with pytest.raises(_Killed):
        E77.input_oracle(run, 'B', tag='t', checkpoint=ck, **ORACLE_KW)
    assert ck.exists() == (kill_at > 2)
    fake_c2.update(n=0, kill_at=None)
    E77.input_oracle(run, 'B', tag='t', checkpoint=ck, **ORACLE_KW)
    for a, b in zip(_read_outputs(run, 'B', 't'), ref_out):
        pd.testing.assert_frame_equal(a, b)
    # a checkpoint of another search is refused
    with pytest.raises(ValueError, match='another search'):
        E77.input_oracle(run, 'B', tag='t', checkpoint=ck, **dict(ORACLE_KW, seed=21))


@pytest.mark.skipif(not os.environ.get('E86_SIMULATE'), reason='simulates a generation of 32 candidates (minutes); '
                    'set E86_SIMULATE=1, or run --part C-check')
def test_first_generation_reproduces_e77(tmp_path, monkeypatch):
    monkeypatch.setattr(E86, '_search_dir', lambda: tmp_path)
    res = E86.check_search(tmp_path)
    assert res['reproduced']


# ---------------------------------------------------------------------------
# part A's statistics

def test_change_stats():
    rng = np.random.default_rng(0)
    idx = E86.targets()['z8'].index
    m0 = pd.Series(rng.normal(0.1, 0.05, len(idx)), index=idx)
    m1 = m0 + rng.normal(0.0, 0.02, len(idx))
    t = pd.Series(rng.normal(0.18, 0.06, len(idx)), index=idx)
    s = E86.change_stats(m0, m1, t)
    assert s['n_ester'] == 15 and s['n_alcohol'] == 6 and s['n_other'] == 99
    assert s['sq_change_ester'] + s['sq_change_alcohol'] + s['sq_change_other'] == pytest.approx(s['sq_change_total'])
    assert 1.0 + s['sq_change_total'] == pytest.approx((s['changed'] / s['unchanged']) ** 2)
    # a change that delivers exactly what the target asks takes the pattern part to 0
    e0 = (m0 - m0.mean()) - (t - t.mean())
    s2 = E86.change_stats(m0, m0 - e0, t)
    assert s2['changed'] == pytest.approx(0, abs=1e-12) and s2['share_asked'] == pytest.approx(1.0)
    assert s2['corr_asked'] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# the verdict on made-up tables

FLIES = 0.178


def _row(change, group, searched_on, g, singles=3, raw=-0.004, mean_effect=0.181, sq_all=0.2, sq_six=0.15, flagged=0):
    return dict(change=change, group=group, searched_on=searched_on, unseen=E86.UNSEEN[searched_on], g=g,
                n_singles_lower=singles, n_singles=4, raw_loss_change=raw, mean_effect_cal=mean_effect,
                sq_reduction_all=sq_all, sq_reduction_six=sq_six, n_flagged=flagged, n_points=4)


def _table(cand=(0.0036, 0.0031), b1=(-0.003, -0.004, -0.002, -0.005), b2=(-0.012, -0.015), b3=None,
           c1=(0.0025, 0.002), c2=(0.003, 0.0032), **kw):
    rows = [_row('cand_A', 'candidate', 'A', cand[0], **kw), _row('cand_B', 'candidate', 'B', cand[1], **kw)]
    rows += [_row(f'B1_perm{s}_{X}', 'B1', X, g) for (X, s), g in zip([('A', 1), ('A', 2), ('B', 1), ('B', 2)], b1)]
    rows += [_row(f'B2_inv_{X}', 'B2', X, g) for X, g in zip('AB', b2)]
    if b3 is not None:
        rows += [_row(f'B3_perm{s}_{X}', 'B3', X, g) for (X, s), g in zip([('A', 3), ('A', 4), ('B', 3), ('B', 4)], b3)]
    rows += [_row(f'C1_{X}', 'C1', X, g, **kw) for X, g in zip('AB', c1)]
    rows += [_row(f'C2_{X}', 'C2', X, g, **kw) for X, g in zip('CD', c2)]
    return pd.DataFrame(rows)


def test_verdict_passes():
    v = E86.verdict(_table(), FLIES)
    assert v['outcome'] == 'passes', v
    assert v['controls']['hold'] and v['gate']['holds'] and not v['controls']['b3_needed']
    assert v['replication']['mean_g'] == pytest.approx((0.0025 + 0.002 + 0.003 + 0.0032) / 4)
    assert v['gate']['mean_g'] == pytest.approx((0.0036 + 0.0031 + 0.003 + 0.0032) / 4)
    assert v['gate']['n_singles_lower'] == 12 and v['gate']['n_singles'] == 16


def test_verdict_pass_needs_three_positive_and_the_mean():
    # a mean of 0.002 or more with three positive passes (C2 kept high, so the gate holds); two positive does not even
    # carry
    v = E86.verdict(_table(c1=(0.0005, 0.0), c2=(0.004, 0.0039)), FLIES)
    assert v['outcome'] == 'passes' and v['replication']['n_positive'] == 3
    assert E86.verdict(_table(c1=(0.006, 0.004), c2=(-0.001, -0.001)), FLIES)['outcome'] == 'falls'
    # below 0.002 but 0.001 or more, three positive: carries
    assert E86.verdict(_table(c1=(0.0016, 0.001), c2=(0.0015, 0.0)), FLIES)['outcome'] == 'carries'
    # below 0.001: falls (the controls holding: outcome (c))
    v = E86.verdict(_table(c1=(0.0008, 0.0004), c2=(0.0006, -0.0002)), FLIES)
    assert v['outcome'] == 'falls' and v['controls']['hold']


@pytest.mark.parametrize('change, expect', [
    (dict(c2=(0.0005, 0.0012)), 'carries'),                 # (i) mean g of the gate's four 0.0021 < 0.003
    (dict(kw_singles=2), 'carries'),                        # (ii) 2 + 2 + 2 + 2 of 16
    (dict(kw_flagged=1), 'carries'),                        # (iii) a point outside the bands
    (dict(kw_raw=0.006), 'carries'),                        # (iii) the raw-count loss up by more than 0.005
    (dict(kw_mean_effect=0.19), 'carries'),                 # (iii) mean APL effect 0.012 from the flies'
    (dict(kw_sq_all=0.14), 'carries'),                      # (iv) 14 % against all recordings
    (dict(kw_sq_six=0.07), 'carries'),                      # (iv) 7 % against the six typical
])
def test_verdict_gate_conditions(change, expect):
    kw = {k[3:]: v for k, v in change.items() if k.startswith('kw_')}
    rest = {k: v for k, v in change.items() if not k.startswith('kw_')}
    v = E86.verdict(_table(**rest, **kw), FLIES)
    assert not v['gate']['holds']
    assert v['outcome'] == expect, v['gate']


def test_verdict_controls():
    # (v): a permuted input at half the candidate's g or more fails the controls (outcome (d))
    v = E86.verdict(_table(b1=(-0.003, 0.0018, -0.002, -0.005)), FLIES)
    assert v['controls']['v'] is False and v['outcome'] == 'incomplete'      # it is above 0: B3 is required
    v = E86.verdict(_table(b1=(-0.003, 0.0018, -0.002, -0.005), b3=(-0.001, -0.002, -0.003, -0.001)), FLIES)
    assert v['controls']['v'] is False and v['controls']['hold'] is False and v['outcome'] == 'falls'
    # in half B's direction the threshold is half of 0.0031
    v = E86.verdict(_table(b1=(-0.003, -0.002, 0.0016, -0.005), b3=(-0.001, -0.002, -0.003, -0.001)), FLIES)
    assert v['controls']['v'] is False and v['outcome'] == 'falls'
    # (vi): an inverted input above 0 fails; exactly 0 holds
    assert E86.verdict(_table(b2=(-0.01, 0.0001)), FLIES)['outcome'] == 'falls'
    v = E86.verdict(_table(b2=(-0.01, 0.0)), FLIES)
    assert v['controls']['vi'] and v['outcome'] == 'passes'


def test_verdict_b3():
    b1 = (-0.003, 0.0005, -0.002, -0.005)
    # a B1 input above 0: B3 is required
    v = E86.verdict(_table(b1=b1), FLIES)
    assert v['controls']['b3_needed'] and v['outcome'] == 'incomplete' and any('B3' in m for m in v['missing'])
    # B3 run: the candidate above all four of its direction and their mean below a quarter of it
    v = E86.verdict(_table(b1=b1, b3=(-0.001, -0.002, -0.003, -0.001)), FLIES)
    assert v['controls']['b3'] and v['controls']['hold'] and v['outcome'] == 'passes'
    # a B3 input above the candidate's g fails
    v = E86.verdict(_table(b1=b1, b3=(0.0037, -0.002, -0.003, -0.001)), FLIES)
    assert v['controls']['b3'] is False and v['outcome'] == 'falls'
    # the four of direction A with a mean of a quarter of 0.0036 or more fail (each below the candidate's)
    v = E86.verdict(_table(b1=(0.0017, 0.0017, -0.002, -0.005), b3=(0.0008, 0.0, -0.003, -0.001)), FLIES)
    assert v['controls']['v'] and v['controls']['b3_rows']['A']['mean_below'] is False and v['outcome'] == 'falls'


def test_verdict_incomplete_and_replicate_kinds():
    t = _table()
    v = E86.verdict(t[t.group != 'C2'], FLIES)
    assert v['outcome'] == 'incomplete' and any(m.startswith('C2') for m in v['missing'])
    v = E86.verdict(_table(c1=(0.003, 0.002), c2=(-0.001, 0.0004)), FLIES)
    assert v['replicate_kinds']['C1']['carries'] and not v['replicate_kinds']['C2']['carries']
