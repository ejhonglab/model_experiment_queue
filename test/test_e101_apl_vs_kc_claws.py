"""Tests of E101's rules on made-up data (experiments/E101_apl_vs_kc_claws.py): the groups, the stimulus design and the date
notes' remarks, the claw totals and the family's line and ratios, the leak's slope and its removal, the place test (a
planted local coupling passes, a leak alone fails at b_G), X3's verdict, the held rule, X1's classes, X2's count rule and Q
over the indicator combinations, the readers' release and raw-only refusals, the cell grid, and (where the shared directory
is mounted) the two-channel reader against E100's on a released 2021 recording."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from experiments import E101_apl_vs_kc_claws as M


def _shared():
    try:
        from kcq import harness as H
        return H.shared_dir()
    except FileNotFoundError:
        return None


# ---------------------------------------------------------------------------
# the groups

def test_groups_as_the_entry_counts_them():
    n = {g: len(v) for g, v in M.GROUPS.items()}
    assert n == {'L': 20, 'L+': 9, 'V0': 23, 'VT': 33, 'G': 6, 'R': 13, 'M': 1}
    stems = [s for v in M.GROUPS.values() for s in v] + list(M.LEFT_OUT)
    assert len(stems) == len(set(stems))
    assert len({M.fly_of(s) for s in M.GROUP_L}) == 7 and len({M.calyx_of(s) for s in M.GROUP_L}) == 7
    assert len({M.fly_of(s) for s in M.GROUP_V0}) == 5 and len({M.calyx_of(s) for s in M.GROUP_V0}) == 8
    assert {M.fly_of(s) for s in M.GROUP_L_PLUS} <= {M.fly_of(s) for s in M.GROUP_L}
    assert len({M.fly_of(s) for s in M.GROUP_L_PLUS}) == 5
    assert set(M.UNPROBED_X12) <= {M.fly_of(s) for s in M.GROUP_L} and len(M.UNPROBED_X12) == 5
    assert set(M.UNPROBED_X3) <= {M.fly_of(s) for s in M.GROUP_L} and len(M.UNPROBED_X3) == 3
    assert M.group_of('20221207_1FLfun2') == 'L' and M.group_of('20220719_2FRfun1') == 'left out'
    assert M.group_of('20250115_1FRfun1') is None
    assert M.fly_of('20221207_1FLfun2') == '20221207_1' and M.calyx_of('20221207_1FLfun2') == '20221207_1L'
    assert M.fly_of('20230126_1FLfun1_000') == '20230126_1' and M.fly_of('20220908_1Ffun1') == '20220908_1'
    assert M.X1_BOUND == pytest.approx(0.585, abs=1e-3)
    assert M.B_MULTIPLES == (0.0, 1.0, 2.0)


# ---------------------------------------------------------------------------
# A1: the design and the date notes

def _frames(blocks, nstim=5):
    """Presentations from blocks given as the stimulus of each block ('pfo' for an all-solvent block): each block one
    solvent then nstim - 1 repeats, as set L."""
    rows = []
    k = 0
    for b, st in enumerate(blocks):
        for pos in range(1, nstim + 1):
            if st == 'pfo' or pos == 1:
                o, c = 'solvent', np.nan
            else:
                o, c = st
            rows.append(dict(recording='r', frame=k, block=b, position=pos, odor=o, log10_conc=c))
            k += 1
    return pd.DataFrame(rows)


def test_block_design():
    f = _frames(['pfo', ('paa', -6.0), ('pa', -6.214), ('paa', -5.0), ('pa', -2.0), 'pfo'])
    d = M.block_design(f)
    assert d['design'] == 'solvent_then_repeats' and d['order'] == 'ascending'
    assert d['stimuli'] == ['paa @ -6', 'pa @ -6.214', 'paa @ -5', 'pa @ -2']
    assert d['presentations']['paa @ -6'] == 4 and d['n_solvent'] == 14
    f = _frames(['pfo', ('paa', -5.0), ('paa', -6.0)])
    assert M.block_design(f)['order'] == 'shuffled'
    f.loc[7, 'odor'] = 'solvent'
    f.loc[7, 'log10_conc'] = np.nan
    assert M.block_design(f)['design'] == 'other'


NOTE = """20221208
MaiTai 1000@0.87

1FRfun1 256 x 224
pockels: 20
gain: A: 100 B:15
z at 12 um
const int channelA[] =
{37,37,37,
37,49,49}
seq11

warpped a lot maybe ~9um
re-adjust

Change to 925 for testing!!!!
1FRfun2 256 x 224
pockels: 25
const int channelA[] = {
37,44,44}
re-run as 1FRfun3

2FLfun1 256 x 224 Maitai at 1040nm
const int channelA[] = {37}
"""


def test_date_note_runs_and_warp():
    r = M.date_note_runs(NOTE)
    assert r['1FRfun1']['wavelength_nm'] == 1000
    assert M.warp_um(' | '.join(r['1FRfun1']['remarks'])) == 9.0
    assert r['1FRfun2']['wavelength_nm'] == 925          # the change written before the run applies to it
    assert r['1FRfun3']['wavelength_nm'] == 925 and 're-run of 1FRfun2' in r['1FRfun3']['remarks']
    assert r['2FLfun1']['wavelength_nm'] == 1040
    assert M.warp_um('warpped for about <3um') == 3.0
    assert np.isnan(M.warp_um('significant warpping'))
    assert M.warp_um('readjust z. run again at 15um') is None


# ---------------------------------------------------------------------------
# the model's claws, line and ratios

def _stage(n_claws=6, stims=('pa', 'ms', 'paa', 'blank @ 0'), scale=1.0, seed=0):
    rng = np.random.default_rng(seed)
    kc = pd.MultiIndex.from_arrays([np.repeat([1, 2, 3], 2), np.arange(n_claws)], names=['kc_id', 'claw'])
    pre = pd.DataFrame(np.full((n_claws, len(stims)), 2.0), index=kc, columns=list(stims))
    drive = {'pa': 3.0, 'ms': 1.0, 'paa': 2.0, 'blank @ 0': 0.0}
    win = pre + pd.DataFrame(np.outer(np.ones(n_claws), [drive[s] * scale for s in stims]), index=kc, columns=list(stims))
    counts = pd.DataFrame(rng.integers(0, 4, size=(3, len(stims))).astype(float),
                          index=pd.Index([1, 2, 3], name='kc_id'), columns=list(stims))
    return dict(claw_win=win, claw_pre=pre, spike_counts=counts)


def test_claw_totals_drive_and_spike_versions():
    on = _stage()
    gl = np.array(['DM1', 'DM1', 'VA2', 'VA2', 'X', 'X'])
    t = M.claw_totals(on, None, gl, np.repeat([1, 2, 3], 2), ['DM1', 'VA2'], ['pa', 'ms', 'paa'])
    t = t.set_index('stimulus')
    assert t.loc['pa', M.claw_col('input', 0.0, 'change')] == pytest.approx(4 * 3.0)
    assert t.loc['pa', M.claw_col('all', 0.0, 'change')] == pytest.approx(6 * 3.0)
    assert t.loc['blank @ 0', M.claw_col('all', 0.0, 'change')] == pytest.approx(0.0)
    # the spike version: (1 - s) drive / scale0 + s counts / scale1, the spike signal without a pre-odor part
    s0, s1 = t['claw_scale_drive'].iloc[0], t['claw_scale_spikes'].iloc[0]
    k = on['spike_counts'].loc[[1, 1, 2, 2, 3, 3], 'pa'].to_numpy()
    want = (0.5 * 3.0 * 6 / s0) + 0.5 * k.sum() / s1
    assert t.loc['pa', M.claw_col('all', 0.5, 'change')] == pytest.approx(want)
    assert s0 == pytest.approx(np.sqrt(np.mean(np.array([3.0, 1.0, 2.0]) ** 2)))


def _fam(exponent=1.0, odor_factor=None):
    """A made-up family table: APL's evoked activity a power of the claws' change, optionally times an odor factor."""
    odor_factor = odor_factor or {}
    rows = []
    for w in ('real', 'P1'):
        for i, o in enumerate(M.H.APL16_ODORS):
            c = 100.0 * (1 + i)
            ev = (c / 100.0) ** exponent * odor_factor.get(o, 1.0)
            for cond in ('on', 'off'):
                rows.append(dict(wiring=w, cond=cond, stimulus=o, apl_win=8.0 + ev, apl_pre=8.0, apl_pre_settled=8.5,
                                 apl_evoked=ev, pn_total=10.0 + i, kc_rate=0.05 + 0.001 * i,
                                 **{M.claw_col(sub, sh, k): (c + 50 if k == 'win' else 50.0 if k == 'pre' else c)
                                    for sub in M.CLAW_SUBSETS for sh in M.CLAW_SHARES for k in ('win', 'pre', 'change')}))
    return pd.DataFrame(rows)


def test_family_line_and_ratios():
    ln = M.family_line(_fam(0.5))
    assert ln['slope_rma'].iloc[0] == pytest.approx(0.5) and np.allclose(ln['residual_log2'], 0, atol=1e-12)
    ln = M.family_line(_fam(1.0, {'pa': 2.0}))
    res = ln.set_index(['wiring', 'odor'])['residual_log2']
    assert res[('real', 'pa')] - res[('real', 'ms')] == pytest.approx(1.0, abs=0.2)
    r = M.family_ratios(_fam(1.0))
    i_pa, i_ms = M.H.APL16_ODORS.index('pa'), M.H.APL16_ODORS.index('ms')
    want = (1 + i_pa) / (1 + i_ms)
    apl = r[(r['wiring'] == 'real') & (r['pair'] == 'pa/ms') & (r['quantity'] == 'APL')]
    claws = r[(r['wiring'] == 'real') & (r['pair'] == 'pa/ms') & (r['quantity'] == 'claws') & (r['share'] == 0.0)]
    assert apl['linear'].iloc[0] == pytest.approx(want) and claws['linear'].iloc[0] == pytest.approx(want)
    assert apl['indicator_min'].iloc[0] <= want <= apl['indicator_max'].iloc[0]


def _fam_blank():
    f = _fam(1.0)
    return pd.concat([f, f[f['stimulus'] == 'pa'].assign(stimulus=M.BLANK_LABEL, apl_win=8.6, apl_evoked=0.6,
                                                         **{M.claw_col(sub, sh, 'win'): 60.0 for sub in M.CLAW_SUBSETS
                                                            for sh in M.CLAW_SHARES})])


def test_ratio_sets_baselines():
    f = _fam_blank()
    rs = M.ratio_sets(f[f['wiring'] == 'real'], 'pa', 'ms', baselines=('e90', 'settled', 'blank'))
    lin = rs[(rs['n'] == 1.0) & (rs['form'] == 'evoked')]
    a = lin[lin['quantity'] == 'APL'].set_index('baseline')['ratio']
    pa, ms = f[(f.wiring == 'real') & (f.stimulus == 'pa')].iloc[0], f[(f.wiring == 'real') & (f.stimulus == 'ms')].iloc[0]
    assert a['e90'] == pytest.approx((pa.apl_win - 8.0) / (ms.apl_win - 8.0))
    assert a['settled'] == pytest.approx((pa.apl_win - 8.5) / (ms.apl_win - 8.5))
    assert a['blank'] == pytest.approx((pa.apl_win - 8.6) / (ms.apl_win - 8.6))
    # the claws: at E90's pre-odor mean under 'e90', at the blank's window under 'blank' and under the partner 'settled'
    c = lin[(lin['quantity'] == 'claws') & (lin['subset'] == 'input') & (lin['share'] == 0.0)].set_index('baseline')['ratio']
    col = M.claw_col('input', 0.0, 'win')
    assert c['e90'] == pytest.approx((pa[col] - 50.0) / (ms[col] - 50.0))
    assert c['blank'] == pytest.approx((pa[col] - 60.0) / (ms[col] - 60.0)) == pytest.approx(c['settled'])


def test_blank_level_zero_and_by_hand():
    f = _fam_blank()
    g = f[f['wiring'] == 'real']
    # under the blank's level the blank's response is zero at both stages, in both forms
    z = M.ratio_sets(g, M.BLANK_LABEL, 'ms', baselines=('blank',))
    assert np.allclose(z['ratio'], 0.0)
    # both forms give the ratio computed by hand
    rs = M.ratio_sets(g, 'pa', 'ms', baselines=('blank',)).set_index(['quantity', 'subset', 'share', 'n', 'form'])['ratio']
    pa, ms = g[g.stimulus == 'pa'].iloc[0], g[g.stimulus == 'ms'].iloc[0]
    a1, a2, a0 = pa.apl_win, ms.apl_win, 8.6
    apl = M.ratio_sets(g, 'pa', 'ms', baselines=('blank',))
    apl = apl[apl['quantity'] == 'APL'].set_index(['n', 'form'])['ratio']
    assert apl[(2.0, 'evoked')] == pytest.approx(((a1 - a0) / (a2 - a0)) ** 2)
    assert apl[(2.0, 'baseline')] == pytest.approx((a1 ** 2 - a0 ** 2) / (a2 ** 2 - a0 ** 2))
    col = M.claw_col('input', 0.25, 'win')
    c1, c2, c0 = pa[col], ms[col], 60.0
    assert rs[('claws', 'input', 0.25, 1.5, 'evoked')] == pytest.approx(((c1 - c0) / (c2 - c0)) ** 1.5)
    assert rs[('claws', 'input', 0.25, 1.5, 'baseline')] == pytest.approx((c1 ** 1.5 - c0 ** 1.5) / (c2 ** 1.5 - c0 ** 1.5))


def test_run_without_a_blank_refused():
    g = _fam(1.0)
    g = g[g['wiring'] == 'real']
    for bl in ('blank', 'settled'):
        with pytest.raises(ValueError, match='blank'):
            M.ratio_sets(g, 'pa', 'ms', baselines=(bl,))
    assert len(M.ratio_sets(g, 'pa', 'ms', baselines=('e90',)))


def test_partner_clause():
    # APL's or Q's class that the partner level does not give is reading-dependent: not decided
    assert M.partner_held('inside', ['inside', 'inside'], 'APL') == 'inside'
    assert M.partner_held('inside', ['inside', 'flatter'], 'APL') == 'not decided'
    assert M.partner_held('flatter', ['steeper', 'flatter'], 'Q') == 'not decided'
    assert M.partner_held('flatter', ['not decided', 'flatter'], 'Q') == 'not decided'
    # the claws have no partner level
    assert M.partner_held('flatter', ['steeper', 'steeper'], 'claws') == 'flatter'
    assert M.partner_held('not decided', ['inside', 'inside'], 'APL') == 'not decided'


def _a3_run(blank_win=8.78, settled=8.79, rate=0.0, vary=0.0):
    stims = ['pa', 'ms', M.BLANK_LABEL]
    kc = pd.MultiIndex.from_arrays([np.repeat([1, 2], 2), np.arange(4)], names=['kc_id', 'claw'])
    run = dict(meta=dict(claw_glomerulus=np.array(['DM1', 'DM1', 'VA2', 'X']), kc_of_claw=np.repeat([1, 2], 2)))
    for _, ct in M.A3_CONDS:
        apl = pd.DataFrame(dict(win=[20.0, 19.0, blank_win], pre=8.27, pre_settled=[settled, settled + vary, settled]),
                           index=pd.Index(stims, name='stimulus'))
        pre = pd.DataFrame(np.full((4, 3), 100.0), index=kc, columns=stims)
        pre.iloc[0, 1] += vary * 1000
        resp = pd.DataFrame(np.zeros((2, 3)), columns=stims)
        resp.iloc[0, 2] = rate
        run[ct] = dict(apl=apl, claw_pre=pre, responses=resp)
    return run


def test_a3_blank_checks():
    rows = M.a3_blank_checks('real', _a3_run())
    assert len(rows) == 2 and all(r['passed'] for r in rows)
    assert not any(r['passed'] for r in M.a3_blank_checks('real', _a3_run(blank_win=8.9)))      # 1.3 % off
    assert not any(r['passed'] for r in M.a3_blank_checks('real', _a3_run(vary=1e-3)))          # levels differ
    assert not any(r['passed'] for r in M.a3_blank_checks('real', _a3_run(rate=1.0)))           # the blank recruits KCs


# ---------------------------------------------------------------------------
# the leak, its removal and the place test

def test_planted_leak_measured_and_removed():
    rng = np.random.default_rng(1)
    apl = rng.uniform(0.2, 1.0, size=(16, 1))[:, [0] * 50] * (1 + 0.3 * rng.standard_normal((16, 50)))
    kc_green_only = 0.034 * apl + rng.normal(0, 0.002, apl.shape)
    b = np.median(M.leak_slopes(apl, kc_green_only))
    assert b == pytest.approx(0.034, abs=0.003)
    kc = 0.2 * rng.uniform(0.2, 1.0, size=(16, 50)) + 0.034 * apl
    assert np.allclose(M.remove_leak(kc, apl, 0.034), kc - 0.034 * apl)
    assert np.allclose(M.remove_leak(kc, apl, 0.0), kc)


def _x3_flies(true_corr, leak, n_flies=7, r=0.9, seed=0):
    rng = np.random.default_rng(seed)
    rec = M.made_up_recordings(rng, n_flies * 2, r, r, true_corr, leak)
    rows = []
    for i in range(n_flies * 2):
        for bm in M.B_MULTIPLES:
            ha, hk = M.halves(rec['apl'][i]), M.halves(rec['kc'][i])
            hk = tuple(M.remove_leak(k, a, bm * leak) for k, a in zip(hk, ha))
            x = M.x3_recording(ha, hk, rec['odors'], n_perm=200, seed=i)
            rows.append(dict(fly=f'f{i // 2}', b_mult=bm, **x))
    return pd.DataFrame(rows)


def test_local_coupling_passes_and_leak_alone_fails():
    fx = M.fly_x3(_x3_flies(0.9, M.A4_B_G))
    assert (fx[fx['b_mult'] == 1.0]['cls'] == 'pass').all()
    fx = M.fly_x3(_x3_flies(0.0, M.A4_B_G, seed=3))
    assert (fx[fx['b_mult'] == 1.0]['cls'] == 'fail').all()        # a leak alone does not pass once removed
    assert (fx[fx['b_mult'] == 1.0]['cls'] != 'pass').all()


def test_x3_recording_statistics():
    rng = np.random.default_rng(2)
    rec = M.made_up_recordings(rng, 1, 0.999, 0.999, 1.0, 0.0)
    ha, hk = M.halves(rec['apl'][0]), M.halves(rec['kc'][0])
    x = M.x3_recording(ha, hk, rec['odors'], n_perm=100)
    assert x['A'] > 0.95 and x['ceiling'] > 0.95 and x['null_pct'] < x['A']
    assert x['n_stimuli'] == sum(M.A4_ODORS) and x['n_odors'] == len(M.A4_ODORS)
    perms = M.other_odor_assignments(rec['odors'], 50)
    od = np.asarray(rec['odors'])
    assert (od[perms] != od[None, :]).all()
    with pytest.raises(ValueError):
        M.other_odor_assignments(['a', 'a'], 5)


def _classes(cls_by_b, flies=('20221016_3', '20221016_4', '20221207_1', '20221207_2', '20221208_1', '20221213_2',
                               '20221217_2')):
    return pd.DataFrame([dict(fly=f, b_mult=b, cls=c[i] if isinstance(c, (list, tuple)) else c)
                         for b, c in cls_by_b.items() for i, f in enumerate(flies)])


def test_x3_verdict():
    allpass = {0.0: 'pass', 1.0: 'pass', 2.0: 'pass'}
    assert M.x3_verdict(_classes(allpass), True)['verdict'] == 'follows'
    assert M.x3_verdict(_classes(allpass), False)['verdict'] == 'not decided'       # B1's cross-check failed
    leak = {0.0: 'pass', 1.0: 'pass', 2.0: 'neither'}
    assert M.x3_verdict(_classes(leak), True)['verdict'] == "follows at the controls' leak"
    assert M.x3_verdict(_classes(leak), True, b_g_alone=True)['verdict'] == 'follows'
    assert M.x3_verdict(_classes({0.0: 'fail', 1.0: 'fail', 2.0: 'fail'}), True)['verdict'] == 'does not follow'
    # five pass, but only one of the three unopened flies (2022-10-16 fly 4, 2022-12-07 fly 2, 2022-12-08 fly 1)
    c = ['pass', 'fail', 'pass', 'fail', 'pass', 'pass', 'pass']
    assert M.x3_verdict(_classes({0.0: c, 1.0: c, 2.0: c}), True)['verdict'] == 'not decided'
    c = ['pass', 'pass', 'fail', 'pass', 'fail', 'pass', 'pass']
    assert M.x3_verdict(_classes({0.0: c, 1.0: c, 2.0: c}), True)['verdict'] == 'follows'
    few = _classes(allpass, flies=('20221016_4', '20221207_2', '20221208_1', '20221213_2'))
    assert M.x3_verdict(few, True)['verdict'] == 'not decided'


def test_fly_x3_rule():
    pr = pd.DataFrame([dict(fly='a', b_mult=0.0, A=0.35, null_pct=0.1, repeat_r_apl=0.8, repeat_r_kc=0.8),
                       dict(fly='b', b_mult=0.0, A=0.1, null_pct=0.0, repeat_r_apl=0.8, repeat_r_kc=0.8),
                       dict(fly='c', b_mult=0.0, A=0.1, null_pct=0.0, repeat_r_apl=0.4, repeat_r_kc=0.8),
                       dict(fly='d', b_mult=0.0, A=0.32, null_pct=0.4, repeat_r_apl=0.8, repeat_r_kc=0.8),
                       dict(fly='e', b_mult=0.0, A=0.2, null_pct=0.1, repeat_r_apl=0.8, repeat_r_kc=0.8)])
    c = M.fly_x3(pr).set_index('fly')['cls']
    assert c.to_dict() == {'a': 'pass', 'b': 'fail', 'c': 'neither', 'd': 'fail', 'e': 'neither'}


def test_held_rule():
    k = [(b, r) for b in (0.0, 1.0) for r in M.READINGS]
    assert M.held_class({x: 'inside' for x in k}, 'inside', 5, M.HELD_X12_MIN) == 'inside'
    assert M.held_class({x: 'inside' for x in k}, 'inside', 3, M.HELD_X12_MIN) == 'not decided'
    assert M.held_class({x: 'inside' for x in k}, 'flatter', 5, M.HELD_X12_MIN) == 'not decided'
    d = {x: 'inside' for x in k}
    d[(1.0, 'block')] = 'flatter'
    assert M.held_class(d, 'inside', 5, M.HELD_X12_MIN) == 'not decided'
    assert M.held_class({x: 'not decided' for x in k}, 'not decided', 5, 4) == 'not decided'
    assert M.held_class({'carrying': 'follows', 'block': 'follows'}, None, 0, 0) == 'follows'


# ---------------------------------------------------------------------------
# X1 and X2

def _resp(factor_pa=1.0, slope=1.0, seed=0):
    rng = np.random.default_rng(seed)
    stims = [f'pa @ {c:g}' for c in (-5.0, -4.4, -3.8, -2.0)] + [f'ms @ {c:g}' for c in (-4.3, -3.6, -2.8, -2.0)] + \
            [f'paa @ {c:g}' for c in (-5.0, -4.0, -3.0, -2.0)]
    kc = pd.Series(rng.uniform(0.1, 1.0, len(stims)), index=stims)
    apl = kc ** slope * np.where([s.startswith('pa @') for s in stims], factor_pa, 1.0)
    return apl, kc


def test_x1_recording_and_verdict():
    apl, kc = _resp(1.0)
    x = M.x1_recording(apl, kc)
    assert all(abs(v) < 1e-9 for v in x['d'].values()) and x['slope'] == pytest.approx(1.0)
    apl, kc = _resp(2.0)
    x = M.x1_recording(apl, kc)
    assert x['d']['pa/ms'] == pytest.approx(1.0, abs=0.35) and x['d']['paa/pa'] < 0
    # a recording with one kept ms stimulus enters no pair with ms
    keep = [s for s in apl.index if not s.startswith('ms') or s == 'ms @ -2']
    x = M.x1_recording(apl[keep], kc[keep])
    assert np.isnan(x['d']['pa/ms']) and np.isnan(x['d']['paa/ms']) and np.isfinite(x['d']['paa/pa'])
    n5 = {f'{a}/{b}': 5 for a, b in M.X1_PAIRS}
    inside = {p: (-0.3, 0.4) for p in n5}
    assert M.x1_verdict(inside, n5) == 'one function'
    assert M.x1_verdict(dict(inside, **{'pa/ms': (0.7, 1.2)}), n5) == 'depends on the odor'
    assert M.x1_verdict(dict(inside, **{'pa/ms': (-1.2, -0.6)}), n5) == 'depends on the odor'
    assert M.x1_verdict(dict(inside, **{'pa/ms': (0.2, 0.9)}), n5) == 'not decided'
    assert M.x1_verdict(inside, dict(n5, **{'paa/pa': 4})) == 'not decided'


def test_fly_mean_interval_and_ratios():
    m, lo, hi = M.fly_mean_interval(np.array([1.0, 2.0, 3.0]), n_boot=2000)
    assert m == pytest.approx(2.0) and lo < m < hi
    pr = pd.DataFrame(dict(fly=['a', 'a', 'b'], num=[0.4, 0.9, 0.6], den=[0.2, 0.3, 0.2]))
    r = M.flies_ratio(pr, n_boot=200)
    assert r['ratio'] == pytest.approx(np.exp(np.mean([np.mean(np.log([2.0, 3.0])), np.log(3.0)])))
    assert r['n_flies'] == 2 and r['n_left_out'] == 0
    pq = pd.DataFrame(dict(fly=['a', 'b'], apl_num=[0.6, 0.8], apl_den=[0.3, 0.4], kc_num=[0.2, 0.3], kc_den=[0.2, 0.1]))
    q = M.flies_q(pq, n_boot=200)
    assert q['ratio'] == pytest.approx(np.exp(np.mean([np.log(2.0), np.log(2.0 / 3.0)])))


def test_x2_count_rule():
    assert M.x2_quantity_verdict({'r1': 'flatter', 'r2': 'flatter', 'r3': 'inside'}) == 'flatter'
    assert M.x2_quantity_verdict({'r1': 'flatter', 'r2': 'flatter', 'r3': 'steeper'}) == 'not decided'
    assert M.x2_quantity_verdict({'r1': 'inside', 'r2': 'inside', 'r3': 'steeper'}) == 'inside'
    assert M.x2_quantity_verdict({'r1': 'steeper', 'r2': 'steeper', 'r3': 'not decided'}) == 'steeper'
    assert M.x2_quantity_verdict({'r1': 'inside', 'r2': 'flatter', 'r3': 'steeper'}) == 'not decided'


def test_q_over_indicator_combinations():
    rows = []
    for w, apl_lin, claw_lin in (('real', 1.2, 1.0), ('P1', 3.0, 1.0)):
        for form in M.FORMS:
            for n in M.APL_NS:
                rows.append(dict(wiring=w, estimate='ms', baseline='e90', quantity='APL', subset='', share=np.nan, n=n,
                                 form=form, ratio=apl_lin ** n))
            for sh in M.CLAW_SHARES:
                for n in M.CLAW_NS:
                    rows.append(dict(wiring=w, estimate='ms', baseline='e90', quantity='claws', subset='input', share=sh,
                                     n=n, form=form, ratio=claw_lin))
    sets = pd.DataFrame(rows)
    q = M.q_classes(sets, 1.5, 2.0).set_index('wiring')
    assert q.loc['real', 'cls'] == 'inside'          # 1.2^2.5 = 1.58 lies in [1.5, 2.0]: one combination is enough
    assert q.loc['real', 'q_linear'] == pytest.approx(1.2)
    assert q.loc['P1', 'cls'] == 'steeper'           # 3.0 and above: every combination beyond the interval
    # an interval below 1 with every Q above 1: on the other side of 1, which E98's classify calls flatter
    assert M.q_classes(sets, 0.5, 0.9).set_index('wiring').loc['real', 'cls'] == 'flatter'


# ---------------------------------------------------------------------------
# the stimulus responses, the kept stimuli, the cell grid

def test_stimulus_responses_readings_and_kept():
    f = _frames(['pfo', ('pa', -2.0), ('ms', -2.0), 'pfo'])
    f['is_solvent'] = f['odor'] == 'solvent'
    f['stim'] = ['solvent' if s else f'{o} @ {c:g}' for s, o, c in zip(f['is_solvent'], f['odor'], f['log10_conc'])]
    v = np.zeros(len(f))
    v[(f['block'] == 1) & (f['position'] == 1)] = 0.1          # the pa block's own solvent
    v[(f['stim'] == 'pa @ -2')] = [1.0, 1.2, 1.4, 1.6]
    v[(f['stim'] == 'ms @ -2')] = 0.5
    c = M.stimulus_responses(f, v, 'carrying')['value']
    base = np.mean([0.0, 0.1, 0.0, 0.0])                       # the first-in-block solvents of the four blocks
    assert c['pa @ -2'] == pytest.approx(1.3 - base)
    assert M.stimulus_responses(f, v, 'block')['value']['pa @ -2'] == pytest.approx(1.3 - 0.1)
    assert M.stimulus_responses(f, v, 'first_two')['value']['pa @ -2'] == pytest.approx(1.1 - base)
    maps = np.repeat(v[:, None], 3, axis=1)
    assert np.allclose(M.stimulus_responses(f, maps, 'carrying').loc['pa @ -2', [0, 1, 2]], 1.3 - base)
    sol = np.array([0.0, 0.02, -0.02, 0.01])
    kept = M.kept_stimuli(pd.Series({'pa @ -2': 1.0, 'ms @ -2': 0.03}), pd.Series({'pa @ -2': 0.2, 'ms @ -2': 0.2}),
                          sol, sol)
    assert kept == ['pa @ -2']


def test_cell_grid_and_maps():
    mask = np.zeros((40, 40), bool)
    mask[4:36, 8:32] = True
    g = M.cell_grid(mask, pixel_um=0.5, cell_um=4.0)
    assert g['cell_px'] == 8 and g['shape'] == (4, 3) and g['n_cells'] == 12
    stack = np.stack([np.indices((40, 40))[0].astype(float)])          # each pixel's row
    cm = M.cell_maps(stack, g)
    assert cm.shape == (1, 12)
    assert cm[0, 0] == pytest.approx(np.arange(4, 12).mean())
    grid = M.cells_to_grid(cm[0], g)
    assert grid.shape == (4, 3) and np.isfinite(grid).all()


def test_experiment_xml_path_other_sex_letter(tmp_path):
    (tmp_path / '20221005_1FRfun2_Experiment.xml').write_text('<x/>')
    (tmp_path / '20221005_2FLfun1_Experiment.xml').write_text('<x/>')
    x, src = M.experiment_xml_path(tmp_path, '20221005_1MRfun2')
    assert x.name == '20221005_1FRfun2_Experiment.xml' and src == '20221005_1FRfun2'
    x, src = M.experiment_xml_path(tmp_path, '20221005_2FLfun1')
    assert src == '20221005_2FLfun1'
    assert M.experiment_xml_path(tmp_path, '20221005_1MRfun3') == (None, '')


def test_cross_correlogram_and_common_share():
    rng = np.random.default_rng(0)
    a = M.smooth_field(rng, 15, 1.5).reshape(15, 15)
    x = M.cross_correlogram(a, a, 4.0)
    assert x['peak_lag_cells'] == (0, 0) and x['zero_lag'] == pytest.approx(1.0) and x['length_um'] > 0
    y = M.cross_correlogram(np.roll(a, 3, axis=1), a, 4.0)
    assert abs(y['peak_lag_cells'][1]) == 3 and not y['peak_within_one_cell']
    P = np.tile(1 + 0.3 * a.ravel(), (5, 1))
    assert M.common_map_share(P)['common_share'] == pytest.approx(1.0)
    Q = 1 + 0.3 * M.smooth_field(rng, 15, 1.5, (5,))
    c = M.common_map_share(Q, Q)
    assert c['claws_share_of_rest'] == pytest.approx(1.0) and c['common_share'] < 0.5


def test_a4_rule():
    rows = [dict(r_apl=r, r_kc=r, true_share=s, b_mult=b, power_pass=p)
            for r in (0.3, 0.9) for s in (0.5,) for b, p in ((0.0, 0.9), (1.0, 0.85), (2.0, 0.6 if r == 0.3 else 0.9))]
    rule = M.a4_rule(pd.DataFrame(rows))
    assert rule['power_2bg_half_ceiling'] == pytest.approx(0.75) and rule['b_g_alone']
    rows = [dict(r, power_pass=0.9) for r in rows]
    assert not M.a4_rule(pd.DataFrame(rows))['b_g_alone']


# ---------------------------------------------------------------------------
# the refusals and the readers

RAW = "raw only: not used until further notice (request in viki's terminal 2026-11-15 15:28)"
OK = "processed: usable (Yang's processed products only)"


@pytest.mark.skipif(_shared() is None, reason='shared directory not mounted')
def test_refusals(monkeypatch):
    # the mechanism, whatever the module's current release state
    monkeypatch.setattr(M, 'YANG_RELEASED', False)
    monkeypatch.setattr(M, 'CHANNEL_MODE', 'waiting')
    row = dict(stem='20221207_1FLfun2', use_until_further_notice=OK, folder='/nonexistent')
    with pytest.raises(PermissionError, match='not released'):
        M.check_readable(row)
    with pytest.raises(PermissionError, match='not released'):
        M.read_two_channel(row)
    with pytest.raises(PermissionError, match='raw only'):
        M.check_readable(dict(row, use_until_further_notice=RAW))
    from experiments import E100_apl_weak_odors as E100
    with pytest.raises(PermissionError, match='raw only'):
        M.check_readable(dict(stem=E100.released_2021()[0], use_until_further_notice=RAW))
    for part in (M.part_b, M.part_c, M.part_d, M.part_e, M.score):
        with pytest.raises(NotImplementedError, match='release'):
            part()
    # released, but the channels neither confirmed nor read as inferred: still waiting
    monkeypatch.setattr(M, 'YANG_RELEASED', True)
    for part in (M.part_b, M.part_c, M.score):
        with pytest.raises(NotImplementedError, match='CHANNEL_MODE'):
            part()


def test_gate_needs_the_channel_check(tmp_path, monkeypatch):
    import json
    monkeypatch.setattr(M, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(M, 'YANG_RELEASED', True)
    monkeypatch.setattr(M, 'CHANNEL_MODE', M.INFERRED_CHECKED)
    M._gate('B')                                                   # part B runs before the check
    for part in ('C', 'D', 'E', 'score'):
        with pytest.raises(M.ChannelCheckNotPassed, match='not run'):
            M._gate(part)
    for v in ('swapped', 'failed'):
        (tmp_path / 'channel_check.json').write_text(json.dumps(dict(verdict=v)))
        with pytest.raises(M.ChannelCheckNotPassed, match=v):
            M._gate('C')
    (tmp_path / 'channel_check.json').write_text(json.dumps(dict(verdict='passed')))
    M._gate('C')
    monkeypatch.setattr(M, 'CHANNEL_MODE', 'confirmed by Yang')
    (tmp_path / 'channel_check.json').unlink()
    M._gate('D')                                                   # confirmed by Yang: no check needed


@pytest.mark.skipif(_shared() is None or not (_shared() / 'data' / 'yang_apl_recordings' / 'recordings.csv').exists(),
                    reason="viki's copy not there")
def test_two_channel_reader_on_a_released_recording():
    res = M.reader_identity_released()
    assert res['passed'], res


@pytest.mark.skipif(_shared() is None or not (_shared() / 'data' / 'yang_apl_recordings' / 'recordings.csv').exists(),
                    reason="viki's copy not there")
def test_a1_rows_from_recordings_csv():
    recs = pd.read_csv(M.copy_dir() / 'recordings.csv')
    t = M.a1_rows(recs)
    n = t.groupby('group').size()
    assert n['L'] == 20 and n['L+'] == 9 and n['V0'] == 23 and n['VT'] == 33 and n['G'] == 6 and n['R'] == 13
    assert n['M'] == 1 and n['left out'] == 1
    assert not t[t['group'] == 'raw only']['usable'].any()
    assert t[t['group'].isin(list(M.GROUPS))]['usable'].all()


# ---------------------------------------------------------------------------
# parts B to E and the score, end to end on made-up recordings (the readers replaced, the outputs in tmp_path)

L_CONCS = {'pa': (-6.214, -5.612, -5.01, -4.408, -3.806, -2.0), 'paa': (-6.0, -5.0, -4.0, -3.0, -2.0),
           'ms': (-5.112, -4.334, -3.556, -2.778, -2.0)}


def _made_up_recording(stem, group, rng, leak=0.034, coupling=0.9, green_only=False, red_only=False, n_px=60,
                       reverse_leak=0.01):
    """A made-up two-channel recording of set L's design: APL (C1) and the KCs (C2) with odor-specific patterns coupled
    by `coupling`, the leak of green into red `leak`; `green_only` (group G) without a KC response, `red_only` (group R)
    without an APL response, its C1 carrying `reverse_leak` of the KC channel."""
    stims = [(o, c) for o in L_CONCS for c in L_CONCS[o]]
    f = _frames(['pfo'] + stims + ['pfo'])
    f['is_solvent'] = f['odor'] == 'solvent'
    f['stim'] = ['solvent' if s else f'{o} @ {c:g}' for s, o, c in zip(f['is_solvent'], f['odor'], f['log10_conc'])]
    U = {o: M.smooth_field(rng, n_px, 6.0).reshape(n_px, n_px) for o in L_CONCS}
    V = {o: coupling * U[o] + np.sqrt(1 - coupling ** 2) * M.smooth_field(rng, n_px, 6.0).reshape(n_px, n_px)
         for o in L_CONCS}
    apl, kc = [], []
    for st, o, c in zip(f['is_solvent'], f['odor'], f['log10_conc']):
        if st:
            a = rng.normal(0, 0.02, (n_px, n_px))
            k = rng.normal(0, 0.005, (n_px, n_px))
        else:
            ks = 0.05 + 0.25 * (c - min(L_CONCS[o])) / (max(L_CONCS[o]) - min(L_CONCS[o]))
            a = (0.0 if red_only else 5 * ks * (1 + 0.3 * U[o])) + rng.normal(0, 0.03 if not red_only else 0.005,
                                                                               (n_px, n_px))
            k = (0.0 if green_only else ks * (1 + 0.3 * V[o])) + rng.normal(0, 0.005, (n_px, n_px))
        apl.append(a + (reverse_leak * k if red_only else 0.0))
        kc.append(k + leak * a)
    apl, kc = np.stack(apl), np.stack(kc)
    f['C1'], f['C2'] = apl.mean(axis=(1, 2)), kc.mean(axis=(1, 2))
    return dict(stem=stem, fly=M.fly_of(stem), calyx=M.calyx_of(stem), group=group, plane_um=15.0, frames=f,
                mask=np.ones((n_px, n_px), bool), pixel_um=1.0, xml=None, registration_matched=True,
                stacks={'C1': apl, 'C2': kc}, channels=['C1', 'C2'], truncated=False)


def _made_up_model_sets():
    rows = []
    stims = ['pa', 'ms', 'paa', 'paa [validation2/diag] @ -5', 'paa [kiwi-control/diag] @ -5', 'paa [megamat/diag] @ -5',
             'ms [kiwi-control/diag] @ -3', 'ms [megamat/diag] @ -3', 'ms [validation2/diag] @ -3',
             'ms [kiwi-control/control] @ -3', M.BLANK_LABEL]
    level = {'pa': 1.3, 'ms': 1.0, 'paa': 1.1, M.BLANK_LABEL: 0.0}
    for w in M.ALL_WIRINGS:
        for s in stims:
            v = level.get(s, 0.6 if s.startswith('paa') else 1.0)
            rows.append(dict(wiring=w, cond='on', stimulus=s, apl_win=8 + 11 * v, apl_pre=8.0, apl_pre_settled=8.5,
                             apl_evoked=11 * v, pn_total=1.0, kc_rate=0.05,
                             **{M.claw_col(sub, sh, k): (1000 + 5000 * v if k == 'win' else 1000.0 if k == 'pre'
                                                         else 5000 * v)
                                for sub in M.CLAW_SUBSETS for sh in M.CLAW_SHARES for k in ('win', 'pre', 'change')}))
    return M.model_x2(pd.DataFrame(rows))


R_TEST = ('20220715_4FRfun2', '20221005_2FLfun1', '20221005_2FLfun2', '20221005_2FRfun1')


def _swap(rec):
    """The same recording with its two channels exchanged (a wrong channel assignment)."""
    f = rec['frames'].copy()
    f['C1'], f['C2'] = rec['frames']['C2'].to_numpy(), rec['frames']['C1'].to_numpy()
    return dict(rec, frames=f, stacks={'C1': rec['stacks']['C2'], 'C2': rec['stacks']['C1']})


def _controls(rng, swap=False, dual_in_g=False):
    g = [_made_up_recording(s, 'G', rng, green_only=not dual_in_g, n_px=20) for s in M.GROUP_G]
    r = [_made_up_recording(s, 'R', rng, red_only=True, n_px=20) for s in R_TEST]
    if swap:
        g, r = [_swap(x) for x in g], [_swap(x) for x in r]
    return g, r


def test_channel_check_pass_swap_fail(tmp_path, monkeypatch):
    rng = np.random.default_rng(11)
    # a pass: the green-only flies respond in C1 with C2 at the leak, the red-only flies in C2 with C1 near zero
    g, r = _controls(rng)
    t, res = M.channel_check(g + r)
    assert res['verdict'] == 'passed', res
    assert (t['cls'] == 'as inferred').all()
    assert t[t['group'] == 'G']['expected_other_fraction'].max() < 0.06   # the planted leak, 0.034
    assert {f['fly'] for f in res['flies']} == {M.fly_of(s) for s in M.GROUP_G + R_TEST}
    # a swap: every control's channels exchanged
    g, r = _controls(rng, swap=True)
    t, res = M.channel_check(g + r)
    assert res['verdict'] == 'swapped' and (t['cls'] == 'swapped').all()
    # a fail: the 'green-only' flies carry a KC response in C2 (0.2 of APL's), so C2 is more than the leak
    g, r = _controls(rng, dual_in_g=True)
    t, res = M.channel_check(g + r)
    assert res['verdict'] == 'failed' and (t[t['group'] == 'G']['cls'] == 'neither').all()
    assert (t[t['group'] == 'R']['cls'] == 'as inferred').all()
    # a fail: no red-only fly with two channels
    g, _ = _controls(rng)
    one = [dict(x, frames=x['frames'].drop(columns='C1')) for x in _controls(rng)[1]]
    t, res = M.channel_check(g + one)
    assert res['verdict'] == 'failed' and (t[t['group'] == 'R']['cls'] == 'not evaluable').all()
    # one recording of three swapped in a green-only fly still leaves it 'as inferred' (2 of 3); two of three do not
    g, r = _controls(rng)
    gs = [_swap(x) if x['stem'] == '20221005_1MRfun2' else x for x in g]
    assert M.channel_check(gs + r)[1]['verdict'] == 'passed'
    gs = [_swap(x) if x['stem'] in ('20221005_1MRfun2', '20221005_1MRfun3') else x for x in g]
    assert M.channel_check(gs + r)[1]['verdict'] == 'failed'
    # part B writes the verdict and parts C to E then refuse
    recs = {'L': [_made_up_recording(s, 'L', rng, n_px=20) for s in M.GROUP_L]}
    g, r = _controls(rng, swap=True)
    recs.update(G=g, R=r)
    monkeypatch.setattr(M, 'read_groups', lambda groups, maps=True: [x for gg in groups for x in recs.get(gg, [])])
    monkeypatch.setattr(M, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(M, 'YANG_RELEASED', True)
    monkeypatch.setattr(M, 'CHANNEL_MODE', M.INFERRED_CHECKED)
    b = M.part_b(n_perm=50)
    assert b['channel_check'] == 'swapped'
    import json
    assert json.loads((tmp_path / 'channel_check.json').read_text())['verdict'] == 'swapped'
    assert (tmp_path / 'channel_check.csv').exists() and (tmp_path / 'b_leak.json').exists()
    for part in (M.part_c, M.part_d, M.part_e, M.score):
        with pytest.raises(M.ChannelCheckNotPassed, match='swapped'):
            part()


def test_channel_class_numbers():
    resp = {'C1': pd.Series({'pa @ -2': 1.0, 'ms @ -2': 0.5}), 'C2': pd.Series({'pa @ -2': 0.1, 'ms @ -2': 0.4}),
            'C1_solvent_sd': 0.01, 'C2_solvent_sd': 0.01}
    c = M.channel_class(resp, 'C1', 'C2')
    assert c['cls'] == 'as inferred' and c['expected_strongest'] == 'pa @ -2'
    assert c['expected_other_fraction'] == pytest.approx(0.1)
    assert c['swapped_strongest'] == 'ms @ -2' and c['swapped_other_fraction'] == pytest.approx(0.5 / 0.4)
    # the other channel at 0.2 of the expected one's: not the leak
    c = M.channel_class(dict(resp, C2=pd.Series({'pa @ -2': 0.2, 'ms @ -2': 0.1})), 'C1', 'C2')
    assert c['cls'] == 'neither'
    # the expected channel below twice its solvent SD, or below 0.02 dF/F: does not respond
    c = M.channel_class(dict(resp, C1_solvent_sd=0.6), 'C1', 'C2')
    assert not c['expected_responds'] and c['cls'] == 'neither'
    small = dict(resp, C1=resp['C1'] * 0.015, C2=resp['C2'] * 0.0)
    assert M.channel_class(dict(small, C1_solvent_sd=0.001), 'C1', 'C2')['cls'] == 'neither'


def test_parts_b_to_e_and_score_on_made_up_recordings(tmp_path, monkeypatch):
    import json
    rng = np.random.default_rng(5)
    recs = {'L': [_made_up_recording(s, 'L', rng) for s in M.GROUP_L],
            'G': [_made_up_recording(s, 'G', rng, green_only=True) for s in M.GROUP_G],
            'R': [_made_up_recording(s, 'R', rng, red_only=True, n_px=20) for s in R_TEST]}
    monkeypatch.setattr(M, 'read_groups', lambda groups, maps=True: [r for g in groups for r in recs.get(g, [])])
    monkeypatch.setattr(M, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(M, 'load_a3_run', lambda label, d=None: None)
    monkeypatch.setattr(M, 'e_remy_rows', lambda pairs=(('pa', 'ms'),): [])
    monkeypatch.setattr(M, 'YANG_RELEASED', False)
    with pytest.raises(NotImplementedError, match='release'):
        M.part_b()
    monkeypatch.setattr(M, 'YANG_RELEASED', True)
    monkeypatch.setattr(M, 'CHANNEL_MODE', M.INFERRED_CHECKED)
    (tmp_path / 'a4_rule.json').write_text(json.dumps(dict(b_g_alone=True)))
    b = M.part_b(n_perm=100)
    assert b['channel_check'] == 'passed'
    assert b['b_G'] == pytest.approx(0.034, abs=0.006)           # the planted leak, measured in the green-only flies
    assert b['cross_check_passed'] is True and not b['outcome_h'] and not b['outcome_h_b5']
    b2 = M.read_table(tmp_path / 'b2_reverse.csv')
    assert len(b2) and b2['share_of_set_l'].abs().max() < 0.1      # red-only flies: APL's channel far below set L's
    c = M.part_c(n_boot=500)
    x1 = json.loads((tmp_path / 'c_x1.json').read_text())
    assert x1['all']['carrying b1']['verdict'] == 'one function'
    M.part_d(n_perm=100, cell_ums=(4.0,))
    d = M.read_table(tmp_path / 'd_x3.csv')
    assert set(d['reading']) == set(M.READINGS) and set(d['b_mult']) == set(M.B_MULTIPLES)
    M.part_e()
    e = M.read_table(tmp_path / 'e_claw_distances.csv')
    assert len(e) and e['distance_noise_removed'].notna().any()
    pd.DataFrame(dict(stem=list(M.GROUP_L), out_of_sensitivity_readings=['warp 9 um'] + [''] * (len(M.GROUP_L) - 1))
                 ).to_csv(tmp_path / 'a1_recordings.csv', index=False)
    v = M.score(model_sets=_made_up_model_sets())
    assert v['X1']['held'] == 'one function'
    assert v['X3']['by_reading']['carrying']['verdict'] == 'follows' and v['X3']['held'] == 'follows'
    assert v['model_levels'] == dict(carrying='blank', partner='settled', shown_deciding_nothing='e90')
    assert v['caveat'] == M.QUALITY_CAVEAT
    for q in M.X2_QUANTITIES:
        assert v['X2']['held'][q] in ('flatter', 'inside', 'steeper', 'not decided')
        assert v['X2']['lines'][q].startswith(f'{q}: ') and "E90's level (decides nothing)" in v['X2']['lines'][q]
        qq = v['X2']['quantities'][q]
        assert qq['held'] == M.partner_held(qq['held_before_partner'], list(qq['partner_classes'].values()), q)
        assert qq['n_unprobed_flies_min'] >= 0 and isinstance(qq['unprobed_classes'], list)
    assert set(v['X2']['quantities']['APL']['partner_classes']) == {'carrying b0', 'carrying b1'}
    assert (tmp_path / 'verdict.json').exists()
    with open(tmp_path / 'score_x2.csv') as f:
        assert f.readline().strip() == f'# {M.QUALITY_CAVEAT}'
    t = M.read_table(tmp_path / 'score_x2.csv')
    assert set(t['a0']) == set(M.A0_CHOICES)
    assert set(t[t['a0'] == 'settled']['reading']) == {'carrying'}                  # the partner: the carrying reading only
    assert dict(t[['a0', 'level_role']].drop_duplicates().values.tolist()) == M.A0_ROLE
    for name in ('channel_check.csv', 'b1_slopes.csv', 'c_responses.csv', 'c_x2_flies.csv', 'd_x3.csv',
                 'e_claw_distances.csv', 'e_model_distances.csv', 'score_x3_flies.csv'):
        with open(tmp_path / name) as f:
            assert f.readline().strip() == f'# {M.QUALITY_CAVEAT}', name
    for name in ('channel_check.json', 'b_leak.json', 'c_x1.json', 'verdict.json'):
        assert json.loads((tmp_path / name).read_text())['caveat'] == M.QUALITY_CAVEAT, name
    assert (tmp_path / 'score_x3_flies.csv').exists()
    assert (tmp_path / 'e_model_distances.csv').exists()
    assert v['sensitivity']['left_out'] == [M.GROUP_L[0]] and v['sensitivity']['X1']['held'] == 'one function'
    assert v['sensitivity']['X3']['held'] == 'follows' and set(v['sensitivity']['planes']) == {'15 um'}
    assert v['V0']['n_flies'] == 0 and v['V0']['not_repeated'] == [] and v['V0']['not_decided_in_v0']
    assert 'pixel for pixel' in v['alignment']['note'] and v['alignment']['set_l_not_registered_together'] == []
    assert v['cross_check_rule'] is None
    # B1's cross-check failed: X1 and X2's claws and Q at b_G not decided, so not held; APL's X2 untouched; no 'follows'
    bl = json.loads((tmp_path / 'b_leak.json').read_text())
    (tmp_path / 'b_leak.json').write_text(json.dumps(dict(bl, cross_check_passed=False)))
    w = M.score(model_sets=_made_up_model_sets())
    assert w['X1']['held'] == 'not decided' and w['X1']['held_if_slope_carried'] == 'one function'
    assert all(x == 'not decided' for kk, x in w['X1']['by_reading'].items() if kk.endswith('b1'))
    assert w['X2']['held']['APL'] == v['X2']['held']['APL']
    for q in M.KC_QUANTITIES:
        assert w['X2']['held'][q] == 'not decided'
        assert w['X2']['held_if_slope_carried'][q] == v['X2']['held'][q]
    assert w['X3']['held'] != 'follows' and w['cross_check_rule']
    assert w['sensitivity']['X1']['held'] == 'not decided'


def test_a3_apl_tables_made_up(tmp_path, monkeypatch):
    d = tmp_path / 'A3'
    d.mkdir()
    stims = ['pa', 'ms', 'paa', 'paa [validation2/diag] @ -5', 'ms [megamat/diag] @ -3', M.BLANK_LABEL]
    win = {'pa': 21.0, 'ms': 19.0, 'paa': 20.0, 'paa [validation2/diag] @ -5': 17.0, 'ms [megamat/diag] @ -3': 18.0,
           M.BLANK_LABEL: 8.8}
    for ct in ('on', 'acute_off'):
        pd.DataFrame(dict(win=[win[s] for s in stims], pre=8.27, pre_settled=8.79),
                     index=pd.Index(stims, name='stimulus')).to_csv(d / f'real_{ct}_apl.csv')
        r = pd.DataFrame(np.zeros((10, len(stims))), columns=stims)
        r.iloc[0, 0] = 1.0
        for k in ('responses', 'spike_counts', 'pn_win', 'pn_pre'):
            r.to_parquet(d / f'real_{ct}_{k}.parquet')
    monkeypatch.setattr(M, 'a3_dir', lambda: d)
    monkeypatch.setattr(M, 'out_dir', lambda: tmp_path)
    res = M.a3_apl_tables(('real',))
    x = pd.read_csv(tmp_path / 'a3_apl_x2.csv')
    r1 = x[(x['ratio_id'] == 'r1') & (x['estimate'] == 'ms') & (x['a0'] == 'e90')].iloc[0]
    assert r1['apl_linear'] == pytest.approx((21.0 - 8.27) / (19.0 - 8.27))
    r1b = x[(x['ratio_id'] == 'r1') & (x['estimate'] == 'ms') & (x['a0'] == 'blank')].iloc[0]
    assert r1b['apl_linear'] == pytest.approx((21.0 - 8.8) / (19.0 - 8.8))
    r3 = x[(x['ratio_id'] == 'r3') & (x['a0'] == 'settled')].iloc[0]
    assert r3['apl_linear'] == pytest.approx((17.0 - 8.79) / (20.0 - 8.79)) and r3['kc_rate_num'] == 0.0
    assert set(x['ratio_id']) == {'r1', 'r2', 'r3'} and res['n_rows'] == 2 * len(stims)


# ---------------------------------------------------------------------------
# A2 again at the blank's level, C4's model side, the flies' strengths and the archive's variants

def test_at_blank_and_ratio_table():
    f = _fam_blank()
    f.loc[f['stimulus'] == M.BLANK_LABEL, 'pn_total'] = 5.0
    b = M.at_blank(f, 'on').set_index(['wiring', 'stimulus'])
    assert b.loc[('real', M.BLANK_LABEL), 'apl_blank'] == 0.0
    assert b.loc[('real', M.BLANK_LABEL), M.blank_col('input', 0.25)] == 0.0
    pa = f[(f.wiring == 'real') & (f.stimulus == 'pa') & (f.cond == 'on')].iloc[0]
    ms = f[(f.wiring == 'real') & (f.stimulus == 'ms') & (f.cond == 'on')].iloc[0]
    col = M.claw_col('input', 0.0, 'win')
    assert b.loc[('real', 'pa'), 'apl_blank'] == pytest.approx(pa.apl_win - 8.6)
    assert b.loc[('real', 'pa'), M.blank_col('input', 0.0)] == pytest.approx(pa[col] - 60.0)
    assert b.loc[('real', 'pa'), 'pn_blank'] == pytest.approx(pa.pn_total - 5.0)
    with pytest.raises(ValueError, match='blank'):
        M.at_blank(_fam(1.0), 'on')
    with pytest.raises(ValueError, match='level'):
        M.ratio_table(f, [dict(numerator='pa', denominator='ms')], 'settled')
    # APL acutely off: its own run's blank, not APL on's
    m = f['cond'] == 'off'
    f.loc[m, 'apl_win'] += 3.0
    f.loc[m & (f['stimulus'] == M.BLANK_LABEL), 'apl_win'] -= 2.0
    pairs = [dict(pair='pa/ms', numerator='pa', denominator='ms'), dict(pair='x/ms', numerator='x', denominator='ms')]
    t = pd.concat([M.ratio_table(f, pairs, 'blank', c) for c in ('on', 'off')])
    t = t[t['wiring'] == 'real']
    assert set(t['pair']) == {'pa/ms'}                                  # a pair the run lacks is skipped
    apl_on = (pa.apl_win - 8.6) / (ms.apl_win - 8.6)
    claws = (pa[col] - 60.0) / (ms[col] - 60.0)
    assert t[(t.cond == 'on') & (t.quantity == 'APL')]['linear'].iloc[0] == pytest.approx(apl_on)
    assert t[(t.cond == 'off') & (t.quantity == 'APL')]['linear'].iloc[0] == \
        pytest.approx((pa.apl_win + 3 - 9.6) / (ms.apl_win + 3 - 9.6))
    assert t[(t.cond == 'on') & (t.quantity == 'claws') & (t.subset == 'input') & (t.share == 0.0)]['linear'].iloc[0] == \
        pytest.approx(claws)
    q = t[(t.cond == 'on') & (t.quantity == 'Q')].iloc[0]
    assert q['linear'] == pytest.approx(apl_on / claws) and q['indicator_min'] <= q['linear'] <= q['indicator_max']
    pn = t[(t.cond == 'on') & (t.quantity == 'PN total')].iloc[0]
    assert pn['linear'] == pytest.approx((pa.pn_total - 5.0) / (ms.pn_total - 5.0))
    kc = t[(t.cond == 'on') & (t.quantity == 'KC response rate')].iloc[0]
    assert kc['linear'] == pytest.approx(pa.kc_rate / ms.kc_rate) and kc['kc_rate_num'] == pa.kc_rate
    e = M.ratio_table(f, pairs[:1], 'e90', 'on')
    assert e[e.quantity == 'APL']['linear'].iloc[0] == pytest.approx((pa.apl_win - 8.0) / (ms.apl_win - 8.0))
    assert e[e.quantity == 'PN total']['linear'].iloc[0] == pytest.approx(pa.pn_total / ms.pn_total)


def _fam_full(stims, wirings=('real', 'P1')):
    """A made-up A3 family table: per stimulus a claw drive c (0 for the blank); at the blank's level the claws read c and
    APL sqrt(c / 100) with APL on, c / 100 acutely off; E90's pre-odor means sit below the blank's windows."""
    rows = []
    for w in wirings:
        for i, s in enumerate(stims):
            c = 0.0 if s == M.BLANK_LABEL else 100.0 * (1 + i)
            for cond in ('on', 'off'):
                win = (8.6 + np.sqrt(c / 100.0)) if cond == 'on' else (9.0 + c / 100.0)
                rows.append(dict(wiring=w, cond=cond, stimulus=s, apl_win=win, apl_pre=8.0, apl_pre_settled=8.6,
                                 apl_evoked=win - 8.0, pn_total=10.0 + c, kc_rate=0.05 + c / 1e5,
                                 **{M.claw_col(sub, sh, k): (60.0 + c if k == 'win' else 50.0 if k == 'pre' else 10.0 + c)
                                    for sub in M.CLAW_SUBSETS for sh in M.CLAW_SHARES for k in ('win', 'pre', 'change')}))
    return pd.DataFrame(rows)


def test_part_a2_blank_and_c4_model_made_up(tmp_path, monkeypatch):
    import json
    hs = [M.hallem_label(o, c) for o in M.HALLEM_ODORS for c in M.HALLEM_CONCS]
    stims = list(M.H.APL16_ODORS) + hs + ['paa [validation2/diag] @ -5', 'ms [megamat/diag] @ -3', M.BLANK_LABEL]
    fam = _fam_full(stims)
    fam.to_csv(tmp_path / 'a3_family.csv', index=False)
    monkeypatch.setattr(M, 'out_dir', lambda: tmp_path)
    a2 = M.part_a2_blank(('real', 'P1'))
    ln = pd.read_csv(tmp_path / 'a2_blank_lines.csv')
    b = ln[(ln.level == 'blank') & (ln.subset == 'input') & (ln.share == 0.0)]
    assert np.allclose(b['slope_rma'], 0.5) and np.allclose(b['residual_log2'], 0.0, atol=1e-12)
    assert len(b) == 2 * 16 and set(ln['level']) == {'blank', 'e90'}
    e = ln[(ln.level == 'e90') & (ln.subset == 'input') & (ln.share == 0.0)]
    want = M.family_line(fam[fam['stimulus'].isin(list(M.H.APL16_ODORS) + [M.BLANK_LABEL])])
    assert np.allclose(e['slope_rma'].to_numpy(), want['slope_rma'].to_numpy())
    r = pd.read_csv(tmp_path / 'a2_blank_ratios.csv')
    assert set(r['pair']) == {'pa/ms', 'paa/ms', 'paa/pa'} and set(r['cond']) == {'on', 'off'}
    i_pa, i_ms = stims.index('pa'), stims.index('ms')
    apl = r[(r.baseline == 'blank') & (r.cond == 'on') & (r.pair == 'pa/ms') & (r.quantity == 'APL')]
    assert np.allclose(apl['linear'], np.sqrt((1 + i_pa) / (1 + i_ms)))
    apl_off = r[(r.baseline == 'blank') & (r.cond == 'off') & (r.pair == 'pa/ms') & (r.quantity == 'APL')]
    assert np.allclose(apl_off['linear'], (1 + i_pa) / (1 + i_ms))
    s = json.loads((tmp_path / 'a2_blank.json').read_text())
    assert any(x['level'] == 'blank' and x['subset'] == 'input' and x['share'] == 0.0
               and x['slope_rma_min'] == pytest.approx(0.5) for x in s['lines'])
    assert a2['lines_input_drive']
    c4 = M.part_c4_model(('real', 'P1'))
    h = pd.read_csv(tmp_path / 'c4_hallem_ratios.csv')
    i4p, i4m = stims.index(M.hallem_label('pa', -4)), stims.index(M.hallem_label('ms', -4))
    x = h[(h.series == 'pa over ms') & (h.log10_conc == -4) & (h.cond == 'on')]
    assert np.allclose(x[x.quantity == 'APL']['linear'], np.sqrt((1 + i4p) / (1 + i4m)))
    assert np.allclose(x[(x.quantity == 'claws') & (x.share == 0.0) & (x.subset == 'input')]['linear'], (1 + i4p) / (1 + i4m))
    i2p = stims.index(M.hallem_label('pa', -2))
    y = h[(h.series == 'pa over pa at 10^-2') & (h.log10_conc == -4) & (h.cond == 'off') & (h.quantity == 'APL')]
    assert np.allclose(y['linear'], (1 + i4p) / (1 + i2p))
    lv = pd.read_csv(tmp_path / 'c4_hallem_levels.csv')
    assert set(lv['stimulus']) == set(hs) | {'pa', 'ms'} and set(lv['cond']) == {'on', 'off'}
    off = pd.read_csv(tmp_path / 'c4_apl_off.csv')
    assert set(off['ratio_id']) == {'r1', 'r2', 'r3'} and set(off['cond']) == {'on', 'off'}
    r3 = off[(off.ratio_id == 'r3') & (off.cond == 'off') & (off.quantity == 'APL')]
    i5, ipaa = stims.index('paa [validation2/diag] @ -5'), stims.index('paa')
    assert np.allclose(r3['linear'], (1 + i5) / (1 + ipaa))
    assert c4['not_positive'] == [] and (tmp_path / 'c4_model.json').exists()


def test_flies_strengths():
    rows = []
    for fly, k in (('f1', 1.0), ('f2', 2.0)):
        for rec in (f'{fly}_r1', f'{fly}_r2'):
            for conc, a in ((-2.0, 0.8), (-4.0, 0.2)):
                for bm in (0.0, 1.0):
                    rows.append(dict(recording=rec, fly=fly, group='L', reading='carrying', b_mult=bm, stim=f'pa @ {conc:g}',
                                     odor='pa', log10_conc=conc, apl=a * k, kc=0.1 * k * (1 if conc == -2 else 0.5),
                                     kept=True))
    t = M.flies_strengths(pd.DataFrame(rows), groups=('L',), n_boot=200).set_index(['b_mult', 'stim'])
    assert t.loc[(0.0, 'pa @ -4'), 'apl_mean'] == pytest.approx((0.2 + 0.4) / 2)
    assert t.loc[(0.0, 'pa @ -4'), 'apl_over_2'] == pytest.approx(0.25)
    assert t.loc[(0.0, 'pa @ -4'), 'kc_over_2'] == pytest.approx(0.5) and t.loc[(0.0, 'pa @ -2'), 'apl_over_2'] == 1.0
    assert t.loc[(1.0, 'pa @ -4'), 'n_flies'] == 2 and t.loc[(1.0, 'pa @ -4'), 'n_recordings'] == 4


def test_write_variants(tmp_path, monkeypatch):
    import json
    from experiments import E83_cross_wiring_fit as E83
    monkeypatch.setattr(M, 'out_dir', lambda: tmp_path)
    monkeypatch.setattr(E83, 'members', lambda: {'real': ('E67_apl_effect_fits', 'fit_real'), 'P1': ('E72_x', 'fit_P1')})
    runs = pd.DataFrame([dict(wiring=w, cond=c, n_stimuli=52, input_sha256='abc', param_fixed_thr=100.0 + i,
                              param_wAPLKC_scale=0.1 if c == 'on' else 0.0, queue_commit='q', pyolfsysm_commit='p')
                         for i, w in enumerate(('real', 'P1')) for c in ('on', 'acute_off')])
    runs.to_csv(tmp_path / 'a3_runs.csv', index=False)
    _fam_full(list(M.H.APL16_ODORS) + [M.BLANK_LABEL]).to_csv(tmp_path / 'a3_family.csv', index=False)
    (tmp_path / 'a3_identity.json').write_text(json.dumps({'real': dict(reproduces_fit=True, reproduces_b1=True)}))
    v = M.write_variants(('real', 'P1')).set_index('tag')
    assert v.loc['real', 'fixed_thr'] == 100.0 and v.loc['P1', 'fixed_thr'] == 101.0
    assert v.loc['real', 'wAPLKC_scale'] == 0.1                       # the APL-on run's parameters
    assert v.loc['P1', 'ref_experiment'] == 'E72_x' and v.loc['real', 'ref_fit'] == 'fit_real'
    assert v.loc['real', 'reproduces_fit'] == True and pd.isna(v.loc['P1', 'reproduces_fit'])   # noqa: E712
    fam = _fam_full(list(M.H.APL16_ODORS))
    assert v.loc['real', 'rate_on'] == pytest.approx(fam[(fam.wiring == 'real') & (fam.cond == 'on')]['kc_rate'].mean())
    assert (tmp_path / 'variants.csv').exists()
