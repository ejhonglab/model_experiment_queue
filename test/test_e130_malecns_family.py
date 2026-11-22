#!/usr/bin/env python3
"""Tests for E130 (experiments/E130_malecns_family.py) and E10_joint_fit.py's --connectome (model_experiment_queue): the
connectome's kwargs and tag, the default unchanged (base_kws, fit_one's call, bouton_base's caches), main's refusals;
E130's KC-KC matrices, the factors on a base, the held counts, the fits' tags, commands and starts, and the classes and
the verdict on made-up readings. No simulations. Needs the shared directory (E10 reads the data-target registry at
import); skipped without it."""

import sys

import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
J = pytest.importorskip('experiments.E10_joint_fit')
E = pytest.importorskip('experiments.E130_malecns_family')
from experiments.E3_drive_dependent_gain import BOUTON  # noqa: E402


# ---------------------------------------------------------------------------
# E10_joint_fit.py's --connectome

def test_connectome_kws_and_tag():
    assert J.CONNECTOMES[0] == J.DEFAULT_CONNECTOME == 'hemibrain-first'
    assert J.connectome_kws() == {} and J.connectome_tag() == ''
    assert J.connectome_kws('hemibrain-argmax') == dict(multibouton_claws='argmax')
    for c in ('malecns-left', 'malecns-right'):
        assert J.connectome_kws(c) == dict(pn2kc_connections=c, multibouton_claws='argmax')
        assert J.connectome_tag(c) == f'cx-{c}_'
    with pytest.raises(ValueError):
        J.connectome_kws('flywire')
    with pytest.raises(ValueError):
        J.connectome_tag('hemibrain')


def test_base_kws_default_unchanged():
    """The default is E3's BOUTON (multibouton_claws 'first') exactly; the others change only what they select."""
    plain = dict(BOUTON, make_plots=False, silent=True)
    assert J.base_kws(None) == plain == J.base_kws(None, connectome='hemibrain-first')
    assert plain['multibouton_claws'] == 'first' and 'pn2kc_connections' not in plain
    a = J.base_kws(None, connectome='hemibrain-argmax')
    assert a == dict(plain, multibouton_claws='argmax')
    m = J.base_kws(None, connectome='malecns-left')
    assert m == dict(plain, multibouton_claws='argmax', pn2kc_connections='malecns-left')
    assert J.base_kws(40.0, 'spikes') == dict(plain, claw_thr=40.0, pn_claw_to_apl=False)


def _fake_base(n=12):
    idx = pd.Index([f'kc{i}' for i in range(n)], name=H.KC_ID_LEVEL)
    wnorm = {k: pd.Series(np.linspace(1, 2, n), index=idx) for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL')}
    breadth = pd.Series(np.arange(n, dtype=float), index=idx)
    initial = dict(fixed_thr=100.0, wAPLKC_scale=0.1, wKCAPL_scale=1e-4, wAPLPN_scale=0.1, wPNAPL_scale=1e-3)
    return initial, wnorm, breadth


@pytest.mark.parametrize('connectome', [None, 'hemibrain-first', 'hemibrain-argmax', 'malecns-right'])
def test_fit_one_passes_the_connectome(monkeypatch, connectome):
    import pyolfsysm.mb_model as mb
    seen = {}

    def fake(orn, **kw):
        seen.clear()
        seen.update(kw)
        return dict(history=pd.DataFrame(), params={}, loss=0.0)

    monkeypatch.setattr(mb, 'optimize_mb_model_params', fake)
    initial, wnorm, breadth = _fake_base()
    kw = dict(kckc0=60.0, maxfev=1, method='Nelder-Mead', thr0_factor=0.7)

    def loss_fn(o):
        return {}
    J.fit_one(2.0, None, initial, wnorm, breadth, loss_fn, **kw)
    plain = dict(seen)
    J.fit_one(2.0, None, initial, wnorm, breadth, loss_fn, **kw,
              **({} if connectome is None else dict(connectome=connectome)))
    extra = J.connectome_kws(connectome or 'hemibrain-first')
    assert plain['multibouton_claws'] == 'first' and 'pn2kc_connections' not in plain
    assert set(seen) == set(plain) | set(extra)
    for k, v in seen.items():
        if k in extra:
            assert v == extra[k]
        else:
            assert v is plain[k] or (isinstance(v, (pd.Series, pd.DataFrame)) and v.equals(plain[k])) or v == plain[k], k


def _fake_fit_mb_model(calls, n=6):
    idx = pd.MultiIndex.from_arrays([[f'kc{i}' for i in range(n)]], names=[H.KC_ID_LEVEL])

    def fake(orn, **kw):
        calls.append(kw)
        s = pd.DataFrame(np.arange(n * 2, dtype=float).reshape(n, 2), index=idx, columns=['a', 'b'])
        pd0 = {}
        for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL'):
            pd0[k] = pd.Series(np.ones(n), index=idx)
            pd0[f'{k}_scale'] = 1.0
        return (s >= 2).astype(float), s, pd.DataFrame(np.ones((n, 1)), index=idx, columns=['g']), pd0
    return fake


def test_bouton_base_connectome(monkeypatch, tmp_path):
    """The default: E3's base and hal's saved breadth (one call, no breadth run, no connectome kwarg). Another
    connectome: its base from --base-params (or connectome_base), its kwargs in both calls, the breadth run cached with
    the cx-<connectome>_ prefix (hal's saved run not used)."""
    import json
    import pyolfsysm.mb_model as mb
    calls = []
    monkeypatch.setattr(mb, 'fit_mb_model', _fake_fit_mb_model(calls))
    orn = pd.DataFrame(np.ones((1, 2)), index=['g'], columns=['a', 'b'])
    J.bouton_base(orn, tmp_path, None)
    assert len(calls) == 1 and calls[0]['multibouton_claws'] == 'first' and 'pn2kc_connections' not in calls[0]
    assert calls[0]['fixed_thr'] == J.e3_bouton_params()['fixed_thr']
    calls.clear()
    base = dict(fixed_thr=200.0, wAPLKC_scale=0.2, wKCAPL_scale=2e-4, wAPLPN_scale=0.2, wPNAPL_scale=5e-4)
    (tmp_path / 'b.json').write_text(json.dumps({'malecns-left': dict(base, source='x')}))
    init, _, _, breadth = J.bouton_base(orn, tmp_path, None, connectome='malecns-left',
                                        base_params=f'{tmp_path}/b.json:malecns-left')
    assert init == base and len(calls) == 2
    assert all(c['pn2kc_connections'] == 'malecns-left' and c['multibouton_claws'] == 'argmax' for c in calls)
    assert calls[0]['fixed_thr'] == 200.0 and calls[1]['wAPLKC'] == 0.0
    assert (tmp_path / 'cx-malecns-left_base_acute_off_spike_counts.parquet').exists()
    # connectome_base: E3's for the default; otherwise tuned once and cached
    calls.clear()
    assert J.connectome_base(orn, 'hemibrain-first', tmp_path) == J.e3_bouton_params()
    with pytest.raises(NotImplementedError):
        J.bouton_base(orn, tmp_path, 40.0, connectome='malecns-left')


def _main(monkeypatch, tmp_path, argv):
    orig_r, orig_o = H.results_dir, H.outputs_dir
    monkeypatch.setattr(H, 'results_dir', lambda e: (tmp_path / 'results') if e == 'E130_check' else orig_r(e))
    monkeypatch.setattr(H, 'outputs_dir', lambda e: (tmp_path / 'outputs') if e == 'E130_check' else orig_o(e))
    monkeypatch.setattr(sys, 'argv', ['E10_joint_fit.py', '--experiment', 'E130_check', *argv])
    J.main()


@pytest.mark.parametrize('argv, message', [
    (['--connectome', 'malecns-left'], 'needs --kckc-matrix'),
    (['--connectome', 'hemibrain-argmax', '--claw-thr', '40'], 'not supported'),
    (['--connectome', 'malecns-right', '--apl-drive', 'spikes', '--kckc-matrix', 'x.parquet', '--kckc-matrix-tag', 'x'],
     'not supported'),
    (['--connectome', 'hemibrain'], 'invalid choice'),
], ids=['male_without_matrix', 'claw_thr', 'spikes', 'unknown'])
def test_main_refuses(monkeypatch, tmp_path, capsys, argv, message):
    with pytest.raises(SystemExit) as e:
        _main(monkeypatch, tmp_path, argv)
    assert e.value.code == 2
    assert message in capsys.readouterr().err


# ---------------------------------------------------------------------------
# E130: the KC-KC matrices

def test_kckc_counts_and_normalization():
    t = pd.DataFrame(dict(bodyId_pre=[1, 1, 2, 3, 3, 9, 2], bodyId_post=[2, 2, 1, 3, 1, 1, 9],
                          roi_pre=['CA(L)'] * 6 + ['SLP(L)'], roi_post=['CA(L)'] * 7))
    M, info = E.kckc_counts(t, [3, 1, 2])
    # post x pre in the given order; the autapse 3->3 and the synapses with KC 9 dropped
    assert list(M.index) == [3, 1, 2] and list(M.columns) == [3, 1, 2]
    assert M.loc[2, 1] == 2 and M.loc[1, 2] == 1 and M.loc[1, 3] == 1 and M.to_numpy().sum() == 4
    assert info == dict(n_synapses=4, n_autapses_dropped=1, roi_pre={'CA(L)': 4})
    N = E.normalize_kckc(M)
    assert np.isclose(N.to_numpy().sum(axis=1).mean(), 1.0)
    assert np.allclose(N.to_numpy(), M.to_numpy() * 3 / 4)
    with pytest.raises(ValueError):
        E.normalize_kckc(M * 0)
    with pytest.raises(ValueError):
        E.kckc_counts(t, [1, 1, 2])


def test_malecns_kckc_uses_the_table(monkeypatch):
    t = pd.DataFrame(dict(bodyId_pre=[1, 2, 2], bodyId_post=[2, 1, 3], roi_pre=['CA(R)'] * 3, roi_post=['CA(R)'] * 3))
    monkeypatch.setitem(E._CACHE, 'kckc_table', t)
    M = E.malecns_kckc('malecns-right', [1, 2, 3])
    assert np.isclose(M.to_numpy().sum(), 3.0) and M.loc[2, 1] == 1.0 and M.loc[3, 2] == 1.0
    with pytest.raises(ValueError):
        E.malecns_kckc('hemibrain-argmax', [1, 2, 3])


def test_step0_kckc_and_apl():
    idx = pd.MultiIndex.from_arrays([[1, 1, 2, 3, 3, 3], [0, 1, 0, 0, 1, 2], ['g', 'g', 'ab', "a'b'", "a'b'", "a'b'"]],
                                    names=[H.KC_ID_LEVEL, 'claw_id', H.KC_TYPE_LEVEL])
    W = pd.DataFrame(np.eye(6, 2), index=idx, columns=pd.MultiIndex.from_arrays([['A', 'B'], [10, 11], [0, 0]],
                                                                                 names=['glomerulus', 'pn_id', 'bouton_id']))
    types = E.kc_types_of_wiring(W)
    assert types.to_dict() == {1: 'g', 2: 'ab', 3: "a'b'"}
    M = pd.DataFrame([[1, 2, 0], [0, 0, 3], [0, 0, 0]], index=[1, 2, 3], columns=[1, 2, 3], dtype=float)
    r = E.step0_kckc('x', M, dict(n_autapses_dropped=0), types)
    assert r['synapses'] == 6 and np.isclose(r['within_type_share'], 1 / 6)
    assert np.isclose(r['kcs_no_input'], 1 / 3) and r['kcs_no_output'] == 0.0 and r['kcs_in_table'] == 1.0
    assert np.isclose(r['incoming_g'], 3 / 2) and np.isclose(r['incoming_ab'], 3 / 2) and r['incoming_apbp'] == 0
    a = dict(apl_kc=pd.Series([0, 2, 3, 0, 1, 1.0]), kc_apl=pd.Series([1, 1, 0, 0, 0, 1.0]),
             apl_pn=pd.Series([1.0, 1.0]), pn_apl=pd.Series([2.0, 0.0]))
    s = E.step0_apl('x', W, a)
    assert np.isclose(s['claws_no_apl_kc'], 2 / 6) and np.isclose(s['bouton_share_apl_output'], 2 / 9)
    assert np.isclose(s['bouton_share_apl_input'], 2 / 5) and np.isclose(s['apl_kc_per_kc_g'], 2.0)
    wr = E.step0_wiring('x', W).set_index('kc_type')
    assert wr.loc['all', 'n_kcs'] == 3 and wr.loc["a'b'", 'claws_per_kc_mean'] == 3 and wr.loc['all', 'n_boutons'] == 2


# ---------------------------------------------------------------------------
# E130: step 1's point and the fits

def test_transfer_params():
    e3 = dict(fixed_thr=250.0, wAPLKC_scale=0.13, wKCAPL_scale=7e-5, wAPLPN_scale=0.13, wPNAPL_scale=3e-4)
    c2 = dict(fixed_thr=100.0, wAPLKC_scale=0.12, wKCAPL_scale=7.5e-5, wAPLPN_scale=0.37, wPNAPL_scale=9e-4, wKCKC_scale=26.4)
    p = E.transfer_params(e3, True, c2, e3)
    assert all(np.isclose(p[k], c2[k], rtol=1e-15) for k in c2)
    base = {k: 2 * v for k, v in e3.items()}
    q = E.transfer_params(base, False, c2, e3)
    assert all(np.isclose(q[k], 2 * c2[k]) for k in E.PARAM5) and q['wKCKC_scale'] == 0.0


def test_held_targets():
    assert E.held_targets(1732, 16, 1852, 4067, 1732 * 16) == (1852, 4067)
    on, off = E.held_targets(1864, 16, 1852, 4067, 1732 * 16)
    assert on == round(1852 * 1864 / 1732) and off == round(4067 * 1864 / 1732)


def test_names():
    assert E.calyx_of('malecns-left-kckc0') == 'malecns-left' and E.is_kckc0('malecns-left-kckc0')
    assert E.calyx_of('hemibrain-argmax') == 'hemibrain-argmax' and not E.is_kckc0('hemibrain-argmax')
    with pytest.raises(ValueError):
        E.calyx_of('flywire-left')
    assert set(E.FIT_TARGETS) == {t for ts in E.LANES.values() for t in ts}
    assert len(E.S1_RUNS) == 6


def test_fit_specs_tags_and_commands():
    specs = {s.name: s for s in E.fit_specs()}
    assert len(specs) == 10
    r = float(H.data_target('ctrl_frac'))
    s = specs['malecns-left']
    assert (s.start, s.suffix, s.kckc0) == ('factors', 'cap3', False)
    assert E.e10_tag(s) == f'cx-malecns-left_rate{r:.3f}_aplboth_kcm-malecns-left_k2_cap3'
    s0 = specs['malecns-right-kckc0']
    assert (s0.start, s0.kckc0) == ('held', True)
    assert E.e10_tag(s0) == f'cx-malecns-right_rate{r:.3f}_aplboth_kckc0_kcm-malecns-right_k2_cap3'
    assert specs['malecns-right-kckc0_s2'].start == 'factors' and specs['malecns-left_s2'].start == 'held'
    ha = specs['hemibrain-argmax']
    assert E.e10_tag(ha) == f'cx-hemibrain-argmax_rate{r:.3f}_aplboth_k2_cap3'
    # the tag is the one E10_joint_fit.py builds from the same options (its prefix rules)
    assert E.e10_tag(s0) == (J.connectome_tag('malecns-right') + f'rate{r:.3f}_' + 'aplboth_' + 'kckc0_'
                             + 'kcm-malecns-right_' + 'k2' + '_cap3')
    c = E.fit_command(s, python='py')
    assert c[:3] == ['py', '-u', 'experiments/E10_joint_fit.py']
    i = c.index('--connectome')
    assert c[i + 1] == 'malecns-left'
    assert c[c.index('--base-params') + 1] == 'results/E130_malecns_family/bases.json:malecns-left'
    assert c[c.index('--kckc-matrix') + 1] == 'outputs/E130_malecns_family/matrices/kckc_malecns-left.parquet'
    assert c[c.index('--init-from') + 1] == 'rate0.068_aplboth_k2_cap3' and '--kckc-fixed' not in c
    assert c[c.index('--experiment') + 1] == 'E130_malecns_family/fits/malecns-left'
    for flag, v in (('-k', '2'), ('--aplpn-max', '3'), ('--rate-target', 'ctrl_frac'), ('--loss', 'apl'),
                    ('--apl-readout', 'both'), ('--tag', '_cap3'), ('--maxfev', '150')):
        assert c[c.index(flag) + 1] == v, flag
    c0 = E.fit_command(s0, python='py')
    assert c0[c0.index('--kckc-fixed') + 1] == '0'
    assert c0[c0.index('--init-params') + 1] == 'results/E130_malecns_family/starts.json:S1_malecns-right-kckc0'
    ca = E.fit_command(ha, python='py')
    assert '--kckc-matrix' not in ca and '--kckc-fixed' not in ca


def test_second_start_and_counted(monkeypatch):
    rows = {}
    monkeypatch.setattr(E, 'fit_row', lambda spec: rows.get(spec.name))
    assert E.counted('malecns-left') is None
    rows['malecns-left'] = pd.Series({E.CHOOSE_COL: 0.060, 'loss': 4.0})
    c = E.counted('malecns-left')
    assert c['spec'].name == 'malecns-left' and c['n_starts'] == 1
    rows['malecns-left'] = pd.Series({E.CHOOSE_COL: 0.070, 'loss': 4.0})
    assert E.needs_second(rows['malecns-left']) and E.counted('malecns-left') is None
    rows['malecns-left_s2'] = pd.Series({E.CHOOSE_COL: 0.065, 'loss': 4.1})
    c = E.counted('malecns-left')
    assert c['spec'].name == 'malecns-left_s2' and c['n_starts'] == 2
    assert not E.needs_second(pd.Series({E.CHOOSE_COL: 0.066}))


def test_bound_flags():
    base = dict(fixed_thr=100.0, wAPLPN_scale=0.1)
    f = E.bound_flags(dict(fixed_thr=25.1, wAPLKC_scale=1.0, wAPLPN_scale=0.3, wKCKC_scale=299.0), dict(wAPLKC_scale=1.0), base)
    assert f['thr_at_bound'] and not f['aplkc_at_bound'] and f['aplpn_at_cap'] and f['kckc_at_bound']
    g = E.bound_flags(dict(fixed_thr=50.0, wAPLKC_scale=29.9, wAPLPN_scale=0.2, wKCKC_scale=20.0), dict(wAPLKC_scale=1.0), base)
    assert not g['thr_at_bound'] and g['aplkc_at_bound'] and not g['aplpn_at_cap'] and not g['kckc_at_bound']


# ---------------------------------------------------------------------------
# E130: the classes and the verdict on made-up readings

T = E.TARGETS_READ


def _same(v):
    return {t: v for t in T}


def test_class_new_default():
    fam = _same(0.056)
    assert E.class_new_default(_same(0.058), fam, True)['cls'] == 'unchanged'
    assert E.class_new_default(dict(_same(0.056), C=0.0595), fam, True)['cls'] == 'changed'
    assert E.class_new_default(_same(0.0531), fam, True)['cls'] == 'unchanged'
    assert E.class_new_default(_same(0.056), fam, False)['cls'] is None


def test_class_inside():
    rng = {t: (0.047, 0.064) for t in T}
    assert E.class_inside(_same(0.05), rng, True)['cls'] == 'inside'
    assert E.class_inside(dict(_same(0.05), B=0.07), rng, True)['cls'] == 'outside (worse side)'
    assert E.class_inside(dict(_same(0.05), A=0.04), rng, True)['cls'] == 'outside (better side)'
    assert E.class_inside(dict(_same(0.05), A=0.04, D=0.07), rng, True)['cls'] == 'outside (both sides)'
    assert E.class_inside(_same(0.05), rng, False)['cls'] is None


def test_class_population():
    tri = {t: list(np.linspace(0.050, 0.060, 21)) for t in T}
    assert E.class_population(_same(0.049), tri, True)['cls'] == 'real calyces average better'
    # below the lowest but not 0.003 below the median on one target
    t2 = dict(tri, A=list(np.linspace(0.0495, 0.0515, 21)))
    assert E.class_population(_same(0.049), t2, True)['cls'] == 'as the re-draws'
    assert E.class_population(_same(0.061), tri, True)['cls'] == 'worse'
    assert E.class_population(_same(0.055), tri, True)['cls'] == 'as the re-draws'
    assert E.class_population(_same(0.049), tri, False)['cls'] is None
    assert E.class_population(_same(0.049), tri, True)['targets']['all']['rank'] == 0


def test_class_kckc():
    free = {c: _same(0.055) for c in E.MALE}
    ok = {c: True for c in E.MALE}
    assert E.class_kckc(free, {c: _same(0.061) for c in E.MALE}, ok, True)['cls'] == 'replicated'
    assert E.class_kckc(free, {c: _same(0.057) for c in E.MALE}, ok, True)['cls'] == 'not replicated'
    assert E.class_kckc(free, {'malecns-left': _same(0.061), 'malecns-right': _same(0.056)}, ok, True)['cls'] == 'mixed'
    assert E.class_kckc(free, {c: _same(0.061) for c in E.MALE}, dict(ok, **{'malecns-left': False}), True)['cls'] == 'not testable'
    assert E.class_kckc(free, {c: _same(0.061) for c in E.MALE}, ok, False)['cls'] is None


def _made_up_readings():
    """Readings and runs of every unit the verdict reads (calibrated readout), with the family's eight between 0.047
    and 0.064."""
    fam = dict(zip(['real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3'], [0.056, 0.064, 0.052, 0.055, 0.047, 0.060, 0.061, 0.058]))
    rows, runs = [], []

    def add(unit, v, kind=None, u=None, in_bands=True):
        rows.append(dict(unit=unit, readout='cal', **{f'{t}_pattern': v for t in T}))
        if kind is not None:
            runs.append(dict(kind=kind, unit=u, in_bands=in_bands, realism_flag='' if in_bands else 'rate_on',
                             band_rate_on=in_bands, band_ratio=True, band_rate_rho_on=True, band_opk_on=True, band_opk_off=True))
    for w, v in fam.items():
        add(f'family:{w}', v, 'family', w)
        add(f'E110_B:{w}', v + 0.014, 'E110_B', w)
    add('family8', 0.048)
    add('E110_B8', 0.062)
    for t, v in {'hemibrain-argmax': 0.0575, 'malecns-left': 0.053, 'malecns-right': 0.066,
                 'malecns-left-kckc0': 0.060, 'malecns-right-kckc0': 0.073}.items():
        add(f'fit:{t}', v, 'fit', t)
    for i, (_, x, y) in enumerate(E.triples()):
        add(f'tri_{x}_{y}', 0.050 + 0.0005 * i)
    add('real3', 0.049)
    return pd.DataFrame(rows), pd.DataFrame(runs)


def test_verdict_on_made_up_readings():
    rd, runs = _made_up_readings()
    v = E.verdict(rd, runs, {c: True for c in E.MALE})
    assert v['new_default']['cls'] == 'unchanged'
    assert v['male']['malecns-left']['cls'] == 'inside'
    assert v['male']['malecns-right']['cls'] == 'outside (worse side)'
    assert v['population']['cls'] == 'real calyces average better'     # below 0.050, the lowest, and the median 0.055
    assert v['kckc']['cls'] == 'replicated'
    assert np.isclose(v['kckc_hemibrain_eight']['population']['all'], 0.014)
    assert v['realism_flags'] == [] and v['no_class_both_starts'] == []
    rd2 = rd.copy()
    rd2.loc[rd2['unit'] == 'real3', [f'{t}_pattern' for t in T]] = 0.0505
    assert E.verdict(rd2, runs, {c: True for c in E.MALE})['population']['cls'] == 'as the re-draws'
    # a run outside the bands: no class for it, and a realism flag
    runs.loc[(runs['kind'] == 'fit') & (runs['unit'] == 'malecns-left'), ['in_bands', 'band_rate_on']] = False
    v2 = E.verdict(rd, runs, {c: True for c in E.MALE})
    assert v2['male']['malecns-left']['cls'] is None and v2['kckc']['cls'] is None and v2['population']['cls'] is None
    assert [f['unit'] for f in v2['realism_flags']] == ['malecns-left']


def test_comparison_table_on_made_up_readings():
    rd, _ = _made_up_readings()
    md = E.comparison_table(rd, ['family:real', 'fit:malecns-left', 'no-such-unit'])
    lines = md.splitlines()
    assert lines[0] == '| statistic | family:real | fit:malecns-left |'
    assert '| pattern part, readout | 0.0560 | 0.0530 |' in lines
    assert all(len(line.split('|')) == 5 for line in lines)
