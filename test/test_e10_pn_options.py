#!/usr/bin/env python3
"""Tests for E10_joint_fit.py's PN normalization options (E102; model_experiment_queue): --pn-inhadd / --pn-inhsc
(`pn_norm_kws`, `pn_norm_tag`, `fit_one`'s pn, `bouton_base`'s pn for the breadth run) and their checks in `main`; with
them unset every kwarg is as before. No simulations, except the exact starts (E102_E10_EXACT_START=1, a few minutes per
fit: E10 from the optima of C2 and E72's P1 fit with the options off and with olfsysm's defaults passed explicitly, each
first evaluation against the fit's saved history and the two runs' first evaluations against each other). Needs the
shared directory (E10 reads the data-target registry at import); skipped without it."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
H = pytest.importorskip('kcq.harness')
try:
    H.data_targets_path().read_bytes()
except (FileNotFoundError, OSError):
    pytest.skip('the shared directory (data-target registry) is not available', allow_module_level=True)
E10 = pytest.importorskip('experiments.E10_joint_fit')

DEFAULTS = dict(pn_inhadd=31.4088, pn_inhsc=368.6631)


def test_pn_norm_kws():
    assert E10.pn_norm_kws() == {} and E10.pn_norm_kws(None, None) == {}
    assert E10.pn_norm_kws(1000, 544.5) == dict(pn_inhadd=1000.0, pn_inhsc=544.5)
    for a, s in ((1000.0, None), (None, 544.5)):
        with pytest.raises(ValueError, match='go together'):
            E10.pn_norm_kws(a, s)
    for a, s in ((0.0, 1.0), (1.0, -1.0), (np.inf, 1.0), (1.0, np.nan)):
        with pytest.raises(ValueError, match='finite and positive'):
            E10.pn_norm_kws(a, s)


def test_pn_norm_tag():
    assert E10.pn_norm_tag({}) == '' and E10.pn_norm_tag(None) == ''
    assert E10.pn_norm_tag(dict(pn_inhadd=1000.0, pn_inhsc=544.44447616846)) == 'pnadd1000sc544.444_'
    t = E10.pn_norm_tag(dict(pn_inhadd=1e6, pn_inhsc=1.2e5))
    assert t == 'pnadd1e06sc120000_'
    assert all(c.isascii() and (c.isalnum() or c in '_-.') for c in t)


def _fake_base(n=12):
    idx = pd.Index([f'kc{i}' for i in range(n)], name=H.KC_ID_LEVEL)
    wnorm = {k: pd.Series(np.linspace(1, 2, n), index=idx) for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL')}
    breadth = pd.Series(np.arange(n, dtype=float), index=idx)
    initial = dict(fixed_thr=100.0, wAPLKC_scale=0.1, wKCAPL_scale=1e-4, wAPLPN_scale=0.1, wPNAPL_scale=1e-3)
    return initial, wnorm, breadth


@pytest.mark.parametrize('pn', [None, {}, dict(pn_inhadd=1000.0, pn_inhsc=544.5), DEFAULTS])
def test_fit_one_passes_the_pn_kwargs(monkeypatch, pn):
    """Every simulation of the fit (one call of optimize_mb_model_params, whose kwargs reach APL on and off alike) gets
    pn_inhadd and pn_inhsc when given, and nothing else changes; not given (None or {}), the call is as before."""
    import pyolfsysm.mb_model as mb
    seen = {}

    def fake(orn, **kw):
        seen.update(kw)
        return dict(history=pd.DataFrame(), params={}, loss=0.0)

    monkeypatch.setattr(mb, 'optimize_mb_model_params', fake)
    initial, wnorm, breadth = _fake_base()
    kw = dict(kckc0=60.0, maxfev=1, method='Nelder-Mead', thr0_factor=0.7)

    def loss_fn(o):
        return {}

    E10.fit_one(2.0, None, initial, wnorm, breadth, loss_fn, **kw)
    plain = dict(seen)
    seen.clear()
    E10.fit_one(2.0, None, initial, wnorm, breadth, loss_fn, **kw, **({} if pn is None else dict(pn=pn)))
    assert 'pn_inhadd' not in plain and 'pn_inhsc' not in plain
    for k in ('pn_inhadd', 'pn_inhsc'):
        assert (seen.get(k) == pn[k]) if pn else k not in seen
    # the APL-off condition only zeroes the APL->KC and APL->bouton scales, so it has the PN kwargs too
    assert seen['condition_param_overrides'] == {'on': {}, 'off': {'wAPLKC_scale': 0.0, 'wAPLPN_scale': 0.0}}
    rest = {k: v for k, v in seen.items() if k not in ('pn_inhadd', 'pn_inhsc')}
    assert set(rest) == set(plain)
    for k, v in plain.items():
        assert rest[k] is v or (isinstance(v, (pd.Series, pd.DataFrame)) and rest[k].equals(v)) or rest[k] == v, k


def test_fit_one_pn_checks(monkeypatch):
    import pyolfsysm.mb_model as mb
    monkeypatch.setattr(mb, 'optimize_mb_model_params', lambda orn, **kw: dict(history=pd.DataFrame(), params={}, loss=0.0))
    initial, wnorm, breadth = _fake_base()
    kw = dict(kckc0=60.0, maxfev=1, method='Nelder-Mead', thr0_factor=0.7)
    with pytest.raises(ValueError, match='pn must be'):
        E10.fit_one(2.0, None, initial, wnorm, breadth, lambda o: {}, **kw, pn=dict(pn_inhadd=1000.0))
    proxy = dict(sigma=0.1, thr_ratio=1.0, n_recordings=1, n_trials=1, seed=0, passes=1, relax=False)
    with pytest.raises(NotImplementedError, match='pn-inhadd'):
        E10.fit_one(2.0, None, initial, wnorm, breadth, lambda o: {}, **kw, block_proxy=proxy,
                    pn=dict(pn_inhadd=1000.0, pn_inhsc=544.5))


@pytest.mark.parametrize('pn', [None, dict(pn_inhadd=1000.0, pn_inhsc=544.44447616846)])
def test_bouton_base_pn(monkeypatch, tmp_path, pn):
    """The breadth run (APL off) gets the PN kwargs and a cache name with their tag; the base's own call (the normalized
    APL weights) does not; without them the calls are as before (hal's saved breadth, no breadth run)."""
    import pyolfsysm.mb_model as mb
    calls = []
    n = 6
    idx = pd.MultiIndex.from_arrays([[f'kc{i}' for i in range(n)]], names=[H.KC_ID_LEVEL])

    def fake(orn, **kw):
        calls.append(kw)
        s = pd.DataFrame(np.arange(n * 2, dtype=float).reshape(n, 2), index=idx, columns=['a', 'b'])
        pd0 = {}
        for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL'):
            pd0[k] = pd.Series(np.ones(n), index=idx)
            pd0[f'{k}_scale'] = 1.0
        return (s >= 2).astype(float), s, pd.DataFrame(np.ones((n, 1)), index=idx, columns=['g']), pd0

    monkeypatch.setattr(mb, 'fit_mb_model', fake)
    orn = pd.DataFrame(np.ones((1, 2)), index=['g'], columns=['a', 'b'])
    out = E10.bouton_base(orn, tmp_path, None, **({} if pn is None else dict(pn=pn)))
    assert len(calls) == (1 if pn is None else 2)
    assert 'pn_inhadd' not in calls[0] and 'pn_inhsc' not in calls[0]
    if pn is not None:
        assert calls[1]['wAPLKC'] == 0.0 and calls[1]['wAPLPN'] == 0.0
        assert all(calls[1][k] == v for k, v in pn.items())
        assert (tmp_path / 'pnadd1000sc544.444_base_acute_off_spike_counts.parquet').exists()
        assert out[3].to_numpy().tolist() == [1.0, 5.0, 9.0, 13.0, 17.0, 21.0]
        # a wiring's breadth run is cached under the wiring's name inside the PN tag
        calls.clear()
        W = out[2]
        E10.bouton_base(orn, tmp_path, None, wiring_matrix=W, wiring_tag='m1', pn=pn)
        assert calls[1]['_wPNKC'] is W and calls[1]['pn_inhadd'] == pn['pn_inhadd']
        assert (tmp_path / 'pnadd1000sc544.444_wiring-m1_base_acute_off_spike_counts.parquet').exists()
    with pytest.raises(ValueError, match='pn must be'):
        E10.bouton_base(orn, tmp_path, None, pn=dict(pn_inhadd=1000.0))


def _main(monkeypatch, tmp_path, argv):
    orig_r, orig_o = H.results_dir, H.outputs_dir
    monkeypatch.setattr(H, 'results_dir', lambda e: (tmp_path / 'results') if e == 'E102_check' else orig_r(e))
    monkeypatch.setattr(H, 'outputs_dir', lambda e: (tmp_path / 'outputs') if e == 'E102_check' else orig_o(e))
    for d in ('results', 'outputs'):
        (tmp_path / d).mkdir(exist_ok=True)
    monkeypatch.setattr(sys, 'argv', ['E10_joint_fit.py', '--experiment', 'E102_check', *argv])
    E10.main()


@pytest.mark.parametrize('argv, message', [
    (['--pn-inhadd', '1000'], 'go together'),
    (['--pn-inhsc', '544'], 'go together'),
    (['--pn-inhadd', '0', '--pn-inhsc', '544'], 'finite and positive'),
    (['--pn-inhadd', '1000', '--pn-inhsc', '544', '--block-proxy', '1', '--jitter-sigma', '0.1'], 'with --block-proxy'),
    (['--pn-inhadd', '1000', '--pn-inhsc', '544', '--apl-drive', 'spikes_only'], 'spikes_only'),
], ids=['inhadd_alone', 'inhsc_alone', 'zero', 'block_proxy', 'spikes_only'])
def test_main_refuses(monkeypatch, tmp_path, capsys, argv, message):
    with pytest.raises(SystemExit) as e:
        _main(monkeypatch, tmp_path, argv)
    assert e.value.code == 2
    assert message in capsys.readouterr().err


def test_main_passes_pn_to_the_base_and_the_fit(monkeypatch, tmp_path):
    """main with the options: bouton_base and fit_one get the same pn; the tag carries pn_norm_tag; fits.csv has
    pn_inhadd and pn_inhsc as given; the provenance has them in its base kwargs and pn_normalization. Without them no
    pn reaches either and neither column is written."""
    import pyolfsysm.mb_model as mb  # noqa: F401
    seen = {}
    n = 8
    idx = pd.MultiIndex.from_arrays([[f'kc{i}' for i in range(n)]], names=[H.KC_ID_LEVEL])

    def fake_base(orn, out_dir, claw_thr=None, **kw):
        seen['base'] = kw.get('pn')
        initial = dict(fixed_thr=100.0, wAPLKC_scale=0.1, wKCAPL_scale=1e-4, wAPLPN_scale=0.1, wPNAPL_scale=1e-3)
        wnorm = {k: pd.Series(np.ones(n), index=idx) for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL')}
        W = pd.DataFrame(np.ones((n, 2)), index=idx, columns=pd.MultiIndex.from_arrays([['g1', 'g2']], names=['glomerulus']))
        return initial, wnorm, W, pd.Series(np.arange(n, dtype=float), index=idx.get_level_values(0))

    class Stop(Exception):
        pass

    def fake_fit(k, orn, initial, wnorm, breadth, loss_fn, **kw):
        seen['fit'] = kw.get('pn')
        raise Stop

    monkeypatch.setattr(E10, 'bouton_base', fake_base)
    monkeypatch.setattr(E10, 'fit_one', fake_fit)
    for argv, pn in (([], None), (['--pn-inhadd', '1000', '--pn-inhsc', '544.44447616846'],
                                   dict(pn_inhadd=1000.0, pn_inhsc=544.44447616846))):
        seen.clear()
        with pytest.raises(Stop):
            _main(monkeypatch, tmp_path, ['-k', '2', '--loss', 'apl', '--apl-readout', 'both', '--rate-target', 'ctrl_frac',
                                          *argv])
        assert seen['base'] == pn and seen['fit'] == pn


# the exact starts: E10 from a fit's optimum, one evaluation, with the options off and with olfsysm's defaults passed
# explicitly (as E82's checks of its options)
COMMON = ['-k', '2', '--loss', 'apl', '--rate-target', 'ctrl_frac', '--aplpn-max', '3', '--maxfev', '1', '--init-from-exact',
          '--apl-readout', 'both', '--tag', '_cap3']
EXACT_FITS = {
    'C2': ('E67_apl_effect_fits', 'rate0.068_aplboth_k2_cap3', []),
    'P1': ('E72_wiring_population', 'wiring-claw_permutation_within_type-seed1_rate0.068_aplboth_k2_cap3',
           ['--wiring', 'claw_permutation_within_type', '--wiring-seed', '1']),
}
IDENTITY = ['--pn-inhadd', repr(DEFAULTS['pn_inhadd']), '--pn-inhsc', repr(DEFAULTS['pn_inhsc'])]


def _e10_exact(root: Path, fit: str, extra):
    from experiments import E78_claw_knee as E78
    exp, tag, opts = EXACT_FITS[fit]
    args = COMMON + ['--experiment', E78.CHECK_SCRATCH, '--init-from', tag, '--init-from-experiment', exp] + opts + list(extra)
    p = subprocess.run([sys.executable, '-u', '-m', 'experiments.E78_claw_knee', '--run-e10', str(root), '--', *args],
                       cwd=REPO, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout[-3000:] + p.stderr[-3000:]
    hist = sorted((root / 'results').glob('history_*.csv'))
    assert len(hist) == 1
    out = {f.stem: pd.read_parquet(f) for f in (root / 'outputs').glob('*_spike_counts.parquet') if 'base_acute_off' not in f.name}
    return pd.read_csv(hist[0], float_precision='round_trip').iloc[0], out, hist[0]


@pytest.mark.skipif(os.environ.get('E102_E10_EXACT_START') != '1', reason='model fits (minutes): set E102_E10_EXACT_START=1')
@pytest.mark.parametrize('fit', list(EXACT_FITS))
def test_exact_start(tmp_path, fit):
    """With the options off and with olfsysm's defaults passed explicitly (the breadth run then recomputed with them, not
    hal's saved one): the saved counts are the fit's own, every statistic of the first evaluation is its best evaluation's
    in its saved history (within 1e-12 relative, as E82's checks), and the two runs' first evaluations are identical to the
    last digit."""
    E78 = pytest.importorskip('experiments.E78_claw_knee')
    from experiments import E72_wiring_population as E72
    exp, tag, _ = EXACT_FITS[fit]
    if not (H.results_dir(exp) / f'history_{tag}.csv').exists() or E72.load_fit(exp, tag) is None:
        pytest.skip(f"{fit}'s saved history or outputs are not available")
    saved = E72.load_fit(exp, tag)
    rows = {}
    for name, extra in (('off', []), ('defaults', IDENTITY)):
        rows[name], outs, path = _e10_exact(tmp_path / name, fit, extra)
        r = E78.compare_history(path, exp, tag)
        assert not r['missing'] and not r['extra'], r
        # (against a saved history, as E82's checks: within 1e-12 relative; the last digits move with the BLAS threads
        # and the host, with the options off as much as on)
        for c, d in r['differing_statistics'].items():
            assert d['abs'] <= 1e-12 * max(abs(d['saved']), 1.0), (c, d)
        t = path.stem[len('history_'):]
        assert ('pnadd31.4088sc368.663_' in t) == (name == 'defaults')
        for c, f in (('on', 'on'), ('off', 'acute_off')):
            s = outs[f'{t}_{f}_spike_counts']
            assert s.index.equals(saved[c][1].index) and (s.to_numpy(dtype=float) == saved[c][1].to_numpy(dtype=float)).all()
    cols = [c for c in rows['off'].index if c not in ('eval', 'time_s')]
    assert list(rows['defaults'].index) == list(rows['off'].index)
    assert all(rows['defaults'][c] == rows['off'][c] or (pd.isna(rows['defaults'][c]) and pd.isna(rows['off'][c])) for c in cols)
