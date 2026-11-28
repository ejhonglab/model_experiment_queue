"""Tests of E145's write-up tables (analyses/E145_writeup_tables.py) on planted score files; no shared-directory reads
and no simulation."""
import json

import numpy as np
import pandas as pd

from analyses import E145_writeup_tables as T


def _cmp(change, spread=0.05, worse=False, ref=1.0, se=0.01):
    return dict(complete=True, change=change, spread=spread, worse_than_spread=worse, ref=ref, arm=ref + change, se=se)


def _cls(cls, **kw):
    r = dict(cls=cls, megamat=_cmp(-0.1, 0.02), stages={s: _cmp(-0.05) for s in T.STAGES}, stage_worse=[],
             guard_worse=[], panel_worse=[], reasons=[], guard={s: _cmp(0.0) for s in T.STAGES},
             panels={k: _cmp(0.01) for k in T.PANEL_KEYS}, panels_beside={k: _cmp(0.2, 0.1) for k in T.BESIDE_KEYS})
    r.update(kw)
    return r


def _fit(z, rises, at_bound=()):
    det = {k: dict(rise=v, determined=v >= 0.01) for k, v in rises.items()}
    cz = {k: (0.0 if not det[k]['determined'] else z[k]) for k in z}
    return dict(n_evals=126, best_eval=116, best_z=z, at_bound=list(at_bound), determined=det, classed_z=cz,
                classed_loss=0.912)


def _d():
    z = dict(aplpn=float(np.log(0.427)), kckc=float(np.log(2.53)), inhadd=float(np.log(52.0)), inhsc=float(np.log(0.71)))
    fit = _fit(z, dict(aplpn=0.007, kckc=0.0139, inhadd=0.766, inhsc=0.0516))
    hist = pd.DataFrame(dict(eval=[0, 116], loss_sim=[1.0396, 0.9052]))
    cons = pd.DataFrame([dict(measure='strength', statistic='stage_error', n_evals=126, rho=0.802, yardstick=0.8)])
    rep = dict(complete=True, aplpn=dict(value=1.0, determined=False, repeats=False),
               kckc=dict(value=2.53, determined=True, repeats=True), inhadd=dict(value=1633.0, determined=True,
                                                                                 repeats=True),
               inhsc_gain=dict(value=0.71, determined=True, read=False))
    rec_f = dict(cell='F', variant='E145_F', cls=_cls("one stage at another's cost", stage_worse=['kc_claw'],
                                                     panel_worse=['validation2/validation2/on']),
                 price=dict(complete=True, price=True, price_value=0.0137), repeats=rep, fit={})
    rec_x = dict(cell='X', variant='E145_X', cls=_cls('no class', reasons=['at a bound: kckc high']),
                 price=dict(complete=True, price=True, price_value=0.0381), repeats=rep, fit={})
    ex = dict(complete=True, cls=dict(X='no class', F="one stage at another's cost", W='holds'), fall_repeats=True,
              megamat_fixes=dict(inhadd=False, aplpn=False, kckc=True), hold_repeats=False, the_flies=False,
              the_wirings=False, no_class_cells=['X'], n_cells_hold=1, outcome=['d', 'f'])
    empty = dict(fit={}, hist=pd.DataFrame(), det=pd.DataFrame(), cons=pd.DataFrame(), rescore=pd.DataFrame(),
                 reads=pd.DataFrame())
    cells = {c: dict(empty) for c in T.CELLS}
    cells['F'] = dict(empty, fit=fit, hist=hist, cons=cons)
    shares = pd.DataFrame([dict(cell='F', share_scored=0.588), dict(cell='X', share_scored=1.0)])
    e143 = dict(derived=dict(fall_repeats=dict(share_scored=0.931), outcome_d=dict(share_scored=0.995)),
                shortcut_check=dict(tol=0.02))
    shortcut = pd.DataFrame([dict(cell='F', statistic='stage', diff=0.0427, within_tol=False),
                             dict(cell='F', statistic='pn_dendrite', diff=-0.01, within_tol=True)])
    return dict(classes=dict(cells=dict(X=rec_x, F=rec_f), experiment=ex, control_fall=0.04), cells=cells,
                shares=shares, e143=e143, shortcut=shortcut, prediction={}, step0=pd.DataFrame(), mega=pd.DataFrame())


def test_helpers():
    assert T._rng(pd.Series([0.16, 0.1601])) == '0.160'
    assert T._rng(pd.Series([0.1, 0.2]), '.2f') == '0.10 to 0.20'
    assert T._names(['apl16/apl16/on', 'kc_claw']) == 'APL16, APL on; KC claws'
    assert T._names([]) == 'none'
    assert np.isclose(T.factor_value('inhadd', float(np.log(52.0))), 52.0 * 31.4088, rtol=1e-4)
    assert np.isclose(T.factor_value('kckc', float(np.log(2.53))), 2.53)
    assert T._fmt_factor('kckc', 2.53) == '×2.530' and T._fmt_factor('inhadd', 1632.6) == '1633'


def test_fits_and_returns_tables():
    d = _d()
    V = T.W.Values()
    rows = [r for r in T.fits_table(d, V) if r.startswith('| F,')]
    assert len(rows) == 1
    assert '| 126 | 1.040 | 0.905 (116) | ×0.427 | ×2.530 | 1633 | ×0.710 | none | 0.802 |' in rows[0]
    ret = [r for r in T.returns_table(d, V) if r.startswith('| F,')]
    assert len(ret) == 4
    assert '| APL→bouton | +0.0070 | no | ×1.000 | no |' in ret[0]          # undetermined: back to the start
    assert '| KC–KC | +0.0139 | yes | ×2.530 | yes |' in ret[1]
    assert ret[3].endswith('| not read |')
    v = V.frame().set_index(['arm', 'statistic'])['value']
    assert v[('F', 'consistency_rho')] == 0.802 and v[('F', 'classed_aplpn')] == 1.0 and v[('F', 'repeats_kckc')] == 1


def test_classes_and_experiment_tables():
    d = _d()
    V = T.W.Values()
    rows = {r.split(' |')[0]: r for r in T.classes_table(d, V) if r.startswith('| X,') or r.startswith('| F,')}
    x = next(r for k, r in rows.items() if k.startswith('| X'))
    f = next(r for k, r in rows.items() if k.startswith('| F'))
    assert 'no class (at a bound: kckc high)' in x and 'yes (+0.0381)' in x and x.endswith('| 1.000 |')
    assert "one stage at another's cost" in f and '-0.100 (0.020)' in f and '| KC claws |' in f
    assert '| validation2 |' in f and f.endswith('| 0.588 |')
    ex = T.experiment_table(d, V)
    assert any(r.startswith('| the fall repeats') and r.endswith('| yes | 0.931 |') for r in ex)
    assert any(r.startswith('| megamat fixes KC–KC | yes |') for r in ex)
    assert any(r.startswith('| (d) the cells disagree | yes | 0.995 |') for r in ex)
    assert any(r.startswith('| (a) ') and '| no |' in r for r in ex)
    v = V.frame().set_index(['arm', 'statistic'])['value']
    assert v[('experiment', 'd')] == 1.0 and v[('experiment', 'a')] == 0.0 and v[('F', 'price_value')] == 0.0137


def test_shortcut_table():
    d = _d()
    V = T.W.Values()
    rows = [r for r in T.shortcut_table(d, V) if r.startswith('| F,')]
    assert rows == ['| F, the other flies: S, half b | 0.0427 | 0.0427 | 1 of 2 |']


def test_kckc_table_and_class_without_bound():
    import numpy as np
    hi = float(np.log(7.3))
    rng = np.random.default_rng(0)
    n = 40
    z = dict(aplpn=rng.normal(-1.7, 0.1, n), kckc=np.r_[np.full(3, hi), rng.uniform(0.0, 1.9, n - 3)],
             inhadd=rng.normal(4.2, 0.05, n), inhsc=rng.normal(-1.0, 0.05, n))
    loss = 0.90 - 0.008 * z['kckc'] + rng.normal(0, 0.001, n)
    loss[0] = 1.03
    hist = pd.DataFrame(dict(eval=np.arange(n), loss_sim=loss, **{f'z_{k}': v for k, v in z.items()}))
    hist.loc[0, 'z_kckc'] = 0.0
    resc = pd.concat([pd.DataFrame(dict(eval=np.arange(n), measure='strength', half=h, stage_error=loss,
                                        kc_claw_error=c + 0.05 * hist.z_kckc)) for h, c in (('fit', 0.97), ('held', 0.79))])
    d = _d()
    d['cells']['X'] = dict(d['cells']['X'], hist=hist, rescore=resc, fit=dict(best_eval=int(np.argmin(loss))))
    d['classes']['cells']['X']['cls'].update(stage_worse=['kc_claw'], improves=False)
    d['classes']['cells']['X']['cls']['megamat']['improves'] = True
    V = T.W.Values()
    rows = T.kckc_table(d, V)
    assert any('Evaluations at the KC–KC bound (×7.3): 2 of 40' in r for r in rows)   # eval 0 moved off it
    sl = T.local_slope(hist, 0.05)
    assert sl['n'] == n - 1 and abs(sl['slope'] + 0.008) < 0.002                       # all but the start
    assert any(r.startswith('| at ×7.3 (the bound) | 2 |') for r in rows)
    assert rows[-2] == "- X's class with the bound rule left out: one stage at another's cost."
    assert T.kckc_bin(7.3, 7.3).startswith('at ×7.3') and T.kckc_bin(6.0, 7.3) == '×5.83 to ×7.3'
    assert T.kckc_bin(1.0, 7.3) == 'below ×1.8' and T.kckc_bin(1.8, 7.3) == '×1.8 to ×3.24'
    rec = dict(cls=dict(megamat=dict(improves=True), reasons=['at a bound: kckc high'], improves=True, panel_worse=[],
                        panel_missing=[]))
    assert T.class_without_bound(rec) == 'holds'
    rec['cls']['reasons'].append('leaves a realism band on P3')
    assert T.class_without_bound(rec) == 'no class'
