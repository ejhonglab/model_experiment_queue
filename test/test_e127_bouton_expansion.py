"""Tests of E127 step 0's reader (no simulation, no saved data)."""
import numpy as np
import pandas as pd
import pytest

from experiments import E127_bouton_expansion as E


def test_expand_identity_and_reference():
    r = np.array([-3.0, 0.0, 5.0, 10.0, 40.0])
    assert np.allclose(E.expand(r, 10.0, 1.0), np.clip(r, 0, None))
    y = E.expand(r, 10.0, 2.0)
    assert y[3] == pytest.approx(10.0) and y[4] == pytest.approx(160.0) and y[2] == pytest.approx(2.5)
    # r_ref is a common factor: r_ref^(1 - q) r^q
    assert np.allclose(E.expand(r, 20.0, 2.5), (20.0 / 10.0) ** (1 - 2.5) * E.expand(r, 10.0, 2.5))


def test_expanded_boutons_off_and_on():
    pn_w, pn_p = np.array([[30.0], [10.0]]), np.array([[10.0], [10.0]])
    b_w, b_p = np.array([[25.0], [4.0]]), np.array([[9.0], [9.0]])
    w, p = E.expanded_boutons(pn_w, pn_p, b_w, b_p, 10.0, 1.0, 'on')
    assert np.allclose(w, b_w) and np.allclose(p, b_p)            # q 1 reproduces the saved boutons
    w, p = E.expanded_boutons(pn_w, pn_p, b_w, b_p, 10.0, 2.0, 'off')
    assert np.allclose(w, [[90.0], [10.0]]) and np.allclose(p, [[10.0], [10.0]])
    w, _ = E.expanded_boutons(pn_w, pn_p, b_w, b_p, 10.0, 2.0, 'on')
    assert np.allclose(w, [[85.0], [4.0]])                         # the held inhibition (5, 6) subtracted
    w, _ = E.expanded_boutons(pn_w, pn_p, b_w, b_p, 20.0, 2.0, 'on')
    assert np.allclose(w, [[40.0], [0.0]])                         # 45 - 5; 5 - 6 clipped at 0


def test_claws_changed():
    W = np.array([[1.0, 0.0], [0.0, 1.0]])
    c = np.array([[3.0], [1.0]])
    assert np.allclose(E.claws_changed(c, W, np.array([[2.0], [-5.0]])), [[5.0], [0.0]])


def _summ(rows):
    return pd.DataFrame([dict(criterion=2.0, reading='rate', cond='on', factor=1.0, q=q, totals_cv=t,
                              totals_cv_in=ti, resp_cv=0.5, resp_cv_in=ri, resp_cv_above=ra, top1=0.06, top1_in=pi,
                              top1_above=pa) for q, t, ti, ri, ra, pi, pa in rows])


def test_gate_rules():
    closed = _summ([(1.0, 0.1, 0, 0, 0, 0, 0), (2.0, 0.3, 0, 8, 0, 8, 0), (3.0, 0.4, 7, 0, 8, 0, 8)])
    v = E.gate(closed, 2.0)
    assert not v['opens'] and v['q'] is None and v['q_in']['resp_cv'] == [2.0, 2.0]
    opened = _summ([(1.0, 0.1, 0, 0, 0, 0, 0), (2.0, 0.4, 6, 8, 0, 3, 5), (2.5, 0.4, 8, 0, 2, 0, 1)])
    v = E.gate(opened, 2.0)
    assert v['opens'] and v['q'] == 2.0 and not v['opens_strict']
