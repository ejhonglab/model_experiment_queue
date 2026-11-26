"""Tests of experiments/E139_normalization_ladder.py (E139, the weaker antennal-lobe normalization without the calyx
refit, along its strength): the points' coordinates in E136's A1 and the planned reads. No model is run."""
import numpy as np
import pytest

from experiments import E136_stage_population_fit as F
from experiments import E139_normalization_ladder as L


def test_point_z():
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    z = L.point_z('B', 1000.0)
    assert z['aplpn'] == pytest.approx(np.log(0.155)) and z['kckc'] == 0.0 and z['inhsc'] == 0.0
    assert PN_INHADD_DEFAULT * np.exp(z['inhadd']) == pytest.approx(1000.0)
    kw = F.pn_norm_kws(z)
    from experiments.E77_input_screen import inhsc_for
    assert kw['pn_inhadd'] == pytest.approx(1000.0) and kw['pn_inhsc'] == pytest.approx(inhsc_for(1000.0))
    assert L.point_z('K', 314.0)['kckc'] == pytest.approx(np.log(3.24))


def test_planned_points():
    p1, p2 = L.planned_points(['1']), L.planned_points(['2'])
    assert len(p1) == 96 and len(p2) == 96
    assert {(p.setting, p.rung) for p in p2} == {('B', 314.0), ('B', 1000.0), ('K', 1000.0)}
    assert {p.panel for p in p1} == {'megamat', 'validation2', 'mixtures'} and 'apl16' in {p.panel for p in p2}
    assert len({L.point_tag(p) for p in p1 + p2}) == 192
    assert L.point_tag(L.Point('F', 314.0, 'megamat', 'P3')) == 'F_inhadd314_megamat_P3'
