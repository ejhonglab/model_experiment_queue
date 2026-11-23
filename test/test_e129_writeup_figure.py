"""Tests of E129's write-up figure tables (analyses/E129_writeup_figure.py); a synthetic verdict, no shared-directory reads."""
import pytest

from analyses import E129_writeup_figure as F


def _verdict():
    def arm(v, cls):
        tg = {t: {'populations': {p: dict(arm=v + 0.001 * i + (0.002 if p == 'T4' else 0), family=0.05, E110_A=0.06,
                                         half_way=0.055) for p in F.POPULATIONS},
                  'd_family_per_wiring': {w: v - 0.05 + 0.0001 * k for k, w in
                                          enumerate(['real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3'])}}
              for i, t in enumerate(F.TARGETS)}
        return dict(cls=cls, targets=tg)
    return {'classes': {'Gpre': arm(0.07, 'neither'), 'G': arm(0.07, 'neither'), 'T': arm(0.057, 'neither'),
                        'P': arm(0.054, 'neither'), 'S': arm(0.050, 'keeps the gain')}}


def test_class_margins():
    M = F.class_margins(_verdict())
    assert len(M) == (len(F.ARMS) + 2) * len(F.TARGETS) * len(F.POPULATIONS)
    r = M[(M.kind == 'T') & (M.target == 'B') & (M.population == 'T4')].iloc[0]
    assert r.value == pytest.approx(0.057 + 0.002 + 0.002)  # target B is the third (0.001 per target) and r.margin == pytest.approx(r.value - 0.055)
    fam = M[M.kind == 'family']
    assert list(fam.margin) == pytest.approx([-0.005] * len(fam)) and (fam.cls == '').all()
    ref = M[M.kind == 'E110_A']
    assert list(ref.margin) == pytest.approx([0.005] * len(ref))
    assert set(M[M.kind == 'S'].cls) == {'keeps the gain'}


def test_wiring_diffs():
    D = F.wiring_diffs(_verdict())
    assert len(D) == len(F.ARMS) * 8
    assert set(D[D.population == 'S4'].wiring) == set(F.S_WIRINGS)
    r = D[(D.kind == 'S') & (D.wiring == 'M3')].iloc[0]
    assert r.population == 'T4' and r.d_family == pytest.approx(0.0007)


def test_draw(tmp_path):
    v = _verdict()
    p = tmp_path / 'f.png'
    F.draw(F.class_margins(v), F.wiring_diffs(v), p)
    assert p.stat().st_size > 1000
