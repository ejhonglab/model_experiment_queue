"""The registry's like-for-like level entries (kcq/data_targets.py; Yang's question of 2026-11-16, unresolved_issues.md):
the ORN cosine distances come from the same fly-mean vectors whose Pearson distances are the ORN reference, and the
cosine of made-up vectors is computed as 1 - cosine similarity over the 120 pairs."""
import numpy as np
import pandas as pd
import pytest

H = pytest.importorskip('kcq.harness')
DT = pytest.importorskip('kcq.data_targets')

ODORS = list(H.APL16_ODORS)


def test_orn_cosine_pairs_on_made_up_vectors(tmp_path):
    rng = np.random.default_rng(0)
    V = pd.DataFrame(rng.normal(size=(36, 16)), index=[f'g{i}' for i in range(36)], columns=ODORS)
    p = tmp_path / 'vectors.csv'
    V[ODORS[::-1]].to_csv(p)
    d = DT.orn_cosine_pairs(p)
    assert len(d) == 120 and list(d.index[:2]) == [(ODORS[0], ODORS[1]), (ODORS[0], ODORS[2])]
    a, b = V[ODORS[3]].to_numpy(), V[ODORS[9]].to_numpy()
    assert d[(ODORS[3], ODORS[9])] == pytest.approx(1 - a @ b / np.linalg.norm(a) / np.linalg.norm(b), rel=1e-12)
    assert DT.orn_cosine_pairs(p)[(ODORS[0], ODORS[1])] == pytest.approx(d[(ODORS[0], ODORS[1])])


def test_orn_vectors_are_the_reference_matrix_vectors():
    path = DT.INPUTS['orn_vectors']
    if not path.exists():
        pytest.skip(f'{path} not on this host')
    V = pd.read_csv(path, index_col=0)[ODORS]
    assert V.shape == (36, 16) and not V.isna().any().any()
    ref = H.orn16_consensus_pairs()
    for (a, b), r in ref.items():
        assert 1 - np.corrcoef(V[a], V[b])[0, 1] == pytest.approx(r, abs=1e-12)
    cos = DT.orn_cosine_pairs()
    assert 0 < cos.mean() < ref.mean()
