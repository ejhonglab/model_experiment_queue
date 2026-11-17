#!/usr/bin/env python3
"""The natural-source alignment over all 120 pairs of Yang's 16 odors, not only the 66 pairs of the 12 megamat
odors (Yang's request on her version of the summary page, 2026-11-16; this script is the source of that answer's
numbers).

The page's natural-source alignment is the ρ (Spearman) of the residual of a KC matrix around its straight-line
prediction from the ORN matrix with the natural-source distances of the 66 megamat pairs (`kcq.harness.
reference_pairs()`, column vcf_hyp: the hyperbolic embedding of natural-source co-occurrence, D = 3, the mean of ten
seeds; the registry's ctrl_vcf / tnt_vcf). That table has the 12 megamat odors only. Anoop's embedding of the
natural-source statistics of the Hallem odors (D = 4; olfsysm's pyolfsysm/data/from_anoop/hallem_embeddings/
Hallem_VCF_D4_emb.csv, a square matrix whose header holds InChIs without "InChI=") covers all 16. Each genotype's KC
matrix (cosine distance per recording, mean over recordings) is read seven ways:

- kc_d4_120, res_d4_120: ρ of the KC distances / of their residual around the 120-pair ORN prediction with D4;
- kc_d4_66, res_d4_66: the same on the 66 megamat pairs (residual around the 66-pair ORN prediction);
- kc_d3_66, res_d3_66: the same with the page's embedding (res_d3_66 is the page's statistic; checked against the
  registry);
- kc_orn_120: ρ of the KC distances with the ORN distances over the 120 pairs;

each with a 90 % interval from 1000 resamples of the recordings (seed 5, control first), its mean over single
recordings and, for res_d4_120, kc_d4_120 and res_d3_66, a permutation test of the odor labels of the embedding
(10000 permutations; seed 7 for the 16 odors, 8 for the 12; one-sided p). The embeddings themselves are compared on
the 66 megamat pairs (D4 against the page's embedding and its raw distances, the page's ten seeds against each other)
and D4 with the ORN distances.

The odors' InChIs are olfdb's (below); (-)-alpha-pinene is in Anoop's matrix without its stereo layers and is matched
by its InChI cut before "/t" (recorded in provenance.json). Writes <shared>/results/yang_natural_source_120/:
alignment.csv, pairs.csv, embeddings.csv and provenance.json.
Usage: ./venv/bin/python analyses/yang_natural_source_120.py
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from kcq import harness as H

# olfdb's odors table, field `inchi` (http://hal:8080/api/odors, read 2026-11-16); not fetched at run time
ODOR_INCHI: Dict[str, str] = {
    '2h': 'InChI=1S/C7H14O/c1-3-4-5-6-7(2)8/h3-6H2,1-2H3',
    'IaA': 'InChI=1S/C7H14O2/c1-6(2)4-5-9-7(3)8/h6H,4-5H2,1-3H3',
    'pa': 'InChI=1S/C7H14O2/c1-3-4-5-6-9-7(2)8/h3-6H2,1-2H3',
    '2-but': 'InChI=1S/C4H8O/c1-3-4(2)5/h3H2,1-2H3',
    'eb': 'InChI=1S/C6H12O2/c1-3-5-6(7)8-4-2/h3-5H2,1-2H3',
    'ep': 'InChI=1S/C5H10O2/c1-3-5(6)7-4-2/h3-4H2,1-2H3',
    '6al': 'InChI=1S/C6H12O/c1-2-3-4-5-6-7/h6H,2-5H2,1H3',
    't2h': 'InChI=1S/C6H10O/c1-2-3-4-5-6-7/h4-6H,2-3H2,1H3/b5-4+',
    '1-8ol': 'InChI=1S/C8H18O/c1-2-3-4-5-6-7-8-9/h9H,2-8H2,1H3',
    '1-5ol': 'InChI=1S/C5H12O/c1-2-3-4-5-6/h6H,2-5H2,1H3',
    '1-6ol': 'InChI=1S/C6H14O/c1-2-3-4-5-6-7/h7H,2-6H2,1H3',
    'ms': 'InChI=1S/C8H8O3/c1-11-8(10)6-4-2-3-5-7(6)9/h2-5,9H,1H3',
    'paa': 'InChI=1S/C8H8O/c9-7-6-8-4-2-1-3-5-8/h1-5,7H,6H2',
    '-aPine': 'InChI=1S/C10H16/c1-7-4-5-8-6-9(7)10(8,2)3/h4,8-9H,5-6H2,1-3H3/t8-,9-/m0/s1',
    'euc': 'InChI=1S/C10H18O/c1-9(2)8-4-6-10(3,11-9)7-5-8/h8H,4-7H2,1-3H3',
    'B-myr': 'InChI=1S/C10H16/c1-5-10(4)8-6-7-9(2)3/h5,7H,1,4,6,8H2,2-3H3',
}
STATS = ('kc_d4_120', 'res_d4_120', 'kc_d4_66', 'res_d4_66', 'kc_d3_66', 'res_d3_66', 'kc_orn_120')
GENOTYPES = ('TNT3c', 'TNT')   # control, APL output blocked (TNT)
SEEDS = dict(bootstrap=5, permutations_16=7, permutations_12=8)
N_BOOT, N_PERM = 1000, 10000
PAGE_SEEDS = [f'vcf5q_ND030_D03_seed{k:02d}' for k in range(1, 11)]   # the page's embedding is their mean


def embedding_path() -> Path:
    """Anoop's D = 4 natural-source distances of the Hallem odors, in the installed pyolfsysm's data folder."""
    return H.pyolfsysm_data_dir() / 'from_anoop' / 'hallem_embeddings' / 'Hallem_VCF_D4_emb.csv'


def tidy_path() -> Path:
    """The lab's megamat distances in other spaces (al_analysis, as apl_pair_analysis/apl_analysis.py reads it)."""
    return H.al_analysis_dir() / 'scripts' / 'corr_minus_orn_vs_vcf_repro' / 'data' / 'tidy_abbrev_dists.tsv'


def load_embedding(path: Path) -> pd.DataFrame:
    """A square distance matrix whose header holds the InChIs (without "InChI=") and whose rows are in its order."""
    emb = pd.read_csv(path)
    assert emb.shape[0] == emb.shape[1], f'{path} is not square'
    emb.index = emb.columns
    return emb


def odor_keys(columns: Sequence[str], inchis: Dict[str, str] = ODOR_INCHI) -> Dict[str, Tuple[str, str]]:
    """odor -> (its key in the embedding's header, 'exact' or 'stereo layers dropped'), the key being the InChI
    without "InChI=", or failing that the InChI cut before "/t".

    Raises:
        KeyError: an odor that is in neither form.
    """
    cols, out = set(columns), {}
    for o, inchi in inchis.items():
        key = inchi[len('InChI='):] if inchi.startswith('InChI=') else inchi
        if key in cols:
            out[o] = (key, 'exact')
        elif key.split('/t')[0] in cols:
            out[o] = (key.split('/t')[0], 'stereo layers dropped')
        else:
            raise KeyError(f'{o} ({inchi}) is not in the embedding')
    return out


def odor_matrix(emb: pd.DataFrame, keys: Dict[str, Tuple[str, str]], odors: Sequence[str]) -> np.ndarray:
    k = [keys[o][0] for o in odors]
    return emb.loc[k, k].to_numpy(dtype=float)


def upper(m: np.ndarray) -> np.ndarray:
    """A square matrix's upper triangle in the order of itertools.combinations of its rows."""
    return m[np.triu_indices(len(m), 1)]


def lookup(table, pairs: Sequence[Tuple[str, str]]) -> np.ndarray:
    """The values of a pair-indexed Series for `pairs`, each pair looked up in either orientation (NaN if in
    neither)."""
    out = []
    for a, b in pairs:
        if (a, b) in table.index:
            out.append(table.loc[(a, b)])
        elif (b, a) in table.index:
            out.append(table.loc[(b, a)])
        else:
            out.append(np.nan)
    return np.asarray(out, dtype=float)


def resid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """`y` minus its least-squares straight line on `x` (the residual around the linear ORN prediction)."""
    b, a = np.polyfit(x, y, 1)
    return y - (a + b * x)


def rho(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman's ρ; NaN when either side is constant (a residual of zeros, say)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if np.ptp(a) <= 1e-12 * max(1.0, np.abs(a).max()) or np.ptp(b) <= 1e-12 * max(1.0, np.abs(b).max()):
        return float('nan')
    return float(spearmanr(a, b)[0])


def make_stats(d4: np.ndarray, orn: np.ndarray, is66: np.ndarray, orn66: np.ndarray, hyp66: np.ndarray
        ) -> Callable[[np.ndarray], Dict[str, float]]:
    """The seven statistics of a vector of 120 KC distances (pairs in combinations order of H.APL16_ODORS)."""
    def stats(kc: np.ndarray) -> Dict[str, float]:
        return dict(kc_d4_120=rho(kc, d4), res_d4_120=rho(resid(kc, orn), d4),
            kc_d4_66=rho(kc[is66], d4[is66]), res_d4_66=rho(resid(kc[is66], orn66), d4[is66]),
            kc_d3_66=rho(kc[is66], hyp66), res_d3_66=rho(resid(kc[is66], orn66), hyp66),
            kc_orn_120=rho(kc, orn))
    return stats


def bootstrap(X: Dict[str, np.ndarray], stats, n: int = N_BOOT, seed: int = SEEDS['bootstrap']) -> Dict[str, pd.DataFrame]:
    """genotype -> resamples x statistics, the recordings of each genotype drawn with replacement from one
    generator, the genotypes in the order of `X`."""
    rng = np.random.default_rng(seed)
    out = {}
    for g, x in X.items():
        out[g] = pd.DataFrame([stats(x[rng.integers(0, len(x), len(x))].mean(0)) for _ in range(n)])
    return out


def permutations(n_odors: int, n: int, seed: int) -> np.ndarray:
    """n random permutations of the odors (n x n_odors), from `np.random.default_rng(seed)`."""
    rng = np.random.default_rng(seed)
    return np.array([rng.permutation(n_odors) for _ in range(n)])


def perm_test(fixed: np.ndarray, Dm: np.ndarray, perms: np.ndarray) -> Tuple[float, float, float, float]:
    """Spearman ρ of `fixed` (a vector over the upper triangle) with the matrix Dm, against the null of Dm with its
    odors permuted (Dm[np.ix_(p, p)] for each p of `perms`).

    Returns:
        (observed ρ, one-sided p = (number of null values >= observed + 1) / (n + 1), the null's SD, its 95th
        percentile).
    """
    iu = np.triu_indices(len(Dm), 1)
    observed = rho(fixed, upper(Dm))
    # Spearman = Pearson of ranks; a permutation only moves the triangle's values, so their ranks move with them
    rmat = np.zeros_like(Dm, dtype=float)
    rmat[iu] = rankdata(upper(Dm))
    rmat = rmat + rmat.T
    null_ranks = rmat[perms[:, iu[0]], perms[:, iu[1]]]
    f = rankdata(fixed)
    f = (f - f.mean()) / np.linalg.norm(f - f.mean())
    nr = null_ranks - null_ranks.mean(1, keepdims=True)
    null = (nr @ f) / np.linalg.norm(nr, axis=1)
    p = (np.sum(null >= observed - 1e-12) + 1) / (len(null) + 1)
    return observed, float(p), float(null.std()), float(np.percentile(null, 95))


def seed_agreement(tidy: pd.DataFrame, pairs: Sequence[Tuple[str, str]], cols: Sequence[str] = PAGE_SEEDS) -> float:
    """The mean Spearman ρ between every two of the page's embedding seeds over `pairs`."""
    vals = {c: lookup(tidy[c], pairs) for c in cols}
    return float(np.mean([rho(vals[a], vals[b]) for a, b in combinations(cols, 2)]))


def _sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main() -> None:
    odors = list(H.APL16_ODORS)
    pairs = list(combinations(odors, 2))
    rec = {g: H.yang_recording_rdms(g, 'cosine') for g in GENOTYPES}   # recordings x 120 pairs, cosine distance
    for g in GENOTYPES:
        assert list(rec[g].columns) == pairs, f'{g}: the pairs are not in combinations order'
    X = {g: rec[g].to_numpy(dtype=float) for g in GENOTYPES}
    emb = load_embedding(embedding_path())
    keys = odor_keys(emb.columns)
    assert set(keys) == set(odors), 'an odor of the panel has no InChI'
    d4 = np.array([emb.loc[keys[a][0], keys[b][0]] for a, b in pairs])
    D16 = odor_matrix(emb, keys, odors)
    assert np.allclose(D16, D16.T) and np.allclose(upper(D16), d4)
    orn = H.orn16_consensus_pairs().reindex(pd.MultiIndex.from_tuples(pairs)).to_numpy(dtype=float)
    assert not np.isnan(orn).any(), 'an ORN pair is missing'
    ref = H.reference_pairs()
    is66 = np.array([a in H.MEGAMAT12 and b in H.MEGAMAT12 for a, b in pairs])
    p66 = [p for p, m in zip(pairs, is66) if m]
    orn66, hyp66, raw66 = lookup(ref['orn'], p66), lookup(ref['vcf_hyp'], p66), lookup(ref['vcf_raw'], p66)
    assert not np.isnan(orn66).any() and not np.isnan(hyp66).any()
    stats = make_stats(d4, orn, is66, orn66, hyp66)
    full = {g: stats(X[g].mean(0)) for g in GENOTYPES}
    reg = json.loads((H.shared_dir() / 'results' / 'data_targets' / 'data_targets.json').read_text())
    for g, name in (('TNT3c', 'ctrl_vcf'), ('TNT', 'tnt_vcf')):
        assert abs(full[g]['res_d3_66'] - reg[name]['value']) < 1e-6, f'res_d3_66 of {g} is not the registry\'s {name}'
    boot = bootstrap(X, stats)
    per_rec = {g: pd.DataFrame([stats(x) for x in X[g]]).mean() for g in GENOTYPES}
    perms16 = permutations(16, N_PERM, SEEDS['permutations_16'])
    D12 = np.zeros((12, 12))
    D12[np.triu_indices(12, 1)] = hyp66   # the 66 megamat pairs are the combinations of H.MEGAMAT12 in order
    D12 = D12 + D12.T
    assert [o for o in odors if o in H.MEGAMAT12] == list(H.MEGAMAT12)
    perms12 = permutations(12, N_PERM, SEEDS['permutations_12'])
    rows = []
    for g in GENOTYPES:
        kc = X[g].mean(0)
        tests = {'res_d4_120': perm_test(resid(kc, orn), D16, perms16), 'kc_d4_120': perm_test(kc, D16, perms16),
            'res_d3_66': perm_test(resid(kc[is66], orn66), D12, perms12)}
        lo, hi = boot[g].quantile(0.05), boot[g].quantile(0.95)
        for s in STATS:
            t = tests.get(s)
            if t is not None:
                assert abs(t[0] - full[g][s]) < 1e-9
            rows.append(dict(genotype=g, statistic=s, n_pairs=120 if s.endswith('120') else 66,
                embedding='hallem_d4' if '_d4_' in s else 'page_d3' if '_d3_' in s else '',
                on='residual' if s.startswith('res_') else 'kc', value=full[g][s],
                lo=np.percentile(boot[g][s], 5), hi=np.percentile(boot[g][s], 95), per_recording_mean=per_rec[g][s],
                perm_p=t[1] if t else None, null_sd=t[2] if t else None, null_q95=t[3] if t else None))
    alignment = pd.DataFrame(rows)
    kc_c, kc_t = X['TNT3c'].mean(0), X['TNT'].mean(0)
    pairs_tab = pd.DataFrame(dict(odor_a=[a for a, _ in pairs], odor_b=[b for _, b in pairs], megamat66=is66, d4=d4,
        orn16=orn, kc_ctrl=kc_c, kc_tnt=kc_t, resid_ctrl=resid(kc_c, orn), resid_tnt=resid(kc_t, orn)))
    tidy = pd.read_csv(tidy_path(), sep='\t', index_col=[0, 1])
    has_raw = ~np.isnan(raw66)
    d466 = d4[is66]
    embeddings = pd.DataFrame([dict(statistic=k, value=v) for k, v in dict(
        d4_vs_page_embedding_66=rho(d466, hyp66),
        d4_vs_page_raw_66=rho(d466[has_raw], raw66[has_raw]),
        page_embedding_vs_raw_66=rho(hyp66[has_raw], raw66[has_raw]),
        page_seed_to_seed_mean_66=seed_agreement(tidy, p66),
        orn_vs_d4_120=rho(orn, d4),
        n_megamat_pairs_without_raw=int((~has_raw).sum())).items()])
    out = H.shared_dir() / 'results' / 'yang_natural_source_120'
    out.mkdir(parents=True, exist_ok=True)
    alignment.to_csv(out / 'alignment.csv', index=False)
    pairs_tab.to_csv(out / 'pairs.csv', index=False)
    embeddings.to_csv(out / 'embeddings.csv', index=False)
    inputs = dict(hallem_d4_embedding=embedding_path(), page_reference_pairs=H.apl_pair_dir() / 'shared12_pairs.csv',
        page_embedding_seeds=tidy_path(), per_recording_rdms=H.yang_data_dir() / 'computed' / 'rdm_per_recording_long.csv',
        orn16=H.shared_dir() / 'data' / 'model_outputs' / 'inputs' / 'orn16_corrdist_olfdb_36glom_flymean.csv',
        data_targets=H.shared_dir() / 'results' / 'data_targets' / 'data_targets.json')
    prov = dict(script='analyses/yang_natural_source_120.py', **H.code_commits(),
        inputs={k: dict(path=str(p), sha256=_sha256(p)) for k, p in inputs.items()},
        odor_inchi_source='olfdb odors table, field inchi (http://hal:8080/api/odors), read 2026-11-16',
        odor_keys={o: dict(inchi=ODOR_INCHI[o], embedding_key=k, match=m) for o, (k, m) in keys.items()},
        stereo_layers_dropped=[o for o, (_, m) in keys.items() if m != 'exact'],
        seeds=SEEDS, n_bootstrap=N_BOOT, n_permutations=N_PERM, page_embedding_seed_columns=PAGE_SEEDS,
        generated=datetime.now(timezone.utc).isoformat(timespec='seconds'),
        host=subprocess.run(['hostname'], capture_output=True, text=True).stdout.strip())
    (out / 'provenance.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    pd.set_option('display.width', 220)
    print(alignment.round(4).to_string(index=False))
    print(embeddings.round(4).to_string(index=False))
    print(f'written to {out}')


if __name__ == '__main__':
    main()
