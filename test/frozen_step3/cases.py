"""The frozen cases on the queue's side: every function that moved to hong2p or pyolfsysm (2026-11-12), called through
the queue's names (kcq.harness, experiments.E69 / E77 / E72) with the queue's defaults, on the frozen inputs
(`frozen_outputs.objects`). hong2p's and olfsysm's copies of these cases (their test_frozen_from_queue.py) call the
moved functions directly, with the queue's parameters, and must give the same outputs.

Each case is name -> (function of the inputs, solver): `solver` marks outputs of root-finding or fitting, compared to a
looser tolerance only when the library versions differ from the frozen files' (`frozen_outputs.check`).
"""
from typing import Callable, Dict, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

Cases = Dict[str, Tuple[Callable[[dict], object], bool]]
RATE = 59.273 / 20


def _pairs(O):
    return H.pair_distances(O['counts_on'])


def _model(O):
    return H.pair_distances(O['counts_on']), H.pair_distances(O['counts_off'])


def hong2p_cases() -> Cases:
    """hong2p.rdm, hong2p.stats, hong2p.imaging and hong2p.viz."""
    def boot_fn(m):
        return {'mean_on': float(m['TNT3c'].mean()), 'rho': H.spearman(m['TNT3c'], m['TNT'])}

    def reversed_pairs(O):
        p = H.pairs_from_square(O['sq_ctrl'])
        rev = pd.Series(p.to_numpy()[::-1], index=pd.MultiIndex.from_tuples([(b, a) for a, b in p.index[::-1]]))
        return rev

    def with_nan(O):
        p = H.pairs_from_square(O['sq_ctrl'])
        a = p.copy()
        a.iloc[[3, 50]] = np.nan
        b = H.subset_pairs(p, H.MEGAMAT12 + ['paa']) * 0.5
        return a, b

    c: Cases = {
        'odor_order': (lambda O: H.odor_order(['B-myr', 'zz', '2h', 'aa', 'pa', 'euc']), False),
        'pair_list': (lambda O: H.pair_list(['ms', '2h', 'zz', 'IaA', 'aa']), False),
        'pairs_from_square': (lambda O: H.pairs_from_square(O['sq_ctrl']), False),
        'square_from_pairs': (lambda O: H.square_from_pairs(H.pairs_from_square(O['sq_ctrl'])), False),
        'square_from_pairs_subset': (lambda O: H.square_from_pairs(H.subset_pairs(H.pairs_from_square(O['sq_ctrl']),
            ['eb', '2h', 'paa', 'ms'])), False),
        'canonical_pairs': (lambda O: H.canonical_pairs(reversed_pairs(O)), False),
        'align_pairs': (lambda O: H.align_pairs(*with_nan(O)), False),
        'subset_pairs': (lambda O: H.subset_pairs(H.pairs_from_square(O['sq_ctrl']), H.MEGAMAT12), False),
        'pair_distances_cosine': (lambda O: H.pair_distances(O['counts_on']), False),
        'pair_distances_cosine_off': (lambda O: H.pair_distances(O['counts_off']), False),
        'pair_distances_pearson': (lambda O: H.pair_distances(O['orn'], measure='pearson'), False),
        'pair_distances_pearson_nosilent': (lambda O: H.pair_distances(O['counts_on'], 'pearson', drop_silent=True), False),
        'pair_distances_tm': (lambda O: H.pair_distances(O['tm_mat']), False),
        'block_label': (lambda O: [H.block_label(a, b) for a, b in (('2h', 'IaA'), ('2h', 'eb'), ('ms', 'paa'),
            ('euc', '-aPine'), ('zz', 'ms'))], False),
        'block_means': (lambda O: H.block_means(H.pairs_from_square(O['sq_ctrl'])), False),
        'per_odor_means': (lambda O: H.per_odor_means(H.pairs_from_square(O['sq_ctrl'])), False),
        'pair_design': (lambda O: H.pair_design(H.pairs_from_square(O['sq_ctrl']).index), False),
        'pair_design_odors': (lambda O: H.pair_design(H.subset_pairs(_pairs(O), ['2h', 'eb', 'ms']).index,
            ['ms', 'eb', '2h', 'paa']), False),
        'odor_shares': (lambda O: H.odor_shares(O['data_effect']), False),
        'density_null_similarity': (lambda O: H.density_null_similarity(O['tm_mat']), False),
        'shuffled_pair_distances': (lambda O: H.shuffled_pair_distances(O['tm_mat'], np.random.default_rng(0)), False),
        'shuffled_pair_distances_thr': (lambda O: H.shuffled_pair_distances(O['tm_mat'], np.random.default_rng(1),
            measure='pearson', values='thresholded', cells='responsive', thresh=0.8), False),
        'transform_matrix': (lambda O: [H.transform_matrix(O['tm_mat'].iloc[:60], values=v, cells=cl)
            for v in ('graded', 'binary', 'thresholded') for cl in ('all', 'responsive')], False),
        'transform_matrix_thresh': (lambda O: H.transform_matrix(O['tm_mat'], values='binary', thresh=0.75), False),
        'linear_resid': (lambda O: H.linear_resid(*_model(O)), False),
        'isotonic_resid': (lambda O: H.isotonic_resid(*_model(O)), False),
        'spearman_pearson': (lambda O: [H.spearman(*_model(O)), H.pearson(*_model(O)),
            H.spearman(*with_nan(O)), H.pearson(*with_nan(O))], False),
        'resample_recording_means': (lambda O: H.resample_recording_means({'TNT3c': O['rec_ctrl'], 'TNT': O['rec_tnt']},
            np.random.default_rng(1)), False),
        'resample_recording_means_fly': (lambda O: H.resample_recording_means({'TNT3c': O['rec_ctrl'],
            'TNT': O['rec_tnt']}, np.random.default_rng(2), unit='fly'), False),
        'bootstrap_recordings': (lambda O: H.bootstrap_recordings(boot_fn, {'TNT3c': O['rec_ctrl'], 'TNT': O['rec_tnt']},
            n_boot=30), False),
        'bootstrap_recordings_fly': (lambda O: H.bootstrap_recordings(boot_fn, {'TNT3c': O['rec_ctrl'],
            'TNT': O['rec_tnt']}, n_boot=25, seed=3, ci=(10, 90), unit='fly'), False),
        'split_half_reliability': (lambda O: H.split_half_reliability(O['rec_ctrl'], n_splits=40), False),
        'split_half_reliability_pearson': (lambda O: H.split_half_reliability(O['rec_tnt'], n_splits=30, seed=4,
            method='pearson'), False),
        'oklab_diverging_cmap': (lambda O: H.oklab_diverging_cmap(), False),
        'oklab_diverging_cmap_args': (lambda O: H.oklab_diverging_cmap(11, l_end=0.5, c_end=0.2, hue_neg=200.0,
            hue_pos=40.0), False),
        'window_volumes': (lambda O: [H.window_volumes(w, RATE) for w in ((2, 8), (2.0, 8.1), (0, 2), (-2, -1))], False),
        'trace_excerpts': (lambda O: H.trace_excerpts(O['traces'], O['onsets'], 15, 60), False),
        'windowed_responses': (lambda O: H.windowed_responses(H.trace_excerpts(O['traces'], O['onsets'], 15, 60),
            [(6, 24), (0, 6), (-15, 0)], first_volume=-15), False),
        'windowed_responses_seconds': (lambda O: H.windowed_responses(H.trace_excerpts(O['traces'], O['onsets'], 15, 60),
            [(2.0, 8.1), (0.5, 3.0)], first_volume=-15, volume_rate=RATE), False),
        'indicator_kernel': (lambda O: [H.indicator_kernel(np.linspace(-1, 30, 311), r, d)
            for r, d in ((1.3, 20.0), (0.0, np.inf), (0.0, 5.0), (0.7, np.inf))], False),
        'indicator_kernel_integral': (lambda O: [H.indicator_kernel_integral(np.linspace(-1, 30, 311), r, d)
            for r, d in ((1.3, 20.0), (0.0, np.inf), (0.0, 5.0), (0.7, np.inf))], False),
        'kernel_signal': (lambda O: [H.kernel_signal(np.arange(-15, 119) / RATE, 0.8, 10.0, 'onset'),
            H.kernel_signal(np.arange(-15, 119) / RATE, 0.8, 10.0, 'steady', pulse=2.0),
            H.kernel_signal(np.arange(-15, 119) / RATE, 0.8, 10.0, 'histogram', edges=[0.0, 1.0, 2.5],
                weights=[2.0, 1.0])], False),
        'fit_indicator_kernel': (lambda O: H.fit_indicator_kernel(O['synth_t'], O['synth_y'], 'onset'), True),
        'fit_indicator_kernel_steady': (lambda O: H.fit_indicator_kernel(O['synth_t'], O['synth_y'], 'steady',
            pulse=2.0), True),
        'fit_indicator_kernel_fixed': (lambda O: H.fit_indicator_kernel(O['synth_t'], O['synth_y'], 'onset',
            fixed=(1.0, 8.0)), False),
        'trial_values': (lambda O: [H.trial_values(O['trials_ctrl'][:, :40], v) for v in ('graded', 'thresholded',
            'binary')] + [H.trial_values(O['trials_ctrl'][:, :40], 'thresholded', thresh=0.5)], False),
        'trial_cosines': (lambda O: H.trial_cosines(O['trials_ctrl']), False),
        'noise_corrected_pair_distances': (lambda O: H.noise_corrected_pair_distances(O['trials_ctrl'], O['odors']),
            False),
        'noise_corrected_pair_distances_rel': (lambda O: H.noise_corrected_pair_distances(
            H.trial_values(O['trials_tnt'], 'thresholded'), O['odors'], min_reliability=0.3), False),
        'noise_free_pair_distances': (lambda O: [H.noise_free_pair_distances(O['trials_ctrl'], O['odors'], repeats=r,
            values=v) for r in ('different', 'adjacent', 'first_third') for v in ('graded', 'thresholded')], False),
        'noise_free_pair_distances_binary': (lambda O: H.noise_free_pair_distances(O['trials_tnt'], O['odors'],
            values='binary', thresh=0.5), False),
    }
    return c


def queue_cases() -> Cases:
    """Wrappers whose moved core needs the queue's data (only the queue's test runs these)."""
    return {
        'kc_trial_noise_sd': (lambda O: [H.kc_trial_noise_sd('TNT3c'), H.kc_trial_noise_sd('TNT'),
            H.kc_noise_floor()], False),
        'code_commits_keys': (lambda O: sorted(k for k in H.code_commits() if k.startswith(('queue', 'pyolfsysm',
            'al_analysis'))), False),
    }


def pyolfsysm_cases() -> Cases:
    """pyolfsysm.readout, pyolfsysm.kc_stats, pyolfsysm.effects and pyolfsysm.tuning."""
    from experiments import E69_kc_firing_ceiling as E69
    from experiments import E72_wiring_population as E72
    import frozen_step3.frozen_outputs as F

    def lut(O):
        return H.rank_match_map(O['counts_on'], O['ref_amp'])

    def z(O, **kw):
        kw = dict(dict(q=0.128, noise=2.0), **kw)
        return H.session_zscore(O['counts_on'].iloc[:200], **kw)

    def jitter(O):
        rng = np.random.default_rng(5)
        c = O['counts_on'].iloc[:100].to_numpy(dtype=float)
        return np.maximum(c[None, None] + rng.integers(-1, 2, (2, 3) + c.shape), 0.0)

    def hold(O, **kw):
        from experiments import E77_input_screen as S
        size = 200 * 16
        sim = F.toy_simulator()
        stub = {'p': dict(F.TOY_PARAMS), 'target_on': 160, 'target_off': 420, 'size': size, 'wnorm': None}
        orig = S.c2_parts
        S.c2_parts = lambda: stub
        try:
            res = S.hold_rates(lambda: sim, **kw)
        finally:
            S.c2_parts = orig
        res.pop('seconds')
        return res

    def population(O):
        sims = {lab: F.toy_simulator(seed=s) for lab, s in (('real', 0), ('P1', 1), ('M1', 2))}
        ratios = {'real': (0.8, 1.2), 'P1': (0.7, 1.0), 'M1': (0.9, 1.5)}
        theta = dict(fixed_thr=25.0, wAPLKC_scale=0.6, wAPLPN_scale=1.2, wKCKC_scale=0.2)
        outs = E72.simulate_population(sims, ratios, theta)
        ds = [{c: H.pair_distances(o[c][1].set_axis(O['odors'], axis=1)) for c in ('on', 'off')} for o in outs]
        return outs, E72.mean_distances(ds)

    def bisect(O, target, bracket, step):
        sim = F.toy_simulator(seed=7)
        p = dict(F.TOY_PARAMS)
        return E69.bisect_factor(lambda f: sim(dict(p, fixed_thr=p['fixed_thr'] * f)), target, 200 * 16, bracket, step,
            'frozen')

    def by_class(O):
        on = [O['trials_ctrl'].mean(axis=0), O['trials_ctrl2'].mean(axis=0)]
        off = [O['trials_tnt'].mean(axis=0)]
        return H.apl_effect_by_class(on, off)

    c: Cases = {
        'readout': (lambda O: [H.readout(O['counts_on'], k, p) for k, p in (('raw', None), ('sqrt', None),
            ('power', 0.5), ('saturating', 3.0), ('log1p', 2.0), ('binary', None), ('binary', 3.0))]
            + [H.readout(O['counts_on'], 'saturating', 3.0, floor=0.5), H.readout(O['counts_on'], 'saturating', 4.0,
                floor=-1.0), H.readout(O['counts_on'], 'rank_match', lut(O))], False),
        'readout_pair_distances': (lambda O: [H.readout_pair_distances(O['counts_on'], 'raw'),
            H.readout_pair_distances(O['counts_off'], 'saturating', 3.0),
            H.readout_pair_distances(O['counts_on'], 'saturating', 3.0, floor=0.5),
            H.readout_pair_distances(O['counts_on'], 'binary', 2.0)], False),
        'readout_stat': (lambda O: [H.readout_stat(O['counts_on'], 'saturating', 3.0),
            H.readout_stat(O['counts_on'], 'saturating', 3.0, 'top1_share'),
            H.readout_stat(O['counts_on'], 'saturating', 3.0, responders_min=2.0),
            H.readout_stat(O['counts_on'], 'power', 0.5, 'responder_cv'),
            H.readout_stat(O['counts_on'], 'saturating', 3.0, floor=0.5, responders_min=2.0)], False),
        'calibrate_readout': (lambda O: [H.calibrate_readout(O['counts_on'], 'saturating', O['cal_target'],
            responders_min=2.0), H.calibrate_readout(O['counts_on'], 'power', 0.05),
            H.calibrate_readout(O['counts_on'], 'log1p', 0.04, stat='participation_ratio'),
            H.calibrate_readout(O['counts_on'], 'saturating', 0.5, stat='responder_cv'),
            H.calibrate_readout(O['counts_on'], 'saturating', O['cal_target'], responders_min=2.0, floor=0.3)], True),
        'rank_match_map': (lambda O: lut(O), False),
        'kernel_window_weights': (lambda O: [H.kernel_window_weights(np.linspace(-1, 10, 45), r, d, **kw)
            for r, d, kw in ((0.5, 3.0, {}), (0.0, np.inf, {}), (1.0, 5.0, dict(window=(2.0, 6.0))),
                (0.3, 20.0, dict(ref_time=0.5)), (0.0, 4.0, {}), (2.0, np.inf, {}))], False),
        'kernel_weighted_counts': (lambda O: [H.kernel_weighted_counts(O['spikes'], pd.Index(range(O['n_kc'])),
            pd.Index(O['odors']), 1.0, 5.0), H.kernel_weighted_counts(O['spikes'], pd.Index(range(O['n_kc'])),
            pd.Index(O['odors']), 0.0, np.inf, stim_onset=0.5, window=(2.0, 8.1))], False),
        'expected_pair_distances_noisy': (lambda O: [H.expected_pair_distances_noisy(O['counts_on'], 0.0),
            H.expected_pair_distances_noisy(O['counts_on'], 0.7, 3), H.expected_pair_distances_noisy(O['tm_mat'], 0.4)],
            False),
        'sample_pseudo_recordings': (lambda O: [H.sample_pseudo_recordings(O['counts_on'].iloc[:50], 0.4,
            n_recordings=3, seed=1), H.sample_pseudo_recordings(O['counts_on'].iloc[:40].to_numpy(), 0.3, n_recordings=2,
            n_trials=2, seed=2, return_trials=True), H.sample_pseudo_recordings(jitter(O)[:, :, :30], 0.5, seed=3)],
            False),
        'session_zscore': (lambda O: [z(O), z(O, amp='saturating', K=5.0), z(O, center=True, n_blank=2),
            z(O, divisor=3.0), z(O, drift=0.5, noise_ratio=0.3), H.session_zscore(jitter(O), q=0.128, noise=1.5)],
            False),
        'zcell_pair_distances': (lambda O: [H.zcell_pair_distances(z(O), r) for r in ('graded', 'thresholded',
            'binary')] + [H.zcell_pair_distances(H.session_zscore(jitter(O), q=0.128, noise=1.5), 'thresholded',
                thresh=0.8)], False),
        'zcell_stats': (lambda O: [H.zcell_stats(z(O)), H.zcell_stats(H.session_zscore(jitter(O), q=0.128, noise=1.5)),
            H.zcell_stats(z(O, amp='saturating', K=5.0), thresh=0.8)], False),
        'calibrate_session_zscore': (lambda O: H.calibrate_session_zscore(O['counts_on'], q=0.128, amp='linear',
            target_rate=0.068), True),
        'calibrate_session_zscore_sat': (lambda O: H.calibrate_session_zscore(O['counts_on'].iloc[:150], q=0.128,
            amp='saturating', target_rate=0.068, target_amp=2.5, target_silent_sd=0.4, n_K_grid=9, max_iter=4), True),
        'response_fractions': (lambda O: H.response_fractions(O['resp_on']), False),
        'model_rdms': (lambda O: H.model_rdms(O['resp_on'], O['counts_on']), False),
        'amplitude_stats': (lambda O: [H.amplitude_stats(O['counts_on']), H.amplitude_stats(O['tm_mat'])], False),
        'responder_cv': (lambda O: [H.responder_cv(O['counts_on'].to_numpy()[O['counts_on'].to_numpy() >= 2]),
            H.responder_cv(np.array([np.nan, 1.0, 3.0])), H.responder_cv(np.array([]))], False),
        'poisson_binomial_bins': (lambda O: [H.poisson_binomial_bins(np.clip(O['tm_mat'].to_numpy() / 6, 0, 1)),
            H.poisson_binomial_bins(np.clip(O['tm_mat'].to_numpy() / 6, 0, 1), edges=((0, 1), (2, None)))], False),
        'opk_class': (lambda O: [H.opk_class_labels(), H.opk_class_labels(((0, 2), (3, None))),
            H.opk_class_shares(O['resp_on']), H.opk_class_shares(O['resp_off'].to_numpy() > 0)], False),
        'class_cosine_terms_array': (lambda O: H.class_cosine_terms_array(O['tm_mat'].to_numpy(),
            O['tm_mat'].to_numpy() > 1.0), False),
        'class_cosine_terms': (lambda O: [H.class_cosine_terms(O['tm_mat']),
            H.class_cosine_terms(H.readout(O['counts_on'], 'saturating', 3.0), resp=O['counts_on'] >= 2),
            H.class_cosine_terms(O['trials_ctrl'].mean(axis=0), thresh=0.8, edges=((0, 1), (2, None)))], False),
        'cosine_terms': (lambda O: [H.cosine_terms(O['tm_mat']), H.cosine_terms(O['counts_on'], thresh=2.5)], False),
        'apl_effect': (lambda O: H.apl_effect(*_model(O)), False),
        'apl_effect_stats': (lambda O: H.apl_effect_stats(*_model(O), O['data_effect']), False),
        'apl_effect_pair_parts': (lambda O: H.apl_effect_pair_parts(*_model(O), O['data_effect']), False),
        'apl_effect_decomposition': (lambda O: H.apl_effect_decomposition(*_model(O), O['data_effect']), False),
        'apl_effect_by_class': (lambda O: by_class(O), False),
        'apl_effect_by_class_resp': (lambda O: H.apl_effect_by_class(H.readout(O['counts_on'], 'saturating', 3.0),
            H.readout(O['counts_off'], 'saturating', 3.0), resp_on=O['counts_on'] >= 2, resp_off=O['counts_off'] >= 2),
            False),
        'bisect_factor': (lambda O: [bisect(O, 500, (0.1, 10.0), 1.1), bisect(O, 2900, (0.5, 2.0), 1.25)], False),
        'hold_rates': (lambda O: hold(O, label='frozen'), False),
        'hold_rates_targets': (lambda O: hold(O, target_off=600, target_on=250, params=dict(F.TOY_PARAMS,
            fixed_thr=40.0), label='frozen2', apl_bracket=(0.5, 4.0)), False),
        'simulate_population': (lambda O: population(O), False),
    }
    return c
