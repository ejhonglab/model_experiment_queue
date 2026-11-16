# Yang's questions on the summary page

Each entry is a question Yang asked about one sentence (or table row) of the all-experiments summary, with the
answer given. `reports/yang_summary.py build` underlines the sentence in Yang's version of the page and shows the
question and the answer in a floating window when it is clicked. `sentence` is the sentence as it reads in Yang's
version (correlations), `sentence_main` as it reads on the main page (distances); numbers in either may be out of
date, every number is matched as a wildcard. Answers are written in correlations, with the main page's values beside
them where that helps, and take their numbers from placeholders.

## noise-floor-method
asked: 2026-11-16
sentence: From the three trials of each odor, the correlation between two independent trial means of the same odor is the noise ceiling (0.77 in control flies, 0.81 with APL output blocked (TNT); between single trials 0.56 and 0.61).
sentence_main: From the three trials of each odor, the distance between two independent trial means of the same odor is the noise floor (0.23 in control flies, 0.19 with APL output blocked (TNT); between single trials 0.44 and 0.39).
question: Is the method to calculate the cosine distance between every pair of trials of the same odor, average across the 3 pairs (since there are 3 repeats), then average across odors?

Not quite. That is how the other number in the sentence is made, the one between single trials ({{inv:data:ctrl_single_trial_floor ~ .2f}} in control flies and {{inv:data:tnt_single_trial_floor ~ .2f}} with APL output blocked (TNT)). The {{inv:data:ctrl_trial_mean_floor ~ .2f}} and {{inv:data:tnt_trial_mean_floor ~ .2f}} are not measured from pairs of trials: they are an estimate of the correlation that two independent three-trial means of one odor would have.

**The value between trial means**, per recording and odor (X is 3 trials × KCs):

1. Take the trial mean x̄ over the 3 trials, and its squared norm |x̄|².
2. Take the noise power s²: each KC's variance across the 3 trials (n − 1 in the denominator), summed over KCs.
3. Averaging 3 trials with independent noise leaves s²/3 of noise power in |x̄|², so the signal power is |a|² = |x̄|² − s²/3.
4. Two independent three-trial means share the signal a and each carries its own noise of power s²/3, so their expected cosine similarity is |a|² / |x̄|². The main page gives the same thing as a distance, (s²/3) / |x̄|².
5. Average over the 16 odors within a recording, then over the recordings (11 control, 8 with APL output blocked (TNT)).

It has to be an estimate: only 3 trials exist, and they cannot be split into two separate three-trial means.

**The value between single trials** is the one described in the question: the cosine similarity of each of the 3 pairs of trials of an odor, averaged over the 3 pairs, then over the odors of a recording, then over recordings.

Two cautions on the trial-mean value. It assumes that the noise of the 3 trials of one odor is independent, and it is the ratio of two expected values, an approximation of the expected cosine. Code: `analyse` in kcq/noise_floor.py; the registry entries are made in kcq/data_targets.py.

On the main page these are distances: {{data:ctrl_trial_mean_floor ~ .2f}} and {{data:tnt_trial_mean_floor ~ .2f}} between trial means, {{data:ctrl_single_trial_floor ~ .2f}} and {{data:tnt_single_trial_floor ~ .2f}} between single trials.

## earlier-correction-too-low
asked: 2026-11-16
sentence: the same by the earlier correction, which takes the noise of two odors to be independent (too high)
sentence_main: the same by the earlier correction, which takes the noise of two odors to be independent (too low)
question: "Too low" means the corrected pair distance is too low? And the noise is not independent because?

Yes. On the main page the row is a mean pair distance, and the earlier correction put it too low. Here the same row is a mean pair correlation, so the earlier correction put it too high: it removed more than the trial noise.

| statistic | control flies | APL output blocked (TNT) |
|---|---|---|
| mean pair correlation as measured, on the three-trial means | {{inv:data:ctrl_mean_dist ~ .2f}} | {{inv:data:tnt_mean_dist ~ .2f}} |
| without the trial noise, by the estimate across repeats | {{inv:data:noise_free_mean_dist_ctrl ~ .2f}} | {{inv:data:noise_free_mean_dist_tnt ~ .2f}} |
| without the trial noise, by the earlier correction (too high) | {{inv:data:ctrl_mean_dist_noise_corrected ~ .2f}} | {{inv:data:tnt_mean_dist_noise_corrected ~ .2f}} |

**Why the noise of two odors is not independent.** In the recordings, the responses to two different odors presented in the same repeat of the panel are more alike than the responses to the same two odors in different repeats. Between single trials the correlation of two different odors is {{inv:E81/distances_by_kind.csv:units=z;values=graded;genotype=TNT3c;level=genotype;pairset=all;kind=same_repeat:value ~ .3f}} in the same repeat against {{inv:E81/distances_by_kind.csv:units=z;values=graded;genotype=TNT3c;level=genotype;pairset=all;kind=different_repeat:value ~ .3f}} in different repeats in control flies, and {{inv:E81/distances_by_kind.csv:units=z;values=graded;genotype=TNT;level=genotype;pairset=all;kind=same_repeat:value ~ .3f}} against {{inv:E81/distances_by_kind.csv:units=z;values=graded;genotype=TNT;level=genotype;pairset=all;kind=different_repeat:value ~ .3f}} with APL output blocked (TNT). With independent noise the two would be equal. The difference, the same-repeat gap, is {{data:repeat_gap_ctrl ~ .3f}} and {{data:repeat_gap_tnt ~ .3f}}, and it is above zero in every recording.

So part of a KC's trial-to-trial variability is common to the odors of a repeat: {{E81/variance_parts.csv:units=z;genotype=TNT3c;level=genotype;cls=responders;statistic=share_shared:value ~ .0%}} of the trial-to-trial variance of a responding KC's response in control flies and {{E81/variance_parts.csv:units=z;genotype=TNT;level=genotype;cls=responders;statistic=share_shared:value ~ .0%}} with APL output blocked (TNT). About three quarters of the gap follows the time between the two presentations ({{E81/order_gap_decomposition.csv:units=z;basis=positions;genotype=TNT3c;statistic=separation_share:value ~ .0%}} and {{E81/order_gap_decomposition.csv:units=z;basis=positions;genotype=TNT;statistic=separation_share:value ~ .0%}}), mostly a slow drift over the session, and the rest is a step tied to the repeat itself. The page measures this and does not say what in the fly or the preparation causes it.

**Why that biases the earlier correction.** It divided the product of the two odors' trial means, x̄ₐ·x̄ᵦ, by the two signal norms (|x̄|² − s²/3 for each odor), on the assumption that noise adds nothing to x̄ₐ·x̄ᵦ. The two trial means are taken over the same three repeats, so the noise the two odors share within each repeat does add to that product. The correction counted it as similarity between the odors while taking all the noise out of the norms. The corrected correlation comes out too high, most with APL output blocked (TNT), where the shared part is larger, and the mean APL effect too high with it.

**What replaced it.** The estimate across repeats uses only products of trials from different repeats, which the shared part does not enter: the cosine between single trials of the two odors from different repeats, divided by the geometric mean of the two odors' cosines between their own trials across the same repeats. The section on the flies' variability across the three repeats has the details ([there](#repeat-structure)).

## orn-kc-level-comparison
asked: 2026-11-16
sentence: Twice as many KCs respond to each odor (0.150 instead of 0.068 of the KCs, with the per-odor pattern conserved), and the mean pair correlation rises from 0.16 to 0.34, above the ORN input's own 0.29.
sentence_main: Twice as many KCs respond to each odor (0.150 instead of 0.068 of the KCs, with the per-odor pattern conserved), and the mean pair distance falls from 0.84 to 0.66, below the ORN input's own 0.71.
question: I think this might not be a fair comparison: the ORNs use a dense code with about 50 units, whereas the KC space is sparse and high dimensional. Using cosine to measure the correlation of the ORN representations and comparing between the two spaces doesn't seem fair. Am I wrong?

No, the comparison with the ORN input is not a fair one, and it does not hold once both sides are measured the same way. Two details differ from the question, and they make the point stronger. The ORN values on this page are Pearson correlations across the 36 glomeruli imaged in both ORN panels, not cosine similarities. The KC values are cosine similarities across all KCs of a recording. The clause therefore sets two different measures of two different codes side by side.

| mean pair correlation (120 pairs) | ORNs | KCs, control flies | KCs, APL output blocked (TNT) |
|---|---|---|---|
| as on this page: Pearson for the ORNs, cosine similarity for the KCs | {{inv:data:orn_mean_dist ~ .2f}} | {{inv:data:ctrl_mean_dist ~ .2f}} | {{inv:data:tnt_mean_dist ~ .2f}} |
| Pearson correlation for all three | {{inv:data:orn_mean_dist ~ .2f}} | {{lit:0.14 # one minus the mean over the 120 pairs of H.yang_mean_rdm('TNT3c', 'pearson') (0.858), computed 2026-11-16 in viki's session for Yang's question; not a registry entry}} | {{lit:0.24 # one minus the mean over the 120 pairs of H.yang_mean_rdm('TNT', 'pearson') (0.759), computed 2026-11-16 in viki's session for Yang's question; not a registry entry}} |
| cosine similarity for all three | {{lit:0.65 # the mean over the 120 pairs of the cosine similarity between the two odors' fly-mean ORN vectors (36 glomeruli, the vectors of data_prep/orn16_distances_from_olfdb.py; mean cosine distance 0.351), computed 2026-11-16 in viki's session for Yang's question; not a registry entry}} | {{inv:data:ctrl_mean_dist ~ .2f}} | {{inv:data:tnt_mean_dist ~ .2f}} |

**With one measure on both sides, the KCs stay less correlated than the ORN input even with APL output blocked (TNT)**: 0.24 against 0.29 by Pearson, 0.34 against 0.65 by cosine. "Above the ORN input's own" is true only of the mixed pairing in the first row.

**Why the measure matters so much for the ORNs.** The ORN code is dense and positive: 96 % of the glomerulus responses are above zero. Every pair of odors then shares a common positive part, which the cosine counts as similarity and the Pearson correlation removes. The KC responses are sparse and z-scored, centred near zero, so the two measures nearly agree in control flies and differ more with APL output blocked (TNT), where twice as many KCs respond.

**Levels are not comparable across the two stages in any case.** The level of a correlation depends on how dense the code is, on the number of units and on noise. Trial noise lowers the KC values (without it, by the estimate across repeats, {{inv:data:noise_free_mean_dist_ctrl ~ .2f}} and {{inv:data:noise_free_mean_dist_tnt ~ .2f}}), and the ORN vectors are means over flies, which carry less of it.

**What can be compared.** Control flies against APL output blocked (TNT): the same kind of cells and the same measure, and the mean pair correlation about doubles. And the order of the pairs, which no change of scale touches: ρ on odor-odor correlations with the ORNs is {{data:rho_ctrl_orn ~ .2f}} in control flies and {{data:rho_tnt_orn ~ .2f}} with APL output blocked (TNT). The page's statistics are of these two kinds.

The Pearson values for the KCs and the cosine value for the ORNs are not on the main page; they were computed for this answer from the same recordings and the same ORN vectors.
