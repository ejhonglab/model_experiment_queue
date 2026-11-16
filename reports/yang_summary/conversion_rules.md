# Converting the summary's distances to correlations (Yang's version)

Yang's version of the all-experiments summary gives every odor-odor distance as a correlation,
correlation = 1 − distance, and changes nothing else. These are the rules for rewording template lines; they are
followed by whoever converts a chunk written by `reports/yang_summary.py export` (a language-model agent, or a person).
The numbers themselves are never typed: a value is a token, and a token is marked to be shown as 1 − value.

## What you get and what you return

A chunk file lists template lines as `@@ <line number> <flag>` followed by the line. A placeholder appears as a token
`⟨k|VALUE|HINT⟩`: k numbers the placeholders of that line, VALUE is what the main page prints, HINT says what the
value is. Flag `*`: the line must be answered. Flag `=`: already converted, shown as it now reads (context; shows the
wording in use). No flag: context that probably needs nothing.

Write the answers to `<chunk file name without .md>.out.md` in the same folder, one entry per line answered:

```
@@ 168
<the whole converted line, on one line>
@@ 171 KEEP
```

- Every starred line gets an entry: the converted line, or `KEEP` when it needs no change.
- An unstarred line gets an entry only if it needs conversion (a row of numbers under a converted header, a sentence
  about pairs being "apart" that the star rule missed).
- Tokens: use every k of the line exactly once. Copy a token unchanged (or write `⟨k⟩`) to keep the value;
  `⟨k|inv⟩` shows 1 − value (a distance becomes a correlation); `⟨k|neg⟩` shows −value (a signed change of a distance
  becomes a change of a correlation); `⟨k|lit|new text⟩` replaces a hand-entered value (empty HINT) that cannot be
  inverted by arithmetic. Never type a value that a token carries.
- First line of your final message: your exact model ID. Then one line per chunk: how many lines converted, how many
  KEEP, and any line you were unsure of with the reason. Edit no other file.

## The one change

An odor-odor distance d, between two odors' response vectors at any stage (ORNs, PNs, boutons, claws, KCs, a model's
KCs, single trials or trial means), becomes the correlation c = 1 − d.

- The word is "correlation": pair correlation, mean pair correlation, odor-odor correlations, correlation matrix,
  KC correlations. Where the text names the function: "cosine distance" becomes "cosine similarity", "one minus the
  cosine similarity" becomes "the cosine similarity", "one minus the Pearson correlation" becomes "the Pearson
  correlation", "Pearson correlation distance" becomes "Pearson correlation".
- A token that is a distance gets `|inv`. A distance typed in the text is rewritten by hand as 1 − d with the same
  number of decimals (0.88 becomes 0.12; "more than 0.8 apart" becomes "with a correlation below 0.2").
- Every statement of direction about a level reverses with it: higher, rises, above, the top, largest, far apart,
  pushed apart, near-maximal distance become lower, falls, below, the bottom, smallest, nearly uncorrelated,
  decorrelated, near-zero correlation, and the reverse. Check each sentence: "falls from ⟨1⟩ to ⟨2⟩, below the ORN
  input's own ⟨3⟩" becomes "rises from ⟨1|inv⟩ to ⟨2|inv⟩, above the ORN input's own ⟨3|inv⟩".
- A range swaps its ends: "⟨1⟩–⟨2⟩" (smallest to largest distance) becomes "⟨2|inv⟩–⟨1|inv⟩".
- Spatial words: "apart" becomes "decorrelated" or "less correlated"; "close", "closer", "near" (of two odors'
  responses) become "correlated", "more correlated", "alike"; "moves two odors apart" becomes "decorrelates two
  odors"; "brings them closer" becomes "makes them more correlated"; "too far apart" becomes "too weakly correlated";
  "too close" becomes "too strongly correlated".
- The noise floor of the distances (the distance between two presentations, or two trial means, of the same odor)
  becomes the "noise ceiling" of the correlations, with `|inv`; "noise inflates every cosine distance toward 1"
  becomes "noise pulls every correlation toward 0".

## What keeps its number (no `|inv`, no `|neg`)

- **A ρ, a Pearson r or an R² between two sets of pair values**, both of which are converted: the same on
  correlations as on distances. Only the word changes ("ρ on odor-odor correlations"). The same holds for the isotonic
  R² and for a root mean squared difference between two sets of distances (a loss).
- **A difference of two distances keeps its size and sign; say it with the two terms in the other order.** "A's
  distance minus B's distance" becomes "B's correlation minus A's correlation". So, with their numbers untouched:
  - the per-pair APL effect, "the distance with APL on minus the distance with APL off", becomes "the correlation with
    APL off minus the correlation with APL on" (in the KC data, with APL output blocked (TNT) minus the control
    flies'); "positive where APL moves two odors apart" becomes "positive where APL decorrelates two odors". The mean
    APL effect, the APL-effect loss, its mean, pattern, amplitude and shape parts, the spread, ρ and Pearson r on the
    per-pair APL effect, the floor and the ceilings are all unchanged; "in units of cosine distance" becomes "in
    units of correlation";
  - the block: the nine fruity × small-ester pairs "sit above the curve" becomes "sit below the curve", and "the block
    is the mean residual of the nine pairs" becomes "the block is how far the nine pairs sit below that curve on
    average (the fit's value minus the measured correlation)". Every block value stays as it is, positive;
  - the same-repeat gap, "the distance in different repeats minus the distance in the same repeat", becomes "the
    correlation in the same repeat minus the correlation in different repeats";
  - the excess shared responders, "the KC distance actually measured minus the distance expected", becomes "the
    correlation expected minus the correlation actually measured"; the sentence that reads its sign stays as written;
  - a `diff(...)` token that is one of these (the mean APL effect from two mean distances, control minus TNT) is
    unchanged.
- **A change described with a direction word** keeps its number and reverses the word: "lowered the mean distance by
  ⟨1⟩" becomes "raised the mean correlation by ⟨1⟩".
- **The structure beyond the ORN input** (a residual of one pair matrix from a fit on another) and the ρ's computed
  on it, the structure vs wild-type somata and the natural-source alignment, keep their values. Where the natural-source
  alignment is defined, say that it is the ρ between the residual and the odors' natural-source similarity (the
  natural-source distances in reverse order), the same number as on the main page; elsewhere "natural-source
  distances" and "VCF" stay as names of that reference.
- The standard deviation or spread of a set of distances is the standard deviation of the correlations.

## What needs `|neg`

- A signed change of a level that the text or a table gives with its sign and no direction word: a column "change of
  the mean distance" with cells like +0.03 becomes "change of the mean correlation" with `|neg`.
- A residual from a monotone or linear fit quoted with its sign for a pair or an odor (not the block): either `|neg`
  with "above" and "below" reversed, or keep the number and say "below the fit by"; choose whichever leaves the
  sentence plain, and do the same for every value of that table or paragraph.
- A ρ between a set of pair distances and something that is not converted (a response rate, a count of shared
  glomeruli, a separation in time): its sign reverses on correlations, so `|neg`, unless the other quantity is one of
  the unchanged differences above.

## What is not an odor-odor distance (leave it)

The total variation distance and "the distance of one recording from the mean of the others" (distributions of the
number of odors per KC); distances in space or time (µm, seconds, presentations apart, "how far apart the two odors
were presented"); a KC's distance or margin from its threshold; "the floor of the loss", "the measurement ceiling",
"the ceiling on ρ", "the firing ceiling"; "near-identical ORN input"; "decorrelation level" where it names the model's
operating level. A starred line that holds only these is `KEEP`.

## When a value cannot be said as a correlation

A ratio or a percentage change of distances, or a distance on another scale, has no correlation form. Keep the number
and say plainly that it is of the distance: "the mean distance (1 − correlation) falls by ⟨1⟩ %". Do not invent an
equivalent.

## Figures and buttons

The images are the main page's and are drawn in distances. A caption stays true to its image: keep what it says the
image shows and the directions as drawn ("the red circles sit above the curve" stays in a caption), and if the caption
mentions distances add once, at its end, "Drawn as distances (distance = 1 − correlation)." A caption that mentions no
distance is `KEEP`. `!button[...](...)` cells and figure paths never change.

## What never changes

- The claims. Nothing added, nothing dropped, the same order, the same number of sentences in a line, the same
  emphasis. Do not shorten, tidy or correct anything else, even an error; note it in your report instead.
- The structure: a heading's level and `{#id}`, a list marker and its indentation, the number of cells of a table row
  and its `{span=N}`, every link target, a hover definition's term before the colon, `??? ` at the start of a fold's
  title.
- The page's wording rules: third person, the character ρ, the genotype labels "control flies" and "APL output
  blocked (TNT)", odors by their abbreviations, no experiment IDs.

## Examples

Source (a list item):

```
- **Noise floor and the distance without the trial noise:** two noisy response vectors are never parallel, so noise inflates every cosine distance toward 1. From the three trials of each odor, the distance between two independent trial means of the same odor is the noise floor (⟨1|0.23|ctrl_trial_mean_floor⟩ in control flies, ⟨2|0.19|tnt_trial_mean_floor⟩ with APL output blocked (TNT); between single trials ⟨3|0.44|ctrl_single_trial_floor⟩ and ⟨4|0.39|tnt_single_trial_floor⟩).
```

Converted:

```
- **Noise ceiling and the correlation without the trial noise:** two noisy response vectors are never parallel, so noise pulls every correlation toward 0. From the three trials of each odor, the correlation between two independent trial means of the same odor is the noise ceiling (⟨1|inv⟩ in control flies, ⟨2|inv⟩ with APL output blocked (TNT); between single trials ⟨3|inv⟩ and ⟨4|inv⟩).
```

Source (table rows):

```
| mean pair distance (120 pairs) | ⟨1|0.48|orn_mean_dist⟩ | ⟨2|0.84|ctrl_mean_dist⟩ | ⟨3|0.66|tnt_mean_dist⟩ |
| range of pair distances | ⟨1|0.05|orn_min_dist⟩–⟨2|1.16|orn_max_dist⟩ | ⟨3|0.34|ctrl_min_dist⟩–⟨4|1.00|ctrl_max_dist⟩ | ⟨5|0.18|tnt_min_dist⟩–⟨6|0.89|tnt_max_dist⟩ |
| pairs more than 0.8 apart | ⟨1|55|orn_n_gt08⟩ | ⟨2|92|ctrl_n_gt08⟩ | ⟨3|17|tnt_n_gt08⟩ |
| ρ on odor-odor correlations with the ORNs | — | ⟨1|0.65|rho_ctrl_orn⟩ | ⟨2|0.73|rho_tnt_orn⟩ |
```

Converted, as the entries would be returned if these were lines 173, 177, 178 and 181 (the last row needs nothing):

```
@@ 173
| mean pair correlation (120 pairs) | ⟨1|inv⟩ | ⟨2|inv⟩ | ⟨3|inv⟩ |
@@ 177
| range of pair correlations | ⟨2|inv⟩–⟨1|inv⟩ | ⟨4|inv⟩–⟨3|inv⟩ | ⟨6|inv⟩–⟨5|inv⟩ |
@@ 178
| pairs with a correlation below 0.2 | ⟨1⟩ | ⟨2⟩ | ⟨3⟩ |
@@ 181 KEEP
```

Source (a paragraph with typed numbers and a named difference):

```
**The one pair-specific exception: the block.** The nine fruity × small-ester pairs, ⟨1|0.58|orn_fruity_small_mean⟩ apart in the ORN input and 0.88 in control flies, fall to 0.60 when APL output is blocked (TNT), well below the monotone trend. That residual, ⟨2|+0.071|block⟩ (90 % confidence interval ⟨3|+0.045|block.lo⟩ to ⟨4|+0.098|block.hi⟩), is the block referred to throughout; it means APL removes a population of KCs that respond to both ester groups.
```

Converted:

```
**The one pair-specific exception: the block.** The nine fruity × small-ester pairs, with a correlation of ⟨1|inv⟩ in the ORN input and 0.12 in control flies, rise to 0.40 when APL output is blocked (TNT), well above the monotone trend. That shortfall of the control flies' correlation below the trend, ⟨2⟩ (90 % confidence interval ⟨3⟩ to ⟨4⟩), is the block referred to throughout; it means APL removes a population of KCs that respond to both ester groups.
```

Source and converted (the main target's definition):

```
For each of the 120 pairs, the per-pair APL effect is the pair's distance with APL on minus its distance with APL off: in the KC data, the control flies' distance minus the distance with APL output blocked (TNT), from the two KC matrices of the first figure (cosine distance, mean over recordings); in a model, the distance with APL on minus the distance of the same model with APL removed and nothing re-tuned (the acute knockout), as cosine distance on spike counts. It is positive where APL moves two odors apart.
```

```
For each of the 120 pairs, the per-pair APL effect is the pair's correlation with APL off minus its correlation with APL on: in the KC data, the correlation with APL output blocked (TNT) minus the control flies' correlation, from the two KC matrices of the first figure (cosine similarity, mean over recordings); in a model, the correlation of the same model with APL removed and nothing re-tuned (the acute knockout) minus its correlation with APL on, as cosine similarity on spike counts. It is positive where APL decorrelates two odors.
```
