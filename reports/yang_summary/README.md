# Yang's version of the all-experiments summary

Yang asked (2026-11-16, in viki's terminal) for her own version of the summary page with two changes:

1. every odor-odor distance given as a correlation (correlation = 1 − distance: the cosine similarity for KCs and
   models, the Pearson correlation for ORNs);
2. every sentence she has asked a question about underlined, the question and its answer opening in a floating
   window when the sentence is clicked.

It is built from the same template as the main page (`reports/summary_template.md`), which is never edited for it.

## Files here

| file | what it is |
|---|---|
| `questions.md` | Yang's questions and the answers, one entry per question, anchored to a sentence of the page |
| `intro.md` | the note under the page's title that says what this version changes |
| `conversion_rules.md` | how a line of the template is reworded from distances to correlations |
| `correlation_lines.jsonl` | the converted lines, keyed by the sha1 of the template line each one replaces |
| `inhibition_short.md` | the short report on how implementations of inhibition change the KC code (its own page) |

The tool is `reports/yang_summary.py` (tests in `test/test_yang_summary.py`).

## Building the page

```
./venv/bin/python reports/yang_summary.py build
```

writes `outputs/_artifacts/summary_yang/index.html` and `build_report.json` beside it (lines converted, lines not yet
converted, which questions found their sentence). A template line that changed since it was converted falls back to
the main page's wording and is tagged "distances" on the page until it is converted again.

## Adding a question

Add an entry to `questions.md` (the format is at the top of that file and in the tool's docstring): the sentence as
it reads in Yang's version, the sentence as it reads on the main page, the question in Yang's words, and the answer
in correlations with its numbers as placeholders. Then build. A question Yang withdraws is removed from the file.

## Converting lines after the template changed

```
./venv/bin/python reports/yang_summary.py status     # what is not converted
./venv/bin/python reports/yang_summary.py export     # chunks under outputs/_artifacts/summary_yang/convert/
# reword the starred lines of each chunk into chunk_NN.out.md, by conversion_rules.md
./venv/bin/python reports/yang_summary.py import --by "<who converted them>"
./venv/bin/python reports/yang_summary.py build
```

When template lines were only split (between sentences, into paragraphs, a bullet list or sub-bullets, no word
changed), carry their entries over to the pieces first, so that only what changed is converted anew:

```
./venv/bin/python reports/yang_summary.py carry --dry-run   # what would be carried, what not
./venv/bin/python reports/yang_summary.py carry             # writes the entries; the rest in convert/carry_failures.md
```

Rewording is summary text, so on a queue host it is a Fable subagent's work (CLAUDE.md, "Permissions by model and
effort"); the tool's code is Opus work.

## What this version does not change

Figures are the main page's and are drawn in distances, except the first figure, which is drawn again in correlations
(`results/summary_yang/summary_background_correlation.png` in the shared directory, by the build). The progress plot
keeps its own labels ("ρ on odor-odor distances").

## The scheme of the baseline model

Yang asked (2026-11-17) for a scheme of the original model: its connectivity and the transformation at each step.
`reports/yang_model_scheme.py` (tests in `test/test_yang_model_scheme.py`) draws it to
`results/summary_yang/model_scheme.png` in the shared directory, with a provenance file beside it. Its counts are
read from the results tables and its constants are checked against the simulator's source at every run. The answer
`baseline-model-scheme` in `questions.md` and the short report embed the figure; run the script again when the
baseline's wiring or the simulator's defaults change.

## The short report on inhibition

Yang asked (2026-11-17) for a short version of the report on one question: how different implementations of
inhibition change the KC code, against the flies. Its text is `inhibition_short.md`, with the same placeholders as
the answers in `questions.md`. `reports/yang_short_report.py build` (tests in `test/test_yang_short_report.py`) fills
and converts it to `outputs/_artifacts/inhibition_short_yang/index.html`. At her request it does not dwell on the
fruity × small-ester block.

## Provenance

- The tool and its tests: written by a code subagent, claude-opus-5-5 at xhigh, spawned from viki 2026-11-16.
- The conversion rules, the questions' answers and the note under the title: written by viki's session
  (claude-fable-5-2) at Yang's request, 2026-11-16.
- The scheme's script and the short report's build tool, with their tests: written by code subagents, claude-opus-5-5
  at xhigh, spawned from viki 2026-11-17; the scheme's wording and the short report's text by viki's session
  (claude-fable-5-2).
- The converted lines: reworded by Fable subagents (claude-fable-5-2, spawned from viki 2026-11-16) following
  `conversion_rules.md`; each entry of `correlation_lines.jsonl` records who converted it and when.
