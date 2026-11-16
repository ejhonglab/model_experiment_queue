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

Rewording is summary text, so on a queue host it is a Fable subagent's work (CLAUDE.md, "Permissions by model and
effort"); the tool's code is Opus work.

## What this version does not change

Figures are the main page's and are drawn in distances, except the first figure, which is drawn again in correlations
(`results/summary_yang/summary_background_correlation.png` in the shared directory, by the build). The progress plot
keeps its own labels ("ρ on odor-odor distances").

## Provenance

- The tool and its tests: written by a code subagent, claude-opus-5-5 at xhigh, spawned from viki 2026-11-16.
- The conversion rules, the questions' answers and the note under the title: written by viki's session
  (claude-fable-5-2) at Yang's request, 2026-11-16.
- The converted lines: reworded by Fable subagents (claude-fable-5-2, spawned from viki 2026-11-16) following
  `conversion_rules.md`; each entry of `correlation_lines.jsonl` records who converted it and when.
