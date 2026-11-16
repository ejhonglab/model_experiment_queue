#!/usr/bin/env python3
"""Builds the Questions for Yang page: the questions markdown as a page, with a form to answer them.

The queue writes its questions about Yang's recordings to the shared directory's
questions_for_yang/questions_<YYYY-MM-DD_HHMM>.md. This converts the newest of them (or the one given)
with md_to_artifact.py's converter (numbered questions keep their sub-bullets inside them) and adds an
answer form below the questions: an optional name, one box per numbered question, labelled "Question
N: <its bold title>" (the titles as the queue monitor reads them, question_titles in
coordination/monitor/server.py), and the questions file's name, sent to the monitor's /api/answers.
The monitor saves each submission as a new text file next to the questions,
answers_<questions stamp>__<YYYY-MM-DD_HHMMSS>.txt. The form works only where the monitor serves the
page (its Report tab, "Questions for Yang"); a copy opened anywhere else says so when sent.

Then save the copy the monitor shows (not published as an artifact):
    ./venv/bin/python reports/save_report_copy.py outputs/_artifacts/questions_for_yang/index.html \\
        questions_for_yang --current current-questions-for-yang.html
Usage: ./venv/bin/python reports/build_questions_page.py [QUESTIONS_MD] [--out HTML]
"""
from __future__ import annotations

import argparse
import html as html_mod
import json
import os
import re
import sys
from pathlib import Path

from coordination.monitor import server as monitor
from reports import md_to_artifact as M

REPO = Path(__file__).resolve().parents[1]
OUT_HTML = REPO / 'outputs' / '_artifacts' / 'questions_for_yang' / 'index.html'
# where the monitor serves the page (coordination/monitor/README.md)
REPORT_SITE = 'http://hal:{}'.format(monitor.DEFAULT_PORT)

FORM_CSS = """
  :root { --ok: #1a7f37; --bad: #b3261e; }
  @media (prefers-color-scheme: dark) { :root { --ok: #5cc27a; --bad: #f0766c; } }
  form.answers .field { margin: 16px 0; }
  form.answers label { display: block; font-weight: 600; margin-bottom: 4px; }
  form.answers input[type=text], form.answers textarea { display: block; width: 100%; box-sizing: border-box;
    font: inherit; color: var(--ink); background: var(--bg); border: 1px solid var(--line); border-radius: 6px;
    padding: 6px 8px; }
  form.answers input[type=text] { max-width: 24em; }
  form.answers textarea { min-height: 6em; resize: vertical; }
  form.answers input:focus, form.answers textarea:focus { outline: 2px solid var(--accent); outline-offset: -1px; }
  form.answers .actions { display: flex; flex-wrap: wrap; align-items: center; gap: 8px 14px; margin: 20px 0 6px; }
  form.answers button { font: inherit; font-weight: 600; padding: 7px 18px; border-radius: 6px;
    border: 1px solid var(--accent); background: var(--accent); color: #fff; cursor: pointer; }
  form.answers button:disabled { opacity: 0.6; cursor: progress; }
  #answers-status.ok { color: var(--ok); } #answers-status.bad { color: var(--bad); }
  p.formnote { color: var(--muted); font-size: 0.9em; }
"""

# sends the form with fetch (asking for JSON) so the page stays as it is, the text in the boxes
# included, and says what was saved; without it the form posts as a plain form and the monitor
# replies with a page
FORM_JS = """<script>
(function () {
  var f = document.getElementById('answers-form');
  if (!f || !window.fetch || !window.FormData || !window.URLSearchParams) return;
  var st = document.getElementById('answers-status'), btn = f.querySelector('button[type=submit]');
  function say(text, cls) { st.textContent = text; st.className = cls || ''; }
  f.addEventListener('submit', function (e) {
    e.preventDefault();
    var boxes = f.querySelectorAll('textarea'), any = false;
    for (var k = 0; k < boxes.length; k++) if (boxes[k].value.trim()) any = true;
    if (!any) { say('Nothing to send: every answer box is empty.', 'bad'); return; }
    btn.disabled = true; say('Sending...');
    fetch(f.getAttribute('action'), {method: 'POST', headers: {'Accept': 'application/json',
        'Content-Type': 'application/x-www-form-urlencoded;charset=UTF-8'},
        body: new URLSearchParams(new FormData(f)).toString()})
      .then(function (r) { return r.json().catch(function () { throw new Error('HTTP ' + r.status); }); })
      .then(function (d) {
        if (d.ok) say('Saved: ' + d.answers + (d.answers === 1 ? ' answer' : ' answers') + ' in ' + d.file
          + '. The text stays in the boxes; sending again saves a new file.', 'ok');
        else say('Not saved: ' + d.error, 'bad');
      })
      .catch(function () { say(@UNREACHABLE@, 'bad'); })
      .then(function () { btn.disabled = false; });
  });
})();
</script>"""


def shared_dir() -> Path:
    """The shared directory, found as kcq.harness.shared_dir finds it (KC_REMAPPING_SHARED, else
    Tom/claude_kc_remapping under /local, /mnt or /media tensor), without importing the model's packages."""
    for c in [os.environ.get('KC_REMAPPING_SHARED')] + [
            f'/{root}/tensor/Tom/claude_kc_remapping' for root in ('local', 'mnt', 'media')]:
        if c and Path(c).is_dir():
            return Path(c)
    raise FileNotFoundError('shared directory not found; set KC_REMAPPING_SHARED')


def newest_questions(folder: Path) -> Path:
    """The newest questions_<YYYY-MM-DD_HHMM>.md in `folder`, by the time in its name."""
    files = sorted(p for p in folder.iterdir() if monitor.QUESTIONS_FILE.fullmatch(p.name))
    if not files:
        raise FileNotFoundError(f'no questions_<YYYY-MM-DD_HHMM>.md in {folder}')
    return files[-1]


def answer_form(questions_name: str, titles: dict) -> str:
    """The answer form's HTML (heading, notes, fields, button and its script) for the questions
    `titles` ({number: title}, as monitor.question_titles gives them) of the file `questions_name`."""
    e = html_mod.escape
    site = f'<a href="{REPORT_SITE}/#{monitor.ANSWERS_VIEW}" target="_top">{REPORT_SITE}</a>'
    out = ['<h2 id="answers">Answers</h2>',
        '<p>Write in the box of any question you want to answer, then press "Send answers"; empty boxes '
        'are left out. The answers are saved as a text file in the queue\'s shared directory, in '
        f'{e(monitor.ANSWERS_DIR)}/ next to this file of questions, named with the time of the questions '
        'and the time the answers were sent. They can be sent more than once: each submission is saved '
        'as a new file, so an answer can be added or corrected later.</p>',
        f'<p class="formnote">The form works only on the queue\'s report site ({site}, Report tab '
        '"Questions for Yang"); a copy of this page opened anywhere else cannot send it.</p>',
        f'<form id="answers-form" class="answers" method="post" action="/api/answers" accept-charset="utf-8">',
        f'<input type="hidden" name="questions_file" value="{e(questions_name)}">',
        '<div class="field"><label for="answer-name">Your name (optional)</label>'
        f'<input type="text" id="answer-name" name="name" maxlength="{monitor.ANSWER_NAME_MAX_CHARS}" '
        'autocomplete="name"></div>']
    for n, title in titles.items():
        label = f'Question {n}' + (f': {M.inline(title)}' if title else '')
        out.append(f'<div class="field"><label for="answer_{n}">{label}</label>'
            f'<textarea id="answer_{n}" name="answer_{n}" rows="4" '
            f'maxlength="{monitor.ANSWER_MAX_CHARS}"></textarea></div>')
    unreachable = ('Not saved: the report site could not be reached. The form works only on the queue\'s '
        f'report site, {REPORT_SITE} (Report tab "Questions for Yang").')
    out += ['<div class="actions"><button type="submit">Send answers</button>'
        '<span id="answers-status" role="status" aria-live="polite"></span></div>', '</form>',
        FORM_JS.replace('@UNREACHABLE@', json.dumps(unreachable))]
    return '\n'.join(out)


def build_page(md: str, questions_name: str, figdir: Path | None = None) -> str:
    """The questions page: `md` converted (numbered questions with their sub-bullets inside them) and
    the answer form below it, for the questions file named `questions_name`."""
    titles = monitor.question_titles(md)
    if not titles:
        raise ValueError(f'{questions_name} has no numbered questions')
    page = M.convert_markdown(md, figdir, nested_items=True)
    head_end, main_end = '</style>\n</head>', '\n</main>\n'
    if page.count(head_end) != 1 or page.count(main_end) != 1:
        raise RuntimeError('the converted page does not have the layout md_to_artifact.py gives')
    page = page.replace(head_end, FORM_CSS + head_end)
    return page.replace(main_end, '\n' + answer_form(questions_name, titles) + main_end)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('questions', type=Path, nargs='?', default=None,
        help="the questions markdown (default: the newest questions_*.md in the shared directory's "
        f'{monitor.ANSWERS_DIR}/)')
    ap.add_argument('--out', type=Path, default=OUT_HTML, help='the page (default: %(default)s)')
    args = ap.parse_args()
    questions = args.questions or newest_questions(shared_dir() / monitor.ANSWERS_DIR)
    if not monitor.QUESTIONS_FILE.fullmatch(questions.name):
        sys.exit(f'{questions.name}: the monitor takes answers only to files named questions_<YYYY-MM-DD_HHMM>.md')
    page = build_page(questions.read_text(), questions.name, questions.parent)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(page)
    n_q = len(monitor.question_titles(questions.read_text()))
    n_box = len(re.findall(r'<textarea ', page))
    print(f'wrote {args.out} ({len(page) // 1024} KB) from {questions}: {n_q} questions, {n_box} answer boxes')
    probs = M.lint(page)
    if probs:
        print('LINT PROBLEMS:'); print('\n'.join(' - ' + p for p in probs)); sys.exit(1)
    print('lint clean')


if __name__ == '__main__':
    main()
