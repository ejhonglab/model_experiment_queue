"""Tests of reports/build_questions_page.py on a scratch questions file, with a scratch shared directory and every
output in tmp_path, never in the shared directory or in outputs/_artifacts."""
import os
import re
import subprocess
import sys
from pathlib import Path

from coordination.monitor import server as monitor
from reports import build_questions_page as B
from reports import md_to_artifact as M

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / 'reports' / 'build_questions_page.py'

QUESTIONS = """# Questions for Yang

Asked 2026-11-15 (file `questions_for_yang/questions_2026-11-15_1614.md`).

Answer any of them.

## The APL-only recordings

1. **Use.** May the queue use your files? The sets are:
   - (a) January 2025;
   - (b) 2022-10-31.

   Which should be left out?
2. **Channels.** Is `_C1` always ChanA?

## The calyx table

10. **roi_peak_dff.** How is it computed?
    - the ROI's own trace;
    - or the map's mean.
11. A question without a bold title.
"""


def test_build_page():
    page = B.build_page(QUESTIONS, 'questions_2026-11-15_1614.md')
    assert M.lint(page) == []
    # the questions: one list per section, sub-bullets and the indented paragraph inside their question
    assert re.findall(r'<ol start="(\d+)">', page) == ['1', '10']
    assert ('<li><b>Use.</b> May the queue use your files? The sets are:<ul><li>(a) January 2025;</li>'
        '<li>(b) 2022-10-31.</li></ul>\n<p>Which should be left out?</p></li><li><b>Channels.</b>') in page
    assert '<li><b>roi_peak_dff.</b> How is it computed?<ul><li>the ROI' in page
    # the form below the questions: name, one box per question labelled with its title, the hidden file name
    form = page[page.index('<form'):page.index('</form>')]
    assert page.index('</ol>', page.index('<ol start="10">')) < page.index('<h2 id="answers">') < page.index('<form')
    assert 'method="post" action="/api/answers"' in form
    assert '<input type="hidden" name="questions_file" value="questions_2026-11-15_1614.md">' in form
    assert 'name="name"' in form
    assert re.findall(r'<label for="answer_\d+">(.*?)</label>', form) == [
        'Question 1: Use', 'Question 2: Channels', 'Question 10: roi_peak_dff', 'Question 11']
    assert re.findall(r'<textarea id="answer_(\d+)" name="answer_\1"', form) == ['1', '2', '10', '11']
    assert 'maxlength="{}"'.format(monitor.ANSWER_MAX_CHARS) in form
    # self-contained: no external script or stylesheet; the notes on the report site and on resending
    assert not re.search(r'<script[^>]+src=|<link ', page)
    assert 'works only on the queue\'s report site' in page and 'more than once' in page
    # the boxes' numbers are the ones the monitor accepts
    assert set(int(n) for n in re.findall(r'name="answer_(\d+)"', form)) == set(monitor.question_titles(QUESTIONS))


def test_cli_newest_questions(tmp_path):
    shared = tmp_path / 'shared'
    q = shared / 'questions_for_yang'
    q.mkdir(parents=True)
    (q / 'questions_2026-11-15_1614.md').write_text(QUESTIONS)
    (q / 'questions_2026-11-16_0900.md').write_text(QUESTIONS.replace('**Use.**', '**Later.**'))
    (q / 'answers_2026-11-15_1614__2026-11-15_164012.txt').write_text('not a questions file')
    (q / 'questions_draft.md').write_text('1. **Not this one.** x\n')
    out = tmp_path / 'html' / 'index.html'
    r = subprocess.run([sys.executable, str(SCRIPT), '--out', str(out)], cwd=str(REPO), capture_output=True,
        text=True, env=dict(os.environ, KC_REMAPPING_SHARED=str(shared)))
    assert r.returncode == 0, r.stdout + r.stderr
    assert '4 questions, 4 answer boxes' in r.stdout and 'lint clean' in r.stdout
    page = out.read_text()
    assert 'value="questions_2026-11-16_0900.md"' in page and 'Question 1: Later' in page
    assert sorted(p.name for p in q.iterdir()) == ['answers_2026-11-15_1614__2026-11-15_164012.txt',
        'questions_2026-11-15_1614.md', 'questions_2026-11-16_0900.md', 'questions_draft.md']
    # a file whose name the monitor would refuse is refused here too
    bad = tmp_path / 'my_questions.md'
    bad.write_text(QUESTIONS)
    r = subprocess.run([sys.executable, str(SCRIPT), str(bad), '--out', str(out)], cwd=str(REPO),
        capture_output=True, text=True, env=dict(os.environ, KC_REMAPPING_SHARED=str(shared)))
    assert r.returncode != 0 and 'questions_<YYYY-MM-DD_HHMM>.md' in r.stderr
