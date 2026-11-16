"""Tests of reports/build_report_artifact.py's hover definitions (the `%%% hover` blocks of md_to_artifact.py),
on scratch markdown and on the report markdown itself; nothing is written."""
import html
import re

import pytest

from reports import build_report_artifact as B

REPORT = """# APL and the re-drawn wiring

%%% hover
re-drawn: Of the wiring: a "randomized" version <of> the hemibrain's.
model: what is simulated
%%%

Nov 3, 2026 · @Tom

## The re-drawn wiring

A paragraph on the re-drawn wiring, wrapped
over two lines, and a [link on re-drawn](#the-re-drawn-wiring); pre-re-drawn and re-drawnish are other words.

| statistic | re-drawn |
|---|---|
| model on re-drawn | 0.5 |

1. **E1, re-drawn wirings (model)** goal.

- a re-drawn bullet
- [ ] a re-drawn task
"""


def _visible(page: str) -> str:
    """The page's text outside tags (attributes, so the hover titles, removed)."""
    return html.unescape(re.sub(r'<[^>]+>', '', page))


def test_hover_defs():
    md, defs = B.hover_defs(REPORT)
    assert defs == {'re-drawn': 'Of the wiring: a "randomized" version <of> the hemibrain\'s.',
        'model': 'what is simulated'}
    assert '%%%' not in md and 'randomized' not in md and 'what is simulated' not in md
    assert md.startswith('# APL and the re-drawn wiring\n\n\nNov 3, 2026')


def test_hover_blocks_rendered():
    page = B.convert(REPORT)
    assert 'span.def {' in page
    # the block is never text in the page: no marker, no definition outside the hover titles
    text = _visible(page)
    assert '%%%' not in text and 'randomized' not in text and 'what is simulated' not in text
    title = html.escape(B.hover_defs(REPORT)[1]['re-drawn'], quote=True)
    span = f'<span class="def" title="{title}">re-drawn</span>'
    model = '<span class="def" title="what is simulated">model</span>'
    # headings (and so the page title) are left alone
    assert '<h1>APL and the re-drawn wiring</h1>' in page and '<title>APL and the re-drawn wiring</title>' in page
    assert '<h2>The re-drawn wiring</h2>' in page
    # paragraphs (joined over lines), link text but not the href, whole words only
    assert (f'<p>A paragraph on the {span} wiring, wrapped over two lines, and a '
        f'<a href="#the-re-drawn-wiring">link on {span}</a>; pre-re-drawn and re-drawnish are other words.</p>') in page
    # table cells, header and body
    assert f'<th>{span}</th>' in page and f'<td class="txt">{model} on {span}</td>' in page
    # numbered, bullet and task list items
    assert f'<li><b>E1, {span} wirings ({model})</b> goal.</li>' in page
    assert f'<li>a {span} bullet</li>' in page and f'<li class="task">&#9744; a {span} task</li>' in page
    assert page.count('class="def"') == 9


def test_hover_block_errors():
    with pytest.raises(ValueError, match='not closed'):
        B.convert('# T\n\n%%% hover\nre-drawn: x\n\ntext\n')
    with pytest.raises(ValueError, match='unsupported block marker'):
        B.convert('# T\n\n%%% html plot.html\n\ntext\n')


def test_report_markdown():
    """The report itself: every hover block removed and its terms wrapped outside headings."""
    md = B.MD.read_text()
    defs = B.hover_defs(md)[1]
    page = B.convert(md)
    text = _visible(page)
    assert '%%%' not in text
    for term, definition in defs.items():
        assert definition not in text
        body = re.sub(r'<h[1-6][^>]*>.*?</h[1-6]>|<title>.*?</title>', '', page, flags=re.S)
        assert f'title="{html.escape(definition, quote=True)}">{term}</span>' in body
