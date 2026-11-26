"""reports.report_values.format_value: signs, rounding to zero and the specs the templates use."""
import numpy as np
import pytest

from reports.report_values import format_value


@pytest.mark.parametrize('value, spec, expected', [
    (0.0036390849, '+.4f', '+0.0036'),
    (-0.0022, '+.4f', '−0.0022'),
    (-0.00001, '.4f', '0.0000'),
    (-0.00001, '+.4f', '0.0000'),
    (-0.234, '.0%', '−23%'),
    (-0.001, '.0%', '0%'),
    (0.235, '.0%', '24%'),
    (-1234.5, ',.1f', '−1,234.5'),
    (-0.4, ',.0f', '0'),
    (-1.5e-7, '.2e', '−1.50e−07'),
    (-0.0, '.2e', '0.00e+00'),
    (3.7, '.0f', '4'),
    (3.6, 'd', '4'),
    ('-0.0004733129370330814', '+.4f', '−0.0005'),
    ('four', '+.4f', 'four'),
    (2, 'w', 'two'),
    (6.0, 'w', 'six'),
    (np.int64(3), 'w', 'three'),
    (0, 'w', 'zero'),
    (20, 'w', 'twenty'),
    (21, 'w', '21'),
    (-1, 'w', '−1'),
    ('3', 'w', 'three'),
    ('as written', None, 'as written'),
])
def test_format_value(value, spec, expected):
    assert format_value(value, spec) == expected


def test_table_skips_leading_caveat_lines(tmp_path):
    from reports.report_values import _table
    f = tmp_path / 't.csv'
    f.write_text('# a caveat line\n# another\na,b\n1,x#y\n2,z\n')
    df = _table(str(f))
    assert list(df.columns) == ['a', 'b'] and df['b'].tolist() == ['x#y', 'z']
    g = tmp_path / 'u.csv'
    g.write_text('a,b\n1,2\n')
    assert _table(str(g)).shape == (1, 2)


def _shared(monkeypatch, tmp_path):
    from kcq import harness as H
    from reports import report_values as R
    monkeypatch.setattr(H, 'shared_dir', lambda: tmp_path)
    monkeypatch.setattr(R, 'RESULTS', tmp_path / 'results')
    return R


def test_json_key(monkeypatch, tmp_path):
    """`json:<file>:<path>` reads a value, `len(<path>)` a length; a dict or list without len() is refused."""
    import json
    R = _shared(monkeypatch, tmp_path)
    d = tmp_path / 'results' / 'E1_x'
    d.mkdir(parents=True)
    (d / 'lock.json').write_text(json.dumps(dict(a=dict(b=[3, 4], c='inside', **{'six/e134': 'far'}))))
    assert R.resolve('json:results/E1_x/lock.json:a.b.1').value == 4
    assert R.resolve('json:results/E1_x/lock.json:len(a.b)').value == 2
    r = R.resolve('json:results/E1_x/lock.json:a.six/e134')
    assert r.value == 'far' and r.source['kind'] == 'json' and 'lock.json, a.six/e134' in R.tooltip(r)
    with pytest.raises(KeyError):
        R.resolve('json:results/E1_x/lock.json:a.b')
    with pytest.raises(KeyError):
        R.resolve('json:results/E1_x/lock.json:a.missing')


def test_json_key_with_dots_in_a_key(monkeypatch, tmp_path):
    """A dict key that contains dots (E141's variant 'D1.8_inhadd1000') is found by joining the path's parts; a plain
    key still wins, and a path that matches nothing is refused."""
    import json
    R = _shared(monkeypatch, tmp_path)
    d = tmp_path / 'results' / 'E1_x'
    d.mkdir(parents=True)
    (d / 'classes.json').write_text(json.dumps(dict(classes={'D1.8_inhadd1000': dict(cls='holds'),
                                                             'D1': {'8_inhadd1000': dict(cls='plain')},
                                                             'D5.8.x': dict(cls='three parts')})))
    assert R.resolve('json:results/E1_x/classes.json:classes.D1.8_inhadd1000.cls').value == 'plain'
    assert R.resolve('json:results/E1_x/classes.json:classes.D5.8.x.cls').value == 'three parts'
    (d / 'classes.json').write_text(json.dumps(dict(classes={'D1.8_inhadd1000': dict(cls='holds')})))
    assert R.resolve('json:results/E1_x/classes.json:classes.D1.8_inhadd1000.cls').value == 'holds'
    with pytest.raises(KeyError):
        R.resolve('json:results/E1_x/classes.json:classes.D1.9_inhadd1000.cls')


def test_const_key(monkeypatch, tmp_path):
    """`const:<file>:<NAME>[.<index>]` reads a module-level constant without importing the file."""
    from kcq import harness as H
    from reports import report_values as R
    monkeypatch.setattr(H, 'HERE', tmp_path)
    (tmp_path / 'm.py').write_text("import os\nA = (0.045, 0.085)\nB = 0.47 - 0.10\nC: int = 1000\nD = os.sep\n")
    assert R.resolve('const:m.py:A.1').value == 0.085
    assert abs(R.resolve('const:m.py:B').value - 0.37) < 1e-12
    r = R.resolve('const:m.py:C')
    assert r.value == 1000 and r.source['line'] == 4 and format_value(r.value, ',d') == '1,000'
    for bad in ('const:m.py:D', 'const:m.py:E', 'const:missing.py:A'):
        with pytest.raises(KeyError):
            R.resolve(bad)


def test_min_max_keys(monkeypatch, tmp_path):
    """`min:` / `max:` take the extreme over the rows a table key selects; `col!=value` excludes and `*` selects all."""
    R = _shared(monkeypatch, tmp_path)
    d = tmp_path / 'results' / 'E2_y'
    d.mkdir(parents=True)
    (d / 't.csv').write_text('wiring,g,v\nreal,a,0.5\nP1,a,0.2\nP2,b,0.9\nmean,a,0.4\n')
    lo = R.resolve('min:E2/t.csv:g=a;wiring!=mean:v')
    assert lo.value == 0.2 and lo.source['picked'] == dict(wiring='P1')
    assert R.resolve('max:E2/t.csv:g=a;wiring!=mean:v').value == 0.5
    assert R.resolve('max:E2/t.csv:*:v').value == 0.9
    assert 'the lowest of 2 rows' in R.tooltip(lo)
    with pytest.raises(KeyError):
        R.resolve('min:E2:g=a:v')
