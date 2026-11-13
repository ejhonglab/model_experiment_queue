"""reports.report_values.format_value: signs, rounding to zero and the specs the templates use."""
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
    ('as written', None, 'as written'),
])
def test_format_value(value, spec, expected):
    assert format_value(value, spec) == expected
