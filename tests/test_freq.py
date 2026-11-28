import math
import pytest
from ionotec import freq


def test_gps_l1():
    assert freq.getFrequency("G01", "C1C") == pytest.approx(1575.42e6)


def test_glonass_uses_channel_number():
    # R02 is channel -4, L1 step is 9/16 MHz
    assert freq.getFrequency("R02", "C1C") == pytest.approx(1602e6 - 4 * 0.5625e6)


def test_unknown_glonass_sat_is_nan():
    assert math.isnan(freq.getFrequency("R99", "C1C"))


def test_alpha_gps_matches_module_constant():
    assert freq.getAlpha("G01", "C1C", "C2W") == pytest.approx(freq.gps_alpha)


@pytest.mark.parametrize("sv, c1, c2", [
    ("G01", "C1C", "C1W"),   # same frequency twice
    ("X01", "C1C", "C2W"),   # unknown constellation
    ("R99", "C1C", "C2P"),   # GLONASS sat without channel
])
def test_alpha_is_nan_when_undefined(sv, c1, c2):
    assert math.isnan(freq.getAlpha(sv, c1, c2))
