from ionotec import clean_signal


def test_interval_is_half_open():
    assert clean_signal.is_in_interval(0, (0, 10))
    assert not clean_signal.is_in_interval(10, (0, 10))
