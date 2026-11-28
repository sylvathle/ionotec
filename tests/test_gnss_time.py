from datetime import datetime, timedelta
import pandas as pd
from ionotec import gnss
import pytest


def test_gps_week_zero_is_gps_epoch():
    assert gnss.getDateFromSat(0, 0) == datetime(1980, 1, 6)


def test_beidou_uses_its_own_epoch():
    assert gnss.getDateFromSat(1, 60, system="C") == datetime(2006, 1, 8, 0, 1)


def test_remove_outsiders_respects_thres():
    idx = pd.to_datetime(["2024-05-09 00:00", "2024-05-09 00:30",
                          "2024-05-09 01:00", "2024-05-09 05:00"])
    df = pd.DataFrame({"x": range(4)}, index=idx)

    out = gnss.removeOutsiders(df, thres=timedelta(hours=1))

    assert list(out["x"]) == [0, 1, 2]


def test_remove_outsiders_does_not_modify_input():
    idx = pd.date_range("2024-05-09", periods=3, freq="30min")
    df = pd.DataFrame({"x": range(3)}, index=idx)
    gnss.removeOutsiders(df)
    assert list(df.columns) == ["x"]


H = timedelta(hours=1)

def _df(times):
    return pd.DataFrame({"x": range(len(times))}, index=pd.to_datetime(times))


def test_remove_outsiders_all_gaps_too_big_returns_empty():
    df = _df(["2024-05-09 00:00", "2024-05-09 03:00", "2024-05-09 06:00"])
    assert gnss.removeOutsiders(df, thres=H).empty


@pytest.mark.parametrize("times", [[], ["2024-05-09 00:00"]])
def test_remove_outsiders_short_input_is_returned_unchanged(times):
    df = _df(times)
    pd.testing.assert_frame_equal(gnss.removeOutsiders(df, thres=H), df)


def test_remove_outsiders_drops_isolated_first_point():
    df = _df(["2024-05-09 00:00", "2024-05-09 05:00", "2024-05-09 05:30"])
    assert list(gnss.removeOutsiders(df, thres=H)["x"]) == [1, 2]
