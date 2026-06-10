import argparse
import logging
from datetime import datetime, timedelta
import sys
sys.path.append('..')

from tadc import tides as tf


def test_lowest_returns_minimum():
    #Issue #2: Lowest() never updated its running minimum and returned the
    #last tide in the window instead of the lowest
    t0 = datetime(2023, 1, 1)
    dts = [t0, t0 + timedelta(hours=12), t0 + timedelta(hours=24)]
    vals = [-1.0, -0.2, -0.5]
    idx = tf.Lowest(dts, vals, t0 - timedelta(hours=1), t0 + timedelta(hours=25))
    assert idx == 0, "Lowest() did not return the index of the lowest value"
    #Highest() was already correct - make sure it stays that way
    hidx = tf.Highest(dts, [1.0, 0.2, 0.5], t0 - timedelta(hours=1), t0 + timedelta(hours=25))
    assert hidx == 0, "Highest() did not return the index of the highest value"


def test_nearest_tide_before_first_tide():
    #Issue #9: a query time before the first tide returned -1, which callers
    #use to index the LAST tide in the record
    t0 = datetime(2023, 1, 10)
    tide_times = [t0 + timedelta(hours=12 * k) for k in range(20)]
    idx = tf.Nearest_Tide(tide_times, t0 - timedelta(hours=3))
    assert idx == 0, "Nearest_Tide() did not return the first tide for an early query"


def test_nearest_tide_normal_picks():
    #Nearest_Tide() should keep returning the closest tide for in-range queries
    t0 = datetime(2023, 1, 10)
    tide_times = [t0 + timedelta(hours=12 * k) for k in range(20)]
    assert tf.Nearest_Tide(tide_times, t0 + timedelta(hours=5)) == 0
    assert tf.Nearest_Tide(tide_times, t0 + timedelta(hours=7)) == 1
    assert tf.Nearest_Tide(tide_times, t0 + timedelta(hours=300)) == 19


def test_calc_expected_diff_constant_offset():
    #Issue #4: a constant control offset puts every pair exactly at the mean,
    #both buckets came up empty and the error path crashed (NameError/exit)
    t0 = datetime(2023, 1, 1)
    types = ['HH', 'LL', 'H', 'L'] * 10
    HL_Sub = []
    HL_Con = []
    for k in range(40):
        ts = t0 + timedelta(hours=6.21 * k)
        HL_Sub.append([ts, 1.0, types[k]])
        HL_Con.append([ts + timedelta(minutes=30), 1.2, types[k]])
    etd = tf.Calc_Expected_Diff(HL_Sub, HL_Con)
    assert etd == 30, "Expected the constant 30-minute offset to be recovered"


def test_calc_expected_diff_early_control_tides():
    #Issue #4: a few early control tides pull the mean below all regular pairs,
    #the below-mean bucket is legitimately empty and the function used to call
    #exit(-1), killing the host process
    t0 = datetime(2023, 1, 1)
    types = ['HH', 'LL', 'H', 'L'] * 10
    HL_Sub = []
    HL_Con = []
    n_high = -1
    for k in range(40):
        ts = t0 + timedelta(hours=6.21 * k)
        shift = 30
        if types[k][0] == 'H':
            n_high = n_high + 1
            if n_high in (3, 7):
                shift = 30 - 150  #this control tide is 2.5 h early
        HL_Sub.append([ts, 1.0, types[k]])
        HL_Con.append([ts + timedelta(minutes=shift), 1.2, types[k]])
    etd = tf.Calc_Expected_Diff(HL_Sub, HL_Con)
    assert isinstance(etd, int), "Expected an integer time difference, not an abort"


def test_calc_expected_diff_no_pairs_raises_runtimeerror():
    #Issue #4: with no pairable tides at all the function should raise a
    #catchable RuntimeError instead of exiting or dividing by zero
    t0 = datetime(2023, 1, 1)
    HL_Sub = [[t0, 1.0, 'H'], [t0 + timedelta(hours=6.21), -1.0, 'L']]
    HL_Con = [[t0 + timedelta(minutes=900), 1.2, 'H'],
              [t0 + timedelta(hours=6.21) + timedelta(minutes=900), -0.8, 'L']]
    broke = False
    try:
        tf.Calc_Expected_Diff(HL_Sub, HL_Con)
    except RuntimeError:
        broke = True
    assert broke is True, "Expected RuntimeError for unpairable tides"


def test_check_tide_order_warning_formats(caplog):
    #Issue #10: the out-of-order warning passed the datetime as a printf
    #argument with no placeholder, so the logging call itself errored and the
    #diagnostic never reached the user
    t0 = datetime(2023, 1, 1)
    dt = [t0 + timedelta(hours=3 * k) for k in range(6)]
    with caplog.at_level(logging.WARNING, logger='tadc.tides'):
        result = tf.Check_Tide_Order(dt, [0, 2], [4])
    assert result == -1, "Out-of-order tides were not detected"
    messages = [record.getMessage() for record in caplog.records]
    assert 'Tides are out of order at: 2023-01-01 06:00:00' in messages


def test_get_gmt_offset_multidigit():
    #Issue #7: trailing digits were collected walking backwards and appended in
    #that order, reversing multi-digit offsets ('UST10' -> 1, 'UST12' -> 21)
    from tadc.run import Get_GMT_Offset
    assert Get_GMT_Offset('GMT') == 0
    assert Get_GMT_Offset('UST5') == 5
    assert Get_GMT_Offset('UST10') == 10
    assert Get_GMT_Offset('UST11') == 11
    assert Get_GMT_Offset('UST12') == 12


def test_str2bool_make_plots():
    #Issue #10: --make_plots used type=bool, and bool('False') is True, so any
    #value on the command line enabled plotting
    from tadc.run import str2bool
    parser = argparse.ArgumentParser()
    parser.add_argument('--make_plots', type=str2bool, default=False)
    assert parser.parse_args([]).make_plots is False
    assert parser.parse_args(['--make_plots', 'True']).make_plots is True
    assert parser.parse_args(['--make_plots', 'False']).make_plots is False
    assert parser.parse_args(['--make_plots', '0']).make_plots is False
    assert str2bool('yes') is True
    assert str2bool('no') is False
