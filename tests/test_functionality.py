import os
import pandas as pd
import pytest
import subprocess
import sys
sys.path.append('..')

import tadc

          
def test_simple_timeseries():
    fname = 'data/test_timeseries_simple.csv'
    lat = 32.78
    lon = -79.92
    broke = False
    try:
        out = tadc.run(fname=fname, Subordinate_Lat=lat, Subordinate_Lon=lon)
        assert '***Error***' not in out.readme, "TADC failed silently"
    except (ValueError, RuntimeError):
        broke = True
    assert broke is False , "TADC broke"


def test_only_highs_timeseries():
    fname = 'data/test_timeseries_onlyhighs.csv'
    lat = 32.78
    lon = -79.92
    broke = False
    try:
        out = tadc.run(fname=fname, Subordinate_Lat=lat, Subordinate_Lon=lon)
        assert '***Error***' not in out.readme, "TADC failed silently"
    except (ValueError, RuntimeError):
        broke = True
    assert broke is False , "TADC broke"


def test_small_signal_with_flatline_timeseries():
    fname = 'data/test_timeseries_smallsignal_flatline.csv'
    lat = 32.78
    lon = -79.92
    broke = False
    try:
        out = tadc.run(fname=fname, Subordinate_Lat=lat, Subordinate_Lon=lon)
        assert '***Error***' not in out.readme, "TADC failed silently"
    except (ValueError, RuntimeError):
        broke = True
    assert broke is False , "TADC broke"


def test_big_gap_timeseries():
    fname = 'data/test_timeseries_biggap.csv'
    lat = 32.78
    lon = -79.92
    broke = False
    try:
        out = tadc.run(fname=fname, Subordinate_Lat=lat, Subordinate_Lon=lon)
        assert '***Error***' not in out.readme, "TADC failed silently"
    except (ValueError, RuntimeError):
        broke = True
    assert broke is False , "TADC broke"


def test_uneven_time_sampling_timeseries():
    fname = 'data/test_timeseries_uneventimesampling.csv'
    lat = 32.78
    lon = -79.92
    broke = False
    try:
        out = tadc.run(fname=fname, Subordinate_Lat=lat, Subordinate_Lon=lon)
        assert '***Error***' not in out.readme, "TADC failed silently"
    except (ValueError, RuntimeError):
        broke = True
    assert broke is False , "TADC broke"


def test_nontidal_with_spikes_timeseries():
    fname = 'data/test_timeseries_nontidal_spikes.csv'
    lat = 32.78
    lon = -79.92
    broke = False
    try:
        out = tadc.run(fname=fname, Subordinate_Lat=lat, Subordinate_Lon=lon)
        assert '***Error***' not in out.readme, "TADC failed silently"
    except (ValueError, RuntimeError):
        broke = True
    assert broke is False , "TADC broke"
