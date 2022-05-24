import os
import pickle as pkl

import numpy as np

from prereise.gather.flexibilitydata.bus.bus_data import get_bus_fips, get_bus_zip


def test_get_bus_fips():
    """Test the FCC Area API using constant dataframe"""
    bus_pos = np.array(
        [[1, 29.7404, -95.3698], [2, 38.8977, -77.0365], [3, 30.6066, -96.3568]]
    )

    # query for fips data, stored in same folder
    get_bus_fips(bus_pos, "")

    # check result
    with open("bus_fips.pkl", "rb") as fh:
        bus_fips = pkl.load(fh)

    assert bus_fips["fips"] == [48201, 11001, 48041]

    # delete file
    os.remove("bus_fips.pkl")


def test_get_bus_zip():
    """Test the geopy OSM query using constant dataframe"""
    bus_pos = np.array(
        [[1, 29.7404, -95.3698], [2, 38.8977, -77.0365], [3, 30.6066, -96.3568]]
    )

    # query for fips data, stored in same folder
    get_bus_zip(bus_pos, "")

    # check result
    with open("bus_zip.pkl", "rb") as fh:
        bus_zip = pkl.load(fh)

    assert bus_zip["zip"] == [77004, 20500, 77845]

    # delete file
    os.remove("bus_zip.pkl")
