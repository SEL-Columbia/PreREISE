import os
import sys

import pytest

sys.path.append("..")

from doe import doe_data
from geo import geo_data
from geo.batch_process import *


def test_batch_download():
    """Test the downloader from all raw data sources, check if file exist"""

    rel_download_path = "..\\raw"
    collect_all_raw_data(rel_download_path)

    # check downloaded files
    assert os.path.isfile(rel_download_path + "\\county_fips_master.csv")
    assert os.path.isfile(rel_download_path + "\\county_population.csv")
    assert os.path.isfile(rel_download_path + "\\county_to_zip.csv")
    assert os.path.isfile(rel_download_path + "\\iou_zipcodes_2019.csv")
    assert os.path.isfile(rel_download_path + "\\non_iou_zipcodes_2019.csv")


def test_cache_production():
    """Test the functions that produce cached files"""

    rel_raw_path = "..\\raw"
    rel_cache_path = "..\\cache"

    create_geo_cache_files(rel_raw_path, rel_cache_path)

    # check cache files
    assert os.path.isfile(rel_cache_path + "\\eiaid2fips.pkl")
    assert os.path.isfile(rel_cache_path + "\\eiaid2zip.pkl")
    assert os.path.isfile(rel_cache_path + "\\fips2zip.pkl")
    assert os.path.isfile(rel_cache_path + "\\fips_population.pkl")
    assert os.path.isfile(rel_cache_path + "\\zip_population.pkl")
    assert os.path.isfile(rel_cache_path + "\\zip2fips.pkl")
