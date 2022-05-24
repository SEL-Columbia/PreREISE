import os
import urllib.request

import prereise.gather.flexibilitydata.geo.geo_data


def get_cache_from_blob(rel_cache_path="..\\cache"):
    """Download previously uploaded cached files from BLOB storage

    :param str rel_cache_path: folder to store downloaded cache files
    """
    file_dir = os.path.abspath(os.path.dirname(__file__))
    cache_path = os.path.normpath(os.path.join(file_dir, rel_cache_path))

    cache_names = [
        "eiaid2fips.pkl",
        "eiaid2zip.pkl",
        "fips_population.pkl",
        "fips2zip.pkl",
        "zip_population.pkl",
        "zip2fips.pkl",
        "bus_fips.pkl",
        "bus_zip.pkl",
    ]

    blob_path = (
        "https://besciences.blob.core.windows.net/datasets/demand_flexibility_doe/"
    )
    for f in cache_names:
        urllib.request.urlretrieve(blob_path + f, os.path.join(cache_path, f))


def collect_all_raw_data(rel_download_path="..\\raw"):
    """Download all required raw data needed for producing cached files

    :param str rel_download_path: folder to store the downloaded file
    """

    file_dir = os.path.abspath(os.path.dirname(__file__))
    download_path = os.path.normpath(os.path.join(file_dir, rel_download_path))

    geo_data.get_census_data(download_path)
    geo_data.get_crosswalk_data(download_path)
    geo_data.get_lse_region_data(download_path)
    geo_data.get_county_fips_data(download_path)


def create_geo_cache_files(rel_raw_path="..\\raw", rel_cache_path="..\\cache"):
    """Process downloaded raw files and create cached intermediate files

    :param str rel_raw_data_path: folder that contains downloaded raw data
    :param str rel_cache_path: folder to store processed cache files
    """

    file_dir = os.path.abspath(os.path.dirname(__file__))
    raw_path = os.path.normpath(os.path.join(file_dir, rel_raw_path))
    cache_path = os.path.normpath(os.path.join(file_dir, rel_cache_path))

    os.makedirs(cache_path, exist_ok=True)
    geo_data.fips_zip_conversion(raw_path, cache_path)
    geo_data.get_fips_population(raw_path, cache_path)
    geo_data.get_zip_population(raw_path, cache_path)
    geo_data.eiaid_to_zip(raw_path, cache_path)
    geo_data.eiaid_to_fips(raw_path, cache_path)
