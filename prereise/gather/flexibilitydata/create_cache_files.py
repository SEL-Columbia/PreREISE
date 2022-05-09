import os

from doe import doe_data
from geo import geo_data


def collect_all_raw_data(relative_download_path = "raw"):
    """download all required raw data needed for producing cached files

    param string download_path: folder to store the downloaded file
    :return: (*None*)
    """
    
    file_dir = os.path.abspath(os.path.dirname(__file__))
    download_path = os.path.join(file_dir, relative_download_path)
    
    geo_data.get_census_data(download_path)
    geo_data.get_crosswalk_data(download_path)
    geo_data.get_LSE_region_data(download_path)
    geo_data.get_county_fips_data(download_path)

def create_geo_cache_files(relative_raw_path = "raw", relative_cache_path = "cache"):
    """process downloaded raw files and create cached intermediate files

    param string raw_data_path: folder that contains downloaded raw data
    param string cache_path: folder to store processed cache files
    :return: (*None*)
    """
    
    file_dir = os.path.abspath(os.path.dirname(__file__))
    raw_path = os.path.join(file_dir, relative_raw_path)
    cache_path = os.path.join(file_dir, relative_cache_path)
    
    os.makedirs(cache_path, exist_ok=True)
    geo_data.fips_zip_conversion(raw_path, cache_path)
    geo_data.get_fips_population(raw_path, cache_path)
    geo_data.get_zip_population(raw_path, cache_path)
    geo_data.eiaid_to_zip(raw_path, cache_path)
    geo_data.eiaid_to_fips(raw_path, cache_path)    
