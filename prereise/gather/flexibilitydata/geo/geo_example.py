import batch_process

download_path = "..\\raw"
cache_path = "..\\cache"

# download raw data from sources
batch_process.collect_all_raw_data(download_path)

# create cache files after downloading raw data
batch_process.create_geo_cache_files(download_path, cache_path)

# download cache files directly from BLOB storage
batch_process.get_cache_from_blob(cache_path)
