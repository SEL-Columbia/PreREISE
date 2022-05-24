import os
import pickle as pkl

import bus_data

case_path = "case.mat"
cache_path = "..//cache"

# find coordinates of buses in a .mat case
bus_data.get_bus_pos(case_path, cache_path)

with open(os.path.join(cache_path, "bus_pos.pkl"), "rb") as fh:
    bus_pos = pkl.load(fh)

# find the FIPS of all buses and store to cache
bus_data.get_bus_fips(bus_pos, cache_path)

# find the ZIP of all buses and store to cache
bus_data.get_bus_zip(bus_pos, cache_path)

# use the cached file to identify EIA ID of buses and add to bus.csv
bus_csv_path = "bus.csv"
bus_csv_out_path = "bus_with_lse.csv"
bus_data.get_all_bus_eiaid(bus_csv_path, cache_path, bus_csv_out_path)
