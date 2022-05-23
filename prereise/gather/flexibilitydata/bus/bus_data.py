import time

import geopy
import numpy as np
import pandas as pd
import requests
import scipy.io as spio
from geopy.extra.rate_limiter import RateLimiter


def get_bus_pos(case_path):
    """Read a .mat case and extract the lat/lon coordinate of all buses

    param string case_path: path to a .mat case file representing a power network
    :return: (*numpy.ndarray*) -- a list of (busID, lat, lon) stored as a numpy array
    """

    case = spio.loadmat(case_path)
    mpc = case["mpc"]

    bus = mpc["bus"][0][0]
    bus_num = bus.shape[0]
    bus2sub = mpc["bus2sub"][0][0]
    sub = mpc["sub"][0][0]
    subid = mpc["subid"][0][0]

    bus_pos = np.zeros((bus_num, 3))

    # [id, lat, long]
    subid_list = subid[:, 0].tolist()
    for i in range(bus_num):

        bus_sub_idx = subid_list.index(bus2sub[i, 0])
        bus_pos[i, 0] = bus[i, 0]
        bus_pos[i, 1] = sub[bus_sub_idx, 2]
        bus_pos[i, 2] = sub[bus_sub_idx, 3]

    return bus_pos


def get_bus_fips(case_path, cache_path, start_idx=0):
    """Try to get FIPS of each bus in a case mat using FCC AREA API
        Can take hours to run, save to cache file for future use

    param string case_path: path to a .mat case file representing a power network
    param string cache_path: folder to store processed cache files
    param int start_idx: pointer to the index of a bus to start query from
    :return: (*None*)
    """
    bus_pos = get_bus_pos(case_path)
    bus_fips_dict = {
        "busid": bus_pos[:, 0],
        "latitude": bus_pos[:, 1],
        "longitude": bus_pos[:, 2],
        "fips": [0] * bus_num,
    }

    url = "https://geo.fcc.gov/api/census/area"

    for i in range(start_idx, bus_num):
        if i % 1000 == 0:
            with open(os.path.join(file_dir, "../cache/bus_fips_usa.pkl"), "wb") as fh:
                pkl.dump(bus_fips_dict, fh)

        pos = bus_pos[i, :]
        params = {"latitude": pos[1], "longitude": pos[2], "format": "json"}
        r = requests.get(url, params=params)

        assert r.status_code == 200

        res = r.json()
        if res["County"]["FIPS"] is not None:
            bus_fips_dict["fips"][i] = int(res["County"]["FIPS"])
        else:
            bus_fips_dict["fips"][i] = -1

        time.sleep(0.02)

    with open(cache_path + "/bus_fips.pkl", "wb") as fh:
        pkl.dump(bus_fips_dict, fh)


def cleanup_zip(zipdict):
    """Try to cleanup a zip dictionary obtained using online query by converting
        to 5-digit integers. Several possible mis-format are considered

    param dicitonary zipdict: a dictionary containing raw zip-code of buses
    :return: (*dictionary*) -- a dictionary containing 5-digit zip codes
    """
    for i in range(len(zipdict["zip"])):
        try:
            zipdict["zip"][i] = int(zipdict["zip"][i])
        except Exception:
            if ":" in zipdict["zip"][i]:
                zipdict["zip"][i] = int(zipdict["zip"][i].split(":")[0])

            elif "-" in zipdict["zip"][i]:
                zipdict["zip"][i] = int(zipdict["zip"][i].split("-")[0])

            elif "‑" in zipdict["zip"][i]:
                zipdict["zip"][i] = int(zipdict["zip"][i].split("‑")[0])

            elif " " in zipdict["zip"][i]:
                success = 0
                for j in zipdict["zip"][i].split(" "):
                    if len(j) == 5:
                        zipdict["zip"][i] = int(j)
                        success = 1
                if success == 0:
                    zipdict["zip"][i] = -1

            else:
                try:
                    zipdict["zip"][i] = int(zipdict["zip"][i][0:5])
                except Exception:
                    zipdict["zip"][i] = -1

    return zipdict


def get_bus_zip(case_path, cache_path, start_idx=0):
    """Try to get ZIP of each bus in a case mat using geopy
        Can take hours to run, save to cache file for future use

    param string case_path: path to a .mat case file representing a power network
    param string cache_path: folder to store processed cache files
    param int start_idx: pointer to the index of a bus to start query from
    :return: (*None*)
    """
    bus_pos = get_bus_pos(case_path)
    bus_zip_dict = {
        "busid": bus_pos[:, 0],
        "latitude": bus_pos[:, 1],
        "longitude": bus_pos[:, 2],
        "zip": [0] * bus_num,
    }

    geocoder = geopy.Nominatim(user_agent="BES")
    geocode = RateLimiter(
        geocoder.geocode, min_delay_seconds=0.05, return_value_on_exception=None
    )

    def get_zip_code(lat, lon):
        location = geocoder.reverse("{}, {}".format(lat, lon))
        if location is not None:
            address = location.raw["address"]
        else:
            return -1

        if "postcode" in address.keys():
            return address["postcode"]
        else:
            return -1

    for i in range(start_idx, bus_num):
        bus_zip_dict["zip"][i] = int(
            get_zip_code(bus_zip_dict["latitude"][i], bus_zip_dict["longitude"][i])
        )

    bus_zip_dict = cleanup_zip(bus_zip_dict)

    with open(cache_path + "/bus_zip.pkl", "wb") as fh:
        pkl.dump(bus_zip_dict, fh)


def get_all_bus_eiaid(bus_csv_path, cache_path, out_path):
    """Compute the EIA ID of each bus in bus.csv from powersimdata using cached files

    param string bus_csv_path: bus.csv in a powersimdata network model
    param string cache_path: folder to store processed cache files
    param string out_path: output path to store the bus.csv with EIA ID
    :return: (*None*)
    """

    # check all required files
    assert os.path.isfile(
        bus_csv_path
    ), "Incorrect path for network data files bus.csv."
    assert os.path.isfile(
        cache_path + "/bus_fips.pkl"
    ), "Cached file bus_fips.pkl does not exist."
    assert os.path.isfile(
        cache_path + "/bus_zip.pkl"
    ), "Cached file bus_zip.pkl does not exist."
    assert os.path.isfile(
        cache_path + "/eiaid2fips.pkl"
    ), "Cached file eiaid2fips.pkl does not exist."
    assert os.path.isfile(
        cache_path + "/eiaid2zip.pkl"
    ), "Cached file eiaid2zip.pkl does not exist."

    bus_df = pd.read_csv(bus_csv_path)

    with open(cache_path + "/bus_fips.pkl", "rb") as fh:
        bus_fips = pkl.load(fh)

    with open(cache_path + "/bus_zip.pkl", "rb") as fh:
        bus_zip = pkl.load(fh)

    with open(cache_path + "/eiaid2fips.pkl", "rb") as fh:
        eiaid2fips = pkl.load(fh)

    with open(cache_path + "/eiaid2zip.pkl", "rb") as fh:
        eiaid2zip = pkl.load(fh)

    all_eias = list(eia2zip.keys())
    bus_num = bus_df["bus_id"].shape[0]
    bus_df["eia_id"] = np.zeros(bus_num, dtype=int)

    # match with zip
    for i in all_eias:
        e_zips = eia2zip[i]

        # all zips that belong to this eia
        for z in e_zips:
            # find the buses that belong to this zip
            if z in bus_zip["zip"]:
                tar_idxs = [i for i, x in enumerate(bus_zip["zip"]) if x == z]

                # assign the buses to this eia
                for t in tar_idxs:
                    bus_tar = bus_zip["busid"][t]
                    bus_df.loc[bus_df["bus_id"] == bus_tar, ["eia_id"]] = i

    # match again with fips
    for i in all_eias:
        e_fips = eia2fips[i][0]

        # all zips that belong to this eia
        for f in e_fips:
            # find the buses that belong to this zip
            if f in bus_fips["fips"]:
                tar_idxs = [i for i, x in enumerate(bus_fips["fips"]) if x == f]

                # assign the buses to this eia if zip match failed
                for t in tar_idxs:
                    bus_tar = bus_fips["busid"][t]
                    if bus_df.loc[bus_df["bus_id"] == bus_tar, ["eia_id"]].values == 0:
                        bus_df.loc[bus_df["bus_id"] == bus_tar, ["eia_id"]] = i

    bus_df.to_csv(out_path)
