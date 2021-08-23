# This script develops puma-level data directly from census and aggregated from census tract data
import os

import geopandas as gpd
import pandas as pd

from prereise.gather.demanddata.bldg_electrification import const


def aggregate_puma_df(puma_fuel_2010, tract_puma_mapping):
    """Scale census tract data up to puma areas
    :param pandas.DataFrame puma_fuel_2010: household fuel type by puma
    :param pandas.DataFrame tract_puma_mapping: tract to puma mapping
    :return: pandas.DataFrame puma_df: population, residential and commercial
    areas, heating degree days, cooling degree days, residential space heating
    fuel usage fractions
    """
    # Set up puma_df data frame
    puma_df = puma_fuel_2010["state"].to_frame()

    # Initialize columns that will be created via summing/averaging
    sum_columns = [
        "pop_2010",
        "res_area_gbs_m2",
        "com_area_gbs_m2",
    ]
    weighted_sum_columns = [
        "hdd65_normals_2010",
        "cdd65_normals_2010",
    ]

    # Load RECS and CBECS area scales for res and com
    resscales = pd.read_csv(os.path.join(data_dir, "area_scale_res.csv"))
    comscales = pd.read_csv(os.path.join(data_dir, "area_scale_com.csv"))

    # Interpolate a 2010 area to scale model area to corresponding RECS/CBECS area
    resscales["2010_scalar"] = (
        resscales["RECS2009"]
        + (resscales["RECS2015"] - resscales["RECS2009"])
        * (
            (const.target_year - const.recs_date_1)
            / (const.recs_date_2 - const.recs_date_1)
        )
    ) / resscales["Model2010"]
    comscales["2010_scalar"] = (
        comscales["CBECS2003"]
        + (comscales["CBECS2012"] - comscales["CBECS2003"])
        * (
            (const.target_year - const.cbecs_date_1)
            / (const.cbecs_date_2 - const.cbecs_date_1)
        )
    ) / comscales["Model2010"]

    # Collect all state data into one data frame
    tract_data = pd.concat(
        [
            pd.read_csv(os.path.join(data_dir, f"tract_data/tract_data_{state}.csv"))
            for state in const.state_list
        ]
    )

    # Rename tract IDs to match all dataframe naming conventions
    tract_data["id"] = ["tract_" + str(i).rjust(11, "0") for i in tract_data["id"]]
    tract_data["id"] = tract_data["id"].astype("str")
    tract_data.set_index("id", inplace=True)

    # Sum population, housing units, and areas
    grouped_tracts = tract_data.groupby(tract_puma_mapping["puma"])
    for col in sum_columns:
        col_to_sum = col.replace("_gbs", "").replace("_", ".")
        puma_df.loc[grouped_tracts.groups.keys(), col] = grouped_tracts[
            col_to_sum
        ].sum()

    # Population-weighted average hdd, cdd, and acpen
    for col in weighted_sum_columns:
        col_to_sum = col.replace("_normals_2010", "").replace("_res_2010", ".res")
        weighted_elements = tract_data[col_to_sum] * tract_data["pop.2010"]
        puma_df[col] = (
            weighted_elements.groupby(tract_puma_mapping["puma"]).sum()
            / tract_data["pop.2010"].groupby(tract_puma_mapping["puma"]).sum()
        )

    # Scale puma area from gbs to 2010 RECS/CBECS
    for state in const.state_list:
        state_row_scale_res = resscales[resscales.eq(state).any(1)].reset_index()
        state_row_scale_com = comscales[comscales.eq(state).any(1)].reset_index()
        res_scalar = state_row_scale_res["2010_scalar"][0]
        com_scalar = state_row_scale_com["2010_scalar"][0]
        puma_df.loc[puma_df["state"] == state, "res_area_2010_m2"] = (
            puma_df[puma_df["state"] == state]["res_area_gbs_m2"] * res_scalar
        )
        puma_df.loc[puma_df["state"] == state, "com_area_2010_m2"] = (
            puma_df[puma_df["state"] == state]["com_area_gbs_m2"] * com_scalar
        )

    # Calculate res fractions of fuel usage based off puma_fuel_2010 household data
    puma_df["frac_sh_res_natgas"] = (
        puma_fuel_2010["hh_utilgas"] / puma_fuel_2010["hh_total"]
    )
    puma_df["frac_sh_res_fok"] = puma_fuel_2010["hh_fok"] / puma_fuel_2010["hh_total"]
    puma_df["frac_sh_res_othergas"] = (
        puma_fuel_2010["hh_othergas"] / puma_fuel_2010["hh_total"]
    )
    puma_df["frac_sh_res_coal"] = puma_fuel_2010["hh_coal"] / puma_fuel_2010["hh_total"]
    puma_df["frac_sh_res_wood"] = puma_fuel_2010["hh_wood"] / puma_fuel_2010["hh_total"]
    puma_df["frac_sh_res_solar"] = (
        puma_fuel_2010["hh_solar"] / puma_fuel_2010["hh_total"]
    )
    puma_df["frac_sh_res_elec"] = puma_fuel_2010["hh_elec"] / puma_fuel_2010["hh_total"]
    puma_df["frac_sh_res_other"] = (
        puma_fuel_2010["hh_other"] / puma_fuel_2010["hh_total"]
    )
    puma_df["frac_sh_res_none"] = puma_fuel_2010["hh_none"] / puma_fuel_2010["hh_total"]

    return puma_df


def scale_fuel_fractions(puma_df, regions, fuel):
    """Scale census tract data up to puma areas
    :param pandas.DataFrame puma_df: output of aggregate_puma_df()
    :param list of lists regions: state regions used to scale fuel fractions
    :param list fuel: types of fuel
    :return: pandas.DataFrame puma_df_frac_ff: fractions of natural gas, fuel
    oil and kerosone, propane, and electricity used for space heating, hot 
    water, cooking, and other in residential and commercial buildings
    """
    
    for c in const.classes:
        if c == "res":
            uselist = ["dhw", "other"]
        else:
            uselist = ["sh", "dhw", "cook"]
        for u in uselist:
            frac_area = pd.DataFrame(columns=fuel)
    
            # Compute frac_area for each fuel type in each region
            for i in regions:
                fuellist = []
                for j in fuel:
                    region_df = puma_df[puma_df["state"].isin(i)].reset_index()
                    fuellist.append(
                        sum(
                            region_df[f"frac_sh_res_{j}"]
                            * region_df[f"{c}_area_2010_m2"]
                        )
                        / sum(region_df[f"{c}_area_2010_m2"])
                    )
                df_i = len(frac_area)
                frac_area.loc[df_i] = fuellist
    
            # Values calculated externally
            frac_scale = pd.read_csv(os.path.join(data_dir, f"frac_target_{u}_{c}.csv"))
    
            downscalar = frac_scale / frac_area
    
            upscalar = (frac_scale - frac_area) / (1 - frac_area)
    
            # Scale frac_hh_fuel to frac_com_fuel
            for f in fuel:
                scalar = 1
                fraccom = []
                for i in range(len(puma_df)):
                    
                    for j in range(len(regions)):
                        if puma_df["state"][i] in regions[j]:
                            region_index = j
                    if downscalar[f][region_index] <= 1:
                        scalar = downscalar[f][region_index]
                        fraccom.append(puma_df[f"frac_sh_res_{f}"][i] * scalar)
                    else:
                        scalar = upscalar[f][region_index]
                        fraccom.append(
                            (1 - puma_df[f"frac_sh_res_{f}"][i]) * scalar
                            + puma_df[f"frac_sh_res_{f}"][i]
                        )
                puma_df[f"frac_{u}_{c}_{f}"] = fraccom
    
    
    # Sum coal, wood, solar and other fractions for frac_com_other
    puma_df["frac_sh_com_other"] = puma_df[
        ["frac_sh_res_coal", "frac_sh_res_wood", "frac_sh_res_solar", "frac_sh_res_other"]
    ].sum(axis=1)
    
    for c in const.classes:
        if c == "res":
            uselist = ["sh", "dhw", "other"]
        else:
            uselist = ["sh", "dhw", "cook"]
        for u in uselist:
            puma_df[f"frac_ff_{u}_{c}_2010"] = puma_df[
                [
                    f"frac_{u}_{c}_natgas",
                    f"frac_{u}_{c}_othergas",
                    f"frac_{u}_{c}_fok",
                ]
            ].sum(axis=1)
            puma_df[f"frac_elec_{u}_{c}_2010"] = puma_df[
                f"frac_{u}_{c}_elec"
            ]
    return puma_df


def puma_timezone_join(timezones, pumas):
    """Assign timezone to puma
    :param geopandas.DataFrame timezones: US timezones
    :param geopandas.DataFrame pumas: US pumas
    :return: pandas.Series puma_timezone["TZID"]: timezone for every puma
    """
    
    puma_timezone = gpd.overlay(pumas, timezones.to_crs("EPSG:4269"))
    
    puma_timezone["area"] = puma_timezone.area
    puma_timezone.sort_values("area", ascending=False, inplace=True)
    puma_timezone = puma_timezone.drop_duplicates(subset="GEOID10", keep="first")
    puma_timezone.sort_values("GEOID10", ascending=True, inplace=True)
    
    return puma_timezone["TZID"]


if __name__ == "__main__":
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

    # Load ACS fuel data for 2010
    puma_fuel_2010 = pd.read_csv(
        os.path.join(data_dir, "puma_fuel_2010.csv"), index_col="puma"
    )

    # Load tract_puma_mapping
    tract_puma_mapping = pd.read_csv(
        os.path.join(data_dir, "tract_puma_mapping.csv"), index_col="tract"
    )

    puma_df = aggregate_puma_df(puma_fuel_2010, tract_puma_mapping)

    puma_df_frac_ff = scale_fuel_fractions(puma_df, const.regions, const.fuel)

    # Add time zone information
    timezones_shp = gpd.GeoDataFrame(gpd.read_file(os.path.join(data_dir, "tz_us.shp")))
    pumas_shp = gpd.GeoDataFrame(gpd.read_file(os.path.join(data_dir, "pumas.shp")))
    puma_timezones = pd.read_csv(
    os.path.join(data_dir, "puma_timezone.csv"), index_col="puma"
    )
    puma_df_frac_ff["timezone"] = puma_timezones["timezone"]

    puma_df_frac_ff.to_csv(os.path.join(data_dir, "puma_data.csv"))
