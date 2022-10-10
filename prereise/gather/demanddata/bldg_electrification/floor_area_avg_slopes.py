import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from prereise.gather.demanddata.bldg_electrification import const
from prereise.gather.demanddata.bldg_electrification.helper import (
    read_shapefile,
    state_shp_overlay,
    zone_shp_overlay,
)


def get_zone_floor_area(iso, zone_shape, pumas_shp):
    """Computes the zone floor area for each ISO.

    :param str iso: abbrev. name of ISO.
    :param geopandas.GeoDataFrame zone_shape: geo data frame of zone(BA) shape file
    :param geopandas.GeoDataFrame pumas_shp: geo data frame of pumas shape file
    :return: (*pandas.DataFrame*) -- Floor area in square meters for all the zones
            with breakdowns of residential, commercial, total heated and total cooled

    .. note:: zone floor area in square meters saved as csv into Profiles/result_stats
    """
    zone_floor_area = pd.DataFrame()
    for zone in zone_names[iso]:
        puma_data_zone = zone_shp_overlay(
            zone_name_shps[iso][zone_names[iso].index(zone)], zone_shape, pumas_shp
        )
        puma_data_zone = puma_data_zone[~(puma_data_zone["frac_in_zone"] < 0.05)]

        total_area_zone_cool_res = (
            (
                puma_data_zone[f"res_area_{base_year}_m2"]
                * puma_data_zone["AC_penetration"]
            )
            * puma_data_zone["frac_in_zone"]
        ).sum()
        total_area_zone_cool_com = (
            (
                puma_data_zone[f"com_area_{base_year}_m2"]
                * puma_data_zone["AC_penetration"]
            )
            * puma_data_zone["frac_in_zone"]
        ).sum()
        total_area_zone_heat_res = (
            (
                puma_data_zone[f"res_area_{base_year}_m2"]
                * puma_data_zone[f"frac_elec_sh_res_{base_year}"]
            )
            * puma_data_zone["frac_in_zone"]
        ).sum()
        total_area_zone_heat_com = (
            (
                puma_data_zone[f"com_area_{base_year}_m2"]
                * puma_data_zone[f"frac_elec_sh_com_{base_year}"]
            )
            * puma_data_zone["frac_in_zone"]
        ).sum()
        total_area_res = (
            puma_data_zone[f"res_area_{base_year}_m2"] * puma_data_zone["frac_in_zone"]
        ).sum()
        total_area_com = (
            puma_data_zone[f"com_area_{base_year}_m2"] * puma_data_zone["frac_in_zone"]
        ).sum()
        zone_floor_area = pd.concat(
            [
                zone_floor_area,
                pd.DataFrame(
                    {
                        "total res": total_area_res,
                        "total com": total_area_com,
                        "heat": total_area_zone_heat_res + total_area_zone_heat_com,
                        "% Heated": (
                            (total_area_zone_heat_res + total_area_zone_heat_com)
                            / (total_area_res + total_area_com)
                        )
                        * 100,
                        "res heat": total_area_zone_heat_res,
                        "com heat": total_area_zone_heat_com,
                        "cool": total_area_zone_cool_res + total_area_zone_cool_com,
                        "res cool": total_area_zone_cool_res,
                        "com cool": total_area_zone_cool_com,
                    },
                    index=[zone],
                ),
            ]
        )

    zone_floor_area.to_csv(
        f"./Profiles/result_stats/{iso_name[iso]}_zone_floor_area_m2.csv"
    )
    return zone_floor_area


def main_plots(iso, zone_shape, pumas_shp, state_shp, country_shp, plot_boolean, size):
    """Creats floor area avraged slopes for all zones within the ISO for one year.

    :param str iso: abbrev. name of ISO.
    :param geopandas.GeoDataFrame zone_shape: geo data frame of zone(BA) shape file
    :param geopandas.GeoDataFrame pumas_shp: geo data frame of pumas shape file
    :param geopandas.GeoDataFrame state_shp: geo data frame of state shape file
    :param geopandas.GeoDataFrame country_shp: geo data frame of nation shape file
    :param boolean plot_boolean: show the plot or not.
    :param int size: defining the image size of plots in dpi.

    .. note:: Floor area averaged heating and cooling slope, error and map plots for
        all zones in each ISO saved as png and csv into Profiles/result_stats/hourly_plots
    """

    zone_floor_area = get_zone_floor_area(iso, zone_shape, pumas_shp)

    # Slope plots for all zones in each ISO in btu/m2/C

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    for zone in zone_names[iso]:
        dayhour_fits = pd.read_csv(
            f"https://besciences.blob.core.windows.net/datasets/bldg_el/dayhour_fits/{zone}_dayhour_fits_{base_year}.csv",
            index_col=0,
        )
        for wk_wknd in ["wk", "wknd"]:
            dayhour_fits[f"s.heat.{wk_wknd}"] = (
                abs(dayhour_fits[f"s.heat.{wk_wknd}"])
                / zone_floor_area.loc[zone, "heat"]
                * const.conv_mw_to_btu
            )
            dayhour_fits[f"s.heat.stderr.{wk_wknd}"] = (
                dayhour_fits[f"s.heat.stderr.{wk_wknd}"]
                / zone_floor_area.loc[zone, "heat"]
                * const.conv_mw_to_btu
            )
            dayhour_fits[f"s.cool.{wk_wknd}.db"] = (
                dayhour_fits[f"s.cool.{wk_wknd}.db"]
                / zone_floor_area.loc[zone, "cool"]
                * const.conv_mw_to_btu
            )
            dayhour_fits[f"s.cool.db.stderr.{wk_wknd}"] = (
                dayhour_fits[f"s.cool.db.stderr.{wk_wknd}"]
                / zone_floor_area.loc[zone, "cool"]
                * const.conv_mw_to_btu
            )

        zone_dayhour = pd.DataFrame()
        zone_dayhour["s.heat"] = (
            dayhour_fits["s.heat.wk"] * 5 + dayhour_fits["s.heat.wknd"] * 2
        ) / 7
        zone_dayhour["s.heat.stderr"] = (
            np.sqrt(
                dayhour_fits["s.heat.stderr.wk"] ** 2 * 5**2
                + dayhour_fits["s.heat.stderr.wknd"] ** 2 * 2**2
            )
            / 7
        )
        zone_dayhour["s.cool.db"] = (
            dayhour_fits["s.cool.wk.db"] * 5 + dayhour_fits["s.cool.wknd.db"] * 2
        ) / 7
        zone_dayhour["s.cool.db.stderr"] = (
            np.sqrt(
                dayhour_fits["s.cool.db.stderr.wk"] ** 2 * 5**2
                + dayhour_fits["s.cool.db.stderr.wknd"] ** 2 * 2**2
            )
            / 7
        )
        ax1.plot(zone_dayhour.index, zone_dayhour["s.heat"], label=zone)
        ax1.set_ylabel("$Btu_e/m^2/^oC$")
        ax1.set_xlabel("hour of day")
        ax1.set_title(f"heating slopes, {iso_name[iso]}")
        ax2.plot(zone_dayhour.index, zone_dayhour["s.cool.db"], label=zone)
        ax2.set_ylabel("$Btu_e/m^2/^oC$")
        ax2.set_xlabel("hour of day")
        ax2.set_title(f"cooling slopes, {iso_name[iso]}")
    ax1.legend()
    ax2.legend()
    if plot_boolean:
        fig1.savefig(
            f"./Profiles/result_stats/hourly_plots/zone_plots/{iso_name[iso]}_heating.png",
            dpi=size,
        )
        fig2.savefig(
            f"./Profiles/result_stats/hourly_plots/zone_plots/{iso_name[iso]}_cooling.png",
            dpi=size,
        )

    # Create hourly slope plots with error bars for each zone within each ISO

    iso_slope = pd.DataFrame()

    iso_dayhour_fits = pd.DataFrame(index=np.arange(0, 24))
    for wk_wknd in ["wk", "wknd"]:

        iso_floor_area = zone_floor_area.loc[zone_names[iso]].sum()

        # read hourly slopes
        dayhour_fits = {
            i: pd.read_csv(
                f"https://besciences.blob.core.windows.net/datasets/bldg_el/dayhour_fits/{zone}_dayhour_fits_{base_year}.csv",
                index_col=0,
            )
            for i, zone in enumerate(zone_names[iso])
        }

        iso_dayhour_fits[f"s.heat.{wk_wknd}"] = (
            np.sum(
                [
                    dayhour_fits[i][f"s.heat.{wk_wknd}"].to_list()
                    for i in range(len(zone_names[iso]))
                ],
                axis=0,
            )
            / iso_floor_area["heat"]
            * const.conv_mw_to_btu
        )

        iso_dayhour_fits[f"s.cool.{wk_wknd}.db"] = (
            np.sum(
                [
                    dayhour_fits[i][f"s.cool.{wk_wknd}.db"].to_list()
                    for i in range(len(zone_names[iso]))
                ],
                axis=0,
            )
            / iso_floor_area["cool"]
            * const.conv_mw_to_btu
        )

        iso_dayhour_fits[f"s.heat.stderr.{wk_wknd}"] = (
            np.sqrt(
                np.sum(
                    [
                        np.square(dayhour_fits[i][f"s.heat.stderr.{wk_wknd}"].to_list())
                        for i in range(len(zone_names[iso]))
                    ],
                    axis=0,
                )
            )
            / iso_floor_area["heat"]
            * const.conv_mw_to_btu
        )

        iso_dayhour_fits[f"s.cool.db.stderr.{wk_wknd}"] = (
            np.sqrt(
                np.sum(
                    [
                        np.square(
                            dayhour_fits[i][f"s.cool.db.stderr.{wk_wknd}"].to_list()
                        )
                        for i in range(len(zone_names[iso]))
                    ],
                    axis=0,
                )
            )
            / iso_floor_area["cool"]
            * const.conv_mw_to_btu
        )
        if plot_boolean:
            fig, ax1 = plt.subplots()
            ax1.errorbar(
                iso_dayhour_fits.index,
                np.abs(iso_dayhour_fits[f"s.heat.{wk_wknd}"]),
                yerr=iso_dayhour_fits[f"s.heat.stderr.{wk_wknd}"] * 1.96,
                fmt="x",
                capsize=3,
                label="Heating Slopes",
            )
            ax1.errorbar(
                iso_dayhour_fits.index + 0.3,
                np.abs(iso_dayhour_fits[f"s.cool.{wk_wknd}.db"]),
                yerr=iso_dayhour_fits[f"s.cool.db.stderr.{wk_wknd}"] * 1.96,
                fmt="x",
                capsize=3,
                color="tab:orange",
                label="Cooling Slopes",
            )
            ax1.set_ylabel("$Btu_e/h/m^2/^oC$", fontsize=14)
            ax1.set_xlabel("Hour of the Day", fontsize=14)
            ax1.legend(fontsize=14)
            text_wk_wknd = "Weekday" if wk_wknd == "wk" else "Weekend"
            plt.title(
                f"{text_wk_wknd} Hourly Profile of Heating and Cooling slopes, {iso_name[iso]}",
                fontsize=12,
            )
            plt.savefig(
                f"./Profiles/result_stats/hourly_plots/iso plot with error/{iso}_{wk_wknd}.png",
                dpi=size,
            )

    iso_dayhour = pd.DataFrame()
    iso_dayhour["s.heat"] = (
        iso_dayhour_fits["s.heat.wk"] * 5 + iso_dayhour_fits["s.heat.wknd"] * 2
    ) / 7
    iso_dayhour["s.heat.stderr"] = (
        np.sqrt(
            iso_dayhour_fits["s.heat.stderr.wk"] ** 2 * 5**2
            + iso_dayhour_fits["s.heat.stderr.wknd"] ** 2 * 2**2
        )
        / 7
    )
    iso_dayhour["s.cool.db"] = (
        iso_dayhour_fits["s.cool.wk.db"] * 5 + iso_dayhour_fits["s.cool.wknd.db"] * 2
    ) / 7
    iso_dayhour["s.cool.db.stderr"] = (
        np.sqrt(
            iso_dayhour_fits["s.cool.db.stderr.wk"] ** 2 * 5**2
            + iso_dayhour_fits["s.cool.db.stderr.wknd"] ** 2 * 2**2
        )
        / 7
    )

    if plot_boolean:
        fig, ax1 = plt.subplots()
        ax1.errorbar(
            iso_dayhour.index,
            np.abs(iso_dayhour["s.heat"]),
            yerr=iso_dayhour["s.heat.stderr"] * 1.96,
            fmt="x",
            capsize=3,
            label="Heating Slopes",
        )
        ax1.errorbar(
            iso_dayhour.index + 0.3,
            np.abs(iso_dayhour["s.cool.db"]),
            yerr=iso_dayhour["s.cool.db.stderr"] * 1.96,
            fmt="x",
            capsize=3,
            color="tab:orange",
            label="Cooling Slopes",
        )
        ax1.set_ylabel("$Btu_e/h/m^2/^oC$", fontsize=14)
        ax1.set_xlabel("Hour of the Day", fontsize=14)
        ax1.legend(fontsize=14)
        plt.title(
            f"Hourly Profile of Heating and Cooling slopes, {iso_name[iso]}",
            fontsize=14,
        )
        plt.savefig(
            f"./Profiles/result_stats/hourly_plots/iso plot with error/{iso}.png",
            dpi=size,
        )

    iso_slope.loc[iso, "s.heat"] = np.mean(iso_dayhour["s.heat"])
    iso_slope.loc[iso, "s.heat.stderr"] = (
        np.sqrt(np.sum(iso_dayhour["s.heat.stderr"] ** 2)) / 24
    )
    iso_slope.loc[iso, "s.cool.db"] = np.mean(iso_dayhour["s.cool.db"])
    iso_slope.loc[iso, "s.cool.db.stderr"] = (
        np.sqrt(np.sum(iso_dayhour["s.cool.db.stderr"] ** 2)) / 24
    )

    iso_slope.to_csv(
        f"./Profiles/result_stats/hourly_plots/iso plot with error/{iso_name[iso]}_iso_slope_with_error.csv"
    )

    # Creating the zone and iso level heating and cooling load in mw/C and btu/m2/C

    zone_slope_df = pd.DataFrame()
    for zone in zone_names[iso]:
        hourly_data = pd.read_csv(
            f"https://besciences.blob.core.windows.net/datasets/bldg_el/dayhour_fits/{zone}_dayhour_fits_{base_year}.csv",
            index_col=0,
        )

        htg_mean = (
            np.mean(hourly_data.loc[:, "s.heat.wk"]) * 5
            + np.mean(hourly_data.loc[:, "s.heat.wknd"]) * 2
        ) / 7
        clg_mean = (
            np.mean(hourly_data.loc[:, "s.cool.wk.db"]) * 5
            + np.mean(hourly_data.loc[:, "s.cool.wknd.db"]) * 2
        ) / 7
        zone_slope_df = pd.concat(
            [
                zone_slope_df,
                pd.DataFrame({"Heating": htg_mean, "Cooling": clg_mean}, index=[zone]),
            ]
        )

    zone_slope_df.to_csv(f"./Profiles/result_stats/{iso_name[iso]}_zone_elec_mw_c.csv")
    zone_elec_btu_m2_c = pd.DataFrame()
    zone_elec_btu_m2_c["Heating"] = (
        zone_slope_df["Heating"] / zone_floor_area["heat"] * const.conv_mw_to_btu
    )
    zone_elec_btu_m2_c["Cooling"] = (
        zone_slope_df["Cooling"] / zone_floor_area["cool"] * const.conv_mw_to_btu
    )
    zone_elec_btu_m2_c.to_csv(
        f"./Profiles/result_stats/{iso_name[iso]}_zone_elec_btu_m2_c.csv"
    )

    # Creating ISO map plots containing load zones with the heating and cooling slope in btu/m2/C
    if (
        iso == "NY"
        or iso == "CA"
        or iso == "PJM"
        or iso == "SPP"
        or iso == "NW"
        or iso == "SE"
    ):
        zone_elec_btu_m2_c[
            (zone_elec_btu_m2_c >= 10) | (zone_elec_btu_m2_c <= -10)
        ] = np.nan
        zone_elec_btu_m2_c.loc["NYIS-ZOND", :] = [np.nan, np.nan]
        zone_elec_btu_m2_c.loc["NYIS-ZONH", :] = [np.nan, np.nan]
        zone_elec_btu_m2_c.loc["WALC", :] = [np.nan, np.nan]
        zone_elec_btu_m2_c.loc["PJM-RECO", :] = [np.nan, np.nan]
        zone_elec_btu_m2_c.loc["SWPP-WFEC", :] = [np.nan, np.nan]
        zone_elec_btu_m2_c.loc["PSEI", :] = [np.nan, np.nan]
        zone_elec_btu_m2_c.loc["AECI", :] = [np.nan, np.nan]

    zone_shp = pd.DataFrame()

    if iso == "NE":
        for i in [
            zone_shp_ma,
            zone_shp_me,
            zone_shp_nh,
            zone_shp_vt,
            zone_shp_ct,
            zone_shp_ri,
        ]:
            zone_shp = zone_shp.append(i)
    elif iso == "PJM":
        for i in [
            zone_shp_de,
            zone_shp_il,
            zone_shp_in,
            zone_shp_ky,
            zone_shp_md,
            zone_shp_nc,
            zone_shp_mi,
            zone_shp_nj,
            zone_shp_oh,
            zone_shp_pa,
            zone_shp_va,
            zone_shp_wva,
            zone_shp_tn,
            zone_shp_dc,
        ]:
            zone_shp = zone_shp.append(i)
    elif iso == "SPP":
        for i in [
            zone_shp_ks,
            zone_shp_ok,
            zone_shp_nm,
            zone_shp_tx,
            zone_shp_ar,
            zone_shp_la,
            zone_shp_mo,
            zone_shp_sd,
            zone_shp_nd,
            zone_shp_mt,
            zone_shp_mn,
            zone_shp_ia,
            zone_shp_wy,
            zone_shp_ne,
        ]:
            zone_shp = zone_shp.append(i)
    elif iso == "MISO":
        for i in [
            zone_shp_la,
            zone_shp_ar,
            zone_shp_ms,
            zone_shp_mi,
            zone_shp_mo,
            zone_shp_ky,
            zone_shp_in,
            zone_shp_il,
            zone_shp_ia,
            zone_shp_mn,
            zone_shp_wi,
            zone_shp_nd,
            zone_shp_sd,
            zone_shp_tx,
            zone_shp_mt,
        ]:
            zone_shp = zone_shp.append(i)
    elif iso == "SW":
        for i in [
            zone_shp_az,
            zone_shp_nm,
            zone_shp_co,
            zone_shp_nv,
            zone_shp_wy,
            zone_shp_sd,
            zone_shp_ne,
        ]:
            zone_shp = zone_shp.append(i)
    elif iso == "NW":
        for i in [
            zone_shp_ca,
            zone_shp_wa,
            zone_shp_or,
            zone_shp_id,
            zone_shp_nv,
            zone_shp_ut,
            zone_shp_wy,
            zone_shp_mt,
        ]:
            zone_shp = zone_shp.append(i)
    elif iso == "SE":
        for i in [
            zone_shp_mo,
            zone_shp_ky,
            zone_shp_ms,
            zone_shp_tn,
            zone_shp_al,
            zone_shp_ga,
            zone_shp_nc,
            zone_shp_sc,
            zone_shp_fl,
            zone_shp_va,
        ]:
            zone_shp = zone_shp.append(i)
    elif iso == "USA" or iso == "Outliers":
        zone_shp = state_shp_overlay("United States", country_shp, zone_shape)

    else:
        zone_shp = state_shp_overlay(iso, state_shp, zone_shape)

    zone_shp.index = zone_shp["BA"]

    for i in range(len(zone_name_shps[iso])):
        zone_shp.loc[zone_name_shps[iso][i], "BA"] = zone_names[iso][i]

    for use in ["Heating", "Cooling"]:
        if iso == "Outliers" and use == "Heating":
            vmax = 250
        elif iso == "Outliers" and use == "Cooling":
            vmax = 100
        else:
            vmax = 10
        zone_shp.index = zone_shp["BA"]
        zone_shp[use] = abs(zone_elec_btu_m2_c[use])
        if plot_boolean:
            fig, ax = plt.subplots(1, 1)
            zone_shp.plot(
                column=use,
                ax=ax,
                vmin=0,
                vmax=vmax,
                legend=True,
                cmap=plt.cm.hot_r,
                edgecolor="black",
                missing_kwds={"color": "lightgrey"},
                legend_kwds={
                    "label": "$Btu_e/m^2/^oC$",
                },
            )
            plt.title(f"Load Zone {use} slopes, {iso_name[iso]}")
            plt.tick_params(
                axis="both",
                which="both",
                bottom=False,
                top=False,
                right=False,
                left=False,
                labelbottom=False,
                labelleft=False,
            )

            plt.savefig(
                f"./Profiles/result_stats/hourly_plots/zone_plots/zone model/{iso}_{use}_zone_model.png",
                dpi=size,
            )


if __name__ == "__main__":

    os.makedirs(
        os.path.join(
            os.path.dirname(__file__),
            "Profiles/result_stats/hourly_plots/zone_plots/zone model",
        ),
        exist_ok=True,
    )

    os.makedirs(
        os.path.join(
            os.path.dirname(__file__),
            "Profiles/result_stats/hourly_plots/iso plot with error",
        ),
        exist_ok=True,
    )

    # Reading Balancing Authority, Pumas, state and country Boundary shapefiles for overlaying
    zone_shape = read_shapefile(
        "https://besciences.blob.core.windows.net/shapefiles/USA/balancing-authorities/ba_area/ba_area.zip"
    )
    pumas_shp = read_shapefile(
        "https://besciences.blob.core.windows.net/shapefiles/USA/pumas-overlay/pumas_overlay.zip"
    )
    state_shp = read_shapefile(
        "https://besciences.blob.core.windows.net/shapefiles/USA/state-outlines/cb_2018_us_state_20m.zip"
    )
    country_shp = read_shapefile(
        "https://besciences.blob.core.windows.net/shapefiles/USA/nation-outlines/cb_2018_us_nation_20m.zip"
    )

    zone_shp_ma = state_shp_overlay("MA", state_shp, zone_shape)
    zone_shp_me = state_shp_overlay("ME", state_shp, zone_shape)
    zone_shp_nh = state_shp_overlay("NH", state_shp, zone_shape)
    zone_shp_vt = state_shp_overlay("VT", state_shp, zone_shape)
    zone_shp_ct = state_shp_overlay("CT", state_shp, zone_shape)
    zone_shp_ri = state_shp_overlay("RI", state_shp, zone_shape)
    zone_shp_de = state_shp_overlay("DE", state_shp, zone_shape)
    zone_shp_il = state_shp_overlay("IL", state_shp, zone_shape)
    zone_shp_in = state_shp_overlay("IN", state_shp, zone_shape)
    zone_shp_ky = state_shp_overlay("KY", state_shp, zone_shape)
    zone_shp_md = state_shp_overlay("MD", state_shp, zone_shape)
    zone_shp_tn = state_shp_overlay("TN", state_shp, zone_shape)
    zone_shp_nc = state_shp_overlay("NC", state_shp, zone_shape)
    zone_shp_mi = state_shp_overlay("MI", state_shp, zone_shape)
    zone_shp_nj = state_shp_overlay("NJ", state_shp, zone_shape)
    zone_shp_oh = state_shp_overlay("OH", state_shp, zone_shape)
    zone_shp_pa = state_shp_overlay("PA", state_shp, zone_shape)
    zone_shp_va = state_shp_overlay("VA", state_shp, zone_shape)
    zone_shp_wva = state_shp_overlay("WV", state_shp, zone_shape)
    zone_shp_dc = state_shp_overlay("DC", state_shp, zone_shape)
    zone_shp_ks = state_shp_overlay("KS", state_shp, zone_shape)
    zone_shp_ok = state_shp_overlay("OK", state_shp, zone_shape)
    zone_shp_nm = state_shp_overlay("NM", state_shp, zone_shape)
    zone_shp_tx = state_shp_overlay("TX", state_shp, zone_shape)
    zone_shp_ar = state_shp_overlay("AR", state_shp, zone_shape)
    zone_shp_la = state_shp_overlay("LA", state_shp, zone_shape)
    zone_shp_mo = state_shp_overlay("MO", state_shp, zone_shape)
    zone_shp_sd = state_shp_overlay("SD", state_shp, zone_shape)
    zone_shp_nd = state_shp_overlay("ND", state_shp, zone_shape)
    zone_shp_mt = state_shp_overlay("MT", state_shp, zone_shape)
    zone_shp_mn = state_shp_overlay("MN", state_shp, zone_shape)
    zone_shp_ia = state_shp_overlay("IA", state_shp, zone_shape)
    zone_shp_wy = state_shp_overlay("WY", state_shp, zone_shape)
    zone_shp_ne = state_shp_overlay("NE", state_shp, zone_shape)
    zone_shp_ms = state_shp_overlay("MS", state_shp, zone_shape)
    zone_shp_wi = state_shp_overlay("WI", state_shp, zone_shape)
    zone_shp_az = state_shp_overlay("AZ", state_shp, zone_shape)
    zone_shp_co = state_shp_overlay("CO", state_shp, zone_shape)
    zone_shp_nv = state_shp_overlay("NV", state_shp, zone_shape)
    zone_shp_ca = state_shp_overlay("CA", state_shp, zone_shape)
    zone_shp_wa = state_shp_overlay("WA", state_shp, zone_shape)
    zone_shp_or = state_shp_overlay("OR", state_shp, zone_shape)
    zone_shp_id = state_shp_overlay("ID", state_shp, zone_shape)
    zone_shp_ut = state_shp_overlay("UT", state_shp, zone_shape)
    zone_shp_al = state_shp_overlay("AL", state_shp, zone_shape)
    zone_shp_ga = state_shp_overlay("GA", state_shp, zone_shape)
    zone_shp_sc = state_shp_overlay("SC", state_shp, zone_shape)
    zone_shp_fl = state_shp_overlay("FL", state_shp, zone_shape)

    zone_names = {
        "NY": [
            "NYIS-ZONA",
            "NYIS-ZONB",
            "NYIS-ZONC",
            "NYIS-ZOND",
            "NYIS-ZONE",
            "NYIS-ZONF",
            "NYIS-ZONG",
            "NYIS-ZONH",
            "NYIS-ZONI",
            "NYIS-ZONJ",
            "NYIS-ZONK",
        ],
        "TX": [
            "ERCO-C",
            "ERCO-E",
            "ERCO-FW",
            "ERCO-N",
            "ERCO-NC",
            "ERCO-S",
            "ERCO-SC",
            "ERCO-W",
        ],
        "CA": [
            "CISO-PGAE",
            "CISO-SCE",
            "CISO-SDGE",
            "IID",
            "WALC",
            "LDWP",
            "TIDC",
            "BANC",
        ],
        "NE": [
            "ISNE-4000",
            "ISNE-4001",
            "ISNE-4002",
            "ISNE-4003",
            "ISNE-4004",
            "ISNE-4005",
        ],
        "PJM": [
            "PJM-AE",
            "PJM-AEP",
            "PJM-AP",
            "PJM-ATSI",
            "PJM-BC",
            "PJM-CE",
            "PJM-DAY",
            "PJM-DEOK",
            "PJM-DUQ",
            "PJM-DPL",
            "PJM-DOM",
            "PJM-EKPC",
            "PJM-JC",
            "PJM-ME",
            "PJM-PE",
            "PJM-PN",
            "PJM-PEP",
            "PJM-PL",
            "PJM-PS",
            "PJM-RECO",
        ],
        "SPP": [
            "SWPP-CSWS",
            "SWPP-EDE",
            "SWPP-GRDA",
            "SWPP-KACY",
            "SWPP-KCPL",
            "SWPP-LES",
            "SWPP-MPS",
            "SWPP-NPPD",
            "SWPP-OKGE",
            "SWPP-OPPD",
            "SWPP-SECI",
            "SWPP-SPRM",
            "SWPP-SPS",
            "SWPP-WAUE",
            "SWPP-WFEC",
            "SWPP-WR",
        ],
        "MISO": [
            "MISO-0001",
            "MISO-0004",
            "MISO-0006",
            "MISO-0027",
            "MISO-0035",
            "MISO-8910",
        ],
        "SW": [
            "PNM",
            "EPE",
            "TEPC",
            "AZPS",
            "WACM",
            "PACE",
            "PSCO",
        ],
        "NW": [
            "PSEI",
            "DOPD",
            "AVA",
            "CHPD",
            "GCPD",
            "BPAT",
            "PGE",
            "PACW",
            "IPCO",
            "NWMT",
            "NEVP",
        ],
        "SE": [
            "AECI",
            "SOCO",
            "AEC",
            "TVA",
            "Carolina",
            "Florida",
        ],
        "USA": [
            "NYIS-ZONA",
            "NYIS-ZONB",
            "NYIS-ZONC",
            "NYIS-ZOND",
            "NYIS-ZONE",
            "NYIS-ZONF",
            "NYIS-ZONG",
            "NYIS-ZONH",
            "NYIS-ZONI",
            "NYIS-ZONJ",
            "NYIS-ZONK",
            "ERCO-C",
            "ERCO-E",
            "ERCO-FW",
            "ERCO-N",
            "ERCO-NC",
            "ERCO-S",
            "ERCO-SC",
            "ERCO-W",
            "CISO-PGAE",
            "CISO-SCE",
            "CISO-SDGE",
            "IID",
            "WALC",
            "LDWP",
            "TIDC",
            "BANC",
            "ISNE-4000",
            "ISNE-4001",
            "ISNE-4002",
            "ISNE-4003",
            "ISNE-4004",
            "ISNE-4005",
            "PJM-AE",
            "PJM-AEP",
            "PJM-AP",
            "PJM-ATSI",
            "PJM-BC",
            "PJM-CE",
            "PJM-DAY",
            "PJM-DEOK",
            "PJM-DUQ",
            "PJM-DPL",
            "PJM-DOM",
            "PJM-EKPC",
            "PJM-JC",
            "PJM-ME",
            "PJM-PE",
            "PJM-PN",
            "PJM-PEP",
            "PJM-PL",
            "PJM-PS",
            "PJM-RECO",
            "SWPP-CSWS",
            "SWPP-EDE",
            "SWPP-GRDA",
            "SWPP-KACY",
            "SWPP-KCPL",
            "SWPP-LES",
            "SWPP-MPS",
            "SWPP-NPPD",
            "SWPP-OKGE",
            "SWPP-OPPD",
            "SWPP-SECI",
            "SWPP-SPRM",
            "SWPP-SPS",
            "SWPP-WAUE",
            "SWPP-WFEC",
            "SWPP-WR",
            "MISO-0001",
            "MISO-0004",
            "MISO-0006",
            "MISO-0027",
            "MISO-0035",
            "MISO-8910",
            "PNM",
            "EPE",
            "TEPC",
            "AZPS",
            "WACM",
            "PACE",
            "PSCO",
            "PSEI",
            "DOPD",
            "AVA",
            "CHPD",
            "GCPD",
            "BPAT",
            "PGE",
            "PACW",
            "IPCO",
            "NWMT",
            "NEVP",
            "AECI",
            "SOCO",
            "AEC",
            "TVA",
            "Carolina",
            "Florida",
        ],
        "Outliers": [
            "NYIS-ZOND",
            "NYIS-ZONH",
            "WALC",
            "PJM-RECO",
            "SWPP-WFEC",
            "PSEI",
            "AECI",
        ],
    }

    zone_name_shps = {
        "NY": [
            "NYISO-A",
            "NYISO-B",
            "NYISO-C",
            "NYISO-D",
            "NYISO-E",
            "NYISO-F",
            "NYISO-G",
            "NYISO-H",
            "NYISO-I",
            "NYISO-J",
            "NYISO-K",
        ],
        "TX": [
            "ERCO-C",
            "ERCO-E",
            "ERCO-FW",
            "ERCO-N",
            "ERCO-NC",
            "ERCO-S",
            "ERCO-SC",
            "ERCO-W",
        ],
        "CA": [
            "CISO-PGAE",
            "CISO-SCE",
            "CISO-SDGE",
            "IID",
            "WALC",
            "LADWP",
            "TID",
            "BANC",
        ],
        "NE": [
            "ISONE-Massachusetts",
            "ISONE-Maine",
            "ISONE-New Hampshire",
            "ISONE-Vermont",
            "ISONE-Connecticut",
            "ISONE-Rhode Island",
        ],
        "PJM": [
            "PJM_AE",
            "PJM_AEP",
            "PJM_AP",
            "PJM_ATSI",
            "PJM_BGE",
            "PJM_ComEd",
            "PJM_DAY",
            "PJM_DEO&K",
            "PJM_DLCO",
            "PJM_DP&L",
            "PJM_Dominion",
            "PJM_EKPC",
            "PJM_JCP&L",
            "PJM_METED",
            "PJM_PECO",
            "PJM_PENELEC",
            "PJM_PEPCO",
            "PJM_PPL",
            "PJM_PSEG",
            "PJM_RECO",
        ],
        "SPP": [
            "SPP-CSWS",
            "SPP-EDE",
            "SPP-GRDA",
            "SPP-KACY",
            "SPP-KCPL",
            "SPP-LES",
            "SPP-MPS",
            "SPP-NPPD",
            "SPP-OKGE",
            "SPP-OPPD",
            "SPP-SECI",
            "SPP-SPRM",
            "SPP-SPS",
            "SPP-WAUE",
            "SPP-WFEC",
            "SPP-WR",
        ],
        "MISO": [
            "MISO-0001",
            "MISO-0004",
            "MISO-0006",
            "MISO-0027",
            "MISO-0035",
            "MISO-8910",
        ],
        "SW": [
            "PNM",
            "EPE",
            "TEPC",
            "Arizona",
            "WACM",
            "PACE",
            "PSCO",
        ],
        "NW": [
            "PSEI",
            "DOPD",
            "AVA",
            "CHPD",
            "GCPD",
            "BPAT",
            "PGE",
            "PACW",
            "IPCO",
            "MT_west",
            "NEVP",
        ],
        "SE": [
            "AECI",
            "SOCO",
            "AEC",
            "TVA",
            "Carolina",
            "Florida",
        ],
        "USA": [
            "NYISO-A",
            "NYISO-B",
            "NYISO-C",
            "NYISO-D",
            "NYISO-E",
            "NYISO-F",
            "NYISO-G",
            "NYISO-H",
            "NYISO-I",
            "NYISO-J",
            "NYISO-K",
            "ERCO-C",
            "ERCO-E",
            "ERCO-FW",
            "ERCO-N",
            "ERCO-NC",
            "ERCO-S",
            "ERCO-SC",
            "ERCO-W",
            "CISO-PGAE",
            "CISO-SCE",
            "CISO-SDGE",
            "IID",
            "WALC",
            "LADWP",
            "TID",
            "BANC",
            "ISONE-Massachusetts",
            "ISONE-Maine",
            "ISONE-New Hampshire",
            "ISONE-Vermont",
            "ISONE-Connecticut",
            "ISONE-Rhode Island",
            "PJM_AE",
            "PJM_AEP",
            "PJM_AP",
            "PJM_ATSI",
            "PJM_BGE",
            "PJM_ComEd",
            "PJM_DAY",
            "PJM_DEO&K",
            "PJM_DLCO",
            "PJM_DP&L",
            "PJM_Dominion",
            "PJM_EKPC",
            "PJM_JCP&L",
            "PJM_METED",
            "PJM_PECO",
            "PJM_PENELEC",
            "PJM_PEPCO",
            "PJM_PPL",
            "PJM_PSEG",
            "PJM_RECO",
            "SPP-CSWS",
            "SPP-EDE",
            "SPP-GRDA",
            "SPP-KACY",
            "SPP-KCPL",
            "SPP-LES",
            "SPP-MPS",
            "SPP-NPPD",
            "SPP-OKGE",
            "SPP-OPPD",
            "SPP-SECI",
            "SPP-SPRM",
            "SPP-SPS",
            "SPP-WAUE",
            "SPP-WFEC",
            "SPP-WR",
            "MISO-0001",
            "MISO-0004",
            "MISO-0006",
            "MISO-0027",
            "MISO-0035",
            "MISO-8910",
            "PNM",
            "EPE",
            "TEPC",
            "Arizona",
            "WACM",
            "PACE",
            "PSCO",
            "PSEI",
            "DOPD",
            "AVA",
            "CHPD",
            "GCPD",
            "BPAT",
            "PGE",
            "PACW",
            "IPCO",
            "MT_west",
            "NEVP",
            "AECI",
            "SOCO",
            "AEC",
            "TVA",
            "Carolina",
            "Florida",
        ],
        "Outliers": [
            "NYISO-D",
            "NYISO-H",
            "WALC",
            "PJM_RECO",
            "SPP-WFEC",
            "PSEI",
            "AECI",
        ],
    }

    iso_name = {
        "NY": "New York",
        "TX": "Texas",
        "CA": "California",
        "NE": "New England",
        "PJM": "PJM Interconnection",
        "SPP": "Southwest Power Pool",
        "MISO": "Midcontinent ISO",
        "SW": "Southwest",
        "NW": "Northwest",
        "SE": "Southeast",
        "USA": "United States",
        "Outliers": "Model Outliers",
    }

    # Use base_year for model results
    base_year = const.base_year

    # If produce profile plots
    plot_boolean = True

    # Plot size in dpi
    size = 700

    for iso in [
        "NY",
        "TX",
        "CA",
        "NE",
        "PJM",
        "SPP",
        "MISO",
        "SW",
        "NW",
        "SE",
        "USA",
        "Outliers",
    ]:
        main_plots(
            iso, zone_shape, pumas_shp, state_shp, country_shp, plot_boolean, size
        )

    # Delete the tmp folder that holds the shapefiles localy after the script is run to completion
    shutil.rmtree(os.path.join("tmp"), ignore_errors=False, onerror=None)
