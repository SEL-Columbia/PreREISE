import os

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from zone_profile_generator import zone_shp_overlay


def get_zone_floor_area(iso_list):
    """Computes the zone floor area for each ISO.

    :param list of str iso_list: list of abbrev. names of ISO.

    :return: (*pandas.DataFrame*) zone_floor_area -- Floor area for all the zones with breakdowns of res, com and total heating and total cooling
    :return: csv -- Saves the csv of the floor area into result_stats
    """
    zone_floor_area = pd.DataFrame()
    for iso in iso_list:
        for zone in zone_names[iso]:
            puma_data_zone = zone_shp_overlay(
                zone_name_shps[iso][zone_names[iso].index(zone)]
            )
            puma_data_zone = puma_data_zone[~(puma_data_zone["frac_in_zone"] < 0.05)]

            total_area_zone_cool_res = (
                (puma_data_zone["res_area_2019_m2"] * puma_data_zone["AC_penetration"])
                * puma_data_zone["frac_in_zone"]
            ).sum()
            total_area_zone_cool_com = (
                (puma_data_zone["com_area_2019_m2"] * puma_data_zone["AC_penetration"])
                * puma_data_zone["frac_in_zone"]
            ).sum()
            total_area_zone_heat_res = (
                (
                    puma_data_zone["res_area_2019_m2"]
                    * puma_data_zone["frac_elec_sh_res_2019"]
                )
                * puma_data_zone["frac_in_zone"]
            ).sum()
            total_area_zone_heat_com = (
                (
                    puma_data_zone["com_area_2019_m2"]
                    * puma_data_zone["frac_elec_sh_com_2019"]
                )
                * puma_data_zone["frac_in_zone"]
            ).sum()
            total_area_res = (
                puma_data_zone["res_area_2019_m2"] * puma_data_zone["frac_in_zone"]
            ).sum()
            total_area_com = (
                puma_data_zone["com_area_2019_m2"] * puma_data_zone["frac_in_zone"]
            ).sum()
            zone_floor_area = pd.concat(
                [
                    zone_floor_area,
                    pd.DataFrame(
                        {
                            "total res": total_area_res,
                            "total com": total_area_com,
                            "heat": total_area_zone_heat_res + total_area_zone_heat_com,
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
    zone_floor_area.to_csv("./Profiles/result_stats/zone_floor_area.csv")
    return zone_floor_area


def state_shp_overlay(state):
    """Select load zones within a state

    :param str state: abbrev. of state

    :return: (*geopandas.GeoDataFrame*) state boundaries and load zones within it
    """
    state_shapefile = gpd.GeoDataFrame(
        gpd.read_file(
            os.path.join(
                os.path.dirname(__file__), "shapefiles", "cb_2018_us_state_20m.shp"
            )
        )
    )
    state_shp = state_shapefile[state_shapefile["STUSPS"] == state]
    zone_shp = gpd.GeoDataFrame(
        gpd.read_file(
            os.path.join(os.path.dirname(__file__), "shapefiles", "ba_area.shp")
        )
    ).to_crs("EPSG:4269")
    zone_shp["area"] = zone_shp["geometry"].to_crs({"proj": "cea"}).area

    zone_state = gpd.overlay(zone_shp, state_shp.to_crs("EPSG:4269"))
    zone_state["area"] = zone_state["geometry"].to_crs({"proj": "cea"}).area
    zone_state["area_frac"] = [
        zone_state["area"][i]
        / list(zone_shp[zone_shp["BA"] == zone_state["BA"][i]]["area"])[0]
        for i in range(len(zone_state))
    ]
    zone_state.loc[zone_state["area_frac"] >= 0.99, "area_frac"] = 1
    zone_state = zone_state.drop(zone_state[zone_state["area_frac"] <= 0.00001].index)

    return zone_state


def main_plots(iso_list):
    """Creats floor area avraged slopes for all zones within the ISO for one year.

    :param list of str iso_list: list of abbrev. names of ISO.

    :return: Errors plots and Map plots of floor area averaged heating and cooling slopes for all zones in each ISO
    """

    zone_floor_area = get_zone_floor_area(iso_list)

    # Create hourly slope plots with error bars for each zone within each ISO

    iso_slope = pd.DataFrame()

    for iso in iso_list:
        iso_dayhour_fits = pd.DataFrame(index=np.arange(0, 24))
        for wk_wknd in ["wk", "wknd"]:
            iso_floor_area = zone_floor_area[
                zone_floor_area.index.isin(zone_names[iso])
            ].sum(axis=0)
            # read hourly slopes
            dayhour_fits = {}
            i = 0
            for zone in zone_names[iso]:
                dayhour_fits[i] = pd.read_csv(
                    f"./dayhour_fits/{zone}_dayhour_fits_2019.csv", index_col=0
                )
                i += 1

            iso_dayhour_fits[f"s.heat.{wk_wknd}"] = (
                np.sum(
                    [
                        dayhour_fits[i][f"s.heat.{wk_wknd}"].to_list()
                        for i in range(len(zone_names[iso]))
                    ],
                    axis=0,
                )
                / iso_floor_area["heat"]
                * 3412141
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
                * 3412141
            )

            iso_dayhour_fits[f"s.heat.stderr.{wk_wknd}"] = (
                np.sqrt(
                    np.sum(
                        [
                            np.square(
                                dayhour_fits[i][f"s.heat.stderr.{wk_wknd}"].to_list()
                            )
                            for i in range(len(zone_names[iso]))
                        ],
                        axis=0,
                    )
                )
                / iso_floor_area["heat"]
                * 3412141
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
                * 3412141
            )

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
                dpi=200,
            )

        iso_dayhour = pd.DataFrame()
        iso_dayhour["s.heat"] = (
            iso_dayhour_fits[f"s.heat.wk"] * 5 + iso_dayhour_fits[f"s.heat.wknd"] * 2
        ) / 7
        iso_dayhour["s.heat.stderr"] = (
            np.sqrt(
                iso_dayhour_fits[f"s.heat.stderr.wk"] ** 2 * 5**2
                + iso_dayhour_fits[f"s.heat.stderr.wknd"] ** 2 * 2**2
            )
            / 7
        )
        iso_dayhour["s.cool.db"] = (
            iso_dayhour_fits["s.cool.wk.db"] * 5
            + iso_dayhour_fits["s.cool.wknd.db"] * 2
        ) / 7
        iso_dayhour["s.cool.db.stderr"] = (
            np.sqrt(
                iso_dayhour_fits[f"s.cool.db.stderr.wk"] ** 2 * 5**2
                + iso_dayhour_fits[f"s.cool.db.stderr.wknd"] ** 2 * 2**2
            )
            / 7
        )

        fig, ax1 = plt.subplots()
        ax1.errorbar(
            iso_dayhour.index,
            np.abs(iso_dayhour[f"s.heat"]),
            yerr=iso_dayhour[f"s.heat.stderr"] * 1.96,
            fmt="x",
            capsize=3,
            label="Heating Slopes",
        )
        ax1.errorbar(
            iso_dayhour.index + 0.3,
            np.abs(iso_dayhour[f"s.cool.db"]),
            yerr=iso_dayhour[f"s.cool.db.stderr"] * 1.96,
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
            dpi=200,
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
        "./Profiles/result_stats/hourly_plots/iso plot with error/iso_slope_with_error.csv"
    )

    # Creating the zone and iso level heating and cooling load in mw/C and btu/m2/C

    zone_slope_df = pd.DataFrame()
    for iso in iso_list:
        for zone in zone_names[iso]:
            hourly_data = pd.read_csv(
                f"./dayhour_fits/{zone}_dayhour_fits_2019.csv", index_col=0
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
                    pd.DataFrame(
                        {"Heating": htg_mean, "Cooling": clg_mean}, index=[zone]
                    ),
                ]
            )

    zone_slope_df.to_csv("./Profiles/result_stats/zone_elec_mw_c.csv")
    zone_elec_btu_m2_c = pd.DataFrame()
    zone_elec_btu_m2_c["Heating"] = (
        zone_slope_df["Heating"] / zone_floor_area["heat"] * 3412141
    )
    zone_elec_btu_m2_c["Cooling"] = (
        zone_slope_df["Cooling"] / zone_floor_area["cool"] * 3412141
    )
    zone_elec_btu_m2_c.to_csv("./Profiles/result_stats/zone_elec_btu_m2_c.csv")

    # Creating ISO map plots containing load zones with the heating and cooling slope in btu/m2/C

    zone_elec_btu_m2_c[(zone_elec_btu_m2_c >= 9) | (zone_elec_btu_m2_c <= -9)] = np.nan
    zone_elec_btu_m2_c.loc["NYIS-ZOND", :] = [np.nan, np.nan]
    zone_elec_btu_m2_c.loc["NYIS-ZONH", :] = [np.nan, np.nan]
    zone_elec_btu_m2_c.loc["WALC", :] = [np.nan, np.nan]

    for iso in iso_list:
        if iso == "NE":
            zone_shp = pd.DataFrame()
            zone_shp_ma = state_shp_overlay("MA")
            zone_shp_me = state_shp_overlay("ME")
            zone_shp_nh = state_shp_overlay("NH")
            zone_shp_vt = state_shp_overlay("VT")
            zone_shp_ct = state_shp_overlay("CT")
            zone_shp_ri = state_shp_overlay("RI")
            for i in [
                zone_shp_ma,
                zone_shp_me,
                zone_shp_nh,
                zone_shp_vt,
                zone_shp_ct,
                zone_shp_ri,
            ]:
                zone_shp = zone_shp.append(i)
        else:
            zone_shp = state_shp_overlay(iso)

        zone_shp.index = zone_shp["BA"]

        for i in range(len(zone_name_shps[iso])):
            zone_shp.loc[zone_name_shps[iso][i], "BA"] = zone_names[iso][i]

        for use in ["Heating", "Cooling"]:
            zone_shp.index = zone_shp["BA"]
            zone_shp[use] = abs(zone_elec_btu_m2_c[use])
            vmax = np.max(np.max(zone_elec_btu_m2_c.loc[zone_names[iso]])) + 1

            fig, ax = plt.subplots(1, 1)
            heatmap = zone_shp.plot(
                column=use,
                ax=ax,
                vmin=0,
                vmax=6,
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
                dpi=200,
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
    }

    iso_name = {
        "NY": "New York",
        "TX": "Texas",
        "CA": "California",
        "NE": "New England",
    }

    main_plots(["NY", "TX", "CA", "NE"])
