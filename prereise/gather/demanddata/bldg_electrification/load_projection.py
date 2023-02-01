import os
from itertools import product

import numpy as np
import pandas as pd
from pandas.tseries.holiday import USFederalHolidayCalendar as calendar  # noqa: N813

from prereise.gather.demanddata.bldg_electrification import const
from prereise.gather.demanddata.bldg_electrification.ff2elec_profile_generator_cook import (
    generate_cook_profiles,
)
from prereise.gather.demanddata.bldg_electrification.ff2elec_profile_generator_dhw import (
    generate_dhw_profiles,
)
from prereise.gather.demanddata.bldg_electrification.ff2elec_profile_generator_htg import (
    generate_htg_profiles,
)
from prereise.gather.demanddata.bldg_electrification.helper import (
    read_shapefile,
    zone_shp_overlay,
)
from prereise.gather.demanddata.bldg_electrification.zone_profile_generator import (
    zonal_data,
)


class scenarios:
    # read and process a scenario (reflects either model base year or a projection year) from user defined parameters
    def __init__(self, id, input_series, other=None):
        """define scenarios.
            For a base scenario, read and save building stock inputs
            For a projection scenario, some fields are projected from base scenario based on the projection scenario inputs

        #param str id: id of scenario, 'base' refers to modeled base year scenario
        #param Pandas Series input_series: contains information of building stock energy usages from scenario inputs
        #param class scenarios other: base scenario. When other == None, then define a base scenario from scenario inputs
        """
        self.id = id
        self.year = int(input_series["year"])
        self.hp_type_heat = input_series.pop("heat_hp_type")
        self.hp_type_dhw = input_series.pop("dhw_hp_type")
        self.cook_efficiency = input_series.pop("cook_eff")
        input_series = input_series.astype("float64")
        self.cool_energy_intensity = input_series["cool_energy_intensity(relative)"]
        self.stats = input_series
        if other == None:
            self.stats = (
                input_series.dropna()
            )  # dropna is supposed to drop rows that only used for defining future scenarios
            self._compute_base_scenario()

        else:
            self._compute_project_scenario(other)

    def _compute_base_scenario(self):
        self.pop = self.stats["pop"]
        self.floor_area_m2 = {
            "res": self.stats["res_area_m2"],
            "com": self.stats["com_area_m2"],
        }
        self.hp_heat_area_m2 = (
            self.stats["frac_hp_res_heat"] * self.floor_area_m2["res"]
            + self.stats["frac_hp_com_heat"] * self.floor_area_m2["com"]
        )
        self.hp_heat_frac = self.hp_heat_area_m2 / (
            self.floor_area_m2["res"] + self.floor_area_m2["com"]
        )
        self.resist_heat_area_m2 = (
            self.stats["frac_resist_res_heat"] * self.floor_area_m2["res"]
            + self.stats["frac_resist_com_heat"] * self.floor_area_m2["com"]
        )
        self.resist_heat_frac = self.resist_heat_area_m2 / (
            self.floor_area_m2["res"] + self.floor_area_m2["com"]
        )
        self.elec_cool_m2 = (
            self.stats["frac_elec_res_cool"] * self.floor_area_m2["res"]
            + self.stats["frac_elec_com_cool"] * self.floor_area_m2["com"]
        )
        return

    def _compute_project_scenario(self, other):
        self.pop = other.pop * (1 + self.stats["pop_ann_grow_rate"]) ** (
            self.year - other.year
        )
        self.floor_area_m2 = {}
        for clas in const.classes:
            if np.isnan(self.stats[f"{clas}_area_ann_grow_rate"]):
                # floor area growth rate is assumed to be proportional to population growth by default, unless user defines a growth rate for floor area
                self.floor_area_m2[clas] = other.floor_area_m2[clas] * (
                    1 + self.stats["pop_ann_grow_rate"]
                ) ** (self.year - other.year)
            else:
                self.floor_area_m2[clas] = other.floor_area_m2[clas] * (
                    1 + self.stats[f"{clas}_area_ann_grow_rate"]
                ) ** (self.year - other.year)

            if not np.isnan(self.stats[f"frac_hp_{clas}_heat"]):
                # calculate fraction of floor area using heat pump, resistance heat or fossil fuel furnace as major heating appliances
                # Users can either define the fractions directly or the shift in major heating fuels from the base scenario or assume BAU cases
                self.stats[f"ff2hp_{clas}"] = (
                    1
                    - self.stats[f"frac_ff_{clas}_heat"]
                    / other.stats[f"frac_ff_{clas}_heat"]
                )
            elif not np.isnan(self.stats[f"ff2hp_{clas}"]):
                self.stats[f"frac_hp_{clas}_heat"] = (
                    other.stats[f"frac_hp_{clas}_heat"]
                    + other.stats[f"frac_resist_{clas}_heat"]
                    * self.stats[f"resist2hp_{clas}"]
                    + other.stats[f"frac_ff_{clas}_heat"] * self.stats[f"ff2hp_{clas}"]
                )
                self.stats[f"frac_resist_{clas}_heat"] = other.stats[
                    f"frac_resist_{clas}_heat"
                ] * (1 - self.stats[f"resist2hp_{clas}"])
                self.stats[f"frac_ff_{clas}_heat"] = other.stats[
                    f"frac_ff_{clas}_heat"
                ] * (1 - self.stats[f"ff2hp_{clas}"])
            else:
                self.stats[f"frac_hp_{clas}_heat"] = other.stats[f"frac_hp_{clas}_heat"]
                self.stats[f"frac_resist_{clas}_heat"] = other.stats[
                    f"frac_resist_{clas}_heat"
                ]
                self.stats[f"frac_ff_{clas}_heat"] = other.stats[f"frac_ff_{clas}_heat"]
                self.stats[f"ff2hp_{clas}"] = 0
                self.stats[f"resist2hp_{clas}"] = 0

            if np.isnan(self.stats[f"frac_elec_{clas}_cool"]):
                self.stats[f"frac_elec_{clas}_cool"] = other.stats[
                    f"frac_elec_{clas}_cool"
                ]
            if np.isnan(self.stats[f"frac_ff_dhw_{clas}"]):
                self.stats[f"frac_ff_dhw_{clas}"] = other.stats[f"frac_ff_dhw_{clas}"]
            cook_other = "cook" if clas == "com" else "other"
            if np.isnan(self.stats[f"frac_ff_{cook_other}_{clas}"]):
                self.stats[f"frac_ff_{cook_other}_{clas}"] = other.stats[
                    f"frac_ff_{cook_other}_{clas}"
                ]

        self.hp_heat_area_m2 = (
            self.stats["frac_hp_res_heat"] * self.floor_area_m2["res"]
            + self.stats["frac_hp_com_heat"] * self.floor_area_m2["com"]
        )

        self.hp_heat_frac = self.hp_heat_area_m2 / (
            self.floor_area_m2["res"] + self.floor_area_m2["com"]
        )
        self.resist_heat_area_m2 = (
            self.stats["frac_resist_res_heat"] * self.floor_area_m2["res"]
            + self.stats["frac_resist_com_heat"] * self.floor_area_m2["com"]
        )
        self.resist_heat_frac = self.resist_heat_area_m2 / (
            self.floor_area_m2["res"] + self.floor_area_m2["com"]
        )
        self.elec_cool_m2 = (
            self.stats["frac_elec_res_cool"] * self.floor_area_m2["res"]
            + self.stats["frac_elec_com_cool"] * self.floor_area_m2["com"]
        )
        return

    def floor_area_growth(self, other):
        # compound floor area growth
        return (self.floor_area_m2["res"] + self.floor_area_m2["com"]) / (
            other.floor_area_m2["res"] + other.floor_area_m2["com"]
        )

    def floor_area_growth_type(self, other, clas):
        # compound floor area growth by building type
        return self.floor_area_m2[clas] / other.floor_area_m2[clas]

    def frac_hp_growth(self, other):
        # floor area growth ratio that use hp as main heating appliance
        return self.hp_heat_area_m2 / other.hp_heat_area_m2

    def frac_resist_growth(self, other):
        # floor area growth ratio that use resistance heat as main heating source
        return self.resist_heat_area_m2 / other.resist_heat_area_m2

    def frac_cool_growth(self, other):
        # floor area growth ratio that have electric air conditioning
        return self.elec_cool_m2 / other.elec_cool_m2

    def frac_htg_ff2hp(self, other, clas):
        # fraction of floor area electrified for heating
        return other.stats[f"frac_ff_{clas}_heat"] - self.stats[f"frac_ff_{clas}_heat"]

    def frac_dhw_ff2hp(self, other, clas):
        # fraction of floor area electrified for dhw
        return other.stats[f"frac_ff_dhw_{clas}"] - self.stats[f"frac_ff_dhw_{clas}"]

    def frac_cook_ff2hp(self, other, clas):
        # fraction of floor area electrified for cooking
        cook_other = "cook" if clas == "com" else "other"
        return (
            other.stats[f"frac_ff_{cook_other}_{clas}"]
            - self.stats[f"frac_ff_{cook_other}_{clas}"]
        )

    def frac_cooling_eff_change(self, other):
        # account for future cooling efficiency improve
        return self.cool_energy_intensity / other.cool_energy_intensity

    def compare_hp_heat_type(self, other):
        # return True if the heat pump type for projection scenario is the same as that of base scenario
        return self.hp_type_heat == other.hp_type_heat


def temp_to_energy(temp_series, hourly_fits_df, db_wb_fit, base_scen, hp_heat_cop):
    """Compute baseload, heating, and cooling electricity for a certain hour of year under model base year scenario

    :param pandas.Series temp_series: data for the given hour.
    :param pandas.DataFrame hourly_fits_df: hourly and week/weekend breakpoints and
        coefficients for electricity use equations.
    :param class scenarios base_scen: reference scenario class
    :return: (*list*) -- [baseload, heating, cooling]
    """
    temp = temp_series["temp_c"]
    temp_wb = temp_series["temp_c_wb"]
    dark_frac = temp_series["hourly_dark_frac"]
    zone_hour = temp_series["hour_local"]

    heat_eng = 0
    mid_cool_eng = 0
    cool_eng = 0
    hp_eng = 0
    resist_eng = 0

    wk_wknd = (
        "wk"
        if temp_series["weekday"] < 5 and ~temp_series["holiday"]  # boolean value
        else "wknd"
    )

    (
        t_bpc,
        t_bph,
        i_heat,
        s_heat,
        s_dark,
        i_cool,
        s_cool_db,
        s_cool_wb,
    ) = (
        hourly_fits_df.at[zone_hour, f"t.bpc.{wk_wknd}.c"],
        hourly_fits_df.at[zone_hour, f"t.bph.{wk_wknd}.c"],
        hourly_fits_df.at[zone_hour, f"i.heat.{wk_wknd}"],
        hourly_fits_df.at[zone_hour, f"s.heat.{wk_wknd}"],
        hourly_fits_df.at[zone_hour, f"s.dark.{wk_wknd}"],
        hourly_fits_df.at[zone_hour, f"i.cool.{wk_wknd}"],
        hourly_fits_df.at[zone_hour, f"s.cool.{wk_wknd}.db"],
        hourly_fits_df.at[zone_hour, f"s.cool.{wk_wknd}.wb"],
    )

    base_eng = s_heat * t_bph + s_dark * dark_frac + i_heat

    if temp <= t_bph:
        heat_eng = -s_heat * (t_bph - temp)
        # seperate resistance heat and heat pump energy by COP. HP assumed to be midium performance HP that's used in Mike's Joule paper
        # ! Since we now have our estimation of current heat pump adoption rate, this function can be merged to the same one in zone_profile_generator.py
        cop_hp = hp_heat_cop.loc[round(temp, 1), "cop"]
        hp_eng = (
            heat_eng
            * (base_scen.hp_heat_frac / cop_hp)
            / (base_scen.resist_heat_frac + base_scen.hp_heat_frac / cop_hp)
        )
        resist_eng = (
            heat_eng
            * (base_scen.resist_heat_frac)
            / (base_scen.resist_heat_frac + base_scen.hp_heat_frac / cop_hp)
        )

    if temp >= t_bph:
        cool_eng = (
            s_cool_db * temp
            + s_cool_wb
            * (
                temp_wb
                - (db_wb_fit[0] * temp**2 + db_wb_fit[1] * temp + db_wb_fit[2])
            )
            + i_cool
        )

    if temp > t_bpc and temp < t_bph:
        mid_cool_eng = ((temp - t_bpc) / (t_bph - t_bpc)) ** 2 * (
            s_cool_db * t_bph
            + s_cool_wb
            * (
                temp_wb
                - (db_wb_fit[0] * temp**2 + db_wb_fit[1] * temp + db_wb_fit[2])
            )
            + i_cool
        )

    return [base_eng, hp_eng, resist_eng, max(cool_eng, 0) + max(mid_cool_eng, 0)]


def scale_energy(
    base_energy, temp_df, base_scen, new_scen, midperfhp_cop, advperfhp_cop
):
    """project energy consumption for each projection scenarios from the base scenario

    :param pandas.DataFrame base_energy: dataframe of disaggregated electricity consumptions for all weather years
    :param pandas.DataFrame temp_df: weather records the given hours

    :return (*pandas.DataFrame*) scen_load_MWh -- hourly electricity consumption induced by heat pump heating, resistance heating, cooling, and 'base' loads for a projection scenario
    """
    base_load_scaler = new_scen.floor_area_growth(base_scen)
    cool_load_scaler = new_scen.frac_cooling_eff_change(
        base_scen
    ) * new_scen.frac_cool_growth(base_scen)
    resist_load_scaler = new_scen.frac_resist_growth(base_scen)
    if new_scen.compare_hp_heat_type(base_scen):
        hp_load_scaler = new_scen.frac_hp_growth(base_scen)
    else:
        if base_scen.hp_type_heat == "midperfhp":
            base_hp_heat_cop = midperfhp_cop
        elif base_scen.hp_type_heat == "advperfhp":
            base_hp_heat_cop = advperfhp_cop
        if new_scen.hp_type_heat == "midperfhp":
            new_hp_heat_cop = midperfhp_cop
        elif new_scen.hp_type_heat == "advperfhp":
            new_hp_heat_cop = advperfhp_cop

        hp_cop_adv = new_hp_heat_cop / base_hp_heat_cop  # <1

        hp_cop_scaler = pd.Series(dtype="float64")
        for i in temp_df.index:
            hp_cop_scaler.loc[i] = hp_cop_adv.loc[
                round(temp_df.loc[i, "temp_c"], 1), "cop"
            ]
        hp_load_scaler = hp_cop_scaler * new_scen.frac_hp_growth(base_scen)

    scen_load_MWh = pd.DataFrame(index=base_energy.index)
    scen_load_MWh["base_load_mw"] = base_energy["base_load_mw"] * base_load_scaler
    if (
        new_hp_profile == "elec"
    ):  # if user select to use the current electricity heat pump consumption profiles for heating load projection
        scen_load_MWh["heat_hp_load_mw"] = (
            base_energy["heat_hp_load_mw"] * hp_load_scaler
        )
    else:
        scen_load_MWh["heat_existing_hp_load_mw"] = base_energy["heat_hp_load_mw"]
    scen_load_MWh["heat_resist_load_mw"] = (
        base_energy["heat_resist_load_mw"] * resist_load_scaler
    )
    scen_load_MWh["cool_load_mw"] = base_energy["cool_load_mw"] * cool_load_scaler

    return scen_load_MWh


def ff_electrify_profiles(weather_years, puma_data, new_scen, base_scen):
    """return hourly electricity loads for a projection scenario from converting fossil fuel heating, dhw and cooking to electric ones
    :param list weather_years: user defined year(s) of weather profile for load projection
    :param pandas.DataFrame puma_data: puma data within zone, output of zone_shp_overlay()
    :param class scenarios new_scen: projection scenario class
    :param class scenarios base_scen: reference scenario class

    :return (*pandas.DataFrame*) ff2hp_load_MWh -- hourly projection load from converting fossil fuel consumption to electricity for projection scenarios given weather conditions from selected weather years
    """

    def ff2hp_dhw_profiles(clas):
        ff2hp_dhw_pumas_yrs = pd.DataFrame()
        for weather_year, state in product(weather_years, zone_states):
            if not os.path.isfile(
                os.path.join(
                    os.path.dirname(__file__),
                    "Profiles",
                    f"elec_dhw_ff2hp_{clas}_{state}_{weather_year}_{hp_type_dhw}_mw.csv",
                )
            ):
                print(f"generating hot water profiles for {state}...")
                generate_dhw_profiles(weather_year, [state], clas, hp_type_dhw)

        for weather_year in weather_years:
            ff2hp_dhw_pumas = pd.concat(
                list(
                    pd.Series(data=zone_states).apply(
                        lambda x: pd.read_csv(
                            os.path.join(
                                os.path.dirname(__file__),
                                "Profiles",
                                f"elec_dhw_ff2hp_{clas}_{x}_{weather_year}_{hp_type_dhw}_mw.csv",
                            )
                        )
                    )
                ),
                axis=1,
            )
            ff2hp_dhw_pumas_yrs = pd.concat(
                [ff2hp_dhw_pumas_yrs, ff2hp_dhw_pumas], ignore_index=True
            )
        ff2hp_dhw_yrs = (
            ff2hp_dhw_pumas_yrs[puma_data.index]
            .mul(puma_data["frac_in_zone"])
            .sum(axis=1)
        )
        return ff2hp_dhw_yrs

    def ff2hp_htg_profiles(clas):
        ff2hp_htg_pumas_yrs = pd.DataFrame()
        for weather_year, state in product(weather_years, zone_states):
            if not os.path.isfile(
                os.path.join(
                    os.path.dirname(__file__),
                    "Profiles",
                    f"elec_htg_ff2hp_{clas}_{state}_{weather_year}_{hp_type_heat}_mw.csv",
                )
            ):
                print(f"generating ff heating profiles for {state}...")
                generate_htg_profiles(weather_year, [state], clas, hp_type_heat)

        for weather_year in weather_years:
            ff2hp_htg_pumas = pd.concat(
                list(
                    pd.Series(data=zone_states).apply(
                        lambda x: pd.read_csv(
                            os.path.join(
                                os.path.dirname(__file__),
                                "Profiles",
                                f"elec_htg_ff2hp_{clas}_{x}_{weather_year}_{hp_type_heat}_mw.csv",
                            )
                        )
                    )
                ),
                axis=1,
            )
            ff2hp_htg_pumas_yrs = pd.concat(
                [ff2hp_htg_pumas_yrs, ff2hp_htg_pumas], ignore_index=True
            )
        ff2hp_htg_yrs = (
            ff2hp_htg_pumas_yrs[puma_data.index]
            .mul(puma_data["frac_in_zone"])
            .sum(axis=1)
        )
        return ff2hp_htg_yrs

    def ff2hp_cook_profiles(clas):
        for state in zone_states:
            if not os.path.isfile(
                os.path.join(
                    os.path.dirname(__file__),
                    "Profiles",
                    f"elec_cook_ff2hp_{clas}_{state}_{base_year}_{cook_eff}_mw.csv",
                )
            ):
                print(f"generating ff cooking profiles for {state}...")
                generate_cook_profiles(base_year, [state], clas, cook_eff)

        ff2hp_cook_pumas = pd.concat(
            list(
                pd.Series(data=zone_states).apply(
                    lambda x: pd.read_csv(
                        os.path.join(
                            os.path.dirname(__file__),
                            "Profiles",
                            f"elec_cook_ff2hp_{clas}_{x}_{base_year}_{cook_eff}_mw.csv",
                        ),
                        index_col=0,
                    )
                )
            )
        ).iloc[:, 0]
        ff2hp_cook = (
            ff2hp_cook_pumas.loc[puma_data.index].mul(puma_data["frac_in_zone"]).sum()
        )
        return ff2hp_cook

    zone_states = list(set(puma_data["state"]))
    hours_utc_weather_years = pd.date_range(
        start=f"{weather_years[0]}-01-01",
        end=f"{weather_years[-1]+1}-01-01",
        freq="H",
        tz="UTC",
    )[:-1]
    hp_type_dhw = new_scen.hp_type_dhw
    hp_type_heat = new_scen.hp_type_heat
    cook_eff = new_scen.cook_efficiency
    ff2hp_load_MWh = pd.DataFrame(index=hours_utc_weather_years)
    for clas in const.classes:
        frac_dhw_ff2hp = new_scen.frac_dhw_ff2hp(base_scen, clas)
        if frac_dhw_ff2hp != 0:
            ff2hp_load_MWh[f"dhw_{clas}"] = ff2hp_dhw_profiles(clas).to_list()
            ff2hp_load_MWh[f"dhw_{clas}"] = (
                ff2hp_load_MWh[f"dhw_{clas}"]
                * new_scen.floor_area_growth_type(base_scen, clas)
                * frac_dhw_ff2hp
            )  # scale energy consumption by floor area information

        frac_htg_ff2hp = new_scen.frac_htg_ff2hp(base_scen, clas)
        if new_hp_profile == "ff" and frac_htg_ff2hp != 0:
            ff2hp_load_MWh[f"htg_{clas}"] = ff2hp_htg_profiles(clas).to_list()
            ff2hp_load_MWh[f"htg_{clas}"] = (
                ff2hp_load_MWh[f"htg_{clas}"]
                * new_scen.floor_area_growth_type(base_scen, clas)
                * frac_htg_ff2hp
            )

        frac_cook_ff2hp = new_scen.frac_cook_ff2hp(base_scen, clas)
        if frac_cook_ff2hp != 0:
            ff2hp_load_MWh[f"cook_{clas}"] = ff2hp_cook_profiles(clas)
            ff2hp_load_MWh[f"cook_{clas}"] = (
                ff2hp_load_MWh[f"cook_{clas}"]
                * new_scen.floor_area_growth_type(base_scen, clas)
                * new_scen.frac_cook_ff2hp(base_scen, clas)
            )
    return ff2hp_load_MWh


def predict_scenario(zone_name, zone_name_shp, base_scen, new_scens, weather_years):
    """load projection for one zone for all selected weather years.

    :param str zone_name: name of load zone used to save profile.
    :param str zone_name_shp: name of load zone within shapefile.
    :param class scenarios new_scen: projection scenario class
    :param class scenarios base_scen: reference scenario class
    :param list weather_years: user defined year(s) of weather profile for load projection
    """
    hours_utc_weather_years = pd.date_range(
        start=f"{weather_years[0]}-01-01",
        end=f"{weather_years[-1]+1}-01-01",
        freq="H",
        tz="UTC",
    )[:-1]
    # prepare weather dataframe
    puma_data_zone = zone_shp_overlay(zone_name_shp, zone_shp, pumas_shp)
    temp_df = pd.DataFrame()
    for year in weather_years:
        hours_utc_year = pd.date_range(
            start=f"{year}-01-01", end=f"{year+1}-01-01", freq="H", tz="UTC"
        )[:-1]
        temp_df_yr, stats = zonal_data(puma_data_zone, hours_utc_year, year)
        temp_df = pd.concat([temp_df, temp_df_yr], ignore_index=True)
    temp_df.index = hours_utc_weather_years

    # compute least-square estimator for relation between WBT and DBT
    t_bpc_start = 10
    db_wb_regr_df = temp_df[
        (temp_df.index.year == base_year) & (temp_df["temp_c"] >= t_bpc_start)
    ]
    db_wb_fit = np.polyfit(db_wb_regr_df["temp_c"], db_wb_regr_df["temp_c_wb"], 2)

    hourly_fits_df = pd.read_csv(
        f"https://besciences.blob.core.windows.net/datasets/bldg_el/dayhour_fits/{zone_name}_dayhour_fits_{base_year}.csv",
        index_col=0,
    )

    midperfhp_cop = pd.read_csv(f"./data/cop_temp_htg_midperfhp.csv")
    advperfhp_cop = pd.read_csv(f"./data/cop_temp_htg_advperfhp.csv")
    midperfhp_cop.index = midperfhp_cop["temp"]
    advperfhp_cop.index = advperfhp_cop["temp"]

    if base_scen.hp_type_heat == "midperfhp":
        base_hp_heat_cop = midperfhp_cop
    elif base_scen.hp_type_heat == "advperfhp":
        base_hp_heat_cop = advperfhp_cop
    else:
        raise KeyError("hp type not defined. Choose from 'midperfhp' and 'advperfhp'")

    # run reference scenario
    zone_profile_refload_MWh = pd.DataFrame(  # noqa: N806
        {"hour_utc": hours_utc_weather_years}
    )

    energy_list = zone_profile_refload_MWh.hour_utc.apply(
        lambda x: temp_to_energy(
            temp_df.loc[x], hourly_fits_df, db_wb_fit, base_scen, base_hp_heat_cop
        )
    )

    (
        zone_profile_refload_MWh["base_load_mw"],
        zone_profile_refload_MWh["heat_hp_load_mw"],
        zone_profile_refload_MWh["heat_resist_load_mw"],
        zone_profile_refload_MWh["cool_load_mw"],
    ) = (
        energy_list.apply(lambda x: x[0]),
        energy_list.apply(lambda x: x[1]),
        energy_list.apply(lambda x: x[2]),
        energy_list.apply(lambda x: x[3]),
    )
    zone_profile_refload_MWh.set_index("hour_utc", inplace=True)

    # scale energy to each projection scenarios
    elec_profile_load_MWh = {}
    zone_profile_load_MWh = {}
    zone_profile_load_MWh["base"] = zone_profile_refload_MWh
    ff2hp_profile_load_MWh = {}
    for id, scenario in new_scens.items():
        elec_profile_load_MWh[id] = scale_energy(
            zone_profile_refload_MWh,
            temp_df,
            base_scen,
            scenario,
            midperfhp_cop,
            advperfhp_cop,
        )
        ff2hp_profile_load_MWh[id] = ff_electrify_profiles(
            weather_years, puma_data_zone, scenario, base_scen
        )
        zone_profile_load_MWh[id] = pd.concat(
            [elec_profile_load_MWh[id], ff2hp_profile_load_MWh[id]], axis=1
        )

    return zone_profile_load_MWh


if __name__ == "__main__":
    # Use base_year for model fitting
    base_year = const.base_year

    # Weather year to produce load profiles. If multiple years, then time series result will show for more than one year
    weather_years = [2018, 2019]

    # new heat pump load profile assumption. User can select whether to use electric profile or fossil fuel profile to estimate
    # electrified fossil fuel consumption for heating
    new_hp_profile = "elec"  # "elec" or "ff"

    # Reading Balancing Authority and Pumas shapefiles for overlaying
    zone_shp = read_shapefile(
        "https://besciences.blob.core.windows.net/shapefiles/USA/balancing-authorities/ba_area/ba_area.zip"
    )
    pumas_shp = read_shapefile(
        "https://besciences.blob.core.windows.net/shapefiles/USA/pumas-overlay/pumas_overlay.zip"
    )

    zone_names = [
        "NYIS-ZONA",
    ]

    zone_name_shps = [
        "NYISO-A",
    ]

    for i in range(len(zone_names)):
        zone_name, zone_name_shp = zone_names[i], zone_name_shps[i]

        scen_data = pd.read_csv(
            os.path.join(
                os.path.dirname(__file__),
                "projection",
                "scenario_inputs",
                f"{zone_name}_stats.csv",
            ),
            index_col=0,
        )

        base_scenarios = scenarios("base", scen_data.pop("yr2019"))
        print(f"base scenario: year {base_year}, weather year: {weather_years}")

        proj_scenarios = {}
        for name, values in scen_data.iteritems():
            proj_scenarios[name] = scenarios(name, values, base_scenarios)
            print(f"projection scenario {name}, year {proj_scenarios[name].year}")

        os.makedirs(
            os.path.join(os.path.dirname(__file__), "Profiles"),
            exist_ok=True,
        )
        zone_profile_load_MWh = predict_scenario(
            zone_name, zone_name_shp, base_scenarios, proj_scenarios, weather_years
        )

        os.makedirs(
            os.path.join(os.path.dirname(__file__), "projection", "results"),
            exist_ok=True,
        )

        for name, values in zone_profile_load_MWh.items():
            zone_profile_load_MWh[name].to_csv(
                os.path.join(
                    os.path.dirname(__file__),
                    "projection",
                    "results",
                    f"{zone_name}_{name}_mwh.csv",
                )
            )
