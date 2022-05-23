import os
import pickle as pkl
import urllib.request

import numpy as np
import pandas as pd


def get_census_data(download_path):
    """download county population data from USA Census website

    param string download_path: folder to store the downloaded file
    :return: (*None*)
    """

    os.makedirs(download_path, exist_ok=True)
    census_file_link = "https://www2.census.gov/programs-surveys/popest/datasets/2010-2020/counties/totals/co-est2020-alldata.csv"
    urllib.request.urlretrieve(census_file_link, download_path+'/county_population.csv')


def get_crosswalk_data(download_path):
    """download FIPS-ZIP crosswalk data from USPS and convert to csv

    param string download_path: folder to store the downloaded file
    :return: (*None*)
    """
    
    # download
    os.makedirs(download_path, exist_ok=True)
    usps_file_link="https://www.huduser.gov/portal/datasets/usps/COUNTY_ZIP_122021.xlsx"
    urllib.request.urlretrieve(usps_file_link, download_path+'/county_to_zip.xlsx')

    # convert to csv
    usps_df = pd.read_excel(download_path+"/county_to_zip.xlsx")
    usps_df = usps_df.sort_values(by=["county", "zip"])
    usps_df.to_csv(download_path+"/county_to_zip.csv", index=False)
    os.remove(download_path+"/county_to_zip.xlsx")


def get_LSE_region_data(download_path):
    """download LSE service region data

    param string download_path: folder to store the downloaded file
    :return: (*None*)
    """
    
    os.makedirs(download_path, exist_ok=True)
    iou_file_link = "https://data.openei.org/files/4042/iou_zipcodes_2019.csv"
    niou_file_link = "https://data.openei.org/files/4042/non_iou_zipcodes_2019.csv"

    urllib.request.urlretrieve(iou_file_link, download_path+"/iou_zipcodes_2019.csv")
    urllib.request.urlretrieve(niou_file_link, download_path+"/non_iou_zipcodes_2019.csv")
    

def get_county_fips_data(download_path):
    """download county FIPS data

    param string download_path: folder to store the downloaded file
    :return: (*None*)
    """
    
    os.makedirs(download_path, exist_ok=True)
    county_fips_link = "https://github.com/kjhealy/fips-codes/raw/master/county_fips_master.csv"
    urllib.request.urlretrieve(county_fips_link, download_path+"/county_fips_master.csv")


def eiaid_to_zip(raw_data_path, cache_path):
    """Find the service region (list of ZIP codes) for every LSE identified by their EIA ID
        Create a dictionary with EIA ID as keys for list of zip codes in the cache folder
        
    param string raw_data_path: folder that contains downloaded raw data
    param string cache_path: folder to store processed cache files        
    :return: (*None*)
    """

    assert os.path.isfile(raw_data_path + "/iou_zipcodes_2019.csv"), "Input file iou_zipcodes_2019.csv does not exist."
    assert os.path.isfile(raw_data_path + "/non_iou_zipcodes_2019.csv"), "Input file non_iou_zipcodes_2019.csv does not exist."

    iou_df = pd.read_csv(raw_data_path + "/iou_zipcodes_2019.csv")
    niou_df = pd.read_csv(raw_data_path + "/non_iou_zipcodes_2019.csv")

    all_ids = list(set(iou_df['eiaid'].to_list() + niou_df['eiaid'].to_list()))

    id2zip = {}

    for i in all_ids:
        iouzips = iou_df.loc[iou_df['eiaid']==i]['zip'].to_list()
        niouzips = niou_df.loc[niou_df['eiaid']==i]['zip'].to_list()
        allzips = sorted(list(set(iouzips + niouzips)))
        id2zip[i] = allzips

    with open(cache_path + "/eiaid2zip.pkl", 'wb') as fh:
        pkl.dump(id2zip, fh)


def eiaid_to_fips(raw_data_path, cache_path):
    """Find the service region (list of FIPS codes) for every LSE identified by their EIA ID
        Create a dictionary with EIA ID as keys for list of FIPS codes in the cache folder
        
    param string raw_data_path: folder that contains downloaded raw data
    param string cache_path: folder to store processed cache files        
    :return: (*None*)
    """

    assert os.path.isfile(cache_path + "/eiaid2zip.pkl"), "Cached file eiaid2zip.pkl does not exist."
    assert os.path.isfile(cache_path + "/zip2fips.pkl"), "Cached file zip2fips.pkl does not exist."
    assert os.path.isfile(cache_path + "/fips_population.pkl"), "Cached file fips_population.pkl does not exist."
    
    # load cached data
    with open(cache_path+"/eiaid2zip.pkl", 'rb') as fh:
        eiaid2zip = pkl.load(fh)

    with open(cache_path+"/zip2fips.pkl", 'rb') as fh:
        zip2fips = pkl.load(fh)

    with open(cache_path+"/fips_population.pkl", 'rb') as fh:
        fips_pop_df = pkl.load(fh)

    eia_keys = list(eiaid2zip.keys())
    eia_num = len(eia_keys)

    eiaid2fips = {}
    for i in range(eia_num):
        key = eia_keys[i]
        zips = eiaid2zip[key]
        all_fips = []
        all_pops = []

        # aggregate all zip codes 
        for j in zips:
            try:
                fipss, weights = zip2fips[j]
            except:
                continue
            fips_num = len(fipss)
        
            for k in range(fips_num):
                # total population in this fips
                total_pops = fips_pop_df.loc[fips_pop_df["fips"]==fipss[k],"population"].to_list()[0]
                if fipss[k] in all_fips:
                    all_pops[all_fips.index(fipss[k])] += weights[k]
                else:
                    all_fips.append(fipss[k])
                    all_pops.append(weights[k]*total_pops)
                    
        eiaid2fips.update({key: [all_fips, all_pops]})

    with open(cache_path+"/eiaid2fips.pkl", 'wb') as fh:
        pkl.dump(eiaid2fips, fh)
    


def fips_zip_conversion(raw_data_path, cache_path):
    """Create a two-way mapping for all ZIP and FIPS in the crosswalk data
        save to dictionary files for future use

    param string raw_data_path: folder that contains downloaded raw data
    param string cache_path: folder to store processed cache files
    :return: (*None*)
    """
    assert os.path.isfile(raw_data_path + "/county_to_zip.csv"), "Input file county_to_zip.csv does not exist."

    df_raw = pd.read_csv(raw_data_path+"/county_to_zip.csv")
    
    all_fips = df_raw['county'].astype('int32')
    all_zip = df_raw['zip'].astype('int32')
    all_weights = df_raw['tot_ratio']

    # create zip -> counties mapping
    zip2fips = {}
    
    for i in pd.unique(all_zip):
        idx = all_zip.index[all_zip == i]
        
        cty = all_fips[idx].tolist()
        wgt = all_weights[idx].tolist()
        
        zip2fips.update({i: (cty, wgt)})

    with open(cache_path+"/zip2fips.pkl", 'wb') as fh:
        pkl.dump(zip2fips, fh)

    # create county -> zips mapping
    fips2zip = {}
    
    for i in pd.unique(all_fips):
        idx = all_fips.index[all_fips == i]
        
        zips = all_zip[idx].tolist()
        wgt = all_weights[idx].tolist()
        
        fips2zip.update({i: (zips, wgt)})


    with open(cache_path+"/fips2zip.pkl", 'wb') as fh:
        pkl.dump(fips2zip, fh)


def get_fips_population(raw_data_path, cache_path):
    """Match county population and county FIPS data to produce concise FIPS population
        save to a dictonary in cache folder with key being 5-digit FIPS codes.

    param string raw_data_path: folder that contains downloaded raw data
    param string cache_path: folder to store processed cache files
    :return: (*None*)
    """

    assert os.path.isfile(raw_data_path + "/county_population.csv"), "Input file county_population.csv does not exist."
    assert os.path.isfile(raw_data_path + "/county_fips_master.csv"), "Input file county_fips_master.csv does not exist."

    cty_pop_df = pd.read_csv(raw_data_path+"/county_population.csv", encoding='cp1252')
    cty_name_df = pd.read_csv(raw_data_path+"/county_fips_master.csv", encoding='cp1252')

    pops = np.zeros(len(cty_name_df.index))
    for cty in cty_name_df.index:
        name = cty_name_df['county_name'][cty]
        state = cty_name_df['state_name'][cty]

        tmp = cty_pop_df.loc[(cty_pop_df["CTYNAME"]==name) & (cty_pop_df["STNAME"]==state), "POPESTIMATE2019"]
        # if no data
        if len(tmp) == 0:
            pops[cty] = 0
        elif len(tmp) == 1:
            pops[cty] = tmp
        else:
            pops[cty] = tmp.to_list()[0]
            
    cty_name_df["population"] = pops

    new_df = cty_name_df.loc[:, ["fips", "county_name", "population"]]

    with open(cache_path+"/fips_population.pkl", 'wb') as fh:
        pkl.dump(new_df, fh)



def get_zip_population(raw_data_path, cache_path):
    """Compute population of each ZIP code using percentage share and FIPS population
        save to a dictonary in cache folder with key being zip codes.
        
    param string raw_data_path: folder that contains downloaded raw data
    param string cache_path: folder to store processed cache files
    :return: (*None*)
    """

    assert os.path.isfile(cache_path + "/fips_population.pkl"), "Cached fips_population.pkl does not exist."
    assert os.path.isfile(cache_path + "/zip2fips.pkl"), "Cached file zip2fips.pkl does not exist."
    
    with open(cache_path+"/fips_population.pkl", 'rb') as fh:
        fips_pop_df = pkl.load(fh)
    
    with open(cache_path+"/zip2fips.pkl", 'rb') as fh:
        zip2fips = pkl.load(fh)

    zips = sorted(list(zip2fips.keys()))

    zip_num = len(zips)
    zip_pops = {}

    for i in range(zip_num):
        z = zips[i]

        fipss = zip2fips[z][0]
        pcts = zip2fips[z][1]
        zip_pops[z] = 0
        for f in range(len(fipss)):
            tmp = fips_pop_df.loc[fips_pop_df['fips']==fipss[f], 'population'].values
            if len(tmp) > 0:
                zip_pops[z] += tmp[0] * pcts[f]


    with open(cache_path+"/zip_population.pkl", 'wb') as fh:
        pkl.dump(zip_pops, fh)

