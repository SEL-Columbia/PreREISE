import pickle as pkl
import pandas as pd
import numpy as np
import scipy.io as spio
import geopy
import requests
import time
import os
import urllib.request
from geopy.extra.rate_limiter import RateLimiter

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


def get_bus_pos(case_path):
    """Read a .mat case and extract the lat/lon coordinate of all buses 

    param string case_path: path to a .mat case file representing a power network
    :return: (*numpy.ndarray*) -- a list of (busID, lat, lon) stored as a numpy array
    """

    case = spio.loadmat(case_path)
    mpc = case['mpc']

    bus = mpc['bus'][0][0]
    bus_num = bus.shape[0]
    bus2sub = mpc['bus2sub'][0][0]
    sub = mpc['sub'][0][0]
    subid = mpc['subid'][0][0]

    bus_pos = np.zeros((bus_num, 3))

    # [id, lat, long]
    subid_list = subid[:, 0].tolist()
    for i in range(bus_num):
            
        bus_sub_idx = subid_list.index(bus2sub[i, 0])
        bus_pos[i, 0] = bus[i, 0]
        bus_pos[i, 1] = sub[bus_sub_idx, 2]
        bus_pos[i, 2] = sub[bus_sub_idx, 3]

    return bus_pos



def get_bus_fips(case_path, start_idx = 0):
    """Try to get FIPS of each bus in a case mat using FCC AREA API
        Can take hours to run, save to cache file for future use

    param string case_path: path to a .mat case file representing a power network
    param int start_idx: pointer to the index of a bus to start query from
    :return: (*None*)
    """
    file_dir = os.path.abspath(os.path.dirname(__file__))

    bus_pos = get_bus_pos(case_path)
    bus_fips_dict = {'busid': bus_pos[:, 0],
                    'latitude': bus_pos[:, 1],
                    'longitude': bus_pos[:, 2],
                    'fips': [0]*bus_num}
        
        
    url = "https://geo.fcc.gov/api/census/area"

    for i in range(start_idx, bus_num):
        if i % 1000 == 0:
            with open(os.path.join(file_dir, "../cache/bus_fips_usa.pkl"), 'wb') as fh:
                pkl.dump(bus_fips_dict, fh)
                
        pos = bus_pos[i, :]
        params = {"latitude" : pos[1],
              "longitude" : pos[2],
              "format": "json"
        }
        r = requests.get(url, params=params)

        assert r.status_code == 200
        
        res = r.json()
        if res['County']['FIPS'] is not None:
            bus_fips_dict['fips'][i] = int(res['County']['FIPS'])
        else:
            bus_fips_dict['fips'][i] = -1
            
        time.sleep(0.02)

    with open(os.path.join(file_dir, "../cache/bus_fips_usa.pkl"), 'wb') as fh:
        pkl.dump(bus_fips_dict, fh)
        

def cleanup_zip(zipdict):
    """Try to cleanup a zip dictionary obtained using online query by converting 
        to 5-digit integers. Several possible mis-format are considered
        
    param dicitonary zipdict: a dictionary containing raw zip-code of buses
    :return: (*dictionary*) -- a dictionary containing 5-digit zip codes
    """
    for i in range(len(zipdict['zip'])):
        try:
            zipdict['zip'][i] = int(zipdict['zip'][i])
        except Exception:
            if ':' in zipdict['zip'][i]:
                zipdict['zip'][i] = int(zipdict['zip'][i].split(':')[0])

            elif '-' in zipdict['zip'][i]:
                zipdict['zip'][i] = int(zipdict['zip'][i].split('-')[0])

            elif '‑' in zipdict['zip'][i]:
                zipdict['zip'][i] = int(zipdict['zip'][i].split('‑')[0])
                
            elif ' ' in zipdict['zip'][i]:
                success = 0
                for j in zipdict['zip'][i].split(' '):
                    if len(j) == 5:
                        zipdict['zip'][i] = int(j)
                        success = 1
                if success == 0:
                    zipdict['zip'][i] = -1
                    
            else:    
                try:
                    zipdict['zip'][i] = int(zipdict['zip'][i][0:5])
                except Exception:
                    zipdict['zip'][i] = -1

    return zipdict

def get_bus_zip(case_path, start_idx = 0):
    """Try to get ZIP of each bus in a case mat using geopy
        Can take hours to run, save to cache file for future use

    param string case_path: path to a .mat case file representing a power network
    param int start_idx: pointer to the index of a bus to start query from
    :return: (*None*)
    """
    bus_pos = get_bus_pos(case_path)
    bus_zip_dict = {'busid': bus_pos[:, 0],
                'latitude': bus_pos[:, 1],
                'longitude': bus_pos[:, 2],
                'zip': [0]*bus_num}

    geocoder = geopy.Nominatim(user_agent="BES")
    geocode = RateLimiter(geocoder.geocode, min_delay_seconds = 0.05, return_value_on_exception = None) 

    def get_zip_code(lat, lon):
        location = geocoder.reverse("{}, {}".format(lat, lon))
        if location is not None:
            address = location.raw['address']
        else:
            return -1

        if 'postcode' in address.keys():
            return address['postcode']
        else:
            return -1

    for i in range(start_idx, bus_num):
        bus_zip_dict['zip'][i] = int(get_zip_code(bus_zip_dict['latitude'][i], bus_zip_dict['longitude'][i]))

    bus_zip_dict = cleanup_zip(bus_zip_dict)
    
    with open("../cache/bus_zip_usa.pkl", 'wb') as fh:
        pkl.dump(bus_zip_dict, fh)


def eiaid_to_zip(raw_data_path, cache_path):
    """Find the service region (list of ZIP codes) for every LSE identified by their EIA ID
        Create a dictionary with EIA ID as keys for list of zip codes in the cache folder
        
    param string raw_data_path: folder that contains downloaded raw data
    param string cache_path: folder to store processed cache files        
    :return: (*None*)
    """

    iou_df = pd.read_csv(raw_data_path + "/iou_zipcodes_2019.csv")
    niou_df = pd.read_csv(raw_data_path + "/non_iou_zipcodes_2019.csv")

    all_ids = list(set(iou_df['eiaid'].to_list() + niou_df['eiaid'].to_list()))

    id2zip = {}

    for i in all_ids:
        iouzips = iou_df.loc[iou_df['eiaid']==i]['zip'].to_list()
        niouzips = niou_df.loc[niou_df['eiaid']==i]['zip'].to_list()
        allzips = list(set(iouzips + niouzips))
        allzips.sort()
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
    
    # load cached data
    with open(cache_path+"/eiaid2zip.pkl", 'rb') as fh:
        eiaid2zip = pkl.load(fh)

    with open(cache_path+"/zip2fips.pkl", 'rb') as fh:
        zip2fips = pkl.load(fh)

    fips_pop_df = pd.read_csv(cache_path+"/fips_population.csv")

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
                total_pops = fips_pop_df.loc[fips_pop_df["fips"]==fipss[k],"population"]
                total_pops = total_pops.to_list()[0]
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
    fip2zips = {}
    
    for i in pd.unique(all_fips):
        idx = all_fips.index[all_fips == i]
        
        zips = all_zip[idx].tolist()
        wgt = all_weights[idx].tolist()
        
        fip2zips.update({i: (zips, wgt)})


    with open(cache_path+"/fip2zips.pkl", 'wb') as fh:
        pkl.dump(fip2zips, fh)


def get_fips_population(raw_data_path, cache_path):
    """Match county population and county FIPS data to produce concise FIPS population
        save to a dictonary in cache folder with key being 5-digit FIPS codes.

    param string raw_data_path: folder that contains downloaded raw data
    param string cache_path: folder to store processed cache files
    :return: (*None*)
    """

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
    new_df.to_csv(cache_path+"/fips_population.csv", index=False)



def get_zip_population(raw_data_path, cache_path):
    """Compute population of each ZIP code using percentage share and FIPS population
        save to a dictonary in cache folder with key being zip codes.
        
    param string raw_data_path: folder that contains downloaded raw data
    param string cache_path: folder to store processed cache files
    :return: (*None*)
    """
    
    fips_pop_df = pd.read_csv(cache_path+"/fips_population.csv")
    
    with open(cache_path+"/zip2fips.pkl", 'rb') as fh:
        zip2fips = pkl.load(fh)

    zips = list(zip2fips.keys())
    zips.sort()

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
