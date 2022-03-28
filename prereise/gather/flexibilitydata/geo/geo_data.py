import pickle as pkl
import pandas as pd
import numpy as np
import scips as sp
import geopy
import requets
import time
from geopy.extra.rate_limiter import RateLimiter

"""
    eiaid_to_zip()

Find the service region (list of ZIP codes) for every LSE identified by their EIA ID
Create and save a dictionary with EIA ID as keys for list of zip codes 
"""
def eiaid_to_zip():
    iou_path = "../raw/iou_zipcodes_2019.csv"
    niou_path = "../raw/non_iou_zipcodes_2019.csv"
    dictpath = "../cache/eiaid2zip.pkl"

    iou_df = pd.read_csv(iou_path)
    niou_df = pd.read_csv(niou_path)

    all_ids = list(set(iou_df['eiaid'].to_list() + niou_df['eiaid'].to_list()))

    id2zip = {}

    for i in all_ids:
        iouzips = iou_df.loc[iou_df['eiaid']==i]['zip'].to_list()
        niouzips = niou_df.loc[niou_df['eiaid']==i]['zip'].to_list()
        allzips = list(set(iouzips + niouzips))
        allzips.sort()
        id2zip[i] = allzips

    with open(dictpath, 'wb') as fh:
        pickle.dump(id2zip, fh)

"""
    get_bus_pos(case_path)

Read a .mat case and extract the lat/lon coordinate of all buses
return a list of (busID, lat, lon)
"""
def get_bus_pos(case_path):
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


"""
    get_bus_fips(case_path)

Try to get FIPS of each bus in a case mat using FCC AREA API
Can take hours to run, save to cache file for future use
"""
def get_bus_fips(case_path):

    bus_pos = get_bus_pos(case_path)
    bus_fips_dict = {'busid': bus_pos[:, 0],
                    'latitude': bus_pos[:, 1],
                    'longitude': bus_pos[:, 2],
                    'fips': [0]*bus_num}
        
        
    url = "https://geo.fcc.gov/api/census/area"
    start_idx = 0

    for i in range(start_idx, bus_num):
        if i % 1000 == 0:
            print(i)
            with open("../cache/bus_fips_usa.pkl", 'wb') as fh:
                pkl.dump(bus_fips_dict, fh)
                
        pos = bus_pos[i, :]
        params = {'latitude' : pos[1],
              'longitude' : pos[2],
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

    with open("../cache/bus_fips_usa.pkl", 'wb') as fh:
        pkl.dump(bus_fips_dict, fh)

"""
    get_bus_zip(case_path)

Try to get ZIP of each bus in a case mat using geopy
Can take hours to run, save to cache file for future use
"""
def get_bus_zip(case_path):
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
        bus_zip_dict['zip'][i] = get_zip_code(bus_zip_dict['latitude'][i], bus_zip_dict['longitude'][i])

    with open("../cache/bus_zip_usa.pkl", 'wb') as fh:
        pkl.dump(bus_zip_dict, fh)

"""
    fips_zip_conversion()

Create a two-way mapping for all buses in the synthetic grid.
"""
def fips_zip_conversion():
    with open('../cache/bus_fips.pkl', 'rb') as fh:
        bf = pickle.load(fh)

    df_raw = pd.read_excel("../raw/county_to_zip.xlsx")
    all_fips = df_raw['COUNTY'].astype('int32')
    all_zip = df_raw['ZIP'].astype('int32')
    all_weights = df_raw['TOT_RATIO']

    row_num = len(all_zip)

    # create zip -> counties mapping
    zip2fips = {}
    
    for i in pd.unique(all_zip):
        idx = all_zip.index[all_zip == i]
        
        cty = all_fips[idx].tolist()
        wgt = all_weights[idx].tolist()
        
        zip2fips.update({i: (cty, wgt)})

    with open('../cache/zip2fips.pkl', 'wb') as fh:
        pickle.dump(zip2fips, fh)

    # create county -> zips mapping
    fip2zips = {}
    
    for i in pd.unique(all_fips):
        idx = all_fips.index[all_fips == i]
        
        zips = all_zip[idx].tolist()
        wgt = all_weights[idx].tolist()
        
        fip2zips.update({i: (zips, wgt)})


    with open('../cache/fip2zips.pkl', 'wb') as fh:
        pickle.dump(fip2zips, fh)

"""
    get_fips_population()

Match county population and county FIPS data to produce concise FIPS population
"""
def get_fips_population():

    cty_pop_df = pd.read_csv('../raw/county_population.csv', encoding='cp1252')
    cty_name_df = pd.read_csv('../raw/county_fips_master.csv', encoding='cp1252')

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
            print(name,state,tmp,len(tmp), tmp.to_list())
            pops[cty] = tmp.to_list()[0]
            
    cty_name_df["population"] = pops

    new_df = cty_name_df.loc[:, ["fips", "county_name", "population"]]
    new_df.to_csv("../cache/fips_population.csv", index=False)


"""
    get_zip_population()

Compute population of each ZIP code using percentage share and FIPS population
"""
def get_zip_population():
    fips_pop_df = pd.read_csv("../cache/fips_population.csv")
    
    with open("../cache/zip2fips.pkl", 'rb') as fh:
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


    with open("../cache/zip_population.pkl", 'wb') as fh:
        pkl.dump(zip_pops, fh)
