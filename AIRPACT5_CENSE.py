# -*- coding: utf-8 -*-
"""
Python functions related to the Cardiopulmonary Events from Smoke Estimator (CENSE)
project at Washington State University. This project is led by Dr. Joe Vaughn
with collaboration from Dr. Matt Kadlec from the Washington State Department
of Ecology.

Python code written by Von P. Walden, Washington State University
"""

def relativeHealthRisk(date, latitude, longitude, age, condition):
    """
    This function determines the relative increase in health risk for 
    various conditions as a function of the atmospheric PM2.5 levels.
    
    Input:
        date      - desired data for determining health risk
                     (e.g., date = pd.to_datetime('now') - pd.Timedelta('1D'))
        lat       - desired latitude (North = +, South = -)
        lon       - desired longitude (East = +, West = -)
        age       - age of person for determining health risk (0-110)
        condition - medical condition for determining health risk
                        One of the following:
                            'Non-trauma EDV CRFs [ICD 10: A00-R99]',
                            'Asthma HA CRFs [ICD 10: J45]',
                            'WA baseline incidence rate (per person per day)',
                            'Asthma symptoms onset CRFs [ICD 10: J45]',
                            'COPD Emergency HA  CRFs [ICD 10: J41,J42,J43,J44]',
                            'Acute Otitis Media (ear ache) Outpatient healthcare CRFs [ICD-10: H65, H66, H67, H68, H69, J10, J11, J12]',
                            'Respiratory Disease HA CRFs [ICD 10: J00-J99]',
                            'Ischemic heart disease EDV CRFs [ICD 10: I20-I25, I46, I49]'

    Output:
        personalRiskIncrease - Health risk increase for a person with given
                                latitude, longitude, age and condition for 
                                given date.
        riskIncrease         - Map of health risk increase for a person with 
                                desired age and condition for given date.
    """
    import numpy  as np
    import pandas as pd
    import xarray as xr
    import requests
    
    def find_WRF_pixel(latvar,lonvar,lat0,lon0):
        # Read latitude and longitude from file into numpy arrays
        # Renamed findWRFpixel from original function, naive_fast, written by Vikram Ravi.
        latvals = latvar[:]
        lonvals = lonvar[:]
        dist_sq = (latvals-lat0)**2 + (lonvals-lon0)**2
        minindex_flattened = dist_sq.argmin()  # 1D index of min element
        iy_min, ix_min = np.unravel_index(minindex_flattened, latvals.shape)
        
        return int(iy_min),int(ix_min)

    # ....Determine date.    now   = pd.to_datetime('now') - pd.Timedelta('1D')
    year  = str(date.year)
    month = str(date.month).zfill(2)
    day   = str(date.day).zfill(2)
    
    # ....Download PM2.5 forecast from AIRPACT5 for desire date.
    webAddress = 'http://lar.wsu.edu/airpact/airraid/' + year + '/' + 'PM25_24hr_' + year + month + day + '.json'
    airpact    = pd.read_json('[' + requests.get(webAddress).content.decode('utf-8')[:-2] + ']')
    r = airpact.ROW.values.max()
    c = airpact.COL.values.max()
    pm25 = airpact.RollingA24_PM25.values.reshape(c,r).T
    
    # ....Read the latitudes and longitudes for the AIRPACT forecast
    grid       = xr.open_dataset('GRIDCRO2D')
    latitudes  = np.reshape(grid.LAT.values,(258,285))
    longitudes = np.reshape(grid.LON.values,(258,285))
    ilat, ilon = find_WRF_pixel(latitudes, longitudes, latitude, longitude)

    # ....Determine row and column for desired latitude and longitude.
    
    # ....Read the Concentration Response Functions (CRFs)
    #     These were obtained from Matt Kadlec, WA Department of Ecology
    CRFs = pd.read_csv('CRFs by age year-Table 1.csv', header=[1], na_values=['NaN', '-'])
    
    # ....Calculate the relative increase in health risk
    CRF = CRFs[condition]
    riskIncrease = 1. - np.exp(-CRF[age] * pm25)
    personalRiskIncrease = riskIncrease[ilat, ilon]
    
    return personalRiskIncrease, riskIncrease, CRF

