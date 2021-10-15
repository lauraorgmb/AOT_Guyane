# -*- coding: utf-8 -*-
"""
Created on Thu May 27 09:29:19 2021

@author: Laura Orgambide

Desription : The following script is meant to download Gridded AOT EDR data
            from VIIRS sensor on board of the Suomi-NPP satelitte. The data is
            accessible from an FTP server :
            ftp://ftp.star.nesdis.noaa.gov/pub/smcd/VIIRS_Aerosol/viirs_aod_
            gridded/idps/snpp/edraot550 .
            They are saved into 4 folders corresponding to four classes of PM10
            thresholds ([0:50[ // [50:100[ // [100:150[ // [150:] mg.m-3).
            To run the script, don't forget to donwload PM10 dataframe!
"""

import os
import shutil
from datetime import date, timedelta
import urllib.request as request
from urllib.error import URLError
from contextlib import closing
import pandas as pd

def date_range(date_start, date_end):
    """
    Create a range of datetime with 1 day of interval.

    Parameters
    ----------
    date_start : datetime.date
        START OF THE PERIOD OF INTEREST.
    date_end : datetime.date
        END OF THE PERIOD OF INTEREST.

    Returns
    -------
    df_dates : dataframe
        DATAFRAME CONTAINING ALL DATES OF THE PERIOD OF INTEREST.

    """
    df_dates = pd.date_range(date_start,date_end-timedelta(days=1),freq='d')

    return df_dates


#%%Batch download

if __name__ == "__main__":
    
    #set working directory
    WD = os.getcwd()
    VIIRS_WD = os.path.join(WD, r"Guyane\Gridded_AOTEDR")
    if not os.path.exists(VIIRS_WD):
        os.makedirs(VIIRS_WD)

    #set url from where to download files
    PRE_URL = "ftp://ftp.star.nesdis.noaa.gov/pub/smcd/VIIRS_Aerosol"
    URL = os.path.join(PRE_URL, "viirs_aod_gridded/idps/snpp/edraot550")

    #create beforehand 4 folders to sort AOT data according to PM10 values
    class_AOT = ["AOT_inf_50", "AOT_sup_50", "AOT_sup_100", "AOT_sup_150"]

    #sort downloaded images according to PM10 concentration
    for folder in class_AOT:
        #set class working directory
        WD_class = os.path.join(VIIRS_WD, folder)
        if not os.path.exists(WD_class):
            os.makedirs(WD_class)
        os.chdir(WD_class)

        #set PM10 dataframe path
        TAB_PATH = os.path.join(WD, r"Guyane\tableau_pics_mean_%s.xlsx"%(folder))

        #open table to obtain PM10 concentration
        PM10_df = pd.read_excel(TAB_PATH)
        PM10_array = PM10_df.to_numpy()
        for date_PM10 in PM10_array[:,1]:
            date_str = date_PM10.strftime("%Y%m%d") #convert to string format

            #set url from where to download files
            file = "npp_aot550_edr_gridded_0.25_%s.high.bin.gz"%date_str
            link = URL + date_str[:4] + "/" + file

            try :
                with closing(request.urlopen(link)) as r:
                    with open(file, 'wb') as f:
                        shutil.copyfileobj(r, f)
                        print("An image was found on %s"%(date_str))
                        continue
            except URLError:
                print("No image were found on %s."%(date_str))
                continue

#%%Download unique tile for precise date

    #set working directory
    FOLDER = "2012_all"
    WD_YEAR = os.path.join(VIIRS_WD, FOLDER)
    if not os.path.exists(WD_YEAR):
        os.makedirs(WD_YEAR)
    os.chdir(WD_YEAR)

    ##give the dates of the tiles to download as str
    #create a range of datetime
    date_df = date_range(date(2012,1,1), date(2012,12,31))

    #donwload files according to date of interest
    for i in range(date_df.shape[0]):

        date_str = date_df[i].strftime("%Y%m%d") #convert to string format
        file = "npp_aot550_edr_gridded_0.25_%s.high.bin.gz"%date_str
        link = URL + date_str[:4] + "/" + file

        try:
            with closing(request.urlopen(link)) as r:
                with open(file, 'wb') as f:
                    shutil.copyfileobj(r, f)
                    print("An image was found for this date : %s"%(date_str))
        except URLError:
            print("No image were found on %s"%(date_str))
