# -*- coding: utf-8 -*-
"""
Created on Thu May 27 09:29:19 2021

@author: Laura Orgambide

Desription : The following script is meant to download Gridded AOT EDR data
            from VIIRS sensor on board of the Suomi-NPP satelitte. The data is
            accessible from an FTP server : 
            ftp://ftp.star.nesdis.noaa.gov/pub/smcd/VIIRS_Aerosol/viirs_aod_gridded/idps/snpp/edraot550 .
            They are saved into 4 folders corresponding to four classes of PM10
            thresholds ([0:50[ // [50:100[ // [100:150[ // [150:] mg.m-3). 
"""

import os, shutil
import pandas as pd
import numpy as np
import urllib.request as request
from contextlib import closing
from datetime import date, datetime, timedelta

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
    df : dataframe
        DATAFRAME CONTAINING ALL DATES OF THE PERIOD OF INTEREST.

    """
    df = pd.date_range(date_start,date_end-timedelta(days=1),freq='d')
    
    return df
    

#%%Batch download

class_AOT = ["AOT_inf_50", "AOT_sup_50", "AOT_sup_100", "AOT_sup_150"] # create beforehand 4 folders to sort AOT data according to PM10 values

i = 0 # counter for images not available in the FTP sever
    
for folder in class_AOT:
    
    #set working directory
    wd = r"C:\Users\cath\Desktop\Laura\Données\AOT_VIIRS\Gridded_AOTEDR"
    wd = os.path.join(wd, folder)
    os.chdir(wd)
    
    #open table for looking for when PM10 was < 50 mg.m-3
    tab_path = r"C:\Users\cath\Desktop\Laura\Données\PM10\tableau_pics_mean_%s.xlsx"%(folder)
    df = pd.read_excel(tab_path)
    np = df.to_numpy()
    for date in np[:,1]:
        date_str = date.strftime("%Y%m%d") # convert datetime format to date string
        
        #set url from where to download files 
        url = "ftp://ftp.star.nesdis.noaa.gov/pub/smcd/VIIRS_Aerosol/viirs_aod_gridded/idps/snpp/edraot550/"
        file = "npp_aot550_edr_gridded_0.25_%s.high.bin.gz"%date_str
        link = url + date_str[:4] + "/" + file
       
        try : 
            with closing(request.urlopen(link)) as r:
                with open(file, 'wb') as f:
                    """
                    Create condition of cloud cover and data availability over ROI
                    """
                    # shutil.copyfileobj(r, f)
                    # print("An image was found for this date : %s"%(date_str))
                    continue
        except:
            i +=1
            # print("No image were found for this date : %s. \n Image number %s."%(date_str,i))
            continue
        
print(i)
            
#%%Download unique tile for precise date
    
#set working directory
wd = r"C:\Users\cath\Desktop\Laura\Données\AOT_VIIRS\Gridded_AOTEDR"
folder = "2012_all"
wd = os.path.join(wd, folder)
os.chdir(wd)

#give the dates of the tiles to download as str
# date_list = ["20120502", "20120512", "20120522", "20120601", "20120611", "20120621", "20120701", "20120711", "20120721", "20120731"]
df = date_range(date(2012,1,1), date(2012,12,31)) # create a range of datetime
date_list = df.to_numpy()


for i in range(df.shape[0]):
    
    date_str = df[i].strftime("%Y%m%d") # convert datetime format to date string
    # print(date_str)
    #set url from where to download files
    url = "ftp://ftp.star.nesdis.noaa.gov/pub/smcd/VIIRS_Aerosol/viirs_aod_gridded/idps/snpp/edraot550/"
    file = "npp_aot550_edr_gridded_0.25_%s.high.bin.gz"%date_str
    link = url + date_str[:4] + "/" + file
    
    try: 
        with closing(request.urlopen(link)) as r:
            with open(file, 'wb') as f:
                """
                Create condition of cloud cover and data availability over ROI
                """
                shutil.copyfileobj(r, f)
                print("An image was found for this date : %s"%(date_str))
    except:
        print("No image were found for this date : %s"%(date_str))
        pass
    