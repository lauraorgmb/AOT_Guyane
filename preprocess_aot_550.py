# -*- coding: utf-8 -*-
"""
Created on Mon May 31 15:24:41 2021

@author: Laura Orgambide, Université Rennes 2 et Laboratoire Ecologie, Environ-
nement, Interactions des Systèmes Amazoniens, CNRS (Cayenne).

This script has been developed during an internship at Cayenne at CNRS in 2021.
The aim of the study is to create a modelisation of plumes of Saharian aerosols
flying over French Guiana.
ATMO-Guyane produced a dataset of PM10 concentration in Cayenne, Matoury and
Kourou from 2010 to 2020. A dataset of VIIRS AOT EDR 550 nm product has been
downloaded via a FTP (File Protocol Transfer) during the same period. The prin-
cipal step to modelize the presence of desertic aerosols over French Guiana is
to compare PM10 and VIIRS AOT EDR 550 nm measurements. In order to evaluate
VIIRS AOT product, via a linear regression, a coefficient of determination will
be produced and will give an idea of the quality of the VIIRS measurements.
"""

#%% IMPORTS AND FUNCTIONS

import os
import sys
import glob
import time
import datetime as dt
import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

def lons_lats():
    """
    Create longitudes and latitudes arrays of the Gridded AOT EDR 550 nm data.

    Returns
    -------
    lons : array
        ARRAY OF LONGITUDES OF SHAPE (720, 1440).
    lats : array
        ARRAY OF LATITUDES OF SHAPE (720, 1440).
    """
    resolution = 0.25
    min_lon = -179.875
    max_lon =  179.875
    min_lat =  -89.875
    max_lat =   89.875

    #create an horizontal matrix between -180 to 180 each 0.25
    h_lon = np.arange(min_lon, max_lon+0.25, resolution)
    lons = np.tile(h_lon, (720,1)) #clone column 720 times
    #create an vertical matrix between -90 to 90 each 0.25
    v_lat = np.arange(min_lat, max_lat+0.25, resolution)[:, np.newaxis]
    lats = np.tile(v_lat, (1,1440)) #clone line 1440 times

    return lons, lats

def retrieve_coordinates(lons, lats, lon_lat):
    """
    Retrieve the position x, y in the matrix for the given geographic
    coordinates.

    Parameters
    ----------
    lons : array
        ARRAY OF LONGITUDES OF THE WHOLE GRIDDED AOT EDR 550 NM PRODUCTS.
    lats : array
        ARRAY OF LATITUDES OF THE WHOLE GRIDDED AOT EDR 550 NM PRODUCTS.
    lon_lat : list
        GEOGRAPHIC LONGITUDE AND LATITUDE OF THE STATION.

    Returns
    -------
    x : int
        POSITION FOR GIVEN LATITUDE POINT IN LATS ARRAY.
    y : int
        POSITION FOR GIVEN LONGITUDE POINT IN LONS ARRAY.
    """
    #retreive possible x positions in array
    x_list = []
    for i_lat in range(lats.shape[0]-1):
        if lon_lat[1]-0.25 < lats[i_lat,0] < lon_lat[1]+0.25:
            x_list.append(i_lat)
        else:
            pass

    #retrieve possible y positions in array
    y_list = []
    for j_lon in range(lons.shape[1]-1):
        if lon_lat[0]-0.25 < lons[0,j_lon] < lon_lat[0]+0.25:
            y_list.append(j_lon)
        else:
            pass

    #select array coordinates (x,y) most closed to geographic coordinates
    for i_lat in x_list :
        x_pos = min(x_list, key=lambda i_lat:abs(i_lat-lon_lat[1]))
    for j_lon in y_list:
        y_pos = min(y_list, key=lambda j_lon:abs(j_lon-lon_lat[0]))

    return x_pos, y_pos

def open_bin(bin_file):
    """
    Open .high.bin file after being decompressed (gz) with PowerISO software.

    Parameters
    ----------

    bin_file : str
        FILENAME CONDUCTING TO THE .BIN FILE

    Returns
    -------
    aot_edr : numpy array
        ARRAY CONTAINING AOT VALUES.
    nAOT : numpy array
        ARRAY CONTAINING NUMBER OF PIXELS USED TO COMPUTE AOT MEAN IN THE GRID BOX.
    """

    #open .high.bin file
    imnp = np.fromfile(bin_file, dtype = np.single) #type : single precision float
    imnp2 = np.fromfile(bin_file, dtype = np.int_) #type : long integer
    ##NB : the data is contained in the same matrix array use to open the data
    ##NB2 : note that there are 2 073 600 cells which, divided by 2,
    ##corresponds to the 1 036 800 cells as described in the readme text.

    # reshape arrays
    array = np.reshape(imnp, [int(len(imnp)/1440), 1440])
    # array = np.flipud(array)
    array2 = np.reshape(imnp2, [int(len(imnp)/1440), 1440])
    # array2 = np.flipud(array2)
    aot = array[:720, :1440] # cut the principal arrays
    n_aot = array2[720:, :1440] # cut the principal array

    # mask invalid values
    aot = ma.masked_where(aot > 5, aot) # mask aberrant values
    aot = ma.masked_where(aot < -100, aot) # mask invalid values
    n_aot = ma.masked_where(n_aot > 100, n_aot) # mask invalid values

    return aot, n_aot

def projection(array2project, date, lons, lats, subregion):
    """
    Project AOT data over studied areas.

    Parameters
    ----------
    array2project : numpy array.
        AOT EDR 550 NM DATA.
    date : str
        DATE OF THE AOT MEASUREMENTS.
    lons : numpy array
        ARRAY OF SHAPE (720,1440) CONTAINING WORLDWIDE LONGITUDES.
    lats : numpy array
        ARRAY OF SHAPE (720,1440) CONTAINING WORLDWIDE LATITUDES.
    subregion : str
        NAME OF THE STUDIED SUBREGION.

    Returns
    -------
    None.

    """
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.OCEAN) #add ocean
    ax.add_feature(cartopy.feature.RIVERS) #add rivers
    ax.add_feature(cartopy.feature.BORDERS) #add borders
    ax.coastlines()
    levels = np.linspace(0.0, 1.0) #set limits of the colorbar
    #define colorbar and apply array to project to the plot
    cbar = plt.contourf(lons, lats, array2project, 60, levels=levels,
                        cmap = 'gist_rainbow_r', transform=ccrs.PlateCarree())

    #zoom over a specific region and set crs of the given coordinates
    if subregion in ["Plateau des Guyanes", "plateau des guyanes"]:
        ax.set_extent([-60, -50, 1, 7], crs=ccrs.PlateCarree())

    elif subregion in ["Cayenne", "cayenne"]:
        ax.set_extent([-53.25, -52.25, 4.5, 5.5], crs=ccrs.PlateCarree())

    elif subregion in ["Kourou", "kourou"]:
        ax.set_extent([-52.95, -52.35, 4.85, 5.45], crs=ccrs.PlateCarree())

    elif subregion in ["Matoury", "matoury"]:
        ax.set_extent([-52.45, -52.25, 4.85, 4.6], crs=ccrs.PlateCarree())

    elif subregion in ["Saint-Georges", "saint-georges", "St-Georges",
                       "st-georges"]:
        ax.set_extent([-52.5, -52, 3.25, 3.75], crs=ccrs.PlateCarree())

    elif subregion in ["Guadeloupe", "guadeloupe"]:
        ax.set_extent([-61, -62, 15.5, 16.5], crs=ccrs.PlateCarree())

    else:
        print("You did not select a particular region of interest.")
        print("""Please select your subregion parameter in the following list :
              \nPlateau des Guyanes, plateau des guyanes, Cayenne, cayenne,
              Kourou, kourou, Saint-Georges, saint-georges, St-Georges,
              st-georges, Guadeloupe, guadeloupe.""")

    plt.colorbar(cbar, orientation='vertical') #add colorbar to plot
    plt.title("VIIRS AOT EDR 550 nm /n in %s on %s"%(subregion, date))
    plt.savefig('Fig_%s_%s.png'%(subregion, date)) #save projection
    plt.show() #show projection

def axes3_plot(x, y_1, y_2, site):
    """
    Plot two datasets (y1,y2) according to x abscisse and save automatically
    the figure.

    Parameters
    ----------
    x : list or array
        ABSCISSE OF THE SCATTER PLOT.
    y_1 : list or array
        DATA TO PLOT.
    y_2 : list or array
        SECOND DATA TO SCATTER PLOT.
    site : str
        NAME OF THE STUDIED STATION.

    Returns
    -------
    None.

    """
    fig, ax1 = plt.subplots() #add an axis to the figure
    ax2 = ax1.twinx() #create a twin x axis
    ax1.scatter(x, y_1, c='b') #scatter first dataset
    ax2.scatter(x, y_2, c='r', marker = "x") #scatter second dataset

    ax1.set_ylabel("Gridded AOT EDR 550 nm") #define axis labels
    ax2.set_ylabel("n AOT EDR 550 nm")
    ax1.legend("AOT") #define legend
    ax2.legend("n AOT")
    ax1.tick_params(axis='x', rotation=45) #set a rotation to tick labels
    plt.title("Gridded AOT EDR at 550 nm over %s"%(site)) #set a title
    plt.show() #show figure

    fig.savefig("AOT_nAOT_%s.png"%(site)) #save figure


def linear_reg(data1, data2):
    """
    Compute linear regression.
    Source : https://www.it-swarm-fr.com/fr/python/comment-calculer-r-
    squared-avec-python-et-numpy/957531669/

    Parameters
    ----------
    data1 : numpy array
        DATASET TO PUT IN ABSCISSE AXIS.
    data2 : numpy array
        DATASET TO PUT IN ORDINATE AXIS.

    Returns
    -------
    polynomial : int
        Coefficients of the regression line.
    r_squared : int
        Determination coefficient between both datasets.

    """
    #compute polynomial coefficients and put in list
    coeffs = np.polyfit(data1, data2, 1)
    polynomial = coeffs.tolist()

    #r-squared
    p_r2 = np.poly1d(coeffs)

    #fit values, and mean
    yhat = p_r2(data1)
    ybar = np.sum(data2)/len(data2)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((data2 - ybar)**2)
    r_squared = ssreg / sstot

    return polynomial, r_squared

def scatter_plot(x_data, y_data, r_squared, coeffs_reg):
    """
    Plot x and y datasets and display results of linear regression
    (r_squared and polynomial).

    Parameters
    ----------
    x_data : numpy array
        DESCRIPTION.
    y_data : numpy array
        DESCRIPTION.
    r_squared : int
        R SQUARED COMPUTED FROM X AND Y DATASETS.
    coeffs_reg : int
        COEFFICIENTS OF REGRESSION LINE COMPUTED FROM X AND Y DATASETS.

    Returns
    -------
    None.

    """
    #scatter AOT EDR and in-situ measurements of PM10
    plt.scatter(x_data,y_data)
    plt.plot(x_data, coeffs_reg[0]*x_data+coeffs_reg[1], "g--") #add line
    plt.title("AOT EDR 550 nm correlation with PM10")
    plt.xlabel("in-situ PM10 \n R2: %s"%(str(r_squared)))
    plt.ylabel("Gridded AOT EDR 550 nm")
    plt.show() #display the graph
    plt.savefig("Fig1_scatter_plot.png") #save the figure in the new working directory

def df_concat(df1,df2,df3,df4):
    """
    Concatenate four dataframes.

    Parameters
    ----------
    df1 : dataframe
        DATAFRAME OF DATASETS 1.
    df2 : dataframe
        DATAFRAME OF DATASETS 2.
    df3 : dataframe
        DATAFRAME OF DATASETS 3.
    df4 : dataframe
        DATAFRAME OF DATASETS 4.

    Returns
    -------
    concat_df : dataframe
        RESULT OF THE CONCATENATION OF FOUR DATAFRAMES.

    """
    arrays2add = np.append(df4, df3, axis=0) #append dataframe to an array
    arrays2add  = np.append(arrays2add, df2, axis=0)
    #name array containing all AOT EDR values (from all classes)
    concat_array = np.append(arrays2add, df1, axis=0)
    concat_df = pd.DataFrame(concat_array) #convert the array to dataframe

    return concat_df

#%% DEFINITIONS OF VARIABLES

if __name__ == "__main__":

    #define parameters for AOT study
    YEAR = "2020"
    classes = ["AOT_sup_150", "AOT_sup_100", "AOT_inf_50", "AOT_sup_50"]
    stations_ATMO = ["Cayenne"]
    stations_AERONET = ["Paracou","Nouragues","Surinam"]

    #set working directories
    WD = r'C:\Users\cath\Desktop\Laura\Données'
    VIIRS_WD = os.path.join(WD, r'AOT_VIIRS\Gridded_AOTEDR\%s'%(YEAR))
    result_dir = os.path.join(VIIRS_WD, r'%s\Results'%(YEAR))
    PM10_WD = os.path.join(WD, 'PM10')

    #define the area of interest according to the name of the station given
    for station in stations_ATMO:
        if station == "Cayenne":
            geo_lonlat = [-52.310793, 4.933351]
            COL_NB = 3
        elif station == "Kourou":
            geo_lonlat = [-52.645549, 5.163618]
            COL_NB = 4
        elif station == "Matoury":
            geo_lonlat = [-52.323817, 4.848194]
            COL_NB = 5
        else:
            print("The name of the station is probably wrong.")

    #predefined objects
    LONS, LATS = lons_lats() #create lon and lat information arrays
    #decoment following lines to resampple arrays
    # LONS = cv2.resize(LONS, (2880,1440), interpolation = cv2.INTER_NEAREST)
    # LATS = cv2.resize(LATS, (2880,1440), interpolation = cv2.INTER_NEAREST)
    x_arr, y_arr = retrieve_coordinates(LONS, LATS, geo_lonlat)
    print("""Latitude %s corresponds to x = %s in the matrix."""
         """ Longitude %s corresponds to y = %s in the matrix."""
          %(geo_lonlat[1], x_arr, geo_lonlat[0], y_arr))

#%% 1ST PART : EXTRACTION OF AOT, SCATTER AOT AND NAOT, PROJECT DATA

start = time.time() #start the timer to compute algorithm execution time

for station in stations_ATMO:
    #create a Dataframe containing the date, the aot_edr value and n_AOT
    df_aot_class = pd.DataFrame(columns = ["Date", "AOT EDR 550 nm", "nAOT"])

    for c in classes:
        #set direction to the folder
        class_dir = os.path.join(VIIRS_WD, c)
        if not os.path.exists(class_dir):
            os.makedirs(class_dir)
        os.chdir(class_dir)
        f_list = glob.glob("*.bin") #set a list of files to open
        #define the folder where the results will be saved
        result_folder = os.path.join(class_dir, "Results")

        #append to dataframe date, AOT and nAOT values
        for filename in f_list:
            #get file's date format in datetime and in string
            file_date = dt.datetime.strptime(filename[28:36], "%Y%m%d").date()
            file_datestr = dt.datetime.strftime(file_date, "%Y-%m-%d")

            if file_datestr[5:7] in ["01", "02", "03", "04", "05", "06", "06",
                                     "07", "08", "09", "10", "11", "12"]:
                try:
                    aot_edr, n_aot_edr = open_bin(filename) #open high.bin file
                    # aot_edr = cv2.resize(aot_edr, (2880,1440),
                    #                      interpolation = cv2.INTER_NEAREST)
                    # n_aot_edr = cv2.resize(n_aot_edr, (2880,1440),
                    #                        interpolation = cv2.INTER_NEAREST)
                    # aot_edr = ma.masked_where(aot_edr == -9999, aot_edr)
                    # n_aot_edr = ma.masked_where(n_aot_edr < -100, n_aot_edr)

                    #research date, AOT and nAOT values for the second pair of
                    #coordinates given by retrieve_PointCoordinates function
                    df_aot_class = df_aot_class.append({"Date" : file_date,
                                    "AOT EDR 550 nm" :
                                        np.ma.mean(aot_edr[x_arr,y_arr]),
                                    "nAOT" :
                                        np.ma.mean(n_aot_edr[x_arr,y_arr])},
                                   ignore_index = True)
                    #project AOT data over area of interest
                    os.chdir(result_folder)
                    projection(aot_edr, file_datestr, LONS,
                               LATS, subregion = "Plateau des Guyanes")

                except ValueError:
                    print("The image of %s is corrupted."%(file_datestr))


        df_aot_class.to_excel("aot_edr_df_%s_%s.xlsx"%(station, c))
        #after saving the dataframe, you have to open it in excel and convert
        ##it to the right format (put write separator : ",") if it is not
        ##In excel : go to Data>Convert and choose "," separator

    #plot 3 axes
    np_aot_class = df_aot_class.to_numpy() # convert dataframe to numpy array
    axes3_plot(np_aot_class[:,0], np_aot_class[:,1], np_aot_class[:,2],
               station)
    del df_aot_class, aot_edr, n_aot_edr, f_list

end = time.time()

print("The execution time of the proccess spent : ", end-start)

#%% 2ND PART : SCATTER AOT AND NAOT FOR 4 CLASSES AND SCATTER AOT
#              AND ATMO GUYANE STATION MEASUREMENTS

#open the dataframe containing mean PM10 for each station from 2010 to 2020
os.chdir(PM10_WD)
list_PM10 = glob.glob("seven_low_season.xlsx")
for PM10 in list_PM10:
    mean_PM10 = pd.read_excel(PM10)
    for station in stations_ATMO:
        D = 0 #counter
        for c in classes:
            D +=1
            os.chdir(os.path.join(VIIRS_WD, c))
            path_name = glob.glob("aot_edr_df_%s*.xlsx"%(station))
            globals()['df%s'%(D)] = pd.read_excel(path_name[0])

        #concat all the VIIRS AOT EDR data into one array for each station
        if len(classes) == 4:
            globals()['concat_%s'%(station)] = df_concat(
                df1,df2,df3,df4).to_numpy()
            globals()['total_AOT_EDR_%s'%(station)] = pd.DataFrame(
                globals()['concat_%s'%(station)],
                columns = ['Days', 'Date', 'VIIRS AOT', 'nAOT'])
            globals()['total_AOT_EDR_%s'%(station)].to_excel(
                "total_AOT_EDR_%s.xlsx"%(station)) #save the dataframe
            globals()['total_AOT_EDR_%s'%(station)].drop_duplicates(
                subset = "Date", keep = 'first', inplace = True)
            globals()['concat_%s'%(station)] = globals(
                )['total_AOT_EDR_%s'%(station)].to_numpy()
        else:
            sys.exit("The following part of the script need to be run with 4 \
                  classes in order to have the most data included in the \
                      computation of r squared. Please add more classes in \
                          the parameters.")

        #create new dataframe cointaining VIIRS AOT and AMTO-Guyane PM10 values
        df = pd.DataFrame(columns = ["Date", "AOT", "PM10", "nAOT"])
        mean_PM10 = mean_PM10.to_numpy()

        #create dataframe containing corresponding PM10 and AOT EDR 550 nm
        for a in range(globals()['concat_%s'%(station)].shape[0]):
            for b in range(mean_PM10.shape[0]):
                if mean_PM10[b,1].date() == globals()['concat_%s'%(
                        station)][a,1].date() and mean_PM10[b,COL_NB] >= 50:
                    aot_pm10 = df.append({"Date" : mean_PM10[b,1].date(),
                                    "AOT" : globals()['concat_%s'%(
                                        station)][a,2],
                                    "PM10" : mean_PM10[b, COL_NB],
                                    "nAOT" :  globals()['concat_%s'%(
                                        station)][a,3]},
                                   ignore_index = True)
                else: #images dates do not correspond to PM10 date
                    pass

        df_nAOT_r2 = pd.DataFrame(columns=["nAOT", "r2", "nb_img", "tot_img"])

        #plot the previous dataframe
        aot_pm10.drop_duplicates(subset = "Date",
                                  keep = 'first',
                                  inplace = True)
        tot_img = aot_pm10.shape[0]
        print("Nombre de jours où PM10 > 50 µg.m-3 : ", tot_img)

        #delete nan values
        aot_pm10 = aot_pm10.replace('--', np.nan)
        aot_pm10 = aot_pm10.dropna()

        quarter_list = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5,
        2.75, 3, 3.25, 3.5, 3.75, 4]

        for quarter_n_aot in quarter_list:
            aot_pm10 = aot_pm10.drop(
                aot_pm10[aot_pm10.nAOT <= quarter_n_aot].index)
            if aot_pm10.shape[0] != 0:
                array2scatter = aot_pm10.to_numpy() #convert to numpy array
                #define PM10 data
                PM10 = np.asarray(array2scatter[:,2], dtype= np.float)
                #define VIIRS AOT EDR data
                viirs_aot = np.asarray(array2scatter[:,1], dtype= np.float)
                #compute coefficient of regression line and r squared
                coeff_line, r2 = linear_reg(PM10, viirs_aot)
                #stock values in a dataframe
                df_nAOT_r2 = df_nAOT_r2.append({"nAOT": quarter_n_aot,
                                                "r2" : r2,
                                                "nb_img": aot_pm10.shape[0],
                                                "tot_img": tot_img},
                                               ignore_index = True)
            else: #break the loop when there is no more images
                break

        os.chdir(WD)
        df_nAOT_r2.to_excel("list_R2_%s_%s.xlsx"%(YEAR, station))
        scatter_plot(PM10, viirs_aot, r2, coeff_line)
