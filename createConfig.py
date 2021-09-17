#!/usr/bin/env python
from __future__ import print_function, division
from argparse import ArgumentParser
import numpy as np
from datetime import datetime, timedelta
from skyfield.api import Topos, load, EarthSatellite
from astropy.io import fits
from astropy.wcs import WCS
import yaml
from tqdm import tqdm
import spacetracktool as st
from spacetracktool import operations as ops
import matplotlib.pyplot as plt


def dayOfYear(utc):
    """
    converts utc to day of year
    """
    time_delta = utc - datetime.strptime(str(utc.year)+"-01-01 00:00:00", '%Y-%m-%d %H:%M:%S')
    return time_delta.total_seconds()/(24.0*60*60) +1
    

def tle_query(utc, args):
    """
    returns the historic tle for the object
    """
    query = st.SpaceTrackClient(args.user, args.passwd)
    start_date = utc + timedelta(hours=-24*args.pastDaysRange)
    end_date = utc 
    date_range = ops.make_range_string(str(start_date.year) + "-" + str(start_date.month).zfill(2) + "-"+ str(start_date.day).zfill(2),
        str(end_date.year) + "-" + str(end_date.month).zfill(2) + "-"+ str(end_date.day).zfill(2))
    result = query.tle_query(epoch=date_range, norad_cat_id=args.noradid)
    return result


def getIntialGuess(result, obs_dayofyear):

    ## get historic data
    e_array, i_array, ra_array = [], [], []
    aop_array, ma_array, mm_array = [], [], []
    doy_array, line1_array, line2_array = [], [], []
    date_array = []

    for i in tqdm(range(len(result.json()[:]))):
        line1 = result.json()[i]["TLE_LINE1"]
        line2 = result.json()[i]["TLE_LINE2"]
        utc = datetime.strptime(result.json()[i]["EPOCH"], '%Y-%m-%d %H:%M:%S')
        days = dayOfYear(utc)
        date_array.append(days)
        line1_array.append(line1)
        line2_array.append(line2)
        e = float(line2[26:33])/10000000
        i =  float(line2[8:16])
        ra = float(line2[17:25])
        aop = float(line2[34:42])
        ma = float(line2[43:51])
        mm = float(line2[52:63])
        doy = float(line1[20:32])
        e_array.append(e)
        i_array.append(i)
        ra_array.append(ra)
        aop_array.append(aop)
        ma_array.append(ma)
        mm_array.append(mm)
        doy_array.append(doy)

    ## typecast to numpy arrays
    e_array = np.array(e_array)
    i_array = np.array(i_array)
    ra_array = np.array(ra_array)
    aop_array = np.array(aop_array)
    ma_array = np.array(ma_array)
    mm_array = np.array(mm_array)
    doy_array = np.array(doy_array)
    date_array = np.array(date_array)

    ## unwrap values
    ra_array = np.degrees(np.unwrap(np.radians(ra_array), discont=np.pi))
    aop_array = np.degrees(np.unwrap(np.radians(aop_array), discont=np.pi))
    i_array = np.degrees(np.unwrap(np.radians(i_array), discont=np.pi))

    ## get the most likely change in value per day (used for generating boundary conditions)
    e_per_day = np.nanmedian(np.divide(np.abs(e_array[1:] - e_array[:-1]), np.abs(doy_array[1:] - doy_array[:-1])))
    i_per_day = np.nanmedian(np.divide(np.abs(i_array[1:] - i_array[:-1]), np.abs(doy_array[1:] - doy_array[:-1])))
    ra_per_day = np.nanmedian(np.divide(np.abs(ra_array[1:] - ra_array[:-1]), np.abs(doy_array[1:] - doy_array[:-1])))
    aop_per_day = np.nanmedian(np.divide(np.abs(aop_array[1:] - aop_array[:-1]), np.abs(doy_array[1:] - doy_array[:-1])))
    mm_per_day = np.nanmedian(np.divide(np.abs(mm_array[1:] - mm_array[:-1]), np.abs(doy_array[1:] - doy_array[:-1])))

    ## get the most recent tle value and get intial guess
    ref_dates = date_array - obs_dayofyear
    ## since we only query the tles from the past, the max ref date is the intial gues
    guess_arg = int(np.where(ref_dates == np.max(ref_dates))[0])
    
    ## intial guess
    guess_i = i_array[guess_arg]
    guess_ra = ra_array[guess_arg]
    guess_e = e_array[guess_arg]
    guess_aop = aop_array[guess_arg] 
    guess_ma = ma_array[guess_arg]
    guess_mm = mm_array[guess_arg]
    guess_line1 = line1_array[guess_arg]
    guess_line2 = line2_array[guess_arg]
    guess_date = date_array[guess_arg]
    
    ## wrap the values to keep them between 0 and 360
    while guess_ra > 360:
        guess_ra -= 360
    
    while guess_ra < 0:
        guess_ra += 360

    while guess_aop > 360:
        guess_aop -= 360
    
    while guess_aop < 0:
        guess_aop += 360

    ## time b/w intial guess and observation
    time_diff = abs(guess_date- obs_dayofyear)

    ## create output yaml config
    yaml_file = [{'initialGuess': {'i': float(guess_i),
                                   'ra': float(guess_ra),
                                   'e': float(guess_e),
                                   'aop': float(guess_aop),
                                   'mm': float(guess_mm)}}, 
                {'maxLimits': { 'i': float(guess_i + i_per_day*time_diff),
                                'ra': float(guess_ra + ra_per_day*time_diff),
                                'e': float(guess_e + e_per_day*time_diff),
                                'aop': float(guess_aop + aop_per_day*time_diff),
                                'mm': float(guess_mm + mm_per_day*time_diff)}},
                 {'minLimits': { 'i': float(guess_i - i_per_day*time_diff),
                                'ra': float(guess_ra - ra_per_day*time_diff),
                                'e': float(guess_e - e_per_day*time_diff),
                                'aop': float(guess_aop - aop_per_day*time_diff),
                                'mm': float(guess_mm - mm_per_day*time_diff)}}, 
                {'niter': {'major': 300, 'minor': 80, 'step': 0.1}},
                {'tle': {'line2': guess_line1, 'line3': guess_line2}}]

    return yaml_file


def main(args):

    ## obtain the wcs object and the datetime of the observation
    hdu = fits.open(args.wcsFile)
    obs_date = datetime.strptime(hdu[0].header["DATE-OBS"],'%Y-%m-%dT%H:%M:%S.%f')
    obs_dayofyear = dayOfYear(obs_date)

    ## perform the tle query
    result = tle_query(obs_date,  args)

    yaml_file = getIntialGuess(result, obs_dayofyear)

    ## write to file
    with open("auto_created_config" + str(args.noradid) + ".yaml" , "w") as output_file:
        yaml.dump(yaml_file, output_file)




if __name__ == "__main__":
    parser = ArgumentParser("create_config_file", description="creates the config file for orbit determination")
    parser.add_argument("--noradid", type=int, required=True, help="the norad id of the satellite")
    parser.add_argument("--wcsFile", required=True, help="a dummy file that has the wcs object and the date of the observation")
    parser.add_argument("--pastDaysRange", type=int, default=60, help="the number of past days to search for historic data")
    parser.add_argument("--debug", type=bool, default=False, help="runs the script in debug mode")
    parser.add_argument("--user", required=True, help="the user name of space-track.org account")
    parser.add_argument("--passwd", required=True, help="the password of space-track.org account")
    args = parser.parse_args()

    global debug
    debug = args.debug

    if debug:
        print("running in debug mode")

    main(args)
