#!/usr/bin/env python
from __future__ import division, print_function
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from argparse import ArgumentParser
from tqdm import tqdm
from astropy.nddata import Cutout2D
from datetime import datetime, timedelta
from spacetracktool import operations as ops
from skyfield.api import EarthSatellite
from skyfield.api import Topos, load
import matplotlib.pyplot as plt
import spacetracktool as st
import scipy.optimize as opt
from os import path
import json
import csv
import sys


def obtainTimeSteps(noradid):
    """
    returns the timeSteps to search for signal
    by searching the json file made by satSearch.py
    Parameters
    ----------
    noradid         : the norad id
    Returns
    -------
    headTimeSteps   : the list of head timeSteps
    tailTimeSteps   : the list of tail timeSteps
    midTimeSteps    : the list of mid timeSteps  
    """
    with open("filtered_summary.json") as f:
        filtered_summary = json.load(f)

    for line in filtered_summary:
        if filtered_summary[line]['norad'] == str(noradid):
            timeSteps = list(map(int,filtered_summary[line]["timeSteps"].split("-")))
        
    return timeSteps, timeSteps, timeSteps

def twoD_Gaussian(data_tuple, amplitude, xo, yo, sigma):
    """
    model of a 2D circular beam
    Parameters
    ----------
    data_tupe   : the x and y mesh grid
    amplitude   : the amp at the peak
    xo          : the peak x location
    yo          : the peak y location
    sigma       : the variance of circular beam
    Returns
    -------
    returns a linear array of the beam
    """

    (x, y) = data_tuple
    xo = float(xo)
    yo = float(yo)    
    a = (1**2)/(2*sigma**2) 
    b =  0
    c =  (1)/(2*sigma**2)
    g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()


def obtainTLE(noradid):
    """
    for a given norad id, obtains the tle for the relevant date range

    Parameters
    ----------
    noradid     : the norad id of the satellite

    Returns
    -------
    returns the tle of the satellite (tle = line1 + line2 + line3)
    """

    time2 = refUTC + timedelta(hours=24)
    day1, month1, year1 = str(refUTC.day).zfill(2), str(refUTC.month).zfill(2), str(refUTC.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)

    if path.exists("{}.txt".format(noradid)) == False:

        if debug:
            print("requesting file from server")

        result = query.tle_query(epoch=date_range, norad_cat_id=noradid)

        with open("{}.txt".format(noradid), "w") as outfile:
            json.dump(result.json(), outfile)

        line1 = result.json()[0]['OBJECT_NAME']
        line2 = result.json()[0]['TLE_LINE1']
        line3 = result.json()[0]['TLE_LINE2']

    else:

        if debug:
            print('TLE file found on disk. Not downloading.')

        with open("{}.txt".format(noradid)) as json_file:
            result = json.load(json_file)

        line1 = result[0]['OBJECT_NAME']
        line2 = result[0]['TLE_LINE1']
        line3 = result[0]['TLE_LINE2']

    return line1, line2, line3

def getSatMWA(line1, line2, line3):
    """
    creates a satellite object and observer object
    
    Paramters
    ---------
    line1   : tle line 1
    line2   : tle line 2
    line3   : tle line 3
    Returns
    -------
    satellite   : the sat object
    mwa         : the observer object
    ts          : the time object
    """

    ts = load.timescale(builtin=True)
    satellite = EarthSatellite(line2, line3, line1, ts)
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827)

    return satellite, mwa, ts

def getHeadPeaks(headTimeSteps, noradid, mwa, sat, ts, args):
    """
    finds all the peak pixels in head using gaussian fit

    Parameters
    ----------
    headTimeSteps   : the list of head time steps
    noradid         : the norad id of the satellite
    mwa             : observer object
    sat             : the satellite object
    ts              : time object

    Returns
    -------
    x_array         : array of x pixels of head
    y_array         : array of y pixels of head
    ra_array        : array of ra of head
    dec_array       : array of dec of head
    """
    x_array, y_array = [], []
    x_err_array, y_err_array = [], []

    y = np.linspace(0, 4, 5)
    x = np.linspace(0, 4, 5)
    x, y = np.meshgrid(x, y)

    for t in headTimeSteps:
        hdu = fits.open("testData/"+args.filePrefix + "-t" + str(t).zfill(4) + ".fits")
        UTCTime = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        data = hdu[0].data
        data = maskData(data, UTCTime, mwa, sat, wcs, ts)

        if np.all(data == 0):
            print("file full of zeros")
            print("aborting....")
            sys.exit(0)
        
        row, col = np.where(data == np.max(data))
        beam_cutout = Cutout2D(beam, (col, row), (5 ,5 ))
        cutout = Cutout2D(data, (col, row), (5,5))
        temp = cutout.data / beam_cutout.data
        temp /= np.nansum(temp)
        
        initial_guess = (temp[1,1], 1, 1, 2)
        popt, pconv = opt.curve_fit(twoD_Gaussian, (x, y), temp.ravel(), p0=initial_guess)
        perr = np.sqrt(np.diag(pconv))
        mu_x, mu_y = popt[1], popt[2]
        sigma_x, sigma_y = perr[1], perr[2]
        solX = float(col + (mu_x - 2))
        solY = float(row + (mu_y - 2))
        
        x_array.append(solX)
        y_array.append(solY)
        x_err_array.append(sigma_x)
        y_err_array.append(sigma_y)

    return x_array, y_array, x_err_array, y_err_array


def getTailPeaks(tailTimeSteps, noradid, mwa, sat, ts, args):
    """
    finds all the peak pixels in tail using gaussian fit

    Parameters
    ----------
    tailTimeSteps   : the list of tail time steps
    noradid         : the norad id of the satellite
    mwa             : observer object
    sat             : the satellite object
    ts              : time object

    Returns
    -------
    x_array         : array of x pixels of tail
    y_array         : array of y pixels of tail
    ra_array        : array of ra of tail
    dec_array       : array of dec of tail
    """
    x_array, y_array = [], []
    x_err_array, y_err_array = [], []

    y = np.linspace(0, 4, 5)
    x = np.linspace(0, 4, 5)
    x, y = np.meshgrid(x, y)

    for t in tailTimeSteps:
        hdu = fits.open("testData/"+"Neg" + args.filePrefix + "-t" + str(t).zfill(4) + ".fits")
        UTCTime = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        data = hdu[0].data
        data = maskData(data, UTCTime, mwa, sat, wcs, ts)

        if np.all(data == 0):
            print("file full of zeros")
            print("aborting....")
            sys.exit(0)
        
        row, col = np.where(data == np.min(data))
        beam_cutout = Cutout2D(beam, (col, row), (5 ,5 ))
        cutout = Cutout2D(data, (col, row), (5,5))
        temp = cutout.data / beam_cutout.data
        temp /= np.nansum(temp)
        
        initial_guess = (temp[1,1], 1, 1, 2)
        popt, pconv = opt.curve_fit(twoD_Gaussian, (x, y), temp.ravel(), p0=initial_guess)
        perr = np.sqrt(np.diag(pconv))
        mu_x, mu_y = popt[1], popt[2]
        sigma_x, sigma_y = perr[1], perr[2]
        solX = float(col + (mu_x - 2))
        solY = float(row + (mu_y - 2))
        
        x_array.append(solX)
        y_array.append(solY)
        x_err_array.append(sigma_x)
        y_err_array.append(sigma_y)

    return x_array, y_array, x_err_array, y_err_array

        




def maskData(data,  UT, mwa, sat, wcs, ts):
    """
    this functions masks everying outside a 10 km cone from satellite
    Parameters
    ----------
    data    : the output data from RFISeeker
    UT      : the UTC of the timeStep
    mwa     : the observer object
    sat     : the satellite object
    wcs     : world coordinate system
    ts      : timescale object (required by skyfield)

    Returns
    -------
    masked version of the input data variable
    """

    ## propogate the satellite to the utc
    time = ts.utc(UT.year, UT.month,  UT.day, UT.hour, UT.minute, UT.second)
    sat.at(time)
    difference = sat - mwa
    topocentric = difference.at(time)
    ra, dec, distance = topocentric.radec()
    ra = np.degrees(ra.radians)
    dec = np.degrees(dec.radians)
    px, py = wcs.all_world2pix([[ra, dec]], 1)[0]
    LOS_range = distance.m
    
    ## mask the data
    radius = np.degrees(70000/LOS_range) # cone radius (s = r * theta)
    number_of_pixels = radius/imgScale ### check this line

    y = np.linspace(0, (imgSize - 1), imgSize)
    x = np.linspace(0, (imgSize - 1), imgSize)
    x, y = np.meshgrid(x, y)

    array_radius = np.sqrt((x - px)**2 + (y-py)**2 )
    array_radius[array_radius > number_of_pixels] = -10
    array_radius[array_radius != -10] = 1
    array_radius[array_radius != 1] = 0
   
    return data * array_radius


def fitLine(head_x, head_y, tail_x, tail_y,order=2):

    combined_x = head_x + tail_x
    combined_y = head_y + tail_y
    
    if motion == "E-W":
        ## check if x is going min-max or max-min

        if head_x[0] - head_x[-1] > 0:
            # motion x max to min
            # sort the values max to min
            combined_y = [x for _,x in sorted(zip(combined_x, combined_y),reverse=True)]
            combined_x.sort(reverse=True)
            
        elif head_x[0] - head_x[-1] < 0:
            # motion x min to max
            # sort the values min to max
            combined_y = [x for _,x in sorted(zip(combined_x, combined_y))]
            combined_x.sort()
        
        else:
            print("bug")
            print("cannot classify motion as min to max for E-W")
            print("aborting...")
            sys.exit(0)

        x_array = np.linspace(combined_x[0], combined_x[-1], 100)
        combined_x = np.array(combined_x)
        combined_y = np.array(combined_y)
        z = np.polyfit(combined_x.ravel(), combined_y.ravel(), order)
        f = np.poly1d(z)
        output = [x_array, f(x_array), f]
    
    if motion == "N-S":
        
        ## check if y is going min-max or max-min

        if head_y[0] - head_y[-1] > 0:
            #motion y max to min
            # sort the values max to min
            combined_x = [x for _,x in sorted(zip(combined_y, combined_x), reverse=True)]
            combined_y.sort(reverse=True)

        elif head_y[0] - head_y[-1] < 0:
            # motion y min to max
            # sort values min to max
            combined_x = [x for _,x in sorted(zip(combined_y, combined_x))]
            combined_y.sort()

        else:
            print("bug")
            print("cannot classify motion as min to max for N-W")
            print("aborting")
            sys.exit(0)


        y_array = np.linspace(combined_y[0], combined_y[-1], 100)
        combined_x = np.array(combined_x)
        combined_y = np.array(combined_y)
        z = np.polyfit(combined_y.ravel(), combined_x.ravel(), order)
        f = np.poly1d(z)
        output = [f(y_array), y_array, f]

    return output

def UTC2s(time):
    """
    converts utc to s with respect to refUTC
    
    Paramters
    ---------
    time    : time to convert into ref seconds
    Returns
    -------
    time converted to ref sec
    """
    
    return (time - refUTC).total_seconds()


def getMidPoints(midTimeSteps, mwa, sat, ts, x_fit, y_fit, function, args):
    
    x_array, y_array, time_array = [], [], []

    for t in midTimeSteps:
        hdu_head = fits.open("testData/"+args.filePrefix + "-t" + str(t).zfill(4) + ".fits")
        hdu_tail = fits.open("testData/"+"Neg" + args.filePrefix + "-t" + str(t).zfill(4) + ".fits")
        head_data = hdu_head[0].data
        head_data /= np.abs(np.sum(head_data))
        tail_data = hdu_tail[0].data
        tail_data /= np.abs(np.sum(tail_data))

        UTCTime = datetime.strptime(hdu_head[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        time_array.append(UTC2s(UTCTime))

        head_data = maskData(head_data, UTCTime, mwa, sat, wcs, ts)
        tail_data = maskData(tail_data, UTCTime, mwa, sat, wcs, ts)

        streak = head_data + tail_data

        ## iterate throug fit and find boundarys
        fit_val = []
        #x_prev = 
        x_inv, y_inv = [], []
        for x, y in zip(x_fit, y_fit):
            fit_val.append(streak[int(y),int(x)])

        if motion == "E-W":
            min_arg = np.where(fit_val == np.min(fit_val))[0]
            max_arg = np.where(fit_val == np.max(fit_val))[0]


            ### check dimension
            if float(min_arg.shape[0]) > 1:
                min_arg = min_arg[0]
            elif float(max_arg.shape[0]) >1:
                max_arg = max_arg[0]

            x_min = x_fit[min_arg]
            x_max = x_fit[max_arg]

            x_streak = np.linspace(x_min, x_max, 100)
            y_streak = function(x_streak)


            fit_val = []
            #x_inv, y_inv = [], []
            for x, y in zip(x_streak, y_streak):
                fit_val.append(streak[int(y), int(x)])
  
            
            mid_arg = np.where(np.abs(fit_val) == np.min(np.abs(fit_val)))[0]
            if len(mid_arg) == 1:
                x_mid = x_streak[mid_arg]
                y_mid = y_streak[mid_arg]
            elif len(mid_arg) > 1:
                x_mid = x_streak[mid_arg[0]]
                y_mid = y_streak[mid_arg[0]]


        if motion == "N-S":
            min_arg = np.where(fit_val == np.min(fit_val))[0]
            max_arg = np.where(fit_val == np.max(fit_val))[0]

            ### check dimension
            if float(min_arg.shape[0]) > 1:
                min_arg = min_arg[0]
            elif float(max_arg.shape[0]) >1:
                max_arg = max_arg[0]

            y_min = y_fit[min_arg]
            y_max = y_fit[max_arg]

            y_streak = np.linspace(y_min, y_max, 100)
            x_streak = function(y_streak)

            fit_val = []
            for x,y in zip(x_streak, y_streak):
                fit_val.append(streak[int(y), int(x)])

            mid_arg = np.where(np.abs(fit_val) == np.min(np.abs(fit_val)))[0]
            if len(mid_arg) == 1:
                x_mid = x_streak[mid_arg]
                y_mid = y_streak[mid_arg]
            elif len(mid_arg) > 1:
                x_mid = x_streak[mid_arg[0]]
                y_mid = y_streak[mid_arg[0]]

        #print(int(x_mid))
        x_mid = float(x_mid)
        y_mid = float(y_mid)
        x_array.append(x_mid)
        y_array.append(y_mid)

    return x_array, y_array, time_array
       


def XY2RADEC(x,y, xerr, yerr):
    """
    converts pixel coords to ra dec

    Parameters
    ----------
    x           : x pixel
    y           : y pixel
    xerr        : x error in pixel 
    yerr        : y error in pixel

    Returns
    -------
    RA          : right ascension (deg)
    DEC         : declination (deg)
    RAerr       : error in right ascension (deg)
    DECerr      : error in declination (deg)
    """
    pixcrd = np.array([[0 ,0 ],[x, y]], dtype=np.float64)
    world = wcs.wcs_pix2world(pixcrd, 0)
    ra, dec = world[1]

    ## calculate errors
    pixcrd1 = np.array([[0 ,0 ],[x+xerr, y+yerr]], dtype=np.float64)
    world1 = wcs.wcs_pix2world(pixcrd1, 0)
    ra1, dec1 = world1[1]

    pixcrd2 = np.array([[0 ,0 ],[x-xerr, y-yerr]], dtype=np.float64)
    world2 = wcs.wcs_pix2world(pixcrd2, 0)
    ra2, dec2 = world2[1]

    pixcrd3 = np.array([[0 ,0 ],[x+xerr, y-yerr]], dtype=np.float64)
    world3 = wcs.wcs_pix2world(pixcrd3, 0)
    ra3, dec3 = world3[1]

    pixcrd4 = np.array([[0 ,0 ],[x-xerr, y+yerr]], dtype=np.float64)
    world4 = wcs.wcs_pix2world(pixcrd4, 0)
    ra4, dec4 = world4[1]

    ra_array = [ra1, ra2, ra3, ra4]
    dec_array = [dec1, dec2, dec3, dec4]
    ra_err = (np.max(ra_array) - np.min(ra_array))/2.0
    dec_err = (np.max(dec_array) - np.min(dec_array))/2.0

    return ra, dec, ra_err, dec_err

def s2UTC(s):
    """
    converts ref s to utc
    Parameters
    ----------
    s   : input ref seconds
    Returns
    -------
    converted utc
    """
    
    return refUTC + timedelta(seconds=s)

def  getUTC(mid_time):
    """
    convert ref time in sec to UTC
    Paramters
    ---------
    mid_time    : mid ref time array
    Returns
    -------
    mid_UTC     : mid utc array
    """

    mid_UTC  = []

    for mt in mid_time:
        mid_UTC.append(s2UTC(mt))

    return mid_UTC


def main(args):

    ## get relevant timeSteps
    global headTimeSteps, tailTimeSteps, midTimeSteps
    headTimeSteps, tailTimeSteps, midTimeSteps = obtainTimeSteps(args.norad)

    if debug:
        print("the selected timeSteps " + str(headTimeSteps))

    ## get global variables
    global wcs, imgSize, refUTC, imgScale
    hdu = fits.open("testData/" + args.filePrefix + "-t" + str(headTimeSteps[0]).zfill(4) + ".fits")
    wcs = WCS(hdu[0].header, naxis=2)
    imgSize = hdu[0].header["NAXIS1"]
    imgScale = np.abs(hdu[0].header["CDELT2"])
    refUTC = datetime.strptime(hdu[0].header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')

    ## obtain TLE for the date range
    line1, line2, line3 = obtainTLE(args.norad)

    ## get satellite object
    sat, mwa, ts = getSatMWA(line1, line2, line3)

    if debug:
        print("line1 {}\nline2 {}\nline3 {}".format(line1, line2, line3))

    ## get head data
    head_x, head_y, head_x_err, head_y_err = getHeadPeaks(headTimeSteps, args.norad, mwa, sat, ts, args)
    tail_x, tail_y, tail_x_err, tail_y_err = getTailPeaks(tailTimeSteps, args.norad, mwa, sat, ts, args)

    
    
    ## determine if the motion if E-W or N-S 
    deltaX = abs(np.max(head_x) - np.min(head_x))
    deltaY = abs(np.max(head_y) - np.min(head_y))
    
    global motion
    if deltaX >= deltaY:
        motion = "E-W"
    elif deltaY > deltaX:
        motion = "N-S"

    if debug:
        print("motion classified as " + motion)

    ## monte carlo thro to propogate the errors in head and
    ## tail positions to mid positions

    if debug:
        print("monte carlo sampling of errors")

    mid_x, mid_y, mid_time = [], [], []

    for iter in tqdm(range(100)):
        head_x_tmp, head_y_tmp, tail_x_tmp, tail_y_tmp = [], [], [], []
        for h_x, h_y, x_err, y_err in zip(head_x, head_y, head_x_err, head_y_err):
            sx = np.random.normal(h_x, x_err, 200)
            x = float(np.random.choice(sx, 1)[0])
            sy = np.random.normal(h_y, y_err, 200)
            y = float(np.random.choice(sy, 1)[0])
            head_x_tmp.append(x)
            head_y_tmp.append(y)

        for t_x, t_y, x_err, y_err in zip(tail_x, tail_y, tail_x_err, tail_y_err):
            sx = np.random.normal(t_x, x_err, 200)
            x = float(np.random.choice(sx, 1)[0])
            sy = np.random.normal(t_y, y_err, 200)
            y = float(np.random.choice(sy, 1)[0])
            tail_x_tmp.append(x)
            tail_y_tmp.append(y)

        x_fit, y_fit , function = fitLine(head_x_tmp, head_y_tmp, tail_x_tmp, tail_y_tmp)
        x_array , y_array, t_array = getMidPoints(midTimeSteps, mwa, sat, ts, x_fit, y_fit, function, args)
        mid_x.append(x_array)
        mid_y.append(y_array)
        mid_time = t_array

    mid_x = np.array(mid_x)
    mid_y = np.array(mid_y)
    mean_x = np.average(mid_x, axis=0)
    mean_y = np.average(mid_y, axis=0)
    std_x = np.std(mid_x, axis=0)
    std_y = np.std(mid_y, axis=0)

    ## save data to file
    ra_array, dec_array, ra_err_array, dec_err_array = [], [], [], []
    x_array, y_array, x_err_array, y_err_array = [], [], [], []

    ## get position info
    for x, xerr, y, yerr in zip(mean_x, std_x, mean_y, std_y):
        ra, dec, ra_err, dec_err= XY2RADEC(x,y, xerr, yerr)
        ra_array.append(ra)
        dec_array.append(dec)
        ra_err_array.append(ra_err)
        dec_err_array.append(dec_err)
        x_array.append(x)
        y_array.append(y)
        x_err_array.append(xerr)
        y_err_array.append(yerr)
        
    ## get time info
    mid_UTC = getUTC(t_array)

    with open(str(args.norad) + "_extracted_angular_measurements_" + str(args.obs) + ".csv", "w") as vsc:
        
        if debug:
            print("writting data to file")
        
        thewriter = csv.writer(vsc)
        thewriter.writerow(['x', 'x_err', 'y', 'y_err', 'RA', 'RA_err', 'DEC', 'DEC_err', 'UTC', 'timeStep'])

        for x, x_err, y, y_err, ra, ra_err, dec, dec_err, utc, t in zip(x_array, x_err_array, y_array, y_err_array, ra_array, ra_err_array, dec_array, dec_err_array, mid_UTC, t_array):
            output = [round(x,4), round(x_err,4), round(y,4), round(y_err,4), round(ra, 4), round(ra_err,4), round(dec,4), round(dec_err,4), utc, t]
            thewriter.writerow(output)




    
        

if __name__ == "__main__":
    parser = ArgumentParser("extractAngularMeasurements", description="extracts angular position measurements of the satellite pass")
    parser.add_argument("--obs", required=True, type=int, help="the obsid")
    parser.add_argument("--norad", required=True, type=int, help="the norad id of the satellite")
    parser.add_argument("--beamFile", required=True, help="the beam file")
    parser.add_argument("--user", required=True, help="the user name of space-track.org")
    parser.add_argument("--passwd", required=True, help="the password of space-track.org")
    parser.add_argument("--debug", default=False, type=bool, help="run script in debug mode (default=False)")
    parser.add_argument("--filePrefix", default="6Sigma1FloodfillSigmaRFIBinaryMapPeakFlux", type=str, help="the prefix of the file name")
    args = parser.parse_args()

    global debug, query, beam
    debug = args.debug
    query = st.SpaceTrackClient(args.user, args.passwd)
    hdu = fits.open(args.beamFile)
    beam = hdu[0].data

    if debug:
        print("running extractAngularMeasurements.py in debug mode")

    main(args)
    