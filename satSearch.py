#!/usr/bin/env python
from __future__ import print_function, division
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from argparse import ArgumentParser
import spacetracktool as st
from spacetracktool import operations as ops
from datetime import datetime, timedelta
from os import path
from skyfield.api import Topos, load
from skyfield.api import EarthSatellite
import json
import copy
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def getTLE(t):
    """
    gets the catalog of unique objects for the requested date range

    Parameters
    ----------
    t           : timestep for reference date (can be any from observation)

    Returns
    -------
    catalog     : the catalog of obtained objects
    entries     : number of objects in catalog
    """
    ## obtain tle for objects updated within one week
    time1 = startUTC + timedelta(hours=-24*1)
    time2 = startUTC + timedelta(hours=24*1) 

    day1, month1, year1 = str(time1.day).zfill(2), str(time1.month).zfill(2), str(time1.year)
    day2, month2, year2 = str(time2.day).zfill(2), str(time2.month).zfill(2), str(time2.year)
    custom_name = year1 + "-" + month1 + "-" + day1 + "__" +  year2 + "-" + month2 + "-" + day2
    entries = 0

    if path.exists("TLE_catalog" + custom_name + ".txt") == False:
                
        if debug:
            print("requesting file from server")

        date_range = ops.make_range_string(year1 + "-" + month1 + "-" + day1, year2 + "-" + month2 + "-" + day2)
        result = query.tle_query(epoch=date_range)
        #result = query.tle_publish_query(publish_epoch=date_range)
        
        ## write catalog to file
        with open("TLE_catalog" + custom_name + ".txt", "w") as outfile:
            json.dump(result.json(), outfile)
        
        entries = len(result.json()[:])
        catalog = result.json()

        ## iterate thro and find the unique objects found
        
        norad_array = []
        for i in range(entries):
            line2 = catalog[i]["TLE_LINE1"]
            norad = line2[2:7]
            apogee = float(catalog[i]['APOGEE'])
            perigee = float(catalog[i]["PERIGEE"])
            
            ## filter out nonLEO objects and duplicates
            if int(norad) in norad_array:
                continue

            elif perigee > 2000:
                continue

            else:
                norad_array.append(int(norad))
            

    else:
        if debug:
            print("tle file found on disk. Not downloading.")
            print("using file "+ "TLE_catalog" + custom_name + ".txt")

        with open("TLE_catalog" + custom_name + ".txt") as json_file:
            result = json.load(json_file)

        entries = len(result[:])
        catalog = result

        ## iterate thro and find the unique objects found
        
        norad_array = []
        for i in range(entries):
            line2 = catalog[i]["TLE_LINE1"]
            norad = line2[2:7]
            apogee = float(catalog[i]['APOGEE'])
            perigee = float(catalog[i]["PERIGEE"])

            ## filter out non leo objects and duplicates
            if int(norad) in norad_array:
                continue

            elif perigee > 2000:
                continue

            else:
                norad_array.append(int(norad))

    ## sort out catalog to have only the closest epoch TLE
    sorted_catalog = []

    if debug:
        print("sorting catalog to have only the closest epoch TLE")

    for norad in tqdm(norad_array):
        ref_time = []
        index_array = []
        for i in range(entries):
            line2 = catalog[i]["TLE_LINE1"]
            local_norad = int(line2[2:7])
            if local_norad != norad:
                continue
            epoch = datetime.strptime(catalog[i]['EPOCH'], '%Y-%m-%d %H:%M:%S')
            diff = abs((startUTC - epoch).total_seconds())
            ref_time.append(diff)
            index_array.append(i)
        try:
            min_index = int(np.where(ref_time == np.min(ref_time))[0])
        except:
            min_index = int(np.where(ref_time == np.min(ref_time))[0][0])
    
        sorted_catalog.append(catalog[index_array[min_index]])
    
    return sorted_catalog, norad_array 


def getSatXY(line1, line2, line3, UT, mwa, ts):
    """
    determines if the satellite is within fov and if so returns 
    the x y pixel location of the satellite

    Paramters
    ---------
    line1   : line1 of tle
    line2   : line2 of tle
    line3   : line3 of tle
    UT      : UTC time
    mwa     : skyfield observer object
    ts      : skyfield time object

    Returns
    -------
    px                  : x location of satellite
    py                  : y location of satellite
    number_of_pixels    : search cone radius
    visible             : visible? (bool)
    """

    visible = False
    number_of_pixels = 0
    sat = EarthSatellite(line2, line3, line1, ts)

    # propogate satellite to utc
    time = ts.utc(UT.year, UT.month,  UT.day, UT.hour, UT.minute, UT.second)
    sat.at(time)
    difference = sat - mwa
    topocentric = difference.at(time)

    # determine angular and pixel space location of satellite
    ra, dec, distance = topocentric.radec()
    ra, dec = np.degrees(ra.radians), np.degrees(dec.radians)
    px, py = wcs.all_world2pix([[ra, dec]], 1)[0]

    ## check if satellite within image

    if 0 < px < imgSize and 0 < py < imgSize:

        if np.sqrt((px - (imgSize/2))**2 + (py -(imgSize/2))**2) < imgSize/3:
            LOS_range = distance.m

            if LOS_range < 2000000:
                visible = True

            radius = np.degrees(70000/LOS_range) # 25 km offset cone in orbit (can be made smaller) s=rxtheta
            number_of_pixels = radius/pixel_scale 

    return px, py, number_of_pixels, visible

def checkDetected(px,py,radius,data):
    """
    creates a mask of given radius and checks if we a 
    event match
    
    Paramters
    ---------
    px      : the x pixel location
    py      : the y pixel location
    radius  : the radius of search cone in pixels
    data    : the data from rfiseeker

    Returns
    -------
    detected    : detected?? (bool)
    """
    x = np.linspace(0, (imgSize-1), imgSize)
    y = np.linspace(0, (imgSize-1), imgSize)
    x, y = np.meshgrid(x, y)

    array_radius = np.sqrt((x - px)**2 + (y - py)**2)
    array_radius[array_radius > radius] = -10
    array_radius[array_radius != -10] = 1
    array_radius[array_radius != 1] = 0

    masked_data = array_radius * data
    
    if np.sum(masked_data) > 0:
        detected = True

    else:
        detected = False

    return detected


def filterDetections(detectionSummary):
    """
    filters events that do not fit the below critera
    -appear more than three times
    -have atleast 2 consec envents
    -no detect timespacing less than 3

    Paramters
    ---------
    detectionSummary   : the dictionary of the detected events

    Returns
    -------
    filteredSummary    : the output detection summary
    """



    ## remove objects will less than 4 occurances and has no consec detections
    filteredSummary = dict()
    for obj in detectionSummary:

        noradid = detectionSummary[obj]["norad"]
        total = detectionSummary[obj]["total"]
        timeSteps = detectionSummary[obj]["timeSteps"]

        if float(total) > 4:

            timeStep_array = timeSteps.split("-")

            for i in timeStep_array:

                if str(int(i) + 1) in timeSteps:
                    
                    filteredSummary[obj] = dict(timeSteps=timeSteps, total=total, norad=noradid)

    

    ## removes detection with more that 3 timestep diff
    all_timeSteps = []
    for obj in filteredSummary:

        total = filteredSummary[obj]["total"]
        timeSteps = filteredSummary[obj]["timeSteps"]

        output_timeSteps = []
        output_total = 0
        timeStep_array = timeSteps.split("-")

        for i in timeStep_array:

            if str(int(i)+1) in timeStep_array or str(int(i)-1) in timeStep_array or \
                str(int(i)+2) in timeStep_array or str(int(i)-2) in timeStep_array or \
                    str(int(i)+3) in timeStep_array or str(int(i)-3) in timeStep_array:
                
                output_timeSteps.append(str(i))
                all_timeSteps.append(int(i))
                output_total += 1

        filteredSummary[obj]["total"] = str(output_total)

        filteredSummary[obj]["timeSteps"] = "-".join(output_timeSteps)

    return filteredSummary



def main(args):

    timeSteps = np.arange(args.t1, args.t2)

    if debug:
        print("the selected timeSteps are " + str(timeSteps))

    ## get tle catalog
    catalog, unique_norads = getTLE(timeSteps[0])

    if debug:
        print("obtained TLE for {0} objects".format(len(unique_norads)))

    ## begin search
    ts = load.timescale()
    mwa = Topos("26.703319405555554 S", "116.91558083333334 E", elevation_m= 377.827)
    globalData = np.zeros((imgSize, imgSize))

    ## the below is used for plotting
    sat_x = []
    sat_y = []

    detectionSummary = dict()    

    if debug:
        print("starting search...")

    for t in tqdm(timeSteps):

        hdu = fits.open("testData/6Sigma1FloodfillSigmaRFIBinaryMap-t" + str(t).zfill(4) + ".fits")
        hduSeed = fits.open("testData/6Sigma1FloodfillSigmaRFIBinaryMapSeed-t" + str(t).zfill(4) + ".fits")
        NeghduSeed = fits.open("testData/Neg6Sigma1FloodfillSigmaRFIBinaryMapSeed-t" + str(t).zfill(4) + ".fits")
        dataSeed = hduSeed[0].data
        NegdataSeed = NeghduSeed[0].data
        data = hdu[0].data
        time = datetime.strptime(hdu[0].header['DATE-OBS'][:-2], '%Y-%m-%dT%H:%M:%S')
        globalData += data

        ## plot
        ax = plt.subplot(1,1,1, projection=wcs)
        tempData = globalData
        tempData = np.ma.masked_where(tempData == 0 , tempData)
        cmap= copy.copy(plt.cm.Purples)
        cmap.set_bad(color="white")
        plt.imshow(tempData, origin="lower", vmax=1, vmin=-1, cmap= cmap)
        plt.title("UTC " + str(time))
        plt.xlabel("RA (Degrees)")
        plt.ylabel("DEC (Degrees)")
        plt.grid()

        
        for satNo in range(len(unique_norads)):
            line1 = "sat"
            line2 = catalog[satNo]["TLE_LINE1"]
            line3 = catalog[satNo]["TLE_LINE2"]
            norad = line2[2:7]

            ## check if within FOV
            x, y, radius, visible = getSatXY(line1, line2, line3, time, mwa, ts)

            if visible:
                
                ## check if head detected
                headDetected = checkDetected(x, y , radius, dataSeed)
                tailDetected = checkDetected(x, y , radius, NegdataSeed)
                
                sat_x.append(x)
                sat_y.append(y)

                ## proceed only if head and tail is detected
                if headDetected and tailDetected:
                    
                    circle = plt.Circle((x, y), radius, fill=False, edgecolor="lime")

                    if norad in detectionSummary:

                        ## obtain  prev values
                        prev_timeSteps = detectionSummary[str(norad)]["timeSteps"]
                        prev_total = detectionSummary[str(norad)]["total"]
                        current_timeSteps = prev_timeSteps + "-" + str(t)
                        current_total = str(int(prev_total) + 1)

                        ## update values
                        detectionSummary[str(norad)]["timeSteps"] = current_timeSteps
                        detectionSummary[str(norad)]["total"] = current_total
                          

                    else:

                        entry = dict(timeSteps=str(t), total="1", norad=norad)
                        detectionSummary[str(norad)] = entry

                else:
                    circle = plt.Circle((x, y), radius, fill=False, edgecolor="red")

                ax.add_artist(circle)

        plt.scatter(sat_x, sat_y, marker="x", color="magenta",s=1)
        plt.savefig("img" + str(t).zfill(2) + ".png", dpi=300)
        plt.clf()

    if debug:
        print("all detected events\n" + str(detectionSummary))
    
    output_summary = filterDetections(detectionSummary)

    if debug:
        print("filtered events\n" + str(output_summary))
    
    ## write to file
    json.dump(output_summary, open("filtered_summary.json", "w"))
    json.dump(detectionSummary, open("detection_summary.json", "w"))

    #### "That's all Folks!"" ###        
    if debug:
        print(" That's all Folks! ")
    

if __name__ == "__main__":
    parser = ArgumentParser("SatSearch", description="searches for satellite passes within the observation")
    parser.add_argument("--t1", required=True, type=int, help="the start timeStep")
    parser.add_argument("--t2", required=True, type=int, help="the end timeStep")
    parser.add_argument("--user", required=True, help="User name for space-track.org")
    parser.add_argument("--passwd", required=True, help="the password for space-track.org")
    parser.add_argument("--debug", default=False, type=bool, help="run the script in debug mode")
    args = parser.parse_args()

    global debug
    debug = args.debug

    global query
    query = st.SpaceTrackClient(args.user, args.passwd)

    ## get header info and make them global
    hdu = fits.open("testData/6Sigma1FloodfillSigmaRFIBinaryMap-t" + str(args.t1).zfill(4) + ".fits")
    
    global wcs, imgSize, pixel_scale, startUTC
    wcs = WCS(hdu[0].header, naxis=2)
    imgSize = hdu[0].header["NAXIS1"]
    pixel_scale = hdu[0].header["CDELT2"]
    startUTC = datetime.strptime(hdu[0].header['DATE-OBS'][:-2], '%Y-%m-%dT%H:%M:%S')

    if debug:
        print("running SatSearch.py in debug mode")

    main(args)

