#!/usr/bin/env python
from __future__ import print_function, division
import csv
from argparse import ArgumentParser
import numpy as np
from datetime import datetime, timedelta
from skyfield.api import Topos, load, EarthSatellite
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import minimize, basinhopping, curve_fit
import matplotlib.pyplot as plt
import yaml
from tqdm import tqdm


def samplerFunc(x, y, x_err, y_err):
    """
    function used to calculate the uncertainities related to the orbit determination.
    """
    s1 = np.random.normal(x, x_err, 200) 
    s2 = np.random.normal(y, y_err, 200) 

    return np.random.choice(s1, 1)[0], np.random.choice(s2, 1)[0]



def readDetectionFile(args):
    """
    extracts the x, xerr, y, yerr, and time information from csv file
    """
    x_array, y_array = [], []
    utc_array, t_array = [], []
    with open("{}_extracted_angular_measurements_{}.csv".format(args.norad, args.obs)) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")

        for row in csv_reader:
            x,x_err,y,y_err,RA,RA_err,DEC,DEC_err,UTC,timeStep = row

            if x == "x":
                continue
            else:
                ## perturb x, y based on error (used for boostrap sampling)
                x, y = samplerFunc(float(x), float(y), float(x_err), float(y_err))

                x_array.append(x)
                y_array.append(y)
                
                utc_array.append(UTC)
                t_array.append(float(timeStep))

    return x_array, y_array, utc_array, t_array


def str2datetime(str_date):
    """
    converts string datetime to datetime object (two possible formats)
    """

    try:
        UT = datetime.strptime(str_date, '%Y-%m-%d %H:%M:%S')

    except:
        UT = datetime.strptime(str_date, '%Y-%m-%d %H:%M:%S.%f')

    return UT


def createTsTimeVector(utc_array):
    """
    create skyfield datetime object from datetime time object
    """
    
    year_array, month_array, day_array = [], [], []
    hour_array, minute_array, second_array = [], [], []

    for utc in utc_array:
        start_utc = str2datetime(utc)
        start_date, start_time = str(start_utc).split(" ")
        start_year, start_month, start_day = start_date.split("-")
        start_hour, start_minute, start_second = start_time.split(":")
        year_array.append(int(start_year))
        month_array.append(int(start_month))
        day_array.append(int(start_day))
        hour_array.append(int(start_hour))
        minute_array.append(int(start_minute))
        second_array.append(float(start_second)) ## note float and not int cos second can be decimal

    times = ts.utc(year_array, month_array, day_array, hour_array, minute_array, second_array)

    return times


def getIntialGuess(cfg):
    """
    parse the yaml file and output the relevant orbital elements
    """
    output = cfg[0]["initialGuess"]["aop"], cfg[0]["initialGuess"]["e"], cfg[0]["initialGuess"]["i"], cfg[0]["initialGuess"]["mm"], cfg[0]["initialGuess"]["ra"]
    
    return output

def dayOfYear(utc):
    """
    converts utc to day of year
    """
    time_delta = utc - datetime.strptime(str(utc.year)+"-01-01 00:00:00", '%Y-%m-%d %H:%M:%S')
    return time_delta.total_seconds()/(24.0*60*60) +1


def getSatXY(line1, line2, line3, time_vector):
    """
    propogate tle elements to epoch and determine the x,y position using wcs object
    """

    sat = EarthSatellite(line2, line3)
    sat.at(time_vector)
    difference = sat - mwa
    topocentric = difference.at(time_vector)
    ra, dec, distance = topocentric.radec()
    ra = np.degrees(ra.radians)
    dec = np.degrees(dec.radians)
    pix_coords = wcs.all_world2pix(np.array([ra, dec]).T, 1)

    px, py = pix_coords.T

    return px, py


def updateTLE_LINE2(old_line2, utc_array):
    """
    update tle line2 with the new day of year epoch
    """
    start_utc = str2datetime(utc_array[0])
    end_utc = str2datetime(utc_array[-1])

    deltaTime = (end_utc - start_utc).total_seconds()
    observation_epoch = dayOfYear(start_utc + timedelta(seconds=deltaTime/2.0))

    line2_part1 = old_line2[0:20]
    line2_part2 = "{:12.08f}".format(observation_epoch)
    line2_part3 = old_line2[32:-1]
    temp = line2_part1 + line2_part2 + line2_part3
    line2 = temp + sumCheck(temp)

    return line2, start_utc + timedelta(seconds=deltaTime/2.0)

def sumCheck(digit):
    """
    sum check for tle line 2 (used to account for missing bits during transmission)
    """

    digit = digit.replace("-", "1")    
    a = sum(int(x) for x in digit if x.isdigit())

    return str(a%10)

def fine_range_mean_anomaly_using_time_elapsed(cfg, utc):
    """
    calculate mean anomaly using the equation from paper.
    """

    prev_ma = float(cfg[4]["tle"]["line3"][43:51])
    prev_mm = float(cfg[4]["tle"]["line3"][52:63])
    prev_doy = float(cfg[4]["tle"]["line2"][20:32])

    degrees_elapsed = (dayOfYear(utc) - prev_doy)*prev_mm*360 + prev_ma

    while degrees_elapsed > 360:
        degrees_elapsed -= 360

    while degrees_elapsed  < 0:
        degrees_elapsed += 360
    
    return degrees_elapsed

def constructTLE(i, ra, e, aop, ma, mm, sid):
    """
    convert orbital elements and noradid into tle format
    """

    ## note that curve_fit/basinhopping can acidentally create negative
    ## values for orbital elements. The below section wraps it arround
    ## to a positive value

    if ra < 0:
        ra += 360
    
    if aop < 0:
        aop += 360

    if ma < 0:
        ma += 360

    line_template_part1 = "2 " + str(sid).ljust(5) + " "
    line_template_part2 = str("%08.4f" %i) + " " + str("%08.4f" %ra) +\
     " " +  str(str("%09.7f" %e).split(".")[1]) + " " + str("%08.4f" %aop) \
     + " " + str("%08.4f" %ma) + " " + str("%011.8f" %mm)  + "14849"
    line_template_part3 = sumCheck(line_template_part1 + line_template_part2)

    tle = line_template_part1 + line_template_part2 + line_template_part3

    return tle

def orbitDetermination(cfg, args):
    """
    does orbit determination
    """

    ## load angular position measurements with errors
    x_array, y_array, utc_array, t_array = readDetectionFile(args)
    time_vector = createTsTimeVector(utc_array)

    ## load intial guess
    guess_aop, guess_e, guess_i, guess_mm, guess_ra = getIntialGuess(cfg)
    
    ## determine the current day of year for observation and determin the approx. mean anomaly
    old_line2 = cfg[4]["tle"]["line2"]
    new_line2, epoch_for_ma = updateTLE_LINE2(old_line2, utc_array)
    guess_ma = fine_range_mean_anomaly_using_time_elapsed(cfg, epoch_for_ma)
    line3 = constructTLE(guess_i, guess_ra, guess_e, guess_aop, guess_ma, guess_mm, args.norad)

    ## define boundary limits
    i_min, i_max = cfg[2]['minLimits']['i'], cfg[1]['maxLimits']['i']
    ra_min, ra_max = cfg[2]['minLimits']['ra'], cfg[1]['maxLimits']['ra']
    aop_min, aop_max = cfg[2]['minLimits']['aop'], cfg[1]['maxLimits']['aop']
    ma_min, ma_max = (guess_ma - 8), (guess_ma + 8)
    e_min, e_max = cfg[2]['minLimits']['e'], cfg[1]['maxLimits']['e']
    mm_min, mm_max = cfg[2]['minLimits']['mm'], cfg[1]['maxLimits']['mm']

    constr = ({'type': 'ineq',
                'fun': lambda x: x[0] - aop_min},
                {'type': 'ineq',
                'fun': lambda x: aop_max -x[0]},
                {'type': 'ineq',
                'fun': lambda x: x[1] - ra_min},
                {'type': 'ineq',
                'fun': lambda x: ra_max - x[1]},
                {'type': 'ineq',
                'fun': lambda x: x[2] - ma_min},
                {'type': 'ineq',
                'fun': lambda x: ma_max - x[2]})

    minimizer_kwargs = {"method":"COBYLA", 
                         "constraints": constr, 
                         'options': ({'maxiter': cfg[3]['niter']['minor']})}


    def accept_test(f_new, x_new, f_old, x_old):
        output = True
        aop, ra, ma = x_new
        if aop < aop_min or aop > aop_max:
            output = False
        if ra < ra_min or ra > ra_max:
            output = False
        if ma < ma_min or ma > ma_max:
            output = False

        return output

    def func_2_minimise(i, ra, e, aop, ma, mm):

        line3 = constructTLE(i, ra, e, aop, ma, mm, args.norad)
        loc_x_array, loc_y_array = getSatXY("sat" ,new_line2, line3, time_vector)

        offset = np.mean((loc_x_array - x_array)**2 + (loc_y_array - y_array)**2)

        return offset

    def function_wraper2(x):

        aop , ra, ma = x
    
        return func_2_minimise(guess_i, ra, guess_e, aop, ma, guess_mm)

    ## basinhopping for orb elements with large dynamic change (3 of the 6 elements)
    result = basinhopping(function_wraper2, x0=(guess_aop, guess_ra, guess_ma), niter=cfg[3]["niter"]["major"],
            minimizer_kwargs=minimizer_kwargs,
            accept_test=accept_test)

    basin_aop, basin_ra, basin_ma = result.x
    line3 = constructTLE(guess_i, basin_ra, guess_e, basin_aop, basin_ma, guess_mm, args.norad)

    ## curve_fit using all 6 orbital elements (fine adjustments)
    def curve_func(time_vector_temp, i, ra, e, aop, ma, mm):
            
        line3 = constructTLE(i*guess_i, ra*basin_ra, e*guess_e, aop*basin_aop, ma*basin_ma, mm*guess_mm, args.norad)
        loc_x_array, loc_y_array = getSatXY("sat" ,new_line2, line3, time_vector)
           
        return list(loc_x_array) + list(loc_y_array)

    curve_fit_bounds = ((i_min/guess_i, ra_min/basin_ra, e_min/guess_e, aop_min/basin_aop, ma_min/basin_ma, mm_min/guess_mm),
                         (i_max/guess_i, ra_max/basin_ra, e_max/guess_e, aop_max/basin_aop, ma_max/basin_ma, mm_max/guess_mm))

    diffValue = 1 ## if curve_fit does not converge, try reducing or inceasing diffValue
    popt, pcov = curve_fit(curve_func, time_vector, list(x_array) + list(y_array), p0=(1,1,1,1,1,1), 
                        method="dogbox", diff_step=diffValue, bounds=curve_fit_bounds)

    new_i, new_ra, new_e, new_aop, new_ma, new_mm = popt[0]*guess_i, popt[1]*basin_ra, popt[2]*guess_e, popt[3]*basin_aop, popt[4]*basin_ma, popt[5]*guess_mm

    line3 = constructTLE(new_i, new_ra, new_e, new_aop, new_ma, new_mm, args.norad)
    fit_x, fit_y = getSatXY("sat" ,new_line2, line3, time_vector)
    return new_i, new_ra, new_e, new_aop, new_ma, new_mm, new_line2, line3, fit_x, fit_y




def main(args):

    ## load initial guess
    if debug:
        print("loading config file...", end="")

    with open(args.config) as yamlfile:
        cfg = yaml.load(yamlfile, Loader=yaml.FullLoader)

    if debug:
        print("done")

    global wcs, mwa, ts
    mwa = Topos("26.701276778 S", "116.670846137 E", elevation_m= 377.827)
    hdu = fits.open(args.wcsFile)
    wcs = WCS(hdu[0].header, naxis=2)
    imgSize = hdu[0].header["NAXIS1"]
    ts = load.timescale()

    ## boostrap to determine the "most-likely" value and error
    bootstrap_e, bootstrap_i, bootstrap_ra = [], [], []
    bootstrap_ma, bootstrap_mm, bootstrap_aop = [], [], []

    for i in tqdm(range(args.niter)):
        i, ra, e, aop, ma, mm, new_tle_line2 , _, _, _ = orbitDetermination(cfg, args)
        bootstrap_i.append(i)
        bootstrap_ra.append(ra)
        bootstrap_e.append(e)
        bootstrap_aop.append(aop)
        bootstrap_ma.append(ma)
        bootstrap_mm.append(mm)


    err_i, err_ra, err_e, err_aop, err_ma, err_mm = np.nanstd(bootstrap_i),np.nanstd(bootstrap_ra),np.nanstd(bootstrap_e), \
                                        np.nanstd(bootstrap_aop), np.nanstd(bootstrap_ma), np.nanstd(bootstrap_mm)
    
    new_i, new_ra, new_e, new_aop, new_ma, new_mm = np.nanmedian(bootstrap_i),np.nanmedian(bootstrap_ra),np.nanmedian(bootstrap_e), \
                                        np.nanmedian(bootstrap_aop), np.nanmedian(bootstrap_ma), np.nanmedian(bootstrap_mm)

    
    print(" i {} err {} ra {} err {} e {} err {} aop {} err {} ma {} err {} mm {} err {}".format(new_i, err_i, new_ra, err_ra, new_e, err_e, new_aop, err_aop,
                                                                                new_ma, err_ma, new_mm, err_mm))

    ## get new tle
    new_tle_line3 = constructTLE(new_i, new_ra, new_e, new_aop, new_ma, new_mm, args.norad)

    ## save data to file
    with open(str(args.obs) + "n" + str(args.norad) + ".txt", "w") as the_file:
        the_file.write("Orbit Determination Solution\n")
        the_file.write("i {0} ra {1} e {2} aop {3} ma {4} mm {5}\n".format(new_i, new_ra, new_e, new_aop, new_ma, new_mm))
        the_file.write("and the corresponding error\n")
        the_file.write("i {0} ra {1} e {2} aop {3} ma {4} mm {5}\n".format(err_i, err_ra, err_e, err_aop, err_ma, err_mm))
        the_file.write("below is the new tle\n")
        the_file.write("{}\n".format(new_tle_line2))
        the_file.write("{}".format(new_tle_line3))


        


    
    


    


if __name__ == "__main__":
    parser =  ArgumentParser("Orbit determination", description="performs orbit determination")
    parser.add_argument("--obs", required=True, type=int, help="the obs id")
    parser.add_argument("--norad", required=True, type=int, help="the norad id of the satellite")
    parser.add_argument("--config", required=True, help="the config file that contains the intial guess")
    parser.add_argument("--wcsFile", required=True, help="any file that contains the relevant wcs object")
    parser.add_argument("--niter", default=10, type=int, help="the no. of bootstrap iterations")
    parser.add_argument("--debug", default=False, type=bool, help="run the script in debug mode (default=False)")
    args = parser.parse_args()

    global debug
    debug = args.debug

    if debug:
        print("running in debug mode")

    main(args)
    