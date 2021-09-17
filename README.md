# MWA-OrbitDetermination

Below are the scipts used for performing orbit determination in the paper 
"Orbit Determination and Catalogue Maintenance for LEO Objects using Murchison Widefield Array"

The demo pipeline does orbit determination from real data in 4 steps

1) identify autonomously all the objects detected in the observation (satSearch.py)
2) extract angular position measurements of the objects detected (extractAngularMeasurements.py)
3) create a config file for orbit determination using past published TLE data with boundary conditions (createConfig.py)
4) perfrom orbit determination to estimate the orbital elements along withs associated errors (orbitFit.py)

## Step 1 Autonomous cross-matching of all detected events in the observation

The positive and negative events detected by [RFISeeker](https://github.com/StevePrabu/RFISeeker) is the input data required by satSearch.py. 
More information on how RFISeeker works can be found [Prabu et al 2020b](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/lowfrequency-blind-survey-of-the-low-earth-orbit-environment-using-noncoherent-passive-radar-with-the-murchison-widefield-array/BF1BFD69F15D72D65514E95868F21DBA)

The output of satSearch.py 

![output](https://github.com/PhD-Misc/MWASSA/blob/master/image1.gif)

In above animation, the blue contours are all the 6 sigma events detected by RFISeeker. SatSearch.py does an API query for the epoch to identify all 
the objects within the FOV of the observation, followed by a cone search for the event (with 7km positional error). In the above animation, the search cone is green if the event is detected by the script, else red. The detected events (along with other required information) are saved to disk for further analysis.

## Step 2 Extraction of angular position measurements of the pass

extractAngularMeasurements.py extracts the angular position measurements by trying to extract the mid-point of the streak (more info in paper)

This steps writes to disk the following information to disk (x(pixel), x_err, y(pixel), y_err, RA, RA_err, DEC, DEC_err, UTC, time-step). 

## Step 3 Create config file with intial guess and boundary conditions for orbit determinaiton

This step, looks at historically published TLE for the object, and uses the most recent past TLE as the intial guess. It also generates boundary condition
for each orbital element by inspecting the historic evolution of the TLE of the object. 

## Step 4 Orbit Determination

Inputs the config file created in step 3 and the angular position measurements obtained in step 2 to perform orbit fit to the satellite pass. The script also writes to disk the estimated orbital elements (along with uncertainities) for the object at the epoch of measurement.





 




