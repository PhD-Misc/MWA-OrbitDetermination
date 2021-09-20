# MWA-OrbitDetermination

Below are the scipts used for performing orbit determination in the paper 
"Orbit Determination and Catalogue Maintenance for LEO Objects using Murchison Widefield Array"

The demo pipeline does orbit determination from real data in 4 steps

1) identify autonomously all the objects detected in the observation (satSearch.py)
2) extract angular position measurements of the objects detected (extractAngularMeasurements.py)
3) create a config file for orbit determination using past published TLE data with boundary conditions (createConfig.py)
4) perfrom orbit determination to estimate the orbital elements along with its associated errors (orbitFit.py)

The steps can be executed on the demo data (observation ID 1157468632) by running the ```demoPipeline.sh``` shell script. 
In order to be able to run the shell script, the user will require an account in [space-track.org](https://www.space-track.org/auth/login).



To run the demo pipeline,
```
./demoPipeline.sh -u ${spaceTrackUser} -p ${spaceTrackPassword}
```

The python modules required by the demo pipeline can be obtained from the [docker image](https://hub.docker.com/layers/steveprabu/mypython/second/images/sha256-412b04389dabd0d668102da0f076e6085263d4e25ba5fad0b5f6abdcd4fbb5ca?context=repo) (alternatively, the required python modules with appropriate versions can be installed using the requirements.txt file provided using syntax ```pip install -r requirements.txt```)

The steps performed by the demo pipeline is explain below.

## Step 1 Autonomous cross-matching of all detected events in the observation
```
python satSearch.py --t1 1 --t2 55 --user ${spaceTrackUser} --passwd ${spaceTrackPassword} --debug True
```
The positive and negative events detected by [RFISeeker](https://github.com/StevePrabu/RFISeeker) is the input data required by satSearch.py. 
More information on how RFISeeker works can be found [Prabu et al 2020b](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/lowfrequency-blind-survey-of-the-low-earth-orbit-environment-using-noncoherent-passive-radar-with-the-murchison-widefield-array/BF1BFD69F15D72D65514E95868F21DBA)

SatSearch.py does an API query (using the user name and password provided) to obtain the TLEs of all the objects near the epoch of the observation. Using the obtained TLEs, it does a cone search for every satellite within the FOV. The cone search can identify objects with a position error upto 7km (which determines the size of the search cone which reduces with altitude). The detected events (along with other required information) is saved to disk for further analysis.

Animation of satellite cone search performed by satSearch.py is shown below

![output](https://github.com/PhD-Misc/MWASSA/blob/master/image1.gif)

## Step 2 Extraction of angular position measurements of the pass

```
python extractAngularMeasurements.py --obs 1157468632 --norad 20580 --beamFile testData/beam.fits --user ${spaceTrackUser} --passwd ${spaceTrackPassword} --debug True
```

extractAngularMeasurements.py extracts the angular position measurements by trying to extract the mid-point of the streak (more info in paper)

This step write to disk the (x(pixel), x_err, y(pixel), y_err, RA, RA_err, DEC, DEC_err, UTC, time-step) of the object for every time-step that it was detected
within the observation.


## Step 3 Create config file with intial guess and boundary conditions for orbit determinaiton

```
python createConfig.py --noradid 20580 --wcsFile testData/6Sigma1FloodfillSigmaRFIBinaryMap-t0000.fits --debug True --user ${spaceTrackUser} --passwd ${spaceTrackPassword}
```

This step, looks at historically published TLE for the object, and uses the most recent past TLE as the intial guess. It also generates boundary condition
for each orbital element by inspecting the historic evolution of the TLE of the object. 

## Step 4 Orbit Determination

```
python orbitFit.py --obs 1157468632 --norad 20580 --config auto_created_config20580.yaml --wcsFile testData/6Sigma1FloodfillSigmaRFIBinaryMap-t0000.fits --debug True
```

Inputs the config file created in step 3 and the angular position measurements obtained in step 2 to perform orbit fit to the satellite pass. The script also writes to disk the estimated orbital elements (along with uncertainities) for the object at the epoch of measurement.

Once, the orbitFit.py has sucessfully run, it will generate a file called 1157468632n20580.txt (of the format ${observationID}n${noradID}.txt) that contains the 
determined orbital elements and its errors.



 




