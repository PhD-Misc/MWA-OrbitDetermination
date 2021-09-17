#!/bin/bash

usage()
{
echo "demoPipeline.sh [-u username] [-p password]
    -u  : space-track.org username
    -p  : space-track.org password" 1>&2;
exit 1;
}

spaceTrackUser=
spaceTrackPassword=

while getopts 'u:p:' OPTION
do
    case "$OPTION" in
        u)
            spaceTrackUser=${OPTARG}
            ;;
        p)
            spaceTrackPassword=${OPTARG}
            ;;
        ? | : | h)
            usage
            ;;
    esac
done

## if either arguments not provided, pring for help
if [[ -z ${spaceTrackUser} ]]
then
    echo "user name not provided."
    usage
fi
if [[ -z ${spaceTrackPassword} ]]
then
    echo "password not provided."
    usage
fi



#step 1
python satSearch.py --t1 1 --t2 55 --user ${spaceTrackUser} --passwd ${spaceTrackPassword} --debug True
#step 2
python extractAngularMeasurements.py --obs 1157468632 --norad 20580 --beamFile testData/beam.fits --user ${spaceTrackUser} --passwd ${spaceTrackPassword} --debug True
#step 3
python createConfig.py --noradid 20580 --wcsFile testData/6Sigma1FloodfillSigmaRFIBinaryMap-t0000.fits --debug True --user ${spaceTrackUser} --passwd ${spaceTrackPassword}
#step 4
python orbitFit.py --obs 1157468632 --norad 20580 --config auto_created_config20580.yaml --wcsFile testData/6Sigma1FloodfillSigmaRFIBinaryMap-t0000.fits --debug True

echo "done."



