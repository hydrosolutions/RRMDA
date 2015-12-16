#!/bin/bash
# Adapts the MATLABPATH to include the bin directory, runs matlab with argument/script file name.
# use : ./matlab_batcher.sh <script>.m
# @arg1 absolute path of bin directory or 'same_as_script'
# @arg2 name of matlab script

#echo "<OUTPUT>" # Content of blackbox model is written between <output> by openDA.

# What is the script location? Gives the full path of this script no matter where it is being called from.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

MATLAB_EXECUTABLE=/Applications/MATLAB_R2015a.app/bin/matlab

# Removes first and last quote of first argument.
DIRECTORY=${1}
DIRECTORY="${DIRECTORY#\"}" # Remove first fit at the front of the string.
DIRECTORY="${DIRECTORY%\"}" # Remove first fit at the back of the string.
#echo "$DIRECTORY"

# If the first argument of this script is "same_as_script" then use DIR as bin-path.
if [ "$DIRECTORY" == "same_as_script" ]; then
    DIRECTORY="$DIR"
fi

MFILE=${2}
#echo "$MFILE"
MFILE="${MFILE#\"}" # Remove first fit at the front of the string.
MFILE="${MFILE%\"}" # Remove first fit at the back of the string.
#echo "$MFILE"


# Store current MATLABPATH
old=$MATLABPATH
if [ "$DIRECTORY" == "same_as_script" ]; then
    export MATLABPATH=$DIRECTORY
else
    # Adapt paths in setup.mat.
    echo "cd app/Themi/src/"
    cd app/Themi/src/
    echo "pwd: $PWD"
    MFILE2="oda_adaptSetupMat"
    echo $MATLAB_EXECUTABLE -nojvm -nodisplay -nosplash -r "$MFILE2; exit"
    $MATLAB_EXECUTABLE -nojvm -nodisplay -nosplash -r "$MFILE2; exit"
    echo "cd DIRECTORY: $DIRECTORY"
    cd $DIRECTORY

    # Adapt MATLABPATH
    MPATH2="$DIRECTORY/app/Themi/src:$DIRECTORY/app/Themi/resources/geometry:$DIRECTORY/app/Themi/resources/samples:$DIRECTORY/src:$DIRECTORY/src/data:$DIRECTORY/src/dbase:$DIRECTORY/src/enkf:$DIRECTORY/src/enkf/prm:$DIRECTORY/src/rrm:$DIRECTORY/src/rrm/budyko:$DIRECTORY/src/setup"
    export MATLABPATH=$MPATH2
fi
echo MATLABPATH=$MATLABPATH

## matlab command by parsing the matlabFunction.
echo $MATLAB_EXECUTABLE -nojvm -nodisplay -nosplash -r "$MFILE; exit"
$MATLAB_EXECUTABLE -nojvm -nodisplay -nosplash -r "$MFILE; exit"

# Restore matlab path
export MATLABPATH=$old

#echo "</OUTPUT>"

