#!/bin/sh

if [ -z $2 ]
then
    echo "Must specify <script input> <script output> <optional:GridSetup>"
    exit
fi

INPUT=$1
OUTPUT=$2
NOGRIDSETUP=0

if (( "$#" == "3" ))
    then
    echo "Option grid setup: $3"
    let NOGRIDSETUP=$3
fi

echo "cd $PWD" &> $OUTPUT
if (( "$NOGRIDSETUP" != "1" )); then
    echo "Grid setup is enabled."
    echo "source /vols/cms/grid/setup.sh" >> $OUTPUT
    echo "export SCRAM_ARCH=slc5_amd64_gcc462" >> $OUTPUT
else 
    echo "Grid setup is disabled."
    echo "export SCRAM_ARCH=slc6_amd64_gcc472" >> $OUTPUT
fi
echo "eval \`scramv1 runtime -sh\`" >> $OUTPUT
echo "source $PWD/scripts/setup_libs.sh" >> $OUTPUT
echo "eval $INPUT" >> $OUTPUT
chmod +x $OUTPUT
