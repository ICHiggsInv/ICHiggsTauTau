JOBWRAPPER=./scripts/generate_job.sh
#JOBSUBMIT=true
JOBSUBMIT="./scripts/submit_ic_batch_job.sh hepmedium.q"


CONFIG=scripts/TauID_2012.cfg
echo $CONFIG

FILELIST=filelists/Apr11_Data_53X
PREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/rlane/Apr11/Data_53X/
#PATHS=(
#'SingleMu-2012A-13Jul2012-v1'
#'SingleMu-2012A-recover-06Aug2012-v1'
#'SingleMu-2012B-13Jul2012-v1'
#'SingleMu-2012C-24Aug2012-v1'
#'SingleMu-2012C-PromptReco-v2'
#'SingleMu-2012C-EcalRecover-11Dec2012-v1'
#'SingleMu-2012D-PromptReco-v1'
#)
#for i in "${PATHS[@]}"
#do
#  echo "$i"
#
# JOB="$i"_mt_2012
# $JOBWRAPPER "./bin/TauIDMeasurement --is_data=1 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=mt --output_name=Data_"$JOB".root \
# --input_prefix=$PREFIX &> jobs/$JOB.log" jobs/$JOB.sh
# $JOBSUBMIT jobs/$JOB.sh
#
#  JOB="$i"_mt_2012
#  $JOBWRAPPER "./bin/TauIDMeasurement --is_data=1 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=mt  --output_name=Data_"$JOB".root \
#   --special_mode=3 --input_prefix=$PREFIX &> jobs/Special_3_$JOB.log" jobs/Special_3_$JOB.sh
#  $JOBSUBMIT jobs/Special_3_$JOB.sh
#done


PREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/agilbert/Apr11/Data_53X/
PATHS=(
'Embedded-MTau-2012A-13Jul2012-v1'
'Embedded-MTau-2012A-recover-06Aug2012-v1'
'Embedded-MTau-2012B-13Jul2012-v1'
'Embedded-MTau-2012C-24Aug2012-v1'
'Embedded-MTau-2012C-PromptReco-v2'
'Embedded-MTau-2012D-PromptReco-v1'
)
TSCALE=(
'0'
'1'
'2'
)

for j in "${TSCALE[@]}"
do
  for i in "${PATHS[@]}"
  do
    echo "$i"

    JOB="$i"_mt_2012
    $JOBWRAPPER "./bin/TauIDMeasurement --is_data=1 --is_embedded=1 --tau_scale_mode=$j --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=mt --output_name=Embedded_$JOB.root \
    --input_prefix=$PREFIX &> jobs/$JOB-$j.log" jobs/$JOB-$j.sh
    $JOBSUBMIT jobs/$JOB-$j.sh
  done
done


FILELIST=filelists/Apr11_MC_53X
PREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/agilbert/Apr11/MC_53X/

PATHS=(
'TTJets'
'TT-v1'
'TT-v2'
'T-tW'
'Tbar-tW'
'WWJetsTo2L2Nu'
'WZJetsTo2L2Q'
'WZJetsTo3LNu'
'ZZJetsTo2L2Nu'
'ZZJetsTo2L2Q'
'ZZJetsTo4L'
#'GluGluToHToTauTau_M-125'
#'VBF_HToTauTau_M-125'
#'WH_ZH_TTH_HToTauTau_M-125'
)
for i in "${PATHS[@]}"
do
 echo "$i"

 JOB="$i"_mt_2012
 $JOBWRAPPER "./bin/TauIDMeasurement --is_data=0 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=mt --output_name=$JOB.root \
 --input_prefix=$PREFIX &> jobs/$JOB.log" jobs/$JOB.sh
 $JOBSUBMIT jobs/$JOB.sh
done

PATHS=(
'WJetsToLNu-v1'
'WJetsToLNu-v2'
'W1JetsToLNu'
'W2JetsToLNu'
'W3JetsToLNu'
'W4JetsToLNu'
)

for i in "${PATHS[@]}"
do
 echo "$i"

 JOB=WJetsToLNuSoup_"$i"_mt_2012
 $JOBWRAPPER "./bin/TauIDMeasurement --is_data=0 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=mt --output_name=$JOB.root \
 --input_prefix=$PREFIX &> jobs/$JOB.log" jobs/$JOB.sh
 $JOBSUBMIT jobs/$JOB.sh
done

PATHS=(
'DYJetsToLL'
'DY1JetsToLL'
'DY2JetsToLL'
'DY3JetsToLL'
'DY4JetsToLL'
)
for j in "${TSCALE[@]}"
do
  for i in "${PATHS[@]}"
  do
    echo "$i"

    JOB=DYJetsToTauTauSoup_"$i"_mt_2012
    $JOBWRAPPER "./bin/TauIDMeasurement --ztautau_mode=1 --is_data=0 --tau_scale_mode=$j --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=mt --output_name=$JOB.root \
    --input_prefix=$PREFIX &> jobs/$JOB-$j.log" jobs/$JOB-$j.sh
    $JOBSUBMIT jobs/$JOB-$j.sh

    JOB=DYJetsToLL-LSoup_"$i"_mt_2012
    $JOBWRAPPER "./bin/TauIDMeasurement --faked_tau_selector=1 --ztautau_mode=2 --is_data=0 --tau_scale_mode=$j --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=mt --output_name=$JOB.root \
    --input_prefix=$PREFIX &> jobs/$JOB-$j.log" jobs/$JOB-$j.sh
    $JOBSUBMIT jobs/$JOB-$j.sh

    JOB=DYJetsToLL-JSoup_"$i"_mt_2012
    $JOBWRAPPER "./bin/TauIDMeasurement --faked_tau_selector=2 --ztautau_mode=2 --is_data=0 --tau_scale_mode=$j --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=mt --output_name=$JOB.root \
    --input_prefix=$PREFIX &> jobs/$JOB-$j.log" jobs/$JOB-$j.sh
    $JOBSUBMIT jobs/$JOB-$j.sh
  done
done

