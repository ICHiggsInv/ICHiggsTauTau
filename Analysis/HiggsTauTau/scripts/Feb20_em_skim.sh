JOBWRAPPER=./scripts/generate_job.sh
#JOBSUBMIT=true
JOBSUBMIT="./scripts/submit_ic_batch_job.sh hepmedium.q"


CONFIG=scripts/Paper_2012.cfg
echo $CONFIG

FILELIST=filelists/Feb20_Data_53X
SKIMPATH=$SSD/Feb20/Data_53X
PREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/agilbert/Feb20/Data_53X/


PATHS=(
'MuEG-2012A-13Jul2012-v1'
'MuEG-2012A-recover-06Aug2012-v1'
'MuEG-2012B-13Jul2012-v1'
'MuEG-2012C-24Aug2012-v1'
'MuEG-2012C-PromptReco-v2'
'MuEG-2012C-EcalRecover-11Dec2012-v1'
'MuEG-2012D-PromptReco-v1'
)
for i in "${PATHS[@]}"
do
  echo "$i"

  mkdir -p $SKIMPATH/$i/em_skim/
  mkdir -p $SKIMPATH/Special_20_$i/em_skim/
  mkdir -p $SKIMPATH/Special_21_$i/em_skim/
  mkdir -p $SKIMPATH/Special_22_$i/em_skim/
  mkdir -p $SKIMPATH/Special_23_$i/em_skim/
  mkdir -p $SKIMPATH/Special_24_$i/em_skim/

  JOB="$i"_em_skim
  $JOBWRAPPER "./bin/HiggsTauTau --is_data=1 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=em --output_folder=./ --output_name=Dummy_$JOB.root \
  --do_skim=true --skim_path=$SKIMPATH/"$i"/em_skim/ --input_prefix=$PREFIX &> jobs/$JOB.log ; \
  ./scripts/hadd_and_filelist.sh $SKIMPATH "$i"/em_skim "$FILELIST"_"$JOB".dat" jobs/$JOB.sh
  $JOBSUBMIT jobs/$JOB.sh

  JOB=Special_23_"$i"_em_skim
  $JOBWRAPPER "./bin/HiggsTauTau --is_data=1 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=em --output_folder=./ --output_name=Dummy_$JOB.root \
  --do_skim=true --skim_path=$SKIMPATH/Special_23_"$i"/em_skim/ --special_mode=23 --input_prefix=$PREFIX &> jobs/$JOB.log ; \
  ./scripts/hadd_and_filelist.sh $SKIMPATH Special_23_"$i"/em_skim "$FILELIST"_"$JOB".dat" jobs/$JOB.sh
  $JOBSUBMIT jobs/$JOB.sh
  
  JOB=Special_24_"$i"_em_skim
  $JOBWRAPPER "./bin/HiggsTauTau --is_data=1 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=em --output_folder=./ --output_name=Dummy_$JOB.root \
  --do_skim=true --skim_path=$SKIMPATH/Special_24_"$i"/em_skim/ --special_mode=24 --input_prefix=$PREFIX &> jobs/$JOB.log ; \
  ./scripts/hadd_and_filelist.sh $SKIMPATH Special_24_"$i"/em_skim "$FILELIST"_"$JOB".dat" jobs/$JOB.sh
  $JOBSUBMIT jobs/$JOB.sh

  JOB=Special_20_"$i"_em_skim
  $JOBWRAPPER "./bin/HiggsTauTau --is_data=1 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=em --output_folder=./ --output_name=Dummy_$JOB.root \
  --do_skim=true --skim_path=$SKIMPATH/Special_20_"$i"/em_skim/ --special_mode=20 --input_prefix=$PREFIX &> jobs/$JOB.log ; \
  ./scripts/hadd_and_filelist.sh $SKIMPATH Special_20_"$i"/em_skim "$FILELIST"_"$JOB".dat" jobs/$JOB.sh
  $JOBSUBMIT jobs/$JOB.sh

  JOB=Special_21_"$i"_em_skim
  $JOBWRAPPER "./bin/HiggsTauTau --is_data=1 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=em --output_folder=./ --output_name=Dummy_$JOB.root \
  --do_skim=true --skim_path=$SKIMPATH/Special_21_"$i"/em_skim/ --special_mode=21 --input_prefix=$PREFIX &> jobs/$JOB.log ; \
  ./scripts/hadd_and_filelist.sh $SKIMPATH Special_21_"$i"/em_skim "$FILELIST"_"$JOB".dat" jobs/$JOB.sh
  $JOBSUBMIT jobs/$JOB.sh

  JOB=Special_22_"$i"_em_skim
  $JOBWRAPPER "./bin/HiggsTauTau --is_data=1 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=em --output_folder=./ --output_name=Dummy_$JOB.root \
  --do_skim=true --skim_path=$SKIMPATH/Special_22_"$i"/em_skim/ --special_mode=22 --input_prefix=$PREFIX &> jobs/$JOB.log ; \
  ./scripts/hadd_and_filelist.sh $SKIMPATH Special_22_"$i"/em_skim "$FILELIST"_"$JOB".dat" jobs/$JOB.sh
  $JOBSUBMIT jobs/$JOB.sh
done

PATHS=(
'Embedded-EM-2012A-13Jul2012-v1'
'Embedded-EM-2012A-recover-06Aug2012-v1'
'Embedded-EM-2012B-13Jul2012-v1'
'Embedded-EM-2012C-24Aug2012-v1'
'Embedded-EM-2012C-PromptReco-v2'
'Embedded-EM-2012D-PromptReco-v1'
)
for i in "${PATHS[@]}"
do
  echo "$i"

  mkdir -p $SKIMPATH/$i/em_skim/

  JOB="$i"_em_skim
  $JOBWRAPPER "./bin/HiggsTauTau --is_data=1 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=em --output_folder=./ --output_name=Dummy_$JOB.root \
  --do_skim=true --skim_path=$SKIMPATH/"$i"/em_skim/ --input_prefix=$PREFIX &> jobs/$JOB.log ; \
  ./scripts/hadd_and_filelist.sh $SKIMPATH "$i"/em_skim "$FILELIST"_"$JOB".dat" jobs/$JOB.sh
  $JOBSUBMIT jobs/$JOB.sh
done


FILELIST=filelists/Feb20_MC_53X
SKIMPATH=$SSD/Feb20/MC_53X
PREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/rlane/Feb20/MC_53X/
PATHS=(
'DYJetsToLL'
'DY1JetsToLL'
'DY2JetsToLL'
'DY3JetsToLL'
'DY4JetsToLL'
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
'GluGluToHToTauTau_M-125'
'VBF_HToTauTau_M-125'
'WH_ZH_TTH_HToTauTau_M-125'
'SUSYGluGluToHToTauTau_M-160'
'SUSYBBHToTauTau_M-160'
#'GluGluToHToTauTau_M-90' 
#'GluGluToHToTauTau_M-95' 
#'GluGluToHToTauTau_M-100' 
#'GluGluToHToTauTau_M-105' 
#'GluGluToHToTauTau_M-110' 
#'GluGluToHToTauTau_M-115' 
#'GluGluToHToTauTau_M-120' 
#'GluGluToHToTauTau_M-130' 
#'GluGluToHToTauTau_M-135' 
#'GluGluToHToTauTau_M-140' 
#'GluGluToHToTauTau_M-145' 
#'GluGluToHToTauTau_M-150' 
#'GluGluToHToTauTau_M-155' 
#'GluGluToHToTauTau_M-160' 
#'VBF_HToTauTau_M-90' 
#'VBF_HToTauTau_M-95' 
#'VBF_HToTauTau_M-100' 
#'VBF_HToTauTau_M-105' 
#'VBF_HToTauTau_M-110' 
#'VBF_HToTauTau_M-115' 
#'VBF_HToTauTau_M-120' 
#'VBF_HToTauTau_M-130' 
#'VBF_HToTauTau_M-135' 
#'VBF_HToTauTau_M-140' 
#'VBF_HToTauTau_M-145' 
#'VBF_HToTauTau_M-150' 
#'VBF_HToTauTau_M-155' 
#'VBF_HToTauTau_M-160' 
#'WH_ZH_TTH_HToTauTau_M-90' 
#'WH_ZH_TTH_HToTauTau_M-95' 
#'WH_ZH_TTH_HToTauTau_M-100' 
#'WH_ZH_TTH_HToTauTau_M-105' 
#'WH_ZH_TTH_HToTauTau_M-110' 
#'WH_ZH_TTH_HToTauTau_M-115' 
#'WH_ZH_TTH_HToTauTau_M-120' 
#'WH_ZH_TTH_HToTauTau_M-130' 
#'WH_ZH_TTH_HToTauTau_M-135' 
#'WH_ZH_TTH_HToTauTau_M-140' 
#'WH_ZH_TTH_HToTauTau_M-145'
#'WH_ZH_TTH_HToTauTau_M-150' 
#'WH_ZH_TTH_HToTauTau_M-155' 
#'WH_ZH_TTH_HToTauTau_M-160' 
)
for i in "${PATHS[@]}"
do
  echo "$i"

  mkdir -p $SKIMPATH/$i/em_skim/

  JOB="$i"_em_skim
  $JOBWRAPPER "./bin/HiggsTauTau --is_data=0 --cfg=$CONFIG --filelist="$FILELIST"_"$i".dat --channel=em --output_folder=./ --output_name=Dummy_$JOB.root \
  --do_skim=true --skim_path=$SKIMPATH/"$i"/em_skim/ --input_prefix=$PREFIX &> jobs/$JOB.log ; \
  ./scripts/hadd_and_filelist.sh $SKIMPATH "$i"/em_skim "$FILELIST"_"$JOB".dat" jobs/$JOB.sh
  $JOBSUBMIT jobs/$JOB.sh
done


