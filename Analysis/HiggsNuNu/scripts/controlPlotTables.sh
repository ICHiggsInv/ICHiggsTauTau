#!/bin/sh
PRODUCTION=Apr04
PARAMS=./filelists/$PRODUCTION/Params${PRODUCTION}.dat
#PARAMS=./filelists/$PRODUCTION/Params${PRODUCTION}_noSignal.dat
DOLUMIXSWEIGHT=false #should be true if controlplots is to take care of lumi*xs/events weight

for CHANNEL in nunulowmet #nunuiglep mumu #nunu taunu enu munu #mumu nunuiglep
  do
  for MET in 130 #0 #70
    do
    for SYST in central #JESUP JESDOWN JERBETTER JERWORSE ELEEFFUP ELEEFFDOWN MUEFFUP MUEFFDOWN #PUUP PUDOWN
      do
	FOLDER=output_rereco/$CHANNEL/MET$MET/ # ./oldanalysisruns/090813_nopuidorjer/output/$CHANNEL/MET$MET/
	PLOTDIR=TABLES_rereco/$CHANNEL/MET$MET/
	PLOTDIRQCD=TABLES_rereco/$CHANNEL/MET$MET/QCD/

	if [ "$SYST" != "central" ] #if not doing central
            then
	    if [[ "$SYST" != PU* ]] #For PU syst info is in central root file if not doing PU get output from syst subdirectory
		then
		FOLDER=$FOLDER"/"$SYST"/"
# 	    else #If doing PU syst pass correct options
#  		DOPUSYST="true"
	
#  		if [ "$SYST" = "PUUP" ]
#  		    then
#  		    PUUPORDOWN="true" #note this is false by default so PUDOWN will have correct option already
#  		fi
	    fi
	    PLOTDIR=$PLOTDIR"/"$SYST"/"
	    PLOTDIRQCD=$PLOTDIR"/"$SYST"/QCD/"
        fi

	
	mkdir -p $PLOTDIR
	mkdir -p $PLOTDIRQCD
	mkdir -p $PLOTDIR/wjetsComp/
	BLIND=0
	if [ "$CHANNEL" != "nunu" ] || (( "$MET" != "130" ))
	    then
	    let BLIND=0
	fi
	
 #     echo $CHANNEL " " $MET " " $DOQCD " : BLIND variable is "$BLIND
	
	
###### jpt_1
#./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
#    --plot_name="jpt_1"  --x_axis_label="Leading Jet p_{T} [GeV]" \
#    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=1000 \
#    --y_axis_min=0.01 --extra_pad=1000 \
#    --rebin=20 \
#    --norm_bins=false --verbose=false \
#    --log_y=true \
#    --paramfile=$PARAMS

###### n_jets
./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --folder=$FOLDER --plot_dir=$PLOTDIR \
    --plot_name="n_jets"  --x_axis_label="Number of jets" \
    --blind=$BLIND \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=20 \
    --y_axis_min=0.01 --extra_pad=1000 \
    --rebin=1 \
    --norm_bins=false --verbose=false \
    --plot_qcd=false \
    --log_y=true \
    --dolumixsweight=$DOLUMIXSWEIGHT \
    --paramfile=$PARAMS \
    --use_rereco_data=true

 
./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --folder=$FOLDER --plot_dir=$PLOTDIR \
    --plot_name="n_jets_puUp"  --x_axis_label="Number of jets" \
    --blind=$BLIND \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=20 \
    --y_axis_min=0.01 --extra_pad=1000 \
    --rebin=1 \
    --norm_bins=false --verbose=false \
    --plot_qcd=false \
    --log_y=true \
    --dolumixsweight=$DOLUMIXSWEIGHT \
    --paramfile=$PARAMS\
    --use_rereco_data=true

./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --folder=$FOLDER --plot_dir=$PLOTDIR \
    --plot_name="n_jets_puDown"  --x_axis_label="Number of jets" \
    --blind=$BLIND \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=20 \
    --y_axis_min=0.01 --extra_pad=1000 \
    --rebin=1 \
    --norm_bins=false --verbose=false \
    --plot_qcd=false \
    --log_y=true \
    --dolumixsweight=$DOLUMIXSWEIGHT \
    --paramfile=$PARAMS\
    --use_rereco_data=true

###### dphijj
./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --folder=$FOLDER --plot_dir=$PLOTDIR  \
    --plot_name="dphijj"  --x_axis_label="#Delta#phi_{jj}" \
    --blind=$BLIND  --x_blind_min=0 --x_blind_max=1. \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=3.15 \
    --y_axis_min=0.01 --extra_pad=100000 \
    --rebin=4 \
    --plot_wjets_comp=false \
    --norm_bins=false \
    --plot_qcd=false \
    --log_y=true \
    --dolumixsweight=$DOLUMIXSWEIGHT \
    --paramfile=$PARAMS\
    --use_rereco_data=true

./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --folder=$FOLDER --plot_dir=$PLOTDIRQCD  \
    --plot_name="dphijj"  --x_axis_label="#Delta#phi_{jj}" \
    --blind=$BLIND  --x_blind_min=0 --x_blind_max=1. \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=3.15 \
    --y_axis_min=0.01 --extra_pad=100000 \
    --rebin=4 \
    --plot_wjets_comp=false \
    --norm_bins=false \
    --plot_qcd=true \
    --log_y=true \
    --dolumixsweight=$DOLUMIXSWEIGHT \
    --paramfile=$PARAMS\
    --use_rereco_data=true

./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --folder=$FOLDER --plot_dir=$PLOTDIR  \
    --plot_name="dphijj"  --x_axis_label="#Delta#phi_{jj}" \
    --blind=$BLIND  --x_blind_min=0 --x_blind_max=1. \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=3.15 \
    --y_axis_min=0.01 --extra_pad=2 \
    --rebin=4 \
    --plot_wjets_comp=false \
    --norm_bins=false \
    --plot_qcd=false \
    --log_y=false \
    --dolumixsweight=$DOLUMIXSWEIGHT \
    --paramfile=$PARAMS\
    --use_rereco_data=true

./bin/ControlPlots --cfg=scripts/controlPlot.cfg  \
    --folder=$FOLDER --plot_dir=$PLOTDIR"/wjetsComp/"  \
    --plot_name="dphijj"  --x_axis_label="#Delta#phi_{jj}" \
    --blind=$BLIND  --x_blind_min=0 --x_blind_max=1. \
    --custom_x_axis_range=true --x_axis_min=0 --x_axis_max=3.15 \
    --y_axis_min=0.01 --extra_pad=100000 \
    --rebin=4 \
    --plot_wjets_comp=true \
    --norm_bins=false \
    --plot_qcd=false \
    --log_y=true \
    --dolumixsweight=$DOLUMIXSWEIGHT \
    --paramfile=$PARAMS\
    --use_rereco_data=true

    done
  done
done