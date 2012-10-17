#!/bin/bash


###Note.  We assume at this point that fMRI images have been coregistered to the T1 and FS images.  Technically, these are still preprocessing steps.
###We will use the WM and CSF masks here.  Specifically, we will extract the signal from these masks and enter them as nuisance regressors in a regression analysis.  Currently, the script begins with the bandpass filtered files (output of 3dFourier).

#Extract signal from CSF mask.
if [ ! -e Median_B2_CSF.1D ] ; then 
 3dmaskave -mask B2_CSF+tlrc -median -mrange 1 1 -quiet \
  Rest_bp+tlrc >> Median_B2_CSF.1D
fi

#Extract signal from WM mask.
if [ ! -e Median_WM.1D ] ; then 
 3dmaskave -mask $restdir/WhiteMatter_ROIs+tlrc -median -mrange 1 2 -quiet \
  Rest_bp+tlrc >> Median_WM.1D
fi

#Use 3dDeconvolve run a regression where we covary out CSF, WM, and motion.  Remaining signal is the 'good' signal left over (see 'errts') that we'll use for later analyses.

if [ ! -e Rest_bp_Decon+tlrc.BRIK ] ; then 
3dDeconvolve -input Rest_bp+tlrc \
  -mask t1h+tlrc.nii \
  -GOFORIT 4 \
  -polort 1 -num_stimts 8 \
  -stim_file 1 Median_B2_CSF.1D -stim_label 1 Median_CSF \
  -stim_file 2 Median_WM.1D -stim_label 2 Median_WM \
  -stim_file 3 Rest_mt.1D[1] -stim_label 3 roll -stim_base 3 \
  -stim_file 4 Rest_mt.1D[2] -stim_label 4 pitch -stim_base 4 \
  -stim_file 5 Rest_mt.1D[3] -stim_label 5 yaw -stim_base 5 \
  -stim_file 6 Rest_mt.1D[4] -stim_label 6 dS -stim_base 6 \
  -stim_file 7 Rest_mt.1D[5] -stim_label 7 dL -stim_base 7 \
  -stim_file 8 Rest_mt.1D[6] -stim_label 8 dP -stim_base 8 \
  -full_first -float \
  -tout -rout -fout -bucket Rest_bp_Decon+tlrc -fitts full_fitts_Decon+tlrc -errts errts_Decon+tlrc
fi

#Remove linear and higher order trends not already removed by 3dFourier.  This may not be necessary but doesn't seem to hurt (or does it?!).
if [ ! -e errts_Decon_dt+tlrc.BRIK ] ; then 
3dDetrend -prefix errts_Decon_dt+tlrc -polort 3 errts_Decon+tlrc
fi

#extract the timeseries from ROIs.  Archana wants this from all FS regions as well as all voxels.  (all voxels seems somewhat strange and difficult at the moment so let's just work with the FS regions for now).  
#Here's an example of how to do it for a particular region.  Note, in this example, Left_VS2 is the name of the region.
 
if [ ! -e Median_Left_VS2_dt.1D ] ; then 
 3dmaskave -mask left_VS+tlrc.nii -median -mrange 1 1 -quiet \
  errts_Decon_dt+tlrc >> Median_Left_VS2_dt.1D
fi

done
done


#After repeating this step for all the regions, we will have a time series for each subject for each region.  Using some program (not sure if it will be AFNI or matlab or something), we need to create a NXN matrix of correlations for each subject.
