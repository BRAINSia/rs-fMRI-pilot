#/bin/bash

#Dave.  It would be nice for us to be able to execute this script by providing name and session id at the terminal prompt.

#1.  Import files into Nifti format

Args="-time:zt 32 132 2000 alt+z -epan -orient RAS"
# numOfSlices = 32 (input)
# numOfVols = 132  (from DICOM)
# tr(ms) = 2000    (from DICOM: Repetition Time [0018, 0080])
to3d $Args -prefix /paulsen/Experiments/rsFMRI/${subjid}/${sessionid}/Rest.nii /paulsen/MRx/FMRI_HD_${siteid}/${subjid}/${sessionid}ANONRAW/FMRI_RestingStateConnectivity/*/*.IMA

#Dave: at this point, we need to either add full paths below or cd into the subjects correct directory.

3drefit -deoblique Rest.nii

#Note. Some of the 'args' aboe will change from study to study. Would be nice to have a way to extract the relevant information from the dicom headers.

#2.  Despike (attenuate time series outliers) and remove first four volumes (to allow scanner to get to steady state)

3dDespike -prefix Rest_ds.nii Rest.nii'[4..131]'

#Note.  The number of volumes to be removed depends on many volumes there were.  Would be nice if we can automate this so that first four volumes consistently are removed. Note that  the volumes start at 0.

#3.  Register all volumes to a single volume.  Save a motion file (mt.1D) file.

3dvolreg -prefix Rest_vr.nii -cubic -tshift 0 -zpad 3 -maxite 50  -x_thresh 0.001 -rot_thresh 0.001 \
-delta 0.1 -final Fourier -twopass -twodup -coarse 2 2 -coarserot -base 9 -dfile Rest_mt.1D Rest_ds.nii

#Note.  We are choosing volume number (after getting rid of the first 4 volumes).  This is pretty arbitrary and there are more sophisticated ways of choosing the 'best looking' volume using 3dTqual.

#4.  Scaling all voxels.  We are actually just doing a mean shift (i.e., not taking percent signal change for rsfMRI compared to task fMRI).

3dTstat -prefix Rest_mn.nii -mean Rest_vr.nii
3dcalc -prefix Rest_ms.nii -a Rest_vr.nii -b Rest_mn.nii -expr '(a - b) + 1000'

#Note. Rest_mn.nii will be used later as a registration base (i.e., we will warp this to T1).

#5.  Bandpass filter (.011 to .1 Hz) using 3dFourier.

3dFourier -prefix Rest_bp.nii -highpass .011 -lowpass .1 -retrend Rest_ms.nii

#6.  Spatial smoothing with 6mm smoothing kernel.

3dmerge -prefix Rest_smooth.nii -1blur_fwhm 6 -doall -1noneg -1clip 100 Rest_bp.nii

#7.  Coregister to T1 which is in the same space as the FreeSurfer structural images (we will be extracting signal from FreeSurfer ROIs later).

     a.  convert the FS label image to nifti format

	FreeSurfer image locations

	/paulsen/Experiments/20120722_JOY_DWI/FMRI_HD_120/$subjectid/$sessionid/10_AUTO.NN3Tv20110419/JOY_v51_2011_$subjectid_$sessionid/

	I have copied an excel spreadsheet I got from Joy that has subjectid and sessionid (Note: sessionid is referred to as experiment in the
	spreadsheet) information into our rsfMRI directory (list_JTM_edit_20120815.xlsx).  I
	will need to go over all the details of this spreadsheet with you.

	you can use the program mri_convert (which i think is a freesurfer utility to do this).

	example:
	mri_convert --in_type .mgz --out_type nii.gz subjectsfile.mgz subjectsfile.nii.gz

     b. use an afni program called 3dAllineate to do a rigid registration between the functional image and the T1.  This will happen in two steps.  first, we
     fit the mean fMRI image to the T1.  Then we apply that transform to all fMRI volumes.

# Rest_smooth.nii
     i. 3dAllineate \
         -base /paulsen/Experiments/20120722_JOY_DWI/FMRI_HD_120/$subjectid/$sessionid/10_AUTO.NN3Tv20110419/$sessionid_AVG_T1.nii.gz \
         -source Rest_mn.nii \
         -prefix rRest_mn.nii \
         -warp shr \
         -cost mi \
         -cmass \
         -interp quintic \
         -final quintic \
         -1Dmatrix_save Rest_mn-T1

	  ii. 3dAllineate \
          -master /paulsen/Experiments/20120722_JOY_DWI/FMRI_HD_120/$subjectid/$sessionid/10_AUTO.NN3Tv20110419/$sessionid_AVG_T1.nii.gz  -mast_dxyz 2.0 \
          -1Dmatrix_apply Rest_mn-T1.aff12.1D \
          -input Rest_bp.nii \
          -final quintic \
          -prefix rRest_bp.nii

###############End of first phase of preprocessing#################


##############Begin Phase 2 of preprocessing#############################
#Phase 2 involves regressing out signal of non-interest from white matter and CSF.


#1.  Generate a CSF mask using Autoworkup output.  I have been running it on the command line within BRAINS2 but it should be scripted/automated.

b2 load Image /paulsen/Experiments/20120722_JOY_DWI/FMRI_HD_120/$subjectid/$sessionid/10_AUTO.NN3Tv20110419/$sessionid_ACPC_Class.nii.gz
b2 threshold image i1 2 upperThreshold=  30
b2 load Talairach-Parameters /paulsen/Experiments/20120722_JOY_DWI/FMRI_HD_120/$subjectid/$sessionid/10_AUTO.NN3Tv20110419/???Talairach.bnd
b2 load Talairach-Box /opt/brains2/bin/talairach/vent2_box
b2 convert Talairach-Box to Mask talbox1 talpar1
b2 And masks m1 m2
b2 save mask /IPLlinux/oleary/functional/BD_COGA/142005/61793711/10_AUTO.v020/61793711_CSF.mask brains2 m3
b2 sum masks m3
b2 save image /IPLlinux/oleary/functional/BD_COGA/142005/61793711/10_AUTO.v020/61793711_CSF.nii.gz nifti i2 data-type= unsigned-8bit plane= coronal





