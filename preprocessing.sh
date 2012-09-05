#/bin/bash

#Dave.  It would be nice for us to be able to execute this script by providing name and session id at the terminal prompt.

#1.  Import files into Nifti format

Args="-time:zt 32 132 2000 alt-z -epan -orient RAS"
# numOfSlices = 32 (input [0054, 0081] <-- should be standard, but it's only for PET!)
# numOfVols = 132  (from DICOM)
# tr(ms) = 2000    (from DICOM: Repetition Time [0018, 0080])
to3d $Args -prefix /paulsen/Experiments/rsFMRI/${subjid}/${sessionid}/Rest.nii \
           /paulsen/MRx/FMRI_HD_${siteid}/${subjid}/${sessionid}ANONRAW/FMRI_RestingStateConnectivity/*/*.IMA

#Dave: at this point, we need to either add full paths below or cd into the subjects correct directory.

3drefit -deoblique Rest.nii  ### WHERE DOES THE OUTPUT FROM 3DRefit go???

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
3dcalc -prefix Rest_ms.nii -a Rest_vr.nii -b Rest_mn.nii -expr 'a - b'
# 3dcalc -prefix Rest_ms.nii -a Rest_vr.nii -expr 'a - Rest_mn.nii'

#Note. Rest_mn.nii will be used later as a registration base (i.e., we will warp this to T1).

#5.  Bandpass filter (.011 to .1 Hz) using 3dFourier.

3dFourier -prefix Rest_bp.nii -highpass .011 -lowpass .1 -retrend Rest_zpad.nii

#6.  Spatial smoothing with 6mm smoothing kernel.

3dmerge -prefix Rest_smooth.nii -1blur_fwhm 6 -doall -1noneg -1clip 100 Rest_bp.nii

#7.  Pad with zeros in IS dimension to make sure nothing gets cut off during later coregistration and warping.

3dZeropad -prefix Rest_zpad.nii -IS 44 Rest_ms.nii  ### WAS THIS SUPPOSED TO CONNECT TO 3dFourier!?!

#Note. Not clear if 3dZeropad has to happen at this point but probably better to do it after spatial smoothing.

#8.  Coregister to T1 which is in the same space as the FreeSurfer structural images (we will be extracting signal from FreeSurfer ROIs later).

Dave: can you please try to script this out?  I have provided conceptual directions.

	A.  Convert T1.mgz into nifti format.  (This step can be done earlier.)  There is a FreeSurfer command called mri_convert that will do this.  Joy told me
	about it so if you have questions, you can ask her.
	I believe the actual command is as simple as: mri_convert oldimage newimage.

	FreeSurfer image locations

	/paulsen/Experiments/20120722_JOY_DWI/FMRI_HD_120/$subjectid/$sessionid/10_AUTO.NN3Tv20110419/JOY_v51_2011_$subjectid_$sessionid/

	I have copied an excel spreadsheet I got from Joy that has subjectid and sessionid (Note: sessionid is referred to as experiment in the
	spreadsheet) information into our rsfMRI directory (list_JTM_edit_20120815.xlsx).  I
	will need to go over all the details of this spreadsheet with you.

	B.



