#!/bin/bash

#export PATH=${PATH}:/Users/vince/development/DICOM/dcmtk-3.6.0/dcmdata/apps:/opt/afni
#export DCMDICTPATH=/Users/vince/development/DICOM/dcmtk-3.6.0/dcmdata/data/dicom.dic
export PATH=${PATH}:/opt/afni


dataDir=$1
scanId=$2

cd $dataDir/SCANS

for series in *
do
  snum=`echo $series | awk '{printf "%03d", $1}'`
  cd $dataDir/SCANS/$series/DICOM
  eximage=`ls -1 *.dcm | head -1`
  descr=`dicom_hdr $eximage | grep "Series Description" | awk -F "/" '{print $5}'`
  tr=`dicom_hdr $eximage | grep Repetition | awk -F "/" '{print $5}'`
  nSlices=`strings $eximage | grep "sGroupArray.anMember\[" | wc | awk '{print $1}'`
  
  case $descr in
    *localizer*)
      fmri=0
      mri=0
      label=${snum}_Loc
      ;;
    *Field*)
      fmri=0
      mri=0
      dti=0
      label=${snum}_FieldMap
      ;;
    *MPRAGE*)
      fmri=0
      mri=1
      dti=0
      type=spgr
      label=${snum}_T1
      ;;  
    *T2*)
      fmri=0
      mri=1
      dti=0
      type=fse
      label=${snum}_T2
      ;; 
    *REST*)
      fmri=1
      mri=0
      dti=0
      label=${snum}_REST
      ;;
    *GO*)
      fmri=1
      mri=0
      dti=0
      label=${snum}_GoNoGo
      ;;
    *SIMON*)
      fmri=1
      mri=0
      dti=0
      label=${snum}_Simon
      ;;
    *COLOR*)
      fmri=1
      mri=0
      dti=0
      label=${snum}_ColorCD
      ;;
    *SHAPE*)
      fmri=1
      mri=0
      dti=0
      label=${snum}_ShapeCD
      ;;
    *SWITCH*)
      fmri=1
      mri=0
      dti=0
      label=${snum}_SwitchCD
      ;;
    *MIXED*)
      fmri=1
      mri=0
      dti=0
      label=${snum}_MixedCD
      ;;
    *DCCS*)
      fmri=1
      mri=0
      dti=0
      label=${snum}_DCCS
      ;;
    *)
      fmri=0
      mri=0
      dti=0
      label=${snum}_Unknown
      ;;  
  esac

  if [ "$fmri" == "1" ]; then 
    order=`strings $eximage | grep sSliceArray.ucMode | awk '{print $3}'`
    isOdd=`expr $nSlices % 2`
    
    if [ "order" == "0x4" ]; then
      if [ "$isOdd" == "1" ]; then
        sliceOrder="alt+z"
      else
        sliceOrder="alt+z2"
      fi
    fi
    if [ "order" == "0x2" ]; then
      sliceOrder="seq-z"
    else
      sliceOrder="seq+z"
    fi
    echo "$series ( $snum ) : $label $tr $nSlices $order $isOdd"
  else
    echo "$series ( $snum ) : $label $tr $nSlices"
  fi  

  images=`ls *dcm | sort -t. -k 5 -n`
  nimages=`echo $images | wc | awk '{print $2}'`
  
  
#  echo $images
  echo "==========================================="
  echo "$series $tr"
  echo "$descr"
  echo "$nimages $nSlices $snum"
  echo "==========================================="
  echo "==========================================="
  mkdir -p $dataDir/Analysis/$label
  
  if [ $mri == 1 ]; then
    to3d -$type -session $dataDir/Analysis/${label} -prefix ${scanId}_${label} -view orig $images
    cd $dataDir/Analysis/$label
    3dAFNItoNIFTI -prefix ${scanId}_${label}.nii.gz ${scanId}_${label}+orig
  fi 
  
  if [ $fmri == 1 ]; then
    to3d -session $dataDir/Analysis/${label} -prefix ${scanId}_${label} -time:zt $nSlices $nimages $tr ${sliceOrder} ${images}
  fi 
  
done
exit



