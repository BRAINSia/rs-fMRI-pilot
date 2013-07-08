#!/bin/bash

NO_DICOM_ERROR=4
dicomDir=$1
MODALITY_=('fmri' 'mri' 'dti' 'none')

# Determine the DICOM extension
dicomExt="*.dcm"
echo "Trying dcm..."
files=($(ls -a1 $dicomDir/${dicomExt}))
if [[ ${#files[@]} == 0 ]]; then
    dicomExt=".*.IMA"
    files=($(ls -a1 $dicomDir/${dicomExt}))
    echo "IMA files"
    if [[ ${#files[@]} == 0 ]]; then
        return ${NO_DICOM_ERROR}
    fi
fi
# Get the name of the last dicom file
eximage=`ls -1 $dicomDir/${dicomExt} | head -1`
#   description
descr=`dicom_hdr $eximage | grep "Series Description" | awk -F "/" '{print $5}'`
#   repetition time
tr=`dicom_hdr $eximage | grep Repetition | awk -F "/" '{print $5}'`
# Number of slices in each timepoint
# Expanded grep search avoids error when sGroupArray.anMember line is incomplete and should be ignored
nSlices=`strings $eximage | grep "sGroupArray.anMember\[[0-9]\+\][\t ]\+= [-0-9]\+" | wc | awk '{print $1}'`
# Sort the dicom series
images=`ls $dicomDir/${dicomExt} | sort -t. -k 5 -n`
# Count the number of dicom files
nimages=`echo $images | wc | awk '{print $2}'`

#   echo "$images"
#   echo "==========================================="
#   echo "$tr"
#   echo "$descr"
#   echo "$nimages"
#   echo "$nSlices"
#   echo "==========================================="
#   echo "==========================================="

type=""
case $descr in
    *localizer*)
        modality=${MODALITY_[3]}
        label=${snum}_Loc
        ;;
    *Field*)
        modality=${MODALITY_[3]}
        label=${snum}_FieldMap
        ;;
    *MPRAGE*)
        modality=${MODALITY_[1]}
        type="spgr"
        label=${snum}_T1
        ;;
    *T2*)
        modality=${MODALITY_[1]}
        type="fse"
        label=${snum}_T2
        ;;
    *REST*)
        modality=${MODALITY_[0]}
        label=${snum}_REST
        ;;
    *GO*)
        modality=${MODALITY_[0]}
        label=${snum}_GoNoGo
        ;;
    *SIMON*)
        modality=${MODALITY_[0]}
        label=${snum}_Simon
        ;;
    *COLOR*)
        modality=${MODALITY_[0]}
        label=${snum}_ColorCD
        ;;
    *SHAPE*)
        modality=${MODALITY_[0]}
        label=${snum}_ShapeCD
        ;;
    *SWITCH*)
        modality=${MODALITY_[0]}
        label=${snum}_SwitchCD
        ;;
    *MIXED*)
        modality=${MODALITY_[0]}
        label=${snum}_MixedCD
        ;;
    *DCCS*)
        modality=${MODALITY_[0]}
        label=${snum}_DCCS
        ;;
    CONNECTIVITY)
        modality=${MODALITY_[0]}
        label=${snum}_CONNECTIVITY
        ;;
    Functional*)
        modality=${MODALITY_[0]}
        label=${snum}_CONNECTIVITY
        ;;
    ep2d_pace_pmu_float_conn_132_*)
        modality=${MODALITY_[0]}
        label=${snum}_CONNECTIVITY
        ;;
    *)
        echo "$descr"
        modality=${MODALITY_[3]}
        label=${snum}_Unknown
        ;;
esac

if [ "$modality" == "fmri" ]; then
    # grep the value of sSliceArray.ucMode
    order=`strings $eximage | grep sSliceArray.ucMode | awk '{print $3}'`
    isOdd=`expr $nSlices % 2`
    if [ "$order" == "0x4" ]; then
        # interleaved
        if [ "$isOdd" == "1" ]; then
            sliceOrder="alt+z"
        else
            sliceOrder="alt+z2"
        fi
    elif [ "$order" == "0x2" ]; then
        # descending
        sliceOrder="seq-z"
    elif [ "$order" == "0x1" ]; then
        # ascending
        sliceOrder="seq+z"
    else
        # HACK/GUESS
        # Don't know what the code would be, so using this as a catch-all
        echo "Order: "$order
        if [ "$isOdd" == "1" ]; then
        #     sliceOrder="alt-z"
            echo "Dicom is interleaved"
        else
        #     sliceOrder="alt-z2"
            echo "Dicom is not interleaved"
        fi
        # END HACK
        exit 2
    fi
    echo "$modality $nSlices $nimages $tr $sliceOrder"
elif  [ "$modality" == "mri" ]; then
    echo "$modality $type"
    exit 1
else
    echo "$modality"
    exit 3
fi

#   case $modality in
#     "mri")
#           to3d -$type -session $dataDir/Analysis/${label} -prefix ${scanId}_${label} -view orig $images
#           cd $dataDir/Analysis/$label
#           3dAFNItoNIFTI -prefix ${scanId}_${label}.nii.gz ${scanId}_${label}+orig
#           ;;
#     "fmri")
#           to3d -session $dataDir/Analysis/${label} -prefix ${scanId}_${label} -time:zt $nSlices $nimages $tr ${sliceOrder} ${images}
#           ;;
#     "dti")
#           echo "DTI: nothing to do..."
#           ;;
#     "none"|*)
#           echo "Image is neither MRI or fMRI.  Exiting..."
#           ;;
#   esac
exit 0

####JV NOTES.

#Added the following for fMRI_HD_120 and PHD_120
#    CONNECTIVITY)
#        modality=${MODALITY_[0]}
#        label=${snum}_CONNECTIVITY
#        ;;
#    Functional*)
#        modality=${MODALITY_[0]}
#        label=${snum}_CONNECTIVITY
#        ;;
