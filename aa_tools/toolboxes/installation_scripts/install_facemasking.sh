# This is a script to download FaceMasking from WUSLT NRG. The scripts checks the instalation by looking for the "spm<version>" folder
#
# source install_facemasking.sh <installation_directory>
#
# Input aguments
#   - installation_directory:   This is where FaceMasking will be installed. When configuring aap.directory_conventions.FaceMaskingdir MUST point to '<installation_directory>/mask_face/nrg-improc'

INSTDIR=$1

cd ${INSTDIR}
git clone https://github.com/mmilch01/mask_face

if [[ -d "${INSTDIR}/mask_face/nrg-improc" ]]; then
    echo "FaceMasking has been installed." >&2
else
    echo "FaceMasking has not been installed properly. Exiting..." >&2
    exit -1; 
fi