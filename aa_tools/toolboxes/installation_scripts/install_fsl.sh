# This is a script to download and install FSL. The scripts checks the instalation by looking for the "fsl" folder
#
# source install_fsl.sh <installation_directory> <version> <os> <standard_only>
#
# Input aguments
#   - installation_directory:   This is where freesurfer will be installed. When configuring aap.directory_conventions.fsldir MUST point to '<installation_directory>/fsl'
#   - version:                  Version number of the desired installation according to semantic versioning (x.y.z)
#   - os:                       ID of OS the downloaded version has been compiled for
#   - standard_only:            Only fsl/data/standard is installed. The whole FSL is not installed. 

INSTDIR=$1
VERSION=$2
FSOS=$3
STANDARDONLY=$4

cd ${INSTDIR}
echo "Installing fsl-${VERSION}-${FSOS} in ${INSTDIR}"
wget http://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-${VERSION}-${FSOS}_64.tar.gz
if [[ "x${STANDARDONLY}x" == "x1x" ]]; then
    tar xzf fsl-${VERSION}-${FSOS}_64.tar.gz fsl/data/standard
else
    tar xzf fsl-${VERSION}-${FSOS}_64.tar.gz
    ${INSTDIR}/fsl/etc/fslconf/post_install.sh -f ${INSTDIR}/fsl

    SCRIPTDIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # same folder
    sed 's/\r$//' $SCRIPTDIR/fsl_bash > $SCRIPTDIR/fsl_bash_t; mv $SCRIPTDIR/fsl_bash_t $SCRIPTDIR/fsl_bash
fi
rm fsl-${VERSION}-${FSOS}_64.tar.gz

if [[ -d "${INSTDIR}/fsl" ]]; then
    echo "FSL has been installed." >&2
else
    echo "FSL has not been installed properly. Exiting..." >&2
    exit -1; 
fi