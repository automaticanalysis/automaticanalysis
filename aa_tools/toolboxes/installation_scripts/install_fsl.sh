# This is a script to download and install FSL. The script checks the installation by looking for the "fsl" folder
#
# source install_fsl.sh <installation_directory> <version> <os> <standard_only> <config_file>
#
# Input aguments
#   - installation_directory:   This is where freesurfer will be installed. When configuring aap.directory_conventions.fsldir MUST point to '<installation_directory>/fsl'
#   - version:                  Version number of the desired installation according to semantic versioning (x.y.z)
#   - os:                       ID of OS the downloaded version has been compiled for
#   - standard_only:            Only fsl/data/standard is installed. The whole FSL is not installed.
#   - config_file:              Full path to the FSL configuration script to be generated

INSTDIR=$1
VERSION=$2
FSOS=$3
STANDARDONLY=$4
CONFIGFILE=$5

function fslconfig_bash {
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:INSTDIR/fsl/fslpython/lib:INSTDIR/fsl/lib'
    echo 'export FSLDIR=INSTDIR/fsl'
    echo 'source $FSLDIR/etc/fslconf/fsl.sh'
    echo 'export PATH=$FSLDIR/bin:$PATH'
}

cd ${INSTDIR}
echo "Installing fsl-${VERSION}-${FSOS} in ${INSTDIR}"
wget -q http://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-${VERSION}-${FSOS}_64.tar.gz
if [[ "x${STANDARDONLY}x" == "x1x" ]]; then
    tar xzf fsl-${VERSION}-${FSOS}_64.tar.gz fsl/data/standard
else
    tar xzf fsl-${VERSION}-${FSOS}_64.tar.gz
    ${INSTDIR}/fsl/etc/fslconf/post_install.sh -f ${INSTDIR}/fsl

    # config script
    if [[ $(basename $(echo $0)) == "bash" ]]; then
        CFGSTR=$(fslconfig_bash)
        echo "${CFGSTR//INSTDIR/$INSTDIR}" > $CONFIGFILE
    else
        echo "FSL configuration is implemented only for BASH."
        echo "FSL has not been installed properly. Exiting..." >&2
        exit -1; 
    fi
fi
rm fsl-${VERSION}-${FSOS}_64.tar.gz

if [[ -d "${INSTDIR}/fsl" ]]; then
    echo "FSL has been installed." >&2
else
    echo "FSL has not been installed properly. Exiting..." >&2
    exit -1; 
fi