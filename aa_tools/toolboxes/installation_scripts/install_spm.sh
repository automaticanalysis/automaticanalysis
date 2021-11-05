# This is a script to download SPM. The scripts checks the instalation by looking for the "spm<version>" folder
#
# source install_spm.sh <installation_directory> <version>
#
# Input aguments
#   - installation_directory:   This is where SPM will be installed. When configuring aap.directory_conventions.toolbox.dir MUST point to '<installation_directory>/spm<version>'
#   - version:                  Version number of the desired installation (e.g. 99, 2, 5, 8, 12).

INSTDIR=$1
VERSION=$2

cd ${INSTDIR}
git clone "https://github.com/spm/spm${VERSION}"

if [[ -d "${INSTDIR}/spm${VERSION}" ]]; then
    echo "SPM${VERSION} has been installed." >&2
else
    echo "SPM${VERSION} has not been installed properly. Exiting..." >&2
    exit -1; 
fi