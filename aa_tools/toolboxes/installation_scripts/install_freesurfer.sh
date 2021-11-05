# This is a script to download and install FreeSurfer with mri_deface templates. The scripts checks the instalation by looking for the "freesurfer" folder
#
# source install_freesurfer.sh <installation_directory> <version> <os> <license_string> <template_directory>
#
# Input aguments
#   - installation_directory:   This is where freesurfer will be installed. When configuring aap.directory_conventions.freesurferdir MUST point to '<installation_directory>/freesurfer', 
#                               and aap.directory_conventions.freesurferenvironment MUST point to '<installation_directory>/freesurfer/FreeSurferEnv.sh;' (i.e. with semicolon at the end)
#   - version:                  Version number of the desired installation according to semantic versioning (x.y.z)
#   - os:                       ID of OS the downloaded version has been compiled for
#   - license_string:           The content of the license obtained from FreeSurfer after registration. Please, beaware that it is usually a multiline text, and line breaks (\n) and whitespaces should be also added. 
#                               The example may be: "name@email.com\n0123n *Abcd4c5EFGhI\n FS67jKL8mN.9o"
#   - template_directory:       This is where the FreeSurfer deface templates will be downloaded. It MUST be the same as in aap.directory_conventions.templatedir. 
#                               Please, be aware that template files will be prefixed with "freesurfer_deface_" for better visibility.

INSTDIR=$1
VERSION=$2
FSOS=$3
LIC=$4
TPLDIR=$5

cd ${INSTDIR}
wget -q "https://ftp.nmr.mgh.harvard.edu/pub/dist/freesurfer/${VERSION}/freesurfer-linux-${FSOS}_x86_64-${VERSION}.tar.gz"
tar xzf "freesurfer-linux-${FSOS}_x86_64-${VERSION}.tar.gz"
rm "freesurfer-linux-${FSOS}_x86_64-${VERSION}.tar.gz"
mkdir -p ${TPLDIR}
wget -q -P ${TPLDIR} https://surfer.nmr.mgh.harvard.edu/pub/dist/mri_deface/talairach_mixed_with_skull.gca.gz
wget -q -P ${TPLDIR} https://surfer.nmr.mgh.harvard.edu/pub/dist/mri_deface/face.gca.gz
gunzip ${TPLDIR}/*
mv ${TPLDIR}/talairach_mixed_with_skull.gca ${TPLDIR}/freesurfer_deface_talairach_mixed_with_skull.gca
mv ${TPLDIR}/face.gca ${TPLDIR}/freesurfer_deface_face.gca

if [[ -d "${INSTDIR}/freesurfer" ]]; then
    echo "FreeSurfer has been installed." >&2
else
    echo "FreeSurfer has not been installed properly. Exiting..." >&2
    exit -1; 
fi