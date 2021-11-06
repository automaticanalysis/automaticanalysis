# This is a script to download MATLAB tools also used in the tests. It requires a set of variables (URL_<tool name with lower case>) describing the source.
# See <aa repo>/.github/workflows/tools_urls.sh as examples.
# The script also checks the installation by looking for the corresponding folders.
#
# source install_matlabtools.sh <parameterset> <installation_directory>
#
# Input aguments
#   - parameterset:             XML parameterset containing information of required tools
#   - installation_directory:   This is where FaceMasking will be installed. When configuring aap.directory_conventions.FaceMaskingdir MUST point to '<installation_directory>/mask_face/nrg-improc'

PARAMXML=$1
INSTDIR=$2

function install_tool() {
    name=$1
    url=$2
    folder=$3

    if [[ $url == git+* ]]; then
        git clone ${url:4} 
    fi

    if [[ -d "${folder}" ]]; then
        echo "${name} has been installed." >&2
    else
        echo "${name} has not been installed properly. Exiting..." >&2
        exit -1; 
    fi
}

CWD=$PWD
cd ${INSTDIR}

# Toolboxes
tbxCount=$(xmllint --xpath 'count(//aap/local/directory_conventions/toolbox)' $PARAMXML)
for indTbx in `seq $tbxCount`; do
    tbxName=$(xmllint --xpath '//aap/local/directory_conventions/toolbox['$indTbx']/name/text()' $PARAMXML)
    varURL=URL_${tbxName}
    tbxDir=$(xmllint --xpath '//aap/local/directory_conventions/toolbox['$indTbx']/dir/text()' $PARAMXML)
    install_tool $tbxName ${!varURL} $tbxDir
done

# Dedicated tools
# - FaceMasking from WUSLT NRG
if [[ $(xmllint --xpath 'count(//aap/local/directory_conventions/FaceMaskingdir)' $PARAMXML) > 0 ]]; then
    install_tool FaceMasking $URL_facemasking $(xmllint --xpath '//aap/local/directory_conventions/FaceMaskingdir/text()' $PARAMXML)
fi

cd $CWD