TOOLDIR=$HOME/tools
mkdir $TOOLDIR
mkdir $TOOLDIR/config
TEMPLATEDIR=$TOOLDIR/templates

sudo apt-get update
sudo apt-get install libtinfo5 libtinfo6 dc libxml2-utils
curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o $TOOLDIR/get-pip.py
python $TOOLDIR/get-pip.py --force-reinstall

# All MATLAB tools MUST be installed referred by the parameterset
source $GITHUB_WORKSPACE/.github/workflows/tools_urls.sh
source $GITHUB_WORKSPACE/aa_tools/toolboxes/installation_scripts/install_tools.sh $GITHUB_WORKSPACE/.github/workflows/$PARAMETER_XML $TOOLDIR

echo "FSL: ${LOAD_FSL}; FREESURFER: ${LOAD_FREESURFER}"

if [[ "x${LOAD_FSL}x" == "x1x" ]]; then
    source $GITHUB_WORKSPACE/aa_tools/toolboxes/installation_scripts/install_fsl.sh $TOOLDIR 6.0.5 centos7 0 $TOOLDIR/config/fsl_bash.sh
fi

if [[ "x${LOAD_FREESURFER}x" == "x1x" ]]; then
    source $GITHUB_WORKSPACE/aa_tools/toolboxes/installation_scripts/install_freesurfer.sh $TOOLDIR 7.2.0 centos7 "tibor.auer@gmail.com\n7061\n *Ccpi6x7PAIeQ\n FS96pPK5vW.0g" $TEMPLATEDIR
fi

echo "Free space:"
df -h
mkdir $HOME/.aa
cp $GITHUB_WORKSPACE/.github/workflows/$PARAMETER_XML $HOME/.aa/aap_parameters_user.xml
mkdir $HOME/projects