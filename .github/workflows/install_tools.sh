TOOLDIR=$HOME/tools
mkdir $HOME/tools
TEMPLATEDIR=$HOME/tools/templates

sudo apt-get update
sudo apt-get install libtinfo5 libtinfo6 dc

source $GITHUB_WORKSPACE/aa_tools/toolboxes/installation_scripts/install_spm.sh $TOOLDIR 12

source $GITHUB_WORKSPACE/aa_tools/toolboxes/installation_scripts/install_facemasking.sh $TOOLDIR

echo "FSL: ${LOAD_FSL}; FREESURFER: ${LOAD_FREESURFER}"

if [[ "x${LOAD_FSL}x" == "x1x" ]]; then
    source $GITHUB_WORKSPACE/aa_tools/toolboxes/installation_scripts/install_fsl.sh $TOOLDIR 6.0.5 centos7
else # FSL standards are needed for some modules
    source $GITHUB_WORKSPACE/aa_tools/toolboxes/installation_scripts/install_fsl.sh $TOOLDIR 6.0.5 centos7 1
fi

if [[ "x${LOAD_FREESURFER}x" == "x1x" ]]; then
    source $GITHUB_WORKSPACE/aa_tools/toolboxes/installation_scripts/install_freesurfer.sh $TOOLDIR 7.2.0 centos7 "tibor.auer@gmail.com\n7061\n *Ccpi6x7PAIeQ\n FS96pPK5vW.0g" $TEMPLATEDIR
fi