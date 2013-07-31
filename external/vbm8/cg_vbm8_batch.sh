#! /bin/sh

########################################################
# global parameters
########################################################
# $Id: cg_vbm8_batch.sh 404 2011-04-11 10:03:40Z gaser $

matlab=matlab   # you can use other matlab versions by changing the matlab parameter
writeonly=0
defaults_file=""
CPUINFO=/proc/cpuinfo
ARCH=`uname`

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  check_matlab
  check_files
  get_no_of_cpus
  run_vbm

  exit 0
}


########################################################
# check arguments and files
########################################################

parse_args ()
{
  local optname optarg
  count=0

  if [ $# -lt 1 ]; then
    help
    exit 1
  fi

  while [ $# -gt 0 ]
  do
    optname="`echo $1 | sed 's,=.*,,'`"
    optarg="`echo $2 | sed 's,^[^=]*=,,'`"
    case "$1" in
        --m* | -m*)
            exit_if_empty "$optname" "$optarg"
            matlab=$optarg
            shift
            ;;
        --d* | -d*)
            exit_if_empty "$optname" "$optarg"
            defaults_file=$optarg
            shift
            ;;
        --p* | -p*)
            exit_if_empty "$optname" "$optarg"
            NUMBER_OF_JOBS=$optarg
            shift
            ;;
        --w* | -w*)
            writeonly=1
            ;;
        -h | --help | -v | --version | -V)
            help
            exit 1
            ;;
        -*)
            echo "`basename $0`: ERROR: Unrecognized option \"$1\"" >&2
            ;;
        *)
            ARRAY[$count]=$1
            ((count++))
            ;;
    esac
    shift
  done

}

########################################################
# check arguments
########################################################

exit_if_empty ()
{
  local desc val

  desc="$1"
  shift
  val="$*"

  if [ -z "$val" ]
  then
    echo 'ERROR: No argument given with \"$desc\" command line argument!' >&2
    exit 1
  fi
}

########################################################
# check files
########################################################

check_files ()
{
  
  SIZE_OF_ARRAY="${#ARRAY[@]}"
  if [ "$SIZE_OF_ARRAY" -eq 0 ]
  then
      echo 'ERROR: No files given!' >&2
      help
      exit 1
  fi

  i=0
  while [ "$i" -lt "$SIZE_OF_ARRAY" ]
  do
    if [ ! -f "${ARRAY[$i]}" ]; then
      echo ERROR: File ${ARRAY[$i]} not found
      help
      exit 1
    fi
    ((i++))
  done

}

########################################################
# get # of cpus
########################################################
# modified code from
# PPSS, the Parallel Processing Shell Script
# 
# Copyright (c) 2009, Louwrentius
# All rights reserved.

get_no_of_cpus () {

  if [ -z "$NUMBER_OF_JOBS" ];
  then
    if [ "$ARCH" == "Linux" ]
    then
      NUMBER_OF_JOBS=`grep ^processor $CPUINFO | wc -l`

    elif [ "$ARCH" == "Darwin" ]
    then
      NUMBER_OF_JOBS=`sysctl -a hw | grep -w logicalcpu | awk '{ print $2 }'`

    elif [ "$ARCH" == "FreeBSD" ]
    then
      NUMBER_OF_JOBS=`sysctl hw.ncpu | awk '{ print $2 }'`

    else
      NUMBER_OF_JOBS=`grep ^processor $CPUINFO | wc -l`

    fi
    echo "Found $NUMBER_OF_JOBS processors."

    if [ -z "$NUMBER_OF_JOBS" ]
    then
        echo "$FUNCNAME ERROR - number of CPUs not obtained. Use -p to define number of processes."
        exit 1
    fi
  fi
}

########################################################
# run vbm tool
########################################################

run_vbm ()
{
    cwd=`dirname $0`
    pwd=$PWD
    
    # we have to go into toolbox folder to find matlab files
    cd $cwd
    
    spm8=`dirname $cwd`
    spm8=`dirname $spm8`

    export MATLABPATH=$spm8

    SIZE_OF_ARRAY="${#ARRAY[@]}"
    BLOCK=$((10000* $SIZE_OF_ARRAY / $NUMBER_OF_JOBS ))
    
    # argument empty?
    if [ ! "${defaults_file}" == "" ]; then
        # check wether absolute or relative names were given
        if [ ! -f ${defaults_file} ];  then
            defaults_file=${pwd}/${defaults_file}
        fi
    
        # check whether defaults file exist
        if [ ! -f ${defaults_file} ];  then
            echo $defaults_file not found.
        fi
    fi

    # split files and prepare tmp-file with filenames
    TMP=/tmp/vbm8_$$
    i=0
    while [ "$i" -lt "$SIZE_OF_ARRAY" ]
    do
        count=$((10000* $i / $BLOCK ))
        
        # check wether absolute or relative names were given
        if [ ! -f ${ARRAY[$i]} ];  then
            FILE=${pwd}/${ARRAY[$i]}
        else
            FILE=${ARRAY[$i]}
        fi
        if [ -z "${ARG_LIST[$count]}" ]; then
            ARG_LIST[$count]="$FILE"
        else
            ARG_LIST[$count]="${ARG_LIST[$count]} $FILE"
        fi
        echo ${FILE} >> ${TMP}${count}
        ((i++))
    done
    time=`date "+%Y%b%d_%H%M"`
    vbmlog=${pwd}/vbm8_${time}
    
    i=0
    while [ "$i" -lt "$NUMBER_OF_JOBS" ]
    do
        if [ ! "${ARG_LIST[$i]}" == "" ]; then
            j=$(($i+1))
            COMMAND="cg_vbm8_batch('${TMP}${i}',${writeonly},'${defaults_file}')"
            echo Calculate ${ARG_LIST[$i]}
            echo ---------------------------------- >> ${vbmlog}_${j}.log
            date >> ${vbmlog}_${j}.log
            echo ---------------------------------- >> ${vbmlog}_${j}.log
            echo >> ${vbmlog}_${j}.log
            echo $0 $file >> ${vbmlog}_${j}.log
            echo >> ${vbmlog}_${j}.log
            nohup ${matlab} -nodisplay -nojvm -nosplash -r $COMMAND >> ${vbmlog}_${j}.log 2>&1 &
            echo Check ${vbmlog}_${j}.log for logging information
            echo
        fi
        ((i++))
    done

    exit 0
}

########################################################
# check if matlab exist
########################################################

check_matlab ()
{
  found=`which ${matlab} 2>/dev/null`
  if [ ! -n "$found" ];then
    echo $matlab not found.
    exit 1
  fi
}

########################################################
# help
########################################################

help ()
{
cat <<__EOM__

USAGE:
   cg_vbm8_batch.sh filename|filepattern [-m matlabcommand] [-w]
   
   -m   matlab command
   -p   number of parallel jobs (=number of processors)
   -w		write already segmented images
   -d   optional default file
   
   Only one filename or pattern is allowed. This can be either a single file or a pattern
   with wildcards to process multiple files. Optionally you can set the matlab command 
   with the "-m" option and force to write already estimated segmentations with the "-w" option.

PURPOSE:
   Command line call of VBM8 segmentation

EXAMPLE
   cg_vbm8_batch.sh spm/spm8/canonical/single_subj_T1.nii
   This command will process only the single file single_subj_T1.nii. 
   
   cg_vbm8_batch.sh spm/spm8/canonical/single_subj_T1.nii -d your_vbm8_defaults_file.m
   This command will process only the single file single_subj_T1.nii. The defaults defined
   in your_vbm8_defaults_file.m will be used instead of cg_vbm8_defaults.m.

   cg_vbm8_batch.sh spm/spm8/canonical/*152*.nii
   Using wildcards all files containing the term "152" will be processed. In this case these 
   are the files avg152PD.nii, avg152T1.nii, and avg152T2.nii.

   cg_vbm8_batch.sh spm/spm8/canonical/*152*.nii -m /usr/local/bin/matlab7
   Using wildcards all files containing the term "152" will be processed. In this case these 
   are the files avg152PD.nii, avg152T1.nii, and avg152T2.nii.
   As matlab-command /usr/local/bin/matlab7 will be used.

INPUT:
   analyze or nifti files

OUTPUT:
   segmented images according to settings in cg_vbm8_defaults.m
   vbm8_log_$time.txt for log information

USED FUNCTIONS:
   cg_vbm8_batch.m
   VBM8 toolbox
   SPM8

SETTINGS
   matlab command: $matlab
   
This script was written by Christian Gaser (christian.gaser@uni-jena.de).

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}
