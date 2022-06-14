#!/bin/bash
#
# Plots conservative variables (-c) without Plot3D normalization (-np)
#
LNS3D_DIR=${LNS3D_DIR:=$HOME/git/lns3d}
usage() {
  echo "Usage: lpost3d.sh [-i] [-r] [-f] [-t] [-ij] [-ji] [-h] [-v] filename(s)"
  echo "  -i              output imaginary component"
  echo "  -r              output real component"
  echo "  -t              output at a specified time"
  echo "  -ij             use IJ ordering on input/output"
  echo "  -ji             use JI ordering on input/output"
  echo "  -f              filter the field before output"
  echo "  -h | --help     this usage"
  echo "  -v | --verbose  verbose output"
  exit 2
}
# Note:  this requires that gnu-getopt be active in your path
# on my mac I needed to install using homebrew:  brew install gnu-getopt
# and make sure that /usr/local/opt/gnu-getopt/bin was first in my path
TEST_GETOPT=$(getopt --test)
RET_VALUE=$?
#echo "RET_VALUE = $RET_VALUE"
if [ "$RET_VALUE" != "4" ]; then
  echo "lpost3d.sh Requires GNU getopt, exiting"
  exit 1
fi
#
# Use GNU getopt
#
PARSED_ARGUMENTS=$(getopt -a -o hirftv --long ij,ji,help,verbose -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi
# set default values
COMPLEX="-i"
SWITCH="-ij"
FILTER=
TIME=
#echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -h | --help)  usage         ; shift ;; 
    -i)           COMPLEX="-i"  ; shift ;;
    -r)           COMPLEX=""    ; shift ;;
    --ij)         SWITCH="-ij"  ; shift ;;
    --ji)         SWITCH=""     ; shift ;;
    -f)           FILTER="-f"   ; shift ;;
    -t)           TIME="-t"     ; shift ;;
    -v | --verbose) VERBOSE=1   ; shift ;;
    # -- means the end of the arguments; drop this and break out of loop
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
   esac
done
#echo "COMPLEX     : $COMPLEX"
#echo "SWITCH      : $SWITCH"
#echo "Remaining parameters are: $@"
for file in $@
do
  echo "Processing file $file"
  if [ "$VERBOSE" == 1 ]; then
    echo "Running: $LNS3D_DIR/util/lpost3d $COMPLEX $SWITCH $FILTER "\
         "$TIME $file >> lpost3d.log"
  fi
  $LNS3D_DIR/util/lpost3d $COMPLEX $SWITCH $FILTER $TIME $file >> lpost3d.log 
done
exit 0
