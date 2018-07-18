#!/bin/sh
# This wrapper script is intended to support independent execution.
# 
# This script uses the following environment variables set by the submit MATLAB code:
# MDCE_MATLAB_EXE     - the MATLAB executable to use
# MDCE_MATLAB_ARGS    - the MATLAB args to use
#

# Copyright 2010-2011 The MathWorks, Inc.

export MALLOC_ARENA_MAX=4

echo "Executing: ${MDCE_MATLAB_EXE} ${MDCE_MATLAB_ARGS}"
exec "${MDCE_MATLAB_EXE}" ${MDCE_MATLAB_ARGS}
