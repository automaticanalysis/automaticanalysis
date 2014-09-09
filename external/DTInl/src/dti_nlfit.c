/* Nonlienear fit of the DTI tensor on a single voxel timeseries
 FORMAT [S0 DT] = dti_nlfit(vdata, b, bv)
 vdata - single voxel DWI timeseries
 bvals - [1xN] b-values coresponding to vdata with N volumes
 bvecs - [3xN] DW directions coresponding to vdata with N volumes
 S0    - sphericity index
 DT    - [3x3] diffusion tensor (non-diagonalised)
_______________________________________________________________________

 Copyright (C) 2014 MRC Congition and Brain Sciences Unit
 Marta Correia and Tibor Auer
 $Id: dti_nlfit.c 2014-08-15 16:00:00Z ta02 $ */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Declare variables */
    size_t mrows,ncols,nVol;
    double *data, *bvals, *bvecs, *S0, *tensor;
    
    /* Check for proper number of arguments. */
    if(nrhs!=3)
        mexErrMsgIdAndTxt( "MATLAB:dti_nlfit:invalidNumInputs",
                "Three inputs required.\nType help dti_nlfit!");
    
    /* Check for proper format of arguments. */
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    if( !(mrows==1 < ncols==1) ) {
        mexErrMsgIdAndTxt( "MATLAB:dti_nlfit:inputNotVector",
                "First input must be a [1xN] vector.\nType help dti_nlfit!");
    }
    nVol = mxGetNumberOfElements(prhs[0]);
    
    mrows = mxGetM(prhs[1]);
    ncols = mxGetN(prhs[1]);
    if( (mrows!=1) || (ncols!=nVol) ) {
        mexErrMsgIdAndTxt( "MATLAB:dti_nlfit:invalidInput",
                "Second input must be a [1xN] vector with the same lenght as the first input.\nType help dti_nlfit!");
    }

    mrows = mxGetM(prhs[2]);
    ncols = mxGetN(prhs[2]);
    if( (mrows!=3) || (ncols!=nVol) ) {
        mexErrMsgIdAndTxt( "MATLAB:dti_nlfit:invalidInput",
                "Third input must be a [3xN] vector with the same lenght as the first input.\nType help dti_nlfit!");
    }
    
    /* Get the data */
    data=(double *)mxGetPr(prhs[0]);
    bvals=(double *)mxGetPr(prhs[1]);
    bvecs=(double *)mxGetPr(prhs[2]);
    
    /* Prepare tensor */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);
    S0 = (double *)mxGetPr(plhs[0]);
    tensor = (double *)mxGetPr(plhs[1]);
 
    get_diffusion_tensor(data, nVol, bvals, bvecs, S0, tensor); 
}