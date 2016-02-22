/*
 * Christian Gaser
 * $Id: sanlmMex.c 404 2011-04-11 10:03:40Z gaser $ 
 *
 */

#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include "matrix.h"

extern void anlm(float* ima, int v, int f, int rician, const int* dims);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/* Declarations */
float *ima;
int v,f,ndim,rician;
const int *dims;

/* check inputs */
if (nrhs!=3)
  mexErrMsgTxt("3 inputs required.");
else if (nlhs>0)
  mexErrMsgTxt("No output arguments allowed.");
  
if (!mxIsSingle(prhs[0]))
	mexErrMsgTxt("First argument must be float.");


/* get input image */
ima = (float*)mxGetPr(prhs[0]);

ndim = mxGetNumberOfDimensions(prhs[0]);
if (ndim!=3)
  mexErrMsgTxt("Images does not have 3 dimensions.");
  
dims = mxGetDimensions(prhs[0]);

/* get parameters */
v = (int)(mxGetScalar(prhs[1]));
f = (int)(mxGetScalar(prhs[2]));

rician = 1;

anlm(ima, v, f, rician, dims); 

return;

}

