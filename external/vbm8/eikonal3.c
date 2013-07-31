/*
 * Robert Dahnke
 * $Id: eikonal3.c 419 2011-06-07 14:13:48Z gaser $ 
 *
 */

/* eikonal distance calculation
 * _____________________________________________________________________________
 * Estimates eikonal distance to an object (negative values) with the speed given
 * by the positive values
 *
 *  D = eqdist(F[,vx_vol)
 *
 *  F      (single) object and speed map
 *  vx_vol (double) size of the voxel (default: [1 1 1])
 * 
 * _____________________________________________________________________________
 * Robert Dahnke 2010_01
 * Center of Neuroimaging 
 * University Jena
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

float abs2(float n) { if (n<0) return -n; else return n; }


/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<1)                                       mexErrMsgTxt("ERROR:eikonal3: not enough input elements\n");
  if (nrhs>2)                                       mexErrMsgTxt("ERROR:eikonal3: too many input elements.\n");
  if (nlhs>2)                                       mexErrMsgTxt("ERROR:eikonal3: too many output elements.\n");
  if (mxIsSingle(prhs[0])==0)                       mexErrMsgTxt("ERROR:eikonal3: first  input must be an 3d single matrix\n");
  if (nrhs==2 && mxIsDouble(prhs[1])==0)            mexErrMsgTxt("ERROR:eikonal3: second input must be an double matrix\n");
  if (nrhs==2 && mxGetNumberOfElements(prhs[1])!=3) mexErrMsgTxt("ERROR:eikonal3: second input must have 3 Elements");
  
  int i, j, k, n, ind, x, y, z; 
  int ni, kll=0;
  float diff, maxdiffi, maxdiff=1, Dio, DNi;
  float TH=0.01;

  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]); 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = (int)mxGetNumberOfElements(prhs[0]);
  const int sS[] = {1,3}; 
  mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  double *S = mxGetPr(SS);
  if (nrhs<2) {S[0]=1; S[1]=1; S[2]=1;} else {S = mxGetPr(prhs[1]);}
  
  float s1 = abs2((float)S[0]),s2 = abs2((float)S[1]),s3 = abs2((float)S[2]);
  const float   s12  = sqrt( s1*s1  + s2*s2); /* xy - voxel size */
  const float   s13  = sqrt( s1*s1  + s3*s3); /* xz - voxel size */
  const float   s23  = sqrt( s2*s2  + s3*s3); /* yz - voxel size */
  const float   s123 = sqrt(s12*s12 + s3*s3); /* xyz - voxel size */
        
  /* indices of the euclidean distance NW */
  const float ND[]  = {s123, s12, s123, s13, s1, s13, s123, s12, s123,   s23, s2, s23, s3, 0.0, s3, s23, s2, s23,   s123, s12, s123, s13, s1, s13, s123, s12, s123};

  /* main volumes - actual without memory optimation ... */
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  float *D = (float *)mxGetPr(plhs[0]);  
  float *L = (float *)mxGetPr(plhs[1]);  
    
  /* input variables */
  float *SEG  = (float *)mxGetPr(prhs[0]);

  if ( TH>=0.5 || TH<0.0001 ) mexErrMsgTxt("ERROR:eikonal3: threshhold must be >0.0001 and smaller than 0.5\n");
  
  /* intitialisiation */
  for (i=0;i<nL;i++) {
    if ( SEG[i]<0 ) {D[i]=0; L[i]=SEG[i];} 
    else            {D[i]=FLT_MAX; L[i]=0;} 
  }
  
  while ( ( maxdiff > TH ) && kll<2000 ) {
    maxdiffi=0;
    kll++;
    
    for (z=0;z<sL[2];z++) for (y=0;y<sL[1];y++) for (x=0;x<sL[0];x++) {
      ind = index(x,y,z,sL);
      if ( SEG[ind]>0 && SEG[ind]<FLT_MAX) {
        /* read neighbor values */
        Dio = D[ind];
        
        n = 0;
        /* go through all elements in a 3x3x3 box */
        for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (k=-1;k<=1;k++) {
          ni = index(x+i,y+j,z+k,sL);

          if ( SEG[ind]!=FLT_MAX ) {
            if ( ((x+i)>=0) && ((x+i)<sL[0]) && ((y+j)>=0) && ((y+j)<sL[1]) && ((z+k)>=0) && ((z+k)<sL[2]) ) {
               DNi = D[ni] + ND[n]*SEG[ind];
               if ( DNi<D[ind] ) {D[ind]=DNi; L[ind]=L[ni]; } 
            } 
          }
          n++;
        }
        
        diff  = abs2( Dio - D[ind] );
        if ( maxdiffi<diff ) maxdiffi=diff;  
      }
    }
    
    for (z=sL[2]-1;z>=0;z--) for (y=sL[1]-1;y>=0;y--) for (x=sL[0]-1;x>=0;x--) {
      ind = index(x,y,z,sL);
      if ( SEG[ind]>0 && SEG[ind]<FLT_MAX) {
        /* read neighbor values */
        Dio = D[ind];
        
        n = 0;
        /* go through all elements in a 3x3x3 box */
        for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (k=-1;k<=1;k++) {
          ni = index(x+i,y+j,z+k,sL);

          if ( SEG[ind]!=FLT_MAX ) {
            if ( ((x+i)>=0) && ((x+i)<sL[0]) && ((y+j)>=0) && ((y+j)<sL[1]) && ((z+k)>=0) && ((z+k)<sL[2]) ) {
               DNi = D[ni] + ND[n]*SEG[ind];
               if ( DNi<D[ind] ) {D[ind]=DNi; L[ind]=L[ni]; } 
            } 
          }
          n++;
        }
        
        diff  = abs2( Dio - D[ind] );
        if ( maxdiffi<diff ) maxdiffi=diff;  
      }
    }

    maxdiff = maxdiffi;   
  }
}