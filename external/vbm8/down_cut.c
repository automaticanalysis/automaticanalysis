/*
 * Robert Dahnke
 * $Id: down_cut.c 419 2011-06-07 14:13:48Z gaser $ 
 *
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"

#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

float abs2(float n) { if (n<0) return -n; else return n; }
float sign(float n) { if (n<0) return 1; else return 0; }
float max2(float a, float b) { if (a>b) return a; else return b; }

/* MAINFUNCTION */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<1)                                                      mexErrMsgTxt("ERROR:down_cut01: not enough input elements\n");
  if (nrhs>5)                                                      mexErrMsgTxt("ERROR:down_cut01: too many input elements.\n");
  if (nlhs>2)                                                      mexErrMsgTxt("ERROR:down_cut01: too many output elements.\n");
  if (mxIsSingle(prhs[0])==0)                                      mexErrMsgTxt("ERROR:down_cut01: first input must be an 3d single matrix\n");
  if (mxIsSingle(prhs[1])==0)                                      mexErrMsgTxt("ERROR:down_cut01: second input must be an 3d single matrix\n");
  if (mxIsDouble(prhs[2])==0 || mxGetNumberOfElements(prhs[2])!=1) mexErrMsgTxt("ERROR:down_cut01: third input must one double value\n");
  if (nrhs==4 && mxIsDouble(prhs[3])==0)                           mexErrMsgTxt("ERROR:down_cut01: fourth input must be an double matrix\n");
  if (nrhs==4 && mxGetNumberOfElements(prhs[3])!=3)                mexErrMsgTxt("ERROR:down_cut01: fourth input must have 3 Elements");
  if (nrhs==5 && mxIsDouble(prhs[4])==0)                           mexErrMsgTxt("ERROR:down_cut01: fifth input must be an double matrix\n");
  if (nrhs==5 && mxGetNumberOfElements(prhs[4])!=2)                mexErrMsgTxt("ERROR:down_cut01: fifth input must have 2 Elements");
    
  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = (int)mxGetNumberOfElements(prhs[0]);

  const int sSS[] = {1,3}, sdsv[] = {1,2}; 
  mxArray *SS  = mxCreateNumericArray(2,sSS, mxDOUBLE_CLASS,mxREAL); double*S  = mxGetPr(SS);
  mxArray *dsv = mxCreateNumericArray(2,sdsv,mxDOUBLE_CLASS,mxREAL); double*dd = mxGetPr(dsv);
  float dI = 0.0;
  double *SEGd;
  if (nrhs>=3) {SEGd=mxGetPr(prhs[2]); dI=(float) SEGd[0];}; 
  if (nrhs<4)  {S[0]=1; S[1]=1; S[2]=1;} else {S=mxGetPr(prhs[3]);}
  if (nrhs<5)  {dd[0]=0.1; dd[1]=10;}    else {dd=mxGetPr(prhs[4]);}
  
  float s1 = abs2((float)S[0]),s2 = abs2((float)S[1]),s3 = abs2((float)S[2]);
  const float   s12  = sqrt( s1*s1  + s2*s2); /* xy - voxel size */
  const float   s13  = sqrt( s1*s1  + s3*s3); /* xz - voxel size */
  const float   s23  = sqrt( s2*s2  + s3*s3); /* yz - voxel size */
  const float   s123 = sqrt(s12*s12 + s3*s3); /* xyz - voxel size */
  
  /* indices of the euclidean distance NW */
  const float ND[]  = {s123, s12, s123, s13, s1, s13, s123, s12, s123,   s23, s2, s23, s3, 0.0, s3, s23, s2, s23,   s123, s12, s123, s13, s1, s13, s123, s12, s123};
 
  int ind,i,j,k,x,y,z,n,ni;
    
  /* main volumes - actual without memory optimation ... */
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); /* label map */
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); /* tissue map (speed) */
    
  /* input variables */
  float*ALAB = (float *)mxGetPr(prhs[0]); /* label map */
  float*SEG  = (float *)mxGetPr(prhs[1]); /* tissue map */
  
  /* output variables */
  float*SLAB = (float *)mxGetPr(plhs[0]); /* label map */
  float*DIST = (float *)mxGetPr(plhs[1]); /* distance map */
  
  int   nCV = 0;    /* # voxel of interest (negative voxel that have to processed) */
  int   kll;
  int   kllv = 2000;
  float DISTN;
  
  /* initialisation of parameter volumes */
  for (i=0;i<nL;i++) { 
    SLAB[i] = ALAB[i]; 
    if (isnan(SLAB[i])) SLAB[i] = 0;
    if (SLAB[i]==0)     DIST[i] = FLT_MAX; 
    else {
      if (SLAB[i]==-FLT_MAX) DIST[i] = -FLT_MAX;
      else                  {DIST[i] = 0; nCV++;} 
    } 
  }
  
  /* diffusion */
  int   nC = nCV;
  kll=0;
  while ( nCV>0 && kll<kllv && nC>0 ) {
    kll++; nC=0;

    for (z=0;z<sL[2];z++) for (y=0;y<sL[1];y++) for (x=0;x<sL[0];x++) {
      ind = index(x,y,z,sL);

      if ( (DIST[ind]<=0) && (DIST[ind] != -FLT_MAX ) ) { 
        if (DIST[ind]<0) DIST[ind] = -DIST[ind]; 
        nCV--; /* demark points - also the with zero distance */
        
        n = 0;
        /* go through all elements in a 3x3x3 box */
        for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (k=-1;k<=1;k++) {
          ni = index(x+i,y+j,z+k,sL);
          if ( ((x+i)>=0) && ((x+i)<sL[0]) && ((y+j)>=0) && ((y+j)<sL[1]) && ((z+k)>=0) && 
               ((z+k)<sL[2]) && (SEG[ind]+dI>=SEG[ni]) && ALAB[ni]==0) {
                  
            if (nrhs==5) DISTN = DIST[ind] + dd[0]*ND[n] + dd[1]*max2(0,SEG[ni]);  
            else         DISTN = DIST[ind] + dd[0]*ND[n] + dd[1]*max2(0,4-SEG[ni]);

            if (  (DIST[ni]!=-FLT_MAX) && (abs2(DIST[ni])>abs2(DISTN)) )  { 
              if (DIST[ni]>0) nCV++; nC++;
              DIST[ni] = -DISTN;
              SLAB[ni] = SLAB[ind];
            }
          }
          n++;
        }
          

        if (DIST[ind]==0) DIST[ind] = -FLT_MAX; /* demark start points */   
      }
    }
  }
}