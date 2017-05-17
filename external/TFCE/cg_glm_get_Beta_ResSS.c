#ifndef lint
static char sccsid[]="@(#)cg_glm_get_Beta_ResSS.c  1.01 Christian Gaser 07/09/09";
#endif

#include <math.h>
#include <memory.h>
#include "mex.h"

#include "spm_mapping.h"

#ifdef _OPENMP
#include "omp.h"
#endif

typedef struct{
    int n_subj;
    int n_beta;
    int ini;
    int fin;
    double* Beta;   
    double* ResSS;
    double* X;
    double* pKX;
    double* TH;
    double* W;
    MAPTYPE* maps;    
    MAPTYPE* map_mask;    
} myargument;

void ThreadFunc( myargument arg )
{
  int i, j, k, z, ind1, ind2, n_slices_x_values;
  int n_subj, n_values, n_slices, n_beta, ini, fin;
  double *estimates, *image, *mask;
  double *Beta, *ResSS, *X, *pKX, *W, *TH;
  MAPTYPE *maps, *map_mask;
  double sum, ival;
  double mat[] = {1, 0, 0, 0, 0, 1, 0, 0,  0, 0, 1, 0, 0, 0, 0, 1};

  n_subj = arg.n_subj;
  n_beta = arg.n_beta;   
  Beta = arg.Beta;
  ResSS = arg.ResSS;  
  maps = arg.maps;
  map_mask = arg.map_mask;    
  n_values = maps[0].dim[0]*maps[0].dim[1];
  n_slices = maps[0].dim[2];
  X = arg.X;    
  pKX = arg.pKX;
  TH = arg.TH;  
  W = arg.W;
  ini = arg.ini;
  fin = arg.fin;
  
#pragma omp critical (MALLOC)
{
  estimates = (double *)malloc(n_values*n_subj*sizeof(double));
  image     = (double *)malloc(n_values*sizeof(double));
  mask      = (double *)malloc(n_values*sizeof(double));
}

  n_slices_x_values = n_slices*n_values;

  for(z=ini; z<fin; z++) 
  {
    mat[14] = z + 1.0;  
    ind2 = z*n_values;  

    /* load mask */  
    slice(mat, mask, map_mask[0].dim[0],map_mask[0].dim[1], &map_mask[0], 0, 0.0);
    
    for(i=0; i<n_subj; i++)
    {
      slice(mat, image, maps[i].dim[0],maps[i].dim[1], &maps[i], 0, 0.0);
      for(j=0; j<n_values; j++)
      {
        ival = W[i]*image[j];
        if ((mask[j] > 0) & (ival>TH[i]))
        {
          /* initialize estimates with image values */
          estimates[j + (i*n_values)] = ival;
          /* calculate betas */
          for(k=0; k<n_beta; k++)
            Beta[j + ind2 + (k*n_slices_x_values)] += pKX[k + (i*n_beta)] * ival;
        }
      }
    }

    /* get estimates */
    for(i=0; i<n_subj; i++)
    {
      ind1 = i*n_values;
      for(j=0; j<n_values; j++)
      {
        ival = W[i]*image[j];
        if ((mask[j] > 0) & (ival>TH[i]))
        {
          sum = 0.0;
          /* calculate difference between estimates and original values */
          for(k=0; k<n_beta; k++)
            sum += (X[i + (k*n_subj)] * Beta[j + ind2 + (k*n_slices_x_values)]);
            estimates[j + ind1] -= sum;
        }
      }
    }

    /* calculate sum of residual squares */
    for(j=0; j<n_values; j++)
    {
      if (mask[j] > 0)
      {
        sum = 0.0;
        for(i=0; i<n_subj; i++)
        {
          ind1 = j + (i*n_values);  
             sum += estimates[ind1] * estimates[ind1];
        }
        ResSS[j + ind2] = sum;
      }
    }
  }
  
  free((char *)estimates);
  free((char *)mask);
  free((char *)image);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int n_subj, n_beta, n_slices, n_values, ini, fin;
  int n, i, j, k, z, Nthreads;
  double *Beta, *ResSS, *X, *pKX, *TH, *W;
  MAPTYPE *maps, *map_mask, *get_maps();

  myargument *ThreadArgs;  

  if (nrhs != 6) mexErrMsgTxt("Six input arguments required.");
  if (nlhs != 2) mexErrMsgTxt("Two output arguments required.");
 
  n_subj = mxGetM(prhs[2]);
  n_beta = mxGetN(prhs[2]);

  map_mask = get_maps(prhs[1], &n);
  if (n!=1)
  {
    free_maps(map_mask, n);
    mexErrMsgTxt("Only one single file as mask allowed.");
  }

  maps = get_maps(prhs[0], &n);
  if (n!=n_subj)
  {
    free_maps(maps, n);
    free_maps(map_mask, 1);
    mexErrMsgTxt("Different number of scans in design matrix.");
  }

  for(i=1; i<n_subj; i++)
  {
    if (  maps[i].dim[0] != maps[0].dim[0] ||
      maps[i].dim[1] != maps[0].dim[1] ||
      maps[i].dim[2] != maps[0].dim[2])
      {
        free_maps(maps, n_subj);
        mexErrMsgTxt("Incompatible image dimensions.");
      }
  }

  n_slices = maps[0].dim[2];
  n_values = maps[0].dim[0]*maps[0].dim[1];
  
  if ((n_slices!=map_mask[0].dim[2]) || (n_values!=map_mask[0].dim[0]*map_mask[0].dim[1]))
  {
    free_maps(maps, n);
    free_maps(map_mask, 1);
    mexErrMsgTxt("Incompatible dimensions between mask and images.");
  }
  
  X   = (double*)mxGetPr(prhs[2]);
  pKX = (double*)mxGetPr(prhs[3]);
  TH  = (double*)mxGetPr(prhs[4]);
  W   = (double*)mxGetPr(prhs[5]);
  
  n_subj = mxGetM(prhs[2]);
  if (n_subj!=mxGetM(prhs[4]))
  {
    free_maps(maps, n);
    free_maps(map_mask, 1);
    mexErrMsgTxt("Incompatible dimensions of thresholds.");
  }

  for(i=0; i<n_subj; i++)
    if (mxIsInf(TH[i]))
      TH[i] = -1e15;
  
  plhs[0] = mxCreateDoubleMatrix(n_slices*n_values, n_beta, mxREAL);
  Beta    = (double*)mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(n_slices*n_values, 1, mxREAL);
  ResSS   = (double*)mxGetPr(plhs[1]);
  
  /* initialize Beta and ResSS with zeros */
  memset(ResSS, 0, sizeof(double)*n_slices*n_values);
  memset(Beta, 0, sizeof(double)*n_slices*n_values*n_beta);

  Nthreads = 1;

#ifdef _OPENMP
  Nthreads = omp_get_num_procs();
  omp_set_num_threads(Nthreads);
#endif

/* Reserve room for handles of threads in ThreadList */
ThreadArgs = (myargument*) malloc( Nthreads*sizeof(myargument));

#ifdef _OPENMP
  # pragma omp parallel for default(shared) private(i,ini,fin)
#endif
for (i = 0; i<Nthreads; i++)
{         
    /* Make Thread Structure */
    ini = (int)(i*n_slices)/Nthreads;
    fin = (int)((i+1)*n_slices)/Nthreads;            
    ThreadArgs[i].n_subj = n_subj;
    ThreadArgs[i].n_beta = n_beta;   
    ThreadArgs[i].Beta = Beta;
    ThreadArgs[i].ResSS = ResSS;  
    ThreadArgs[i].maps = maps;
    ThreadArgs[i].map_mask = map_mask;    
    ThreadArgs[i].X = X;    
    ThreadArgs[i].pKX = pKX;
    ThreadArgs[i].TH = TH;  
    ThreadArgs[i].W = W;
    ThreadArgs[i].ini = ini;
    ThreadArgs[i].fin = fin;
    (void)ThreadFunc(ThreadArgs[i]);    
}

  free(ThreadArgs);
  free_maps(maps, n_subj);
  free_maps(map_mask, 1);

}
