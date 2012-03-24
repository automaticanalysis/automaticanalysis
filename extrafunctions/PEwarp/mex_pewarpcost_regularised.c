#include <math.h>
#include "mex.h"
#include "blas.h"
#include "omp.h"

/* Eelke Visser, 2009 (eelke.visser@donders.ru.nl) */

/* #include <time.h> */

/* TODO: Is Ctrl-C from Matlab handled properly? */

void getDef(
    float *def,
    mwSignedIndex *dims,
    mwSignedIndex nbx,
    float *bx,
    mwSignedIndex nby,
    float *by,
    mwSignedIndex nbz,
    float *bz,
    float *x,
    float xscale)
{
    char cN = 'N';
    char cT = 'T';
    float fZero = 0;
    float fOne = 1;
    mwSignedIndex iOne = 1;
    
    mwSignedIndex Alen = dims[0] * nby;
    mwSignedIndex Blen = dims[0] * dims[1];
    int nvox = dims[0] * dims[1] * dims[2];
    
    float *A, *B, *p;
    int i;
    
    for (p = def; p < def + nvox; p++)
        *p = 0;
    
    A = mxMalloc(sizeof(float) * Alen);
    B = mxMalloc(sizeof(float) * Blen);
    
    for (p = A; p < A + Alen; p++)
        *p = 0;
    
    for (p = B; p < B + Blen; p++)
        *p = 0;
    
    for (i = 0; i < nbz; i++)
    {
        float *xi = x + nbx * nby * i;
        float *bzi = bz + dims[2] * i;
        
        /* A_k = scale * bx * x_k */
        sgemm_(&cN, &cN, &dims[0], &nby, &nbx, &xscale, bx, &dims[0], xi, &nbx, &fZero, A, &dims[0]);

        /* B_k = A_k * by^T */
        sgemm_(&cN, &cT, &dims[0], &dims[1], &nby, &fOne, A, &dims[0], by, &dims[1], &fZero, B, &dims[0]);
        
        /* def = def + B_k(linear) * Z_k */
        sger_(&Blen, &dims[2], &fOne, B, &iOne, bzi, &iOne, def, &Blen);
    }
    
    mxFree(A);
    mxFree(B);
}

void peResampleAndApplyJacobian(uint8_T *Y, float *X, float *D, /* float *J, */ uint8_T *mask, const mwSignedIndex *dims)
{
    mwSize sx = dims[0];
    mwSize sy = dims[1];
    mwSize sz = dims[2];
    int k;

    #pragma omp for schedule(dynamic) private(k)
    for (k = 0; k < sz; k++)
    {
        int offset = k * sx * sy;
        
        uint8_T *Yi = Y + offset;
        float *Xi = X + offset;
        float *Di = D + offset;
        /* float *Ji = J + offset; */
        uint8_T *mi = mask + offset;
        
        float *slcstart = X + offset;
        float *slcend = slcstart + sx * sy;
        
        int i, j;
        
        for (j = 0; j < sy; j++)
        {
            for (i = 0; i < sx; i++)
            {
                if (*mi)
                {
                    float intp;

                    /* Compiles and links, but fails with errno set on execution
                     * if math.h is not included. (???) */
                    float frac = modff(*Di, &intp);

                    /* Is this cast always correct? */
                    int iintp = (int)intp;

                    float *nearvox = Xi + iintp * sx;
                    float v;

                    if (nearvox < slcstart + sx)
                        v = *(slcstart + j);
                    else if (nearvox > slcend - sx)
                        v = *(slcend - sx + j);
                    else
                        if (frac < 0)
                            v = ((frac + 1) * *nearvox - frac * *(nearvox - sx)) /* * (*Ji + 1) */;
                        else
                            v = ((-frac + 1) * *nearvox + frac * *(nearvox + sx)) /* * (*Ji + 1) */;

                    if (v > 255)
                        *Yi = 255;
                    else
                        *Yi = (uint8_T)(v + 0.5f);
                }
                else
                    *Yi = 0; /* Not strictly necessary since these voxels are never read. */

                Di++;
                /* Ji++; */
                Yi++;
                Xi++;
                mi++;
            }
        }
    }
}

void hist2(float *H, uint8_T *F, uint8_T *G, uint8_T *mask, mwSize n)
{
    #pragma omp single
    {
        uint8_T *Fend = F + n;
        float *p;
    
        for (p = H; p < H + 65536; p++)
            *p = 0;

        do
        {
            if (*mask++)
                H[256 * *F + *G] += 1;
            
            G++;
        }
        while (++F != Fend);
    }
}

void smoothRows(float *Y, float *X, int rows, int cols, float *skrn, int skrnl)
{
    int r;
    int skrnleft = (skrnl - 1) / 2;
    
    #pragma omp for schedule(static) private(r)
    for (r = 0; r < rows; r++)
    {
        int offset = cols * r;
        
        float *Xstart = X + offset;
        float *Xend = Xstart + cols;
        float *Xi = Xstart;
        float *Yi = Y + offset;
        
        do
        {
            float *s = skrn;
	    /* TODO: Don't move pointer outside allocated range (probably won't cause problems...) */
            float *Xj = Xi - skrnleft;
            int j;
            
            *Yi = 0;
            
            for (j = 0; j < skrnl; j++)
            {
                if (Xj >= Xstart && Xj < Xend)
                    *Yi += *s * *Xj;
                
                s++;
                Xj++;
            }
            
            Yi++;
        } while (++Xi != Xend);
    }
}

void smoothCols(float *Y, float *X, int rows, int cols, float *skrn, int skrnl)
{
    int c;
    int skrnleft = (skrnl - 1) / 2;
    
    #pragma omp for schedule(static) private(c)
    for (c = 0; c < cols; c++)
    {
        float *Xstart = X + c;
        float *Xend = Xstart + rows * cols;
        float *Xi = Xstart;
        float *Yi = Y + c;
        
        do
        {
            float *s = skrn;
	    /* TODO: Don't move pointer outside allocated range (probably won't cause problems...) */
            float *Xj = Xi - skrnleft * cols;
            int j;
            
            *Yi = 0;
            
            for (j = 0; j < skrnl; j++)
            {
                if (Xj >= Xstart && Xj < Xend)
                    *Yi += *s * *Xj;
                
                s++;
                Xj += cols;
            }
            
            Yi += cols;
        } while ((Xi += cols) != Xend);
    }
}

#define EPS 1e-15
double histSumAll(float *H, int length)
{
    float *fi;
    double s = 0;
    
    for (fi = H; fi < H + length; fi++)
    {
        *fi += EPS;
        s += (double)*fi;
    }
    
    return s;
}

double histLogAll(float *H, double s, int length)
{
    float *fi;
    double sL = 0;
    
    for (fi = H; fi < H + length; fi++)
    {
        double t = (double)*fi / s;
        sL += t * log2(t);
    }
    
    return sL;
}

double histLogVert(float *H, double *s1, double s, int rows, int cols)
{
    int i, j;
    double s1L = 0;
    float *fi = H;
    double *di;

    for (di = s1; di < s1 + cols; di++)
        *di = 0;
    
    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            s1[j] += (double)*fi++;

    for (di = s1; di < s1 + cols; di++)
    {
        double t = *di / s;
        s1L += t * log2(t);
    }
    
    return s1L;
}

double histLogHorz(float *H, double *s2, double s, int rows, int cols)
{
    int i, j;
    double s2L = 0;
    float *fi = H;
    double *di;

    for (di = s2; di < s2 + rows; di++)
        *di = 0;    
    
    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            s2[i] += (double)*fi++;
   
    for (di = s2; di < s2 + rows; di++)
    {
        double t = *di / s;
        s2L += t * log2(t);
    }
    
    return s2L;
}

float regularisation(float *J, mwSize n)
{
    /* TODO: Parallelise */
    float r = 0.0;
    float *Jend = J + n;

    do
    {
        r += *J * *J++;
    } while (J != Jend);
    
    r /= n;
    
    return r;
}

double pewarpcost(float *x, float *bX, float *bY, float *dbY, float *bZ,
        mwSignedIndex *bdims, uint8_T *target, float *source, uint8_T *mask,
        mwSignedIndex *ddims, float *skrn, int skrnl, float xscale,
        float regstrength, float *H, uint8_T *warped, float *def)
{
    int nvox = ddims[0] * ddims[1] * ddims[2];
    /* float *def = mxMalloc(sizeof(float) * nvox); */
    float *jac = mxMalloc(sizeof(float) * nvox);
    /* uint8_T *warped = mxMalloc(sizeof(uint8_T) * nvox); */
    /* float *H = mxMalloc(sizeof(float) * 256 * 256); */
    float *Htmp = mxMalloc(sizeof(float) * 256 * 256);

    double *s1 = mxMalloc(sizeof(double) * 256);
    double *s2 = mxMalloc(sizeof(double) * 256);
    double s, sL, s1L, s2L;
    
    double cost;
    
    /* mexPrintf(">>> Start -> %lu (CLOCKS_PER_SEC: %lu)\n", clock(), CLOCKS_PER_SEC); */
            
    getDef(def, ddims, bdims[0], bX, bdims[1], bY, bdims[2], bZ, x, xscale);
    /* The actual Jacobian is 1 + jac, but it is easier to add 1 later. */
    getDef(jac, ddims, bdims[0], bX, bdims[1], dbY, bdims[2], bZ, x, xscale);
    
    /* mexPrintf("getDef done -> %lu\n", clock()); */
    

    /* For some reason, execution is single-threaded without the num_threads clause.
     * TODO: Find out why and change. */
    #pragma omp parallel num_threads(4) shared(s, sL, s1L, s2L)
    {
        peResampleAndApplyJacobian(warped, source, def, /* jac, */ mask, ddims);
        
        /* if (omp_get_thread_num() == 1)
            mexPrintf("peResampleAndApplyJacobian done -> %lu\n", clock()); */
        
        hist2(H, target, warped, mask, nvox);
        
        /* if (omp_get_thread_num() == 1)
            mexPrintf("hist2 done -> %lu\n", clock()); */
        
        smoothRows(Htmp, H, 256, 256, skrn, skrnl);
        smoothCols(H, Htmp, 256, 256, skrn, skrnl);
        
        /* if (omp_get_thread_num() == 1)
            mexPrintf("smooth... done -> %lu\n", clock()); */
        
        #pragma omp single
        {
            s = histSumAll(H, 65536);
        }

        #pragma omp sections
        {
            #pragma omp section
            {
                sL = histLogAll(H, s, 65536);
            }

            #pragma omp section
            {
                s1L = histLogVert(H, s1, s, 256, 256);
            }

            #pragma omp section
            {
                s2L = histLogHorz(H, s2, s, 256, 256);
            }
        }
        
        /* if (omp_get_thread_num() == 1)
            mexPrintf("hist... done -> %lu\n", clock()); */
        
    }
    
    cost = - (s1L + s2L) / sL + regstrength * regularisation(jac, nvox);
    
    /* mexPrintf("regularisation done -> %lu\n", clock());
    mexPrintf("<<<\n"); */
    
    /* mexPrintf("cost: %f  sL: %f  s1L: %f  s2L: %f  #voxels: %d\n", cost, sL, s1L, s2L, nvox); */

    /* mxFree(def); */
    mxFree(jac);
    /* mxFree(warped); */
    /* mxFree(H); */
    mxFree(Htmp);
    mxFree(s1);
    mxFree(s2);

    return cost;
}

void mexFunction(
    int nlhs,
    mxArray *plhs[],
    int nrhs,
    const mxArray *prhs[])
{
  /* Parameters:
   * x - Coefficients
   * bx - Basis vectors in x
   * by - Basis vectors in y
   * dby - Derivatives of basis vectors in y
   * bz - Basis vectors in z
   * target - Source image as uint8s
   * source - Target image as floats
   * mask - Mask
   * skrn - 1D Histogram smoothing kernel
   * xscale - Scaling factor for x
   * regstrength - Regularisation strength
   */
  float *x, *bx, *by, *dby, *bz, *source, *skrn, *xscale, *regstrength;
  uint8_T *target, *mask;
  mwSignedIndex dims[3], bdims[3];
  const mwSize *sdims, *tdims, *mdims;
  mwSize skrnl;
  
  /* Result */
  double *cost;
  float *H;
  uint8_T *warped;
  float *def;
  
  if (nrhs != 11)
      mexErrMsgTxt("Wrong number of arguments.");

  if (!mxIsSingle(prhs[0]) || mxIsComplex(prhs[0])
    || !mxIsSingle(prhs[1]) || mxIsComplex(prhs[1])
    || !mxIsSingle(prhs[2]) || mxIsComplex(prhs[2])
    || !mxIsSingle(prhs[3]) || mxIsComplex(prhs[3])
    || !mxIsSingle(prhs[4]) || mxIsComplex(prhs[4])
    || !mxIsUint8(prhs[5]) || mxIsComplex(prhs[5])
    || !mxIsSingle(prhs[6]) || mxIsComplex(prhs[6])
    || !mxIsUint8(prhs[7]) || mxIsComplex(prhs[7])
    || !mxIsSingle(prhs[8]) || mxIsComplex(prhs[8])
    || !mxIsSingle(prhs[9]) || mxIsComplex(prhs[9])
    || !mxIsSingle(prhs[10]) || mxIsComplex(prhs[10]))
  {
      mexErrMsgTxt("At least one argument is of an incorrect type.");
  }

  dims[0] = (mwSignedIndex)mxGetM(prhs[1]);
  bdims[0] = (mwSignedIndex)mxGetN(prhs[1]);
  bx = (float*)mxGetPr(prhs[1]);
  dims[1] = (mwSignedIndex)mxGetM(prhs[2]);
  bdims[1] = (mwSignedIndex)mxGetN(prhs[2]);
  by = (float*)mxGetPr(prhs[2]);
  
  if (mxGetM(prhs[3]) != dims[1] || mxGetN(prhs[3]) != bdims[1])
      mexErrMsgTxt("Matrices by and dby should have same size.");
  
  dby = (float*)mxGetPr(prhs[3]);
  dims[2] = (mwSignedIndex)mxGetM(prhs[4]);
  bdims[2] = (mwSignedIndex)mxGetN(prhs[4]);
  bz = (float*)mxGetPr(prhs[4]);

  if (mxGetM(prhs[0]) != bdims[0] * bdims[1] * bdims[2])
      mexErrMsgTxt("Wrong number of coefficients.");
  
  x = (float*)mxGetPr(prhs[0]);
  
  if (mxGetNumberOfDimensions(prhs[5]) != 3 || mxGetNumberOfDimensions(prhs[6]) != 3
          || mxGetNumberOfDimensions(prhs[7]) != 3)
      mexErrMsgTxt("All images should be 3D arrays.");
  
  tdims = mxGetDimensions(prhs[5]);
  sdims = mxGetDimensions(prhs[6]);
  mdims = mxGetDimensions(prhs[7]);
    
  if (sdims[0] != tdims[0] || sdims[1] != tdims[1] || sdims[2] != tdims[2]
          || sdims[0] != dims[0] || sdims[1] != dims[1] || sdims[2] != dims[2]
          || sdims[0] != mdims[0] || sdims[1] != mdims[1] || sdims[2] != mdims[2])
      mexErrMsgTxt("Matrix size error.");
  
  target = (uint8_T*)mxGetPr(prhs[5]);
  source = (float*)mxGetPr(prhs[6]);
  mask = (uint8_T*)mxGetPr(prhs[7]);

  if (mxGetM(prhs[8]) != 1)
      mexErrMsgTxt("Histogram smoothing kernel should be 1 x n.");
  
  skrn = (float*)mxGetPr(prhs[8]);
  skrnl = mxGetN(prhs[8]);
  
  if (mxGetNumberOfDimensions(prhs[9]) != 2 || mxGetM(prhs[9]) != 1 || mxGetN(prhs[9]) != 1)
      mexErrMsgTxt("Scaling factor should be a scalar");
  
  xscale = (float*)mxGetPr(prhs[9]);
  
  if (mxGetNumberOfDimensions(prhs[10]) != 2 || mxGetM(prhs[10]) != 1 || mxGetN(prhs[10]) != 1)
      mexErrMsgTxt("Regularisation strength factor should be a scalar");
  
  regstrength = (float*)mxGetPr(prhs[10]);

  plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
  cost = (double*)mxGetPr(plhs[0]);

  plhs[1] = mxCreateNumericMatrix(256, 256, mxSINGLE_CLASS, mxREAL);
  H = (float*)mxGetPr(plhs[1]);
  
  plhs[2] = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
  warped = (uint8_T*)mxGetPr(plhs[2]);
  
  plhs[3] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
  def = (float*)mxGetPr(plhs[3]);

  *cost = pewarpcost(x, bx, by, dby, bz, bdims, target, source, mask, dims, skrn, skrnl, *xscale, *regstrength, H, warped, def);
}
