/* Various functions used by dti_nlfit
_______________________________________________________________________

 Copyright (C) 2014 MRC Congition and Brain Sciences Unit

 Marta Correia
 $Id: Fit_routines.c 2014-08-15 16:00:00Z ta02 $ */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include "values.h"
#include "limits.h"


void invert( double **a, int n, double **out);

void svdcmp(double **a, int m, int n, double *w, double **v);

double pythag(double a, double b);

void fit_linear(double *y, double **x, int n, int N, double **out);

void ROTATE(double **a, double s, double tau, int i, int j, int k, int l);

void jacobi(double **a, int n, double *d, double **v);

void funcs(double *x, double *a, double *ymod, double *dyda, int ma);

void funcs_Behrens(double *x, double *a, double *ymod, double *dyda, int ma);

void mrqcof(double **x, double y[], double sig[], int ndata, int n_var, double *a, int ia[], int ma, double **alpha, double beta[], double *chisq, void (*funcs)(double *, double [], double *, double [], int));

void gaussj(double **a, int n, double **b, int m);

void covsrt(double **covar, int ma, int ia[], int mfit);

void mrqmin(double **x, double y[], double sig[], int ndata, int n_var, double *a, int ia[], int ma, double **covar, double **alpha, double *chisq, void (*funcs)(double *, double [], double *, double [], int), double *alamda);

void exponencial_fit( double **x, int ndata, int ma, int n_var, double *a_inic, double *y_data, double *out);

void Behrens_fit( double **x, int ndata, int ma, int n_var, double *a_inic, double *y_data, double *out);

void tanh_fit( double **x, int ndata, int ma, int n_var, double *a_inic, double *y_data, double *out);

void funcs_tanh(double *x, double *a, double *ymod, double *dyda, int ma);

void swapshort(unsigned char* c);




/* Inverts a matrix a and returns the result in out */

void svdcmp(); 

void invert( double **a, int n, double **out) { 

	int i, j, k; 

	double **bmat, **v, *w; 
	
	double tmp; 

	bmat = (double **) malloc( n * sizeof(double *));  
	v = (double **) malloc( n * sizeof(double *)); 
	w = (double *) malloc( n * sizeof(double)); 
	
	for (i=0; i<n; i++) {
		bmat[i] = (double *) malloc( n * sizeof(double)); 
		v[i] = (double *) malloc( n * sizeof(double)); 
      
		for (j=0; j<n; j++) { 
			bmat[i][j] = a[i][j];
		}
	}
	
	/* printf("INPUT:\n"); 
	for (i=0; i<mdl.n; i++) { 
		for (j=0; j<7; j++) { 
			printf("%g ", bmat[i][j]); 
		}
		printf("\n"); 
	} */ 

	svdcmp(bmat, n, n, w,  v); 

/*
	printf("OUTPUT U:\n"); 
	for (i=0; i<mdl.n; i++) { 
		for (j=0; j<7; j++) { 
			printf("%g ", bmat[i][j]); 
		}
		printf("\n"); 
	}
*//*	
	printf("OUTPUT W:\n"); 
	for (i=0; i<n; i++) { printf("%g ", w[i]); } 
	printf("\n");
*/	/*
	printf("OUTPUT V:\n"); 
	for (i=0; i<7; i++) { 
		for (j=0; j<7; j++) { 
			printf("%g ", v[i][j]); 
		}
		printf("\n"); 
	}
	
	printf("PRODUCT:\n"); 
	for (i=0; i<mdl.n; i++) {
		for (j=0; j<7; j++) {
			tmp = 0.0;
			for (k=0; k<7; k++) { tmp += bmat[i][k] * w[k] * v[j][k]; }
			printf("%g ", tmp); 
		}
		printf("\n"); 
	}
*/
	for (i=0; i<n; i++) { 
		for (j=0; j<n; j++) { 
			out[i][j] = 0.0; 
			for (k=0; k<n; k++) { 
				if ((fabs(w[k]) > FLT_MIN ) && ( fabs(w[k]) > 1e-06 * fabs(w[0])) ) { 
					out[i][j] += v[i][k] * bmat[j][k] / w[k]; 
				} 
			}
		}
	}
	
	for (i=0; i<n; i++) { free(bmat[i]); }
	for (i=0; i<n; i++) { free(v[i]); }
	free (bmat); 
	free (v); 
	free(w); 
}

#define NRANSI
#define SIGN(a,b) (b>=0?(a>=0?a:-a):(a>=0?-a:a))
#define MIN(a,b) (b<a?b:a)
#define MAX(a,b) (b>a?b:a)


void svdcmp(double **a, int m, int n, double *w, double **v)
{
	double pythag(double a, double b);
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1 = (double *) malloc( n * sizeof(double) ); /* rv1=vector(1,n); */
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a[k][i]);
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i != n) {
			for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l-1]=f-g;
				for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=1;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l-1;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) { fprintf(stderr, "no convergence in 30 svdcmp iterations"); exit(1); }
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free(rv1); /* free_vector(rv1,1,n); */
}
#undef NRANSI

/* (C) Copr. 1986-92 Numerical Recipes Software #!-"k$'%)]'V',4$. */


#define NRANSI

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+pow(absb/absa,2.0));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+pow(absa/absb,2.0)));
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software #!-"k$'%)]'V',4$. */

/* ***************************************************************** */
/* ***************************************************************** */


/* Fits the simulated data to the single tensor model */

void fit_linear(double *y, double **x, int n, int N, double **out){

  int i, j, k;
  double sy, sy2, av_y, s2, aux, aux2;
  double *sj, *sj2, *rjy, *sjy2, *av_x, *sum, *temp;
  double  **rjk, **inv_rjk,  **sjk2;

  
  sj=(double *) malloc(n * sizeof(double));
  sj2=(double *) malloc(n * sizeof(double));
  rjy=(double *) malloc(n * sizeof(double));
  sjy2=(double *) malloc(n * sizeof(double));
  av_x=(double *) malloc(n * sizeof(double));
  sum=(double *) malloc(N * sizeof(double));
  temp=(double *) malloc(n * sizeof(double));

  rjk = (double **) malloc( n * sizeof(double *)); 
  inv_rjk = (double **) malloc( n * sizeof(double *)); 
  sjk2 = (double **) malloc( n * sizeof(double *));

  for (i=0; i<n; i++) {
 
   rjk[i] = (double *) malloc( n * sizeof(double)); 
   inv_rjk[i] = (double *) malloc( n * sizeof(double)); 
   sjk2[i] = (double *) malloc( n * sizeof(double));  
  
  }


  /**** Calculates the average of y ****/

  av_y=0.0;

  for(i=0; i<N; i++){

    av_y += y[i]/N;
   
  }

 
  /**** Calculates the averaje of each xj ****/

  for(j=0; j<n; j++){
    
    av_x[j]=0.0;
  }
  
  for(j=0; j<n; j++){
    for(i=0; i<N; i++){
      
      av_x[j] += x[j][i]/N;
   }
  }



  /**** Calculates sy ****/

  sy=0;
  sy2=0;
  for(i=0; i<N; i++){
    
    sy2 += (y[i]-av_y)*(y[i]-av_y)/(N-1);
  }
  sy=sqrt(sy2);



  /**** Calculates each sj ****/

  for(j=0; j<n; j++){
    
    sj[j]=0;
    sj2[j]=0;
  }

  for(j=0; j<n; j++){
    for(i=0; i<N; i++){
      
      sj2[j] += (x[j][i]-av_x[j])*(x[j][i]-av_x[j])/(N-1);
    }
  }
  
  for(j=0; j<n; j++){

    sj[j]=sqrt(sj2[j]);

  }
  


  /**** Calculates each sjy2 ****/
  
  for(j=0; j<n; j++){
    
    sjy2[j]=0;
  }
  
  for(j=0; j<n; j++){
    for(i=0; i<N; i++){

      sjy2[j] += (x[j][i]-av_x[j])*(y[i]-av_y)/(N-1);
    }
  }
  

  
  /**** Calculates each sjk2 ****/
 
  for(j=0; j<n; j++){
    for(k=0; k<n; k++){
     
      sjk2[j][k]=0;
    }
  }
  
 
  for(j=0; j<n; j++){
    for(k=0; k<n; k++){
      for(i=0; i<N; i++){

	sjk2[j][k] += (x[j][i]-av_x[j])*(x[k][i]-av_x[k])/(N-1);
      }
    }
  }
  
 

  /**** Calculates each rjy ****/

  for(j=0; j<n; j++){
   
    rjy[j]=sjy2[j]/(sj[j]*sy);

  }
  
 
 
 /**** Calculates each rjk ****/

  for(j=0; j<n; j++){
    for(k=0; k<n; k++){
      
      rjk[j][k]=sjk2[j][k]/(sj[j]*sj[k]);

    }
  }




  /**** Invert matrix [rjk] ****/

  invert(rjk, n, inv_rjk);
 
  /*
  //for(j=0; j<n; j++){
  //  for(k=0; k<n; k++){
      
    //  fprintf(output, "%lf\t", rjk[j][k]);
   // }
   // fprintf(output, "\n");
 // }

//  for(j=0; j<n; j++){
  //  for(k=0; k<n; k++){
      
  //    printf("inv_rjk=%lf\n",inv_rjk[j][k]);
  //  }
    
//  }
*/

  /**** Calculates the regression coeficients ****/

  for(j=0; j<n; j++){
    
    out[0][j]=0;
  }

  for(j=0; j<n; j++){
    for(k=0; k<n; k++){

      out[0][j] += (sy/sj[j])*(rjy[k]*inv_rjk[j][k]);
    }
  }

  out[0][n]=av_y;

  for(j=0; j<n; j++){
    
    out[0][n] += -out[0][j]*av_x[j];
  }

/*
//  for(j=0; j<n; j++){

//    printf("out[0][%i]=%lf\n", j, out[0][j]);
//  }
*/
  
  /**** Calculates the uncertainties associated with the regression coefficients ***/

  for(i=0; i<N; i++){
   
    sum[i]=0.0;
  }

  for(i=0; i<N; i++){
    for(j=0; j<n; j++){
     
      sum[i] += out[0][j]*(x[j][i]-av_x[j]);
    }
  }

  s2=0;

  for(i=0; i<N; i++){

    s2 += ((y[i]-av_y)-sum[i])*((y[i]-av_y)-sum[i])/(N-n-1);
  }
 

  for(j=0; j<n; j++){

    aux=(s2/sj2[j])*inv_rjk[j][j]/(N-1);
    out[1][j]=sqrt(aux);
  }

 /*
  //for(k=0; k<n; k++){
   // printf("inv_rjk=%lf\n", inv_rjk[k][k]);
 // }
 */

  aux=0;
  out[1][n]=0;
  
  for(j=0; j<n; j++){
 
   temp[j]=0;
  }

  for(j=0; j<n; j++){
    for(k=0; k<n; k++){

      temp[j] += (av_x[j]*av_x[k]*inv_rjk[j][k])/(sj[j]*sj[k]);
    }
  }

  for(j=0; j<n; j++){

    aux += (av_x[j]*av_x[j]*inv_rjk[j][j])/(sj2[j])+temp[j];
  }

  aux2=((1/N)+aux/(N-1))*s2;

  out[1][n]=sqrt(aux2);



  /**** Free memory ****/
 
  free(sj);
  free(sj2);
  free(rjy);
  free(sjy2);
  free(av_x);
  free(sum);
  free(temp);
 
  for (i=0; i<n; i++) { free(rjk[i]); }
  for (i=0; i<n; i++) { free(sjk2[i]); }
  for (i=0; i<n; i++) { free(inv_rjk[i]); } 
  free (rjk);
  free (inv_rjk); 
  free (sjk2); 

}


/* ****************************************************************** */
/* ****************************************************************** */


/* Returns the eighen-values of matrix a in d */
/* Returns the eighrn-vectors of a in v */

void ROTATE(double **a, double s, double tau, int i, int j, int k, int l){
  
  double g, h;  

  g=a[i][j];
  h=a[k][l];
  a[i][j]=g-s*(h+g*tau);
  a[k][l]=h+s*(g-h*tau);

}

void jacobi(double **a, int n, double *d, double **v){

  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  b = (double *) malloc( n * sizeof(double));
  z = (double *) malloc( n * sizeof(double)); 
	
  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }

  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
       
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      free(z);
      free(b);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
	    && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((double)(fabs(h)+g) == (double)fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=0;j<ip;j++) {
	    ROTATE(a,s,tau,j,ip,j,iq);
	  }
	  for (j=ip+1;j<iq;j++) {
	    ROTATE(a,s,tau,ip,j,j,iq);
	  }
	  for (j=iq+1;j<n;j++) {
	    ROTATE(a,s,tau,ip,j,iq,j);
	  }
	  for (j=0;j<n;j++) {
	    ROTATE(v,s,tau,j,ip,j,iq);
	  }
	  
	}
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  printf("\n\nERROR: Too many iterations in routine jacobi\n\n");
  exit(1);
}

/* (C) Copr. 1986-92 Numerical Recipes Software #!-"k$'%)]'V',4$. */

/* ****************************************************************** */
/* ****************************************************************** */

#define NRANSI


#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}





void funcs(double *x, double *a, double *ymod, double *dyda, int ma){

  double arg, ex;

  arg=a[1]*x[1]*x[1]+a[2]*x[2]*x[2]+a[3]*x[3]*x[3]+2*a[4]*x[1]*x[2]+2*a[5]*x[1]*x[3]+2*a[6]*x[2]*x[3];
  ex=exp(-x[0]*arg);

  *ymod=a[0]*ex;

  dyda[0]=ex;
  dyda[1]=(*ymod)*(-x[0]*x[1]*x[1]);
  dyda[2]=(*ymod)*(-x[0]*x[2]*x[2]);
  dyda[3]=(*ymod)*(-x[0]*x[3]*x[3]);
  dyda[4]=(*ymod)*(-x[0]*2*x[1]*x[2]);
  dyda[5]=(*ymod)*(-x[0]*2*x[1]*x[3]);
  dyda[6]=(*ymod)*(-x[0]*2*x[2]*x[3]);

}

void funcs_Behrens(double *x, double *a, double *ymod, double *dyda, int ma){

  double prod, f, d, theta, phi, S0, b, rx, ry, rz;

  f=a[0];
  d=a[1];
  theta=a[2];
  phi=a[3];
  S0=a[4];

  b=x[0];
  rx=x[1];
  ry=x[2];
  rz=x[3];

  prod = sin(theta)*cos(phi)*rx + sin(theta)*sin(phi)*ry + cos(theta)*rz;


  *ymod=S0 * ( (1-f)*exp(-b*d) + f*exp(-b*d*prod*prod) );


  dyda[0]= -S0 * exp(-b*d) + S0 * exp(-b*d*prod*prod);

  dyda[1]=-b*S0*( (1-f)*exp(-b*d) + f*prod*prod*exp(-b*d*prod*prod) );

  dyda[2]=-2*b*d*prod*S0*f*exp(-b*d*prod*prod)*(cos(theta)*cos(phi)*rx + cos(theta)*sin(phi)*ry - sin(theta)*rz);

  dyda[3]=-2*b*d*prod*S0*f*exp(-b*d*prod*prod)*(-sin(theta)*sin(phi)*rx + sin(theta)*cos(phi)*ry);
  
  dyda[4]=(1-f)*exp(-b*d) + f * exp(-b*d*prod*prod);

}


void funcs_tanh(double *x, double *a, double *ymod, double *dyda, int ma){

  double cos_h;

  (*ymod)=(exp(x[0]/a[0])-exp(-x[0]/a[0]))/(exp(x[0]/a[0])+exp(-x[0]/a[0]));

/*  printf("a=%lf x=%lf\n", a[0], x[0]); */

/*  printf("tanh=%lf\n", (*ymod)); */


  cos_h=(exp(x[0]/a[0])+exp(-x[0]/a[0]))/2;


  dyda[0]=1/(a[0]*pow(cos_h, 2));


}





void mrqcof(double **x, double y[], double sig[], int ndata, int n_var, double *a, int ia[],
	int ma, double **alpha, double beta[], double *chisq,
	void (*funcs)(double *, double [], double *, double [], int))
{
  int i,j,k,l,m,mfit=0, aux;
	double *ymod,wt,sig2i,dy,*dyda, *arg;

	dyda = (double *) malloc (ma * sizeof(double)); /*dyda=vector(1,ma);*/
	ymod = (double *) malloc (1 * sizeof(double));
	arg = (double *) malloc (n_var * sizeof(double));

	for (j=0;j<ma;j++)
		if (ia[j]) mfit++;
	for (j=0;j<mfit;j++) {
		for (k=0;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=0;i<ndata;i++) {
	        for(aux=0; aux<n_var; aux++){
	               arg[aux]=x[i][aux];
	        }

		(*funcs)(arg,a,ymod,dyda,ma);


		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-(*ymod);
		for (j=0,l=0;l<ma;l++) {
			if (ia[l]) {
				wt=dyda[l]*sig2i;
				for (k=0,m=0;m<l+1;m++)
					if (ia[m]) alpha[j][k++] += wt*dyda[m];
				beta[j++] += dy*wt;
			}
		}
		*chisq += dy*dy*sig2i;
	}
	for (j=1;j<mfit;j++)
		for (k=0;k<j;k++) alpha[k][j]=alpha[j][k];

	
	free(dyda);
	free(ymod);
	free(arg);
}

/* (C) Copr. 1986-92 Numerical Recipes Software #!-"k$'%)]'V',4$. */




void gaussj(double **a, int n, double **b, int m)
{



	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;



	indxc = (int *) malloc(n * sizeof(int));        /*   indxc=ivector(1,n);*/
	indxr = (int *) malloc(n * sizeof(int));        /*   indxr=ivector(1,n);*/
	ipiv = (int *) malloc(n * sizeof(int));        /*   ipiv=ivector(1,n);*/


	for (j=0;j<n;j++) ipiv[j]=0;



	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} 
				}



		++(ipiv[icol]);



		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0){

		  printf("\n\ngaussj: Singular Matrix\n");
		  exit(1);
		}




		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}


	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}



	free(ipiv);
	free(indxr);
	free(indxc);
}

/* (C) Copr. 1986-92 Numerical Recipes Software #!-"k$'%)]'V',4$. */





void covsrt(double **covar, int ma, int ia[], int mfit)
{
	int i,j,k;
	double temp;

	for (i=mfit;i<ma;i++)
		for (j=0;j<i+1;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit-1;
	for (j=ma-1;j>=0;j--) {
		if (ia[j]) {
			for (i=0;i<ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=0;i<ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		}
	}
}

/* (C) Copr. 1986-92 Numerical Recipes Software #!-"k$'%)]'V',4$. */





void mrqmin(double **x, double y[], double sig[], int ndata, int n_var, double *a, int ia[],
	int ma, double **covar, double **alpha, double *chisq,
	void (*funcs)(double *, double [], double *, double [], int), double *alamda)
{
  void covsrt(double **covar, int ma, int ia[], int mfit);
  void gaussj(double **a, int n, double **b, int m);
  void mrqcof(double **x, double y[], double sig[], int ndata, int n_var, double a[],
	      int ia[], int ma, double **alpha, double beta[], double *chisq,
	      void (*funcs)(double *, double [], double *, double [], int));
  


  
  int j,k,l,m, aux;
  static int mfit;
  static double ochisq,*atry,*beta,*da,**oneda, **temp;
  
  if (*alamda < 0.0) {

    atry = (double *) malloc(ma * sizeof(double));   /* atry=vector(1,ma);*/
    beta = (double *) malloc(ma * sizeof(double));   /* beta=vector(1,ma);*/
    da = (double *) malloc(ma * sizeof(double));   /* da=vector(1,ma);*/
    
    
    mfit=0;
    
    for (j=0;j<ma;j++)
      if (ia[j]) mfit++;


    oneda = (double **) malloc(mfit * sizeof(double*)); /* oneda=matrix(1,mfit,1,1);*/
    for(j=0; j<mfit; j++){
      oneda[j] = (double *) malloc( 1 * sizeof(double));
    }
    
    *alamda=0.001;
    

    
    mrqcof(x,y,sig,ndata,n_var,a,ia,ma,alpha,beta,chisq,funcs);
    
    

    
    
    ochisq=(*chisq);
	  
    for (j=0;j<ma;j++) atry[j]=a[j];
  }
  

  
  temp = (double **) malloc(mfit * sizeof(double *));
  for(aux=0; aux<mfit; aux++){
    temp[aux] = (double *) malloc(mfit * sizeof(double));
  }
  

  
  
  for(j=0; j<mfit; j++){
	  
    for(k=0; k<mfit; k++){
      covar[j][k]=alpha[j][k];
    }
    
    
    
    covar[j][j]=alpha[j][j]*(1.0+(*alamda));
    
    for(k=0; k<mfit; k++){
      temp[j][k]=covar[j][k];
    }
    
    
    oneda[j][0]=beta[j];
    
    
  }
  

  
  
  gaussj(temp, mfit, oneda, 1);
  
  

  
  for(j=0; j<mfit; j++){
    
    for(k=0; k<mfit; k++){
      covar[j][k]=temp[j][k];
    }
    
    da[j]=oneda[j][0];
    
  }
  
  


  if (*alamda == 0.0) {
    covsrt(covar,ma,ia,mfit);
    
    for(j=0; j<mfit; j++){ free(oneda[j]); }  /* free_matrix(oneda,1,mfit,1,1);*/
    free(oneda);
    
    for(j=0; j<mfit; j++){ free(temp[j]); }
    free(temp);
    
    free(da);
    free(beta);
    free(atry);
    
    return;
  }
  
  
  
  for (j=0,l=0;l<ma;l++)
    if (ia[l]) atry[l]=a[l]+da[j++];
  
  
  
  mrqcof(x,y,sig,ndata,n_var,atry,ia,ma,covar,da,chisq,funcs);



  if (*chisq < ochisq) {
    *alamda *= 0.1;
    ochisq=(*chisq);
    
    for(j=0; j<mfit; j++){
      for(k=0; k<mfit; k++){
	alpha[j][k]=covar[j][k];
      }
      beta[j]=da[j];
    }
    for(l=0; l<ma; l++){
	    a[l]=atry[l];
    }
  }
  else{
    (*alamda) *= 10.0;
    (*chisq)=ochisq;
  }
  
  
  for(j=0; j<mfit; j++){ free(oneda[j]); }  /* free_matrix(oneda,1,mfit,1,1);*/
  free(oneda);
  
  for(j=0; j<mfit; j++){ free(temp[j]); }
  free(temp);
  
  free(da);
  free(beta);
  free(atry);
  
}

#undef SWAP
#undef NRANSI

/* (C) Copr. 1986-92 Numerical Recipes Software #!-"k$'%)]'V',4$. */



void exponencial_fit( double **x, int ndata, int ma, int n_var, double *a_inic, double *y_data, double *out){




  int  counter, i, j;

  int ia[ma];

  double **covar, **alpha, *chisq, *alambda, chi_aux, chi_var;
  double y[ndata], sig[ndata];


  covar = (double **) malloc(ma * sizeof(double *));
  alpha = (double **) malloc(ma * sizeof(double *));


  for(i=0; i<ma; i++){

    covar[i] = (double *) malloc(ma * sizeof(double));
    alpha[i] = (double *) malloc(ma * sizeof(double));
  }




  chisq = (double *) malloc(1 * sizeof(double));
  alambda = (double *) malloc(1 * sizeof(double));


  for(i=0; i<ma; i++){

    out[i]=a_inic[i];

    ia[i]=1;
  }

  *alambda = -0.01;


  for(j=0; j<ndata; j++){

    y[j]=y_data[j];

    sig[j]=1;

  }

  counter=1;  
  chi_var=1;



  mrqmin(x, y, sig, ndata, n_var, out, ia, ma, covar, alpha, chisq, funcs, alambda);



  chi_aux=(*chisq);



  while(chi_var > 0.01){

    *alambda = -0.01;
    
    mrqmin(x, y, sig, ndata, n_var, out, ia, ma, covar, alpha, chisq, funcs, alambda);
    


    chi_var=fabs( ((*chisq)-chi_aux)/chi_aux );
      
    chi_aux=(*chisq);
    
  /*  printf("counter=%i\n", counter);*/
    counter++;
  }


/*  printf("\nchisq=%.12lf\n\n",(*chisq));*/

  free(chisq);
  free(alambda);


  for(i=0; i<ma; i++){

    free(covar[i]);
    free(alpha[i]);
  }

  free(covar);
  free(alpha);

}



void Behrens_fit( double **x, int ndata, int ma, int n_var, double *a_inic, double *y_data, double *out){




  int  counter, i, j;

  int ia[ma];

  double **covar, **alpha, *chisq, *alambda, chi_aux, chi_var;
  double y[ndata], sig[ndata];


  covar = (double **) malloc(ma * sizeof(double *));
  alpha = (double **) malloc(ma * sizeof(double *));


  for(i=0; i<ma; i++){

    covar[i] = (double *) malloc(ma * sizeof(double));
    alpha[i] = (double *) malloc(ma * sizeof(double));
  }




  chisq = (double *) malloc(1 * sizeof(double));
  alambda = (double *) malloc(1 * sizeof(double));


  for(i=0; i<ma; i++){

    out[i]=a_inic[i];

    ia[i]=1;
  }

  *alambda = -0.01;


  for(j=0; j<ndata; j++){

    y[j]=y_data[j];

    sig[j]=1;

  }

  counter=1;  
  chi_var=1;



  mrqmin(x, y, sig, ndata, n_var, out, ia, ma, covar, alpha, chisq, funcs_Behrens, alambda);



  chi_aux=(*chisq);



  while(chi_var > 0.01){

    *alambda = -0.01;
    
    mrqmin(x, y, sig, ndata, n_var, out, ia, ma, covar, alpha, chisq, funcs_Behrens, alambda);
    


    chi_var=fabs( ((*chisq)-chi_aux)/chi_aux );
      
    chi_aux=(*chisq);
    
  /*  printf("counter=%i\n", counter);*/
    counter++;
  }


/*  printf("\nchisq=%.12lf\n\n",(*chisq));*/

  free(chisq);
  free(alambda);


  for(i=0; i<ma; i++){

    free(covar[i]);
    free(alpha[i]);
  }

  free(covar);
  free(alpha);

}

/*
//void swapshort(unsigned char* c)
//{
//  unsigned char b1,b2;
//  b2 = *c;
//  b1 = *(c+1);
//  *(c) = b1;
//  *(c+1)=b2;
//}
*/


void tanh_fit( double **x, int ndata, int ma, int n_var, double *a_inic, double *y_data, double *out){




  int  counter, i, j;

  int ia[ma];

  double **covar, **alpha, *chisq, *alambda, chi_aux, chi_var;
  double y[ndata], sig[ndata];


  covar = (double **) malloc(ma * sizeof(double *));
  alpha = (double **) malloc(ma * sizeof(double *));


  for(i=0; i<ma; i++){

    covar[i] = (double *) malloc(ma * sizeof(double));
    alpha[i] = (double *) malloc(ma * sizeof(double));
  }




  chisq = (double *) malloc(1 * sizeof(double));
  alambda = (double *) malloc(1 * sizeof(double));


  for(i=0; i<ma; i++){

    out[i]=a_inic[i];

    ia[i]=1;
  }

  *alambda = -0.01;


  for(j=0; j<ndata; j++){

    y[j]=y_data[j];

    sig[j]=1;

  }

  counter=1;  
  chi_var=1;



  mrqmin(x, y, sig, ndata, n_var, out, ia, ma, covar, alpha, chisq, funcs_tanh, alambda);



  chi_aux=(*chisq);

    printf("\nchisq=%.12lf\n\n",(*chisq));

  while(chi_var > 0.01){

    *alambda = -0.01;
    
    mrqmin(x, y, sig, ndata, n_var, out, ia, ma, covar, alpha, chisq, funcs_tanh, alambda);
    


    chi_var=fabs( ((*chisq)-chi_aux)/chi_aux );
      
    chi_aux=(*chisq);
    
    printf("counter=%i\n", counter);
    counter++;
  }


    printf("\nchisq=%.12lf\n\n",(*chisq));

  free(chisq);
  free(alambda);


  for(i=0; i<ma; i++){

    free(covar[i]);
    free(alpha[i]);
  }

  free(covar);
  free(alpha);

}





