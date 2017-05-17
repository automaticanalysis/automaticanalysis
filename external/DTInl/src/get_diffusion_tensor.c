/* Main function for nonlinear fit of diffusion tensor
_______________________________________________________________________

 Copyright (C) 2014 MRC Congition and Brain Sciences Unit

 Marta Correia and Tibor Auer
 $Id: get_diffusion_tensor.c 2014-08-15 16:00:00Z ta02 $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float get_diffusion_tensor(double *data, int size_d, double *bvals, double *bvecs, double *S0, double *tensor) {
    
    
    const int n=6, nn=3, ma=7, n_var=4;
    
    int N_eff, counter_eff, *b0_id, b0_control;
    
    double *y, *param_inic, *param, **a, **var, **x, *y_eff, *y_data_eff, **var_eff, **x_eff,  aux_var;
    
    double dir_x[size_d], dir_y[size_d], dir_z[size_d];
    
    double r[3]; /* not used?? */
    
    int i, m, l, t, p, j, k, counter, v, zero_control;
    
    /* count b0 in bvals */
    
    int N_b0 = 0;
    for (v=0; v<size_d; v++) {
        
        if(*(bvals+v)==0) N_b0++;
        
    }
    
    b0_id = (int *) malloc(N_b0 * sizeof(int));
    
    
    x = (double **) malloc( n * sizeof(double *));
    
    var = (double **) malloc(size_d * sizeof(double *));
    for(i=0; i<size_d; i++){
        
        var[i] = (double *) malloc(n_var * sizeof(double));
    }
    
    for (l=0; l<n; l++) {
        
        x[l] = (double *) malloc( size_d * sizeof(double));
    }
    
    
    /* get values for x_dir, y_dir, z_dir */
    
    for(i=0; i<size_d; i++){
        
        dir_x[i]=bvecs[i*3];
        dir_y[i]=bvecs[(i*3)+1];
        dir_z[i]=bvecs[(i*3)+2];
        
    }
    
    
    /* set values for var and x */
    
    for(i=0; i<size_d; i++){
        
        x[0][i]=-bvals[i]*dir_x[i]*dir_x[i];
        x[1][i]=-bvals[i]*dir_y[i]*dir_y[i];
        x[2][i]=-bvals[i]*dir_z[i]*dir_z[i];
        x[3][i]=-2*bvals[i]*dir_x[i]*dir_y[i];
        x[4][i]=-2*bvals[i]*dir_x[i]*dir_z[i];
        x[5][i]=-2*bvals[i]*dir_y[i]*dir_z[i];
        
        var[i][0]=bvals[i];
        var[i][1]=dir_x[i];
        var[i][2]=dir_y[i];
        var[i][3]=dir_z[i];
        
    }
    
    
    y = (double *) malloc( size_d * sizeof(double));
    
    param = (double *) malloc(ma * sizeof(double));
    param_inic = (double *) malloc(ma * sizeof(double));
    
    
    a = (double **) malloc( 2 * sizeof(double *));
    
    for(l=0; l<2; l++){
        
        a[l] = (double *) malloc((n+1)*sizeof(double));
    }
    
    
    /* clear b0_id[] */
    
    for(i=0; i<N_b0; i++){
        
        b0_id[i]=0;
        
    }
    
    
    /* populate b0_id with the volume ids for the b=0 volumes */
    
    j=0;
    for(v=0; v<size_d; v++){
        
        if(var[v][0]==0){
            
            b0_id[j]=v;
            j++;
        }
    }
    
    
    /* check whether there are any zero valued data points */
    
    zero_control=0;
    
    for(v=0; v<size_d; v++){
        
        if(data[v]==0) zero_control=1;
        
    }
    
    
    
    if(zero_control==0){
        
        /* if no zeros apply logs and apply linear fit to find initial parameters */
        
        
        for(v=0; v<size_d; v++){
            
            y[v]=log(data[v]);
            
        }
        
        
        fit_linear(y, x, n, size_d, a);
        
        
        for(m=0; m<ma-1; m++){
            
            param_inic[m+1]=a[0][m];
        }
        
        param_inic[0]=exp(a[0][6]);
        
        
        exponencial_fit(var, size_d, ma, n_var, param_inic, data, param);
        
    }
    
    
    else{
        
        
        /* check whether the zero valued data points corresponds to a b=0 volume */
        
        b0_control=0;
        
        for(j=0; j<N_b0; j++){
            
            if(data[b0_id[j]]==0){
                
                b0_control=1;
            }
        }
        
        
        if(b0_control==1){
            
            /* if b0_control==1 fit not performed and return zero for all parameters */
            
            for(j=0; j<7; j++){
                
                param[j]=0;
            }
            
        }
        
        else{
            
            /* if the zero values are elsewhere, exclude those from the 'inital guess' linear fit */
            
            N_eff=0;
            
            for(v=0; v<size_d; v++){
                
                if(data[v] > 0){
                    
                    N_eff++;
                }
            }
            
            y_eff = (double *) malloc( N_eff * sizeof(double));
            
            y_data_eff = (double *) malloc( N_eff * sizeof(double));
            
            
            var_eff = (double **) malloc(N_eff * sizeof(double *));
            for(i=0; i<N_eff; i++){
                
                var_eff[i] = (double *) malloc(n_var * sizeof(double));
            }
            
            x_eff = (double **) malloc( n * sizeof(double *));
            for (l=0; l<n; l++) {
                
                x_eff[l] = (double *) malloc( N_eff * sizeof(double));
            }
            
            
            counter_eff=0;
            
            for(v=0; v<size_d; v++){
                
                if(data[v] > 0){
                    
                    x_eff[0][counter_eff]=x[0][v];
                    x_eff[1][counter_eff]=x[1][v];
                    x_eff[2][counter_eff]=x[2][v];
                    x_eff[3][counter_eff]=x[3][v];
                    x_eff[4][counter_eff]=x[4][v];
                    x_eff[5][counter_eff]=x[5][v];
                    
                    var_eff[counter_eff][0]=var[v][0];
                    var_eff[counter_eff][1]=var[v][1];
                    var_eff[counter_eff][2]=var[v][2];
                    var_eff[counter_eff][3]=var[v][3];
                    
                    y_data_eff[counter_eff]=data[v];
                    
                    y_eff[counter_eff]=log(y_data_eff[counter_eff]);
                    
                    counter_eff++;
                    
                }
            }
            
            aux_var=0;
            
            for(v=0; v<N_eff; v++){
                
                aux_var += y_eff[v];
            }
            
            
            
            if(aux_var==0){
                
                /* if all data is zero, also returns zero for all parameters and exit */
                
                for(j=0; j<7; j++){
                    
                    param[j]=0;
                }
            }
            
            
            else{
                
                fit_linear(y_eff, x_eff, n, N_eff, a);
                
                for(m=0; m<ma-1; m++){
                    
                    param_inic[m+1]=a[0][m];
                    
                }
                
                param_inic[0]=exp(a[0][6]);
                
                
                /* non-linear tensor fit */
                
                exponencial_fit(var, size_d, ma, n_var, param_inic, data, param);
                
                
            }
            
            
            free(y_eff);
            free(y_data_eff);
            
            for(i=0; i<N_eff; i++){
                
                free(var_eff[i]);
            }
            free(var_eff);
            
            for (l=0; l<n; l++) {
                
                free(x_eff[l]);
            }
            free(x_eff);
            
        }
        
    }
    
    /* create results */
    *S0 = param[0];
    
    const int tCol = 3, tRow = 3, index[] = {1,4,5,4,2,6,5,6,3};
    
    l=0;
    for(i=0;i<tRow;i++)
        for(j=0;j<tCol;j++)
            *(tensor+(i*tCol)+j) = param[index[l++]];

    return 0;
    
}
