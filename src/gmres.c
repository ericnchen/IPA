/**
 * gmres.c
 */

#include <stdlib.h>
#include <math.h>
#include "gmres.h"

#include "stdio.h"

double gmres(double *val, int *col_ind, int *row_ptr, int rows, int nnz,
             double *bvec, double *xvec) {

    int m = 10; // restart parameter

    // r0 = b - Ax0
    double *r   = malloc(rows*sizeof(double));
    double *ax0 = malloc(rows*sizeof(double));
    mvmul(ax0, val, col_ind, row_ptr, rows, nnz, xvec);
    for (int i = 0; i < rows; i++) {
        r[i] = bvec[i]-ax0[i];
    }

    // v1 = r0/norm(r0)
    double **v = malloc(rows*sizeof(double));
    for (int i = 0; i < rows; i++) {
        v[i] = malloc((rows+1)*sizeof(double));
    }
    double norm = 0.0;
    for (int i = 0; i < rows; i++) {
        norm += r[i]*r[i];
    }
    norm = sqrt(norm);
    for (int i = 0; i < rows; i++) {
        v[i][0] = r[i]/norm;
    }

    // g = norm(r)*e1
    double *g = malloc((m+1)*sizeof(double));
    g[0] = norm;
    for (int i = 1; i < m+1; i++) {
        g[i] = 0;
    }

    // outer loop allocation
    double **h = malloc((m+1)*sizeof(double));
    double **rm = malloc(m*sizeof(double));
    for (int i = 0; i <= m; i++) {
        h[i] = malloc((m+1)*sizeof(double));
    }
    for (int i = 0; i < m; i++) {
        rm[i] = malloc((m+1)*sizeof(double));
    }
    double tempout;

    // inner loop allocation
    double *c = malloc((m+1)*sizeof(double));
    double *s = malloc((m+1)*sizeof(double));
    double tempin;

    // outer loop
    for (int j = 0; j < m; j++) {
        getkrylov(val, col_ind, row_ptr, rows, nnz, j, v, h);
        rm[0][j] = h[0][j];

        // inner loop
        for (int k = 1; k <= j; k++) {
            tempin = c[k-1]*rm[k-1][j]+s[k-1]*h[k][j];
            rm[k][j] = -s[k-1]*rm[k-1][j]+c[k-1]*h[k][j];
            rm[k-1][j] = tempin;
        }

        // steps 8 and 9
        tempout = sqrt(rm[j][j]*rm[j][j]+h[j+1][j]*h[j+1][j]);
        c[j] = rm[j][j]/tempout;
        s[j] = h[j+1][j]/tempout;
        rm[j][j] = c[j]*rm[j][j]+s[j]*h[j+1][j];
        g[j+1] = -s[j]*g[j];
        g[j] = c[j]*g[j];

    }

    // residual
    double rho = fabs(g[m]);
    for (int i = 0; i <= m; i++) {
        printf("%f\n", fabs(g[i]));
    }

    // xm = x0 + Vm(Rm^-1gm)
    double *irg = malloc(m*sizeof(double));
    for (int i = 9-1; i >= 0; i--) {
        irg[i] = g[i];
        for (int k = i+1; k < 9; k++) {
            irg[i] = irg[i] - rm[i][k]*irg[k];
        }
        irg[i] = irg[i]/rm[i][i];
    }

    // new x vector
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            xvec[i] += v[i][j]*irg[j];
        }
    }

    // deallocate
    free(r);
    free(ax0);
    for (int i = 0; i < rows; i++) {
        free(v[i]);
    }
    free(v);
    free(g);
    for (int i = 0; i <= m; i++) {
        free(h[i]);
    }
    free(h);
    for (int i = 0; i < m; i++) {
        free(rm[i]);
    }
    free(rm);
    free(c);
    free(s);
    free(irg);

    printf("%f\n", rho);
    return rho;
}

void getkrylov(double *val, int *col_ind, int *row_ptr, int rows, int nnz,
               int j, double **v, double **h) {

    // temporary v vector from V
    double *vcurrent = malloc(rows*sizeof(double));
    for (int i = 0; i < rows; i++) {
        vcurrent[i] = v[i][j];
    }

    // w = Avj
    double *w = malloc(rows*sizeof(double));
    mvmul(w, val, col_ind, row_ptr, rows, nnz, vcurrent);

    // loop to get the Hessenberg matrix
    for (int i = 0; i <= j; i++) {
        for (int k = 0; k < rows; k++) {
            h[i][j] += w[k]*v[k][i];
        }
        for (int k = 0; k < rows; k++) {
            w[k] = w[k] - h[i][j]*v[k][i];
        }
    }

    // h_j+1,j = norm(w)
    double norm = 0.0;
    for (int i = 0; i < rows; i++) {
        norm += w[i]*w[i];
    }
    norm = sqrt(norm);
    h[j+1][j] = norm;

    // v_j+1 = w/h_j+1,j
    for (int k = 0; k < rows; k++) {
        v[k][j+1] = w[k]/h[j+1][j];
    }
    
    free(vcurrent);
    free(w);
}

void mvmul(double *out, double *val, int *col_ind, int *row_ptr, int rows,
           int nnz, double *vec) {
    int r1, r2;
    double tmp;
    for (int i = 0; i < rows; i++) {
        r1 = row_ptr[i];
        if (i+1 < rows) {
            r2 = row_ptr[i+1];
        } else {
            r2 = nnz;
        }
        tmp = 0.0;
        for (int j = r1; j < r2; j++) {
            tmp += val[j]*vec[col_ind[j]];
        }
        out[i] = tmp;
    }
}
