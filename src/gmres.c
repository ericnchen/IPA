/**
 * gmres.c
 */

#include <stdlib.h>
//#include <string.h>
#include <math.h>
#include "gmres.h"

double gmres(double *val, int *col_ind, int *row_ptr, int rows, int nnz,
             double *bvec, double *xvec) {

    // r0 = b - Ax0
    double *r   = malloc(rows*sizeof(double));
    double *ax0 = malloc(rows*sizeof(double));
    mvmul(ax0, val, col_ind, row_ptr, rows, nnz, xvec);
    vvsub(r, bvec, ax0, rows);

    // v1 = r0/norm(r0)
    double **v = malloc(rows*sizeof(double));
    for (int i = 0; i < rows; i++) {
        v[i] = malloc((rows+1)*sizeof(double));
    }
    double norm = norm2(r,r);
    for (int i = 0; i < rows; i++) {
        v[i][0] = r[i]/norm;
    }

    // g = norm(r)*e1
    double *g = malloc((rows+1)*sizeof(double));
    g[0] = norm;
    for (int i = 1; i < rows+1; i++) {
        g[i] = 0;
    }

    // outer loop allocation
    double **h = malloc((rows+1)*sizeof(double));
    double **rm = malloc(rows*sizeof(double));
    for (int i = 0; i < rows+1; i++) {
        h[i] = malloc((rows+1)*sizeof(double));
    }
    for (int i = 0; i < rows+1; i++) {
        rm[i] = malloc((rows+1)*sizeof(double));
    }
    double tempout;
    double rho = 1, eps = 1e-08;
    int j_exit;

    // inner loop allocation
    double *c = malloc((rows+1)*sizeof(double));
    double *s = malloc((rows+1)*sizeof(double));
    double tempin;

    // outer loop
    for (int j = 0; j < rows; j++) {
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

        // exit condition
        rho = fabs(g[j+1]);
        if (rho <= eps) {
            j_exit = j;
            break;
        }
    }

    // xm = x0 + Vm(Rm^-1gm)
    double *irg = malloc(rows*sizeof(double));
    for (int i = j_exit-1; i >= 0; i--) {
        irg[i] = g[i];
        for (int k = i+1; k < j_exit; k++) {
            irg[i] = irg[i] - rm[i][k]*irg[k];
        }
        irg[i] = irg[i]/rm[i][i];
    }

    // new x vector
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            xvec[i] += v[i][j]*irg[j];
        }
    }

    // deallocate
    free(r);
    free(ax0);
    free(v);
    free(g);
    free(h);
    free(rm);
    free(c);
    free(s);
    free(irg);

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
    double norm = norm2(w,w,rows);
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

double vvmul(double *vect1, double *vect2, int rows) {
             double out = 0.0;
    for (int i = 0; i < rows; i++) {
        out += vect1[i]*vect2[i];
    }
    return out;
}

double norm2(double *vect1, double *vect2, int rows) {
             double out;
    out = sqrt(vvmul(vect1, vect2, rows));
    return out;
}

void vvsub(double *out, double *vect1, double *vect2, int rows) {
    for (int i = 0; i < rows; i++) {
        out[i] = vect1[i]-vect2[i];
    }
}

void vsmul(double *out, double *vec, double scaler, int rows) {
    for (int i = 0; i < rows; i++) {
        out[i] = vec[i]*scaler;
    }
}
