/**
 * gmres.c
 */

#include <stdlib.h>
//#include <string.h>
#include <math.h>
#include "gmres.h"

void gmres(double *out, double *val, int *col_ind, int *row_ptr, int rows,
           int nnz, double *bvec, double *xvec, double tol) {
//    double *r = malloc(rows*sizeof(double));
//    double *g = malloc((rows+1)*sizeof(double));
//    double *c = malloc((rows+1)*sizeof(double));
//    double *s = malloc((rows+1)*sizeof(double));
//    double *tmp = malloc(rows*sizeof(double));
//    double rho = 28022011; // RIP
//
//    double **v = malloc(rows*sizeof(double));
//    double **h = malloc((rows+1)*sizeof(double));
//    double **rm = malloc(rows*sizeof(double));
//
//    for (int i = 0; i < rows; i++) {
//        g[i] = 0.0;
//        v[i] = malloc((rows+1)*sizeof(double));
//        rm[i] = malloc((rows+1)*sizeof(double));
//    }
//
//    int krylovs = 0;
//    double tmp2, sq;
//
//    do {
//        mvmul(tmp, val, col_ind, row_ptr, rows, nnz, xvec);
//        vvsub(r, bvec, tmp, rows);
//        rho = norm2(r, r);
//        for (int i = 0; i < rows; i++) {
//            v[i][0] = r[i]/rho;
//        }
//        g[0] = rho;
//        for (int j = 0; j < rows; j++) {
//            krylovs++;
//            getkrylov(val, col_ind, row_ptr, rows, nnz, j, v, h);
//            rm[0][j] = h[0][j];
//            for (int k = 1; k <= j; k++) {
//                tmp2 = c[k-1]*rm[k-1][j]+s[k-1]*h[k][j];
//                rm[k][j] = -s[k-1]*rm[k-1][j]+c[k-1]*h[k][j];
//                rm[k-1][j] = tmp2;
//            }
//            sq = sqrt(rm[j][j]*rm[j][j]+h[j+1][j]*h[j+1][j]);
//            c[j] = rm[j][j]/sq;
//            s[j] = h[j+1][j]/sq;
//            rm[j][j] = c[j]*rm[j][j]+s[j]*h[j+1][j];
//            g[j+1] = -s[j]*g[j];
//            g[j] = c[j]*g[j];
//            rho = fabs(g[j+1]);
//        }
//    } while (rho > tol);
//
//    double *irg = malloc(rows*sizeof(double));
//    for (int i = j-1; i >= 0; i--) {
//        irg[i] = g[i];
//        for (k = i+1; k < j; k++) {
//            irg[i] -= rm[i][k]*irg[k];
//        }
//        irg[i] = irg[i]/rm[i][i];
//    }
//
//    for (int i = 0; i < rows; i++) {
//        for (int j = 0; j < rows; j++) {
//            x[i] += v[i][j]*irg[j];
//        }
//    }
//
//    free(r);
//    free(v);
//    free(g);
//    free(tmp);
//    free(v);
}

void getkrylov(double *val, int *col_ind, int *row_ptr, int rows, int nnz,
               int j, double **v, double **h) {

    double *vcurrent = malloc(rows*sizeof(double));
    for (int i = 0; i < rows; i++) {
        vcurrent[i] = v[i][j];
    }

    double *w = malloc(rows*sizeof(double));
    mvmul(w, val, col_ind, row_ptr, rows, nnz, vcurrent);

    for (int i = 0; i <= j; i++) {
        for (int k = 0; k < rows; k++) {
            h[i][j] += w[k]   *v[k][i];
            w[k]    -= h[i][j]*v[k][i];
        }
    }

    double norm;
    norm = norm2(w, w, rows);
    h[j+1][j] = norm;

    for (int k = 0; k < rows; k++) {
        v[k][j+1] = w[k]/h[j+1][j];
    }
    
    free(w);
    free(vcurrent);
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
