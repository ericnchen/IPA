/**
 * gmres.c
 */

#include <stdlib.h>
#include <math.h>
#include "gmres.h"

#include <stdio.h>
#include <errno.h>

void gmres(double *val, int *col_ind, int *row_ptr, int rows, int nnz,
           double *bvec, double *xvec) {

    // solver parameters
    int m = 10; // restart
    double eps = 1e-08; //tolerance

    // r0 = b - Ax0
    double *r   = malloc(rows*sizeof(double));
    double *ax0 = malloc(rows*sizeof(double));
    if (r == NULL || ax0 == NULL) {
        fprintf(stderr, "(-) Could not allocate memory for GMRES.\n");
        exit(errno);
    }
    mvmul(ax0, val, col_ind, row_ptr, rows, nnz, xvec);
    for (int i = 0; i < rows; i++) {
        r[i] = bvec[i]-ax0[i];
    }

    // rho = norm(r)
    double rho = 0.0;
    for (int i = 0; i < rows; i++) {
        rho += r[i]*r[i];
    }
    rho = sqrt(rho);

    // v1 = r0/norm(r0)
    double *v = malloc(rows*(m+1)*sizeof(double));
    if (v == NULL) {
        fprintf(stderr, "(-) Could not allocate memory for GMRES.\n");
        exit(errno);
    }
    double norm = 0.0;
    for (int i = 0; i < rows; i++) {
        norm += r[i]*r[i];
    }
    norm = sqrt(norm);
    for (int i = 0; i < rows; i++) {
        v[i*(m+1)+0] = r[i]/norm;
    }

    // g = norm(r)*e1
    double *g = malloc((m+1)*sizeof(double));
    if (g == NULL) {
        fprintf(stderr, "(-) Could not allocate memory for GMRES.\n");
        exit(errno);
    }
    g[0] = norm;
    for (int i = 1; i <= m; i++) {
        g[i] = 0.0;
    }

    // outer GMRES loop allocation
    double *h = malloc((m+1)*m*sizeof(double));
    double *rm = malloc((m+1)*m*sizeof(double));
    if (h == NULL || rm == NULL) {
        fprintf(stderr, "(-) Could not allocate memory for GMRES.\n");
        exit(errno);
    }
    for (int i = 0; i < (m+1)*m; i ++) h[i] = 0.0;
    double tempout;
    int j_exit = m;

    // inner GMRES loop allocation
    double *c = malloc(m*sizeof(double));
    double *s = malloc(m*sizeof(double));
    if (c == NULL || s == NULL) {
        fprintf(stderr, "(-) Could not allocate memory for GMRES.\n");
        exit(errno);
    }
    double tempin;

    // restart GMRES loop
    while (rho >= eps) {

        // outer loop
        for (int j = 0; j < m; j++) {
            getkrylov(val, col_ind, row_ptr, rows, nnz, j, v, h, m);
            rm[0*m+j] = h[0*m+j];

            // inner loop
            for (int k = 1; k <= j; k++) {
                tempin = c[k-1]*rm[(k-1)*m+j]+s[k-1]*h[k*m+j];
                rm[k*m+j] = -s[k-1]*rm[(k-1)*m+j]+c[k-1]*h[k*m+j];
                rm[(k-1)*m+j] = tempin;
            }
    
            // steps 8 and 9
            tempout = rm[j*m+j]*rm[j*m+j]+h[(j+1)*m+j]*h[(j+1)*m+j];
            tempout = sqrt(tempout);
            c[j] = rm[j*m+j]/tempout;
            s[j] = h[(j+1)*m+j]/tempout;
            rm[j*m+j] = c[j]*rm[j*m+j]+s[j]*h[(j+1)*m+j];
            g[j+1] = -s[j]*g[j];
            g[j] = c[j]*g[j];
    
            // convergence check
            rho = fabs(g[j+1]);
            printf("%d %f\n", j, rho);
    //        if (rho <= 1e-08) {
    //            j_exit = j;
    //            break;
    //        }
        }
        rho = 1e-09;
    }

//
//    // xm = x0 + Vm(Rm^-1gm)
//    double *irg = malloc(m*sizeof(double));
//    if (irg == NULL) {
//        fprintf(stderr, "(-) Could not allocate memory for GMRES.\n");
//        exit(errno);
//    }
//    for (int i = j_exit-1; i >= 0; i--) {
//        irg[i] = g[i];
//        for (int k = i+1; k < j_exit; k++) {
//            irg[i] = irg[i] - rm[i*m+k]*irg[k];
//        }
//        irg[i] = irg[i]/rm[i*m+i];
//    }
//
//    // new x vector
//    for (int i = 0; i < rows; i++) {
//        for (int j = 0; j < m; j++) {
//            xvec[i] += v[i*(m+1)+j]*irg[j];
//        }
//    }
//
    // deallocate
    free(r);
    free(ax0);
    free(v);
    free(g);
    free(h);
    free(rm);
    free(c);
    free(s);
//    free(irg);
//
//    return rho_loc;
}

void getkrylov(double *val, int *col_ind, int *row_ptr, int rows, int nnz,
               int j, double *v, double *h, int m) {

    // temporary v vector from V
    double *vcurrent = malloc(rows*sizeof(double));
    if (vcurrent == NULL) {
        fprintf(stderr, "(-) Could not allocate memory for GMRES.\n");
        exit(errno);
    }
    for (int i = 0; i < rows; i++) {
        vcurrent[i] = v[i*(m+1)+j];
    }

    // w = Avj
    double *w = malloc(rows*sizeof(double));
    if (w == NULL) {
        fprintf(stderr, "(-) Could not allocate memory for GMRES.\n");
        exit(errno);
    }
    mvmul(w, val, col_ind, row_ptr, rows, nnz, vcurrent);

    // loop to get the Hessenberg matrix
    for (int i = 0; i <= j; i++) {
        for (int k = 0; k < rows; k++) {
            h[i*m+j] = h[i*m+j]+ w[k]*v[k*(m+1)+i];
        }
        for (int k = 0; k < rows; k++) {
            w[k] = w[k]-h[i*m+j]*v[k*(m+1)+i];
        }
    }

    // h_j+1,j = norm(w)
    double norm = 0.0;
    for (int i = 0; i < rows; i++) {
        norm += w[i]*w[i];
    }
    norm = sqrt(norm);
    h[(j+1)*m+j] = norm;

    // v_j+1 = w/h_j+1,j
    for (int k = 0; k < rows; k++) {
        v[k*(m+1)+(j+1)] = w[k]/h[(j+1)*m+j];
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
