/**
 * gmres.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gmres.h"


// solves Ax = b for x
void gmres(double *val, int *col_ind, int *row_ptr, int rows, int nnz,
           double *rhs, double *xvec, int max, int m, double tol) {

    // allocate things
    double beta;
    double tempr, tempsqrt, tempsum;
    double *r0 = malloc(rows*sizeof(double));
    double *ax0 = malloc(rows*sizeof(double));
    double *v = malloc(rows*(m+1)*sizeof(double));
    double *h = malloc((m+1)*m*sizeof(double));
    double *w = malloc(rows*sizeof(double));
    double *vtmp = malloc(rows*sizeof(double));
    double *c = malloc(m*sizeof(double));
    double *s = malloc(m*sizeof(double));
    double *r = malloc((m+1)*m*sizeof(double));
    double *g = malloc((m+1)*sizeof(double));
    double *irg = malloc(m*sizeof(double));

    // outer GMRES loop
    for (int outer = 0; outer < max; outer++) {

        // r0 = b - Ax0
        mvmul(val, col_ind, row_ptr, rows, nnz, xvec, ax0);

        for (int i = 0; i < rows; i++) {
            r0[i] = rhs[i]-ax0[i];
        }

        // beta = norm(r0)
        beta = norm2(r0, rows);

        // v1 = r0/beta
        for (int i = 0; i < rows; i++) {
            v[i*(m+1)+0] = r0[i]/beta;
        }

        // g = norm(r0)*e1
        memset(g, 0, (m+1)*sizeof(double));
        g[0] = beta;

        // zero out old values
        memset(h, 0, (m+1)*m*sizeof(double));

        // inner GMRES loop
        for (int j = 0; j < m; j++) {
            
            // wj = Avj
            for (int i = 0; i < rows; i++) {
                vtmp[i] = v[i*(m+1)+j];
            }
            mvmul(val, col_ind, row_ptr, rows, nnz, vtmp, w);

            // compute the Hessenberg and omega matrices
            for (int i = 0; i < j; i++) {
                memset(vtmp, 0, rows*sizeof(double));
                for (int ii = 0; ii < rows; ii++) {
                    vtmp[ii] = v[ii*(m+1)+i];
                }
                h[i*m+j] = vvdot(w,vtmp, rows);
                for (int ii = 0; ii < rows; ii++) {
                    w[ii] = w[ii]-h[i*m+j]*vtmp[ii];
                }
            }
            h[(j+1)*m+j] = norm2(w, rows);

            // get the new Krylov vector
            for (int i = 0; i < rows; i++) {
                v[i*(m+1)+(j+1)] = w[i]/h[(j+1)*m+j];
//                printf("%1.16f\n", w[i]);
            }

            // Givens rotation
            r[0*m+j] = h[0*m+j];
            for (int k = 1; k < j; k++) {
                tempr = c[k-1]*r[(k-1)*m+j]+s[k-1]*h[k*m+j];
                r[k*m+j] = -s[k-1]*r[(k-1)*m+j]+c[k-1]*h[k*m+j];
                r[(k-1)*m+j] = tempr;
            }

            // solve the least-squares problem with Givens rotations
            tempsqrt = sqrt(r[j*m+j]*r[j*m+j]+h[(j+1)*m+j]*h[(j+1)*m+j]);
            c[j] = r[j*m+j]/tempsqrt;
            s[j] = h[(j+1)*m+j]/tempsqrt;
            r[j*m+j] = c[j]*r[j*m+j]+s[j]*h[(j+1)*m+j];
            g[j+1] = -s[j]*g[j];
            g[j] = c[j]*g[j];
printf("%d\t%1.16f\n", j,  fabs(g[j+1]));
        } // inner loop

//        // inverse of upper Hessenberg matrix * g
//        for (int i = m-1; i >= 0; i--) {
//            irg[i] = g[i];
//            for (int k = i+1; k < m; k++) {
//                irg[i] = irg[i]-r[i*m+k]*irg[k];
//            }
//            irg[i] = irg[i]/r[i*m+i];
//        }
//
//        // calculate the new x vector
//        for (int i = 0; i < rows; i++) {
//            tempsum = 0;
//            memset(vtmp, 0, rows*sizeof(double));
//            for (int j = 0; j < rows; j++) {
//                vtmp[j] = v[i*rows+j];
//                tempsum = vvdot(vtmp, irg, rows);
//            }
//            xvec[i] = xvec[i]+tempsum;
//        }
    } // outer loop

//    for (int i = 0; i < rows; i++) {
//      printf("%f\n", xvec[i]);
//    }

    // deallocate
    free(r0);
    free(ax0);
    free(v);
    free(h);
    free(w);
    free(vtmp);
    free(c);
    free(s);
    free(r);
    free(g);
    free(irg);
}

// solves y = Ax
 void mvmul(double *val, int *col_ind, int *row_ptr, int rows, int nnz,
            double *vec, double *out) {
    int start, end;
    double tempproduct;
    for (int i = 0; i < rows; i++) {
        start = row_ptr[i];
        if (i+1 < rows) {
            end = row_ptr[i+1];
        } else {
            end = nnz;
        }
        tempproduct = 0.0;
        for (int j = start; j < end; j++) {
            tempproduct += val[j]*vec[col_ind[j]];
        }
        out[i] = tempproduct;
    }
}

// solves dot(v1,v2)
double vvdot(double *v1, double *v2, int rows) {
    double tempproduct = 0.0;
    for (int i = rows; i < rows; i++) {
        tempproduct += v1[i]*v2[i];
    }
    return tempproduct;
}

// 2-norm of a vector
double norm2(double *val, int rows) {
    double norm = 0.0;
    for (int i = 0; i < rows; i++) {
        norm += val[i]*val[i];
    }
    norm = sqrt(norm);
    return norm;
}
