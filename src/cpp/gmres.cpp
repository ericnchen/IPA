#include <iostream>
#include "solvesystem.h"
using namespace std;

void GMRES(vector<double> &val,
           vector<int> &col_ind,
           vector<int> &row_ptr,
           int rows,
           int nnz,
           vector<double> &rhs,
           vector<double> &xvec,
           int nouter,
           int m,
           double tolerance_outer,
           double toler_inner) {

    // allocate
    double beta;
    double tempr, tempsqrt;
    vector<double> ax0(rows), r0(rows), vtemp(rows), w(rows);
    vector<double> c(m), s(m);
    vector<double> g(m+1);
    vector<vector<double> > v(rows, vector<double>(m+1,0));
    vector<vector<double> > h(m+1, vector<double>(m));
    vector<vector<double> > r(m+1, vector<double>(m));


    for (int outer = 0; outer < nouter; outer++) {
        // r0 = b- ax0
        mvmul(val, col_ind, row_ptr, rows, nnz, xvec, ax0);
        for (int i = 0; i < rows; i++) {
            r0[i] = rhs[i]-ax0[i];
        }

        // beta = norm(r0)
        beta = norm2(r0, rows);

        //v1 = r0/beta
        for (int i = 0; i < rows; i++) {
            v[i][0] = r0[i]/beta;
        }

        // g = norm(r0)*e1
        fill(g.begin(), g.end(), 0);
        g[0] = beta;

        // inner loop
        for (int j = 0; j < m; j++) {

            // wj = Avj
            for (int i = 0; i < rows; i++) {
                vtemp[i] = v[i][j];
            }
            mvmul(val, col_ind, row_ptr, rows, nnz, vtemp, w);

            // compute Hessenberg and omega matrices
            for (int i = 0; i < j; i++) {
                fill(vtemp.begin(), vtemp.end(), 0);
                for (int ii = 0; ii < rows; ii++) {
                    vtemp[ii] = v[ii][i];
                }
                h[i][j] = vvdot(w, vtemp, rows);
                for (int ii = 0; ii < rows; ii++) {
                    w[ii] = w[ii]-h[i][j]*vtemp[ii];
                }
            }
            h[j+1][j] = norm2(w, rows);

            // get the new Krylov vector
            for (int i = 0; i < rows; i++) {
                v[i][j+1] = w[i]/h[j+1][j];
            }

            // givens rotation
            r[0][j] = h[0][j];
            for (int k = 1; k < j; k++) {
                tempr = c[k-1]*r[k-1][j]+s[k-1]*h[k][j];
                r[k][j] = -s[k-1]*r[k-1][j]+c[k-1]*h[k][j];
                r[k-1][j] = tempr;
            }

            // solve th least-squares problem with Givens rotations
            tempsqrt = sqrt(r[j][j]*r[j][j]+h[j+1][j]*h[j+1][j]);
            c[j] = r[j][j]/tempsqrt;
            s[j] = h[j+1][j]/tempsqrt;
            r[j][j] = c[j]*r[j][j]+s[j]*h[j+1][j];
            g[j+1] = -s[j]*g[j];
cout << fabs(g[j+1]) << endl;
            g[j] = c[j]*g[j];
        }

    }
}

void mvmul(vector<double> &val,
           vector<int> &col_ind,
           vector<int> &row_ptr,
           int rows,
           int nnz,
           vector<double> &vec,
           vector<double> &out) {
    double tempproduct;
    for (int i = 0; i < rows; i++) {
        tempproduct = 0.0;
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) {
            tempproduct += val[j]*vec[col_ind[j]];
        }
        out[i] = tempproduct;
    }
}

double vvdot(vector<double> &v1,
           vector<double> &v2,
           int rows) {
    double dotproduct(0.0);
    for (int i = rows; i < rows; i++) {
        dotproduct += v1[i]*v2[i];
    }
    return dotproduct;
}

double norm2(vector<double> &val, int rows) {
    double norm(0.0);
    for (int i = 0; i < rows; i++) {
        norm += val[i]*val[i];
    }
    norm = sqrt(norm);
    return norm;
}
