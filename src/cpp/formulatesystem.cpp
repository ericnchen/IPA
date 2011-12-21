/**
 * formulatesystem.h
 */

#include "formulatesystem.h"

void FormulateSystem::calculate_aptr(int nn) {
    aptr.reserve(nn+1);
    aptr[0] = 0;
    aptr[1] = 2;
    for (int i = 2; i < nn; i++) {
        aptr[i] = aptr[i-1]+3;
    }
    aptr[nn] = aptr[nn-1]+2;
}

void FormulateSystem::calculate_aind(int ne, int nn) {
    int nnz((nn-2)*3+(2)*2);
    aind.reserve(nnz);
    aind[0] = 0;
    aind[1] = 1;
    aind[nnz-2] = ne-1;
    aind[nnz-1] = ne;
    for (int i = 0; i < nn-2; i++) {
        aind[aptr[i+1]  ] = i;
        aind[aptr[i+1]+1] = i+1;
        aind[aptr[i+1]+2] = i+2;
    }
}

void FormulateSystem::calculate_aval(int ne, int nn, double a, double nu) {
    int nnz((nn-2)*3+(2)*2);
    double h(1.0/ne);
    double valupper( a/2-nu/h);
    double vallower(-a/2-nu/h);
    double valdiag (   2*nu/h);
    aval.reserve(nnz);
    aval[0] = -a/2+nu/h;
    aval[1] = valupper;
    aval[nnz-2] = vallower;
    aval[nnz-1] = a/2+nu/h;
    for (int i = 0; i < nn-2; i++) {
        aval[aptr[i+1]  ] = vallower;
        aval[aptr[i+1]+1] = valdiag;
        aval[aptr[i+1]+2] = valupper;
    }
}

void FormulateSystem::calculate_bvec(int ne, int nn) {
    double h(1.0/ne);
    bvec.reserve(nn);
    bvec[0] = h/2;
    bvec[nn-1] = h/2;
    for (int i = 1; i < nn-1; i++) {
        bvec[i] = h;
    }
}

void FormulateSystem::calculate_xvec(int nn) {
    xvec.reserve(nn);
    for (int i = 0; i < nn; i++) {
        xvec[i] = 0.0;
    }
}
