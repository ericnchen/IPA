/**
 * solvesystem.c
 *
 * Solves the system. For now, only (restarted) GMRES is implemented.
 * Preconditioning will come after I've validated my results here, and
 * maybe I'll add other iterative methods like BiCGSTAB, but this is less
 * important.
 */

#include <stdio.h>
#include "structs.h"
#include "gmres.h"

void solvesystem(meshinfo minf, solverinfo sinf,
                 double *aval, int *aind, int *aptr,
                 double *bvec, double *xvec, int debug) {

    for (int i = 0; i < minf.nn; i++) {
        xvec[i] = 0.0;  // initial guess
    }

    if (debug) fprintf(stdout, "(+) Solving the system...\n");

    // call solver
    int nnz = (minf.nn-2)*3+4;
    gmres(aval, aind, aptr, minf.nn, nnz, bvec, xvec, 2, 2, 0.01);

//for (int i = 0; i < minf.nn; i++) {
//    printf("%d %f\n", i, xvec[i]);
//}

}
