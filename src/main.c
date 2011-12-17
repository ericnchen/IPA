/**
 * IPA*
 * 
 * IPA*, the Incompressible Problem Analyzer*, a Finite Element Method (FEM)
 * based solver for fluid flows.
 */

#include <stdlib.h>
#include "main.h"
#include "structs.h"

#include <stdio.h>

int main(void) {

    int debug = 1;
    meshinfo minf;

    parseinput(&minf, debug);

    int *ien = (int *)malloc(minf.ne*minf.nen*sizeof(int));
    int *rng = (int *)malloc(minf.ne*minf.nen*sizeof(int));
    double *xyz = (double *)malloc(minf.nn*sizeof(double));
    genmesh(minf, xyz, ien, rng, debug);

    free(ien);
    free(rng);
    free(xyz);
    return EXIT_SUCCESS;
}
