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

    /** Parse Input
     */
    meshinfo minf;
    solverinfo sinf;
    parseinput(&minf, &sinf, debug);

    /** Generate/Get Mesh
     */
    int *ien = malloc(minf.ne*minf.nen*sizeof(int));
    int *rng = malloc(minf.ne*minf.nen*sizeof(int));
    double *xyz = malloc(minf.nn*sizeof(double));
    genmesh(minf, xyz, ien, rng, debug);

    /** Formulate System
     */
    fluidinfo finf;
    int *aind;
    int *aptr = malloc(minf.nn*sizeof(int));
    gensystem(minf, finf, sinf, aind, aptr, debug);

    for (int i = 0; i < minf.nn; i++) {
        printf("%d\n", aptr[i]);
    }

    /** Cleanup
     */
    free(ien);
    free(rng);
    free(xyz);
    free(aind);
    free(aptr);
    return EXIT_SUCCESS;
}
