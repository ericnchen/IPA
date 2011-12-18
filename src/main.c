/**
 * IPA*
 * 
 * IPA*, the Incompressible Problem Analyzer*, a Finite Element Method (FEM)
 * based solver for fluid flows.
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "structs.h"

int main(int argc, char *argv[]) {

    int debug = 1;

    /** Parse Input
     */
    meshinfo minf;
    solverinfo sinf;
    fluidinfo finf;
    parseinput(&minf, &finf, &sinf, debug);

    /** Generate/Get Mesh
     */
    int *ien = malloc(minf.ne*minf.nen*sizeof(int));
    int *rng = malloc(minf.ne*minf.nen*sizeof(int));
    double *xyz = malloc(minf.nn*sizeof(double));
    if (ien == NULL || rng == NULL || xyz == NULL) {
        fprintf(stderr, "(-) Could not allocate memory for genmesh.\n");
        exit(errno);
    } else {
        genmesh(minf, xyz, ien, rng, debug);
    }

    /** Formulate System
     */
    // I'll probably have to change this malloc somehow when I add quadratic
    // elements, and who knows what I'll have to do when I add SUPG. Look
    // into realloc() for that.
    double *bvec = malloc(minf.nn*sizeof(double));
    double *aval = malloc(((minf.nn-2)*3+(2)*2)*sizeof(double));
    int *aind = malloc(((minf.nn-2)*3+(2)*2)*sizeof(int));
    int *aptr = malloc(minf.nn*sizeof(int));
    if (aval == NULL || aind == NULL || aptr == NULL || bvec == NULL) {
        fprintf(stderr, "(-) Could not allocate memory for gensystem.\n");
    } else {
        gensystem(minf, finf, sinf, aval, aind, aptr, bvec, debug);
    }

    for (int i = 0; i < minf.nn; i++) {
        printf("%f\n", bvec[i]);
    }

    /** Cleanup
     */
    free(ien);
    free(rng);
    free(xyz);
    free(aval);
    free(aind);
    free(aptr);
    free(bvec);
    return EXIT_SUCCESS;
}
