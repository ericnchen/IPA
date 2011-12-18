/**
 * gensystem.c
 *
 * Sets up the necessary components for the system of equations.
 */

#include <stdio.h>
#include <string.h>
#include "structs.h"
#include "assemble.h"

void gensystem(meshinfo minf, fluidinfo finf, solverinfo sinf,
               double *aval, int *aind, int *aptr, double *bvec,
               int debug) {

    /** Choose the proper system to assemble.
     */
    if (minf.nsd == 1) {
        if (debug) fprintf(stdout, "(+) Assembling 1D system...\n");
        if (!strcmp(sinf.element, "linear") &&
            !strcmp(sinf.formulation, "galerkin")) {
            lineargal1d(minf, finf, aval, aind, aptr, bvec);
        } else if (!strcmp(sinf.element, "linear") &&
                   !strcmp(sinf.formulation, "supg")) {
        } else if (!strcmp(sinf.element, "quadratic") &&
                   !strcmp(sinf.formulation, "galerkin")) {
        } else if (!strcmp(sinf.element, "quadratic") &&
                   !strcmp(sinf.formulation, "supg")) {
        }
    } else if (minf.nsd == 2) {
        if (!strcmp(sinf.element, "linear") &&
            !strcmp(sinf.formulation, "galerkin")) {
        } else if (!strcmp(sinf.element, "linear") &&
                   !strcmp(sinf.formulation, "supg")) {
        } else if (!strcmp(sinf.element, "quadratic") &&
                   !strcmp(sinf.formulation, "galerkin")) {
        } else if (!strcmp(sinf.element, "quadratic") &&
                   !strcmp(sinf.formulation, "supg")) {
        }
    }

}
