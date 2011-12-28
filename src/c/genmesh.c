/**
 * genmesh.c
 *
 * In the future this function will read MIXD format mesh files, but for
 * now, a temporary implementation will be to just generate a simple 1D
 * mesh.
 */

#include <stdio.h>
#include <stdlib.h>
#include "structs.h"

void genmesh(meshinfo minf, double *xyz, int *ien, int *rng, int debug) {

    if (debug) fprintf(stdout, "(+) Creating 1D mesh...\n");

    double h = 1.0/minf.ne;

    for (int i = 0; i < minf.nn; i++) {
        xyz[i] = i*h;
    }

    for (int i = 0; i < minf.ne*minf.nen; i+=2) {
        ien[i] = i/2;
        ien[i+1] = i/2+1;
    }

    for (int i = 0; i < minf.ne*minf.nen; i++) {
        rng[i] = 0;
    }
    rng[0] = 1;
    rng[minf.ne*minf.nen-1] = 2;
}
