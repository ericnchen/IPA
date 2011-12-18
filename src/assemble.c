/**
 * assemble.c
 *
 * Generic file to hold all of my system matrix and forcing vector 
 * assembly functions.
 */

#include <stdio.h>
#include "assemble.h"
#include "structs.h"

void lineargal1d(meshinfo minf, fluidinfo finf,
                 double *aval, int *aind, int *aptr) {

    aptr[0] = 0;
    aptr[1] = 2;
    for (int i = 2; i < minf.nn; i++) {
        aptr[i] = aptr[i-1]+3;
    }

    aind[0] = 0;
    aind[1] = 1;
    aind[(minf.nn-2)*3+(2)*2-2] = 9;
    aind[(minf.nn-2)*3+(2)*2-1] = 10;
    for (int i = 0; i < minf.nn-2; i++) {
        aind[aptr[i+1]] = i;
        aind[aptr[i+1]+1] = i+1;
        aind[aptr[i+1]+2] = i+2;
    }

    double h = 1.0/minf.ne;
    double valupper =  finf.a/2-finf.nu/h;
    double vallower = -finf.a/2-finf.nu/h;
    double valdiag  =         2*finf.nu/h;

    aval[0] = -finf.a/2+finf.nu/h;
    aval[1] = valupper;
    aval[(minf.nn-2)*3+(2)*2-2] = vallower;
    aval[(minf.nn-2)*3+(2)*2-1] = finf.a/2+finf.nu/h;
    for (int i = 0; i < minf.nn-2; i++) {
        aval[aptr[i+1]] = vallower;
        aval[aptr[i+1]+1] = valdiag;
        aval[aptr[i+1]+2] = valupper;
    }
}
