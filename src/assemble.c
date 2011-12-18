/**
 * assemble.c
 *
 * Generic file to hold all of my system matrix and forcing vector 
 * assembly functions.
 */

#include <stdio.h>
#include "assemble.h"
#include "structs.h"

void lineargal1d(meshinfo minf, fluidinfo finf, int *aind, int *aptr) {

    aind = malloc(((minf.nn-2)*3+4)*sizeof(int));
   free(aind); 

    aptr[0] = 0;
    aptr[1] = 2;
    for (int i = 2; i < minf.nn; i++) {
        aptr[i] = aptr[i-1]+3;
    }
}
