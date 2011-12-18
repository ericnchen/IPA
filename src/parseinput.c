/**
 * parseinput.c
 *
 * Looks for an input file ('ipa.in') in the current working directory and
 * parses the information.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "structs.h"

void parseinput(meshinfo *minf, fluidinfo *finf, solverinfo *sinf,
                int debug) {

    FILE *pFile;
    pFile = fopen("ipa.in", "r");

    if (debug) fprintf(stdout, "(+) Parsing input...\n");

    if (pFile != NULL) {
        char linebuf[256];
        while(fgets(linebuf, sizeof(linebuf)-1, pFile)) {

            char *key = strtok(linebuf, " \t");
            char *val = strtok(NULL, "\r\n");

            /** Mesh Information
             */
            if (!strcmp(key, "ne")) {
                minf->ne = atoi(val);
            } else if (!strcmp(key, "nn")) {
                minf->nn = atoi(val);
            } else if (!strcmp(key, "nen")) {
                minf->nen = atoi(val);
            } else if (!strcmp(key, "nsd")) {
                if (atoi(val) != 1) {
                    fprintf(stderr, "(-) Only 1D is supported.\n");
                    exit(errno);
                } else {
                    minf->nsd = atoi(val);
                }
            } else if (!strcmp(key, "ndf")) {
                minf->ndf = atoi(val);

            /** Fluid Properties
             */
            } else if (!strcmp(key, "convective_velocity")) {
                finf->a = atof(val);
            } else if (!strcmp(key, "viscosity")) {
                finf->nu = atof(val);

            /** Solver Settings
             */
            } else if (!strcmp(key, "element")) {
                strncpy(sinf->element, val, 255);
            } else if (!strcmp(key, "formulation")) {
                strncpy(sinf->formulation, val, 255);

            /** Catch-All
             */
            } else if (!strcmp(key, "\n")) {
                continue;
            } else if (*key == '#') {
                continue;
            } else {
                fprintf(stderr, "(-) Unrecognized key in input file.\n");
                exit(errno);
            }
        }
    } else {
        fprintf(stderr, "(-) Could not find input file.\n");
        exit(errno);
    }

    fclose(pFile);
}
