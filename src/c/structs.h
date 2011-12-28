/**
 * structs.h
 */

typedef struct {
    int ne, nn, nen, nsd, ndf;
} meshinfo;

typedef struct {
    char element[256];
    char formulation[256];
    char method[256];
} solverinfo;

typedef struct {
    double a;
    double nu;
} fluidinfo;
