/**
 * solvesystem.c
 *
 * Solves the system. For now, only (restarted) GMRES is implemented.
 * Preconditioning will come after I've validated my results here, and
 * maybe I'll add other iterative methods like BiCGSTAB, but this is less
 * important.
 */

#include <iostream>

#include "solvesystem.h"

using namespace std;

void SolveSystem(FormulateSystem &system, MeshSettings minf, SystemSettings sinf, GMRESSettings ginf) {

    int nnz = (minf.nn-2)*3+4;
    GMRES(system.aval, system.aind, system.aptr, minf.nn, nnz,
          system.bvec, system.xvec,
          ginf.nouter, ginf.m, ginf.epsilon_outer, ginf.epsilon_inner);


for (int i = 0; i < minf.nn; i++) {
    cout << system.xvec[i] << endl;
}


}
