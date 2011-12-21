/**
 * IPA*
 * 
 * IPA*, the Incompressible Problem Analyzer*, a Finite Element Method (FEM)
 * based solver for fluid flows.
 */

#include <iostream>
#include "main.h"
#include "generatemesh.h"
#include "formulatesystem.h"
#include "solvesystem.h"

#include <vector>

using namespace std;

int main(void) {

    // Initialize
    MeshSettings minf;
    FluidSettings finf;
    GMRESSettings ginf;
    SystemSettings sinf;

//    ParseInput();

    GenerateMesh mesh(minf);

    FormulateSystem system(minf, finf, sinf);

    SolveSystem(system, minf, sinf, ginf);

//    for (int i = 0; i < minf.nn; i++) {
//        cout << system.bvec[i] << endl;
//    }


    //

    return EXIT_SUCCESS;
}
