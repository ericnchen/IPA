/**
 * formulatesystem.h
 */

#ifndef FORMULATE_SYSTEM
#define FORMULATE_SYSTEM

#include <vector>
#include "main.h"

class FormulateSystem {
    public:
    std::vector<double> aval, bvec, xvec;
    std::vector<int> aind, aptr;

    FormulateSystem(MeshSettings minf, FluidSettings finf, SystemSettings sinf) {
        calculate_aptr(minf.nn);
        calculate_aind(minf.ne, minf.nn);
        calculate_aval(minf.ne, minf.nn, finf.a, finf.nu);
        calculate_bvec(minf.ne, minf.nn);
        calculate_xvec(minf.nn);
    }

    private:
    void calculate_aptr(int nn);
    void calculate_aind(int ne, int nn);
    void calculate_aval(int ne, int nn, double a, double nu);
    void calculate_bvec(int ne, int nn);
    void calculate_xvec(int nn);
};

#endif
