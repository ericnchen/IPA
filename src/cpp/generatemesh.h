#ifndef GENERATEMESH_H
#define GENERATEMESH_H

#include <vector>
#include "main.h"

class GenerateMesh {
    public:
    std::vector<int> ien;
//    std::vector<int> rng;
    std::vector<double> xyz;

    GenerateMesh(MeshSettings minf) {
        calculate_xyz(minf.ne, minf.nn);
        calculate_ien(minf.ne, minf.nen);
//        calculate_rng(minf.ne, minf.nen);
    }

    private:
    void calculate_xyz(int ne, int nn);
    void calculate_ien(int ne, int nen);
//    void calculate_rng(int ne, int nen);
};

#endif
