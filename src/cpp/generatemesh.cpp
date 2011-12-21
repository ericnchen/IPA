/**
 * generatemesh.cpp
 */

#include "generatemesh.h"
using namespace std;

void GenerateMesh::calculate_xyz(int ne, int nn) {
    double h = 1.0/ne;
    xyz.reserve(nn);
    for (int i = 0; i < nn; i++) {
        xyz[i] = i*h;
    }
}

void GenerateMesh::calculate_ien(int ne, int nen) {
    ien.reserve(ne*nen);
    for (int i = 0; i < ne*nen; i+=2) {
        ien[i] = i/2;
        ien[i+1] = i/2+1;
    }
}

//void GenerateMesh::calculate_rng(int ne, int nen) {
//    rng.reserve(ne*nen);
//    for (int i = 0; i < ne*nen; i++) {
//        rng[i] = 0;
//    }
//    rng.front() = 1;
//    rng.back() = 2;
//}
