#ifndef SOLVESYSTEM_H
#define SOLVESYSTEM_H

#include <vector>
#include <algorithm>
#include <cmath>
#include "main.h"
#include "formulatesystem.h"

void SolveSystem(FormulateSystem &system, MeshSettings minf, SystemSettings sinf, GMRESSettings ginf);

void GMRES(std::vector<double> &val,
           std::vector<int> &col_ind,
           std::vector<int> &row_ptr,
           int rows,
           int nnz,
           std::vector<double> &rhs,
           std::vector<double> &xvec,
           int nouter,
           int m,
           double tolerance_outer,
           double toler_inner);

void mvmul(std::vector<double> &val,
           std::vector<int> &col_ind,
           std::vector<int> &row_ptr,
           int rows,
           int nnz,
           std::vector<double> &vec,
           std::vector<double> &out);

double vvdot(std::vector<double> &v1,
           std::vector<double> &v2,
           int rows);

double norm2(std::vector<double> &val, int rows);

#endif
