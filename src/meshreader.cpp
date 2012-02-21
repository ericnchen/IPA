/**
 * meshreader.cpp
 *
 * Reads in a mesh stored in MIXD format. For more information about the
 * MIXD format, see: http://www.cats.rwth-aachen.de/software/formats/mixd.
 *
 * Author   Eric Chen (eric@ericnchen.com)
 * Updated  28 December 2011
 *
 * Released under the MIT License, see included LICENSE file for more info.
 */

#include <string>
#include <fstream>
#include <iostream>
#include "meshreader.h"

//void MeshReader::readXYZ(std::vector<std::vector<double> > xyz,
//                         InputParser parameters) {

void MeshReader::readXYZ(InputParser parameters) {
    std::ifstream file(parameters.mxyz.c_str(),
                       std::ios::in | std::ios::binary);
    if (file.is_open()) {
        std::cout << "mxyz file opened" << std::endl;
    } else {
        std::cerr << "(-) ERROR: mxyz file not found!" << std::endl;
        exit(EXIT_FAILURE);
    }
}
