/**
 * meshreader.h
 *
 * Reads in a mesh stored in MIXD format. For more information about the
 * MIXD format, see: http://www.cats.rwth-aachen.de/software/formats/mixd.
 *
 * Author   Eric Chen (eric@ericnchen.com)
 * Updated  28 December 2011
 *
 * Released under the MIT License, see included LICENSE file for more info.
 */

#ifndef MESHREADER_H
#define MESHREADER_H

#include <vector>
#include "inputparser.h"

class MeshReader {
    public:
//    std::vector<std::vector<double> > xyz;
    // ien
    // rng
    // constructor
    MeshReader(InputParser parameters) {
        readXYZ(parameters);
    }
    private:
//    void readXYZ(std::vector<std::vector<double> > xyz,
//                 InputParser parameters);
    void readXYZ(InputParser parameters);
};

#endif
