/**
 * inputparser.cpp
 *
 * Reads and parses the IPA* input file.
 *
 * Author   Eric Chen (eric@ericnchen.com)
 * Updated  28 December 2011
 *
 * Released under the MIT License, see included LICENSE file for more info.
 */

#include <fstream>
#include <sstream>
#include "inputparser.h"

void InputParser::parse(int &ne, int &nn, int &nen, int &nsd, int &ndf,
                        std::string &mxyz,
                        std::string &mien, std::string &mrng,
                        double &a, double &nu,
                        int &nts, int &nouter, int &ninner) {
    std::ifstream input_file("ipa.in");
    if (input_file.is_open()) {
        std::string line, key, token;
        while (getline(input_file, line)) {
            // read lines
            if (line == "") continue;   // skip blank lines
            std::stringstream iss(line);
            iss >> key >> token;        // tokenize
            if (key == "#") continue;   // skip comment lines
            // mesh parameters
            if      (key == "ne") ne = atoi(token.c_str());
            else if (key == "nn") nn = atoi(token.c_str());
            else if (key == "nen") nen = atoi(token.c_str());
            else if (key == "nsd") nsd = atoi(token.c_str());
            else if (key == "ndf") ndf = atoi(token.c_str());
            // mesh files
            else if (key == "mxyz") mxyz = token;
            else if (key == "mien") mien = token;
            else if (key == "mrng") mrng = token;
            // fluid parameters
            else if (key == "convective_velocity") a = atof(token.c_str());
            else if (key == "viscosity") nu = atof(token.c_str());
            // solver parameters
            else if (key == "nts") nts = atoi(token.c_str());
            else if (key == "nouter") nouter = atoi(token.c_str());
            else if (key == "ninner") ninner = atoi(token.c_str());
        }
    } else {
        std::cerr << "(-) ERROR: Input file not found!" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void InputParser::checkMeshParameters(int ne, int nn, int nen, int nsd) {
    if (nsd == 1) {
        if (nn != ne+1) {
            std::cerr << "(-) ERROR: Incorrect nn for 1D!" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (nen != 2) {
            std::cerr << "(-) ERROR: Incorrect nen for 1D!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else if (nsd == 2) {
        std::cerr << "(-) ERROR: 2D not yet implemented!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else {
        std::cerr << "(-) ERROR: Unsupported nsd!" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void InputParser::checkMeshFiles(std::string mxyz, std::string mien,
                                 std::string mrng) {
    if (mxyz == "") {
        std::cerr << "(-) ERROR: No mxyz file specified!" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (mrng == "") {
        std::cerr << "(-) ERROR: No mrng file specified!" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (mrng == "") {
        std::cerr << "(-) ERROR: No mrng file specified!" << std::endl;
        exit(EXIT_FAILURE);
    }
}
