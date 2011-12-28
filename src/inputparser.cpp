/**
 * inputparser.cpp
 *
 * Parses the IPA* input file.
 *
 * Author   Eric Chen (eric@ericnchen.com)
 * Updated  28 December 2011
 *
 * Released under the MIT License, see included LICENSE file for more info.
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include "inputparser.h"

void InputParser::parse(int &ne, int &nn, int &nen, int &nsd, int &ndf,
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
