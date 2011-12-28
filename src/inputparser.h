/**
 * inputparser.h
 *
 * Reads and parses the IPA* input file.
 *
 * Author   Eric Chen (eric@ericnchen.com)
 * Updated  28 December 2011
 *
 * Released under the MIT License, see included LICENSE file for more info.
 */

#ifndef INPUTPARSER_H 
#define INPUTPARSER_H

class InputParser {
    public:
    // mesh parameters
    int ne;             // number of elements
    int nn;             // number of nodes
    int nen;            // number of elemental nodes
    int nsd;            // number of spatial dimensions
    int ndf;            // number of degrees of freedom
    // mesh files
    std::string mxyz;   // MIXD coordinate file
    std::string mien;   // MIXD connectivity file
    std::string mrng;   // MIXD boundary file
    // fluid parameters
    double a;           // convective velocity
    double nu;          // viscosity
    // solver parameters
    int nts;            // number of time steps
    int nouter;         // number of outer GMRES iterations
    int ninner;         // number of inner GMRES iterations

    // constructor
    InputParser() {
        parse(ne, nn, nen, nsd, ndf, mxyz, mien, mrng,
              a, nu, nts, nouter, ninner);
        checkMeshParameters(ne, nn, nen, nsd);
        checkMeshFiles(mxyz, mien, mrng);
    }

    private:
    void parse(int &ne, int &nn, int &nen, int &nsd, int &ndf,
               std::string &mxyz, std::string &mien, std::string &mrng,
               double &a, double &nu,
               int &nts, int &nouter, int &ninner);
    void checkMeshParameters(int ne, int nn, int nen, int nsd);
    void checkMeshFiles(std::string mxyz, std::string mien,
                        std::string mrng);
};

#endif
