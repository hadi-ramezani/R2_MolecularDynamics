#ifndef OUTPUT_H
#define	OUTPUT_H
#include <fstream>
#include "Initial.h"

using namespace std;

class Output {
public:
    ofstream outf;
    ofstream outf_threebody;
    int width = 16;
    Output(const char *filename, const Configure* conf);
    Output(const Configure* conf);
    void wrap(const Initial *init, const double *box, Vector *const pos);
    void print(const int Step, const double Time, const double Ebond, 
        const double Eangle, const double Edihedral, const double Eimproper,
        const double Evdw, const double Eelec, const double Ekin,
        const double Etot, const double temp);
    void print(const unsigned int frameNum, const double Ebond, 
        const double Eangle, const double Edihedral, const double Eimproper,
        const double Evdw, const double Eelec, const double Emisc, const double Etot);
    void print_threebody(const double rik, const double xik, const double yik, const double zik, const double thetakij,
        const double energy);
    Output();
    virtual ~Output();
private:

};

#endif	/* OUTPUT_H */

