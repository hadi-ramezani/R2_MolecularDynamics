#ifndef OUTPUT_H
#define	OUTPUT_H
#include <fstream>
#include "Initial.h"

using namespace std;

class Output {
public:
    ofstream outf;

    Output(const char *filename, const Configure* conf);
    void wrap(const Initial *init, const double *box, Vector *const pos);
    void Print(const int Step, const double Time, const double Ebond, const double Eangle, const double Edihedral, const double Eimproper, const double Evdw, const double Eelec, const double Ekin, const double Etotal, const double temp);
    void Print(const unsigned int frameNum, const double Ebond, const double Eangle, const double Edihedral, const double Eimproper, const double Evdw, const double Eelec, const double Etotal);
    Output(const Output& orig);
    virtual ~Output();
private:

};

#endif	/* OUTPUT_H */

