#ifndef ANALYSIS_H
#define	ANALYSIS_H

#include "Configure.h"
#include "Initial.h"
#include "Vector.h"
#include "Nonbonded.h"
#include "Bonded.h"
#include "Trajectory.h"
#include "Output.h"
#include "R2.h"

using namespace std;

class Analysis {

public:

    Nonbonded nonbonded;
    Bonded bonded;
    Trajectory dcd;
    Output out;

    Vector* pos;
    Vector* vel;
    Vector* vel2;
    Vector* f;
    double* imass;
    double aBox[3];

    double Etot = 0.0, Ebond = 0.0, Eangle = 0.0, Evdw = 0.0, Eelec = 0.0;
    double Edihedral = 0.0 , Eimproper = 0.0;
    double Emisc = 0.0;

    Analysis(const Configure *conf, const Initial *init, const PDB *pdb, Parameters *params);
    virtual ~Analysis();
};

#endif	/* ANALYSIS_H */
