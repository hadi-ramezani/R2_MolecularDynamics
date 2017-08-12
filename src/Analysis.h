/*
 * Author: Hadi
 *
 */

#include "Configure.h"
#include "Initial.h"
#include "Vector.h"
#include "Nonbonded.h"
#include "Bonded.h"
#include "Trajectory.h"
#include "Output.h"

#ifndef ANALYSIS_H
#define	ANALYSIS_H

class Analysis {

public:

    Nonbonded nonbonded;
    Bonded bonded;
    Trajectory dcd;
    Output out;

    Vector *pos;
    Vector *vel;
    Vector *vel2;
    Vector *ff;
    double *rmass;
    double aBox[3];

    double Etot = 0.0, Ebond = 0.0, Eangle = 0.0, Evdw = 0.0, Eelec = 0.0, Ekin = 0.0;
    const double dtfac = 1.0/TIMEFACTOR;

    Analysis(const Configure *conf, const Initial *init);
    void run(const Configure *conf, const Initial *init);
    virtual ~Analysis();
};

#endif	/* ANALYSIS_H */
