/*
 * File:   Integrator.h
 * Author: amin
 *
 * Created on September 11, 2015, 5:17 PM
 */

#ifndef INTEGRATOR_H
#define	INTEGRATOR_H
#define TIMEFACTOR 48.88821 

#include "Configure.h"
#include "Initial.h"
#include "Vector.h"
#include "Nonbonded.h"
#include "Bonded.h"
#include "Trajectory.h"
#include "Output.h"
#include "Thermostat.h"
#include "PDB.h"

using namespace std;

class Integrator {
public:

    Nonbonded nonbonded;
    Bonded bonded;
    Trajectory dcd;
    Output out;
    Thermostat temp;

    Vector *pos;
    Vector *vel;
    Vector *vel2;
    Vector *f;
    double *imass;

    double Etot = 0.0, Ebond = 0.0, Eangle = 0.0, Evdw = 0.0, Eelec = 0.0, Ekin = 0.0;
    double Edihedral = 0.0 , Eimproper = 0.0;
    double temperature;
    double dt, KB, lambdaDt;
    const double dtfac = 1.0/TIMEFACTOR;

    void Loop(const Configure *conf, const Initial *init);
    void Loop_dpd(const Configure *conf, const Initial *init);
    Integrator(const Configure *conf, const Initial *init, PDB *pdb, Parameters *params);
    virtual ~Integrator();
private:

};

#endif	/* INTEGRATOR_H */



