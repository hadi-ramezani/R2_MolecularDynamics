/*
 * File:   Thermostat.h
 * Author: amin
 *
 * Created on November 26, 2015, 10:32 AM
 */

#ifndef THERMOSTAT_H
#define	THERMOSTAT_H

#include <random>
#include <mkl_vsl.h>
 //#include "/opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/include/mkl_vsl.h"
#include "Initial.h"
#include "Configure.h"
#include "Vector.h"
#include "Nonbonded.h"
#include "Initial.h"

class Thermostat {
public:
        /* initialize random seed: */
    std::default_random_engine generator, ran;

    VSLStreamStatePtr stream;


    int ndof; //number of degree of freedom
    double target;
    double gamma;
    double sigma;
    double dt_s;
    double cut_1;
    double kTs;
    double KB;


    void Initial_vel(const int num, const Initial *init, Vector *const vel);
    void DPD(const int num, const double *dr, const double *dot, double * const dpdF);
    void LA(const int num, const Initial *init, const double *rmass, const Vector *pos, const Nonbonded *nonbonded, Vector * const vel);
    Thermostat(const Configure *conf);
    Thermostat(const Thermostat& orig);
    virtual ~Thermostat();
private:

};

#endif	/* THERMOSTAT_H */

