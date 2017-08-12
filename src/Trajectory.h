/*
 * File:   Trajectory.h
 * Author: amin
 *
 * Created on September 24, 2015, 12:20 PM
 */

#ifndef TRAJECTORY_H
#define	TRAJECTORY_H
#include <fstream>
#include "Vector.h"
#include "Configure.h"
#include "Initial.h"

using namespace std;

class Trajectory {
public:
    fstream dcdf;

    float *X;
    float *Y;
    float *Z;
    double *boxdcd;
    double Etot = 0.0, Ebond = 0.0, Eangle = 0.0, Evdw = 0.0, Eelec = 0.0, Ekin = 0.0;

    void WriteHeader(const char *filename, int natoms);
    void ReadHeader(const char *filename, int natoms);
    void WriteFrame(int natoms, const Vector *coor, const double* box);
    bool ReadFrame(int natoms, Vector *coor, double* box);
    Trajectory(const char *filename, int natoms, const Configure* conf);
    Trajectory(const Trajectory& orig);
    virtual ~Trajectory();
private:

};

#endif	/* TRAJECTORY_H */

