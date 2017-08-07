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

using namespace std;

class Trajectory {
public:
    ofstream dcdf;

    float *X;
    float *Y;
    float *Z;
    double *boxdcd;

    void frames(int natoms, const Vector *coor, const double *box);
    Trajectory(const char *filename,int natoms);
    Trajectory(const Trajectory& orig);
    virtual ~Trajectory();
private:

};

#endif	/* TRAJECTORY_H */

