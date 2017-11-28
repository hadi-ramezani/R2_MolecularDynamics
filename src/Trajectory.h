#ifndef TRAJECTORY_H
#define	TRAJECTORY_H
#include <fstream>
#include "Vector.h"
#include "Configure.h"
#include "Initial.h"
 
using namespace std;

/* Define feature flags for this DCD file */
#define DCD_IS_XPLOR        0x00
#define DCD_IS_CHARMM       0x01
#define DCD_HAS_4DIMS       0x02
#define DCD_HAS_EXTRA_BLOCK 0x04
#define DCD_HAS_64BIT_REC   0x08

#define RECSCALE32BIT 1

class Trajectory {

public:

    fstream dcdf;
// Define variables to read from dcd. These all can be included in a new type
    float* X;
    float* Y;
    float* Z;
    double* boxdcd;
    int natoms;
    int nsets;
    int setsread;
    int istart;
    int nsavc;
    double delta;
    int namnf;
    int nfixed;
    int* freeind;
    float* fixedcoords;
    int reverse;
    int charmm;
    int first;
    int with_unitcell;        

    double Etot = 0.0, Ebond = 0.0, Eangle = 0.0, Evdw = 0.0, Eelec = 0.0, Ekin = 0.0;

    void write_header(const char *filename, int natoms);
    void read_header(int natoms, int nsets, int istart, int nsavc,
     double delta, int namnf, int* freeind, float* fixedcoords, int reverse,
     int* charmm);
    void open_dcd_get_info(const char *filename, int natoms);
    void write_frame(int natoms, const Vector *coor, const double* box);
    bool read_frame(int natoms, Vector *coor, double* box);
    Trajectory(const char *filename, int natoms, const Configure* conf);
    Trajectory(const Trajectory& orig);
    virtual ~Trajectory();
private:

};

#endif	/* TRAJECTORY_H */

