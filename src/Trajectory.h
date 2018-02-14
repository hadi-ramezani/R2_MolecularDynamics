#ifndef TRAJECTORY_H
#define	TRAJECTORY_H
#include <fstream>
#include "Vector.h"
#include "Configure.h"
#include "Initial.h"
 
using namespace std;

/* Define error codes that may be returned by the DCD routines */
#define DCD_SUCCESS      0  /* No problems                     */
#define DCD_EOF         -1  /* Normal EOF                      */
#define DCD_DNE         -2  /* DCD file does not exist         */
#define DCD_OPENFAILED  -3  /* Open of DCD file failed         */
#define DCD_BADREAD     -4  /* read call on DCD file failed    */
#define DCD_BADEOF      -5  /* premature EOF found in DCD file */
#define DCD_BADFORMAT   -6  /* format of DCD file is wrong     */
#define DCD_FILEEXISTS  -7  /* output file already exists      */
#define DCD_BADMALLOC   -8  /* malloc failed                   */
#define DCD_BADWRITE    -9  /* write call on DCD file failed   */

/* Define feature flags for this DCD file */
#define DCD_IS_XPLOR        0x00
#define DCD_IS_CHARMM       0x01
#define DCD_HAS_4DIMS       0x02
#define DCD_HAS_EXTRA_BLOCK 0x04
#define DCD_HAS_64BIT_REC   0x08

#define RECSCALE32BIT 1
#define RECSCALE64BIT   2 
#define RECSCALEMAX 2

/* defines used by write_dcdstep */
#define NFILE_POS 8L
#define NSTEP_POS 20L

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
    int dcdStep;       
    long header_size, filesize, firstframesize, framesize, trajsize;

    double Etot = 0.0, Ebond = 0.0, Eangle = 0.0, Evdw = 0.0, Eelec = 0.0, Ekin = 0.0;

    void write_header(const char *filename, int natoms);
    void read_header(int natoms, int nsets, int istart, int nsavc,
     double delta, int namnf, int* freeind, float* fixedcoords, int reverse,
     int charmm, const Configure *conf);
    void open_dcd_get_info(const char *filename, int natoms, const Configure *conf);
    void write_frame(int natoms, const Vector *coor, const double* box);
    bool read_frame(int natoms, Vector *coor, double* box, const Configure *conf);
    void skip_frames(const int dcdStep, int setsread);
    void read_dcd_step(int natoms, Vector* coor, double* box);
    void print_dcderror(int errcode);
    Trajectory(const char *filename, int natoms, const Configure* conf);
    Trajectory(const Trajectory& orig);

    virtual ~Trajectory();
private:

};

#endif	/* TRAJECTORY_H */

