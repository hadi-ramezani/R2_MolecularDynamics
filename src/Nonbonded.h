/*
 * File:   Nonbonded.h
 * Author: amin
 *
 * Created on September 13, 2015, 4:21 PM
 */

#ifndef NONBONDED_H
#define	NONBONDED_H

#include "iostream"
#include "Vector.h"
#include <vector>
#include "Initial.h"
#include <mkl_vsl.h>
//#include "/opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/include/mkl_vsl.h"



using namespace std;

struct Cell {
    vector<int> atoms;
    int num;
    int nbrlist[14];
};

struct NonbondedAtom {
    vector<int> nbrlist1;
    vector<int> nbrlist2;

    int numnb;
    int cell;
    int type;
};

class Nonbonded  {
public:
    VSLStreamStatePtr stream;

    Cell *cells;
    NonbondedAtom *atoms;
    int ncellx, ncelly, ncellz;
    int *type;

    vector<int>::iterator it;

    float cut, cut2, switch2, pairdist2;
    float cut_1, cut_2, switch_2, pairdist_2;
    float c1, c2, c3, c4;

    float **rmin2_table;
    float **epsi_table;
    float **kqq_table;

    Vector *poshift; // Keep the location of shifted atoms for nonbonded calculation

    void Compute(const Initial *init, const Vector *pos, Vector * const ff, const int num, double& Evdw, double& Eelec);
    void Compute_dpd(const Initial *init, const Vector *pos, const Vector *vel, Vector * const ff, const int num, double& Evdw);
    void Build_cells(const double* box, const double pairdist, const int num);
    void Neighborlist(const double* box, const int num, const Initial *init, const Vector *pos);
    void built_table(const Initial *init, const int ntype, const int num);
    Nonbonded(int num, const double rc, const double switch1, const double pairlistdist, int seed);
    Nonbonded(const Nonbonded& orig);
    virtual ~Nonbonded();
private:

};

#endif	/* NONBONDED_H */

