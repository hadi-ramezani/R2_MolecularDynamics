#ifndef NONBONDED_H
#define	NONBONDED_H

#include "Configure.h"
#include "Vector.h"

class Initial;
class Parameters;
class Vector;
class LJTable;

struct atominfo {
    Vector pos;
    Vector force;
    int ind;
};

struct Pair {
  const Vector *pos1;
  const Vector *pos2;
  double kqq;   // The electrostatic factor
  double vdwA;
  double vdwB;
};

#define COULOMB 332.0636

class Nonbonded {
public:
    Nonbonded(const Initial *, const Parameters *, const Configure *);
    ~Nonbonded();

    void compute(const Initial *init, const Vector *pos, 
                   Vector *f, double& Evdw, double &Eelec, const Configure *conf);
    void Build_cells(const Initial *init, const Vector *pos, Vector *f, const Configure *conf);
    void Build_Neighborlist(const Initial *init, const Vector *pos, Vector *f, const Configure *conf);

    Vector *poshift; // Keep the location of shifted atoms for nonbonded calculation

private:

    int natoms;

    double cut2;
    double switch2;
    double pair2;
  
    // vdW switching
    double c1, c3;
    LJTable *ljTable;

    int xb, yb, zb, xytotb, totb;               // dimensions of decomposition

    atominfo* *boxatom;       // positions, forces, and indicies for each atom  
    int *numinbox, *maxinbox; // Number of atoms in each box
    int **nbrlist;            // List of neighbors for each box

};

#endif	/* NONBONDED_H */

