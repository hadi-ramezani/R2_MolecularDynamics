#ifndef NONBONDED_H
#define	NONBONDED_H

#include "Configure.h"
#include "Vector.h"
#include <vector>

class Initial;
class Parameters;
class Vector;
class LJTable;

using namespace std;

struct atominfo {
    Vector pos;
    Vector force;
    int ind;
};

// defined for a second algorithm for nonbonded calculations
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
    void build_cells(const Initial *init, const Vector *pos, Vector *f, const Configure *conf);
    void build_atomlist(const Initial *init, const Vector *pos, Vector *f, const Configure *conf);

    //second algorithm
    void build_mycells(const Initial *init, const Vector *pos, Vector *f, const Configure *conf);
    void build_neighborlist(const Initial *init, const Vector *pos, Vector *f, const Configure *conf);
    void threebody_neighborlist(const Initial *init, const Vector *pos, Vector *f, const Configure *conf);
    void mycompute(const Initial *init, const Vector *pos, 
                   Vector *f, double& Evdw, double &Eelec, const Configure *conf);

    void compute_threebody(const Initial *init, const Vector *pos, 
                   Vector *f, double& Emisc, const Configure *conf);

    void apply_pbc(const double box[3], const double box_2[3], Vector &dij);

    Vector *poshift; // Keep the location of shifted atoms for nonbonded calculation

private:

    int natoms;

    double cut2;
    double switch2;
    double pair2;
  
    // vdW switching
    double c1, c3;
    LJTable *ljTable;

    int xb, yb, zb, xytotb, totb=0;               // dimensions of decomposition

    atominfo* *boxatom;       // positions, forces, and indicies for each atom  
    int *numinbox, *maxinbox; // Number of atoms in each box
    int **nbrlist;            // List of neighbors for each box

    // for a second algorithm
    Cell *cells;
    NonbondedAtom *atoms;
    int ncellx, ncelly, ncellz;
    int *type;
    vector<int>::iterator it;
    // three body parameters
    double threebody_cut2;
    double threebody_pair2;
    double threebody_ijcut2;
};

#endif	/* NONBONDED_H */

