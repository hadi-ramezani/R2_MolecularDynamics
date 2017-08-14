#ifndef NONBONDED_H
#define	NONBONDED_H

#include "Configure.h"
class Initial;
class Parameters;
class Vector;
class LJTable;

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
                   Vector *f, double& Evdw, double &Eelec);
private:
    int natoms;

    double cut2;
    double switch2;
    double pair2;
  
    // vdW switching
    double c1, c3;
    LJTable *ljTable;
};

#endif	/* NONBONDED_H */

