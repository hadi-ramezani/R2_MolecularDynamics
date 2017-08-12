/*
 * File:   Initial.h
 * Author: Hadi and Amin
 *
 * Created on September 10, 2015, 11:28 AM
 */

#ifndef INITIAL_H
#define	INITIAL_H

#include "Configure.h"
#include "Molecule.h"
#include "Vector.h"
#include "vector"
#include "Parameters.h"

struct BondElem {
  int atom1, atom2;
  double x0;
  double k;
};

struct AngleElem {
  int atom1, atom2, atom3;
  double k, theta0;
  double k_ub, r_ub;
};

struct VdwElem {
    double epsi, rmin;
    double charge;
    int *exclist;
    int numexc;
    int type;
};

struct DpdElem {
    int *exclist;
    int numexc;
    int type;
};

struct Connections {
    int *exclist;
    int numexc;
};

class Initial {
public:

    Vector* pos;
    double* mass;

    int* residue;
    int nresidue;

    BondElem *bonds;
    AngleElem *angles;
    VdwElem *vdw;
    DpdElem *dpd;
    Connections *connections;

    float **aij; // for DPD simulation

    double box[3];

    double gamma, sigma;

    int ntype;
    float tt[3] = {1,2,3};

    bool exclude(const int iatom, const int jatom) const;
    bool exclude_dpd(const int iatom, const int jatom) const;
    void build_bondlist(const Configure *conf, const Molecule *mol, const Parameters *params);
    void build_anglelist(const Configure *conf, const Molecule *mol, const Parameters *params);
    void build_vdwlist(const Configure *conf, const Molecule *mol, const Parameters *params);
    void build_dpdlist(const Configure *conf, const Molecule *mol, const Parameters *params);

    Initial(const Configure *conf, const Molecule *mol, const Parameters *params);
    Initial(const Initial& orig);
    virtual ~Initial();
private:

};

#endif	/* INITIAL_H */
