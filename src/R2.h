#ifndef MINIMIZER_H__
#define MINIMIZER_H__

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define R2VERSION 0.01

#if ( INT_MAX == 2147483647 )
typedef int     int32;
#else
typedef short   int32;
#endif

#define TIMEFACTOR 48.88821

void R2_die(const char *message);

typedef unsigned short Index;

typedef struct atom_name_info {
    char *resname;
    char *atomname;
    char *atomtype;
} AtomNameInfo;

typedef struct atom_constants {
    float mass;
    float charge;
    Index vdw_type;
    int32 status;   // flags telling about this atom
	int32 hydrogenList;
} Atom;

typedef struct bond {
    int32 atom1;
    int32 atom2;
    Index bond_type;
} Bond;

typedef struct angle {
    int32 atom1;
    int32 atom2;
    int32 atom3;
    Index angle_type;
} Angle;

typedef struct dihedral {
    int32 atom1;
    int32 atom2;
    int32 atom3;
    int32 atom4;
    Index dihedral_type;
} Dihedral;

typedef struct improper {
    int32 atom1;
    int32 atom2;
    int32 atom3;
    int32 atom4;
    Index improper_type;
} Improper;

class Exclusion {

public:
    Exclusion(void) : modified(0) {;}
    Exclusion(int a1, int a2, int mod = 0) :
                atom1(a1), atom2(a2), modified(mod) {;}
    int32 atom1;
    int32 atom2;
    Index modified;
    
    int hash(void) const {
                return atom1 + atom2;
    }
    
    int operator==(const Exclusion &o) const {
                return atom1 == o.atom1 && atom2 == o.atom2;
    }
    
    int operator<(const Exclusion &o) const {
            return
            (
                ( atom1 < o.atom1 ) ||
                ( atom1 == o.atom1 && atom2 < o.atom2 ) ||
                ( atom1 == o.atom1 && atom2 == o.atom2 )
            );
    }
};

#endif