/* 
 * File:   Molecule.h
 * Author: amin
 *
 * Created on September 9, 2015, 2:24 PM
 */

#ifndef MOLECULE_H
#define	MOLECULE_H
#include "Configure.h"



typedef struct bond {
    unsigned int atom1 = 0;
    unsigned int atom2 = 0;
    unsigned short bond_type = 0;
} Bond;

typedef struct angle {
    unsigned int atom1 = 0;
    unsigned int atom2 = 0;
    unsigned int atom3 = 0;
    unsigned short angle_type = 0;
} Angle;

typedef struct atom {
        // atoms name information
    int res_ID;
    char resname[5];
    char atomname[5];
    char atomtype[5];
    
    // atoms constant
    float mass;
    float charge;
    float coorx, coory, coorz;

} Atom;

class Molecule {
public:
    
    Atom *atoms;
    Bond *bonds;
    Angle *angles;
    
    Molecule(const Configure *conf);
    void read_psf(const char *filename);
    void read_pdb(const char *filename);
    
    Molecule(const Molecule& orig);
    virtual ~Molecule();
private:

};

#endif	/* MOLECULE_H */
