#ifndef INITIAL_H
#define	INITIAL_H

#include "Configure.h"
#include "Vector.h"
#include "UniqueSet.h"
#include "Hydrogen.h"
#include "Parameters.h"

class Parameters;
class PDB;

// List maintaining the global atom indicies sorted by helix groups.
class Initial {

    typedef struct constraint_params {
        double k;        //  Force constant
        Vector refPos;   //  Reference position for restraint
    } ConstraintParams;

private:
    Atom *atoms;        //  Array of atom structures
    AtomNameInfo *atomNames;//  Array of atom name info.  Only maintained
                //  on node 0 for VMD interface
    Bond *bonds;        //  Array of bond structures
    Angle *angles;      //  Array of angle structures
    Dihedral *dihedrals;    //  Array of dihedral structures
    Improper *impropers;    //  Array of improper structures
    Exclusion *exclusions;  //  Array of exclusion structures
    UniqueSet<Exclusion> exclusionSet;  //  Used for building
    int* consIndexes;   //  Constraint indexes for each atom
    ConstraintParams *consParams;
                //  Parameters for each atom constrained
    double* langevinParams;   //  b values for langevin dynamics
    double* langForceVals;    //  Calculated values for langvin random forces
    int* fixedAtomFlags;    //  1 for fixed, -1 for fixed group, else 0
    double* rigidBondLengths;  //  if H, length to parent or 0. or
                //  if not H, length between children or 0.

    int* *bondsWithAtom;    //  List of bonds involving each atom
    int* *bondsByAtom;  //  List of bonds owned by each atom
    int* *anglesByAtom;     //  List of angles owned by each atom
    int* *dihedralsByAtom;  //  List of dihedrals owned by each atom
    int* *impropersByAtom;  //  List of impropers owned by each atom
    int* *exclusionsByAtom; //  List of exclusions owned by each atom

    int* *all_exclusions;
                //  List of all exclusions, including
                //  explicit exclusions and those calculated
                //  from the bonded structure based on the
                //  exclusion policy
    int* *onefour_exclusions;
                //  List of 1-4 interactions.  This list is
                //  used only if the exclusion policy is 
                //  scaled1-4 to track 1-4 interactions that
                //  need to be handled differently

    void build_lists_by_atom();
                //  Build the list of structures by atom
    

    void read_atoms(FILE *, Parameters *);
                //  Read in atom info from .psf
    void read_bonds(FILE *, Parameters *);
                //  Read in bond info from .psf
    void read_angles(FILE *, Parameters *);
                //  Read in angle info from .psf
    void read_dihedrals(FILE *, Parameters *);
                //  Read in dihedral info from .psf
    void read_impropers(FILE *, Parameters *);
                //  Read in improper info from .psf
    void read_donors(FILE *);
                //  Read in hydrogen bond donors from .psf
    void read_acceptors(FILE *);
                //  Read in hydrogen bond acceptors from .psf
    void read_exclusions(FILE *);
                //  Read in exclusion info from .psf

    void build12excl(void);
    void build13excl(void);
    void build14excl(int);
    void stripHGroupExcl(void);
    void build_exclusions();

    // analyze the atoms, and determine which are oxygen, etc.
    // this is called after a molecule is sent our (or received in)
    void build_atom_status(void);

public:
    int numAtoms;       //  Number of atoms 
    int numBonds;       //  Number of bonds
    int numAngles;      //  Number of angles
    int numDihedrals;   //  Number of dihedrals
    int numImpropers;   //  Number of impropers
    int numExclusions;  //  Number of exclusions
    int numTotalExclusions; //  double Total Number of Exclusions // hack
    int numDonors;
    int numAcceptors;
    int numConstraints; //  Number of atoms constrained
    int numFixedAtoms;  //  Number of fixed atoms
    int numHydrogenGroups;  //  Number of hydrogen groups
    int numRigidBonds;  //  Number of rigid bonds
    int numFixedRigidBonds; //  Number of rigid bonds between fixed atoms

    // The following are needed for error checking because we
    // eliminate bonds, etc. which involve only fixed atoms
    int numCalcBonds;   //  Number of bonds requiring calculation
    int numCalcAngles;  //  Number of angles requiring calculation
    int numCalcDihedrals;   //  Number of dihedrals requiring calculation
    int numCalcImpropers;   //  Number of impropers requiring calculation
    int numCalcExclusions;  //  Number of exclusions requiring calculation

    //  Number of dihedrals with multiple periodicity
    int numMultipleDihedrals; 
    //  Number of impropers with multiple periodicity
    int numMultipleImpropers; 
    // indexes of "atoms" sorted by hydrogen groups
    HydrogenGroup hydrogenGroup;
    int waterIndex;

    Initial(const Configure *conf, Parameters *params);
    ~Initial();        //  Destructor

    void read_psf_file(const char *, Parameters *);
            //  Read in a .psf file given
            //  the filename and the parameter
            //  object to use

#if 0   
    void build_constraint_params(StringList *, StringList *, StringList *,
                     PDB *, char *);
                //  Build the set of harmonic constraint 
                // parameters

    void build_langevin_params(double coupling, int doHydrogen);
    void build_langevin_params(StringList *, StringList *, PDB *, char *);
                //  Build the set of langevin dynamics parameters

    void build_fixed_atoms(StringList *, StringList *, PDB *, char *);
                //  Determine which atoms are fixed (if any)
#endif

    int is_hydrogen(int);     // return true if atom is hydrogen
    int is_oxygen(int);       // return true if atom is oxygen
    int is_hydrogenGroupParent(int); // return true if atom is group parent
    int is_water(int);        // return true if atom is part of water 
    int  get_groupSize(int);     // return # atoms in (hydrogen) group
    int get_mother_atom(int);  // return mother atom of a hydrogen

    //  Get the mass of an atom
    double atommass(int anum) const {
        return(atoms[anum].mass);
    }

    //  Get the residue name of an atom
    char* get_resname(int anum) const {
        return(atomNames[anum].resname);
    }

    //  Get the segment name of an atom
    char* get_segname(int anum) const {
        return(atomNames[anum].segname);
    }

    //  Get the resid from atom index
    int get_resnum(int anum) const {
        return(atomNames[anum].resnum);
    }

    //  Get the charge of an atom
    double atomcharge(int anum) const {
        return(atoms[anum].charge);
    }
    
    //  Get the vdw type of an atom
    Index atomvdwtype(int anum) const {
        return(atoms[anum].vdw_type);
    }

    //  Retrieve a bond structure
    Bond *get_bond(int bnum) const {return (&(bonds[bnum]));}

    //  Retrieve an angle structure
    Angle *get_angle(int anum) const {return (&(angles[anum]));}

    //  Retrieve an improper strutcure
    Improper *get_improper(int inum) const {return (&(impropers[inum]));}

    //  Retrieve a dihedral structure
    Dihedral *get_dihedral(int dnum) const {return (&(dihedrals[dnum]));}

    //  Retrieve an exclusion structure
    Exclusion *get_exclusion(int ex) const {return (&(exclusions[ex]));}

    //  Retrieve an atom type
    const char *get_atomtype(int anum) const {
        if (atomNames == NULL) {
            R2_die("Tried to find atom type on node other than node 0");
        }

        return(atomNames[anum].atomtype);
    }

    //  Lookup atom id from segment, residue, and name
    int get_atom_from_name(const char *segid, int resid, const char *aname) const;

    //  Lookup number of atoms in residue from segment and residue
    int get_residue_size(const char *segid, int resid) const;

    //  Lookup atom id from segment, residue, and index in residue
    int get_atom_from_index_in_residue(const char *segid, int resid, int index) const;

    
    //  The following routines are used to get the list of bonds
    //  for a given atom.  This is used when creating the bond lists
    //  for the force objects
    int *get_bonds_for_atom(int anum) { return bondsByAtom[anum]; }
    int *get_angles_for_atom(int anum) { return anglesByAtom[anum]; }
    int *get_dihedrals_for_atom(int anum) { return dihedralsByAtom[anum]; }
    int *get_impropers_for_atom(int anum) { return impropersByAtom[anum]; }
    int *get_exclusions_for_atom(int anum) { return exclusionsByAtom[anum]; }
    
    //  Check for exclusions, either explicit or bonded.
    //  Inline this funcion since it is called so often
    int checkexcl(int atom1, int atom2) const {
        register int check_int;  //  atom whose array we will search
        int other_int;   //  atom we are looking for

       //  We want to search the array of the smaller atom
       if (atom1<atom2) {
            check_int = atom1;
            other_int = atom2;
       }
       else {
            check_int = atom2;
            other_int = atom1;
       }

       //  Do the search and return the correct value
       register int *list = all_exclusions[check_int];
       check_int = *list;
       while( check_int != other_int && check_int != -1 ) {
          check_int = *(++list);
       }
       return ( check_int != -1 );
    }
    
    //  Check for 1-4 exclusions.  This is only valid when the
    //  exclusion policy is set to scaled1-4. Inline this function
    //  since it will be called so often
    int check14excl(int atom1, int atom2) const {
        register int check_int;  //  atom whose array we will search
        int other_int;   //  atom we are looking for

       //  We want to search the array of the smaller atom
        if (atom1<atom2) {
            check_int = atom1;
            other_int = atom2;
        }
        else {
            check_int = atom2;
            other_int = atom1;
        }

        //  Do the search and return the correct value
        register int *list = onefour_exclusions[check_int];
        check_int = *list;
        while( check_int != other_int && check_int != -1 ) {
            check_int = *(++list);
        }
        return ( check_int != -1 );
    }

    //  Return true or false based on whether the specified atom
    //  is constrained or not.
    int is_atom_constrained(int atomnum) const {
        if (numConstraints) {
            //  Check the index to see if it is constrained
            return(consIndexes[atomnum] != -1);
        }
        else {
            //  No constraints at all, so just return 0
            return(0);
        }
    }

    //  Get the harmonic constraints for a specific atom
    void get_cons_params(double &k, Vector &refPos, int atomnum) const {
        k = consParams[consIndexes[atomnum]].k;
        refPos = consParams[consIndexes[atomnum]].refPos;
    }

    double langevin_param(int atomnum) const {
        return(langevinParams[atomnum]);
    }

    double langevin_force_val(int atomnum) const {
        return(langForceVals[atomnum]);
    }

    int is_atom_fixed(int atomnum) const {
        return (numFixedAtoms && fixedAtomFlags[atomnum]);
    }

    int is_group_fixed(int atomnum) const {
        return (numFixedAtoms && (fixedAtomFlags[atomnum] == -1));
    }

    // 0 if not rigid or length to parent, for parent refers to H-H length
    double rigid_bond_length(int atomnum) const {
        return(rigidBondLengths[atomnum]);
    }
    
    //  Print out list of atoms
    void print_atoms(Parameters *); 
    //  Print out list of bonds    
    void print_bonds(Parameters *); 
    //  Print out list of exclusions
    void print_exclusions();
};

#endif	/* INITIAL_H */
