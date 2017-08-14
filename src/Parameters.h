/* 
 * Author: Hadi
 *
 */

#ifndef PARAMETERS_H
#define	PARAMETERS_H

class StringList;

#include "strlib.h"
#include "Configure.h"
#include "R2.h"

// Check if these are necessary; probably yes
#include <vector>
#include <string>

//  Define the number of Multiples that a Dihedral or Improper
//  bond can have with the Charm22 parameter set
#define MAX_MULTIPLICITY    4

//  Number of characters maximum allowed for storing atom type names
#define MAX_ATOMTYPE_CHARS      6

using namespace std;

typedef struct bond_val {
    double k;       //  Force constant for the bond
    double x0;  //  Rest distance for the bond
} BondValue;

typedef struct angle_val {
    double k;       //  Force constant for angle
    double theta0;  //  Rest angle for angle
    double k_ub;    //  Urey-Bradley force constant
    double r_ub;    //  Urey-Bradley distance
} AngleValue;

typedef struct four_body_consts {
    double k;       //  Force constant
    int n;      //  Periodicity
    double delta;   //  Phase shift
} FourBodyConsts;

typedef struct dihedral_val {
    int multiplicity;
    FourBodyConsts values[MAX_MULTIPLICITY];
} DihedralValue;

typedef struct improper_val {
    int multiplicity;
    FourBodyConsts values[MAX_MULTIPLICITY];
} ImproperValue;

typedef struct nonbondedexcl_val {
    // need to put parameters here...
    // for now, copy bond
    double k;       //  Force constant for the bond
    double x0;  //  Rest distance for the bond
} NonbondedExclValue;

typedef struct vdw_val {
    double sigma;   //  Sigma value
    double epsilon; //  Epsilon value
    double sigma14; //  Sigma value for 1-4 interactions
    double epsilon14; //  Epsilon value for 1-4 interactions
} VdwValue;

//  IndexedVdwPair is used to form a binary search tree that is
//  indexed by vwd_type index.  This is the tree that will be
//  used to search during the actual simulation

typedef struct indexed_vdw_pair {
    Index ind1;     //  Index for first atom type
    Index ind2;     //  Index for second atom type
    double A;           //  Parameter A for this pair
    double A14;     //  Parameter A for 1-4 interactions
    double B;           //  Parameter B for this pair
    double B14;     //  Parameter B for 1-4 interactions
    struct indexed_vdw_pair *right;  //  Right child
    struct indexed_vdw_pair *left;   //  Left child
} IndexedVdwPair;

//  Structures that are defined in Parameters.cpp
struct bond_params;
struct angle_params;
struct improper_params;
struct dihedral_params;
struct vdw_params;
struct vdw_pair_params;

class Parameters {
private:
    char *atomTypeNames; //  Names of atom types
    int AllFilesRead;    //  Flag 1 imples that all
                         //  of the parameter files
                         //  have been read in and
                         //  the arrays have been
                         //  created.
    struct bond_params *bondp;          //  Binary tree of bond params
    struct angle_params *anglep;        //  Binary tree of angle params
    struct improper_params *improperp;  //  Linked list of improper par.
    struct dihedral_params *dihedralp;  //  Linked list of dihedral par.
    struct vdw_params *vdwp;            //  Binary tree of vdw params
    struct vdw_pair_params *vdw_pairp;  //  Binary tree of vdw pairs
    
    BondValue *bond_array;          //  Array of bond params
    AngleValue *angle_array;        //  Array of angle params
    DihedralValue *dihedral_array;  //  Array of dihedral params
    ImproperValue *improper_array;  //  Array of improper params
    VdwValue *vdw_array;            //  Array of vdw params
    IndexedVdwPair *vdw_pair_tree;  //  Tree of vdw pair params
    
    int NumBondParams;          //  Number of bond parameters
    int NumAngleParams;         //  Number of angle parameters
    int NumDihedralParams;      //  Number of dihedral params
    int NumImproperParams;      //  Number of improper params
    int NumVdwParams;           //  Number of vdw parameters
    int NumVdwParamsAssigned;   //  Number actually assigned
    int NumVdwPairParams;       //  Number of vdw_pair params

    int *maxDihedralMults;      //  Max multiplicity for
                                //  dihedral bonds
    int *maxImproperMults;      //  Max multiplicity for
                                //  improper bonds

    void add_bond_param(char *);              //  Add a bond parameter
    struct bond_params *add_to_bond_tree(struct bond_params * , 
                     struct bond_params *);

    void add_angle_param(char *);            //  Add an angle parameter
    struct angle_params *add_to_angle_tree(struct angle_params * , 
                     struct angle_params *);

    void add_dihedral_param(char *, FILE *); //  Add a dihedral parameter
    void add_to_dihedral_list(struct dihedral_params *);
    void add_to_charmm_dihedral_list(struct dihedral_params *);

    void add_improper_param(char *, FILE *); //  Add an improper parameter
    void add_to_improper_list(struct improper_params *);

    void add_vdw_param(char *);             //  Add a vdw parameter
    struct vdw_params *add_to_vdw_tree(struct vdw_params *, 
                     struct vdw_params *);

    void add_vdw_pair_param(char *);        //  Add a vdw pair parameter
    void add_to_vdw_pair_list(struct vdw_pair_params *);

    //  The index_* routines are used to index each of 
    //  the parameters and build the arrays that will be used
    //  for constant time access
    Index index_bonds(struct bond_params *, Index);
    Index index_angles(struct angle_params *, Index);
    Index index_vdw(struct vdw_params *, Index);
    void index_dihedrals();
    void index_impropers();
    
    void convert_vdw_pairs();
    IndexedVdwPair *add_to_indexed_vdw_pairs(IndexedVdwPair *, IndexedVdwPair *);
    
    int vdw_pair_to_arrays(int *, int *, double *, double *, double *, double *, 
                   int, IndexedVdwPair *);

    //  The free_* routines are used by the destructor to deallocate
    //  memory
    void free_bond_tree(struct bond_params *);
    void free_angle_tree(struct angle_params *);
    void free_dihedral_list(struct dihedral_params *);
    void free_improper_list(struct improper_params *);
    void free_vdw_tree(struct vdw_params *);
    void free_vdw_pair_tree(IndexedVdwPair *);
    void free_vdw_pair_list();

    //check if these are necessary; probably not
    vector<string> bond_str;
    vector<string> angle_str;
    vector<string> vdw_str;
    vector<string> dpd_str;

public:

    Parameters(const Configure *conf);
    Parameters(const Parameters& orig);
    virtual ~Parameters();
    // return a string for the Nth atom type.  This can only be
    // called after all the param files have been read and the type
    // names have been indexed.  The Nth atom type refers to the same
    // index of the Nth vdw parameter (i.e. there are NumVdwParams names).
    char *atom_type_name(Index a) {
        return (atomTypeNames + (a * (MAX_ATOMTYPE_CHARS + 1)));
    }

    //  Read a parameter file
    void read_parameter_file(const char *);

    //  Signal the parameter object that all of
    //  the parameter files have been read in
    void done_reading_files();

    //  Signal the parameter object that the
    //  structure file has been read in
    void done_reading_structure();

    //  The assign_*_index routines are used to assign
    //  an index to atoms or bonds.  If an specific atom
    //  or bond type can't be found, then the program 
    //  terminates
    //  Assign a vdw index to an atom
    void assign_vdw_index(char *, Atom *);  
    //  Assign a bond index to a bond
    void assign_bond_index(char *, char *, Bond *); 
    //  Assign an angle index to an angle
    void assign_angle_index(char *, char *, char *, Angle *);
    //  Assign a dihedral index to a dihedral
    void assign_dihedral_index(char *, char*, char*, char *, Dihedral *, int);
    //  Assign an improper index to an improper
    void assign_improper_index(char *, char*, char*, char *, Improper *, int);

    //  The get_*_params routines are the routines that really
    //  do all the work for this object.  Given an index, they
    //  access the parameters and return the relevant information
    void get_bond_params(double *k, double *x0, Index index) const {
        *k = bond_array[index].k;
        *x0 = bond_array[index].x0;
    }

    void get_angle_params(double *k, double *theta0, double *k_ub, double *r_ub,
                  Index index) const {
        *k = angle_array[index].k;
        *theta0 = angle_array[index].theta0;
        *k_ub = angle_array[index].k_ub;
        *r_ub = angle_array[index].r_ub;
    }

    int get_improper_multiplicity(Index index) const {
        return(improper_array[index].multiplicity);
    }

    int get_dihedral_multiplicity(Index index) const {
        return(dihedral_array[index].multiplicity);
    }

    void get_improper_params(double *k, int *n, double *delta, 
                 Index index, int mult) const {
        if ( (mult<0) || (mult>MAX_MULTIPLICITY) ) {
            R2_die("Bad mult index in Parameters::get_improper_params");
        }

        *k = improper_array[index].values[mult].k;
        *n = improper_array[index].values[mult].n;
        *delta = improper_array[index].values[mult].delta;
    }

    void get_dihedral_params(double *k, int *n, double *delta, 
                 Index index, int mult) const {
        if ( (mult<0) || (mult>MAX_MULTIPLICITY) ) {
            R2_die("Bad mult index in Parameters::get_dihedral_params");
        }

        *k = dihedral_array[index].values[mult].k;
        *n = dihedral_array[index].values[mult].n;
        *delta = dihedral_array[index].values[mult].delta;
    }

    void get_vdw_params(double *sigma, double *epsilon, double *sigma14, 
                double *epsilon14, Index index) const {
        *sigma = vdw_array[index].sigma;
        *epsilon = vdw_array[index].epsilon;
        *sigma14 = vdw_array[index].sigma14;
        *epsilon14 = vdw_array[index].epsilon14;
    }

    //  Find a vwd_pair parameter
    int get_vdw_pair_params(Index ind1, Index ind2, double *, double *, double *, double *) const;

    int get_num_vdw_params() const { return NumVdwParamsAssigned; }
};

#endif	/* PARAMETERS_H */

