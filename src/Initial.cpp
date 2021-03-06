#include "Initial.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <algorithm>
#include "strlib.h"
#include "Parameters.h"
#include "PDB.h"
#include "Hydrogen.h"
#include "UniqueSetIter.h"
#include "Configure.h"
#include "Initial.h"

#define UnknownAtom      0x00
#define HydrogenAtom     0x01
#define OxygenAtom       0x02

// Check if these are necessary
#include <iostream>
#include <vector>
#include <algorithm>

Initial::Initial(const Configure *conf, Parameters *param) {

    atoms=NULL;
    atomNames=NULL;
    bonds=NULL;
    angles=NULL;
    dihedrals=NULL;
    impropers=NULL;
    exclusions=NULL;
    bondsWithAtom=NULL;
    bondsByAtom=NULL;
    anglesByAtom=NULL;
    dihedralsByAtom=NULL;
    impropersByAtom=NULL;
    exclusionsByAtom=NULL;
    all_exclusions=NULL;
    onefour_exclusions=NULL;
    langevinParams=NULL;
    langForceVals=NULL;
    fixedAtomFlags=NULL;
    rigidBondLengths=NULL;
    consIndexes=NULL;
    consParams=NULL;

    /*  Initialize counts to 0 */
    numAtoms=0;
    numBonds=0;
    numAngles=0;
    numDihedrals=0;
    numImpropers=0;
    numExclusions=0;
    numConstraints=0;
    numFixedAtoms=0;
    numRigidBonds=0;
    numFixedRigidBonds=0;
    numMultipleDihedrals=0;
    numMultipleImpropers=0;
    numCalcBonds=0;
    numCalcAngles=0;
    numCalcDihedrals=0;
    numCalcImpropers=0;
    numCalcExclusions=0;

    if (param != NULL && conf->psfname != NULL) {
        read_psf_file(conf->psfname, param);
    }
    build_lists_by_atom(); 
}    

Initial::~Initial() {

    delete [] atoms;
    delete [] atomNames;
    delete [] bonds;
    delete [] angles;
    delete [] dihedrals;
    delete [] impropers;
    delete [] exclusions;
    delete [] bondsByAtom;
    delete [] anglesByAtom;
    delete [] dihedralsByAtom;
    delete [] impropersByAtom;
    delete [] exclusionsByAtom;
    delete [] all_exclusions;
    delete [] onefour_exclusions;
    delete [] fixedAtomFlags;
    delete [] rigidBondLengths;
}

void Initial::read_psf_file(const char *fname, Parameters *params) {

    char err_msg[512];  //  Error message for R2_die
    char buffer[512];  //  Buffer for file reading
    int i;      //  Loop counter
    int NumTitle;    //  Number of Title lines in .psf file
    FILE *psf_file;    //  pointer to .psf file
    int ret_code;    //  ret_code from R2_read_line calls

    /* Try and open the .psf file           */
    if ( (psf_file = fopen(fname, "r")) == NULL) {
        sprintf(err_msg, "UNABLE TO OPEN .psf FILE %s", fname);
        R2_die(err_msg);
    }

    /*  Read till we have the first non-blank line of file    */
    ret_code = R2_read_line(psf_file, buffer);

    while ( (ret_code==0) && (R2_blank_string(buffer)) ) {
        ret_code = R2_read_line(psf_file, buffer);
    }

    /*  Check to see if we dropped out of the loop because of a     */
    /*  read error.  This shouldn't happen unless the file is empty */
    if (ret_code!=0) {
        sprintf(err_msg, "EMPTY .psf FILE %s", fname);
        R2_die(err_msg);
    }

    /*  The first non-blank line should contain the word "psf".    */
    /*  If we can't find it, die.               */
    if (!R2_find_word(buffer, "psf")) {
        sprintf(err_msg, "UNABLE TO FIND \"PSF\" STRING IN PSF FILE %s",
        fname);
        R2_die(err_msg);
    }

    /*  Read until we find the next non-blank line      */
    ret_code = R2_read_line(psf_file, buffer);

    while ( (ret_code==0) && (R2_blank_string(buffer)) ) {
        ret_code = R2_read_line(psf_file, buffer);
    }

    /*  Check to see if we dropped out of the loop because of a     */
    /*  read error.  This shouldn't happen unless there is nothing  */
    /*  but the PSF line in the file        */
    if (ret_code!=0) {
        sprintf(err_msg, "MISSING EVERYTHING BUT PSF FROM %s", fname);
        R2_die(err_msg);
    }

    /*  This line should have the word "NTITLE" in it specifying    */
    /*  how many title lines there are        */
    if (!R2_find_word(buffer, "NTITLE")) {
        sprintf(err_msg,"CAN NOT FIND \"NTITLE\" STRING IN PSF FILE %s",
        fname);
        R2_die(err_msg);
    }

    sscanf(buffer, "%d", &NumTitle);

    /*  Now skip the next NTITLE non-blank lines and then read in the*/
    /*  line which should contain NATOM        */
    i=0;

    while ( ((ret_code = R2_read_line(psf_file, buffer)) == 0) && 
        (i<NumTitle) ) {
        if (!R2_blank_string(buffer))
            i++;
    }

    /*  Make sure we didn't exit because of a read error    */
    if (ret_code != 0) {
        sprintf(err_msg, "FOUND EOF INSTEAD OF NATOM IN PSF FILE %s", 
        fname);
        R2_die(err_msg);
    }

    while (R2_blank_string(buffer)) {
        R2_read_line(psf_file, buffer);
    }

    /*  Check to make sure we have the line we want      */
    if (!R2_find_word(buffer, "NATOM")) {
        sprintf(err_msg, "DIDN'T FIND \"NATOM\" IN PSF FILE %s",
        fname);
        R2_die(err_msg);
    }

    /*  Read in the number of atoms, and then the atoms themselves  */
    sscanf(buffer, "%d", &numAtoms);

    read_atoms(psf_file, params);

    /*  Read until we find the next non-blank line      */
    ret_code = R2_read_line(psf_file, buffer);

    while ( (ret_code==0) && (R2_blank_string(buffer)) ) {
        ret_code = R2_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0) {
        R2_die("EOF ENCOUNTERED LOOKING FOR NBONDS IN PSF");
    }

    /*  Look for the string "NBOND"          */
    if (!R2_find_word(buffer, "NBOND")) {
        R2_die("DID NOT FIND NBOND AFTER ATOM LIST IN PSF");
    }

    /*  Read in the number of bonds and then the bonds themselves  */
    sscanf(buffer, "%d", &numBonds);

    if (numBonds)
        read_bonds(psf_file, params);

    /*  Read until we find the next non-blank line      */
    ret_code = R2_read_line(psf_file, buffer);

    while ( (ret_code==0) && (R2_blank_string(buffer)) ) {
        ret_code = R2_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0) {
        R2_die("EOF ENCOUNTERED LOOKING FOR NTHETA IN PSF");
    }

    /*  Look for the string "NTHETA"        */
    if (!R2_find_word(buffer, "NTHETA")) {
        R2_die("DID NOT FIND NTHETA AFTER BOND LIST IN PSF");
    }

    /*  Read in the number of angles and then the angles themselves */
    sscanf(buffer, "%d", &numAngles);

    if (numAngles)
        read_angles(psf_file, params);

    /*  Read until we find the next non-blank line      */
    ret_code = R2_read_line(psf_file, buffer);

    while ( (ret_code==0) && (R2_blank_string(buffer)) ) {
        ret_code = R2_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0) {
        R2_die("EOF ENCOUNTERED LOOKING FOR NPHI IN PSF");
    }

    /*  Look for the string "NPHI"          */
    if (!R2_find_word(buffer, "NPHI")) {
        R2_die("DID NOT FIND NPHI AFTER ANGLE LIST IN PSF");
    }

    /*  Read in the number of dihedrals and then the dihedrals      */
    sscanf(buffer, "%d", &numDihedrals);

    if (numDihedrals)
        read_dihedrals(psf_file, params);

    /*  Read until we find the next non-blank line      */
      ret_code = R2_read_line(psf_file, buffer);

    while ( (ret_code==0) && (R2_blank_string(buffer)) ) {
        ret_code = R2_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0) {
        R2_die("EOF ENCOUNTERED LOOKING FOR NIMPHI IN PSF");
    }

    /*  Look for the string "NIMPHI"        */
    if (!R2_find_word(buffer, "NIMPHI")) {
        R2_die("DID NOT FIND NIMPHI AFTER ATOM LIST IN PSF");
    }

    /*  Read in the number of Impropers and then the impropers  */
    sscanf(buffer, "%d", &numImpropers);

    if (numImpropers)
        read_impropers(psf_file, params);

    /*  Read until we find the next non-blank line      */
    ret_code = R2_read_line(psf_file, buffer);

    while ( (ret_code==0) && (R2_blank_string(buffer)) ) {
        ret_code = R2_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0) {
        R2_die("EOF ENCOUNTERED LOOKING FOR NDON IN PSF");
    }

    /*  Look for the string "NDON"        */
    if (!R2_find_word(buffer, "NDON")) {
        R2_die("DID NOT FIND NDON AFTER ATOM LIST IN PSF");
    }

    /*  Read in the number of hydrogen bond donors and then the donors */
    sscanf(buffer, "%d", &numDonors);

    if (numDonors)
        read_donors(psf_file);

    /*  Read until we find the next non-blank line      */
    ret_code = R2_read_line(psf_file, buffer);

    while ( (ret_code==0) && (R2_blank_string(buffer)) ) {
        ret_code = R2_read_line(psf_file, buffer);
    }

    /*  Check to make sure we didn't hit the EOF      */
    if (ret_code != 0) {
        R2_die("EOF ENCOUNTERED LOOKING FOR NACC IN PSF");
    }

    /*  Look for the string "NACC"        */
    if (!R2_find_word(buffer, "NACC")) {
        R2_die("DID NOT FIND NACC AFTER ATOM LIST IN PSF");
    }

    /*  Read in the number of hydrogen bond donors and then the donors */
    sscanf(buffer, "%d", &numAcceptors);

    if (numAcceptors)
        read_acceptors(psf_file);

    /*  look for the explicit non-bonded exclusion section.     */
    while (!R2_find_word(buffer, "NNB")) {
        ret_code = R2_read_line(psf_file, buffer);

        if (ret_code != 0) {
        R2_die("EOF ENCOUNTERED LOOKING FOR NNB IN PSF FILE");
        }
    }

    /*  Read in the number of exclusions and then the exclusions    */
    sscanf(buffer, "%d", &numExclusions);

    if (numExclusions)
        read_exclusions(psf_file);

    /*  Close the .psf file.  There is a Group section in the .psf  */
    /*  file after the NNB section, but currently, R2 does not    */
    /*  use this section for anything.        */
    fclose(psf_file);

    //  analyze the data and find the status of each atom
    build_atom_status();

    return;
}
/*      END OF FUNCTION read_psf_file      */

/************************************************************************/
/*                  */
/*        FUNCTION read_atoms      */
/*                  */
/*   INPUTS:                */
/*  fd - file pointer to the .psf file        */
/*  params - Parameters object to use for parameters    */
/*                  */
/*  this function reads in the Atoms section of the .psf file.      */
/*   This section consists of numAtoms lines that are of the form:  */
/*     <atom#> <mol> <seg#> <res> <atomname> <atomtype> <charge> <mass> */
/*   Each line is read into the appropriate entry in the atoms array.   */
/*   The parameters object is then used to determine the vdW constants  */
/*   for this atom.              */
/*                  */
/************************************************************************/

void Initial::read_atoms(FILE *fd, Parameters *params) {
    
    char buffer[512];  // Buffer for reading from file
    int atom_number=0;  // Atom number 
    int last_atom_number=0; // Last atom number, used to assure
            // atoms are in order
    char segment_name[11]; // Segment name
    int residue_number; // Residue number
    char residue_name[11];  // Residue name
    char atom_name[11];  // Atom name
    char atom_type[11];  // Atom type
    double charge;    // Charge for the current atom
    double mass;    // Mass for the current atom
    int read_count;    // Number of fields read by sscanf

    /*  Allocate the atom arrays          */
    atoms     = new Atom[numAtoms];
    atomNames = new AtomNameInfo[numAtoms];
    hydrogenGroup.resize(0);

    if (atoms == NULL || atomNames == NULL ) {
        R2_die("memory allocation failed in Initial::read_atoms");
    }

    /*  Loop and read in numAtoms atom lines.      */
    while (atom_number < numAtoms) {
        /*  Get the line from the file        */
        R2_read_line(fd, buffer);

        /*  If its blank or a comment, skip it      */
        if ( (R2_blank_string(buffer)) || (buffer[0] == '!') )
            continue;

        /*  Parse up the line          */
        read_count=sscanf(buffer, "%d %s %i %s %s %s %lf %lf",
            &atom_number, segment_name, &residue_number,
            residue_name, atom_name, atom_type, &charge, &mass);

        /*  Check to make sure we found what we were expecting  */
        if (read_count != 8) {
            char err_msg[128];

            sprintf(err_msg, "BAD ATOM LINE FORMAT IN PSF FILE IN ATOM LINE %d\nLINE=%s",
            last_atom_number+1, buffer);
            R2_die(err_msg);
        }

        /*  Make sure the atoms were in sequence    */
        if (atom_number != last_atom_number+1) {
            char err_msg[128];

            sprintf(err_msg, "ATOM NUMBERS OUT OF ORDER AT ATOM #%d OF PSF FILE",
            last_atom_number+1);
            R2_die(err_msg);
        }

        last_atom_number++;

        /*  Dynamically allocate strings for atom name, atom    */
        /*  type, etc so that we only allocate as much space    */
        /*  for these strings as we really need      */
        int reslength = strlen(residue_name)+1;
        int seglength = strlen(segment_name)+1;
        int namelength = strlen(atom_name)+1;
        int typelength = strlen(atom_type)+1;

        atomNames[atom_number-1].resname = new char[reslength];
        atomNames[atom_number-1].segname = new char[seglength];
        atomNames[atom_number-1].atomname = new char[namelength];
        atomNames[atom_number-1].atomtype = new char[typelength];
  
        if (atomNames[atom_number-1].resname == NULL)
            R2_die("memory allocation failed in Initial::read_atoms");

        /*  Put the values from this atom into the atoms array  */
        strcpy(atomNames[atom_number-1].resname, residue_name);
        strcpy(atomNames[atom_number-1].segname, segment_name);
        strcpy(atomNames[atom_number-1].atomname, atom_name);
        strcpy(atomNames[atom_number-1].atomtype, atom_type);
        atomNames[atom_number-1].resnum = residue_number;

        atoms[atom_number-1].mass = mass;
        atoms[atom_number-1].charge = charge;
        atoms[atom_number-1].status = UnknownAtom;

        /*  Determine the type of the atom (H or O) */
        if (atoms[atom_number-1].mass <=3.5) {
            atoms[atom_number-1].status |= HydrogenAtom;
        } 
        else if ((atomNames[atom_number-1].atomname[0] == 'O') && 
             (atoms[atom_number-1].mass >= 14.0) && 
            (atoms[atom_number-1].mass <= 18.0)) {
            atoms[atom_number-1].status |= OxygenAtom;
        }

        /*  Look up the vdw constants for this atom    */
        params->assign_vdw_index(atomNames[atom_number-1].atomtype, 
            &(atoms[atom_number-1]));
    }

    return;
}
/*      END OF FUNCTION read_atoms      */

/************************************************************************/
/*                  */
/*      FUNCTION read_bonds        */
/*                  */
/*  read_bonds reads in the bond section of the .psf file.  This    */
/*  section contains a list of pairs of numbers where each pair is      */
/*  represents two atoms that are bonded together.  Each atom pair is   */
/*  read in.  Then that parameter object is queried to determine the    */
/*  force constant and rest distance for the bond.      */
/*                  */
/************************************************************************/

void Initial::read_bonds(FILE *fd, Parameters *params) {
    
    int atom_nums[2];  // Atom indexes for the bonded atoms
    char atom1name[11];  // Atom type for atom #1
    char atom2name[11];  // Atom type for atom #2
    register int j;      // Loop counter
    int num_read=0;    // Number of bonds read so far
    int origNumBonds = numBonds;   // number of bonds in file header

    /*  Allocate the array to hold the bonds      */
    bonds=new Bond[numBonds];

    if (bonds == NULL) {
        R2_die("memory allocations failed in Initial::read_bonds");
    }

    /*  Loop through and read in all the bonds      */
    while (num_read < numBonds) {
        /*  Loop and read in the two atom indexes    */
        for (j=0; j<2; j++) {
            /*  Read the atom number from the file.         */
            /*  Subtract 1 to convert the index from the    */
            /*  1 to NumAtoms used in the file to the       */
            /*  0 to NumAtoms-1 that we need    */
            atom_nums[j]=R2_read_int(fd, "BONDS")-1;

            /*  Check to make sure the index isn't too big  */
            if (atom_nums[j] >= numAtoms) {
                char err_msg[128];

                sprintf(err_msg, "BOND INDEX %d GREATER THAN NATOM %d IN BOND # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
                R2_die(err_msg);
            }
        }

        /*  Get the atom type for the two atoms.  When we query */
        /*  the parameter object, we need to send the atom type */
        /*  that is alphabetically first as atom 1.    */
        if (strcasecmp(atomNames[atom_nums[0]].atomtype, 
             atomNames[atom_nums[1]].atomtype) < 0) {
            strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
            strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
        }
        else {
            strcpy(atom2name, atomNames[atom_nums[0]].atomtype);
            strcpy(atom1name, atomNames[atom_nums[1]].atomtype);
        }

        /*  Query the parameter object for the constants for    */
        /*  this bond            */
        Bond *b = &(bonds[num_read]);
        params->assign_bond_index(atom1name, atom2name, b);

        /*  Assign the atom indexes to the array element  */
        b->atom1=atom_nums[0];
        b->atom2=atom_nums[1];

        /*  Make sure this isn't a fake bond meant for shake in x-plor.  */
        double k, x0;
        params->get_bond_params(&k,&x0,b->bond_type);
        if ( k == 0. ) --numBonds;  // fake bond
            else ++num_read;  // real bond
    }

    /*  Tell user about our subterfuge  */
    if ( numBonds != origNumBonds ) {
        cout <<   "Ignored " << origNumBonds - numBonds <<
                " bonds with zero force constants.\n" << endl;
        cout <<  
        "Will get H-H distance in rigid H2O from H-O-H angle.\n" << endl;
    }

    return;
}
/*      END OF FUNCTION read_bonds      */

/************************************************************************/
/*                  */
/*      FUNCTION read_angles        */
/*                  */
/*   INPUTS:                */
/*  fd - File descriptor for .psf file        */
/*  params - Parameters object to query for parameters    */
/*                  */
/*  read_angles reads the angle parameters from the .psf file.      */
/*   This section of the .psf file consists of a list of triplets of    */
/*   atom indexes.  Each triplet represents three atoms connected via   */
/*   an angle bond.  The parameter object is queried to obtain the      */
/*   constants for each bond.            */
/*                  */
/************************************************************************/

void Initial::read_angles(FILE *fd, Parameters *params) {

    int atom_nums[3];  //  Atom numbers for the three atoms
    char atom1name[11];  //  Atom type for atom 1
    char atom2name[11];  //  Atom type for atom 2
    char atom3name[11];  //  Atom type for atom 3
    register int j;      //  Loop counter
    int num_read=0;    //  Number of angles read so far
    int origNumAngles = numAngles;  // Number of angles in file
    /*  Alloc the array of angles          */
    angles=new Angle[numAngles];

    if (angles == NULL) {
        R2_die("memory allocation failed in Initial::read_angles");
    }

    /*  Loop through and read all the angles      */
    while (num_read < numAngles) {
        /*  Loop through the 3 atom indexes in the current angle*/
        for (j=0; j<3; j++) {
            /*  Read the atom number from the file.         */
            /*  Subtract 1 to convert the index from the    */
            /*  1 to NumAtoms used in the file to the       */
            /*  0 to NumAtoms-1 that we need    */
            atom_nums[j]=R2_read_int(fd, "ANGLES")-1;

            /*  Check to make sure the atom index doesn't   */
            /*  exceed the Number of Atoms      */
            if (atom_nums[j] >= numAtoms) {
                char err_msg[128];

                sprintf(err_msg, "ANGLES INDEX %d GREATER THAN NATOM %d IN ANGLES # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
                R2_die(err_msg);
            }
        }

        /*  Place the bond name that is alphabetically first  */
        /*  in the atom1name.  This is OK since the order of    */
        /*  atom1 and atom3 are interchangable.  And to search  */
        /*  the tree of angle parameters, we need the order     */
        /*  to be predictable.          */
        if (strcasecmp(atomNames[atom_nums[0]].atomtype, 
             atomNames[atom_nums[2]].atomtype) < 0) {
            strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
            strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
            strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
        }
        else {
            strcpy(atom1name, atomNames[atom_nums[2]].atomtype);
            strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
            strcpy(atom3name, atomNames[atom_nums[0]].atomtype);
        }

        /*  Get the constant values for this bond from the  */
        /*  parameter object          */
        params->assign_angle_index(atom1name, atom2name, 
        atom3name, &(angles[num_read]));

        /*  Assign the three atom indices      */
        angles[num_read].atom1=atom_nums[0];
        angles[num_read].atom2=atom_nums[1];
        angles[num_read].atom3=atom_nums[2];

        /*  Make sure this isn't a fake angle meant for shake in x-plor.  */
        double k, t0, k_ub, r_ub;
        params->get_angle_params(&k,&t0,&k_ub,&r_ub,angles[num_read].angle_type);
        if ( k == 0. && k_ub == 0. ) --numAngles;  // fake angle
        else ++num_read;  // real angle
    }

    /*  Tell user about our subterfuge  */
    if ( numAngles != origNumAngles ) {
        cout <<   "Ignored " << origNumAngles - numAngles <<
                " angles with zero force constants.\n" << endl;
    }

    return;
}
/*      END OF FUNCTION read_angles      */

/************************************************************************/
/*                  */
/*        FUNCTION read_dihedrals      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for the .psf file        */
/*  params - pointer to parameter object        */
/*                  */
/*  read_dihedreals reads the dihedral section of the .psf file.    */
/*   This section of the file contains a list of quartets of atom       */
/*   numbers.  Each quartet represents a group of atoms that form a     */
/*   dihedral bond.              */
/*                  */
/************************************************************************/

void Initial::read_dihedrals(FILE *fd, Parameters *params) {

    int atom_nums[4];  // The 4 atom indexes
    int last_atom_nums[4];  // Atom numbers from previous bond
    char atom1name[11];  // Atom type for atom 1
    char atom2name[11];  // Atom type for atom 2
    char atom3name[11];  // Atom type for atom 3
    char atom4name[11];  // Atom type for atom 4
    register int j;      // loop counter
    int num_read=0;    // number of dihedrals read so far
    int multiplicity=1;  // multiplicity of the current bond
    int duplicate_bond;  // Is this a duplicate of the last bond
    int num_unique=0;   // Number of unique dihedral bonds

    //  Initialize the array used to check for duplicate dihedrals
    for (j=0; j<4; j++)
        last_atom_nums[j] = -1;

    /*  Allocate an array to hold the Dihedrals      */
    dihedrals = new Dihedral[numDihedrals];

    if (dihedrals == NULL) {
        R2_die("memory allocation failed in Initial::read_dihedrals");
    }

    /*  Loop through and read all the dihedrals      */
    while (num_read < numDihedrals) {
        duplicate_bond = 1;

        /*  Loop through and read the 4 indexes for this bond   */
        for (j=0; j<4; j++) {
            /*  Read the atom number from the file.         */
            /*  Subtract 1 to convert the index from the    */
            /*  1 to NumAtoms used in the file to the       */
            /*  0 to NumAtoms-1 that we need    */
            atom_nums[j]=R2_read_int(fd, "DIHEDRALS")-1;

            /*  Check for an atom index that is too large  */
            if (atom_nums[j] >= numAtoms) {
                char err_msg[128];

                sprintf(err_msg, "DIHEDRALS INDEX %d GREATER THAN NATOM %d IN DIHEDRALS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
                R2_die(err_msg);
            }

            //  Check to see if this atom matches the last bond
            if (atom_nums[j] != last_atom_nums[j]) {
                duplicate_bond = 0;
            }

            last_atom_nums[j] = atom_nums[j];
        }

        /*  Get the atom types for the 4 atoms so we can look  */
        /*  up the constants in the parameter object    */
        strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
        strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
        strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
        strcpy(atom4name, atomNames[atom_nums[3]].atomtype);

        //  Check to see if this is really a new bond or just
        //  a repeat of the last one
        if (duplicate_bond) {
            //  This is a duplicate, so increase the multiplicity
            multiplicity++;

            if (multiplicity == 2) {
                numMultipleDihedrals++;
            }
        }
        else {
            multiplicity=1;
            num_unique++;
        }

        /*  Get the constants for this dihedral bond    */
        params->assign_dihedral_index(atom1name, atom2name, 
        atom3name, atom4name, &(dihedrals[num_unique-1]),
            multiplicity);

        /*  Assign the atom indexes        */
        dihedrals[num_unique-1].atom1=atom_nums[0];
        dihedrals[num_unique-1].atom2=atom_nums[1];
        dihedrals[num_unique-1].atom3=atom_nums[2];
        dihedrals[num_unique-1].atom4=atom_nums[3];

        num_read++;
    }

    numDihedrals = num_unique;

    return;
}
/*      END OF FUNCTION read_dihedral      */

/************************************************************************/
/*                  */
/*        FUNCTION read_impropers      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for .psf file        */
/*  params - parameter object          */
/*                  */
/*  read_impropers reads the improper section of the .psf file.  */
/*   This section is identical to the dihedral section in that it is    */
/*   made up of a list of quartets of atom indexes that define the      */
/*   atoms that are bonded together.          */
/*                  */
/************************************************************************/

void Initial::read_impropers(FILE *fd, Parameters *params) {
    
    int atom_nums[4];  //  Atom indexes for the 4 atoms
    int last_atom_nums[4];  //  Atom indexes from previous bond
    char atom1name[11];  //  Atom type for atom 1
    char atom2name[11];  //  Atom type for atom 2
    char atom3name[11];  //  Atom type for atom 3
    char atom4name[11];  //  Atom type for atom 4
    register int j;      //  Loop counter
    int num_read=0;    //  Number of impropers read so far
    int multiplicity=1;  // multiplicity of the current bond
    int duplicate_bond;  // Is this a duplicate of the last bond
    int num_unique=0;   // Number of unique dihedral bonds

    //  Initialize the array used to look for duplicate improper
    //  entries.  Set them all to -1 so we know nothing will match
    for (j=0; j<4; j++)
        last_atom_nums[j] = -1;

    /*  Allocate the array to hold the impropers      */
    impropers=new Improper[numImpropers];

    if (impropers == NULL) {
        R2_die("memory allocation failed in Initial::read_impropers");
    }

    /*  Loop through and read all the impropers      */
    while (num_read < numImpropers) {
        duplicate_bond = 1;

        /*  Loop through the 4 indexes for this improper  */
        for (j=0; j<4; j++) {
            /*  Read the atom number from the file.         */
            /*  Subtract 1 to convert the index from the    */
            /*  1 to NumAtoms used in the file to the       */
            /*  0 to NumAtoms-1 that we need    */
            atom_nums[j]=R2_read_int(fd, "IMPROPERS")-1;

            /*  Check to make sure the index isn't too big  */
            if (atom_nums[j] >= numAtoms) {
                char err_msg[128];

                sprintf(err_msg, "IMPROPERS INDEX %d GREATER THAN NATOM %d IN IMPROPERS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
                R2_die(err_msg);
            }

            if (atom_nums[j] != last_atom_nums[j]) {
                duplicate_bond = 0;
            }

            last_atom_nums[j] = atom_nums[j];
        }

        /*  Get the atom types so we can look up the parameters */
        strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
        strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
        strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
        strcpy(atom4name, atomNames[atom_nums[3]].atomtype);

        //  Check to see if this is a duplicate improper
        if (duplicate_bond) {
            //  This is a duplicate improper.  So we don't
            //  really count this entry, we just update
            //  the parameters object
            multiplicity++;

            if (multiplicity == 2) {
                //  Count the number of multiples.
                numMultipleImpropers++;
            }
        }
        else {
            //  Not a duplicate
            multiplicity = 1;
            num_unique++;
        }

        /*  Look up the constants for this bond      */
        params->assign_improper_index(atom1name, atom2name, 
        atom3name, atom4name, &(impropers[num_unique-1]),
            multiplicity);

        /*  Assign the atom indexes        */
        impropers[num_unique-1].atom1=atom_nums[0];
        impropers[num_unique-1].atom2=atom_nums[1];
        impropers[num_unique-1].atom3=atom_nums[2];
        impropers[num_unique-1].atom4=atom_nums[3];

        num_read++;
    }

    //  Now reset the numImpropers value to the number of UNIQUE
    //  impropers.  Sure, we waste a few entries in the improper_array
    //  on the master node, but it is very little space . . .
    numImpropers = num_unique;

    return;
}
/*      END OF FUNCTION read_impropers      */

// Read the donors, and ignore them.

void Initial::read_donors(FILE *fd) {
  
    for (int i=0; i<numDonors; i++) {
        R2_read_int(fd, "DONORS")-1;
        R2_read_int(fd, "DONORS")-1;
    }
}

void Initial::read_acceptors(FILE *fd) {

    for (int i=0; i<numAcceptors; i++) {
        R2_read_int(fd, "ACCEPTORS")-1;
        R2_read_int(fd, "ACCEPTORS")-1;
    } 
}

/************************************************************************/
/*                  */
/*      FUNCTION read_exclusions      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for .psf file        */
/*                  */
/*  read_exclusions reads in the explicit non-bonded exclusions     */
/*  from the .psf file.  This section is a little funky, so hang on.    */
/*  Ok, first there is a list of atom indexes that is NumExclusions     */
/*  long.  These are in some sense the atoms that will be exlcuded.     */
/*  Following this list is a list of NumAtoms length that is a list     */
/*  of indexes into the list of excluded atoms.  So an example.  Suppose*/
/*  we have a 5 atom simulation with 3 explicit exclusions.  The .psf   */
/*  file could look like:            */
/*                  */
/*  3!NNB                */
/*  3 4 5                */
/*  0 1 3 3 3              */
/*                  */
/*  This would mean that atom 1 has no explicit exclusions.  Atom 2     */
/*  has an explicit exclusion with atom 3.  Atom 3 has an explicit      */
/*  exclusion with atoms 4 AND 5.  And atoms 4 and 5 have no explicit   */
/*  exclusions.  Got it!?!  I'm not sure who dreamed this up . . .      */
/*                  */
/************************************************************************/

void Initial::read_exclusions(FILE *fd) {
    
    int *exclusion_atoms;  //  Array of indexes of excluded atoms
    register int num_read=0;    //  Number fo exclusions read in
    int current_index;  //  Current index value
    int last_index;    //  the previous index value
    register int insert_index=0;  //  index of where we are in exlcusions array

    /*  Allocate the array of exclusion structures and the array of */
    /*  exlcuded atom indexes          */
    exclusions      = new Exclusion[numExclusions];
    exclusion_atoms = new int[numExclusions];

    if ( (exclusions == NULL) || (exclusion_atoms == NULL) ) {
        R2_die("memory allocation failed in Initial::read_exclusions");
    }

    /*  First, read in the excluded atoms list      */
    for (num_read=0; num_read<numExclusions; num_read++) {
        /*  Read the atom number from the file. Subtract 1 to   */
        /*  convert the index from the 1 to NumAtoms used in the*/
        /*  file to the  0 to NumAtoms-1 that we need    */
        exclusion_atoms[num_read]=R2_read_int(fd, "IMPROPERS")-1;

        /*  Check for an illegal index        */
        if (exclusion_atoms[num_read] >= numAtoms) {
            char err_msg[128];

            sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN PSF FILE", exclusion_atoms[num_read]+1, numAtoms, num_read+1);
            R2_die(err_msg);
        }
    }

    /*  Now, go through and read the list of NumAtoms pointers into */
    /*  the array that we just read in        */
    last_index=0;

    for (num_read=0; num_read<numAtoms; num_read++) {
        /*  Read in the current index value      */
        current_index=R2_read_int(fd, "EXCLUSIONS");

        /*  Check for an illegal pointer      */
        if (current_index>numExclusions) {
            char err_msg[128];

            sprintf(err_msg, "EXCLUSION INDEX %d LARGER THAN NUMBER OF EXLCUSIONS %d IN PSF FILE, EXCLUSION #%d\n", 
            current_index+1, numExclusions, num_read);
            R2_die(err_msg);
        }

        /*  Check to see if it matches the last index.  If so   */
        /*  than this atom has no exclusions.  If not, then     */
        /*  we have to build some exclusions      */
        if (current_index != last_index) {
            /*  This atom has some exlcusions.  Loop from   */
            /*  the last_index to the current index.  This  */
            /*  will include how ever many exclusions this  */
            /*  atom has          */
            for (insert_index=last_index; 
                insert_index<current_index; insert_index++) {
                /*  Assign the two atoms involved.      */
                /*  The first one is our position in    */
                /*  the list, the second is based on    */
                /*  the pointer into the index list     */
                exclusions[insert_index].atom1=num_read;
                exclusions[insert_index].atom2=exclusion_atoms[insert_index];
            }

            last_index=current_index;
        }
    }

    /*  Free our temporary list of indexes        */
    delete [] exclusion_atoms;

    return;
}
/*      END OF FUNCTION read_exclusions      */

/************************************************************************/
/*                  */
/*      FUNCTION print_atoms        */
/*                  */
/*  print_atoms prints out the list of atoms stored in this object. */
/*  It is inteded mainly for debugging purposes.      */
/*                  */
/************************************************************************/

void Initial::print_atoms(Parameters *params) {

    register int i;
    double sigma;
    double epsilon;
    double sigma14;
    double epsilon14;

    cout << "ATOM LIST\n" 
          << "******************************************\n" 
                      << "NUM  NAME TYPE RES  MASS    CHARGE SIGMA   EPSILON SIGMA14 EPSILON14\n" 
        << endl;

    for (i=0; i<numAtoms; i++){
        params->get_vdw_params(&sigma, &epsilon, &sigma14, &epsilon14, 
            atoms[i].vdw_type);

        cout << i+1 << " " << atomNames[i].atomname  
            << " " << atomNames[i].atomtype << " " 
            << atomNames[i].resname  << " " << atoms[i].mass  
            << " " << atoms[i].charge << " " << sigma 
            << " " << epsilon << " " << sigma14 
            << " " << epsilon14 
            << endl;
    }
}
/*      END OF FUNCTION print_atoms      */

/************************************************************************/
/*                  */
/*      FUNCTION print_bonds        */
/*                  */
/*  print_bonds prints out the list of bonds stored in this object. */
/*  It is inteded mainly for debugging purposes.      */
/*                  */
/************************************************************************/

void Initial::print_bonds(Parameters *params) {
    register int i;
    double k;
    double x0;

    cout << "BOND LIST\n" << "********************************\n" 
        << "ATOM1 ATOM2 TYPE1 TYPE2      k        x0" 
        << endl;

    for (i=0; i<numBonds; i++) {
        params->get_bond_params(&k, &x0, bonds[i].bond_type);

        cout << bonds[i].atom1+1 << " " 
            << bonds[i].atom2+1 << " "   
            << atomNames[bonds[i].atom1].atomtype << " "  
            << atomNames[bonds[i].atom2].atomtype << " " << k 
            << " " << x0 << endl;
    }
}
/*      END OF FUNCTION print_bonds      */

/************************************************************************/
/*                  */
/*      FUNCTION print_exclusions      */
/*                  */
/*  print_exlcusions prints out the list of exlcusions stored in    */
/*  this object.  It is inteded mainly for debugging purposes.    */
/*                  */
/************************************************************************/

void Initial::print_exclusions() {
    register int i;

    cout << "EXPLICIT EXCLUSION LIST\n" 
        << "********************************\n" 
        << "ATOM1 ATOM2 " 
        << endl;

    for (i=0; i<numExclusions; i++) {

        cout << exclusions[i].atom1+1 << "  " 
            << exclusions[i].atom2+1 << endl;
    }
}
/*      END OF FUNCTION print_exclusions    */


    /************************************************************************/
    /*                  */
    /*      FUNCTION build_lists_by_atom      */
    /*                  */
    /*  This function builds O(NumAtoms) arrays that store the bonds,   */
    /*  angles, dihedrals, and impropers, that each atom is involved in.    */
    /*  This is a space hog, but VERY fast.  This will certainly have to    */
    /*  change to make things scalable in memory, but for now, speed is the */
    /*  thing!                */
    /*                  */
    /************************************************************************/

void Initial::build_lists_by_atom() {

    register int i;      //  Loop counter
    //register int numFixedAtoms = this->numFixedAtoms;  // many tests
    
    bondsWithAtom = new int *[numAtoms];
    bondsByAtom = new int *[numAtoms];
    anglesByAtom = new int *[numAtoms];
    dihedralsByAtom = new int *[numAtoms];
    impropersByAtom = new int *[numAtoms];
    exclusionsByAtom = new int *[numAtoms];

    int *byAtomSize = new int[numAtoms];

    //  Build the bond lists
    for (i=0; i<numAtoms; i++) {
        byAtomSize[i] = 0;
    }

    for (i=0; i<numBonds; i++) {
        byAtomSize[bonds[i].atom1]++;
        byAtomSize[bonds[i].atom2]++;
    }
    for (i=0; i<numAtoms; i++) {
        bondsWithAtom[i] = new int[byAtomSize[i]+1];
        bondsWithAtom[i][byAtomSize[i]] = -1;
        byAtomSize[i] = 0;
    }
    for (i=0; i<numBonds; i++) {
        int a1 = bonds[i].atom1;
        int a2 = bonds[i].atom2;
        bondsWithAtom[a1][byAtomSize[a1]++] = i;
        bondsWithAtom[a2][byAtomSize[a2]++] = i;
    }
    //  Build the bond lists
    for (i=0; i<numAtoms; i++) {
        byAtomSize[i] = 0;
    }
    numCalcBonds = 0;
    for (i=0; i<numBonds; i++) {
        if ( numFixedAtoms && fixedAtomFlags[bonds[i].atom1]
                            && fixedAtomFlags[bonds[i].atom2] ) continue;
         byAtomSize[bonds[i].atom1]++;
         numCalcBonds++;
    }
    for (i=0; i<numAtoms; i++) {
        bondsByAtom[i] = new int[byAtomSize[i]+1];
        bondsByAtom[i][byAtomSize[i]] = -1;
        byAtomSize[i] = 0;
    }
    for (i=0; i<numBonds; i++) {
        if ( numFixedAtoms && fixedAtomFlags[bonds[i].atom1]
                            && fixedAtomFlags[bonds[i].atom2] ) continue;
        int a1 = bonds[i].atom1;
        bondsByAtom[a1][byAtomSize[a1]++] = i;
    }
    //  Build the angle lists
    for (i=0; i<numAtoms; i++) {
        byAtomSize[i] = 0;
    }
    numCalcAngles = 0;
    for (i=0; i<numAngles; i++) {
        if ( numFixedAtoms && fixedAtomFlags[angles[i].atom1]
                            && fixedAtomFlags[angles[i].atom2]
                            && fixedAtomFlags[angles[i].atom3] ) continue;
        byAtomSize[angles[i].atom1]++;
        numCalcAngles++;
    }
    for (i=0; i<numAtoms; i++) {
        anglesByAtom[i] = new int[byAtomSize[i]+1];
        anglesByAtom[i][byAtomSize[i]] = -1;
        byAtomSize[i] = 0;
    }
    for (i=0; i<numAngles; i++) {
        if ( numFixedAtoms && fixedAtomFlags[angles[i].atom1]
                            && fixedAtomFlags[angles[i].atom2]
                            && fixedAtomFlags[angles[i].atom3] ) continue;
        int a1 = angles[i].atom1;
        anglesByAtom[a1][byAtomSize[a1]++] = i;
    }
    //  Build the improper lists
    for (i=0; i<numAtoms; i++) {
        byAtomSize[i] = 0;
    }
    numCalcImpropers = 0;
    for (i=0; i<numImpropers; i++) {
        if ( numFixedAtoms && fixedAtomFlags[impropers[i].atom1]
                            && fixedAtomFlags[impropers[i].atom2]
                            && fixedAtomFlags[impropers[i].atom3]
                            && fixedAtomFlags[impropers[i].atom4] ) continue;
        byAtomSize[impropers[i].atom1]++;
        numCalcImpropers++;
    }
    for (i=0; i<numAtoms; i++) {
        impropersByAtom[i] = new int[byAtomSize[i]+1];
        impropersByAtom[i][byAtomSize[i]] = -1;
        byAtomSize[i] = 0;
    }
    for (i=0; i<numImpropers; i++) {
        if ( numFixedAtoms && fixedAtomFlags[impropers[i].atom1]
                        && fixedAtomFlags[impropers[i].atom2]
                        && fixedAtomFlags[impropers[i].atom3]
                        && fixedAtomFlags[impropers[i].atom4] ) continue;
        int a1 = impropers[i].atom1;
        impropersByAtom[a1][byAtomSize[a1]++] = i;
    }
    //  Build the dihedral lists
    for (i=0; i<numAtoms; i++) {
        byAtomSize[i] = 0;
    }
    numCalcDihedrals = 0;
    for (i=0; i<numDihedrals; i++) {
        if ( numFixedAtoms && fixedAtomFlags[dihedrals[i].atom1]
                        && fixedAtomFlags[dihedrals[i].atom2]
                        && fixedAtomFlags[dihedrals[i].atom3]
                        && fixedAtomFlags[dihedrals[i].atom4] ) continue;
        byAtomSize[dihedrals[i].atom1]++;
        numCalcDihedrals++;
    }
    for (i=0; i<numAtoms; i++) {
        dihedralsByAtom[i] = new int[byAtomSize[i]+1];
        dihedralsByAtom[i][byAtomSize[i]] = -1;
        byAtomSize[i] = 0;
    }
    for (i=0; i<numDihedrals; i++) {
        if ( numFixedAtoms && fixedAtomFlags[dihedrals[i].atom1]
                        && fixedAtomFlags[dihedrals[i].atom2]
                        && fixedAtomFlags[dihedrals[i].atom3]
                        && fixedAtomFlags[dihedrals[i].atom4] ) continue;
        int a1 = dihedrals[i].atom1;
        dihedralsByAtom[a1][byAtomSize[a1]++] = i;
    }
    //  Build the arrays of exclusions for each atom
    build_exclusions();

    //  Remove temporary structures
    delete [] bondsWithAtom;  bondsWithAtom = 0;

    if (exclusions != NULL)
        delete [] exclusions;

    // 1-4 exclusions which are also fully excluded were eliminated by hash table
    numTotalExclusions = exclusionSet.size();
    exclusions = new Exclusion[numTotalExclusions];
    UniqueSetIter<Exclusion> exclIter(exclusionSet);
    for ( exclIter=exclIter.begin(),i=0; exclIter != exclIter.end(); exclIter++,i++ ) {
        exclusions[i] = *exclIter;
    }
    // Free exclusionSet storage
    // exclusionSet.clear(1);
    exclusionSet.clear();

    for (i=0; i<numAtoms; i++) {
        byAtomSize[i] = 0;
    }
    numCalcExclusions = 0;
    for (i=0; i<numTotalExclusions; i++) {
        if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                            && fixedAtomFlags[exclusions[i].atom2] ) continue;
        byAtomSize[exclusions[i].atom1]++;
        numCalcExclusions++;
    }
    for (i=0; i<numAtoms; i++) {
        exclusionsByAtom[i] = new int[byAtomSize[i]+1];
        exclusionsByAtom[i][byAtomSize[i]] = -1;
        byAtomSize[i] = 0;
    }
    for (i=0; i<numTotalExclusions; i++) {
        if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                        && fixedAtomFlags[exclusions[i].atom2] ) continue;
        int a1 = exclusions[i].atom1;
        exclusionsByAtom[a1][byAtomSize[a1]++] = i;
    }

    //  Allocate an array of int *'s to hold the exclusions for
    //  each atom
    all_exclusions = new int *[numAtoms];

    for (i=0; i<numAtoms; i++) {
        byAtomSize[i] = 0;
    }
    for (i=0; i<numTotalExclusions; i++) {
        // first atom should alway have lower number!
        if ( ! exclusions[i].modified ) {
            if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                            && fixedAtomFlags[exclusions[i].atom2] ) continue;
            byAtomSize[exclusions[i].atom1]++;
        }
    }
    for (i=0; i<numAtoms; i++) {
        all_exclusions[i] = new int[byAtomSize[i]+1];
        all_exclusions[i][byAtomSize[i]] = -1;
        byAtomSize[i] = 0;
    }
    for (i=0; i<numTotalExclusions; i++) {
        if ( ! exclusions[i].modified ) {
            if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                               && fixedAtomFlags[exclusions[i].atom2] ) continue;
            int a1 = exclusions[i].atom1;
            int a2 = exclusions[i].atom2;
            all_exclusions[a1][byAtomSize[a1]++] = a2;
        }
    }

    //  If the exclusion policy is scaled 1-4, then allocate
    //  an array of int *'s to hold the 1-4 interactions
    //  Allocate them all the time and assume they are there! -JCP
    //if (simParams->exclude == SCALED14) { 
            {onefour_exclusions = new int *[numAtoms];
            for (i=0; i<numAtoms; i++) {
                byAtomSize[i] = 0;
            }
            for (i=0; i<numTotalExclusions; i++) {
                // first atom should alway have lower number!
                if ( exclusions[i].modified ) {
                    if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                                     && fixedAtomFlags[exclusions[i].atom2] ) continue;
                byAtomSize[exclusions[i].atom1]++;
                }
            }
            for (i=0; i<numAtoms; i++) {
                onefour_exclusions[i] = new int[byAtomSize[i]+1];
                onefour_exclusions[i][byAtomSize[i]] = -1;
                byAtomSize[i] = 0;
            }
            for (i=0; i<numTotalExclusions; i++) {
                if ( exclusions[i].modified ) {
                    if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                                     && fixedAtomFlags[exclusions[i].atom2] ) continue;
                    int a1 = exclusions[i].atom1;
                    int a2 = exclusions[i].atom2;
                    onefour_exclusions[a1][byAtomSize[a1]++] = a2;
                }
            }
    }

    delete [] byAtomSize;

}
/*    END OF FUNCTION build_lists_by_atom    */

void Initial::build_exclusions() {
    
    for (int i=0; i<numExclusions; i++)
    exclusionSet.add(exclusions[i]);
    
    build12excl();
    build13excl();
    build14excl(1);
    //stripHGroupExcl();
}

void Initial::build12excl(void) {
    
    int *current_val;  //  Current value to check
    register int i;    //  Loop counter to loop through all atoms

    //  Loop through all the atoms marking the bonded interactions for each one
    for (i=0; i<numAtoms; i++) {
        current_val = bondsWithAtom[i];

        //  Loop through all the bonds for this atom
        while (*current_val != -1) {
            if (bonds[*current_val].atom1 == i) {
                if (i<bonds[*current_val].atom2) {
                    exclusionSet.add(Exclusion(i,bonds[*current_val].atom2));
                }
            }
            else {
                if (i<bonds[*current_val].atom1) {
                    exclusionSet.add(Exclusion(i,bonds[*current_val].atom1));
                }
            }

            ++current_val;
        }
    }
}

void Initial::build13excl() {

    int *bond1, *bond2;  //  The two bonds being checked
    int middle_atom;  //  Common third atom
    register int i;    //  Loop counter to loop through all atoms
       
    //  Loop through all the atoms looking at the bonded connections
    //  for each one
    for (i=0; i<numAtoms; i++) {
        bond1 = bondsWithAtom[i];
       
        //  Loop through all the bonds directly connect to atom i
        while (*bond1 != -1) {
            if (bonds[*bond1].atom1 == i) {
                middle_atom=bonds[*bond1].atom2;
            }
            else {
                middle_atom=bonds[*bond1].atom1;
            }

            bond2 = bondsWithAtom[middle_atom];

            //  Now loop through all the bonds connect to the
            //  middle atom
            while (*bond2 != -1) {
                if (bonds[*bond2].atom1 == middle_atom) {
                    if (i < bonds[*bond2].atom2) {
                        exclusionSet.add(Exclusion(i,bonds[*bond2].atom2));
                    }
                }
                else {
                    if (i < bonds[*bond2].atom1) {
                        exclusionSet.add(Exclusion(i,bonds[*bond2].atom1));
                    }
                }

                ++bond2;
            }

            ++bond1;
        }
    }
}
/*      END OF FUNCTION build13excl      */

    /************************************************************************/
    /*                  */
    /*        FUNCTION build14excl      */
    /*                  */
    /*   INPUTS:                */
    /*  lists - Array of IntList objects to put exclusions into    */
    /*                  */
    /*  This function calculates all the 1-4 exclusions (that is,  */
    /*   atoms that are connected via a sequence of three linear bonds) and */
    /*   places these interactions into the array of IntList object passed  */
    /*   in.                */
    /*                  */
    /************************************************************************/


void Initial::build14excl(int modified) {

    int *bond1, *bond2, *bond3;  //  The two bonds being checked
    int mid1, mid2;    //  Middle atoms
    register int i;      //  Counter to loop through all atoms
       
    //  Loop through all the atoms
    for (i=0; i<numAtoms; i++) {  
        // Get all the bonds connect directly to atom i
        bond1 = bondsWithAtom[i];

        while (*bond1 != -1) {
            if (bonds[*bond1].atom1 == i) {
                mid1=bonds[*bond1].atom2;
            }
            else {
                mid1=bonds[*bond1].atom1;
            }

            bond2 = bondsWithAtom[mid1];

            //  Loop through all the bonds connected to atom mid1
            while (*bond2 != -1) {
                if (bonds[*bond2].atom1 == mid1) {
                    mid2 = bonds[*bond2].atom2;
                }
                else {
                    mid2 = bonds[*bond2].atom1;
                }

            //  Make sure that we don't double back to where
            //  we started from.  This causes strange behavior.
            //  Trust me, I've been there . . .
            if (mid2 == i) {
                ++bond2;
                continue;
            }

            bond3=bondsWithAtom[mid2];

            //  Loop through all the bonds connected to mid2
            while (*bond3 != -1) {
                if (bonds[*bond3].atom1 == mid2) {
                    //  Make sure that we don't double back to where
                    //  we started from.  This causes strange behavior.
                    //  Trust me, I've been there . . .
                    //  I added this!!!  Why wasn't it there before?  -JCP
                    if (bonds[*bond3].atom2 != mid1)
                        if (i<bonds[*bond3].atom2) {
                            exclusionSet.add(Exclusion(i,bonds[*bond3].atom2,modified));
                        }
                }
                else {
                    //  Make sure that we don't double back to where
                    //  we started from.  This causes strange behavior.
                    //  Trust me, I've been there . . .
                    //  I added this!!!  Why wasn't it there before?  -JCP
                    if (bonds[*bond3].atom1 != mid1)
                    if (i<bonds[*bond3].atom1) {
                        exclusionSet.add(Exclusion(i,bonds[*bond3].atom1,modified));
                    }
                }

                ++bond3;
            }

            ++bond2;
            }
    
            ++bond1;
        }
    }
}
/*      END OF FUNCTION build14excl      */
    /************************************************************************/
    /*                                                                      */
    /*        FUNCTION stripHGroupExcl                                      */
    /*                                                                      */
    /*  This function removes all exclusions which are entirely             */
    /*  within a single hydrogen group.  This assumes that they all         */
    /*  exist, which should be true for exclusion policy 1-3 or higher.     */
    /*                                                                      */
    /************************************************************************/

void Initial::stripHGroupExcl(void) {

    HydrogenGroup::iterator h_i, h_e, h_j;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();

    for( ; h_i != h_e; ++h_i ) {
        for ( h_j = h_i + 1; h_j != h_e && ! h_j->isGP; ++h_j ) {
            if ( h_i->atomID < h_j->atomID )
                exclusionSet.del(Exclusion(h_i->atomID,h_j->atomID));
            else
                exclusionSet.del(Exclusion(h_j->atomID,h_i->atomID));
        }
    }

}
/*      END OF FUNCTION stripHGroupExcl      */
#if 0
    /************************************************************************/
    /*                  */
    /*      FUNCTION build_constraint_params    */
    /*                  */
    /*   INPUTS:                */
    /*  consref - Value of consref parameter from config file    */
    /*  conskfile - Value of conskfile from config file      */
    /*  conskcol - Value of conskcol from config file      */
    /*  initial_pdb - PDB object that contains initial positions  */
    /*  cwd - Current working directory          */
    /*                  */
    /*  This function builds all the parameters that are necessary  */
    /*   to do harmonic constraints.  This involves looking through    */
    /*   one or more PDB objects to determine which atoms are constrained,  */
    /*   and what the force constant and reference position is force each   */
    /*   atom that is constrained.  This information is then stored    */
    /*   in the arrays consIndexes and consParams.        */
    /*                  */
    /************************************************************************/

void Initial::build_constraint_params(StringList *consref, 
             StringList *conskfile, 
             StringList *conskcol, 
             PDB *initial_pdb,
             char *cwd) {
       
    PDB *refPDB, *kPDB;    //  Pointer to other PDB's if used
    register int i;      //  Loop counter
    int current_index=0;    //  Index into values used
    int kcol = 4;      //  Column to look for force constant in
    double kval = 0;      //  Force constant value retreived
    char filename[129];    //  PDB filename
       
    //  Get the PDB object that contains the reference positions.  If
    //  the user gave another file name, use it.  Otherwise, just use
    //  the PDB file that has the initial coordinates.  i.e., constrain
    //  the atoms around their initial position.  This is the most likely
    //  case anyway
    if (consref == NULL) {
        refPDB = initial_pdb;
    }
    else {
        if (consref->next != NULL) {
            R2_die("Multiple definitions of constraint reference file in configruation file");
        }

        if ( (cwd == NULL) || (consref->data[0] == '/') ) {
        strcpy(filename, consref->data);
        }
        else {
            strcpy(filename, cwd);
            strcat(filename, consref->data);
        }
    
        refPDB = new PDB(filename);
        if ( refPDB == NULL ) {
            R2_die("Memory allocation failed in Initial::build_constraint_params");
        }
        
        if (refPDB->num_atoms() != numAtoms) {
            R2_die("Number of atoms in constraint reference PDB doesn't match coordinate PDB");
        }
    }
       
    //  Get the PDB to read the force constants from.  Again, if the user
    //  gave us another file name, open that one.  Otherwise, just use
    //  the PDB with the initial coordinates
    if (conskfile == NULL) {
        kPDB = initial_pdb;
    }
    else {
        if (conskfile->next != NULL) {
            R2_die("Multiple definitions of constraint constant file in configuration file");
        }

        if ( (consref != NULL) && (strcasecmp(consref->data, conskfile->data) == 0) ) {
            //  Same PDB used for reference positions and force constants
            kPDB = refPDB; 
        }
        else {
            if ( (cwd == NULL) || (conskfile->data[0] == '/') ) {
                strcpy(filename, conskfile->data);
            }
            else {
                strcpy(filename, cwd);
                strcat(filename, conskfile->data);
            }

            kPDB = new PDB(filename);
            if ( kPDB == NULL ) {
                R2_die("Memory allocation failed in Initial::build_constraint_params");
            }
        
            if (kPDB->num_atoms() != numAtoms) {
                R2_die("Number of atoms in constraint constant PDB doesn't match coordinate PDB");
            }
        }
    }
       
    //  Get the column that the force constant is going to be in.  It
    //  can be in any of the 5 floating point fields in the PDB, according
    //  to what the user wants.  The allowable fields are X, Y, Z, O, or
    //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
    //  The default is the 4th field, which is the occupancy
    if (conskcol == NULL) {
        kcol = 4;
    }
    else {
        if (conskcol->next != NULL) {
            R2_die("Multiple definitions of harmonic constraint column in config file");
        }
    
        if (strcasecmp(conskcol->data, "X") == 0) {
            kcol=1;
        }
        else if (strcasecmp(conskcol->data, "Y") == 0) {
            kcol=2;
        }
        else if (strcasecmp(conskcol->data, "Z") == 0) {
            kcol=3;
        }
        else if (strcasecmp(conskcol->data, "O") == 0) {
            kcol=4;
        }
        else if (strcasecmp(conskcol->data, "B") == 0) {
            kcol=5;
        }
        else {
            R2_die("conskcol must have value of X, Y, Z, O, or B");
        }
    }
       
    //  Allocate an array that will store an index into the constraint
    //  parameters for each atom.  If the atom is not constrained, its
    //  value will be set to -1 in this array.
    consIndexes = new int[numAtoms];
       
    if (consIndexes == NULL) {
        R2_die("memory allocation failed in Initial::build_constraint_params()");
    }
       
    //  Loop through all the atoms and find out which ones are constrained
    for (i=0; i<numAtoms; i++) {
        //  Get the k value based on where we were told to find it
        switch (kcol) {
            case 1:
                kval = (kPDB->atom(i))->xcoor();
                break;
            case 2:
                kval = (kPDB->atom(i))->ycoor();
                break;
            case 3:
                kval = (kPDB->atom(i))->zcoor();
                break;
            case 4:
                kval = (kPDB->atom(i))->occupancy();
                break;
            case 5:
                kval = (kPDB->atom(i))->temperaturefactor();
                break;
        }
    
        if (kval > 0.0) {
            //  This atom is constrained
            consIndexes[i] = current_index;
            current_index++;
        }
        else {
            //  This atom is not constrained
            consIndexes[i] = -1;
        }
    }
       
    if (current_index == 0) {
        //  Constraints were turned on, but there weren't really any constrained
        cout <<   "NO CONSTRAINED ATOMS WERE FOUND, BUT CONSTRAINTS ARE ON . . . " << endl;
    }
    else {
        //  Allocate an array to hold the constraint parameters
        consParams = new ConstraintParams[current_index];
    
        if (consParams == NULL) {
            R2_die("memory allocation failed in Initial::build_constraint_params");
        }
    }
       
    numConstraints = current_index;
       
    //  Loop through all the atoms and assign the parameters for those
    //  that are constrained
    for (i=0; i<numAtoms; i++) {
        if (consIndexes[i] != -1) {
            //  This atom is constrained, so get the k value again
            switch (kcol) {
                case 1:
                    consParams[consIndexes[i]].k = (kPDB->atom(i))->xcoor();
                    break;
                case 2:
                    consParams[consIndexes[i]].k = (kPDB->atom(i))->ycoor();
                    break;
                case 3:
                    consParams[consIndexes[i]].k = (kPDB->atom(i))->zcoor();
                    break;
                case 4:
                    consParams[consIndexes[i]].k = (kPDB->atom(i))->occupancy();
                    break;
                case 5:
                    consParams[consIndexes[i]].k = (kPDB->atom(i))->temperaturefactor();
                    break;
            }
       
            //  Get the reference position
            consParams[consIndexes[i]].refPos.x = (refPDB->atom(i))->xcoor();
            consParams[consIndexes[i]].refPos.y = (refPDB->atom(i))->ycoor();
            consParams[consIndexes[i]].refPos.z = (refPDB->atom(i))->zcoor();
        }
    }
       
    //  If we had to create new PDB objects, delete them now
    if (consref != NULL) {
        delete refPDB;
    }
       
    if ((conskfile != NULL) &&
         !((consref != NULL) && 
        (strcasecmp(consref->data, conskfile->data) == 0))) {
        delete kPDB;
    }

}
    /*      END OF FUNCTION build_constraint_params    */

void Initial::build_langevin_params(double coupling, int doHydrogen) {

    //  Allocate the array to hold all the data
    langevinParams = new double[numAtoms];
    langForceVals = new double[numAtoms];

    if ( (langevinParams == NULL) || (langForceVals == NULL) ) {
        R2_die("memory allocation failed in Initial::build_langevin_params()");
    }

    //  Calculate the constant portion of the force values.  Note that
    //  because we need to convert from femtoseconds to picoseconds,
    //  the factor of 0.001 is needed.
    double forceConstant = 0.002*TIMEFACTOR*TIMEFACTOR*BOLTZMAN*
                (simParams->langevinTemp)/(simParams->dt);

    //  Loop through all the atoms and get the b value
    for (int i=0; i<numAtoms; i++) {
        double bval = coupling;

        if ( (! doHydrogen) && is_hydrogen(i) ) bval = 0;

        //  Assign the b value
        langevinParams[i] = bval;

        //  Calculate the random force value
        langForceVals[i] = sqrt(forceConstant*atoms[i].mass*bval);
    }

}

/************************************************************************/
/*                  */
/*      FUNCTION build_langevin_params      */
/*                  */
/*   INPUTS:                */
/*  langfile - Value of langevinfile from config file    */
/*  langcol - Value of langevincol from config file      */
/*  initial_pdb - PDB object that contains initial positions  */
/*      cwd - Current working directory          */
/*                  */
/*  This function builds the array of b values necessary for  */
/*   Langevin dynamics.  It takes the name of the PDB file and the      */
/*   column in the PDB file that contains the b values.  It then  */
/*   builds the array langevinParams for use during the program.  */
/*                  */
/************************************************************************/

void Initial::build_langevin_params(StringList *langfile, 
        StringList *langcol, 
        PDB *initial_pdb,
        char *cwd) {
    
    PDB *bPDB;      //  Pointer to PDB object to use
    int bcol = 4;      //  Column that data is in
    double bval = 0;      //  b value from PDB file
    int i;      //  Loop counter
    double forceConstant;  //  Constant factor in force calc
    char filename[129];    //  Filename
       
    //  Get the PDB object that contains the b values.  If
    //  the user gave another file name, use it.  Otherwise, just use
    //  the PDB file that has the initial coordinates.  
    if (langfile == NULL) {
        bPDB = initial_pdb;
    }
    else {
        if (langfile->next != NULL) {
            R2_die("Multiple definitions of langvein PDB file in configuration file");
        }

        if ( (cwd == NULL) || (langfile->data[0] == '/') ) {
            strcpy(filename, langfile->data);
        }
        else {
            strcpy(filename, cwd);
            strcat(filename, langfile->data);
        }
    
        bPDB = new PDB(filename);
        if ( bPDB == NULL ) {
            R2_die("Memory allocation failed in Initial::build_langevin_params");
        }
        
        if (bPDB->num_atoms() != numAtoms) {
            R2_die("Number of atoms in langevin parameter PDB doesn't match coordinate PDB");
        }
    }
       
    //  Get the column that the b vaules are in.  It
    //  can be in any of the 5 floating point fields in the PDB, according
    //  to what the user wants.  The allowable fields are X, Y, Z, O, or
    //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
    //  The default is the 4th field, which is the occupancy
    if (langcol == NULL) {
        bcol = 4;
    }
    else {
        if (langcol->next != NULL) {
            R2_die("Multiple definitions of langevin parameter column in config file");
        }
    
        if (strcasecmp(langcol->data, "X") == 0) {
            bcol=1;
        }
        else if (strcasecmp(langcol->data, "Y") == 0) {
            bcol=2;
        }
        else if (strcasecmp(langcol->data, "Z") == 0) {
            bcol=3;
        }
        else if (strcasecmp(langcol->data, "O") == 0) {
            bcol=4;
        }
        else if (strcasecmp(langcol->data, "B") == 0) {
            bcol=5;
        }
        else {
            R2_die("langevincol must have value of X, Y, Z, O, or B");
        }
    }
       
//  Allocate the array to hold all the data
    langevinParams = new double[numAtoms];
    langForceVals = new double[numAtoms];
      
    if ( (langevinParams == NULL) || (langForceVals == NULL) ) {
        R2_die("memory allocation failed in Initial::build_langevin_params()");
    }

    //  Calculate the constant portion of the force values.  Note that
    //  because we need to convert from femtoseconds to picoseconds,
    //  the factor of 0.001 is needed.  
    forceConstant = 0.002*TIMEFACTOR*TIMEFACTOR*BOLTZMAN*(simParams->langevinTemp)/(simParams->dt);
       
    //  Loop through all the atoms and get the b value
    for (i=0; i<numAtoms; i++) {
        //  Get the k value based on where we were told to find it
        switch (bcol) {
            case 1:
                bval = (bPDB->atom(i))->xcoor();
                break;
            case 2:
                bval = (bPDB->atom(i))->ycoor();
                break;
            case 3:
                bval = (bPDB->atom(i))->zcoor();
                break;
            case 4:
                bval = (bPDB->atom(i))->occupancy();
                break;
            case 5:
                bval = (bPDB->atom(i))->temperaturefactor();
                break;
        }
    
        //  Assign the b value
        langevinParams[i] = bval;

        //  Calculate the random force value
        langForceVals[i] = sqrt(forceConstant*atoms[i].mass*bval);
    }
       
    //  If we had to create a PDB object, delete it now
    if (langfile != NULL) {
        delete bPDB;
    }
}
/*      END OF FUNCTION build_langevin_params    */

/************************************************************************/
/*                  */
/*      FUNCTION build_fixed_atoms      */
/*                  */
/*   INPUTS:              */
/*  fixedfile - Value of langevinfile from config file    */
/*  fixedcol - Value of langevincol from config file    */
/*  initial_pdb - PDB object that contains initial positions  */
/*      cwd - Current working directory        */
/*                  */
/*  This function builds the list of fixed atoms.      */
/*   It takes the name of the PDB file and the      */
/*   column in the PDB file that contains the flags.  It then  */
/*   builds the array fixedAtomFlags for use during the program.  */
/*                  */
/************************************************************************/

void Initial::build_fixed_atoms(StringList *fixedfile, 
        StringList *fixedcol, 
        PDB *initial_pdb,
        char *cwd) {
    
    PDB *bPDB;      //  Pointer to PDB object to use
    int bcol = 4;      //  Column that data is in
    double bval = 0;      //  b value from PDB file
    int i;      //  Loop counter
    char filename[129];    //  Filename
       
    //  Get the PDB object that contains the b values.  If
    //  the user gave another file name, use it.  Otherwise, just use
    //  the PDB file that has the initial coordinates.  
    if (fixedfile == NULL) {
        bPDB = initial_pdb;
    }
    else {
        if (fixedfile->next != NULL) {
            R2_die("Multiple definitions of fixed atoms PDB file in configuration file");
        }

        if ( (cwd == NULL) || (fixedfile->data[0] == '/') ) {
            strcpy(filename, fixedfile->data);
        }
        else {
            strcpy(filename, cwd);
            strcat(filename, fixedfile->data);
        }
    
        bPDB = new PDB(filename);
        if ( bPDB == NULL ) {
            R2_die("Memory allocation failed in Initial::build_fixed_atoms");
        }
    
        if (bPDB->num_atoms() != numAtoms) {
            R2_die("Number of atoms in fixed atoms PDB doesn't match coordinate PDB");
        }
    }
       
    //  Get the column that the b vaules are in.  It
    //  can be in any of the 5 floating point fields in the PDB, according
    //  to what the user wants.  The allowable fields are X, Y, Z, O, or
    //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
    //  The default is the 4th field, which is the occupancy
    if (fixedcol == NULL) {
        bcol = 4;
    }
    else {
        if (fixedcol->next != NULL) {
            R2_die("Multiple definitions of fixed atoms column in config file");
        }
    
        if (strcasecmp(fixedcol->data, "X") == 0) {
            bcol=1;
        }
        else if (strcasecmp(fixedcol->data, "Y") == 0) {
            bcol=2;
        }
        else if (strcasecmp(fixedcol->data, "Z") == 0) {
            bcol=3;
        }
        else if (strcasecmp(fixedcol->data, "O") == 0) {
            bcol=4;
        }
        else if (strcasecmp(fixedcol->data, "B") == 0) {
            bcol=5;
        }
        else {
            R2_die("fixedatomscol must have value of X, Y, Z, O, or B");
        }
    }
       
    //  Allocate the array to hold all the data
    fixedAtomFlags = new int[numAtoms];
       
    if (fixedAtomFlags == NULL) {
        R2_die("memory allocation failed in Initial::build_fixed_atoms()");
    }
       
    numFixedAtoms = 0;

    //  Loop through all the atoms and get the b value
    for (i=0; i<numAtoms; i++) {
        //  Get the k value based on where we were told to find it
        switch (bcol) {
            case 1:
                bval = (bPDB->atom(i))->xcoor();
                break;
            case 2:
                bval = (bPDB->atom(i))->ycoor();
                break;
            case 3:
                bval = (bPDB->atom(i))->zcoor();
                break;
            case 4:
                bval = (bPDB->atom(i))->occupancy();
                break;
            case 5:
                bval = (bPDB->atom(i))->temperaturefactor();
                break;
        }
    
        //  Assign the b value
        if ( bval != 0 ) {
            fixedAtomFlags[i] = 1;
            numFixedAtoms++;
        }
        else {
            fixedAtomFlags[i] = 0;
        }
    }
       
    //  If we had to create a PDB object, delete it now
    if (fixedfile != NULL) {
        delete bPDB;
    }

    // now figure out how we interact with rigidBonds 
    // this is mainly for degree of freedom counting
    if ( numRigidBonds ) {
        HydrogenGroup::iterator h_i, h_e;
        h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
        int parentIsFixed = 0;
        for( ; h_i != h_e; ++h_i ) {
            if ( h_i->isGP ) {
                parentIsFixed = fixedAtomFlags[h_i->atomID];
                if ( (rigidBondLengths[h_i->atomID] != 0.)  // water
                    && fixedAtomFlags[h_i[1].atomID]
                    && fixedAtomFlags[h_i[2].atomID] ) {
                    ++numFixedRigidBonds;
                }
            } 
            else {
                if ( (rigidBondLengths[h_i->atomID] != 0.)
                    && fixedAtomFlags[h_i->atomID]
                    && parentIsFixed ) {
                    ++numFixedRigidBonds;
                }
            }
        }
    }

}
    /*      END OF FUNCTION build_fixed_atoms    */

#endif

int Initial::is_hydrogen(int anum) {
    return ((atoms[anum].status & HydrogenAtom) != 0);
}

int Initial::is_oxygen(int anum) {
    return ((atoms[anum].status & OxygenAtom) != 0);
}

int Initial::is_hydrogenGroupParent(int anum) {
    return (hydrogenGroup[atoms[anum].hydrogenList].isGP);
}

int Initial::is_water(int anum) {
    return (hydrogenGroup[atoms[anum].hydrogenList].sortVal == 2);
}

int Initial::get_groupSize(int anum) {
    return (hydrogenGroup[atoms[anum].hydrogenList].atomsInGroup);
}

// go through the molecular structure, analyze the status of each atom,
// and save the data in the Atom structures stored for each atom.  This
// could be built up incrementally while the Initial is being read in,
// but doing it all in one shot allows us to just send the basic info
// over the network and have each node calculate the rest of the data on
// it's own.
void Initial::build_atom_status(void) {
    register int i;
    int a1, a2;

    // initialize information for each atom (note that the status has
    // already been initialized during the read/receive phase)
    HydrogenGroupID *hg;
    hg = new HydrogenGroupID[numAtoms];
    for (i=0; i < numAtoms; i++) {
        hg[i].atomID = i;  // currently unsorted
        hg[i].atomsInGroup = 1;  // currently only 1 in group
        hg[i].isGP = 1;  // assume it is a group parent
        hg[i].GPID = i;  // assume it is a group parent
        hg[i].sortVal = 0;  // for group sorting
    }

    // find which atom each hydrogen is bound to
    // also determine number of atoms in each group
    for (i=0; i < numBonds; i++) {
        a1 = bonds[i].atom1;
        a2 = bonds[i].atom2;
        if (is_hydrogen(a1)) {
            // check for hydrogen gas...  For H2, explicitly define the group parent.
            // I have been informed that H3 is not a concern.  This is good since
            // this code will fail for H3.
            if (is_hydrogen(a2)) {
                hg[a1].isGP = 1;
                cout <<   "Found H-H bond - are you sure?" << endl;
            }
            else {
                hg[a2].atomsInGroup++;
                hg[a1].atomsInGroup = 0;
                hg[a1].GPID = a2;
                hg[a1].isGP = 0;
                // check for waters (put them in their own groups: OH or OHH)
                if (is_oxygen(a2))  hg[a2].sortVal++;
            }
        }
        if (is_hydrogen(a2)) {
            hg[a1].atomsInGroup++;
            hg[a2].atomsInGroup = 0;
            hg[a2].GPID = a1;
            hg[a2].isGP = 0;
            // check for waters (put them in their own groups: OH or OHH)
            if (is_oxygen(a1))  hg[a1].sortVal++;
        }
    }

    // sort the hydrogenGroup list and count number of groups
    numHydrogenGroups = 0;
    for(i=0; i<numAtoms; i++) {
        // make H follow their group parents.
        if (!hg[i].isGP)  hg[i].sortVal = hg[hg[i].GPID].sortVal;
        else ++numHydrogenGroups;
        // add to list to sort
        hydrogenGroup.push_back(hg[i]);
    }
    sort(hydrogenGroup.begin(), hydrogenGroup.end());

    // finally, add the indexing from atoms[] to hydrogenGroup[]
    waterIndex = numAtoms;
    for(i=0; i<numAtoms; i++) {
        atoms[hydrogenGroup[i].atomID].hydrogenList = i;
        // identify where waters start
        if ((hydrogenGroup[i].sortVal == 2) && (i < numAtoms))
            waterIndex = i;
    }

    #if 0
    // debugging code for showing sorted atoms
    for(i=0; i<numAtoms; i++)
        cout << i << " atomID=" << hydrogenGroup[i].atomID
        << " isGP=" << hydrogenGroup[i].isGP
        << " parent=" << hydrogenGroup[i].GPID
        << " #" << hydrogenGroup[i].atomsInGroup
        << " sortVal=" << hydrogenGroup[i].sortVal
        << "\n" << endl;
    #endif
    delete [] hg;

    #if 0
  // now deal with rigidBonds
    if ( simParams->rigidBonds != RIGID_NONE ) {
        delete [] rigidBondLengths;
        rigidBondLengths = new double[numAtoms];
        if ( ! rigidBondLengths ) 
            R2_die("Memory allocation failed in Initial::build_atom_status()\n");
            int mode = simParams->rigidBonds;

            // add H-mother lengths or 0 if not constrained
            for (i=0; i < numBonds; i++) {
                a1 = bonds[i].atom1;
                a2 = bonds[i].atom2;
                double dum, x0;
                params->get_bond_params(&dum,&x0,bonds[i].bond_type);
                if (is_hydrogen(a1)) {
                    if ( ! is_hydrogen(a2) && ( is_water(a2) || mode == RIGID_ALL ) ) {
                        rigidBondLengths[a1] = x0;
                        ++numRigidBonds;
                    }
                    else {
                        rigidBondLengths[a1] = 0.;
                    }
                }
                if (is_hydrogen(a2)) {
                    if ( is_water(a1) || mode == RIGID_ALL ) {
                        rigidBondLengths[a2] = x0;
                        ++numRigidBonds;
                    } 
                    else {
                        rigidBondLengths[a2] = 0.;
                    }
                }
            }

        // zero out H-H lengths - water handled below
        HydrogenGroup::iterator h_i, h_e;
        h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
        for( ; h_i != h_e; ++h_i ) {
            if ( h_i->isGP ) rigidBondLengths[h_i->atomID] = 0.;
        }

        // fill in H-H lengths for water by searching angles - yuck
        for (i=0; i < numAngles; i++) {
            a2 = angles[i].atom2;
            if ( ! is_water(a2) ) continue;
            a1 = angles[i].atom1;
            a3 = angles[i].atom3;
            if ( rigidBondLengths[a1] != rigidBondLengths[a3] ) {
                R2_die("Asymmetric water Initial found???  This can't be right.\n");
            }
            double dum, t0;
            params->get_angle_params(&dum,&t0,&dum,&dum,angles[i].angle_type);
            rigidBondLengths[a2] = 2. * rigidBondLengths[a1] * sin(0.5*t0);
            ++numRigidBonds;
        }

    }
    #endif
}
