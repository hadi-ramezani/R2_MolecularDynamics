//-*-c++-*-
/***************************************************************************/
/*       (C) Copyright 1995,1996,1997 The Board of Trustees of the         */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 * PDB Class
 *   Given a PDB file name, read in all the data.
 * As of now, you can only search the file for ATOM and 
 * HETATM information.  The return value is an IntList (of new'ed
 * memory so you have to delete it!) containing the list of all
 * fields that match that criterion, indexed by position in the file.
 * (Hence, 0 is the 1st ATOM or HETATM record, 10 is the eleventh,
 * and so on...).  Note that with these searches there is no
 * way to choose ATOM or HETATM; you have to make that distinguishment
 * yourself.
 ***************************************************************************/



#ifndef PDB_H
#define PDB_H

#include "PDBData.h"
#include "Vector.h"
#include "Lattice.h"
#include "R2.h"

using namespace std;

typedef PDBAtom *PDBAtomPtr ;
typedef struct PAL {
  PDBAtom *data;
  struct PAL *next;
} PDBAtomList;
  
class PDB {
  private:
    PDBAtomList *atomListHead, *atomListTail;
    PDBAtom **atomArray;
      // this doesn't create a copy 
    void add_atom_element(PDBAtom *newAtom); 
    int atomCount;
    
  public:
    PDB( const char *pdbfilename);   // read in PDB from a file
    ~PDB( void);               // clear everything
    void write(const char *outfilename, const char *commentline=NULL); // write the coordinates to a file
       // the following deals only with ATOMs and HETATMs
    int num_atoms( void);

    PDBAtom *atom(int place); // get the nth atom in the PDB file
         // return linked list containing all atoms
    PDBAtomList *atoms(void ) { return atomListHead; }  
         
	// Find the extreme edges of the molecule
    void find_extremes(Vector *low, Vector *high);
    void find_99percent_extremes(Vector *low, Vector *high);

    void set_all_positions(Vector *);	//  Reset all the positions in PDB

    void get_all_positions(Vector *);	//  Get all positions in PDB
};

#endif // PDB_H
