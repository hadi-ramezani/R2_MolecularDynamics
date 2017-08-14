#include <stdio.h>
#include <string.h>
#include <iostream>
#include "PDB.h"


// read in a file and stick all the elements on the appropriate list
PDB::PDB( const char *pdbfilename) {
  FILE *infile;
  char buf[160];

  atomCount = 0;
  atomListHead = atomListTail = NULL;
  infile = fopen(pdbfilename, "r");
  if (! infile) {
     char s[500];
     sprintf(s, "Cannot open file '%s' for input in PDB::PDB.", pdbfilename);
     R2_die(s);
  }
  
    // loop through each line in the file and get a record
  while ( fgets(buf, 150, infile) ) {
   PDBData *newelement;
   char *s;
   for (s=buf; *s && *s!='\n'; s++)  // chop off the '\n'
    ;
   *s = 0;
   if ( s == (buf + 149) ) {
     char s[500];
     sprintf( s, "Input line too long in pdbfile %s.\n", pdbfilename);
     R2_die(s);
   }
   *(s+1) = 0;  // just to be on the safe side

     // I now have a string; make a PDBData element out of it
   newelement = new_PDBData(buf);
   if (!newelement) {
      R2_die("Could not allocate PDBData.\n");
   }
     // I only know how to deal with ATOM and HETATM types; and
     //  I want to throw away the unknown data types; so
   if (newelement -> type() != PDBData::ATOM && 
           newelement -> type() != PDBData::HETATM) {
       delete newelement;
   } else {
       add_atom_element( (PDBAtom *) newelement);
   }
  }  // input while loop
 fclose(infile);
 
 // now I have a linked list, and I know the size.  However,
 // that's got pretty slow access time for ramdom fetches, so
 // I'll turn it into an array
 {
  atomArray = new PDBAtomPtr[atomCount];
  if ( atomArray == NULL )
  {
    R2_die("memory allocation failed in PDB::PDB");
  }
  PDBAtomList *tmp = atomListHead;
  int i=0;                              // just need to copy the pointers
  for (i=0, tmp = atomListHead; tmp != NULL; tmp = tmp -> next, i++) {
    atomArray[i] = tmp -> data;
  }
     // now delete the linked list (w/o deleting the data)
  PDBAtomList *tmp2;
  for (tmp2 = tmp = atomListHead; tmp != NULL; tmp = tmp2) {
    tmp2 = tmp->next;
    delete tmp;
  }
  atomListHead = atomListTail = NULL;
 }  // everything converted
 
}

//  Destructor - delete all the data pointed to by the array
//   and then delete the array
PDB::~PDB( void )
{
	int i;
	for (i=atomCount-1; i>=0; i--)
	   delete atomArray[i];
	delete [] atomArray;
	atomArray = NULL;
	atomCount = 0;
}

// print the PDB file out to a given file name
void PDB::write(const char *outfilename, const char *commentline)
{
	int i;
	char s[200];
	FILE *outfile;
	if ((outfile = fopen(outfilename, "w")) == NULL) {
	   sprintf(s, "Cannot open file '%s' in PDB::write.", outfilename);
	   R2_die(s);
	}

	if (commentline != NULL)
	{
		sprintf(s, "REMARK  %s\n", commentline);
		if (fputs(s, outfile) == EOF)
		{
			R2_die("EOF in PDB::write writing the comment line - file system full?");
		}
	}

	for (i=0; i<atomCount; i++){ // I only contain ATOM/HETATM records
	  atomArray[i]->sprint(s, PDBData::COLUMNS);
	  if ( (fputs(s, outfile)    == EOF) || 
	       (fputc('\n', outfile) == EOF)    ) {
	    sprintf(s, "EOF in PDB::write line %d - file system full?", i);
	    R2_die(s);
	  }
	}
	if (fputs("END\n", outfile) == EOF) {
	   R2_die("EOF in PDB::write while printing 'END' -- file system full?");
	}
	if (fclose(outfile) == EOF) {
	   R2_die("EOF in PDB::write while closing -- file system full?");
	}
	  
}

// store the info on the linked list
void PDB::add_atom_element( PDBAtom *newAtom)
{
  PDBAtomList *tmp = new PDBAtomList;
  if ( tmp == NULL )
  {
    R2_die("memory allocation failed in PDB::add_atom_element");
  }
  tmp -> data = newAtom;
  tmp -> next = NULL;
  
  if (atomListHead == NULL) {        // make the list
    atomListHead = atomListTail = tmp;
  } else {
    atomListTail -> next = tmp;       // add to the tail
    atomListTail = tmp;
  }
  atomCount++;
}


// return the number of atoms found
int PDB::num_atoms( void)
{
  return atomCount;
}


// Reset all the atom positions.  This is used in preparation for
// output in cases like the restart files, etc.
void PDB::set_all_positions(Vector *pos)
{
	int i;
	PDBAtomPtr *atomptr;

	for (i=0, atomptr=atomArray; i<atomCount; atomptr++, i++)
	{
		(*atomptr)->xcoor(pos[i].x);
		(*atomptr)->ycoor(pos[i].y);
		(*atomptr)->zcoor(pos[i].z);
	}
}

//  Get all the atom positions into a list of Vectors
void PDB::get_all_positions(Vector *pos)
{
	int i;
	PDBAtomPtr *atomptr;

	for (i=0, atomptr=atomArray; i<atomCount; atomptr++, i++)
	{
		pos[i].x = (*atomptr)->xcoor();
		pos[i].y = (*atomptr)->ycoor();
		pos[i].z = (*atomptr)->zcoor();
	}
}

//  given an index, return that atom
PDBAtom *PDB::atom(int place)
{
  if (place <0 || place >= atomCount)
    return NULL;
  return atomArray[place];
}


// find the lowest and highest bounds to the atom
void PDB::find_extremes(Vector *low, Vector *high)
{
  PDBAtomPtr *atomptr = atomArray;
  PDBAtom *atom;
  double tmpcoor;
  int i;

  low->x=low->y=low->z=99999;   // larger than a legal PDB file allows -- .1mm!
  high->x=high->y=high->z=-99999;

  // search the array
  // the count down is for speed, I just use i as a counter, and it is
  //  quick to check against 0
  for (i=atomCount ; i>0; i--) {
    atom = *atomptr;
    if ( (low->x) > (tmpcoor=atom->xcoor()) )
	low->x = tmpcoor;
    if ( (high->x) < tmpcoor)
	high->x = tmpcoor;
    
    if ( (low->y) > (tmpcoor=atom->ycoor()) )
	low->y = tmpcoor;
    if ( (high->y) < tmpcoor)
	high->y = tmpcoor;
    
    if ( (low->z) > (tmpcoor=atom->zcoor()) )
	low->z = tmpcoor;
    if ( (high->z) < tmpcoor)
	high->z = tmpcoor;
    atomptr++;  // next!
  }
}
