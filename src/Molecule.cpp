/* 
 * File:   Molecule.cpp
 * Author: amin
 * 
 * Created on September 9, 2015, 2:24 PM
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include "Molecule.h"
#include "Configure.h"


using namespace std;

Molecule::Molecule(const Configure *conf) {
    atoms = new Atom[conf->numAtoms];
    bonds = new Bond[conf->numBonds];
    angles = new Angle[conf->numAngles]; 

    read_psf(conf->psfname);
    read_pdb(conf->pdbname);
}

void Molecule::read_psf(const char *filename){
    
    string str;
    char   seg[256];
    char  *pch;
    size_t found;
    ifstream infile;
    infile.open(filename);

    while(!infile.eof()) {
	getline(infile, str);        

        found=str.find("NATOM");
        // PSF file has 9 columns 
        // atom ID, segment name, residue ID, residue name, atom name, atom type, charge, mass, and an unused 0
        // I only need residue ID, residue name, atom name, atom type, charge, and mass
        if (found!=std::string::npos) {
            Atom *patoms = atoms;
            strncpy(seg, str.c_str(), sizeof(seg));
            seg[sizeof(seg) - 1] = 0;
            pch = strtok (seg," ,:=");
            int numAtoms = stoi(pch); // Read the number of atoms again
            for (int ii=0; ii<numAtoms; ii++){
                getline(infile, str);
                strncpy(seg, str.c_str(), sizeof(seg));
                seg[sizeof(seg) - 1] = 0;
                pch = strtok (seg," ,:="); pch = strtok (NULL, " ,:="); // jump two columns
                pch = strtok (NULL, " ,:="); 
                patoms->res_ID = stoi(pch); // residue ID
                pch = strtok (NULL, " ,:="); strcpy(patoms->resname,pch);
                pch = strtok (NULL, " ,:="); strcpy(patoms->atomname,pch);
                pch = strtok (NULL, " ,:="); strcpy(patoms->atomtype,pch);
                patoms->charge = stof(strtok (NULL, " ,:="));
                patoms->mass = stof(strtok (NULL, " ,:="));
                ++patoms; // next atom
            }     
        }
        
// Read bonds connections                 
        found=str.find("NBOND");
        if (found!=std::string::npos) {
            strncpy(seg, str.c_str(), sizeof(seg));
            seg[sizeof(seg) - 1] = 0;
            pch = strtok (seg," ,:=");
            int numBonds = stoi(pch); // Read the number of bonds again
            Bond *pbonds = bonds; // Pointer to first of bonds array
            int ii=0; 
            getline(infile, str);
            while (str != "" && ii<numBonds){
                strncpy(seg, str.c_str(), sizeof(seg));
                seg[sizeof(seg) - 1] = 0;
                pch = strtok (seg," ,:=");
                while (pch != NULL) { // Each line has Max 4 bonds
                    pbonds->atom1 = stoi(pch) - 1; // Read atom1
                    pch = strtok (NULL," ,:="); 
                    pbonds->atom2 = stoi(pch) - 1; // Read atom2
                    pch = strtok (NULL," ,:=");
                    ++pbonds; ii++;// Increase the pointer to move to next bond
                } 
                getline(infile, str); // Read next line
            }
        }
        
// Read angles connections        
        found=str.find("NTHETA");
        if (found!=std::string::npos) {
            strncpy(seg, str.c_str(), sizeof(seg));
            seg[sizeof(seg) - 1] = 0;
            pch = strtok (seg," ,:=");
            int numAngles = stoi(pch); // Read the number of angles again
            Angle *pangles = angles; // Pointer to the first of angles array
            int ii=0;
            getline(infile, str);
            while (str != "" && ii<numAngles){
                strncpy(seg, str.c_str(), sizeof(seg));
                seg[sizeof(seg) - 1] = 0;
                pch = strtok (seg," ,:=");
                while (pch != NULL) { // Each line has Max 4 bonds
                    pangles->atom1 = stoi(pch) - 1; // Read atom1
                    pch = strtok (NULL," ,:="); 
                    pangles->atom2 = stoi(pch) - 1; // Read atom2
                    pch = strtok (NULL," ,:=");
                    pangles->atom3 = stoi(pch) - 1; // Read atom3
                    pch = strtok (NULL," ,:=");
                    ++pangles; ii++;// Increase the pointer to move to next bond
                } 
                getline(infile, str); // Read next line
            }
           
        }  
    }
    infile.close();
}

void Molecule::read_pdb(const char *filename){
    
    string str;
    char   seg[8];
    char  *pch;
    size_t found;
    Atom *patoms = atoms;
    ifstream infile;
    infile.open(filename);
    
    while(!infile.eof()) {
	getline(infile, str);
        found=str.find("ATOM");
        if (found!=std::string::npos) {
            // Read atoms coordinates
            str.copy(seg,8,30);
            patoms->coorx = stof(seg) ; // X            
            str.copy(seg,8,38);
            patoms->coory = stof(seg) ; // Y
            str.copy(seg,8,46);
            patoms->coorz = stof(seg) ; // Z
            patoms++;
        }
    }
}

Molecule::Molecule(const Molecule& orig) {
}

Molecule::~Molecule() {
}



