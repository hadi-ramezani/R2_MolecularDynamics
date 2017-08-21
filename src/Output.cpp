#include "Output.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>

Output::Output(const char *filename, const Configure* conf) {
    outf.open(filename);
    if ((strncasecmp(conf->analysis, "on", 2) == 0)){
        
        outf << setw(width) << right << "#Frame"
            << setw(width) << right << "Ebond"
            << setw(width) << right << "Eangle"
            << setw(width) << right << "Edihedral"
            << setw(width) << right << "Eimproper"
            << setw(width) << right << "Evdw"
            << setw(width) << right << "Eelec"
            << setw(width) << right << "Emisc"
            << setw(width) << right << "Etotal" << endl;

    } else {
        outf << setw(width) << right << "#Step" 
            << setw(width) << right << "Time" 
            << setw(width) << right << "Ebond"
            << setw(width) << right << "Eangle"
            << setw(width) << right << "Edihedral"
            << setw(width) << right << "Eimproper"
            << setw(width) << right << "Evdw"
            << setw(width) << right <<  "Eelec"
            << setw(width) << right << "Ekin"
            << setw(width) << right << "Etotal"
            << setw(width) << right << "Temp" << endl;
    }
}

// Normal simulation
void Output::print(const int Step, const double Time, const double Ebond, const double Eangle, const double Edihedral, const double Eimproper, 
    const double Evdw, const double Eelec, const double Ekin, const double Etot, const double temp){
    outf << setw(width) << right << Step 
        << setw(width) << right << Time 
        << setw(width) << right << Ebond
        << setw(width) << right << Eangle
        << setw(width) << right << Edihedral
        << setw(width) << right << Eimproper
        << setw(width) << right << Evdw 
        << setw(width) << right << Eelec
        << setw(width) << right << Ekin
        << setw(width) << right << Etot
        << setw(width) << right << temp << endl;
//    printf("%10d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",Step, Time, Ebond, Eangle, Evdw, Eelec, Ekin, Etotal);
}

// Analysis output
void Output::print(const unsigned int frameNum, const double Ebond, const double Eangle, const double Edihedral, const double Eimproper,
    const double Evdw, const double Eelec, const double Emisc, const double Etot) {
    
    outf << setw(width) << right << frameNum
        << setw(width) << right << Ebond 
        << setw(width) << right << Eangle 
        << setw(width) << right << Edihedral
        << setw(width) << right << Eimproper
        << setw(width) << right << Evdw 
        << setw(width) << right << Eelec
        << setw(width) << right << Emisc
        << setw(width) << right << Etot << endl;    
}

void Output::wrap(const Initial *init, const double *box, Vector *const pos){

    /*int iatom = 0;
    for (int mm=0; mm<init->nresidue; mm++){
        bool flag = true;  // the residue is out of the box
        int hatom = iatom;
        while (mm==init->residue[iatom]){
            if ( (pos[iatom].x >= 0 && pos[iatom].x <= box[0]) ){
                 flag = false; // at lease one atom of residue is in the box
            }
            iatom++;
        }
        if (flag) {
            for (int ii=hatom; ii<iatom; ii++){
                pos[ii].x = pos[ii].x - floor(pos[ii].x/box[0])*box[0];
            }
        }
    }
    iatom = 0;
    for (int mm=0; mm<init->nresidue; mm++){
        bool flag = true;  // the residue is out of the box
        int hatom = iatom;
        while (mm==init->residue[iatom]){
            if ( (pos[iatom].y >= 0 && pos[iatom].y <= box[1]) ) {
                 flag = false; // at lease one atom of residue is in the box
            }
            iatom++;
        }
        if (flag) {
            for (int ii=hatom; ii<iatom; ii++){
                pos[ii].y = pos[ii].y - floor(pos[ii].y/box[1])*box[1];
            }
        }
    }
    iatom = 0;
    for (int mm=0; mm<init->nresidue; mm++){
        bool flag = true;  // the residue is out of the box
        int hatom = iatom;
        while (mm==init->residue[iatom]){
            if ( (pos[iatom].z >= 0 && pos[iatom].z <= box[2]) ){
                 flag = false; // at lease one atom of residue is in the box
            }
            iatom++;
        }
        if (flag) {
            for (int ii=hatom; ii<iatom; ii++){
                pos[ii].z = pos[ii].z - floor(pos[ii].z/box[2])*box[2];
            }
        }
    }*/
}

Output::Output(const Output& orig) {
}

Output::~Output() {
    outf.close();
}

