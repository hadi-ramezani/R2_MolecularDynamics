/*
 * File:   Output.cpp
 * Author: amin
 *
 * Created on September 28, 2015, 5:23 PM
 */

#include "Output.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>

Output::Output(const char *filename, const Configure* conf) {
    outf.open(filename);
    if ((strncasecmp(conf->analysis, "on", 2) == 0)){
        outf << "#Frame Ebond Eangle Evdw Eelec Etotal" << endl;
    } else {
        outf << "#Step Time Ebond Eangle Evdw Eelec Ekin Etotal" << endl;
    }
}

void Output::Print(const int Step, const double Time, const double Ebond, const double Eangle, const double Evdw, const double Eelec, const double Ekin, const double Etotal, const double temp){
    outf << Step << " " << Time << " " << Ebond << " " << Eangle << " " << Evdw << " " << Eelec << " " << Ekin << " " << Etotal << " " << temp << endl;
//    printf("%10d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",Step, Time, Ebond, Eangle, Evdw, Eelec, Ekin, Etotal);
}

void Output::Print(const unsigned int frameNum, const double Ebond, const double Eangle, const double Evdw, const double Eelec, const double Etotal) {
    outf << frameNum << " " << Ebond << " " << Eangle << " " << Evdw << " " << Eelec << " " << Etotal << endl;    
}

void Output::wrap(const Initial *init, const double *box, Vector *const pos){

    int iatom = 0;
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
    }
}

Output::Output(const Output& orig) {
}

Output::~Output() {
    outf.close();
}

