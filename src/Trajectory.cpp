/*
 * File:   Trajectory.cpp
 * Author: amin
 *
 * Created on September 24, 2015, 12:20 PM
 */

#include <fstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sstream>
#include "Trajectory.h"

using namespace std;

Trajectory::Trajectory(const char *filename, int natoms) {

    int out;
    char buff[80];
    char title[200];

    dcdf.open(filename, ios::out|ios::binary);

    // write the dcd header
    out = 84;
    dcdf.write ((char*)&out, sizeof (int));
    strcpy(buff, "CORD");
    dcdf.write (buff, 4);
    out=0;
    dcdf.write ((char*)&out, sizeof (int));
    out=1;
    dcdf.write ((char*)&out, sizeof (int));
    out=1000;
    dcdf.write ((char*)&out, sizeof (int));
    out=0;

    for (int ii=0; ii<6; ii++){
        dcdf.write ((char*)&out, sizeof (int));
    }

    float outf = 2.0;
    dcdf.write ((char*)&outf, sizeof (double));

    for (int ii=0; ii<8; ii++){
        dcdf.write((char*)&out, sizeof (int));
    }

    out=24;
    dcdf.write ((char*)&out, sizeof (int));
    out=84;
    dcdf.write ((char*)&out, sizeof (int));
    out=164;
    dcdf.write ((char*)&out, sizeof (int));

    out=2;
    dcdf.write ((char*)&out, sizeof (int));
    sprintf(title,"REMARKS FILENAME= out.dcd CREATED BY NAMD");
    title[79]='\0';
    dcdf.write (title, 80);
    sprintf(title,"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    dcdf.write (buff, 80);

    out=164;
    dcdf.write ((char*)&out, sizeof (int));

    out=4;
    dcdf.write ((char*)&out, sizeof (int));
    dcdf.write ((char*)&natoms, sizeof (int));
    dcdf.write ((char*)&out, sizeof (int));

    X = new float[natoms];
    Y = new float[natoms];
    Z = new float[natoms];
    boxdcd = new double[6];
}

void Trajectory::WriteFrame(int natoms, const Vector *coor, const double *box)  {

    boxdcd[0]=box[0];boxdcd[1]=90.0; boxdcd[2]=box[1];
    boxdcd[3]=90.0;  boxdcd[4]=90.0;   boxdcd[5]=box[2];

    int out = 48;
    dcdf.write ((char*)&out, sizeof (unsigned int));
    dcdf.write ((char*)boxdcd, out);
    dcdf.write ((char*)&out, sizeof (unsigned int));

    for (int ii=0; ii<natoms; ii++){
        X[ii] = coor[ii].x ;//- floor(coor[ii].x/box[0])*box[0];
        Y[ii] = coor[ii].y ;//- floor(coor[ii].y/box[1])*box[1];
        Z[ii] = coor[ii].z ;//- floor(coor[ii].z/box[2])*box[2];
    }

    out = natoms*4;
    dcdf.write ((char*)&out, sizeof (int));
    dcdf.write ((char*)X, natoms*sizeof (float));
    dcdf.write ((char*)&out, sizeof (int));

    dcdf.write ((char*)&out, sizeof (int));
    dcdf.write ((char*)Y, natoms*sizeof (float));
    dcdf.write ((char*)&out, sizeof (int));

    dcdf.write ((char*)&out, sizeof (int));
    dcdf.write ((char*)Z, natoms*sizeof (float));
    dcdf.write ((char*)&out, sizeof (int));

}
Trajectory::Trajectory(const Trajectory& orig) {
}

Trajectory::~Trajectory() {
    dcdf.close();
}

