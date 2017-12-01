#include <fstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sstream>
#include "Trajectory.h"


using namespace std;

Trajectory::Trajectory(const char *filename, int natoms, const Configure *conf) {
    if (!(strncasecmp(conf->analysis, "on", 2) == 0)){
        write_header(filename, natoms);        
    } else {
        open_dcd_get_info(filename, natoms, conf);
        //read_header(filename, natoms);
    }
}

void Trajectory::open_dcd_get_info(const char *filename, int natoms, const Configure *conf) {
    
    // Check if the file is open
    dcdf.open(filename, ios::in|ios::binary);

    if (!dcdf) {
        cout << "Error: could not read the input dcd file " << endl;
        exit(1);
    }


    read_header(natoms, nsets, istart, nsavc, delta, nfixed, freeind,
        fixedcoords, reverse, &charmm, conf);

    int ndims = 3; float newnsets;
    firstframesize = (natoms+4) * ndims * sizeof(float);
    framesize = (natoms-nfixed) * ndims * sizeof(float) + 8 * sizeof(int) + 6 * sizeof(double);

    header_size = dcdf.tellg(); /* save current offset (end of header) */
    dcdf.seekg(0, dcdf.end);
    filesize = dcdf.tellg();
    dcdf.seekg(header_size, dcdf.beg);

    trajsize = filesize - header_size - firstframesize;

    if (trajsize < 0) {
        cout << "the dcd file appears to contain no timesteps." << endl;
        dcdf.close();
        exit(1);
    }

    newnsets = trajsize / framesize + 1;
    cout << "Reading " << newnsets << " frames" << endl;

    if (nsets > 0 && newnsets != nsets) {
        cout << "Warning: DCD header claims " << nsets << 
                " frames, file size indicates there are actually " 
                << newnsets << "frames" << endl;
    }

    if (conf->dcdLast > newnsets || conf->dcdLast < 0){
        nsets = newnsets;
    } else {
        nsets = conf->dcdLast;
    }
    setsread = 0;

    first = 1;

    // Allocate memory for coordinates
    X = new float[natoms];
    Y = new float[natoms];
    Z = new float[natoms];

    if (!X || !Y || !Z) {
        cout << "Unable to allocate space to store coordinates..." << endl;
     if (X) delete[] X;
     if (Y) delete[] Y;
     if (Z) delete[] Z;
     dcdf.close();
   }
}

void Trajectory::read_header(int natoms, int nsets, int istart, int nsavc,
     double delta, int namnf, int* freeind, float* fixedcoords, int reverse,
     int* charmm, const Configure *conf) {

    unsigned int input_integer[2];  /* buffer space */
    int i, ret_val, rec_scale;
    char hdrbuf[84];    /* char buffer used to store header */
    int NTITLE;
    int dcdcordmagic;
    char* corp = (char *) &dcdcordmagic;

    /* coordinate dcd file magic string 'CORD' */
    corp[0] = 'C';
    corp[1] = 'O';
    corp[2] = 'R';
    corp[3] = 'D';

    // First thing in the file should be an 84.
    dcdf.read((char*)input_integer, 2*sizeof(unsigned int));
    if (input_integer[0] == 84 && input_integer[1] == dcdcordmagic) {
        rec_scale=RECSCALE32BIT;
        cout << "detected standard 32-bit DCD file of native endianness\n" << endl;
    }

    /* Buffer the entire header for random access */
    dcdf.read((char *)hdrbuf, 80);

    /* CHARMm-genereate DCD files set the last integer in the     */
    /* header, which is unused by X-PLOR, to its version number.  */
    /* Checking if this is nonzero tells us this is a CHARMm file */
    /* and to look for other CHARMm flags.                        */
    if (hdrbuf[76] != 0) {
        cout << "DCD is charmm";
        (*charmm) = DCD_IS_CHARMM;
    }

    if (*charmm & DCD_IS_CHARMM) {
        /* CHARMM and NAMD versions 2.1b1 and later */
        cout << "CHARMM format DCD file (also NAMD 2.1 and later)" << endl;;
    }

    /* Store the number of sets of coordinates (nsets) */
    nsets = hdrbuf[0];
    /* Store istart, the starting timestep */
    istart = hdrbuf[4];
    /* Store nsvac, the number of timesteps between dcd saves */
    nsavc = hdrbuf[8];
    /* Store namnf, the number of fixed atoms */
    namnf = hdrbuf[32];

    /* Read in the timestep, DELTA */
    /* Note: DELTA is stored as a double with X-PLOR but as a float with CHARMm */
    if ((*charmm) & DCD_IS_CHARMM) {
        float ftmp;
        ftmp = hdrbuf[36];
        delta = (double)ftmp; 
    }

    /* Get the end size of the first block */
    dcdf.read((char *) input_integer, rec_scale*sizeof(int));

    /* Read in the size of the next block */
    input_integer[1] = 0;
    dcdf.read((char *) input_integer, rec_scale*sizeof(int));

     if ((((input_integer[0] + input_integer[1])-4) % 80) == 0) {
        /* Read NTITLE, the number of 80 character title strings there are */
        dcdf.read((char *) &NTITLE, sizeof(int));

        if (NTITLE < 0) {
            cout << "WARNING: Bogus NTITLE value..." << endl;
            exit(1);
        }

        if (NTITLE > 1000) {
            cout << "WARNING: Bogus NTITLE value..." << endl;

            if (NTITLE == 1095062083) {
                cout << "WARNING: Broken Vega ZZ 2.4.0 DCD file detected" << endl;
                cout << "Assuming 2 title lines, good luck..." << endl;
                NTITLE = 2;
            } else {
                cout << "Assuming zero title lines, good luck..." << endl;
                NTITLE = 0;
            }
        }
        for (i=0; i<NTITLE; i++) {
            dcdf.seekg(80, ios::cur);
        }
        /* Get the ending size for this block */
        dcdf.read((char *) input_integer, rec_scale*sizeof(int));
    } else {
            cout << "Bad DCD format..." << endl;
            exit(1);
    }
    /* Read in an integer '4' */
    input_integer[1] = 0;
    dcdf.read((char *) input_integer, rec_scale*sizeof(int));

     /* Read in the number of atoms */
    dcdf.read((char *) input_integer, rec_scale*sizeof(int));
    natoms = *input_integer;

    /* Read in an integer '4' */
    input_integer[1] = 0;
    dcdf.read((char *) input_integer, rec_scale*sizeof(int));

    freeind = NULL;
    fixedcoords = NULL;
    if (namnf != 0) {
        freeind = (int *) calloc(((natoms)-(namnf)), sizeof(int));
        if (freeind == NULL) {
            cout << "Could not allocate memory..." << endl;
            exit(1);
        }
        fixedcoords = (float *) calloc(natoms*4 - namnf, sizeof(float));
        if (fixedcoords == NULL) {
            cout << "Could not allocate memory..." << endl;
            exit(1);
        }
        /* Read in index array size */
        dcdf.read((char *) input_integer, rec_scale*sizeof(int));
        dcdf.read((char *) freeind, natoms-namnf*sizeof(int));
        input_integer[1] = 0;
        dcdf.read((char *) input_integer, rec_scale*sizeof(int));
    }
}

void Trajectory::skip_frames(const int dcdStep){
    dcdf.seekg((dcdStep - 1) * framesize, dcdf.cur);
}

bool Trajectory::read_frame(int natoms, Vector* pos, double* aBox, const Configure *conf)  {

    int out;
    boxdcd = new double[6];
    aBox[0] = aBox[1] = aBox[2] = 1.0;
    aBox[3] = aBox[4] = aBox[5] = 90.0;
    if (setsread == 0){
        skip_frames(conf->dcdFirst);
        setsread+=conf->dcdFirst;
    }
    if (setsread  >= nsets) return false;

    dcdf.read((char*)&out, sizeof (unsigned int));

    // Read box information
    dcdf.read((char*)boxdcd, out);
    // Save box information
    aBox[0] = boxdcd[0]; aBox[1] = boxdcd[2]; aBox[2] = boxdcd[5];
    dcdf.read((char*)&out, sizeof (unsigned int));

    dcdf.read((char*)&out, sizeof (int));
    dcdf.read((char*)X, natoms*sizeof (float));
    dcdf.read((char*)&out, sizeof (int));

    dcdf.read((char*)&out, sizeof (int));
    dcdf.read((char*)Y, natoms*sizeof (float));
    dcdf.read((char*)&out, sizeof (int));

    dcdf.read((char*)&out, sizeof (int));
    dcdf.read((char*)Z, natoms*sizeof (float));
    dcdf.read((char*)&out, sizeof (int));

    for (int ii=0; ii<natoms; ii++){
        pos[ii].x = X[ii]; //- floor(coor[ii].x/box[0])*box[0];
        pos[ii].y = Y[ii]; //- floor(coor[ii].y/box[1])*box[1];
        pos[ii].z = Z[ii]; //- floor(coor[ii].z/box[2])*box[2];
    }

    skip_frames(conf->dcdStep);
    setsread += conf->dcdStep;

    if (!dcdf) {
        return false;
    } else {
        return true;
    }
}

void Trajectory::write_header(const char *filename,int natoms) {

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

void Trajectory::write_frame(int natoms, const Vector *coor, const double *box)  {


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

/*void Trajectory::read_header(const char *filename, int natoms) {

    int out;
    char buff[80];
    char title[200];


    // read the dcd header
    dcdf.read((char*) &out, sizeof (int));

    strcpy(buff, "CORD");
    dcdf.read(buff, 4);

    dcdf.read((char*)&out, sizeof (int));
    dcdf.read((char*)&out, sizeof (int));
    dcdf.read((char*)&out, sizeof (int));

    for (int ii=0; ii<6; ii++){
        dcdf.read((char*)&out, sizeof (int));
    }

    float outf = 2.0;
    dcdf.read((char*)&outf, sizeof (double));

    for (int ii=0; ii<8; ii++){
        dcdf.read((char*)&out, sizeof (int));
    }

    dcdf.read((char*)&out, sizeof (int));
    dcdf.read((char*)&out, sizeof (int));
    dcdf.read((char*)&out, sizeof (int));

    dcdf.read((char*)&out, sizeof (int));
    sprintf(title,"REMARKS FILENAME= out.dcd CREATED BY NAMD");
    title[79]='\0';
    dcdf.read(title, 80);
    sprintf(title,"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    dcdf.read(buff, 80);

    dcdf.read((char*)&out, sizeof (int));

    dcdf.read((char*)&out, sizeof (int));
    dcdf.read((char*)&natoms, sizeof (int));
    dcdf.read((char*)&out, sizeof (int));

    X = new float[natoms];
    Y = new float[natoms];
    Z = new float[natoms];
    boxdcd = new double[6];
}*/

Trajectory::Trajectory(const Trajectory& orig) {
}

Trajectory::~Trajectory() {
    dcdf.close();
}

