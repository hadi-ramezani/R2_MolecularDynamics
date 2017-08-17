#include "Nonbonded.h"
#include "Parameters.h"
#include "Initial.h"
#include "LJTable.h"

Nonbonded::Nonbonded(const Initial *init, 
                                   const Parameters *params,
                                   const Configure *conf) {

    natoms = init->numAtoms;
    cut2 = conf->cutoff;
    cut2 *= cut2;
    switch2 = conf->switchdist;
    switch2 *= switch2;
    pair2 = conf->pairlistdist;
    pair2 *= pair2;
    c1 = 1.0/(cut2-switch2);
    c1 = c1*c1*c1;
    c3 = 4*c1;
    ljTable = new LJTable(params);
    
    poshift = new Vector[natoms];
    memset((void *)poshift, 0, natoms*sizeof(Vector));
}


Nonbonded::~Nonbonded() {
    delete ljTable;
}

void Nonbonded::Build_cells(const Initial *init, const Vector *pos, Vector *f, const Configure *conf) {


    double pairdist = conf->pairlistdist;

    xb = int(conf->box[0] / pairdist) + 1;
    yb = int(conf->box[1] / pairdist) + 1;
    zb = int(conf->box[2] / pairdist) + 1;

    xytotb = yb * xb;
    totb = xytotb * zb;

    nbrlist = new int *[totb];

    //
    // Create neighbor list for each cell
    /* How to construct a cell list? Let assume we have a 3*3*3 grid.
       In our calculations, we only need to count 14 neighbors because we loop
       over all boxes and we don't want to double count things. However, we need
       to be careful to do this corretly. The correct neighbor cells include a 
       L-shaped list (4 neighbors) plus a row above the cell (9 cells). For example,
       the neighbors of cell 13 here are cells 16, 17, 14, 11 and the whole top 
       row (18, 19, 20, 21, 22, 23, 24, 25, 26). In the periodic calculations, we 
       just need to account for the periodic neighbors.
    Bottom row: 6  7  8   Middle row:  15  16  17  Top row:  24  25  26
                3  4  5                12  13  14            21  22  23
                0  1  2                 9  10  11            18  19  20
    */

    int xnb = 0, ynb = 0, znb = 0;
    int aindex = 0, aindexNB = 0;
    for (int zi=0; zi<zb; zi++) {
        for (int yi=0; yi<yb; yi++) {
            for (int xi=0; xi<xb; xi++) {
                int nbrs[14];           // Max possible number of neighbors in 3D
                int n=0;                // Number of neighbors found so far

                aindex = zi * xytotb + yi * xb + xi;
                nbrs[n++] = aindex;
                int counter = 1;
                // Find neighbors and take care of PBC
                for (int zz = 0; zz < 2; zz++){
                    for (int yy = -1; yy < 2; yy++){
                        for (int xx = -1; xx < 2; xx++){

                            znb = zi + zz;
                            //if (znb < 0 ) znb = zb -1;
                            if (znb == zb) znb = 0;

                            ynb = yi + yy;
                            if (ynb < 0 ) ynb = yb -1;
                            else if (ynb == yb) ynb = 0;

                            xnb = xi + xx;
                            if (xnb < 0 ) xnb = xb -1;
                            else if (xnb == xb) xnb = 0;

                            aindexNB = znb * xytotb + ynb * xb + xnb;
                            if (!((zz == 0 && xx == -1 && yy == 1) ||
                                (zz == 0 && xx == -1 && yy == 0) ||
                                (zz == 0 && xx == -1 && yy == -1) ||
                                (zz == 0 && xx == 0 && yy == -1) ||
                                (zz == 0 && xx == 0 && yy == 0))) {
                                    nbrs[n++] = aindexNB;
                            }
                        }
                    }
                }
                nbrlist[aindex] = new int[n+1];
                memcpy((void *)nbrlist[aindex], (void *)nbrs, n*sizeof(int));
                nbrlist[aindex][n] = -1;  // Sentinel for end of neighbors
            }
        }
    }
}

void Nonbonded::Build_Neighborlist(const Initial *init, const Vector *pos, Vector *f, const Configure *conf){
    
    int aindex;
    
    // shift atom coordinates inside the box (i.e. wrap them) for nonbonded calculations
    // Currently we assume the origin is at (0, 0, 0)
    // I'll change this in the future if needed
    for (int i=0; i<natoms; i++){ // For a box with the origin at the lower left vertex
        poshift[i].x = pos[i].x - floor(pos[i].x/conf->box[0]) * conf->box[0];
        poshift[i].y = pos[i].y - floor(pos[i].y/conf->box[1]) * conf->box[1];
        poshift[i].z = pos[i].z - floor(pos[i].z/conf->box[2]) * conf->box[2];
    }
    
    // add small number to make sure all particles are located in the cells
    // avoid the problem for particles located on the border
    double pairdistx = conf->box[0]/xb + 0.001; 
    double pairdisty = conf->box[1]/yb + 0.001; 
    double pairdistz = conf->box[2]/zb + 0.001;

    boxatom = new atominfo*[totb];
    numinbox = new int[totb];
    maxinbox = new int[totb];
    memset((void *)numinbox, 0,totb*sizeof(int));
    memset((void *)maxinbox, 0,totb*sizeof(int));

    //
    // Put all the atoms into their box
    //
    for (int i=0; i<natoms; i++) {
        const Vector *loc = poshift + i;
        const Vector *force = f+i;
        int axb = (int)(loc->x / pairdistx);
        int ayb = (int)(loc->y / pairdisty);
        int azb = (int)(loc->z / pairdistz);
        aindex = azb * xytotb + ayb * xb + axb;
        if (numinbox[aindex] == 0) {   // First atom in the box
            maxinbox[aindex] = 10;
            boxatom[aindex] = new atominfo[10];
        }
        else if (numinbox[aindex] == maxinbox[aindex]) { // Need to resize the box
            atominfo *tmpbox = new atominfo[2*numinbox[aindex]];
            memcpy((void *)tmpbox, (void *)boxatom[aindex], 
            numinbox[aindex]*sizeof(atominfo));
            delete [] boxatom[aindex];
            boxatom[aindex] = tmpbox;
            maxinbox[aindex] *= 2;
        }   
        boxatom[aindex][numinbox[aindex]].pos = *loc;
        boxatom[aindex][numinbox[aindex]].force = *force;
        boxatom[aindex][numinbox[aindex]].ind = i;
        numinbox[aindex]++;
    } 
    delete [] maxinbox;
}

void Nonbonded::compute(const Initial *init, const Vector *pos,
                               Vector *f, double& Evdw, double& Eelec, const Configure *conf) {

    
    double box_2[3];
    box_2[0] = conf->box[0]*0.5; box_2[1] = conf->box[1]*0.5; box_2[2] = conf->box[2]*0.5;

    Build_cells(init, pos, f, conf);
    Build_Neighborlist(init, pos, f, conf);

    // Loop over cells, and compute the interactions between each cell and
    // its neighbors.

    Evdw = Eelec = 0;
    for (int aindex = 0; aindex<totb; aindex++) {  
        atominfo *tmpbox = boxatom[aindex];
        int *tmpnbr = nbrlist[aindex];
        for (int *nbr = tmpnbr; *nbr != -1; nbr++) {  
            atominfo *nbrbox = boxatom[*nbr];
            for (int i=0; i<numinbox[aindex]; i++) {  
                register Vector tmpf;
                register Vector tmppos = tmpbox[i].pos;
                int ind1 = tmpbox[i].ind;
                Index vdwtype1 = init->atomvdwtype(ind1);
                double kq = COULOMB * init->atomcharge(ind1);
                int startj = 0;
                if (aindex == *nbr) startj = i+1;
                int num = numinbox[*nbr];
                for (int j=startj; j<num; j++) {   
                    Vector dr = nbrbox[j].pos - tmppos;
                    // PBC
                    if (dr.x > box_2[0]) dr.x -= conf->box[0];
                    else if (dr.x <= -box_2[0]) dr.x += conf->box[0];
                    if (dr.y > box_2[1]) dr.y -= conf->box[1];
                    else if (dr.y <= -box_2[1]) dr.y += conf->box[1];
                    if (dr.z > box_2[2]) dr.z -= conf->box[2];
                    else if (dr.z <= -box_2[2]) dr.z += conf->box[2];
                    double dist = dr.length2();
                    if(dist > cut2) continue;   
                    int ind2 = nbrbox[j].ind;
                    if (!init->checkexcl(ind1, ind2)) {  // exclusion check 
                        double r = sqrt(dist);
                        double r_1 = 1.0/r; 
                        double r_2 = r_1*r_1;
                        double r_6 = r_2*r_2*r_2;
                        double r_12 = r_6*r_6;
                        double switchVal = 1, dSwitchVal = 0;
                        if (dist > switch2) {
                            // applying the switch 
                            // see http://localscf.com/localscf.com/LJPotential.aspx.html for details
                            double c2 = cut2 - dist;
                            double c4 = c2*(cut2 + 2*dist - 3.0*switch2);
                            switchVal = c2*c4*c1;
                            dSwitchVal = c3*r*(c2*c2-c4);
                        }

                        // get VDW constants
                        Index vdwtype2 = init->atomvdwtype(ind2);
                        const LJTableEntry *entry;
                        if (init->check14excl(ind1,ind2))
                            entry = ljTable->table_val_scaled14(vdwtype1, vdwtype2); 
                        else
                            entry = ljTable->table_val(vdwtype1, vdwtype2);

                        double vdwA = entry->A;
                        double vdwB = entry->B;
                        double AmBterm = (vdwA * r_6 - vdwB)*r_6;
                        Evdw += switchVal*AmBterm;
                        double force_r = ( switchVal * 6.0 * (vdwA*r_12 + AmBterm) *
                        r_1 - AmBterm*dSwitchVal )*r_1;
             
                        // Electrostatics
                        double kqq = kq * init->atomcharge(ind2);
                        double efac = 1.0-dist/cut2;
                        double prefac = kqq * r_1 * efac;
                        Eelec += prefac * efac;
                        force_r += prefac * r_1 * (r_1 + 3.0*r/cut2);
           
                        tmpf -= force_r * dr; 
                        nbrbox[j].force += force_r * dr;

                    } // exclusion check 
                }  
                tmpbox[i].force += tmpf; 
            }     
        }       
    }       

    // 
    // copy forces from atomboxes to the f array
    //
    for (int i = 0; i < totb; i++) {
        for (int j=0; j<numinbox[i]; j++) {
        f[boxatom[i][j].ind] = boxatom[i][j].force;
        }
    }
  
    // free up the storage space allocted for the grid search
    for(int i=0; i < totb; i++) {
        if (numinbox[i])  delete [] boxatom[i];
        delete [] nbrlist[i];
    }
    delete [] nbrlist;
    delete [] boxatom;
    delete [] numinbox;
}
