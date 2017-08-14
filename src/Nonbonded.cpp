#include "Nonbonded.h"
#include "Parameters.h"
#include "Initial.h"
#include "LJTable.h"

struct atominfo {
  Vector pos;
  Vector force;
  int ind;
};

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
}


Nonbonded::~Nonbonded() {
    delete ljTable;
}

void Nonbonded::compute(const Initial *init, const Vector *pos,
                               Vector *f, double& Evdw, double& Eelec) {

    int xb, yb, zb, xytotb, totb;               // dimensions of decomposition

    atominfo **boxatom;       // positions, forces, and indicies for each atom  
    int *numinbox, *maxinbox; // Number of atoms in each box
    int **nbrlist;            // List of neighbors for each box

    int i, j, aindex;

    //
    // find min/max bounds of Initial's coordinates
    //

    double xmin, xmax, ymin, ymax, zmin, zmax;  // extent of atomic coordinates
    for(i=0; i < natoms; i++) {
        const Vector *loc = pos+i;
        if(i==0) {
            xmin = xmax = loc->x; ymin = ymax = loc->y; zmin = zmax = loc->z;
        } 
        else {
            if(loc->x < xmin) xmin = loc->x;
            else if(loc->x > xmax) xmax = loc->x;
            if(loc->y < ymin) ymin = loc->y;
            else if(loc->y > ymax) ymax = loc->y;
            if(loc->z < zmin) zmin = loc->z;
            else if(loc->z > zmax) zmax = loc->z;
        }
    }
 
    // from size of Initial, break up space into boxes of dimension pairlistdist
    // Since I'm recreating the pairlist each time, there's no need to make the
    // boxes any bigger than the cutoff length.

    float pairdist = sqrt(cut2);
    xb = (int)((xmax - xmin) / pairdist) + 1;
    yb = (int)((ymax - ymin) / pairdist) + 1;
    zb = (int)((zmax - zmin) / pairdist) + 1;
    xytotb = yb * xb;
    totb = xytotb * zb;

    boxatom = new atominfo*[totb];
    nbrlist = new int *[totb];
    numinbox = new int[totb];
    maxinbox = new int[totb];
    memset((void *)numinbox, 0,totb*sizeof(int));
    memset((void *)maxinbox, 0,totb*sizeof(int));

    //
    // Put all the atoms into their box
    //
    for (i=0; i<natoms; i++) {
        const Vector *loc = pos+i;
        const Vector *force = f+i;
        int axb = (int)((loc->x - xmin) / pairdist);
        int ayb = (int)((loc->y - ymin) / pairdist);
        int azb = (int)((loc->z - zmin) / pairdist);
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

    //
    // Create neighbor list for each box
    //
    aindex = 0;
    for (int zi=0; zi<zb; zi++) {
        for (int yi=0; yi<yb; yi++) {
            for (int xi=0; xi<xb; xi++) {
                int nbrs[14];           // Max possible number of neighbors in 3D
                int n=0;                // Number of neighbors found so far
                nbrs[n++] = aindex;     // Always include self
                if (xi < xb-1) nbrs[n++] = aindex + 1;
                if (yi < yb-1) nbrs[n++] = aindex + xb;
                if (zi < zb-1) nbrs[n++] = aindex + xytotb;
                if (xi < (xb-1) && yi < (yb-1)) nbrs[n++] = aindex + xb + 1;
                if (xi < (xb-1) && zi < (zb-1)) nbrs[n++] = aindex + xytotb + 1;
                if (yi < (yb-1) && zi < (zb-1)) nbrs[n++] = aindex + xytotb + xb;
                if (xi < (xb-1) && yi > 0)      nbrs[n++] = aindex - xb + 1;
                if (xi > 0 && zi < (zb-1))     nbrs[n++] = aindex + xytotb - 1;
                if (yi > 0 && zi < (zb-1))     nbrs[n++] = aindex + xytotb - xb;
                if (xi < (xb-1) && yi < (yb-1) && zi < (zb-1))
                    nbrs[n++] = aindex + xytotb + xb + 1;
                if (xi > 0 && yi < (yb-1) && zi < (zb-1))
                    nbrs[n++] = aindex + xytotb + xb - 1; 
                if (xi < (xb-1) && yi > 0 && zi < (zb-1))
                    nbrs[n++] = aindex + xytotb - xb + 1;
                if (xi > 0 && yi > 0 && zi < (zb-1))
                    nbrs[n++] = aindex + xytotb - xb - 1;
        
                nbrlist[aindex] = new int[n+1];
                memcpy((void *)nbrlist[aindex], (void *)nbrs, n*sizeof(int));
                nbrlist[aindex][n] = -1;  // Sentinel for end of neighbors
                aindex++;
            }
        }
    }

    //
    // Loop over boxes, and compute the interactions between each box and
    // its neighbors.
    //

    Evdw = Eelec = 0;
    for (aindex = 0; aindex<totb; aindex++) {
        atominfo *tmpbox = boxatom[aindex];
        int *tmpnbr = nbrlist[aindex];
        for (int *nbr = tmpnbr; *nbr != -1; nbr++) {
            atominfo *nbrbox = boxatom[*nbr];

            for (i=0; i<numinbox[aindex]; i++) {
                register Vector tmpf;
                register Vector tmppos = tmpbox[i].pos;
                int ind1 = tmpbox[i].ind;
                Index vdwtype1 = init->atomvdwtype(ind1);
                double kq = COULOMB * init->atomcharge(ind1);
                int startj = 0;
                if (aindex == *nbr) startj = i+1;
                int num = numinbox[*nbr];
                for (j=startj; j<num; j++) {
                    Vector dr = nbrbox[j].pos - tmppos;
                    double dist = dr.length2();
                    if(dist > cut2) continue;   
                    int ind2 = nbrbox[j].ind;
                    if (!init->checkexcl(ind1, ind2)) {
                        double r = sqrt(dist);
                        double r_1 = 1.0/r; 
                        double r_2 = r_1*r_1;
                        double r_6 = r_2*r_2*r_2;
                        double r_12 = r_6*r_6;
                        double switchVal = 1, dSwitchVal = 0;
                        if (dist > switch2) {
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
                }   // Loop over neighbor atoms
                tmpbox[i].force += tmpf; 
            }     // Loop over self atoms
        }       // Loop over neighbor boxes
    }         // Loop over self boxes

    // 
    // copy forces from atomboxes to the f array
    //
    for (i = 0; i < totb; i++) {
        for (j=0; j<numinbox[i]; j++) {
        f[boxatom[i][j].ind] = boxatom[i][j].force;
        }
    }
  
    // free up the storage space allocted for the grid search
    for(i=0; i < totb; i++) {
        if (numinbox[i])  delete [] boxatom[i];
        delete [] nbrlist[i];
    }
    delete [] nbrlist;
    delete [] boxatom;
    delete [] numinbox;
}
