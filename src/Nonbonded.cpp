#include "Nonbonded.h"
#include "Parameters.h"
#include "Initial.h"
#include "LJTable.h"
#include <math.h>
#define PI 3.14159265

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
    
    atoms = new NonbondedAtom[natoms];
    memset((void *)atoms, 0, natoms*sizeof(NonbondedAtom));

    poshift = new Vector[natoms];
    memset((void *)poshift, 0, natoms*sizeof(Vector));

    //three-body cutoff
    threebody_cut2 = conf->threebodyCutoff;
    threebody_cut2 *= threebody_cut2;
    threebody_pair2 = conf->threebodyPairDist;
    threebody_pair2 *= threebody_pair2;
    threebody_ijcut2 = conf->threebodyIJCutoff;
    threebody_ijcut2 *= threebody_ijcut2; 
}

void Nonbonded::build_cells(const Initial *init, const Vector *pos, Vector *f, const Configure *conf) {


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

void Nonbonded::build_mycells(const Initial *init, const Vector *pos, Vector *f, const Configure *conf) {

    int  ncellxy, tcells;
    double pairdist = conf->pairlistdist;

    ncellx = int(conf->box[0]/pairdist) + 1;
    ncelly = int(conf->box[1]/pairdist) + 1;
    ncellz = int(conf->box[2]/pairdist) + 1;
    ncellxy = ncelly * ncellx;
    tcells = ncellxy * ncellz;
    cells = new Cell[tcells];  // Allocate space for total number of cells
    memset((void *)cells, 0, tcells*sizeof(Cell));

        // Create list of neighbor for each box
    int index, nbindex;
    int xnb, ynb, znb;
    int cellindex = 0;
    for (int zi=0; zi<ncellz; zi++) {
        for (int yi=0; yi<ncelly; yi++) {
            for (int xi=0; xi<ncellx; xi++) {
                int nbrs[27];
                int nnb = 0;
                index = zi * ncellxy + yi * ncellx + xi;
                for (int xx=-1; xx<2; xx++){
                    for (int yy=-1; yy<2; yy++){
                        for (int zz=-1; zz<2; zz++){
                            znb = zi + zz;
                            if (znb < 0 ) znb = ncellz -1;
                            else if (znb == ncellz) znb = 0;

                            ynb = yi + yy;
                            if (ynb < 0 ) ynb = ncelly -1;
                            else if (ynb == ncelly) ynb = 0;

                            xnb = xi + xx;
                            if (xnb < 0 ) xnb = ncellx -1;
                            else if (xnb == ncellx) xnb = 0;

                            nbindex = znb * ncellxy + ynb * ncellx + xnb;
                            nbrs[nnb++] = nbindex;

                        }
                    }
                }
            // Copy 14 elemnets of nbrs to cells and delete double elements
            for (int jj=0; jj<14; jj++){ cells[cellindex].nbrlist[jj] = -1; }
            cells[cellindex].num = 0;
            int jj=0;
            bool flag = true;
            for (int ii=0; ii<27; ii++){
                while (cells[cellindex].nbrlist[jj] != -1 && jj<14){
                    if (nbrs[ii] == cells[cellindex].nbrlist[jj]) flag = false;
                    jj++;
                }
                if (flag) {
                    cells[cellindex].nbrlist[cells[cellindex].num++] = nbrs[ii];
                }
                jj=0;
                flag = true;
                if (cells[cellindex].num == 14 || cells[cellindex].nbrlist[cells[cellindex].num - 1] == cellindex) break;
            }
            cellindex++;
            }
        }
    }
}

void Nonbonded::build_atomlist(const Initial *init, const Vector *pos, Vector *f, const Configure *conf){
    
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

void Nonbonded::build_neighborlist(const Initial *init, const Vector *pos, Vector *f, const Configure *conf){
    
    double box_2[3];
    box_2[0] = conf->box[0]*0.5; box_2[1] = conf->box[1]*0.5; box_2[2] = conf->box[2]*0.5;

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
    double pairdistx = conf->box[0]/ncellx + 0.001; 
    double pairdisty = conf->box[1]/ncelly + 0.001; 
    double pairdistz = conf->box[2]/ncellz + 0.001;

    int ncellxy = ncelly * ncellx;
    int tcell   = ncellxy * ncellz;

    for (int ii=0; ii<tcell; ii++) cells[ii].atoms.clear(); // Clear memory to save atoms again
    int cellx, celly, cellz, index;

    for (int ii=0; ii<natoms; ii++){
        const Vector *loc = poshift + ii;
        cellx = int(loc->x/pairdistx);
        celly = int(loc->y/pairdisty);
        cellz = int(loc->z/pairdistz);
        index = cellz * ncellxy + celly * ncellx + cellx;
        cells[index].atoms.push_back(ii);
    }

    for (int iatom=0; iatom<natoms; iatom++) {
        atoms[iatom].nbrlist1.clear(); // Clear the neighbor list
        atoms[iatom].nbrlist2.clear();
    }

    for (int icell=0; icell<tcell; icell++) {
        for (int ii=0; ii<cells[icell].atoms.size(); ii++) {
            int iatom = cells[icell].atoms[ii];
            for (int jj=ii+1; jj<cells[icell].atoms.size(); jj++){ // Atoms i and j are both in one cell
                int jatom = cells[icell].atoms[jj];
                Vector dij = poshift[iatom] - poshift[jatom];
                // PBC
                if (dij.x>box_2[0]) dij.x -= conf->box[0];
                else if (dij.x<=-box_2[0]) dij.x += conf->box[0];
                if (dij.y>box_2[1]) dij.y -= conf->box[1];
                else if (dij.y<=-box_2[1]) dij.y += conf->box[1];
                if (dij.z>box_2[2]) dij.z -= conf->box[2];
                else if (dij.z<=-box_2[2]) dij.z += conf->box[2];
                double dist2 = dij.length2();

                if (dist2 < pair2){
                    if (!init->checkexcl(iatom, jatom)){
                        if (dist2 > cut2) {
                            atoms[iatom].nbrlist1.push_back(jatom);
                        } else {
                            atoms[iatom].nbrlist2.push_back(jatom);
                        }
                    }
                }
            }
            // Atom i is in the cell icell and atom j is in the cell jcell
            for (int nn=0; nn<cells[icell].num-1; nn++){ // Search nearest neighbors; Maximum 13 neighbors;
                int jcell=cells[icell].nbrlist[nn];
                for (int jj=0; jj<cells[jcell].atoms.size(); jj++) { // Loop over all atoms in the jcell
                    int jatom = cells[jcell].atoms[jj];
                    Vector dij = poshift[iatom] - poshift[jatom];
                    // PBC
                    if (dij.x>box_2[0]) dij.x -= conf->box[0];
                    else if (dij.x<=-box_2[0]) dij.x += conf->box[0];
                    if (dij.y>box_2[1]) dij.y -= conf->box[1];
                    else if (dij.y<=-box_2[1]) dij.y += conf->box[1];
                    if (dij.z>box_2[2]) dij.z -= conf->box[2];
                    else if (dij.z<=-box_2[2]) dij.z += conf->box[2];
                    double dist2 = dij.length2();

                    if (dist2 < pair2){
                        if (!init->checkexcl(iatom, jatom)){
                            if (dist2 > cut2) {
                                atoms[iatom].nbrlist1.push_back(jatom);
                            } else {
                                atoms[iatom].nbrlist2.push_back(jatom);
                            }
                        }
                    }
                }
            }

            atoms[iatom].nbrlist1.shrink_to_fit(); atoms[iatom].nbrlist1.reserve(atoms[iatom].nbrlist1.size()+10);
            atoms[iatom].nbrlist2.shrink_to_fit(); atoms[iatom].nbrlist2.reserve(atoms[iatom].nbrlist2.size()+10);

        }
    }
}

void Nonbonded::threebody_neighborlist(const Initial *init, const Vector *pos, Vector *f, const Configure *conf){
    /**
        Makes the neighborlist for three-body calculations.
        The difference between this function and a general neighborlist function
        is that we don't check for exclusion list, meaning that we include 
        bonded atom in the neighborlist.
    */
    double box_2[3];
    box_2[0] = conf->box[0]*0.5; box_2[1] = conf->box[1]*0.5; box_2[2] = conf->box[2]*0.5;

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
    double pairdistx = conf->box[0]/ncellx + 0.001; 
    double pairdisty = conf->box[1]/ncelly + 0.001; 
    double pairdistz = conf->box[2]/ncellz + 0.001;

    int ncellxy = ncelly * ncellx;
    int tcell   = ncellxy * ncellz;

    for (int ii=0; ii<tcell; ii++) cells[ii].atoms.clear(); // Clear memory to save atoms again
    int cellx, celly, cellz, index;

    for (int ii=0; ii<natoms; ii++){
        const Vector *loc = poshift + ii;
        cellx = int(loc->x/pairdistx);
        celly = int(loc->y/pairdisty);
        cellz = int(loc->z/pairdistz);
        index = cellz * ncellxy + celly * ncellx + cellx;
        cells[index].atoms.push_back(ii);
    }

    for (int iatom=0; iatom<natoms; iatom++) {
        atoms[iatom].nbrlist1.clear(); // Clear the neighbor list
        atoms[iatom].nbrlist2.clear();
    }

    for (int icell=0; icell<tcell; icell++) {
        for (int ii=0; ii<cells[icell].atoms.size(); ii++) {
            int iatom = cells[icell].atoms[ii];
            for (int jj=ii+1; jj<cells[icell].atoms.size(); jj++){ // Atoms i and j are both in one cell
                int jatom = cells[icell].atoms[jj];
                Vector dij = poshift[iatom] - poshift[jatom];
                // PBC
                if (dij.x>box_2[0]) dij.x -= conf->box[0];
                else if (dij.x<=-box_2[0]) dij.x += conf->box[0];
                if (dij.y>box_2[1]) dij.y -= conf->box[1];
                else if (dij.y<=-box_2[1]) dij.y += conf->box[1];
                if (dij.z>box_2[2]) dij.z -= conf->box[2];
                else if (dij.z<=-box_2[2]) dij.z += conf->box[2];
                double dist2 = dij.length2();

                if (dist2 < pair2){
                    if (dist2 > cut2) {
                        atoms[iatom].nbrlist1.push_back(jatom);
                    } else {
                        atoms[iatom].nbrlist2.push_back(jatom);
                        }
                }
            }
            // Atom i is in the cell icell and atom j is in the cell jcell
            for (int nn=0; nn<cells[icell].num-1; nn++){ // Search nearest neighbors; Maximum 13 neighbors;
                int jcell=cells[icell].nbrlist[nn];
                for (int jj=0; jj<cells[jcell].atoms.size(); jj++) { // Loop over all atoms in the jcell
                    int jatom = cells[jcell].atoms[jj];
                    Vector dij = poshift[iatom] - poshift[jatom];
                    // PBC
                    if (dij.x>box_2[0]) dij.x -= conf->box[0];
                    else if (dij.x<=-box_2[0]) dij.x += conf->box[0];
                    if (dij.y>box_2[1]) dij.y -= conf->box[1];
                    else if (dij.y<=-box_2[1]) dij.y += conf->box[1];
                    if (dij.z>box_2[2]) dij.z -= conf->box[2];
                    else if (dij.z<=-box_2[2]) dij.z += conf->box[2];
                    double dist2 = dij.length2();

                    if (dist2 < pair2){
                        if (dist2 > cut2) {
                            atoms[iatom].nbrlist1.push_back(jatom);
                        } else {
                            atoms[iatom].nbrlist2.push_back(jatom);
                        }
                    }
                }
            }

            atoms[iatom].nbrlist1.shrink_to_fit(); atoms[iatom].nbrlist1.reserve(atoms[iatom].nbrlist1.size()+10);
            atoms[iatom].nbrlist2.shrink_to_fit(); atoms[iatom].nbrlist2.reserve(atoms[iatom].nbrlist2.size()+10);

        }
    }
}

void Nonbonded::compute(const Initial *init, const Vector *pos,
                               Vector *f, double& Evdw, double& Eelec, const Configure *conf) {

    
    double box_2[3];
    box_2[0] = conf->box[0]*0.5; box_2[1] = conf->box[1]*0.5; box_2[2] = conf->box[2]*0.5;

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

                    }
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
}

void Nonbonded::mycompute(const Initial *init, const Vector *pos,
                               Vector *f, double& Evdw, double& Eelec, const Configure *conf) {

    
    double box_2[3];
    box_2[0] = conf->box[0]*0.5; box_2[1] = conf->box[1]*0.5; box_2[2] = conf->box[2]*0.5;
    Vector loci;

    Evdw = Eelec = 0;
    // loop over all atoms
    for (int iatom=0; iatom<natoms; iatom++){

        // get the vdwtype and charge of the atom i
        Index vdwtype1 = init->atomvdwtype(iatom);
        double kq = COULOMB * init->atomcharge(iatom);
        //get the coordinate of the atom 
        loci = pos[iatom];
        
        // get number of neighbors between the cutoff and pairlistdist
        int nnbr1 = atoms[iatom].nbrlist1.size();
        int i = 0;
        for (int jatom : atoms[iatom].nbrlist1) {
            // get the vector connecting i to j
            Vector dij = loci - pos[jatom];
            //PBC
            if (dij.x > box_2[0]) dij.x -= conf->box[0];
            else if (dij.x<=-box_2[0]) dij.x += conf->box[0];
            if (dij.y > box_2[1]) dij.y -= conf->box[1];
            else if (dij.y <= -box_2[1]) dij.y += conf->box[1];
            if (dij.z > box_2[2]) dij.z -= conf->box[2];
            else if (dij.z <= -box_2[2]) dij.z += conf->box[2];

            // get the squared distance
            double dr2 = dij.length2();

            // check if the j atom moved to the cutoff distance
            if (dr2 < cut2) {
                // if so add the atom to the inner neighborlist
                int pushf = atoms[iatom].nbrlist1[jatom];
                atoms[iatom].nbrlist2.push_back(pushf);
                //now erase that atom from the outer list
                // check this
                atoms[iatom].nbrlist1.erase(atoms[iatom].nbrlist1.begin()+i);
            }
            // check if the atom moved out of the outer cutoff
            else if (dr2 > pair2) {
                // if so remove it from the outer list
                atoms[iatom].nbrlist1.erase(atoms[iatom].nbrlist1.begin()+i);
            }
            i++;
        }
            
        // resize the vector
        atoms[iatom].nbrlist2.shrink_to_fit(); atoms[iatom].nbrlist2.reserve(atoms[iatom].nbrlist2.size()+10);
        // get the number of neighbors within the cutoff 
        int nnbr2 = atoms[iatom].nbrlist2.size();
        // a vector to hold forces on atom i from the neighbors
        register Vector tmpf;
        // loop over all atoms within the cutoff
        for (int jatom : atoms[iatom].nbrlist2) {
            
            // get the vector connecting atom i and j
            Vector dij = loci - pos[jatom];
            //PBC
            if (dij.x > box_2[0]) dij.x -= conf->box[0];
            else if (dij.x <= -box_2[0]) dij.x += conf->box[0];
            if (dij.y > box_2[1]) dij.y -= conf->box[1];
            else if (dij.y <= -box_2[1]) dij.y += conf->box[1];
            if (dij.z > box_2[2]) dij.z -= conf->box[2];
            else if (dij.z <= -box_2[2]) dij.z += conf->box[2];

            // get the squared distance between atom i and j
            double dr2 = dij.length2();
            // if the distance if larger than cutoff skip atom j elae continue
            if(dr2 > cut2) continue;

            // comoute the energy and apply the switch
            double r = sqrt(dr2);
            double r_1 = 1.0/r; 
            double r_2 = r_1*r_1;
            double r_6 = r_2*r_2*r_2;
            double r_12 = r_6*r_6;
            double switchVal = 1, dSwitchVal = 0;
            if (dr2 > switch2) {
                // applying the switch 
                // see http://localscf.com/localscf.com/LJPotential.aspx.html for details
                double c2 = cut2 - dr2;
                double c4 = c2*(cut2 + 2*dr2 - 3.0*switch2);
                switchVal = c2*c4*c1;
                dSwitchVal = c3*r*(c2*c2-c4);
            }

            // get VDW constants
            Index vdwtype2 = init->atomvdwtype(jatom);
            const LJTableEntry *entry;
            if (init->check14excl(iatom,jatom))
                entry = ljTable->table_val_scaled14(vdwtype1, vdwtype2); 
            else
                entry = ljTable->table_val(vdwtype1, vdwtype2);

            double vdwA = entry->A;
            double vdwB = entry->B;
            double AmBterm = (vdwA * r_6 - vdwB)*r_6;
            Evdw += switchVal*AmBterm;
            double lj_force = ( switchVal * 6.0 * (vdwA*r_12 + AmBterm) *
            r_1 - AmBterm*dSwitchVal )*r_1;
             
            // Electrostatics
            double kqq = kq * init->atomcharge(jatom);
            double efac = 1.0 - dr2/cut2;
            double prefac = kqq * r_1 * efac;
            Eelec += prefac * efac;
            double elec_force = prefac * r_1 * (r_1 + 3.0*r/cut2);
            f[jatom] -= (lj_force + elec_force) * dij;
            tmpf -= f[jatom]; 
        }
        f[iatom] += tmpf;
    }
}

void Nonbonded::compute_threebody(const Initial *init, const Vector *pos,
                               Vector *f, double& Emisc, const Configure *conf) {


    double box_2[3];
    box_2[0] = conf->box[0]*0.5; box_2[1] = conf->box[1]*0.5; box_2[2] = conf->box[2]*0.5;
    Vector loci;

    Emisc = 0;
    // loop over all atoms
    for (int iatom=0; iatom<natoms; iatom++){

        //get the coordinate of the atom 
        loci = pos[iatom];
        // get number of neighbors between the cutoff and pairlistdist
        int nnbr1 = atoms[iatom].nbrlist1.size();
        int i = 0;
        for (int jatom : atoms[iatom].nbrlist1) {
            // get the vector connecting i to j
            Vector dij = loci - pos[jatom];
            //PBC
            if (dij.x > box_2[0]) dij.x -= conf->box[0];
            else if (dij.x<=-box_2[0]) dij.x += conf->box[0];
            if (dij.y > box_2[1]) dij.y -= conf->box[1];
            else if (dij.y <= -box_2[1]) dij.y += conf->box[1];
            if (dij.z > box_2[2]) dij.z -= conf->box[2];
            else if (dij.z <= -box_2[2]) dij.z += conf->box[2];

            // get the squared distance
            double dr2 = dij.length2();

            // check if the j atom moved to the cutoff distance
            if (dr2 < threebody_cut2) {
                // if so add the atom to the inner neighborlist
                int pushf = atoms[iatom].nbrlist1[jatom];
                atoms[iatom].nbrlist2.push_back(pushf);
                //now erase that atom from the outer list
                // check this
                atoms[iatom].nbrlist1.erase(atoms[iatom].nbrlist1.begin()+i);
            }
            // check if the atom moved out of the outer cutoff
            else if (dr2 > pair2) {
                // if so remove it from the outer list
                atoms[iatom].nbrlist1.erase(atoms[iatom].nbrlist1.begin()+i);
            }
            i++;
        }
            
        // resize the vector
        atoms[iatom].nbrlist2.shrink_to_fit(); atoms[iatom].nbrlist2.reserve(atoms[iatom].nbrlist2.size()+10);
        // get the number of neighbors within the cutoff 
        int nnbr2 = atoms[iatom].nbrlist2.size();

        // loop over all atoms within the cutoff
        for (int j = 0; j < atoms[iatom].nbrlist2.size(); j++) {
            
            int jatom = atoms[iatom].nbrlist2[j];
            Vector dij = loci - pos[jatom];
            
            //PBC
            apply_pbc(conf->box, box_2, dij);

            double distij = dij.length2();
            // look for k if i and j are bonded
            if (distij < threebody_ijcut2) {
                for (int k = j+1; k < atoms[iatom].nbrlist2.size(); k++) {
                    
                    int katom = atoms[iatom].nbrlist2[k];
                    // check if k and i are different chains
                    int iatom_res = init->get_resnum(iatom);
                    int katom_res = init->get_resnum(katom);

                    if (iatom_res != katom_res) {

                        // get the vector connecting atom i, j, and k
                        Vector dik = loci - pos[katom];
                        Vector djk = pos[jatom] - pos[katom];

                        // PBC
                        apply_pbc(conf->box, box_2, dik);
                        apply_pbc(conf->box, box_2, djk);

                        // get the squared distance between atom i and j
                        double distik = dik.length2();
                        double distjk = djk.length2();

                        // if the distance if larger than cutoff skip
                        if(distik > threebody_cut2 || distjk > threebody_cut2) continue;   

                        double rij = sqrt(distij);
                        double rik = sqrt(distik);
                        double rjk = sqrt(distjk);

                        double rij_1 = 1.0/rij; 
                        double rij_3 = rij_1*rij_1*rij_1;
                        double rik_1 = 1.0/rik; 
                        double rik_3 = rik_1*rik_1*rik_1;
                        double rjk_1 = 1.0/rjk; 
                        double rjk_3 = rjk_1*rjk_1*rjk_1;

                        double cos_thetaijk = (dij*dik)/(rij*rik);
                        double cos_thetajki = -(djk*dij)/(rjk*rij);
                        double cos_thetakij = (dik*djk)/(rik*rjk);

                        double thetaijk = acos(cos_thetaijk)* 180.0 / PI;
                        double thetajki = acos(cos_thetajki)* 180.0 / PI;
                        double thetakij = acos(cos_thetakij)* 180.0 / PI;
                        //cout << thetaijk << " " << thetajki << " " << thetakij << endl;
                        double three_body_ene = (1 + 3 * cos_thetaijk * cos_thetajki * cos_thetakij)*rij_3*rik_3*rjk_3;
                        Emisc += three_body_ene;
                    }
                }
            }
        }
    }
}

void Nonbonded::apply_pbc(const double box[3], const double box_2[3], Vector &dij){
    
    if (dij.x > box_2[0]) dij.x -= box[0];
    else if (dij.x <= -box_2[0]) dij.x += box[0];
    if (dij.y > box_2[1]) dij.y -= box[1];
    else if (dij.y <= -box_2[1]) dij.y += box[1];
    if (dij.z > box_2[2]) dij.z -= box[2];
    else if (dij.z <= -box_2[2]) dij.z += box[2];
}


Nonbonded::~Nonbonded() {
    delete ljTable;
    delete[] atoms;
    delete[] poshift;
    // free up the storage space allocted for the grid search
    if (totb > 0) {
        for(int i=0; i < totb; i++) {
            if (numinbox[i])  delete [] boxatom[i];
            delete [] nbrlist[i];
        }
        delete [] nbrlist;
        delete [] boxatom;
        delete [] numinbox;
    }
}
