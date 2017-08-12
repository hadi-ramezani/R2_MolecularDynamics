/*
 * File:   Nonbonded.cpp
 * Author: amin
 *
 * Created on September 13, 2015, 4:21 PM
 */
#define COULOMB 332.0636 
#define ALIGN __attribute__((aligned(32)))
#define P4 float
#define L1 1024
#define L2 2048
#define nint(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))


#include "Nonbonded.h"
#include "Integrator.h"
#include "Vector.h"
#include <vector>
#include <algorithm>
#include <numeric>

using namespace std;

Nonbonded::Nonbonded(int num, const double rc, const double switch1, const double pairlistdist, int seed) {
    atoms = new NonbondedAtom[num];
    memset((void *)atoms, 0, num*sizeof(NonbondedAtom));

    poshift = new Vector[num];
    memset((void *)poshift, 0, num*sizeof(Vector));

    cut = rc;
    cut2 = rc*rc;
    pairdist2 = pairlistdist*pairlistdist;
    switch2 = switch1*switch1;

    cut_1 = 1.0/rc;
    cut_2 = 1.0f/cut2;
    switch_2 = 1.0f/switch2;
    pairdist_2 = 1.0f/pairdist2;

    c1 = 1.0/(cut2-switch2);
    c1 = c1*c1*c1;
    c3 = 4*c1;

    // Initializing random numbers
    vslNewStream( &stream, VSL_BRNG_SFMT19937, seed );
}

void Nonbonded::built_table(const Initial *init, const int ntype, const int num){
    rmin2_table = new float*[ntype];
    epsi_table = new float*[ntype];
    kqq_table = new float*[ntype];
    for (int ii=0; ii<ntype; ii++){
        rmin2_table[ii] = new float[ntype];
        epsi_table[ii] = new float[ntype];
        kqq_table[ii] = new float[ntype];
    }

    for (int iatom=0; iatom<ntype; iatom++){
        for (int jatom=0; jatom<ntype; jatom++){
            rmin2_table[iatom][jatom] = ((init->vdw[iatom].rmin + init->vdw[jatom].rmin)/2.0)*((init->vdw[iatom].rmin + init->vdw[jatom].rmin)/2.0);
            epsi_table[iatom][jatom] = sqrt(init->vdw[iatom].epsi * init->vdw[jatom].epsi);
            kqq_table[iatom][jatom] = COULOMB * init->vdw[iatom].charge * init->vdw[jatom].charge;
        }
    }

    type = new int[num];
    for (int iatom=0; iatom<num; iatom++){
        atoms[iatom].type = init->vdw[iatom].type;
        type[iatom] = init->vdw[iatom].type;
    }
}
// Divide Box into small cells and put atoms into cells
void Nonbonded::Build_cells(const double* box, const double celldist, const int num){

    int  ncellxy, tcells;
    ncellx = int(box[0]/celldist) + 1;
    ncelly = int(box[1]/celldist) + 1;
    ncellz = int(box[2]/celldist) + 1;
    ncellxy = ncelly * ncellx;
    tcells = ncellxy * ncellz;
    cells = new Cell[tcells];  // Allocate space for total number of cells
    memset((void *)cells, 0, tcells*sizeof(Cell));
    cout << "The simulation box is divided into the " << tcells << " cells for neighbor list calculation." << endl;
    cout << "This works only for NVT and NPT with small volume fluctuation." << endl;

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

void Nonbonded::Neighborlist(const double* box, const int num, const Initial *init, const Vector *pos){

    double box_2[3];
    box_2[0] = box[0]*0.5; box_2[1] = box[1]*0.5; box_2[2] = box[2]*0.5;
    
    // Shift all atoms into the simulation box
    // The origin of the box could differs from NAMD. This is important when we use NAMD coordinates
    for (int ii=0; ii<num; ii++){ // For a box with the origin at the lower left vertex
        poshift[ii].x = pos[ii].x - floor(pos[ii].x/box[0])*box[0];
        poshift[ii].y = pos[ii].y - floor(pos[ii].y/box[1])*box[1];
        poshift[ii].z = pos[ii].z - floor(pos[ii].z/box[2])*box[2];
    }

    double pairdistx = box[0]/ncellx + 0.001; // add small number to make sure all particles are located in the cells
    double pairdisty = box[1]/ncelly + 0.001; // avoid the problem for particles located on the border
    double pairdistz = box[2]/ncellz + 0.001;
    int ncellxy = ncelly * ncellx;
    int tcell   = ncellxy * ncellz;

    for (int ii=0; ii<tcell; ii++) cells[ii].atoms.clear(); // Clear memory to save atoms again
    int cellx, celly, cellz, index;
    for (int ii=0; ii<num; ii++){
        const Vector *loc = poshift + ii;
        cellx = int(loc->x/pairdistx);
        celly = int(loc->y/pairdisty);
        cellz = int(loc->z/pairdistz);
        index = cellz * ncellxy + celly * ncellx + cellx;
        cells[index].atoms.push_back(ii);
    }

    for (int iatom=0; iatom<num; iatom++) {
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
                if (dij.x>box_2[0]) dij.x -= box[0];
                else if (dij.x<=-box_2[0]) dij.x += box[0];
                if (dij.y>box_2[1]) dij.y -= box[1];
                else if (dij.y<=-box_2[1]) dij.y += box[1];
                if (dij.z>box_2[2]) dij.z -= box[2];
                else if (dij.z<=-box_2[2]) dij.z += box[2];
                double dist2 = dij.length2();

                if (dist2 < pairdist2){
                    if (!init->exclude(iatom,jatom)){
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
                    if (dij.x>box_2[0]) dij.x -= box[0];
                    else if (dij.x<=-box_2[0]) dij.x += box[0];
                    if (dij.y>box_2[1]) dij.y -= box[1];
                    else if (dij.y<=-box_2[1]) dij.y += box[1];
                    if (dij.z>box_2[2]) dij.z -= box[2];
                    else if (dij.z<=-box_2[2]) dij.z += box[2];
                    double dist2 = dij.length2();

                    if (dist2 < pairdist2){
                        if (!init->exclude(iatom,jatom)){
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

            if (atoms[iatom].nbrlist2.size()>=L2) {
                cout << "Atom " << iatom << " has " << atoms[iatom].nbrlist2.size() << " neighbors." << endl;;
                cout << "ERROR: The number of atoms in the neighbor list is more than " << L2 << " atoms." << endl;
                cout << "Reduce pairlistdist" << endl;
                exit(1);
            }
        }
    }
}
void Nonbonded::Compute(const Initial *init, const Vector *pos, Vector *const ff, const int num, double &Evdw, double &Eelec){

int ntype = init->ntype;
double rmin2[ntype], epsi[ntype], kqq[ntype];
Vector box, box_2;
box.x = init->box[0]; box.y=init->box[1]; box.z=init->box[2];
box_2 = box*0.5;

Vector loci;
Vector dij[L2];
Vector FF[L2] ;
double dr2[L2], dr_1[L2], dr_2[L2], dr_6[L2], r2_d2[L2] ALIGN;
double ljE[L2], ljF[L2] ALIGN;
double elecE[L2], elecF[L2] ;
double dpdF[L2];
double rmin2_ij[L2], epsi_ij[L2], kqq_ij[L2] ALIGN;
double dot_rv[L2];
int jtypes[L2];
int eraselist[20];

Evdw = 0.0; Eelec = 0.0;
for (int iatom=0; iatom<num; iatom++){
    int itype = atoms[iatom].type;
    for (int ii=0; ii<ntype; ii++){
        rmin2[ii]=rmin2_table[itype][ii];
        epsi[ii]=epsi_table[itype][ii];
        kqq[ii]=kqq_table[itype][ii];
    }

    loci = pos[iatom];

    int nnbr1 = atoms[iatom].nbrlist1.size();
    int index=0;
//#pragma simd
    for (int jatom : atoms[iatom].nbrlist1) {
        dij[index] = loci - pos[jatom];
        index++;
    }

#pragma simd
    for (index=0; index<nnbr1; index++){

        if (dij[index].x>box_2.x) dij[index].x -= box.x;
        else if (dij[index].x<=-box_2.x) dij[index].x += box.x;
        if (dij[index].y>box_2.y) dij[index].y -= box.y;
        else if (dij[index].y<=-box_2.y) dij[index].y += box.y;
        if (dij[index].z>box_2.z) dij[index].z -= box.z;
        else if (dij[index].z<=-box_2.z) dij[index].z += box.z;

        dr2[index] = dij[index].length2();
    }

    int nerase = 0;
    for (index=0; index<nnbr1; index++){
        if (dr2[index] < cut2) {
            int pushf = atoms[iatom].nbrlist1[index];
            atoms[iatom].nbrlist2.push_back(pushf);
            eraselist[nerase] = index;
            nerase++;
        } else if (dr2[index] > pairdist2) {
            eraselist[nerase] = index;
            nerase++;
        }
    }
    for (int ii=0; ii<nerase; ii++){
        atoms[iatom].nbrlist1.erase(atoms[iatom].nbrlist1.begin()+eraselist[ii]-ii);
    }
    atoms[iatom].nbrlist2.shrink_to_fit(); atoms[iatom].nbrlist2.reserve(atoms[iatom].nbrlist2.size()+10);

    int nnbr2 = atoms[iatom].nbrlist2.size();
    index=0;
//#pragma simd
    for (int jatom : atoms[iatom].nbrlist2) {
        dij[index] = loci - pos[jatom];
        jtypes[index] = type[jatom];
        index++;
    }

#pragma simd
    for (index=0; index<nnbr2; index++) {
        if (dij[index].x>box_2.x) dij[index].x -= box.x;
        else if (dij[index].x<=-box_2.x) dij[index].x += box.x;
        if (dij[index].y>box_2.y) dij[index].y -= box.y;
        else if (dij[index].y<=-box_2.y) dij[index].y += box.y;
        if (dij[index].z>box_2.z) dij[index].z -= box.z;
        else if (dij[index].z<=-box_2.z) dij[index].z += box.z;
    }

#pragma simd
    for (index=0; index<nnbr2; index++){
        int jtype = jtypes[index];
        epsi_ij[index] = epsi[jtype];
        kqq_ij[index] = kqq[jtype];

        dr2[index] = dij[index].length2();
        r2_d2[index] = rmin2[jtype]/dr2[index];
    }

#pragma simd
    for (index=0; index<nnbr2; index++){
        dr_2[index] = 1.0/dr2[index];
        dr_1[index] = sqrt(dr_2[index]);

        dr_6[index] = pow(r2_d2[index],3);

        ljE[index] = epsi_ij[index]*dr_6[index]*(dr_6[index]-2.0);
        ljF[index] = 12.0*epsi_ij[index]*dr_6[index]*(dr_6[index]-1.0)*dr_2[index];

        elecE[index] = kqq_ij[index]*dr_1[index]*(1.0 - dr2[index]*cut_2);
        elecF[index] = elecE[index]*(dr_2[index] + 3.0*cut_2);
        elecE[index]*= (1.0 - dr2[index]*cut_2);
    }

    for (index=0; index<nnbr2; index++){
        if (dr2[index] < cut2) {
            FF[index] = (ljF[index] + elecF[index]) * dij[index];
        } else {
            ljE[index] = 0.0;
            FF[index] = 0.0;
        }
    }


    Evdw += accumulate(ljE,ljE+nnbr2,0.0);
    Eelec += accumulate(elecE,elecE+nnbr2,0.0);

    Vector FFt;
    for (index=0; index<nnbr2; index++){
        FFt += FF[index];
    }
    ff[iatom] += FFt;

    index = 0;
    for (int jatom : atoms[iatom].nbrlist2) {
        ff[jatom] -= FF[index];
        index++;
    }
}
}

void Nonbonded::Compute_dpd(const Initial *init, const Vector *pos, const Vector *vel, Vector *const ff, const int num, double &Evdw){

double gamma = init->gamma;
double sigma = init->sigma;
int ntype = init->ntype;
double rmin2[ntype], epsi[ntype], kqq[ntype];
Vector box, box_2;
box.x = init->box[0]; box.y=init->box[1]; box.z=init->box[2];
box_2 = box*0.5;


Vector loci, veli;
Vector dij[L2], norm_dij[L2], vij[L2];
Vector FF[L2] ;
double rand[L2];
double dr[L2], dr2[L2] ALIGN;
double fac[L2] ALIGN;
double omega[L2] ALIGN;
double dot[L2] ALIGN;
double dpdFc[L2], dpdFd[L2], dpdFr[L2] ALIGN;
int jtypes[L2]  ALIGN;
int eraselist[20];
double aij[init->ntype];

Evdw = 0.0;
for (int iatom=0; iatom<num; iatom++){
    for (int ii=0; ii<init->ntype; ii++) {
        aij[ii] = init->aij[init->dpd[iatom].type][ii];
    }

    loci = pos[iatom];
    veli = vel[iatom];

    int nnbr1 = atoms[iatom].nbrlist1.size();
    int index=0;
//#pragma simd
    for (int jatom : atoms[iatom].nbrlist1) {
        dij[index] = loci - pos[jatom];
        index++;
    }

#pragma simd
    for (index=0; index<nnbr1; index++){

        if (dij[index].x>box_2.x) dij[index].x -= box.x;
        else if (dij[index].x<=-box_2.x) dij[index].x += box.x;
        if (dij[index].y>box_2.y) dij[index].y -= box.y;
        else if (dij[index].y<=-box_2.y) dij[index].y += box.y;
        if (dij[index].z>box_2.z) dij[index].z -= box.z;
        else if (dij[index].z<=-box_2.z) dij[index].z += box.z;

        dr2[index] = dij[index].length2();
    }

    int nerase = 0;
    for (index=0; index<nnbr1; index++){
        if (dr2[index] < cut2) {
            int pushf = atoms[iatom].nbrlist1[index];
            atoms[iatom].nbrlist2.push_back(pushf);
            eraselist[nerase] = index;
            nerase++;
        } else if (dr2[index] > pairdist2) {
            eraselist[nerase] = index;
            nerase++;
        }
    }
    for (int ii=0; ii<nerase; ii++){
        atoms[iatom].nbrlist1.erase(atoms[iatom].nbrlist1.begin()+eraselist[ii]-ii);
    }
    atoms[iatom].nbrlist2.shrink_to_fit(); atoms[iatom].nbrlist2.reserve(atoms[iatom].nbrlist2.size()+10);

// Region 2         r < cut
    int nnbr2 = atoms[iatom].nbrlist2.size();
    index=0;
//#pragma simd
    for (int jatom : atoms[iatom].nbrlist2) {
        dij[index] = loci - pos[jatom];
        vij[index] = veli - vel[jatom];
        fac[index] = aij[init->dpd[jatom].type];
        index++;
    }


#pragma simd
    for (index=0; index<nnbr2; index++) {
        if (dij[index].x>box_2.x) dij[index].x -= box.x;
        else if (dij[index].x<=-box_2.x) dij[index].x += box.x;
        if (dij[index].y>box_2.y) dij[index].y -= box.y;
        else if (dij[index].y<=-box_2.y) dij[index].y += box.y;
        if (dij[index].z>box_2.z) dij[index].z -= box.z;
        else if (dij[index].z<=-box_2.z) dij[index].z += box.z;
    }

#pragma simd
    for (index=0; index<nnbr2; index++){
        dr[index] = dij[index].length();
        dot[index] = dij[index].dot(vij[index]);
        omega[index] = 1.0 - dr[index]*cut_1;
    }

//    vdRngUniform ( VSL_RNG_METHOD_UNIFORM_STD , stream, nnbr2, rand, -0.5, 0.5 );
    vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nnbr2, rand, 0.0, 1.0 );

    for (index=0; index<nnbr2; index++){
        dpdFc[index] = fac[index]*omega[index];
        dpdFd[index] = -gamma*omega[index]*omega[index]*dot[index]/dr[index];
        dpdFr[index] = sigma*omega[index]*rand[index];
        }
    
    for (index=0; index<nnbr2; index++){
        if (dr[index] < cut) {
            FF[index] = ((dpdFc[index] + dpdFd[index] + dpdFr[index])/dr[index]) * dij[index];
            Evdw += 0.5*dpdFc[index]*omega[index]*cut;
        } else {
            FF[index] = 0.0;
        }

    }

    Vector FFt;
    for (index=0; index<nnbr2; index++){
        FFt += FF[index];
    }
    ff[iatom] += FFt;

    index = 0;
    for (int jatom : atoms[iatom].nbrlist2) {
        ff[jatom] -= FF[index];
        index++;
    }

}
}

Nonbonded::Nonbonded(const Nonbonded& orig) {
}

Nonbonded::~Nonbonded() {
}
