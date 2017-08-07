/*
 * File:   Initial.cpp
 * Author: amin
 *
 * Created on September 10, 2015, 11:28 AM
 */

#include <iostream>
#include <string.h>
#include <vector>
#include <algorithm>
#include "Initial.h"
#include "Configure.h"
#include "Molecule.h"
#include "Parameters.h"

#define PI 3.14159265

using namespace std;

Initial::Initial(const Configure *conf, const Molecule *mol, const Parameters *params) {

    bonds = new BondElem[conf->numBonds];
    memset((void *)bonds, 0, conf->numBonds*sizeof(BondElem));
    build_bondlist(conf,mol,params);
    angles = new AngleElem[conf->numAngles];
    memset((void *)angles, 0, conf->numAngles*sizeof(AngleElem));
    build_anglelist(conf,mol,params);

    connections = new Connections[conf->numAtoms];
    memset((void *)connections, 0, conf->numAtoms*sizeof(Connections));

    if (strcmp(conf->mode,"lj") == 0) {
        vdw = new VdwElem[conf->numAtoms];
        memset((void *)vdw, 0, conf->numAtoms*sizeof(VdwElem));
        build_vdwlist(conf,mol,params);
    } else if (strcmp(conf->mode,"dpd") == 0) {
        dpd = new DpdElem[conf->numAtoms];
        memset((void *)dpd, 0, conf->numAtoms*sizeof(DpdElem));
        build_dpdlist(conf,mol,params);
    }


    pos = new Vector[conf->numAtoms];
    mass = new double[conf->numAtoms];
    for (int ii=0; ii<conf->numAtoms; ii++){
        pos[ii].x = mol->atoms[ii].coorx;
        pos[ii].y = mol->atoms[ii].coory;
        pos[ii].z = mol->atoms[ii].coorz;

        mass[ii] = mol->atoms[ii].mass;
    }

    residue = new int[conf->numAtoms];
    int res_ID = mol->atoms[0].res_ID;
    char res_name[5];
    strcpy(res_name,mol->atoms[0].resname);
    nresidue = 0;
    for (int ii=0; ii<conf->numAtoms; ii++){
        if (res_ID != mol->atoms[ii].res_ID || strcmp(res_name,mol->atoms[ii].resname) != 0 ) {
            res_ID = mol->atoms[ii].res_ID;
            strcpy(res_name,mol->atoms[ii].resname);
            ++nresidue;
        }
        residue[ii] = nresidue;
    }
    ++nresidue;

    for (int ii=0; ii<3; ii++) { box[ii] = conf->box[ii]; }

    // For DPD simulation
    gamma = conf->gamma;
    sigma = sqrt(2*gamma*conf->temp/conf->timestep);
}

void Initial::build_bondlist(const Configure *conf, const Molecule *mol, const Parameters *params){
    string str;
    char   seg[256];
    char  *pch;
    size_t found1, found2;

    for (int ii=0; ii<conf->numBonds; ii++){
        bonds[ii].atom1 = mol->bonds[ii].atom1;
        bonds[ii].atom2 = mol->bonds[ii].atom2;

        string atype1(mol->atoms[bonds[ii].atom1].atomtype); // Atom1 type
        string atype2(mol->atoms[bonds[ii].atom2].atomtype); // Atom2 type
        for (int jj=0; jj < params->bond_str.size(); jj++){

            str = params->bond_str[jj];
            found1=str.find(atype1); //find atom1 in the str
            if (found1!=std::string::npos) {
                str.replace(found1,atype1.length(),""); //remove the atom1 from str
                found2=str.find(atype2); //find atom1 in the str
                if (found2!=std::string::npos) {
                    strncpy(seg, str.c_str(), sizeof(seg));
                    seg[sizeof(seg) - 1] = 0;
                    pch = strtok (seg," ,:=");pch = strtok (NULL," ,:="); // Jumb two columns
                    bonds[ii].k = stof(strtok (NULL," ,:="));
                    bonds[ii].x0 = stof(strtok (NULL," ,:="));
                }
            }
        } // Loops all bonds parameters
        if (bonds[ii].k==0 && bonds[ii].x0==0) {
            cout << "ERROR : NO BOND PARAMETER FOR BOND BETWEEN ATOM " << bonds[ii].atom1+1 << " AND " << bonds[ii].atom2+1 << endl;
            exit(1);
        }

    } // Loops all bonds
}

void Initial::build_anglelist(const Configure *conf, const Molecule *mol, const Parameters *params){
    string str;
    char   seg[256];
    char  *pch1, *pch2, *pch3;

    for (int ii=0; ii<conf->numAngles; ii++){
        angles[ii].atom1 = mol->angles[ii].atom1;
        angles[ii].atom2 = mol->angles[ii].atom2;
        angles[ii].atom3 = mol->angles[ii].atom3;

        char *atype1 = mol->atoms[angles[ii].atom1].atomtype; // Atom1 type
        char *atype2 = mol->atoms[angles[ii].atom2].atomtype; // Atom2 type
        char *atype3 = mol->atoms[angles[ii].atom3].atomtype; // Atom3 type

        for (int jj=0; jj < params->angle_str.size(); jj++){
            str = params->angle_str[jj];
            strncpy(seg, str.c_str(), sizeof(seg));
            seg[sizeof(seg) - 1] = 0;
            pch1 = strtok (seg," ,:=");pch1 = strtok (NULL," ,:=");
            pch2 = strtok (NULL," ,:=");
            pch3 = strtok (NULL," ,:=");

            if (strcmp(atype2,pch2)==0){
                if ((strcmp(atype1,pch1)==0 && strcmp(atype3,pch3)==0)||(strcmp(atype1,pch3)==0 && strcmp(atype3,pch1)==0)){
                    angles[ii].k = stof(strtok (NULL," ,:="));
                    angles[ii].theta0 = stof(strtok (NULL," ,:=")); // * PI / 180.0 ;
                }
            }
        }
        if (angles[ii].k==0 && angles[ii].theta0==0) {
            cout << "ERROR : NO ANGLE PARAMETER FOR ANGLE BETWEEN ATOM " << angles[ii].atom1+1 << " AND "
                                                                         << angles[ii].atom2+1 << " AND "
                                                                         << angles[ii].atom3+1 << endl;
            exit(1);
        }
    }
}

void Initial::build_vdwlist(const Configure *conf, const Molecule *mol, const Parameters *params){
    string str;
    char   seg[256];
    char  *pch;
    size_t found;
    vector<string> stype;

    // Convert string type to int type
    string atype(mol->atoms[0].atomtype);
    stype.push_back(atype);
    for (int ii=1; ii<conf->numAtoms; ii++){
        bool flag = true;
        string atype(mol->atoms[ii].atomtype);
        for (int jj=0; jj<stype.size(); jj++) {
            if (atype==stype[jj]) {
                vdw[ii].type = jj;
                flag = false;
                break;
            }
        }
        if (flag) {
            vdw[ii].type = stype.size();
            stype.push_back(atype);

        }
    }
    ntype = stype.size();


    // Obtain rmin and epsilon
    for (int ii=0; ii<conf->numAtoms; ii++){
        string atype1(mol->atoms[ii].atomtype);
        for (int jj=0; jj<params->vdw_str.size(); jj++){
            str = params->vdw_str[jj];
            found=str.find(atype1);
            if (found!=std::string::npos) {
                strncpy(seg, str.c_str(), sizeof(seg));
                seg[sizeof(seg) - 1] = 0;
                pch = strtok (seg," ,:="); pch = strtok (NULL," ,:="); // Jump two columns
                vdw[ii].epsi = stof(strtok (NULL," ,:="));
                vdw[ii].rmin = stof(strtok (NULL," ,:="));
                vdw[ii].charge = mol->atoms[ii].charge;
            }
        }
    }

    // find excluded atoms for nonbonded calculation
    int a1, a2, a3;
    int con[conf->numAtoms][20]; // Save the list of all connections
    int ind[conf->numAtoms];     // Save the number of connections
    memset(ind, 0, sizeof ind);
    memset(con, 0, sizeof con);

    for (int ii=0; ii<conf->numAtoms; ii++){
        con[ii][ind[ii]] = ii; ind[ii]++; // exclude itself
    }

    for (int ii=0; ii<conf->numBonds; ii++){ // Loop over all bonds
        a1 = bonds[ii].atom1;
        a2 = bonds[ii].atom2;

        con[a1][ind[a1]] = a2; ind[a1]++; // Save the connection for atom1
        con[a2][ind[a2]] = a1; ind[a2]++;
    }
    for (int ii=0; ii<conf->numAngles; ii++){ //Loop over all angles
        a1 = angles[ii].atom1;
        a2 = angles[ii].atom2;
        a3 = angles[ii].atom3;
        // No need to save atom2 since the bond connections for atom2 was saved in bond section
        con[a1][ind[a1]] = a3; ind[a1]++; // Save the connection for atom1
        con[a3][ind[a3]] = a1; ind[a3]++; // Save the connection for atom1
    }

    for (int ii=0; ii<conf->numAtoms; ii++){ // Loop over all atoms to save connections
        vector<int> temp(con[ii],con[ii]+ind[ii]);

        sort(temp.begin(),temp.end()); // Sort connections
        connections[ii].exclist = new int[ind[ii]]; // Allocate memory to save all connections

        connections[ii].numexc = ind[ii];
        for (int jj=0; jj<ind[ii]; jj++){
            connections[ii].exclist[jj] = temp[jj];
        }
    }
}

void Initial::build_dpdlist(const Configure *conf, const Molecule *mol, const Parameters *params){
    string str;
    char   seg[256];
    char  *pch;
    size_t found1, found2;
    vector<string> stype;

    // Convert string type to int type
    string atype(mol->atoms[0].atomtype);
    stype.push_back(atype);
    for (int ii=1; ii<conf->numAtoms; ii++){
        bool flag = true;
        string atype(mol->atoms[ii].atomtype);
        for (int jj=0; jj<stype.size(); jj++) {
            if (atype==stype[jj]) {
                dpd[ii].type = jj;
                flag = false;
                break;
            }
        }
        if (flag) {
            dpd[ii].type = stype.size();
            stype.push_back(atype);
        }
    }
    ntype = stype.size();
    

    aij = new float*[ntype];
    for (int ii=0; ii<ntype; ii++) {
        aij[ii] = new float[ntype];
        memset((void *)aij[ii], 0, ntype*sizeof(ntype));
    }


    // Obtain aij
    for (int ii=0; ii<ntype; ii++){
        string atype1(stype[ii]);
        for (int line=0; line<params->dpd_str.size(); line++){
            str = params->dpd_str[line];
            found1=str.find(atype1);
            if (found1!=std::string::npos) {
                str.replace(found1,atype1.length(),""); //remove the atom1 from str
                for (int jj=0; jj<ntype; jj++) {
                    string atype2(stype[jj]);
                    found2=str.find(atype2);
                    if (found2!=std::string::npos) {
                        strncpy(seg, str.c_str(), sizeof(seg));
                        seg[sizeof(seg) - 1] = 0;
                        pch = strtok (seg," ,:="); pch = strtok (NULL," ,:="); // Jump two columns
                        aij[ii][jj] = stof(strtok (NULL," ,:="));
                    }
                }
            }
        }
    }

    for (int ii=0; ii<ntype; ii++){
        for (int jj=0; jj<ntype; jj++) {
            if (aij[ii][jj]==0) {
                cout << "There is no interaction parameter between atom type " << stype[ii] << " and " << stype[jj] << endl;
                exit(1);
            }
        }
    }

    // find excluded atoms for nonbonded calculation
    int a1, a2, a3;
    vector< vector<int> > con(conf->numAtoms,vector<int>(20));// Save the list of all connections
    int ind[conf->numAtoms];     // Save the number of connections
    memset(ind, 0, sizeof ind);

    for (int ii=0; ii<conf->numAtoms; ii++){
        con[ii][ind[ii]] = ii; ind[ii]++; // exclude itself
    }

    for (int ii=0; ii<conf->numBonds; ii++){ // Loop over all bonds
        a1 = bonds[ii].atom1;
        a2 = bonds[ii].atom2;
        
        con[a1][ind[a1]] = a2; ind[a1]++; // Save the connection for atom1
        con[a2][ind[a2]] = a1; ind[a2]++;
    }
    for (int ii=0; ii<conf->numAngles; ii++){ //Loop over all angles
        a1 = angles[ii].atom1;
        a2 = angles[ii].atom2;
        a3 = angles[ii].atom3;
        // No need to save atom2 since the bond connections for atom2 was saved in bond section
        con[a1][ind[a1]] = a3; ind[a1]++; // Save the connection for atom1
        con[a3][ind[a3]] = a1; ind[a3]++; // Save the connection for atom1
    }

    for (int ii=0; ii<conf->numAtoms; ii++){ // Loop over all atoms to save connections
        vector<int> temp(con[ii].begin(),con[ii].begin()+ind[ii]);

        sort(temp.begin(),temp.end()); // Sort connections
        connections[ii].exclist = new int[ind[ii]]; // Allocate memory to save all connections

        connections[ii].numexc = ind[ii];
        for (int jj=0; jj<ind[ii]; jj++){
            connections[ii].exclist[jj] = temp[jj];
        }
    }
}
bool Initial::exclude(const int iatom, const int jatom) const {
    int jj=0;
    bool flag = false;
    while (connections[iatom].exclist[jj]<=jatom && jj<connections[iatom].numexc){
        if (connections[iatom].exclist[jj++] == jatom){
            flag = true;
        }
    }
    return flag;
}

Initial::Initial(const Initial& orig) {
}

Initial::~Initial() {
}



