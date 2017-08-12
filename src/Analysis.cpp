/*
 * Author: Hadi
 *
 */

#include <fstream>
#include "Integrator.h"
#include "Configure.h"
#include "Initial.h"
#include "Vector.h"
#include "Nonbonded.h"
#include "Bonded.h"
#include "Analysis.h"

using namespace std;

Analysis::Analysis(const Configure *conf, const Initial *init) : nonbonded(conf->numAtoms, conf->cutoff, conf->switchdist, conf->pairlistdist, conf->seed)
                                                                    , dcd(conf->aDcdName, conf->numAtoms, conf)
                                                                    , out(conf->enename, conf) {
}

void Analysis::run(const Configure *conf, const Initial* init) {
    
    // define again atoms parameters this might help for parallelization
    pos = new Vector[conf->numAtoms];
    rmass = new double[conf->numAtoms];
    for (int ii=0; ii<conf->numAtoms; ii++){
        pos[ii] = init->pos[ii];
        rmass[ii] = 1/init->mass[ii];
    }

    ff = new Vector[conf->numAtoms];
    memset((void *)ff, 0, conf->numAtoms*sizeof(Vector));

    unsigned int frameNum = 0;

    nonbonded.built_table(init, init->ntype, conf->numAtoms);

    // flag to check if we have reached the end of dcd file
    bool flag = true;
    
    while(flag){

        // Read a frame
        flag = dcd.ReadFrame(conf->numAtoms, init->pos, aBox);
        cout << "Reading frame " << frameNum << endl;

        //Compute and print energies

        bonded.Compute_angle(init,pos,ff,conf->numAngles, Eangle);
        bonded.Compute_bond(init, pos, ff, conf->numBonds, Ebond);
        nonbonded.Build_cells(aBox, conf->celldist, conf->numAtoms);
        nonbonded.Neighborlist(aBox, conf->numAtoms, init, pos);
        nonbonded.Compute(init, pos, ff, conf->numAtoms, Evdw, Eelec);
        Etot = Ebond + Eangle + Evdw + Eelec;
        out.Print(frameNum, Ebond, Eangle, Evdw, Eelec, Etot);

        frameNum++;
    }
}

Analysis::~Analysis() {
}