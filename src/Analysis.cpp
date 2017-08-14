#include <fstream>
#include "Integrator.h"
#include "Configure.h"
#include "Initial.h"
#include "Vector.h"
#include "Nonbonded.h"
#include "Bonded.h"
#include "Analysis.h"

using namespace std;

Analysis::Analysis(const Configure *conf, const Initial *init, const PDB *pdb, Parameters *params) : nonbonded(init, params, conf)
                                                                    , bonded(init, params)
                                                                    , dcd(conf->aDcdName, conf->numAtoms, conf)
                                                                    , out(conf->enename, conf) {
}

void Analysis::run(const Configure *conf, const Initial* init, const PDB *pdb) {
    

    // define again atoms parameters this might help for parallelization
    pos = new Vector[conf->numAtoms];

    //pdb->get_all_positions(pos);

    // needed to avoid redifinition of compute functions
    f = new Vector[conf->numAtoms];
    memset((void *)f, 0, conf->numAtoms*sizeof(Vector));

    unsigned int frameNum = 0;

    // flag to check if we have reached the end of dcd file
    bool flag = true;
    
    while(flag){

        // Read a frame
        flag = dcd.ReadFrame(conf->numAtoms, pos, aBox);
        cout << "Reading frame " << frameNum << endl;

        //Compute and print energies
        bonded.compute(pos, f, Ebond, Eangle, Edihedral, Eimproper);
        nonbonded.compute(init, pos, f, Evdw, Eelec);
        Etot = Ebond + Eangle + Edihedral + Eimproper + Evdw + Eelec;
        out.Print(frameNum, Ebond, Eangle, Edihedral, Eimproper, Evdw, Eelec, Etot);
        frameNum++;
    }
}

Analysis::~Analysis() {
}