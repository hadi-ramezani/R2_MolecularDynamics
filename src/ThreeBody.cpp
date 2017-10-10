#include <fstream>
#include "Integrator.h"
#include "Configure.h"
#include "Initial.h"
#include "Vector.h"
#include "Nonbonded.h"
#include "Bonded.h"
#include "Analysis.h"
#include "ThreeBody.h"

using namespace std;

ThreeBody::ThreeBody(const Configure *conf, const Initial *init, const PDB *pdb, Parameters *params): Analysis(conf, init, pdb, params) {

}

void ThreeBody::run(const Configure *conf, const Initial* init, const PDB *pdb) {
    

    // define again atoms parameters this might help for parallelization
    pos = new Vector[conf->numAtoms];

    // needed to avoid redifinition of compute functions
    f = new Vector[conf->numAtoms];
    memset((void *)f, 0, conf->numAtoms*sizeof(Vector));

    unsigned int frameNum = 0;

    // flag to check if we have reached the end of dcd file
    bool flag = true;
    
    while(flag){

        // Read a frame
        flag = dcd.read_frame(conf->numAtoms, pos, aBox);
        cout << "Reading frame " << frameNum << endl;

        //Compute and print energies
        //bonded.compute(pos, f, Ebond, Eangle, Edihedral, Eimproper);
        nonbonded.build_mycells(init, pos, f, conf);
        //nonbonded.build_neighborlist(init, pos, f, conf);
        //nonbonded.mycompute(init, pos, f, Evdw, Eelec, conf);

        nonbonded.threebody_neighborlist(init, pos, f, conf);
        nonbonded.compute_threebody(init, pos, f, Emisc, conf);
        Etot = Ebond + Eangle + Edihedral + Eimproper + Evdw + Eelec + Emisc;
        out.print(frameNum, Ebond, Eangle, Edihedral, Eimproper, Evdw, Eelec, Emisc, Etot);
        frameNum++;
    }
}

ThreeBody::~ThreeBody() {
}