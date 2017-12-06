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
    
    //while(flag){
    for (int ii=0; ii < dcd.nsets; ii++) {
        if (flag) {
            // Read a frame
            if (dcd.setsread == 0){
                cout << "Reading frame " << conf->dcdFirst << endl;
                frameNum = conf->dcdFirst;
            } else if (dcd.setsread <= dcd.nsets && dcd.setsread != 0) {
                cout << "Reading frame " << dcd.setsread << endl;            
                frameNum = dcd.setsread;
            }
            flag = dcd.read_frame(conf->numAtoms, pos, aBox, conf);

            //Compute and print energies
            bonded.compute(pos, f, Ebond, Eangle, Edihedral, Eimproper);

            nonbonded.analysis_cells(init, pos, f, conf, aBox);
            nonbonded.analysis_nbrlist(init, pos, f, aBox);
            nonbonded.analysis_comp(init, pos, f, Evdw, Eelec, aBox);

            nonbonded.analysis_threebody_nbrlist(init, pos, f, aBox);
            nonbonded.analysis_comp_threebody(init, pos, f, Emisc, aBox);
            Etot = Ebond + Eangle + Edihedral + Eimproper + Evdw + Eelec + Emisc;
            if (flag) {
                out.print(frameNum, Ebond, Eangle, Edihedral, Eimproper, Evdw, Eelec, Emisc, Etot, aBox[0], aBox[1], aBox[2]);
            }
        }
    }
}

ThreeBody::~ThreeBody() {
}