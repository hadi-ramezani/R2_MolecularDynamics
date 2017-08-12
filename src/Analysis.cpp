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
                                                                    , out(conf->enename) {
}

void Analysis::run(const Configure *conf, const Initial *init) {
    
    dcd.ReadFrame(conf->numAtoms, init->pos);
}

Analysis::~Analysis() {
}