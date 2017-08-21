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

Analysis::~Analysis() {
}