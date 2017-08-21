#ifndef THREEBODY_H
#define	THREEBODY_H

#include "Configure.h"
#include "Initial.h"
#include "Vector.h"
#include "Nonbonded.h"
#include "Bonded.h"
#include "Trajectory.h"
#include "Output.h"
#include "R2.h"

using namespace std;

class ThreeBody : public Analysis {

public:

    ThreeBody(const Configure *conf, const Initial *init, const PDB *pdb, Parameters *params);
    void run(const Configure *conf, const Initial *init, const PDB *pdb);
    virtual ~ThreeBody();
};

#endif	/* THREEBODY_H */
