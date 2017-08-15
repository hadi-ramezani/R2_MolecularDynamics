#include <fstream>
#include "Integrator.h"
#include "Configure.h"
#include "Initial.h"
#include "Vector.h"
#include "Nonbonded.h"
#include "Bonded.h"

using namespace std;


Integrator::Integrator(const Configure *conf, const Initial *init, PDB *pdb, Parameters *params) : nonbonded(init, params, conf)
                                                                    , bonded(init, params)
                                                                    , dcd(conf->dcdname, conf->numAtoms, conf)
                                                                    , out(conf->enename, conf)
                                                                    , temp(conf){

    // define again atoms parameters this might help for parallelization
    pos = new Vector[conf->numAtoms];
    imass = new double[conf->numAtoms];
    for (int ii=0; ii<conf->numAtoms; ii++){
        imass[ii] = 1/init->atommass(ii);
    }
    
    pdb->get_all_positions(pos);

    vel = new Vector[conf->numAtoms];
    memset((void *)vel, 0, conf->numAtoms*sizeof(Vector));

    vel2 = new Vector[conf->numAtoms];
    memset((void *)vel2, 0, conf->numAtoms*sizeof(Vector));

    f = new Vector[conf->numAtoms];
    memset((void *)f, 0, conf->numAtoms*sizeof(Vector));

    dt=conf->timestep*dtfac;
    KB=0.001987191;
}

void Integrator::Loop(const Configure *conf, const Initial *init){

    bonded.compute(pos, f, Ebond, Eangle, Edihedral, Eimproper);
    nonbonded.compute(init, pos, f, Evdw, Eelec, conf);
    Etot = Ebond + Eangle + Edihedral + Eimproper + Evdw + Eelec + Ekin;
    cout << "Running the simulation..." << endl;

    for (int step=conf->fstep; step < conf->fstep + conf->nsteps + 1; step++){ // main loop

        if (step%conf->dcdFreq == 0) dcd.WriteFrame(conf->numAtoms,pos,conf->box);
        if (step%conf->energyFreq == 0) out.Print(step,step*conf->timestep,Ebond, Eangle, Edihedral, Eimproper, Evdw, Eelec, Ekin, Etot, temperature);

        for (int ii=0; ii<conf->numAtoms; ii++){ // Velocity Verlet
            pos[ii] += dt*vel[ii] + 0.5*dt*dt*f[ii]*imass[ii];
            vel[ii] += 0.5*dt*f[ii]*imass[ii];
        }

        //Compute force at time t+dt
        memset((void *)f, 0, conf->numAtoms*sizeof(Vector));
        bonded.compute(pos, f, Ebond, Eangle, Edihedral, Eimproper);
        nonbonded.compute(init, pos, f, Evdw, Eelec, conf);

        Ekin = 0; 
        for (int ii=0; ii<conf->numAtoms; ii++) {
            vel[ii] += 0.5*dt*f[ii]*imass[ii];
            Ekin += vel[ii]*vel[ii]*init->atommass(ii);
        }
        Ekin *= 0.5;

        Etot = Ebond + Eangle + Edihedral + Eimproper + Evdw + Eelec + Ekin;

        //if (strcmp(conf->thermostat,"on") == 0) temp.LA(conf->numAtoms,init,rmass,pos,&nonbonded,vel);
    }
}


Integrator::~Integrator() {
    delete [] f;
    delete [] vel;
    delete [] pos;
    delete [] imass;
}

