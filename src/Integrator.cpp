/*
 * File:   Integrator.cpp
 * Author: amin
 *
 * Created on September 11, 2015, 5:17 PM
 */

#include <fstream>
#include "Integrator.h"
#include "Configure.h"
#include "Initial.h"
#include "Vector.h"
#include "Nonbonded.h"
#include "Bonded.h"

using namespace std;


Integrator::Integrator(const Configure *conf, const Initial *init)  : nonbonded(conf->numAtoms, conf->cutoff, conf->switchdist, conf->pairlistdist, conf->seed)
                                                                    , dcd(conf->dcdname, conf->numAtoms)
                                                                    , out(conf->enename)
                                                                    , temp(conf){

    // define again atoms parameters this might help me for parallelization
    pos = new Vector[conf->numAtoms];
    rmass = new double[conf->numAtoms];
    for (int ii=0; ii<conf->numAtoms; ii++){
        pos[ii] = init->pos[ii];
        rmass[ii] = 1/init->mass[ii];
    }

    vel = new Vector[conf->numAtoms];
    memset((void *)vel, 0, conf->numAtoms*sizeof(Vector));

    vel2 = new Vector[conf->numAtoms];
    memset((void *)vel2, 0, conf->numAtoms*sizeof(Vector));

    ff = new Vector[conf->numAtoms];
    memset((void *)ff, 0, conf->numAtoms*sizeof(Vector));

    if (strcmp(conf->mode,"dpd") == 0) {
        dt=conf->timestep;
        lambdaDt = dt * conf->lambda;
        KB=1.0;
    } else {
        dt=conf->timestep*dtfac;
        KB=0.001987191;
    }
}

void Integrator::Loop(const Configure *conf, const Initial *init){

    nonbonded.Build_cells(init->box, conf->celldist, conf->numAtoms); // I build cells only once; this might cause problem if the simulation box change too much on NPT ensemble
    nonbonded.built_table(init, init->ntype, conf->numAtoms);
    nonbonded.Neighborlist(init->box, conf->numAtoms, init, pos);

    nonbonded.Compute(init, pos, ff, conf->numAtoms, Evdw, Eelec); // calculate force at time 0

    if (strcmp(conf->rigidBonds,"Yes") !=0 ) {
        bonded.Compute_bond(init,pos,ff,conf->numBonds,Ebond);
    } else {
        cout << "RIGID BONDS" << endl;
    }

    bonded.Compute_angle(init,pos,ff,conf->numAngles, Eangle);

    if (strcmp(conf->settemp,"Yes") == 0 ) temp.Initial_vel(conf->numAtoms,init,vel);

    Ekin =0.0;
    for (int ii=0; ii<conf->numAtoms; ii++){
        Ekin += vel[ii]*vel[ii]*init->mass[ii];
    }
    temperature = Ekin/(temp.ndof*KB); Ekin *= 0.5;

    cout << "Running the simulation..." << endl;

    for (int step=conf->fstep; step < conf->fstep + conf->nsteps + 1; step++){ // main loop

        Etot = Ebond + Eangle + Evdw + Eelec + Ekin;
        if (step%conf->dcdFreq == 0) dcd.frames(conf->numAtoms,pos,init->box);
        if (step%conf->energyFreq == 0) out.Print(step,step*conf->timestep,Ebond, Eangle, Evdw, Eelec, Ekin, Etot, temperature);

        for (int ii=0; ii<conf->numAtoms; ii++){ // Velocity Verlet
            pos[ii] += dt*vel[ii] + 0.5*dt*dt*ff[ii]*rmass[ii];
            vel2[ii] = vel[ii] + 0.5*dt*ff[ii]*rmass[ii];
        }

        //Compute force at time t+dt
        memset((void *)ff, 0, conf->numAtoms*sizeof(Vector));
        if (step%conf->pairlistFreq == 0) nonbonded.Neighborlist(init->box, conf->numAtoms, init, pos);
        if (step%conf->nonbondedFreq == 0) nonbonded.Compute(init, pos, ff, conf->numAtoms, Evdw, Eelec);

        if (strcmp(conf->rigidBonds,"Yes") !=0 ) {
            bonded.Compute_bond(init,pos,ff,conf->numBonds,Ebond);
        } else {
            cout << "RIGID BONDS" << endl;
        }


        bonded.Compute_angle(init,pos,ff,conf->numAngles, Eangle);

        for (int ii=0; ii<conf->numAtoms; ii++){ // Velocity Verlet
            vel[ii] = vel2[ii] + 0.5*dt*ff[ii]*rmass[ii];
        }

        if (strcmp(conf->thermostat,"on") == 0) temp.LA(conf->numAtoms,init,rmass,pos,&nonbonded,vel);

        Ekin = 0;
        for (int ii=0; ii<conf->numAtoms; ii++){ // Velocity Verlet
            Ekin += vel[ii]*vel[ii]*init->mass[ii];
        }

        temperature = Ekin/(temp.ndof*KB); Ekin *= 0.5;
    }
}

void Integrator::Loop_dpd(const Configure *conf, const Initial *init){


    nonbonded.Build_cells(init->box, conf->celldist, conf->numAtoms); // I build cells only once; this might cause problem if the simulation box change too much on NPT ensemble
    nonbonded.Neighborlist(init->box, conf->numAtoms, init, pos);

    nonbonded.Compute_dpd(init, pos, vel2, ff, conf->numAtoms, Evdw ); // calculate force at time 0

    if (strcmp(conf->rigidBonds,"Yes") !=0 ) {
        bonded.Compute_bond(init,pos,ff,conf->numBonds,Ebond);
    } else {
        cout << "RIGID BONDS" << endl;
    }
    bonded.Compute_angle_ub(init,pos,ff,conf->numAngles, Eangle);
    if (strcmp(conf->settemp,"Yes") == 0 ) temp.Initial_vel(conf->numAtoms,init,vel);
    Ekin =0.0;
    for (int ii=0; ii<conf->numAtoms; ii++){
        Ekin += vel[ii]*vel[ii]*init->mass[ii];
    }
    temperature = Ekin/(temp.ndof*KB); Ekin *= 0.5;

    cout << "Running the simulation..." << endl;
    for (int step=conf->fstep; step < conf->fstep + conf->nsteps + 1; step++){ // main loop
        Evdw /= conf->numAtoms;
        Etot = Ebond + Eangle + Evdw + Ekin;
        if (step%conf->wrapFreq ==0 || step%conf->dcdFreq == 0) {
            cout << "Step " << step << " wrap all coordinates around periodic boundaries" << endl;
            out.wrap(init, init->box, pos);
        }
        if (step%conf->dcdFreq == 0) {
            cout << "Step " << step << " write dcd file" << endl;
            dcd.frames(conf->numAtoms,pos,init->box);
        }
        if (step%conf->energyFreq == 0) {
            cout << "Step " << step << " write output file" << endl;
            out.Print(step,step*conf->timestep,Ebond, Eangle, Evdw, Eelec, Ekin, Etot, temperature);
        }

        for (int ii=0; ii<conf->numAtoms; ii++){ // Modified Velocity Verlet
            vel2[ii] = vel[ii] + lambdaDt*ff[ii]*rmass[ii];
            vel[ii]  = vel[ii] + 0.5*dt*ff[ii]*rmass[ii];
            pos[ii]  = pos[ii] + dt*vel[ii];
        }

        //Compute force at time t+dt
        memset((void *)ff, 0, conf->numAtoms*sizeof(Vector));
        if (step%conf->pairlistFreq == 0) nonbonded.Neighborlist(init->box, conf->numAtoms, init, pos);
        if (step%conf->nonbondedFreq == 0) nonbonded.Compute_dpd(init, pos, vel2, ff, conf->numAtoms, Evdw);

        if (strcmp(conf->rigidBonds,"Yes") !=0 ) {
            bonded.Compute_bond(init,pos,ff,conf->numBonds,Ebond);
        } else {
//            cout << "RIGID BONDS" << endl;
        }
        bonded.Compute_angle_ub(init,pos,ff,conf->numAngles, Eangle);


//        ext.Boundary(conf->numAtoms, pos, ff);

        for (int ii=0; ii<conf->numAtoms; ii++){ // Velocity Verlet
            vel[ii] = vel[ii] + 0.5*dt*ff[ii]*rmass[ii];
        }


//        ext.Boundary(conf->numAtoms, pos, vel);
        
        Ekin = 0;
        for (int ii=0; ii<conf->numAtoms; ii++){ // Kinetic energy calculation
            Ekin += vel[ii]*vel[ii]*init->mass[ii];
        }
        temperature = Ekin/(temp.ndof); Ekin *= 0.5;


    }
}
//Integrator::Integrator(const Integrator& orig) : nonbonded() {

//}

Integrator::~Integrator() {
}

