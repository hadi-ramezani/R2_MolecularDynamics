/*
 * File:   Thermostat.cpp
 * Author: amin
 *
 * Created on November 26, 2015, 10:32 AM
 */
#define PI 3.141592653
#define TIMEFACTOR 48.88821
#define L2 2048

#include <random>
#include <iostream>
#include <string.h>
#include "Vector.h"
#include "Thermostat.h"

using namespace std;

Thermostat::Thermostat(const Configure *conf) {
    if (strcmp(conf->mode,"dpd") == 0) {
        KB=1.0;
    } else {
        KB=0.001987191;
    }
    ndof = 3*conf->numAtoms;
    target = conf->temp;
    cout << "The target temperature is : " << target << endl;
    gamma = conf->gamma;
    kTs = sqrt(KB*target);
    dt_s = 1/sqrt(conf->timestep*TIMEFACTOR);
    cut_1 = 1/conf->cutoff;

    // Initializing random numbers
    vslNewStream( &stream, VSL_BRNG_SFMT19937, conf->seed );


}

void Thermostat::Initial_vel(const int num, const Initial *init, Vector *const vel){

    float rand[num];

    vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, num, rand, 0.0, sqrtf(2.0) );
    for ( int ii=0; ii<num; ii++ ){
        vel[ii].x = rand[ii];
    }

    vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, num, rand, 0.0, sqrtf(2.0) );
    for ( int ii=0; ii<num; ii++ ){
        vel[ii].y = rand[ii];
    }

    vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, num, rand, 0.0, sqrtf(2.0) );
    for ( int ii=0; ii<num; ii++ ){
        vel[ii].z = rand[ii];
    }

    double Ek = 0.0;
    for (int ii=0; ii< num; ii++){
        Ek += vel[ii]*vel[ii]*init->mass[ii];
    }
    double temperature = Ek/(ndof*KB);
    double factor = sqrt(target/temperature);

    for (int ii=0; ii< num; ii++){
        vel[ii] *= factor;
    }

}

void Thermostat::DPD(const int num, const double *dr_1, const double *dot, double *const dpdF){
    double omega[num], FD[num], FR[num];
    std::uniform_real_distribution<double> distribution(-0.5,0.5);

    for (int index=0; index<num; index++){
        omega[index] = 1.0 - cut_1/dr_1[index];
        FD[index] = gamma*omega[index]*omega[index]*dot[index];
        FR[index] = sigma*omega[index]*dt_s*distribution(generator);
    }
    for (int index=0; index<num; index++){
        dpdF[index] = FD[index] + FR[index];
    }
}

void Thermostat::LA(const int num, const Initial *init, const double *rmass, const Vector *pos, const Nonbonded *nonbonded, Vector * const vel){

    Vector loci, veli;
    Vector dij[L2], vij[L2];
    Vector vel_sum[num];
    int pairlist[L2];
    memset((void *)vel_sum, 0, num*sizeof(Vector));

    double rmassi, rmassj[L2];
    double dot[L2], miue[L2], deltav[L2];

    Vector box, box_2;
    box.x = init->box[0]; box.y=init->box[1]; box.z=init->box[2];
    box_2 = box*0.5;

    std::uniform_real_distribution<double> distribution_uniform(0,1.0);
    std::normal_distribution<double> distribution_gauss(0.0,1.0);

    for (int iatom=0; iatom<num; iatom++){
            loci = pos[iatom];
            veli = vel[iatom];
            rmassi = rmass[iatom];
            int nnbr = nonbonded->atoms[iatom].nbrlist2.size();

            float rand_uni[nnbr];
            vsRngUniform ( VSL_RNG_METHOD_UNIFORM_STD , stream, nnbr, rand_uni, .0, 1.0 );

            int npair = 0;
            int index = 0;
            for (int jatom : nonbonded->atoms[iatom].nbrlist2) {
                if ( rand_uni[index] < gamma){
                    pairlist[npair] = jatom;
                    dij[npair] = loci - pos[jatom];
                    vij[npair] = veli - vel[jatom];
                    rmassj[npair] = rmass[jatom];
                    npair++;
                }
                index++;
            }


            for (int index=0; index<npair; index++) {
                if (dij[index].x>box_2.x) dij[index].x -= box.x;
                else if (dij[index].x<=-box_2.x) dij[index].x += box.x;
                if (dij[index].y>box_2.y) dij[index].y -= box.y;
                else if (dij[index].y<=-box_2.y) dij[index].y += box.y;
                if (dij[index].z>box_2.z) dij[index].z -= box.z;
                else if (dij[index].z<=-box_2.z) dij[index].z += box.z;
            }

            for (int index=0; index<npair; index++) {
                dij[index] /= dij[index].length();
                dot[index] = dij[index].dot(vij[index]);
            }

            double rand_gauss[npair];
            vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, npair, rand_gauss, 0.0, 1.0 );

            for (int index=0; index<npair; index++) {
                miue[index] = 1.0 / (rmassi + rmassj[index]);
                deltav[index] = rand_gauss[index]*kTs*sqrt(miue[index]);
                deltav[index] -= miue[index]*dot[index];
            }

            for (int index=0; index<npair; index++) {
                int jatom = pairlist[index];
                vel_sum[iatom] += deltav[index]*dij[index]*rmassi;
                vel_sum[jatom] -= deltav[index]*dij[index]*rmassj[index];
            }

    }
    for (int index=0; index<num; index++) {
        vel[index] += vel_sum[index];
    }
}

Thermostat::Thermostat(const Thermostat& orig) {
}

Thermostat::~Thermostat() {
}

