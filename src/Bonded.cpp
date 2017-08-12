/*
 * File:   Bonded.cpp
 * Author: Hadi and Amin
 *
 * Created on September 23, 2015, 9:46 AM
 */

#include "Bonded.h"
#include "Initial.h"


Bonded::Bonded() {
}

void Bonded::Compute_bond(const Initial *init, const Vector *pos,Vector *const ff, const int num, double &Ebond){
    Ebond = 0.0;
    for (int ii=0; ii<num; ii++){
        BondElem *bond = init->bonds + ii;
        Vector r12 = pos[bond->atom1] - pos[bond->atom2];
        double dist = r12.length();
        double diff = dist - bond->x0;
        Vector f12 = -2 * bond->k * diff * (r12/dist);
        Ebond += bond->k * diff * diff;
        ff[bond->atom1] += f12;
        ff[bond->atom2] -= f12;

        bond++;
    }
}

void Bonded::Compute_angle(const Initial *init, const Vector *pos,Vector *const ff, const int num, double &Eangle){
    Eangle = 0;
    for (int ii=0; ii<num; ii++){
        AngleElem *angle = init->angles + ii;
        const Vector *coor1 = pos + angle->atom1;
        const Vector *coor2 = pos + angle->atom2;
        const Vector *coor3 = pos + angle->atom3;

        Vector r12 = *coor1 - *coor2;
        Vector r32 = *coor3 - *coor2;
        double d12 = r12.length();
        double d32 = r32.length();
        double cos_theta = (r12*r32)/(d12*d32);
        if (cos_theta >= 1.0) cos_theta = 1.0;     // use to avoid numerical calculation error
        else if (cos_theta <= -1.0) cos_theta = -1.0;
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        double theta = acos(cos_theta);
        double diff = theta - angle->theta0;
        Eangle += angle->k * diff * diff;

        //Forces
        double d12inv = 1.0/d12;
        double d32inv = 1.0/d32;
        diff *= (-2.0 * angle->k) / sin_theta;
        double c1 = diff * d12inv;
        double c2 = diff * d32inv;
        Vector f12 = c1*(r12*(d12inv*cos_theta) - r32*d32inv);
        Vector f1 = f12;
        Vector f32 = c2*(r32*(d32inv*cos_theta) - r12*d12inv);
        Vector f3 = f32;
        Vector f2 = -f12 - f32;

        if (angle->k_ub > 0.0) {

            Vector r13 = r12 - r32;
            double d13 = r13.length();
            diff = d13 - angle->r_ub;
            Eangle += angle->k_ub * diff * diff;

            // ub forces
            diff *= -2.0*angle->k_ub / d13;
            r13 *= diff;
            f1 += r13;
            f3 -= r13;
    } 

        ff[angle->atom1] += f1;
        ff[angle->atom2] += f2;
        ff[angle->atom3] += f3;
    }
}
void Bonded::Compute_angle_ub(const Initial *init, const Vector *pos,Vector *const ff, const int num, double &Eangle){
    Eangle = 0;
    for (int ii=0; ii<num; ii++){
        AngleElem *angle = init->angles + ii;
        const Vector *coor1 = pos + angle->atom1;
        const Vector *coor3 = pos + angle->atom3;

        Vector r13 = *coor1 - *coor3;
        double d13 = r13.length();

        double diff = d13-angle->theta0;
        Eangle += angle->k * diff * diff;

        diff *= -2.0 * angle->k / d13;
        ff[angle->atom1] += diff * r13;
        ff[angle->atom3] -= diff * r13;
    }
}
Bonded::Bonded(const Bonded& orig) {
}

Bonded::~Bonded() {
}

