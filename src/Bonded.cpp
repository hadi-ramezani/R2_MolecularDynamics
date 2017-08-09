/*
 * File:   Bonded.cpp
 * Author: amin
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
        if (cos_theta >= 1.0) cos_theta = .99999;     // use to avoid numerical calculation error
        else if (cos_theta <= -1.0) cos_theta = -0.99999;
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        double theta = acos(cos_theta);
        double diff = theta - angle->theta0;
        //cout << diff << endl;
        Eangle += angle->k * diff * diff;


        //Force  (HADI please check the equation :D)
        Vector r12n = r12/d12;
        Vector r32n = r32/d32;
        double c12 = -2.0 * angle->k * diff / (sin_theta * d12);
        double c32 = -2.0 * angle->k * diff / (sin_theta * d32);
        Vector f12 = c12*(r12n*cos_theta - r32n);
        Vector f32 = c32*(r32n*cos_theta - r12n);

        ff[angle->atom1] += f12;
        ff[angle->atom2] -= f12+f32;
        ff[angle->atom3] += f32;
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

