#include "Bonded.h"
#include "Initial.h"
#include "Parameters.h"
#include "Vector.h"
#include <math.h>
#include <iostream>

struct BondElem {
  int atom1, atom2;
  double x0;
  double k;
};

struct AngleElem {
  int atom1, atom2, atom3;
  double k, theta0;
  double k_ub, r_ub;
};

struct DihedralElem {
  int atom1, atom2, atom3, atom4;
  int multiplicity;        // Max of 4; 2 would probably suffice 
  double k[4], delta[4];
  int n[4];
};

struct ImproperElem {
  int atom1, atom2, atom3, atom4;
  int multiplicity;        // Max of 4; 2 would probably suffice 
  double k[4], delta[4];
  int n[4];
};

Bonded::Bonded(const Initial *init, Parameters *params) {

    nbonds = init->numBonds;
    nangles = init->numAngles;
    ndihedrals = init->numDihedrals;
    nimpropers = init->numImpropers;

    bonds = new BondElem[nbonds];
    angles = new AngleElem[nangles];
    dihedrals = new DihedralElem[ndihedrals];
    impropers = new ImproperElem[nimpropers];

    build_bondlist(init, params);
    build_anglelist(init, params);
    build_dihedrallist(init, params);
    build_improperlist(init, params);
}

Bonded::~Bonded() {
      
    delete [] bonds;
    delete [] angles;
    delete [] dihedrals;
    delete [] impropers;
}

void Bonded::build_bondlist(const Initial *init, 
                                   const Parameters *params) {
    BondElem *bond = bonds;
    for (int i=0; i<nbonds; i++) {
        bond->atom1 = init->get_bond(i)->atom1;
        bond->atom2 = init->get_bond(i)->atom2;
        params->get_bond_params( &(bond->k), &(bond->x0), init->get_bond(i)->bond_type);
        bond++;
    }
}

void Bonded::build_anglelist(const Initial *init,
                       const Parameters *params) {
    AngleElem *angle = angles;
    for (int i=0; i<nangles; i++) {
        angle->atom1 = init->get_angle(i)->atom1;
        angle->atom2 = init->get_angle(i)->atom2;
        angle->atom3 = init->get_angle(i)->atom3;
        params->get_angle_params( &(angle->k), &(angle->theta0), 
                                  &(angle->k_ub), &(angle->r_ub), 
                                init->get_angle(i)->angle_type);
        angle++;
    }
}

void Bonded::build_dihedrallist(const Initial *init,
                                       const Parameters *params) {
    DihedralElem *dihedral = dihedrals;
    for (int i=0; i<ndihedrals; i++) {
        const Dihedral *dmol = init->get_dihedral(i);
        Index type = dmol->dihedral_type;
        dihedral->atom1 = dmol->atom1;
        dihedral->atom2 = dmol->atom2;
        dihedral->atom3 = dmol->atom3;
        dihedral->atom4 = dmol->atom4;
        dihedral->multiplicity = params->get_dihedral_multiplicity(type); 
    for (int j=0; j<dihedral->multiplicity; j++) 
        params->get_dihedral_params(&(dihedral->k[j]),
                                      &(dihedral->n[j]),
                                    &(dihedral->delta[j]),
                                    type, j);
        dihedral++;
    }
}

void Bonded::build_improperlist(const Initial *init,
                                       const Parameters *params) {
  
    ImproperElem *improper = impropers;
    for (int i=0; i<nimpropers; i++) {
        const Improper *imol = init->get_improper(i);
        Index type = imol->improper_type;
        improper->atom1 = imol->atom1;
        improper->atom2 = imol->atom2;
        improper->atom3 = imol->atom3;
        improper->atom4 = imol->atom4;
        improper->multiplicity = params->get_improper_multiplicity(type); 
    for (int j=0; j<improper->multiplicity; j++) 
          params->get_improper_params(&(improper->k[j]),
                                      &(improper->n[j]),
                                      &(improper->delta[j]),
                                      type, j);
        improper++;
    }
}
void Bonded::compute(const Vector *coords, Vector *f,
                            double& Ebond, double& Eangle, 
                            double &Edihedral, double &Eimproper) const {
  
    Ebond = compute_bonds(coords, f);
    Eangle = compute_angles(coords, f);
    Edihedral = compute_dihedrals(coords, f);
    Eimproper= compute_impropers(coords, f);
}

double Bonded::compute_bonds(const Vector *coords, Vector *f) const {
    
    double energy = 0.0;
    BondElem *bond = bonds;
    for (int i=0; i<nbonds; i++) {
        Vector r12 = coords[bond->atom1] - coords[bond->atom2];
        double r = r12.length();
        double diff = r - bond->x0;
        Vector f12 = r12 * diff * (-2*bond->k) / r;
        energy += bond->k * diff * diff; 
        f[bond->atom1] += f12;
        f[bond->atom2] -= f12; 
        bond++;
    }
    return energy;
}

double Bonded::compute_angles(const Vector *coords, Vector *f) const {
  
    double energy = 0.0;
    AngleElem *angle = angles;
    for (int i=0; i<nangles; i++) {
        Vector f1, f2, f3;
        const Vector *pos1 = coords + angle->atom1;
        const Vector *pos2 = coords + angle->atom2;
        const Vector *pos3 = coords + angle->atom3;
        Vector r12 = *pos1 - *pos2;
        Vector r32 = *pos3 - *pos2;
        double d12 = r12.length();
        double d32 = r32.length();
        double cos_theta = (r12*r32)/(d12*d32);
        if (cos_theta > 1.0) cos_theta = 1.0;
        else if (cos_theta < -1.0) cos_theta = -1.0;
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        double theta = acos(cos_theta);
        double diff = theta-angle->theta0;
        energy += angle->k * diff * diff; 

        // forces
        double d12inv = 1.0/d12;
        double d32inv = 1.0/d32;
        diff *= (-2.0*angle->k) / sin_theta;
        double c1 = diff * d12inv;
        double c2 = diff * d32inv;
        Vector f12 = c1*(r12*(d12inv*cos_theta) - r32*d32inv);
        f1 = f12;
        Vector f32 = c2*(r32*(d32inv*cos_theta) - r12*d12inv);
        f3 = f32;
        f2 = -f12 - f32;

        if (angle->k_ub > 0.0) {
            Vector r13 = r12 - r32;
            double d13 = r13.length();
            diff = d13 - angle->r_ub;
            energy += angle->k_ub * diff * diff;

            // ub forces
            diff *= -2.0*angle->k_ub / d13;
            r13 *= diff;
            f1 += r13;
            f3 -= r13;
        } 
        f[angle->atom1] += f1;
        f[angle->atom2] += f2;
        f[angle->atom3] += f3;

        angle++;
    }
    return energy;
}
 
double Bonded::compute_dihedrals(const Vector *coords, Vector *f) const {
  
    double energy = 0.0;
    DihedralElem *dihedral = dihedrals;
    for (int i=0; i<ndihedrals; i++) {
        const Vector *pos0 = coords + dihedral->atom1;
        const Vector *pos1 = coords + dihedral->atom2;
        const Vector *pos2 = coords + dihedral->atom3;
        const Vector *pos3 = coords + dihedral->atom4;
        const Vector r12 = *pos0 - *pos1;
        const Vector r23 = *pos1 - *pos2;
        const Vector r34 = *pos2 - *pos3;

        Vector dcosdA;
        Vector dcosdB;
        Vector dsindC;
        Vector dsindB;
        Vector f1, f2, f3;

        Vector A, B, C;
        A.cross(r12, r23);
        B.cross(r23, r34);
        C.cross(r23, A);

        double rA = A.length();
        double rB = B.length(); 
        double rC = C.length();

        double cos_phi = (A*B)/(rA*rB);
        double sin_phi = (C*B)/(rC*rB);

        // Normalize B
        rB = 1.0/rB;
        B *= rB;

        double phi = -atan2(sin_phi, cos_phi);

        if (fabs(sin_phi) > 0.1) {
            // Normalize A
            rA = 1.0/rA;
            A *= rA;
            dcosdA = rA*(cos_phi*A-B);
            dcosdB = rB*(cos_phi*B-A);
        }
        else {
            // Normalize C
            rC = 1.0/rC;
            C *= rC;
            dsindC = rC*(sin_phi*C-B);
            dsindB = rB*(sin_phi*B-C);
        }
 
        int mult = dihedral->multiplicity; 
        for (int j=0; j<mult; j++) {
            double k = dihedral->k[j];
            double n = dihedral->n[j];
            double delta = dihedral->delta[j];
            double K, K1;
            if (n) {
                K = k * (1.0+cos(n*phi + delta)); 
                K1 = -n*k*sin(n*phi + delta);
            }
            else {
                double diff = phi-delta;
                if (diff < -M_PI) diff += 2.0*M_PI;
                else if (diff > M_PI) diff -= 2.0*M_PI;
                K = k*diff*diff;
                K1 = 2.0*k*diff;
            }
            energy += K;

            // forces
            if (fabs(sin_phi) > 0.1) {
                K1 = K1/sin_phi;
                f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
                f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
                f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);

                f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
                f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
                f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);

                f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
                        + r34.y*dcosdB.z - r34.z*dcosdB.y);
                f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
                        + r34.z*dcosdB.x - r34.x*dcosdB.z);
                f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
                        + r34.x*dcosdB.y - r34.y*dcosdB.x);
            }
            else {
                //  This angle is closer to 0 or 180 than it is to
                //  90, so use the cos version to avoid 1/sin terms
                K1 = -K1/cos_phi;

                f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
                    - r23.x*r23.y*dsindC.y
                    - r23.x*r23.z*dsindC.z);
                f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
                    - r23.y*r23.z*dsindC.z
                    - r23.y*r23.x*dsindC.x);
                f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
                    - r23.z*r23.x*dsindC.x
                    - r23.z*r23.y*dsindC.y);

                f3 += cross(K1,dsindB,r23);

                f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
                    +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
                    +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
                    +dsindB.z*r34.y - dsindB.y*r34.z);
                f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
                    +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
                    +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
                    +dsindB.x*r34.z - dsindB.z*r34.x);
                f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
                    +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
                    +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
                    +dsindB.y*r34.x - dsindB.x*r34.y);
            }
        }    // end loop over multiplicity
        f[dihedral->atom1] += f1;
        f[dihedral->atom2] += f2-f1;
        f[dihedral->atom3] += f3-f2;
        f[dihedral->atom4] += -f3;

        dihedral++; 
    }
    return energy;
}

double Bonded::compute_impropers(const Vector *coords, Vector *f) const {
  
    double energy = 0.0;
    ImproperElem *improper = impropers;
    for (int i=0; i<nimpropers; i++) {
        const Vector *pos0 = coords + improper->atom1;
        const Vector *pos1 = coords + improper->atom2;
        const Vector *pos2 = coords + improper->atom3;
        const Vector *pos3 = coords + improper->atom4;
        const Vector r12 = *pos0 - *pos1;
        const Vector r23 = *pos1 - *pos2;
        const Vector r34 = *pos2 - *pos3;

        Vector dcosdA;
        Vector dcosdB;
        Vector dsindC;
        Vector dsindB;
        Vector f1, f2, f3;

        Vector A, B, C;
        A.cross(r12, r23);
        B.cross(r23, r34);
        C.cross(r23, A);

        double rA = A.length();
        double rB = B.length(); 
        double rC = C.length();

        double cos_phi = (A*B)/(rA*rB);
        double sin_phi = (C*B)/(rC*rB);

        // Normalize B
        rB = 1.0/rB;
        B *= rB;

        double phi = -atan2(sin_phi, cos_phi);

        if (fabs(sin_phi) > 0.1) {
            // Normalize A
            rA = 1.0/rA;
            A *= rA;
            dcosdA = rA*(cos_phi*A-B);
            dcosdB = rB*(cos_phi*B-A);
        }
        else {
            // Normalize C
            rC = 1.0/rC;
            C *= rC;
            dsindC = rC*(sin_phi*C-B);
            dsindB = rB*(sin_phi*B-C);
        }
 
        int mult = improper->multiplicity; 
        for (int j=0; j<mult; j++) {
            double k = improper->k[j];
            double n = improper->n[j];
            double delta = improper->delta[j];
            double K, K1;
            if (n) {
                K = k * (1.0+cos(n*phi + delta)); 
                K1 = -n*k*sin(n*phi + delta);
            }
            else {
                double diff = phi-delta;
                if (diff < -M_PI) diff += 2.0*M_PI;
                else if (diff > M_PI) diff -= 2.0*M_PI;
                K = k*diff*diff;
                K1 = 2.0*k*diff;
            }
            energy += K;

            // forces
            if (fabs(sin_phi) > 0.1) {
                K1 = K1/sin_phi;
                f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
                f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
                f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);

                f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
                f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
                f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);

                f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
                     + r34.y*dcosdB.z - r34.z*dcosdB.y);
                f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
                    + r34.z*dcosdB.x - r34.x*dcosdB.z);
                f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
                    + r34.x*dcosdB.y - r34.y*dcosdB.x);
            }
            else {
                //  This angle is closer to 0 or 180 than it is to
                //  90, so use the cos version to avoid 1/sin terms
                K1 = -K1/cos_phi;

                f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
                    - r23.x*r23.y*dsindC.y
                    - r23.x*r23.z*dsindC.z);
                f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
                    - r23.y*r23.z*dsindC.z
                    - r23.y*r23.x*dsindC.x);
                f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
                    - r23.z*r23.x*dsindC.x
                    - r23.z*r23.y*dsindC.y);

                f3 += cross(K1,dsindB,r23);

                f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
                    +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
                    +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
                    +dsindB.z*r34.y - dsindB.y*r34.z);
                f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
                    +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
                    +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
                    +dsindB.x*r34.z - dsindB.z*r34.x);
                f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
                    +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
                    +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
                    +dsindB.y*r34.x - dsindB.x*r34.y);
            }
        }    // end loop over multiplicity
        f[improper->atom1] += f1;
        f[improper->atom2] += f2-f1;
        f[improper->atom3] += f3-f2;
        f[improper->atom4] += -f3;

        improper++; 
    }
    return energy;
}
