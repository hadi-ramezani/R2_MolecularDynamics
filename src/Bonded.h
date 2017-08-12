/*
 * File:   Bonded.h
 * Author: Hadi and Amin
 *
 * Created on September 23, 2015, 9:46 AM
 */

#ifndef BONDED_H
#define	BONDED_H

#include "Initial.h"
#include "Vector.h"

class Bonded {
public:
    Bonded();
    void Compute_bond(const Initial *init, const Vector *pos,Vector *const ff, const int num, double &Ebond);
    void Compute_angle(const Initial *init, const Vector *pos,Vector *const ff, const int num, double &Eangle);
    void Compute_angle_ub(const Initial *init, const Vector *pos,Vector *const ff, const int num, double &Eangle);
    Bonded(const Bonded& orig);
    virtual ~Bonded();
private:

};

#endif	/* BONDED_H */

