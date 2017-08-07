/*
 * File:   External.h
 * Author: amin
 *
 * Created on January 18, 2016, 3:39 PM
 */

#ifndef EXTERNAL_H
#define	EXTERNAL_H

#include "iostream"
#include "Vector.h"
#include <vector>
#include "Initial.h"

using namespace std;

class External {
public:

    void Boundary(const int num, const Vector *pos, Vector * const ff);
    void Force(const int num, const Vector *pos, Vector * const ff);
    External();
    External(const External& orig);
    virtual ~External();
private:

};

#endif	/* EXTERNAL_H */

