/*
 * File:   External.cpp
 * Author: amin
 *
 * Created on September 28, 2015, 5:23 PM
 */

#include "External.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>

using namespace std;

External::External() {
    cout << "External force and boundary conditions apply" << endl;
}

void External::Boundary(const int num, const Vector *pos, Vector * const ff){

    // Channel
/*    for (int iatom=0; iatom<num; iatom++){
      if (pos[iatom].y<2.0)  ff[iatom].y = -ff[iatom].y;
      if (pos[iatom].y>10.805)  ff[iatom].y = -ff[iatom].y;
    }
*/
    for (int iatom=0; iatom<num; iatom++){
        if (pos[iatom].y<3.0)  {
          if (pos[iatom].y>2.0) {
            ff[iatom].y += 20*(1-(3.0-pos[iatom].y));
          } else {
            ff[iatom].y += 20;
          }
        }
        if (pos[iatom].y>9.805) {
          if (pos[iatom].y<8.805) {
            ff[iatom].y -= 20*(1-(pos[iatom].y-8.805));
          } else {
            ff[iatom].y -= 20;
          }
        }
    }
}
void External::Force(const int num, const Vector *pos, Vector *const ff){
}

External::External(const External& orig) {
}

External::~External() {
}

