/*
 * File:   Parameters.cpp
 * Author: amin
 *
 * Created on September 9, 2015, 5:16 PM
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include "Parameters.h"
#include "Configure.h"
#include <algorithm>

using namespace std;

Parameters::Parameters(const Configure *conf) {

    read_params(conf->parname);
}

void Parameters::read_params(const char *filename){
    string str;
    string first_word_str;      //  First word of the current line
    size_t found;
    char   first_word[512];
    int skipline; // skip the line
    int skipall = 0; // skip rest of the file
    int  par_type=0;         //  What type of parameter are we currently

    if (!ifstream(filename)) {
        cout << "Unable to open the parameter file" << endl;
        exit(1);
    }
    ifstream infile;
    infile.open(filename);
    while(!infile.eof()) {
        getline(infile, str);
        // Get the first word
        first_word_str = str.substr(0, str.find(" "));
        // convert the first word to char (str won't work with strncmp function)
        strncpy(first_word, first_word_str.c_str(), sizeof(first_word));
        // set this so that the parameters will be stored
        skipline = 0;
        if (str != "" &&
            (strncmp(first_word, "!", 1) != 0) &&
            (strncmp(first_word, "*", 1) != 0) &&
            (strncasecmp(first_word, "END", 3) != 0)) {

            if ( skipall ) {
                cout << "SKIPPING PART OF PARAMETER FILE AFTER RETURN STATEMENT\n" << endl;
                break;
            }
            /*  Now, determine the apropriate parameter type.   */
            if (strncasecmp(first_word, "bond", 4)==0) {
                par_type=1; skipline=1;
            }
            else if (strncasecmp(first_word, "angl", 4)==0) {
                par_type=2; skipline=1;
            }
            else if (strncasecmp(first_word, "thet", 4)==0) {
                par_type=2; skipline=1;
            }
            else if (strncasecmp(first_word, "dihe", 4)==0) {
                par_type=3; skipline=1;
            }
            else if (strncasecmp(first_word, "phi", 3)==0) {
                par_type=3; skipline=1;
            }
            else if (strncasecmp(first_word, "impr", 4)==0) {
                par_type=4; skipline=1;
            }
            else if (strncasecmp(first_word, "imph", 4)==0) {
                par_type=4; skipline=1;
            }
            else if (strncasecmp(first_word, "nbon", 4)==0) {
                par_type=5; skipline=1;
            }
            else if (strncasecmp(first_word, "nonb", 4)==0) {
                par_type=5; skipline=1;
            }
            else if (strncasecmp(first_word, "nbfi", 4)==0) {
                par_type=6; skipline=1;
            }
            else if (strncasecmp(first_word, "hbon", 4)==0) {
                par_type=7; skipline=1;
            }
            else if (strncasecmp(first_word, "cmap", 4)==0) {
                par_type=8; skipline=1;
            }
            else if (strncasecmp(first_word, "nbta", 4)==0) {
                par_type=9; skipline=1;
            }
            else if (strncasecmp(first_word, "thol", 4)==0) {
                par_type=10; skipline=1;
            }
            else if (strncasecmp(first_word, "atom", 4)==0) {
                par_type=11; skipline=1;
            }
            else if (strncasecmp(first_word, "ioformat", 8)==0) {
                skipline=1;
            }
            else if (strncasecmp(first_word, "read", 4)==0) {
                skipline=1;
            }
            else if (strncasecmp(first_word, "return", 4)==0) {
                skipall=1;  skipline=1;
            }
            else if ((strncasecmp(first_word, "nbxm", 4) == 0) ||
                (strncasecmp(first_word, "grou", 4) == 0) ||
                (strncasecmp(first_word, "cdie", 4) == 0) ||
                (strncasecmp(first_word, "shif", 4) == 0) ||
                (strncasecmp(first_word, "vgro", 4) == 0) ||
                (strncasecmp(first_word, "vdis", 4) == 0) ||
                (strncasecmp(first_word, "vswi", 4) == 0) ||
                (strncasecmp(first_word, "cutn", 4) == 0) ||
                (strncasecmp(first_word, "ctof", 4) == 0) ||
                (strncasecmp(first_word, "cton", 4) == 0) ||
                (strncasecmp(first_word, "eps", 3) == 0) ||
                (strncasecmp(first_word, "e14f", 4) == 0) ||
                (strncasecmp(first_word, "wmin", 4) == 0) ||
                (strncasecmp(first_word, "aexp", 4) == 0) ||
                (strncasecmp(first_word, "rexp", 4) == 0) ||
                (strncasecmp(first_word, "haex", 4) == 0) ||
                (strncasecmp(first_word, "aaex", 4) == 0) ||
                (strncasecmp(first_word, "noac", 4) == 0) ||
                (strncasecmp(first_word, "hbno", 4) == 0) ||
                (strncasecmp(first_word, "cuth", 4) == 0) ||
                (strncasecmp(first_word, "ctof", 4) == 0) ||
                (strncasecmp(first_word, "cton", 4) == 0) ||
                (strncasecmp(first_word, "cuth", 4) == 0) ||
                (strncasecmp(first_word, "ctof", 4) == 0) ||
                (strncasecmp(first_word, "cton", 4) == 0) ) {
                if ((par_type != 5) && (par_type != 6) && (par_type != 7) && (par_type != 9)) {
                    cout << "ERROR IN CHARMM PARAMETER FILE " << filename << " in LINE= " << str << endl;
                } 
                else {
                    skipline = 1;
                }
            }        
            else if (par_type == 0) {
                /*  This is an unknown paramter.        */
                cout << "UNKNOWN PARAMETER IN CHARMM PARAMETER FILE " << filename << " in LINE= " << str << endl;
            }
        }
        else {
            skipline = 1;
        }

        if ( (par_type != 0) && (!skipline) ) {
            /*  Now, call the appropriate function based    */
            /*  on the type of parameter we have                */
            /*  I know, this should really be a switch ...  */
            // Note that the code currently only handles bond, angle, and vdw.
            // I included all sections so that the implementation will be easier in the
            // future
            if (par_type == 1) {
                // This is bond parameter
                bond_str.push_back(str);
            }
            else if (par_type == 2) {
                // This is angle parameter
                angle_str.push_back(str);
            }
            else if (par_type == 3) {
                // This is dihedral parameter
                continue;
            }
            else if (par_type == 4) {
                // This is improper parameter
                continue;
            }
            else if (par_type == 5) {
                // This is vdw parameter
                vdw_str.push_back(str);
            }
            else if (par_type == 6) {
                // This is vdw pair parameter
                continue;
            }
            else if (par_type == 7) {
                // This is hb pair parameter
                continue;
            }
            else if (par_type == 8) {
                // This is crossterm parameter
                continue;
            }
            else if (par_type == 9) {

                // This is table pair parameter
                continue;
            }
            else if (par_type == 10) {
                // This is nbthole parameter
                continue;
            }
            else if (par_type == 11) {
                if ( strncasecmp(first_word, "mass", 4) != 0 ) {
                    cout << "UNKNOWN PARAMETER IN CHARMM PARAMETER FILE " << filename << " in LINE= " << str << endl;
                }
            }
        else {
            /*  This really should not occour!      */
            /*  This is an internal error.          */
            /*  This is VERY BAD                        */
            cout << "INTERNAL ERROR IN CHARMM PARAMETER FILE " << filename << "LINE = " << str << endl;
        }
        }
    }
}

Parameters::Parameters(const Parameters& orig) {
}

Parameters::~Parameters() {
}
