/*
 * File:   main.cpp
 * Author: amin
 *
 * Created on September 3, 2015, 1:33 PM
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <new>
#include "Configure.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Initial.h"
#include "Integrator.h"
#include "Nonbonded.h"
#include "Analysis.h"

using namespace std;

int main(int argc, char** argv) {
    time_t t0 = time(0);
    char* cht = ctime(&t0);
    cout << "The local date and time is: " << cht ;

    if (argc != 2){
        cout << "Error: Enter Input file" << endl;
        exit(1);
    }

    Configure conf(argv[1]);  // Read Configure file
    cout << "Input files : " << conf.pdbname << "  " << conf.psfname << "  " << conf.parname << endl;

    // check if the dcd file exists
    /*if (ifstream(conf.dcdname)) {
        cout << "The dcd file exist." << endl;
        cout << "Please backup your dcd file and run your simulation again" << endl;
        exit(1);
    }*/

    if (ifstream(conf.resname)) {
        cout << "Simulation start from restart file" << endl;
    } else {
        cout << "Simulation start from initial configuration with random velocity" << endl;
    }

    cout << "SIMULATION PARAMETERS:" << endl;
    cout << "Number of Atoms   = " << conf.numAtoms << endl;
    cout << "Number of bonds   = " << conf.numBonds << endl;
    cout << "Number of Angles  = " << conf.numAngles << endl;

    Molecule mol(&conf); // Read PSF and PDB files
    Parameters params(&conf); // Read Parameter files

    Initial init(&conf,&mol,&params); // Prepare initial parameteres for run
    // We don't need mol anymore and it would be optimized the code if we delete mol but I don't know how :(

    cout << "Number of Residue = " << init.nresidue << endl;


    Integrator run(&conf,&init); // Main class for MD code

    if (strcmp(conf.mode,"dpd") != 0) {
        if (strncasecmp(conf.analysis, "on", 2) == 0){
            Analysis analysis(&conf,&init);
            analysis.run(&conf, &init);
        }
        else {
            run.Loop(&conf, &init);
        }
    } else {
        cout << "DPD simulation is running" << endl;
        cout << "Cutoff distance must be 1.0" << endl;
        cout << "DPD friction, gamma = " << init.gamma << endl;
        cout << "DPD random,   sigma = " << init.sigma << endl;
        run.Loop_dpd(&conf, &init);
    }

    // Simulation Time
    time_t t1 = time(0);
    cht = ctime(&t1);
    cout << "The local date and time is: " << cht ;
    time_t t2 = t1 - t0;
    cout << "Simulation Time: " << t2 << " seconds" << endl;
    struct tm *tt = gmtime ( &t2 );
    printf ("Simulation Time:  %2d:%02d:%02d\n", (tt->tm_hour)%24, tt->tm_min, tt->tm_sec);
    cout << "END PROGRAM" << endl;
    return 0;
}
