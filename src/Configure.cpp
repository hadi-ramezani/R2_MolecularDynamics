#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "Configure.h"

using namespace std;

Configure::Configure(char *input) {
    string str;
    char   seg[256];
    char  *pch;
    ifstream infile;
    infile.open(input);
    while(!infile.eof()) {
	getline(infile, str);
	if (str!="") {
            strncpy(seg, str.c_str(), sizeof(seg));
            seg[sizeof(seg) - 1] = 0;
            pch = strtok (seg," ,:=");
            string word(pch);
	    if (word == "structure") {
                char *pch = strtok (NULL, " ,:=");
                strcpy(psfname,pch);
        }
        else if (word == "coordinates") {
            char *pch = strtok (NULL, " ,:=");
            strcpy(pdbname,pch);
        }
        else if (word == "parameters") {
            char *pch = strtok (NULL, " ,:=");
            strcpy(parname,pch);
        }
        else if (word == "outputName") {
            char *pch = strtok (NULL, " ,:=");
            strcpy(outname,pch);
            sprintf(dcdname,"%s.dcd",outname);
            sprintf(resname,"%s.res",outname);
            sprintf(enename,"%s.ene",outname);
        }
        else if (word == "inputName") {
            char *pch = strtok (NULL, " ,:=");
            strcpy(inname,pch);
        }
        else if (word == "mode") {
            char *pch = strtok (NULL, " ,:=");
            strcpy(mode,pch);
        }
        else if (word == "thermostat") {
            char *pch = strtok (NULL, " ,:=");
            strcpy(thermostat,pch);
        }
        else if (word == "temperature") {
            char *pch = strtok (NULL, " ,:=");
            temp = stof(pch);
        }
        else if (word == "settemp") {
            char *pch = strtok (NULL, " ,:=");
            strcpy(settemp,pch);
        }
        else if (word == "lambda") {
            char *pch = strtok (NULL, " ,:=");
            lambda = stof(pch);
        }
        else if (word == "friction") {
            char *pch = strtok (NULL, " ,:=");
            gamma = stof(pch);
        }
        else if (word == "seed") {
            char *pch = strtok (NULL, " ,:=");
            seed = stoi(pch);
        }
        else if (word == "timestep") {
            char *pch = strtok (NULL, " ,:=");
            timestep = stof(pch);
        }
        else if (word == "cutoff"){
            char *pch = strtok (NULL, " ,:=");
            cutoff = stof(pch);
        }
        else if (word == "switchdist"){
            char *pch = strtok (NULL, " ,:=");
            switchdist = stof(pch);
            if (switchdist>cutoff) {
                cout << "ERROR: switchdist must be <= cutoff" <<endl;
                exit(1);
            }
        }
        else if (word == "pairlistdist"){
            char *pch = strtok (NULL, " ,:=");
            pairlistdist = stof(pch);
            if (pairlistdist<cutoff) {
                cout << "ERROR: pairlistdist must be >= cutoff" <<endl;
                exit(1);
            }
        }
        else if (word == "celldist"){
            char *pch = strtok (NULL, " ,:=");
            celldist = stof(pch);
            if (celldist<pairlistdist){
                cout << "ERROR: celldist must be >= pairlistdist" <<endl;
                exit(1);
            }
        }
        else if (word == "pairlistFreq"){
            char *pch = strtok (NULL, " ,:=");
            pairlistFreq = stoi(pch);
        }
        else if (word == "nonbondedFreq"){
            char *pch = strtok (NULL, " ,:=");
            nonbondedFreq = stoi(pch);
        }
        else if (word == "rigidBonds") {
            char *pch = strtok (NULL, " ,:=");
            strcpy(rigidBonds,pch);
        }
        else if (word == "box"){
            for (int ii=0; ii<3; ii++){
                char *pch = strtok (NULL, " ,:=");
                box[ii] = stof(pch);
            }
        }
        else if (word == "wrapFreq"){
            char *pch = strtok (NULL, " ,:=");
            wrapFreq = stoi(pch);
        }
        else if (word == "restartFreq"){
            char *pch = strtok (NULL, " ,:=");
            restartFreq = stoi(pch);
        }
        else if (word == "dcdFreq"){
            char *pch = strtok (NULL, " ,:=");
            dcdFreq = stoi(pch);
        }
        else if (word == "energyFreq"){
            char *pch = strtok (NULL, " ,:=");
            energyFreq = stoi(pch);
        }
        else if (word == "firsttimestep"){
            char *pch = strtok (NULL, " ,:=");
            fstep = stoi(pch);
        }
        else if (word == "run"){
            char *pch = strtok (NULL, " ,:=");
            nsteps = stoi(pch);
        }
        else if (word == "analysis"){
            char* pch = strtok(NULL, " ,:=");
            strcpy(analysis,pch);
            }
        else if (word == "aDcdName"){
            char* pch = strtok(NULL, " ,:=");
            strcpy(aDcdName,pch);
            pch = strtok (NULL, " ,:=");
            if (pch != NULL) dcdFirst = stoi(pch); 
            pch = strtok (NULL, " ,:=");
            if (pch != NULL) dcdLast = stoi(pch); 
            pch = strtok (NULL, " ,:=");
            if (pch != NULL) dcdStep = stoi(pch); 
        }
        else if (word == "aMode"){
            char* pch = strtok(NULL, " ,:=");
            strcpy(aMode,pch);                 
        }
        else if (word == "3bodyCutoff"){
            char *pch = strtok (NULL, " ,:=");
            threebodyCutoff = stof(pch);
        }
        else if (word == "3bodyIJCutoff"){
            char *pch = strtok (NULL, " ,:=");
            threebodyIJCutoff = stof(pch);
        }
        else if (word == "3bodyPairDist"){
            char *pch = strtok (NULL, " ,:=");
            threebodyPairDist = stof(pch);
            cout << threebodyPairDist << endl;
        }
        else if (word == "aOutputFilename"){
            char* pch = strtok(NULL, " ,:=");
            strcpy(aOutputFilename,pch); 
        }
    }
}
    infile.close();
    
    // Check the analysis input
    if (strncasecmp(analysis, "on", 2)==0 && (aDcdName[0] == NULL)) {
        cout << "Error: You must define a dcd file with analysis on..." << endl;
        exit(1);
    }

    // Read the number of atoms, bonds, angles from psf file
    infile.open(psfname);
    size_t found;
    while(!infile.eof()) {
	getline(infile, str);

        found=str.find("NATOM");
        if (found!=std::string::npos) {
            strncpy(seg, str.c_str(), sizeof(seg));
            seg[sizeof(seg) - 1] = 0;
            pch = strtok (seg," ,:=");
            numAtoms = stoi(pch);
        }

        found=str.find("NBOND");
        if (found!=std::string::npos) {
            strncpy(seg, str.c_str(), sizeof(seg));
            seg[sizeof(seg) - 1] = 0;
            pch = strtok (seg," ,:=");
            numBonds = stoi(pch);
        }

        found=str.find("NTHETA");
        if (found!=std::string::npos) {
            strncpy(seg, str.c_str(), sizeof(seg));
            seg[sizeof(seg) - 1] = 0;
            pch = strtok (seg," ,:=");
            numAngles = stoi(pch);
        }

    }
    infile.close();
}

Configure::Configure(const Configure& orig) {
}

Configure::~Configure() {
}
