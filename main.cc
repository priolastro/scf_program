// MINIMAL BASIS STO-3G CALCULATION ON HEH+

// THIS IS A LITTLE DUMMY MAIN PROGRAM WHICH CALLS HFCALC
// Attila Szabo and Neil S. Ostlund
// Ed. 2nd (1989) Dover Publications INC.

#include <iostream>
#include <iomanip>
#include <fstream>
#include<vector>

#include "hfcalc.h"
#ifndef __hfcalc__
#define __hfcalc__


using namespace std;

int main(int argc, char** argv) {
    
    int N;   //sto-Ng calculation 
    double R; //bondlength in au
    double Z1,Z2; //Slater orbital exponents
    double ZA,ZB; //Atomic numbers 


    HFCALC(N,R,Z1,Z2,ZA,ZB);
    
    return 0
}