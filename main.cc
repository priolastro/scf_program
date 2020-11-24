// MINIMAL BASIS STO-3G CALCULATION ON HEH+

// THIS IS A LITTLE DUMMY MAIN PROGRAM WHICH CALLS HFCALC
// Attila Szabo and Neil S. Ostlund
//https://github.com/lcb/szabo.py/blob/master/szabo.py


#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include<vector>

using namespace std;


double S_integral(double& alpha, double& beta, double& Rab2) {
    double S=pow((M_PI)/(alpha+beta),1.5)*exp(-(alpha*beta)/(alpha+beta)*Rab2);
    return S;
}

double Integral(int& N, double& R, double& Z1, double& Z2, double& ZA, double& ZB) {
    

    //contraction coefficient and exponent for normalized slater orbital
    // for STO-1G alpha=0.27 and coeff=1, ecc...
    vector<vector<double>> coef{
        {1, 0.678914, 0.444635},
        {0, 0.430129, 0.535328},
        {0, 0, 0.154329},
    };
    vector<vector<double>> expon{
        {0.270950,0.151623, 0.109818},
        {0, 0.851819, 0.405771},
        {0, 0, 2.22766},
    };

    //scaling contraction coefficients
    vector<double> alpha1(N,0);
    vector<double> alpha2(N,0);
    vector<double> cont1(N,0);
    vector<double> cont2(N,0);

    for (int i=0; i<N; i++){
        alpha1[i]=expon[i][N-1]*pow(Z1, 2);
        alpha2[i]=expon[i][N-1]*pow(Z2, 2);
        cont1[i]=coef[i][N-1]*pow(2*alpha1[i]/M_PI,0.75);
        cont2[i]=coef[i][N-1]*pow(2*alpha2[i]/M_PI,0.75);
    }

    double Rab2=pow(R,2);
    double S12=0;
    for (int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            S12+=S_integral(alpha1[i], alpha2[j], Rab2)* cont1[i] * cont2[j];
        }
    }
    return S12;
}

int main(int argc, char** argv) {
    
    cout << "Program for the calculation of scf energy of a biatomic molecule.\n";

    int N=3;   //sto-Ng calculation 
    double R=1.4632; //bondlength in au
    double Z1=2.0925,Z2=1.24; //Slater orbital exponents
    double ZA=2,ZB=1; //Atomic numbers 

    cout << "Input parameters:\n";
    cout << N << " " << R << " " << Z1 << " " << Z2 << " " << ZA << " " << ZB << " \n\n" << endl;

    // HartreeFook_Calculation(N,R,Z1,Z2,ZA,ZB);

    double S12=Integral(N, R, Z1, Z2, ZA, ZB);
    
    cout << "S12= " << S12 << endl;
    
    return 0;
}