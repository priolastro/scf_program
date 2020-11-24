// MINIMAL BASIS STO-3G CALCULATION ON HEH+

// THIS IS alpha LITTLE DUMMY MAIN PROGRAM WHICH CALLS HFCALC
// Attila Szabo and Neil S. Ostlund
//https://github.com/lcb/szabo.py/blob/master/szabo.py


#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include<vector>

using namespace std;

double Function0 (double arg){
    if (arg < 1e-6){
        return 1-arg/3; 
    }
    double F0=0.5*pow(M_PI/arg,0.5)*erf(pow(arg,0.5));
    return F0;
}

double S_integral(double& alpha, double& beta, double& Rab2) {
    double S=pow((M_PI)/(alpha+beta),1.5)*exp(-(alpha*beta)/(alpha+beta)*Rab2);
    return S;
}

double T_integral(double& alpha, double& beta, double Rab2){
    double T = alpha*beta/(alpha+beta)*(3-(2*alpha*beta)/(alpha+beta)*Rab2)*pow(M_PI/(alpha+beta),1.5)*exp(-(alpha*beta)/(alpha+beta)*Rab2);
    return T;
}

double V_integral(double& alpha, double& beta, double Rab2, double Rcp2 , double Z){
    double V = -2*M_PI/(alpha+beta)*Z*exp(-(alpha*beta)/(alpha+beta)*Rab2)*Function0((alpha+beta)*Rcp2);
    return V;
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
    double T11=0;
    double T12=0;
    double T22=0;
    double V11A=0;
    double V12A=0;
    double V22A=0;
    double V11B=0;
    double V12B=0;
    double V22B=0;

    for (int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            double Rap=alpha2[j]*R/(alpha1[i]+alpha2[j]); //perche origine in A che uguale a 0 quindi rimane solo beta
            double Rap2=pow(Rap,2);
            double Rbp2=pow(R-Rap,2);
            S12+=S_integral(alpha1[i], alpha2[j], Rab2)* cont1[i] * cont2[j];
            T11+=T_integral(alpha1[i], alpha1[j], 0)*cont1[i]*cont1[j];
            T12+=T_integral(alpha1[i], alpha2[j], Rab2)*cont1[i]*cont2[j];
            T22+=T_integral(alpha2[i], alpha2[j], 0)*cont2[i]*cont2[j];
            V11A+=V_integral(alpha1[i], alpha1[j], 0, 0, ZA) * cont1[i] * cont1[j];
            V12A+=V_integral(alpha1[i], alpha2[j], Rab2, Rap2, ZA) * cont1[i] * cont2[j];
            V22A+=V_integral(alpha2[i], alpha2[j], 0, Rab2, ZA) * cont2[i] * cont2[j];
            V11B+=V_integral(alpha1[i], alpha1[j], 0, Rab2, ZB) * cont1[i] * cont1[j];
            V12B+=V_integral(alpha1[i], alpha2[j], Rab2, Rbp2, ZB) * cont1[i] * cont2[j];
            V22B+=V_integral(alpha2[i], alpha2[j], 0, 0, ZB) * cont2[i] * cont2[j];
        }
    }
    cout << "S12= " << S12 << endl;
    cout << "T11= " << T11 << endl;
    cout << "T12= " << T12 << endl;
    cout << "T22= " << T22 << endl;
    cout << "V11A= " << V11A << endl;
    cout << "V12A= " << V12A << endl;
    cout << "V22A" << V22A << endl;
    cout << "V11B= " << V11B << endl;
    cout << "V12B" << V12B << endl;
    cout << "V22B" << V22B << endl;

    
    return 0;
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

    Integral(N, R, Z1, Z2, ZA, ZB);


    return 0;
}