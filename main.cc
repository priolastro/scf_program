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
    if (arg < 1.0E-6){
        return 1.0-arg/3.0; 
    }
    double F0 = 0.5*pow(M_PI/arg,0.5)*erf(pow(arg,0.5));
    return F0;
}

double S_integral(double& alpha, double& beta, double& Rab2) {
    double S = pow((M_PI)/(alpha+beta),1.5)*exp(-(alpha*beta)/(alpha+beta)*Rab2);
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

double TE_integral(double A, double B, double C, double D, double RAB2, double RCD2, double RPQ2){
    double TE = 2.0*pow(M_PI,2.5)/((A+B)*(C+D)*pow(A+B+C+D,0.5))*Function0((A+B)*(C+D)*RPQ2/(A+B+C+D))*exp(-A*B*RAB2/(A+B)-C*D*RCD2/(C+D));
    return TE;
}

double Integral(int& N, double& R, double& Z1, double& Z2, double& ZA, double& ZB, double& S12, double& T11, double& T12, double& T22, double& V11A, double& V12A, double& V22A, double& V11B, double& V12B, double& V22B, double& V1111, double& V2111, double& V2121, double& V2211, double& V2221, double& V2222) {
    
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
        //to fit slater function with orbital exponent different than 1 (equation 3.224) 
        alpha1[i]=expon[i][N-1]*pow(Z1, 2);
        alpha2[i]=expon[i][N-1]*pow(Z2, 2);
        cont1[i]=coef[i][N-1]*pow(2*alpha1[i]/M_PI,0.75);
        cont2[i]=coef[i][N-1]*pow(2*alpha2[i]/M_PI,0.75);
    }

    double Rab2=pow(R,2);

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

    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            for (int k=0; k<N; k++){
                for (int l=0; l<N; l++){
                    double Rap=alpha2[i]*R/(alpha2[i]+alpha1[j]);
                    double Rbp=R-Rap;
                    double Raq=alpha2[k]*R/(alpha2[k]+alpha1[l]);
                    double Rbq=R-Raq;
                    double Rpq=Rap-Raq;
                    double Rap2=pow(Rap,2);
                    double Rbp2=pow(Rbp,2);
                    double Raq2=pow(Raq,2);
                    double Rbq2=pow(Rbq,2);
                    double Rpq2=pow(Rpq,2);
                    V1111+=TE_integral(alpha1[i], alpha1[j], alpha1[k], alpha1[l], 0, 0, 0)*cont1[i]*cont1[j]*cont1[k]*cont1[l];
                    V2111+=TE_integral(alpha2[i], alpha1[j], alpha1[k], alpha1[l], Rab2, 0, Rap2)*cont2[i]*cont1[j]*cont1[k]*cont1[l];
                    V2121+=TE_integral(alpha2[i], alpha1[j], alpha2[k], alpha1[l], Rab2, Rab2, Rpq2)*cont2[i]*cont1[j]*cont2[k]*cont1[l];
                    V2211+=TE_integral(alpha2[i], alpha2[j], alpha1[k], alpha1[l], 0, 0, Rab2)*cont2[i]*cont2[j]*cont1[k]*cont1[l];
                    V2221+=TE_integral(alpha2[i], alpha2[j], alpha2[k], alpha1[l], 0, Rab2, Rbq2)*cont2[i]*cont2[j]*cont2[k]*cont1[l];
                    V2222+=TE_integral(alpha2[i], alpha2[j], alpha2[k], alpha2[l], 0, 0, 0)*cont2[i]*cont2[j]*cont2[k]*cont2[l];
                }
            }
        }
    }
    return 0;
}

double Collect(double& S12, double& T11, double& T12, double& T22, double& V11A, double& V12A, double& V22A, double& V11B, double& V12B, double& V22B, double& V1111, double& V2111, double& V2121, double& V2211, double& V2221, double& V2222, vector<vector<double>>& H_core, vector<vector<double>>& S_mat){

    for (int i=0; i<2; i++){
        H_core[i].resize(2);
        S_mat[i].resize(2);
    }
    
    // core Hamiltonian
    H_core[0][0]=T11+V11A+V11B;
    H_core[0][1]=H_core[1][0]=T12+V12A+V12B;
    H_core[1][1]=T22+V22A+V22B;

    // S matrix
    S_mat[0][0]=S_mat[1][1]=(double)1;
    S_mat[0][1]=S_mat[1][0]=S12;
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
    double V1111=0;
    double V2111=0;
    double V2121=0;
    double V2211=0;
    double V2221=0;
    double V2222=0;

    Integral(N, R, Z1, Z2, ZA, ZB, S12, T11, T12, T22, V11A, V12A, V22A, V11B, V12B, V22B, V1111, V2111, V2121, V2211, V2221, V2222);

    vector<vector<double>> H_core(2);
    vector<vector<double>> S_mat(2);


    Collect(S12, T11, T12, T22, V11A, V12A, V22A, V11B, V12B, V22B, V1111, V2111, V2121, V2211, V2221, V2222, H_core, S_mat);

    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
            cout << S_mat[i][j] << endl;
        }
    }

    return 0;
}