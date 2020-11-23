



double HFCALC(N,R,Z1,Z2,ZA,ZB) {

    cout << "STO Calculation for biatomic" << endl;

    //calculation of one and two electron integrals

    vector<double> S12;  
    vector<double> T11;
    vector<double> T12;
    vector<double> T22;
    vector<double> V11A;
    vector<double> V12A;
    vector<double> V22A;
    vector<double> V11B;
    vector<double> V12B;
    vector<double> V22B;
    vector<double> V1111;
    vector<double> V2111;
    vector<double> V2121;
    vector<double> V2211;
    vector<double> V2221;
    vector<double> V2222;

    INTGRL(N, R, Z1, Z2, ZA, ZB);

    vector<double> H;
    vector<double> S;
    vector<double> X;
    vector<double> XT;
    vector<double> TT;

    COLECT(IOP, S12, T11, T12, T22, V11A, V12A, V22A, V11B, V12B, V22B, V1111, V2111, V2121, V2211, V2221, V2222);
    
    SCF(R, ZA, ZB, H, S, X, XT, TT);


}
