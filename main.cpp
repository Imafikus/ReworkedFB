
#include <iostream>
#include "model.h"
#include "unittest.h"
#include <vector>
#include <algorithm>
#include<fstream>
const int MAX_SIZE_X = 1000;

using namespace std;
void initialization(int &observedVars, int &statesZ, int &statesX, int * &X, double * &Pi, double ** &E, double ** &T)
{
    X = new int[MAX_SIZE_X];
    Pi = new double[statesZ];

    cout << "Making array E of size " << statesX << " x " << statesZ << endl;
    E = new double*[statesX];
    T = new double*[statesZ];

    for(int i = 0; i < statesZ; i++)
    {
        T[i] = new double[statesZ];
    }
    for (int i = 0; i<statesX; i++)
        E[i] = new double[statesZ];
}

void getInitialValues(int &observedVars, int &statesZ, int &statesX, int * &X, double * &Pi, double ** &E, double ** &T, int &k)
{
    ifstream inf("input.txt");
    inf >> observedVars;
    inf >> statesZ;
    inf >> statesX;
    cout << "States Z = " << statesZ << " states X = " << statesX << endl;

    initialization(observedVars,statesZ, statesX, X, Pi, E, T);

    cout << "Loading X\n";
    for(int i = 0; i < observedVars; i++) {
        inf >> X[i];
        cout << "X[i] " << X[i] << endl;
	}
    cout << "Loading P\n";
    for(int i = 0; i < statesZ; i++){
        inf >> Pi[i];
        cout << "Pi[i] " << Pi[i] << endl;
    }
    cout << "Loading T" << endl;
    for(int i = 0; i < statesZ; i++) {
        for(int j = 0; j < statesZ; j++) {
            inf >> T[i][j];
            cout << "T[i][j] " << T[i][j] << endl;
        }
    }
    cout << "Loading E\n";
    for(int i = 0; i < statesX; i++)
        for(int j = 0; j < statesZ; j++) {
            cout << i << " " << j << endl;
            inf >> E[i][j];
            cout << "E[i][j] " << E[i][j] << endl;
	}

    inf >> k;
    cout << "k " << k << endl;

}

int main()
{
    int observedVars;
    int statesZ;
    int statesX;
    int * X;
    double * Pi;
    double ** E;
    double ** T;
    int k;


    getInitialValues(observedVars,statesZ, statesX, X, Pi, E, T, k);

    Model model(observedVars, statesZ, statesX, X, Pi, T, E, k);

/*model.testPi();

    /*model.printP();
    cout << "stampam pi" << endl;

    model.printTrans();
    cout << "stampam trrans" << endl;

    model.printEmission();
    cout << "stampam  E" << endl;*/

    model.predict();
    cout << endl;
}


