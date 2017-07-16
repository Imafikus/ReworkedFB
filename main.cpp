
#include <iostream>
#include "model.h"
#include "unittest.h"
#include <vector>
#include <algorithm>
#include<fstream>

using namespace std;
void initialization(int &observedVars, int &statesZ, int &statesX, int * &X, double * &Pi, double ** &E, double ** &T)
{
    X = new int[observedVars];
    Pi = new double[statesZ];

    cout << "Making array E of size " << statesX << " x " << statesZ << endl;
    E = new double*[statesX];
    T = new double*[statesZ];

    for(int i = 0; i < statesZ; i++)
        {

            T[i] = new double[statesZ];
        }
    for (int i = 0; i<statesX; i++) E[i] = new double[statesZ];
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
	cout << X[i] << endl;
	}
    cout << "Loading P\n";
    for(int i = 0; i < statesZ; i++){
        inf >> Pi[i];
	cout << Pi[i] << endl;
    }
    cout << "Loading T\n";
    for(int i = 0; i < statesZ; i++)
        for(int j = 0; j < statesZ; j++) {
            inf >> T[i][j];
	    cout << T[i][j] << endl;
	}
    cout << "Loading E\n";
    for(int i = 0; i < statesX; i++)
        for(int j = 0; j < statesZ; j++) {
	    cout << i << " " << j << endl;
            inf >> E[i][j];
	    cout << E[i][j] << endl;
	}

    inf >> k;
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
    cout << "aaa\n";



    Model model(observedVars, statesZ, statesX, X, Pi, T, E, k);

   /* model.printNumberOfVars();
    cout << endl;
    model.printNumberOfPossibleStatesZ();
    cout << endl;
    model.printNumberOfPossibleStatesX();
    cout << endl;
    model.printX();
    cout << endl;
    model.printP();
    cout << endl;
    model.printEmission();
    cout << endl;
    model.printTrans();
    cout << endl;
    model.printNumberOfIterations();
    cout << endl;*/


    //model.printX();


    model.testPi();
/*
    model.printAlpha();
    cout << endl;
    model.printBeta();
    cout << endl;
    model.printGamma();
    cout << endl;
    model.printP();
    cout << endl;
    model.printKsi();
    cout << endl;
    model.printTrans();
    cout << endl;
    model.printMi();
    cout << endl;
    model.printEmission();
    cout << endl;
    model.printNumberOfIterations();
    cout << endl;
*/
}


