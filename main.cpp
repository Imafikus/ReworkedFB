
#include <iostream>
#include "model.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <time.h>
#include <stdlib.h>

const int MAX_SIZE_X = 10000;

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

void getInitialValues(int &observedVars, int &statesZ, int &statesX, int * &X, double * &Pi, double ** &E, double ** &T, int &iter)
{
    cout << "Give states Z:" << endl;
    cin >> statesZ;

    cout << "Give number of iterations:" << endl;
    cin >> iter;

    ifstream inf("input.txt");
    inf >> observedVars;
    //inf >> statesZ;
    inf >> statesX;
    //cout << "States Z = " << statesZ
    cout << " states X = " << statesX << endl;

    initialization(observedVars,statesZ, statesX, X, Pi, E, T);

    cout << "Loading X\n";
    for(int i = 0; i < observedVars; i++) {
        inf >> X[i];
        cout << "X[i] " << X[i] << endl;
	}
    cout << "Loading P\n";
    for(int i = 0; i < statesZ; i++){
        Pi[i] = 1.0 / statesZ;
    }


    cout << "Loading T" << endl;
    for(int i = 0; i < statesZ; i++)
    {
        for(int j = 0; j < statesZ; j++)
        {
           T[i][j] = rand() % 20;
        }
    }
    cout << "Loading E\n";
    for(int i = 0; i < statesX; i++)
    {
        for(int j = 0; j < statesZ; j++)
        {
            E[i][j] = rand() % 20;
        }
    }

    //inf >> k;
    //cout << "k " << k << endl;

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
    int iter;


    getInitialValues(observedVars,statesZ, statesX, X, Pi, E, T, iter);

    Model model(observedVars, statesZ, statesX, X, Pi, T, E, iter);

    /*X[observedVars] = 100;
    cout << X[observedVars] << endl;*/ //observedVars-1 je poslednji index niza X pre dodavanja
    //
    cout << "Broj iteracija: " << endl;
    model.printNumberOfIterations();

    cout << "Broj stanja X: " << endl;
    model.printNumberOfPossibleStatesX();

    cout << "Broj stanja Z: " << endl;
    model.printNumberOfPossibleStatesZ();

    cout << "Emisiona matrica: " << endl;
    model.printEmission();

    cout << "Tranziciona matrica: " << endl;
    model.printTrans();

    cout << "Niz pi: " << endl;
    model.printP();

    cout << "Niz X" << endl;
    model.printX();


    /*model.testPi();


    model.predict();
    cout << endl;

    model.printP();
    cout << "stampam pi" << endl;

    model.printTrans();
    cout << "stampam trrans" << endl;

    model.printEmission();
    cout << "stampam  E" << endl;*/
}


