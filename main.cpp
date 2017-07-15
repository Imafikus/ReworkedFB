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

    E = new double*[statesX];
    T = new double*[statesZ];

    for(int i = 0; i < statesZ; i++)
        {
            E[i] = new double[statesZ];
            T[i] = new double[statesZ];
        }
}

void getInitialValues(int &observedVars, int &statesZ, int &statesX, int * &X, double * &Pi, double ** &E, double ** &T, int &k)
{

    ifstream inf("input.txt");
    inf >> observedVars;
    inf >> statesZ;
    inf >> statesX;

    initialization(observedVars,statesZ, statesX, X, Pi, E, T);

    for(int i = 0; i < observedVars; i++)
        inf >> X[i];

    for(int i = 0; i < statesZ; i++)
        inf >> Pi[i];

    for(int i = 0; i < statesZ; i++)
        for(int j = 0; j < statesZ; j++)
            inf >> T[i][j];

    for(int i = 0; i < statesX; i++)
        for(int j = 0; j < statesZ; j++)
            inf >> E[i][j];

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



    Model model(observedVars, statesZ, statesX, X, Pi, T, E, k);

    /*model.printNumberOfVars();
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

    model.printAlpha();
    cout << endl;

    model.printBeta();


    double* testedPi = model.getP();

    for (int i = 0; i < statesZ; i++)
    {
        cout << testedPi[i] << endl;
    }

    //UnitTest tester;

    //tester.testP(testedPi, statesZ);
}

