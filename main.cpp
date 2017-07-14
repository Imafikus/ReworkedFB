#include <iostream>
#include "model.h"
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


    cout << observedVars << " " << statesZ << " " << statesX << endl;

    for(int i = 0; i < observedVars; i++)
        cout << X[i] << " ";
    cout << endl;

    for(int i = 0; i < statesZ; i++)
        cout << Pi[i] << " ";
    cout << endl;

    for(int i = 0; i < statesX; i++)
        for(int j = 0; j < statesZ; j++)
            cout << E[i][j] << " ";
        cout << endl;

    for(int i = 0; i < statesZ; i++)
        for(int j = 0; j < statesZ; j++)
            cout << T[i][j] << " ";
        cout << endl;

    cout << k << endl;

}



    /*Model model(observedVars,statesZ, statesX, X, Pi, E, T, k);
    model.train();
    return 0;
}*/
