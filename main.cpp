
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

void getInitialValues(int &observedVars, int &statesZ, int &statesX, int * &X, double * &Pi, double ** &E, double ** &T, int &iter, string &input)
{
    srand(time(NULL));
    cout << "Give states Z:" << endl;
    cin >> statesZ;

    cout << "Give number of iterations:" << endl;
    cin >> iter;

    cout << "Give states X:" << endl;
    cin >> statesX;

    ifstream inf(input);
    inf >> observedVars;
    //inf >> statesX;
    //cout << " states X = " << statesX << endl;

    initialization(observedVars,statesZ, statesX, X, Pi, E, T);

    cout << "Loading X\n";
    for(int i = 0; i < observedVars; i++)
    {
        inf >> X[i];
        cout << "X[i] " << X[i] << endl;
	}

    cout << "Loading P\n";
    for(int i = 0; i < statesZ; i++)
    {
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

    //normalizing T
    double *niz = new double[statesZ];
    for(int i = 0; i < statesZ; i++)
    {
        double zbir = 0;
        for(int j = 0; j < statesZ; j++)
            zbir += T[i][j];
        niz[i] = zbir;
    }

    for(int i = 0; i < statesZ; i++)
    {
        for(int j = 0; j < statesZ; j++)
            T[i][j] /= niz[i];
    }

    cout << "Loading E\n";
    for(int i = 0; i < statesX; i++)
    {
        for(int j = 0; j < statesZ; j++)
        {
            E[i][j] = rand() % 20;
        }
    }
    delete [] niz;

    //normalizing E
    niz = new double[statesX];
    for(int i = 0; i < statesX; i++)
    {
        double zbir = 0;
        for(int j = 0; j < statesZ; j++)
        {
            zbir += E[j][i];
        }
        niz[i] = zbir;
    }

    for(int i = 0; i < statesX; i++)
    {
        for(int j = 0; j < statesZ; j++)
        {
            E[j][i] /= niz[i];
        }
    }




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
    string input = "input.txt";


    getInitialValues(observedVars,statesZ, statesX, X, Pi, E, T, iter, input);

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


    model.testPi();


    model.predict();
    cout << endl;

    /*model.printP();
    cout << "stampam pi" << endl;

    model.printTrans();
    cout << "stampam trrans" << endl;

    model.printEmission();
    cout << "stampam  E" << endl;*/
}


