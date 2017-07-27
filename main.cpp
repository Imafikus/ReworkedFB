
#include <iostream>
#include "model.h"
#include <vector>
#include <sstream>
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

void getInitialValues(int &observedVars, int &statesZ, int &statesX, int * &X, double * &Pi, double ** &E, double ** &T, int &iter, string &input, int &expectedState)
{
    srand(time(NULL));

    statesZ  = 3;
    iter = 100;
    statesX = 5;
    observedVars = 300;

    ifstream inf(input);
    //inf >> observedVars;
    //inf >> statesX;
    //cout << " states X = " << statesX << endl;

    initialization(observedVars,statesZ, statesX, X, Pi, E, T);

    cout << "Loading X\n";
    for(int i = 0; i < observedVars; i++)
    {
        inf >> X[i];
        cout << "X[i] " << X[i] << endl;
	}

	cout << "Loading expectedState\n";
	inf >> expectedState;

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

    cout << "normalizing E\n";
    niz = new double[statesZ];
    for(int i = 0; i < statesZ; i++)
    {
        double zbir = 0;
        for(int j = 0; j < statesX; j++)
        {
            zbir += E[j][i];
        }
        niz[i] = zbir;
    }
    for(int i = 0; i < statesZ; i++)
    {
        for(int j = 0; j < statesX; j++)
        {
            E[j][i] /= niz[i];
        }
    }
}
void runMethodOne()
{
vector<int> v;

    for(int k = 0; k < 400; k++)
    {
        string input = "kretanje";
        int expectedState;
        int observedVars;
        int statesZ;
        int statesX;
        int * X;
        double * Pi;
        double ** E;
        double ** T;
        int iter;

        stringstream ss;
        ss << k;
        input += ss.str() + ".txt";

        getInitialValues(observedVars,statesZ, statesX, X, Pi, E, T, iter, input, expectedState);

        Model model(observedVars, statesZ, statesX, X, Pi, T, E, iter);
        //model.printEmission();

        model.testPi();


        int prediction = model.predict();
        if(prediction == expectedState) v.push_back(1);
        else v.push_back(0);

        delete[] X;
        delete[] Pi;

        for(int i = 0; i < statesZ; i++)
        {
            delete[] T[i];
        }
        delete[] T;

        for(int i = 0; i < statesX; i++)
        {
            delete[] E[i];
        }
        delete[] E;

        cout << k << ". iteracija zavrsena" << endl;
    }
    int good = 0;
    for(int i = 0; i < v.size(); i++)
        if(v.at(i) == 1) good++;
    cout << "pogodjenih: " << good << ":" << v.size() << endl;
}
int main()
{
    runMethodOne();


}


