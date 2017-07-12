#include "model.h"

#include <algorithm>
#include <iostream>

using namespace std;

void Model::initializeT()
{
    transitionProbs = new double*[numberOfPossibleStatesZ];
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        transitionProbs[i] = new double[numberOfPossibleStatesZ];
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            transitionProbs[i][j] = 1.0/numberOfPossibleStatesZ;
    }
}

void Model::initializeE()
{
    emissionProbs = new double*[numberOfPossibleStatesX];
    for(int i = 0; i < numberOfPossibleStatesX; i++)
    {
        emissionProbs[i] = new double[numberOfPossibleStatesZ];
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            emissionProbs[i][j] = 1.0/numberOfPossibleStatesX;
    }
}

void Model::initializeP()
{
    P = new double[numberOfPossibleStatesZ];
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        P[i] = 1.0/numberOfPossibleStatesZ;
}
void Model::initializeAlpha()
{
    alpha = new double[numberOfPossibleStatesZ];
    /*for(int i = 0; i < numberOfPossibleStatesZ; i++)
            alpha[i] = 0;*/
}
void Model::initializeBeta()
{
    beta = new double[numberOfPossibleStatesZ];
    /*for(int i = 0; i < numberOfPossibleStatesZ; i++)
            beta[i] = 0;*/
}
void Model::initializeGamma()
{
    gamma = new double[numberOfPossibleStatesZ];
    /*for(int i = 0; i < numberOfPossibleStatesZ; i++)
        gamma[i] = 0;*/

}
void Model::initializeKsi()
{
    ksi = new double*[numberOfPossibleStatesZ];
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        ksi[i] = new double[numberOfPossibleStatesZ];
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            ksi[i][j] = transitionProbs[i][j] * P[j];
    }
}
void Model::computeAlpha(int currentObservedState)
{
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        alpha[0] += P[i] * emissionProbs[X[currentObservedState]][i];


    for(int i = 1; i < numberOfPossibleStatesZ; i++)
    {
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            {
                alpha[i] += alpha[i-1] * emissionProbs[X[currentObservedState]][i] * transitionProbs[i][j];

            }
    }

}
void Model::computeBeta(int currentObservedState)
{
     beta[numberOfPossibleStatesZ-1] = 1;

        for(int i = numberOfPossibleStatesZ-2; i >=0; i--)
            for(int j = 0; j < numberOfPossibleStatesZ; j++)
                beta[i] += beta[i+1] * emissionProbs[X[currentObservedState]][i] * transitionProbs[i][j];
}
void Model::computeGamma(int currentObservedState)
{

    double zbir = 0;

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        zbir += alpha[i] * beta[i];

    double normCoef = 1 / zbir;

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        gamma[i] = alpha[i] * beta[i] * normCoef;
    }

}

//fillArrayForGamma
//computeKsi
//computeCurrentT
//computeCurrentP
//computeCurrentE

//public functions

void Model::setArrayX(vector<int> &values)
{
    X = new int[values.size()];
    for(int i = 0; i < values.size(); i++)
    {
        X[i] = values.at(i);
    }
}
void Model::printGamma()
{
    for (int j = 0; j < numberOfPossibleStatesZ; j++)
        cout << gamma[j] << " ";
    cout << endl;
}
void Model::printNormalizedMatrix()
{
    for (int i = 0; i < numberOfObservedVars; i++)
        {
            for (int j = 0; j < numberOfPossibleStatesZ; j++)
                cout << normalizedProbs[i][j] << " ";
            cout << endl;
        }

    cout << endl;

}

void Model::train()
{
    int currentObservedState = 0;
    initializeT();
    initializeE();
    initializeP();

    initializeAlpha();
    initializeBeta();
    initializeKsi();


    //for(int currentObservedState = 0; currentObservedState < numberOfPossibleStatesX; currentObservedState++)
      //  {
            computeAlpha(currentObservedState);
            computeBeta(currentObservedState);
            computeGamma(currentObservedState);
            printGamma();
       // }


}

