#include "model.h"

#include <algorithm>
#include <iostream>

using namespace std;

void Model::initializeT()
{
    transitionProbs = new double*[numberOfPossibleStatesZ];
    for(int i = 1; i <= numberOfPossibleStatesZ; i++)
    {
        transitionProbs[i] = new double[numberOfPossibleStatesZ];
        for(int j = 1; j <= numberOfPossibleStatesZ; j++)
            transitionProbs[i][j] = 1 / numberOfPossibleStatesZ;
    }
}

void Model::initializeE()
{
    emissionProbs = new double*[numberOfPossibleStatesX];
    for(int i = 1; i <= numberOfPossibleStatesX; i++)
    {
        emissionProbs[i] = new double[numberOfPossibleStatesZ];
        for(int j = 1; j <= numberOfPossibleStatesZ; j++)
            emissionProbs[i][j] = 1 / numberOfPossibleStatesX;
    }
}

void Model::initializeP()
{
    P = new double[numberOfPossibleStatesZ];
    for(int i = 0; i <= numberOfPossibleStatesZ; i++)
        P[i] = 1 / numberOfPossibleStatesZ;
}
void Model::initializeAlpha()
{
    alpha = new double*[numberOfPossibleStatesZ];
    for(int i = 1; i <= numberOfPossibleStatesZ; i++)
    {
        alpha[i] = new double[numberOfPossibleStatesZ];
        for(int j = 1; j <= numberOfPossibleStatesZ; j++)
            alpha[i][j] = 0;
    }
}
void Model::initializeBeta()
{
    beta = new double*[numberOfPossibleStatesZ];
    for(int i = 1; i <= numberOfPossibleStatesZ; i++)
    {
        beta[i] = new double[numberOfPossibleStatesZ];
        for(int j = 1; j <= numberOfPossibleStatesZ; j++)
            beta[i][j] = 0;
    }
}

void Model::computeAlpha()
{
    for(int i = 1; i <= numberOfPossibleStatesZ; i++)
        alpha[1][i] = P[i]*emissionProbs[X[1]][i];

    for(int k = 2; k <= numberOfObservedVars; k++)
        for(int i = 1; i <= numberOfPossibleStatesZ; i++)
            for(int j = 1; j <= numberOfPossibleStatesZ; j++)
                alpha[k][i] += alpha[k-1][i] * emissionProbs[X[k]][j] * transitionProbs[i][j];
}
void Model::computeBeta()
{
    for(int i = 1; i <= numberOfPossibleStatesZ; i++)
        beta[numberOfObservedVars][i] = 1;

    for(int k = numberOfObservedVars-1; k >= 1; k--)
        for(int i = 1; i <= numberOfPossibleStatesZ; i++)
            for(int j = 1; j <= numberOfPossibleStatesZ; j++)
                beta[k][i] += beta[k+1][i] * emissionProbs[X[k+1]][j] * transitionProbs[i][j];
}
void Model::computeNormalized()
{
    normalizedProbs = new double*[numberOfPossibleStatesZ];
    for (int k = 1; k <= numberOfObservedVars; k++)
        {
            normalizedProbs[k] = new double[numberOfPossibleStatesZ];
            float zbir = 0;

            for(int i = 1; i <= numberOfPossibleStatesZ; i++)
                zbir += alpha[k][i] * beta[k][i];

            float x = 1 / zbir;

            for (int i = 1; i <= numberOfPossibleStatesZ; i++)
                normalizedProbs[k][i] = alpha[k][i] * beta[k][i] * x;
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
    for(int i = 1; i <= values.size(); i++)
    {
        X[i] = values.at(i);
    }
}
void Model::printNormalizedMatrix()
{
    for (int i = 1; i <= numberOfObservedVars; i++)
        {
            for (int j = 1; j <= numberOfPossibleStatesZ; j++)
                cout << normalizedProbs[i][j];
            cout << endl;
        }
}
void Model::train()
{
    double gama[1000];
    double ksi[1000];


    /*setNumberOfpossibleStatesZ();
    getNumberOfObservedVars();
    getNumberOfpossibleStatesX();
    getArrayX();*/

    initializeT();
    initializeE();
    initializeP();
    initializeAlpha();
    initializeBeta();

    computeAlpha();
    computeBeta();
    computeNormalized();

    //printNormalizedMatrix();

}

