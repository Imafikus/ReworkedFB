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
    alpha = new double*[numberOfObservedVars];
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
            alpha[i] = new double[numberOfPossibleStatesZ];
}
void Model::initializeBeta()
{
    beta = new double*[numberOfObservedVars];
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
            beta[i] = new double[numberOfPossibleStatesZ];
}
void Model::initializeGamma()
{
    gamma = new double *[numberOfObservedVars];
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        gamma[i] = new double [numberOfPossibleStatesZ];
}
void Model::initializeKsi()
{
    ksi = new double*[numberOfPossibleStatesZ];
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
            ksi[i] = new double[numberOfPossibleStatesZ];
}
void Model::initializeMi()
{
    mi = new double[numberOfObservedVars];

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        mi[i] = new[numberOfPossibleStatesZ] ;
}
void Model::computeAlpha()
{
    for(int k = 0; k < numberOfObservedVars; k++)
    {
        for(int i = 0; i < numberOfPossibleStatesZ; i++)
            alpha[k][0] += P[i] * emissionProbs[X[k]][i];


        for(int i = 1; i < numberOfPossibleStatesZ; i++)
                for(int j = 0; j < numberOfPossibleStatesZ; j++)
                    alpha[k][i] += alpha[k][i-1] * emissionProbs[X[k]][i] * transitionProbs[i][j];
    }
}
void Model::computeBeta()
{
    for(int k = 0; k < numberOfObservedVars; k++)
    {
        beta[k][numberOfPossibleStatesZ-1] = 1;

        for(int i = numberOfPossibleStatesZ-2; i >=0; i--)
            for(int j = 0; j < numberOfPossibleStatesZ; j++)
                beta[k][i] += beta[k][i+1] * emissionProbs[X[k]][i] * transitionProbs[i][j];
    }
}

void Model::computeGamma()
{
    for(int k = 0; k < numberOfObservedVars; k++)
    {
        double zbir = 0;

        for(int i = 0; i < numberOfPossibleStatesZ; i++)
            zbir += alpha[k][i] * beta[k][i];

        double normCoef = 1 / zbir;

        for(int i = 0; i < numberOfPossibleStatesZ; i++)
        {
            gamma[k][i] = alpha[k][i] * beta[k][i] * normCoef;
        }
    }
}

void Model::computeNextP(int currentObservedState)
{
    double denominator = gamma[currentObservedState];
    double numerator = 0;
    for(int i = 0; i <= currentObservedState; i++)
        numerator += gamma[i];
    P[currentObservedState] = denominator / numerator;
}
void Model::computeKsi(int currentObservedState)
{
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        double zbir = 0;

        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            zbir += alpha[currentObservedState][j] * beta[currentObservedState][j] *
                    emissionProbs[X[currentObservedState]][j] * transitionProbs[i][j];

        double normCoef = 1 / zbir;

        for(int j = 0; j < numberOfPossibleStatesZ; j++)
        {
            ksi[i][j] = alpha[currentObservedState][i] * beta[currentObservedState][j] *
                    emissionProbs[X[currentObservedState]][j] * transitionProbs[i][j] * normCoef;
            cout<<ksi[i][j] <<" ";
        }
        cout <<endl;

    }
}
void Model::computeCurrentT(int currentObservedState)
{
    double br = 0;
    double im1 = 0;
    double im2 = 0;

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        for(int j = 0; j <= currentObservedState; i++)
            im1 += ksi[i][j];

        for(int l = 0; l < currentObservedState; l++)
            for(int j = 0; j <numberOfPossibleStatesZ; j++)
                im2 += ksi[l][j];

        for(int j = 0; j < numberOfPossibleStatesZ; j++)
        {
            br = ksi[i][j];
            transitionProbs[i][j] = br / (im1 + im2);
        }
    }
}
void Model::computeMiForCurrentE(int currentObservedState)
{
    double br = 0;
    double im = 0;
    for(int i = 0; i < numberOfObservedVars)





}

//computeMiForCurrentE
//computeCurrentE

//PUBLIC FUNCTIONS

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
            computeAlpha();
            computeBeta();
            computeGamma(currentObservedState);
            computeKsi(currentObservedState);
            printGamma();
       // }


}

