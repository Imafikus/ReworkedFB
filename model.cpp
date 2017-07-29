#include "model.h"
#include "unittest.h"
#include <algorithm>
#include <iostream>
#include <time.h>
#include <stdlib.h>


using namespace std;

void Model::initializeC()
{
    C = new double[MAX_SIZE_FOR_ARRAYS];
}
void Model:: eraseAlpha()
{
    for(int i = 0; i < numberOfObservedVars; i++)
    {
        delete[] alpha[i];
    }
    delete[] alpha;
}
void Model::eraseC()
{
    delete [] C;
}
void Model::initializeAlpha()
{
    alpha = new double*[MAX_SIZE_ALPHA];
    for(int i = 0; i < MAX_SIZE_ALPHA; i++)
        alpha[i] = new double[numberOfPossibleStatesZ];

    for(int i = 0; i < MAX_SIZE_ALPHA; i++)
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            alpha[i][j] = 0;

}
void Model::initializeBeta()
{
    beta = new double*[numberOfObservedVars];
    for(int i = 0; i < numberOfObservedVars; i++)
            beta[i] = new double[numberOfPossibleStatesZ];

    for(int i = 0; i < numberOfObservedVars; i++)
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            beta[i][j] = 0;
}
void Model::initializeRandomSeed()
{
    randomSeed = new double[3];
}
void Model::computeAlpha()
{
    double sumC = 0;
    for (int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        alpha[0][i] = P[i] * emissionProbs[X[0]][i];
        sumC += alpha[0][i];
    }

    C[0] = 1.0 / sumC;

    for (int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        alpha[0][i] /= sumC;
    }

    for(int k = 1; k < numberOfObservedVars; k++)
    {
        sumC = 0;

        for(int i = 0; i < numberOfPossibleStatesZ; i++)
        {
            double zbir = 0;

                for(int j = 0; j < numberOfPossibleStatesZ; j++)
                    zbir += alpha[k-1][j] * transitionProbs[j][i];
            alpha[k][i] =  zbir * emissionProbs[X[k]][i];
            sumC += alpha[k][i];
        }
        C[k] = 1.0 / sumC;

        for(int i = 0; i < numberOfPossibleStatesZ; i++)
            alpha[k][i] /= sumC;
    }
}
void Model::computeAlphaForPredict()
{
    eraseAlpha();
    eraseC();
    initializeAlpha();
    initializeC();


    double sumC = 0;
    for (int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        alpha[0][i] = P[i] * emissionProbs[X[0]][i];
        sumC += alpha[0][i];
    }

    if(sumC != 0) C[0] = 1.0 / sumC;
    else C[0] = 0;

    for (int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        alpha[0][i] /= sumC;
    }

    for(int k = 1; k < numberOfObservedVars+1; k++)
    {
        sumC = 0;
        for(int i = 0; i < numberOfPossibleStatesZ; i++)
        {
            double zbir = 0;

                for(int j = 0; j < numberOfPossibleStatesZ; j++)
                {
                    zbir += alpha[k-1][j] * transitionProbs[j][i];
                }
            alpha[k][i] = zbir * emissionProbs[X[k]][i];
            sumC += alpha[k][i];
        }
        C[k] = 1.0 / sumC;
        for(int i = 0; i < numberOfPossibleStatesZ; i++)
        {
            if(sumC == 0) alpha[k][i] = 0;
            else alpha[k][i] /= sumC;
        }
    }
}
void Model::computeBeta()
{
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        beta[numberOfObservedVars-1][i] = 1.0;

    for(int k = numberOfObservedVars-2; k >= 0; k--)
    {
        for(int i = 0; i < numberOfPossibleStatesZ; i++)
            {
                beta[k][i] = 0;

                for(int j = 0; j < numberOfPossibleStatesZ; j++)
                    beta[k][i] += beta[k+1][j] * emissionProbs[X[k+1]][j] * transitionProbs[i][j];
                beta[k][i] *= C[k+1];
            }
    }
}

void Model::fit()
{
    double** newT = new double*[numberOfPossibleStatesZ];
    double** newE = new double*[numberOfPossibleStatesX];

    for (int i = 0; i<numberOfPossibleStatesX; i++)
        newE[i] = new double[numberOfPossibleStatesZ];

    for(int j = 0; j < numberOfPossibleStatesZ; j++)
        newT[j] = new double[numberOfPossibleStatesZ];

    for(int j = 0; j < numberOfPossibleStatesZ; j++)
        for(int l = 0; l < numberOfPossibleStatesZ; l++)
            newT[j][l] = 0;

    //making new Pi

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        P[i] = alpha[0][i] * beta[0][i];

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        double DD = 0.0;

        for(int j = 0; j < numberOfObservedVars-1; j++)
           DD += alpha[j][i] * beta[j][i];

        //computes new T

        for(int j = 0; j < numberOfPossibleStatesZ; j++)
        {
            double NN = 0.0;
            for(int t = 0; t < numberOfObservedVars-1; t++)
            {
                NN += C[t+1] * alpha[t][i] * emissionProbs[X[t+1]][j] * beta[t + 1][j];
            }
             newT[i][j] = transitionProbs[i][j] *(NN / DD);
        }

        for(int j = 0; j < numberOfPossibleStatesZ ; j++)
        {
            transitionProbs[i][j] = newT[i][j];
        }

        for(int j = 0; j < numberOfPossibleStatesX; j++)
            newE[j][i] = 0;

        DD += alpha[numberOfObservedVars - 1][i] * beta[numberOfObservedVars - 1][i];

        for(int t = 0; t < numberOfObservedVars; t++)
        {
            newE[X[t]][i] += alpha[t][i] * beta[t][i];
        }

        for(int j = 0; j < numberOfPossibleStatesX; j++)
            newE[j][i] /= DD;
    }
    for (int i = 0; i<numberOfPossibleStatesX; i++)
    {
    	for (int j = 0; j<numberOfPossibleStatesZ; j++) emissionProbs[i][j] = newE[i][j];
    }
    //cleaning up
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        delete[] newT[i];
    }
    delete[] newT;

    for(int i = 0; i < numberOfPossibleStatesX; i++)
    {
        delete[] newE[i];
    }
    delete[] newE;
}
void Model::getRandom()
{
    srand(time(NULL));
    initializeRandomSeed();

    for(int i = 0; i < 3; i++)
    {
        double random = ((double)rand() / RAND_MAX);
        randomSeed[i] = random;
    }
}


//PUBLIC FUNCTIONS

void Model::printNumberOfVars(){cout <<numberOfObservedVars << endl;}

void Model::printNumberOfPossibleStatesZ(){cout << numberOfPossibleStatesZ << endl;}

void Model::printNumberOfPossibleStatesX(){cout << numberOfPossibleStatesX << endl;}

void Model::printX()
{
    for(int i = 0; i < numberOfObservedVars; i++)
        cout << X[i] << " ";
    cout << endl;
}
void Model::printP()
{
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        cout << P[i] << " ";
    cout << endl;
}
void Model::printTrans()
{
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
    {
         for(int j = 0; j < numberOfPossibleStatesZ; j++)
            cout << transitionProbs[i][j] << " ";
        cout<<endl;
    }
}
void Model::printEmission()
{
    for(int i = 0 ; i < numberOfPossibleStatesX; i++)
    {
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            cout << emissionProbs[i][j] << " ";
        cout << endl;
    }
}
void Model::printNumberOfIterations(){cout << numberOfIterations << endl;}
void Model::printAlpha()
{
    for(int i = 0; i < numberOfObservedVars; i++)
    {
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            cout << alpha[i][j] << " ";
        cout << endl;
    }
}
void Model::printBeta()
{
    for(int i = 0; i < numberOfObservedVars; i++)
    {
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            cout << beta[i][j];
        cout << endl;
    }
}
void Model::printRandomSeed()
{
    for(int i = 0; i < 3; i++)
        cout << randomSeed[i] << " ";
}
void Model::printC()
{
    for(int i = 0; i < numberOfObservedVars; i++)
        cout << C[i] << " ";
}

int Model::getNumberOfObservedVars(){return numberOfObservedVars;}
void Model::setNumberOfObservedVars(int m_setNumberOfObservedVars){numberOfObservedVars = m_setNumberOfObservedVars;}


int Model::getNumberOfPossibleStatesZ(){return numberOfPossibleStatesZ;}
void Model::setNumberOfPossibleStatesZ(int m_setNumberOfPossibleStatesZ){numberOfPossibleStatesZ =m_setNumberOfPossibleStatesZ;}


int Model::getNumberOfPossibleStatesX(){return numberOfPossibleStatesX;}
void Model::setNumberOfPossibleStatesX(int m_setNumberOfPossibleStatesX){numberOfPossibleStatesX =m_setNumberOfPossibleStatesX;}


double *Model::getP(){return P;}
void Model::setP(double * m_P){P = m_P;}


double **Model::getTransitonMatrix(){return transitionProbs;}
void Model::setTransitionMatrix(double **m_transitionProbs){transitionProbs = m_transitionProbs;}


double **Model::getEmissionMatrix(){return emissionProbs;}
void Model::setEmissionMatrix (double **m_emissionProbs){emissionProbs = m_emissionProbs;}


int *Model::getArrayX(){return X;}
void Model::setArrayX(int *m_X){X = m_X;}

int Model::getIndexForZ()
{
    printP();
    cout << "stampam Pi" << endl;

    int limit = getNumberOfPossibleStatesZ();
    double *niz = new double[limit];

    double suma = P[0];
    niz[0] = suma;
    int i = 1;
    while((i < numberOfPossibleStatesZ) && (suma < 1.0))
    {
        suma += P[i];
        niz[i] = suma;
        if(1 - suma <= 0.0000000001) break;
        i++;
    }
    //srand(time(NULL));

    double r = randomSeed[0];//((double)rand() / RAND_MAX); // generating radnom value between 0 and 1

    int index;

    cout << "radnom vrednost: " << r << endl;
    i = 0;
    bool foundIt = false;
    while((i < numberOfPossibleStatesZ) && (foundIt == false))
    {
        if(r - niz[i] <= 0)
        {
            index = i;
            foundIt = true;
        }
        i++;
    }
    delete[] niz;
    return index;
}
int Model::getXFromE(int &currentZ)//picks X from col which is determined by Z,
{
    int Z = currentZ;

    for(int i = 0; i < numberOfPossibleStatesX; i++)
        cout << emissionProbs[i][Z] << " ";
    cout << "Stampam za dato Z gore" << endl;

    int limit = getNumberOfPossibleStatesX();
    double *niz = new double[limit];

    double suma = emissionProbs[0][Z];
    niz[0] = suma;
    int i = 1;

    while((i < numberOfPossibleStatesX) && (suma < 1.0))
    {
        suma += emissionProbs[i][Z];
        niz[i] = suma;
        if(1 - suma <= 0.0000000001) break;
        i++;
    }
    for(int i = 0; i < numberOfPossibleStatesX; i++)
        cout << niz[i] << " ";
    cout << endl;
    cout << "Odstampao novi niz " << endl;

    //srand(time(NULL));

    double r = randomSeed[1];//((double)rand() / RAND_MAX); // generating radnom value between 0 and 1

    int index;

    cout << "radnom vrednost: " << r << endl;
    i = 0;
    bool foundIt = false;
    while((i < numberOfPossibleStatesX) && (foundIt == false))
    {
        if(r - niz[i] <= 0)
        {
            index = i;
            foundIt = true;
        }
        i++;
    }
    delete[] niz;
    cout << "trazeni indeks za X je:" << index << endl;
    return index;
}
int Model::getZFromT(int &currentZ)//picks Z from col which is determined by currentZ,
{
    int Z = currentZ;

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        cout << transitionProbs[Z][i] << " ";
    cout << "Stampam za dato Z gore" << endl;

    int limit = getNumberOfPossibleStatesZ();
    double *niz = new double[limit];

    double suma = transitionProbs[Z][0];
    niz[0] = suma;
    int i = 1;

    while((i < numberOfPossibleStatesZ) && (suma < 1.0))
    {
        suma += transitionProbs[Z][i];
        niz[i] = suma;
        if(1 - suma <= 0.0000000001) break;
        i++;
    }
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        cout << niz[i] << " ";
    cout << endl;
    cout << "Odstampao novi niz " << endl;

    //srand(time(NULL));
    double r = randomSeed[2];//((double)rand() / RAND_MAX); // generating radnom value between 0 and 1

    int index;

    cout << "radnom vrednost: " << r << endl;
    i = 0;
    bool foundIt = false;
    while((i < numberOfPossibleStatesZ) && (foundIt == false))
    {
        if(r - niz[i] <= 0)
        {
            index = i;
            foundIt = true;
        }
        i++;
    }
    delete[] niz;
    cout << "trazeni indeks za Z je:" << index << endl;
    return index;
}
void Model::testPi()
{

    initializeAlpha();
    cout << "initializeAlpha();" << endl;
    initializeBeta();
    cout << "initializeBeta();" << endl;
    initializeC();

    for(int i = 0; i < numberOfIterations; i++)
    {
        computeAlpha();

        computeBeta();

        fit();
    }
}
int Model::predict()
{

    double *niz = new double[numberOfObservedVars+1];

    for(int state = 0; state < numberOfPossibleStatesX; state++)
    {
        X[numberOfObservedVars] = state;
        computeAlphaForPredict();
        double suma = 0;
        for(int i = 0; i < numberOfPossibleStatesZ; i++)
        {
            for(int j = 0; j < numberOfPossibleStatesZ; j++)
            {
                suma += emissionProbs[state][i] * transitionProbs[i][j] * alpha[numberOfObservedVars][j];
            }
        }
        niz[state] = suma ;    }

    int maxElem = niz[0];
    int maxIndex = 0;
    for(int i = 1; i < numberOfPossibleStatesX; i++)
        {
            if(niz[i] > maxElem)
            {
                maxElem = niz[i];
                maxIndex = i;
            }
        }
    return maxIndex;

}
