#include "model.h"
#include "unittest.h"
#include <algorithm>
#include <iostream>

using namespace std;

void printMat(double** a, int n, int m)
{
	for (int i =0; i<n; i++)
	{
		for (int j=0; j<m; j++)
		{
			cout << a[i][j] << " ";
		}
		cout << endl;
	}
}
void Model::initializeC()
{
    C = new double[numberOfObservedVars];
}
void Model::initializeAlpha()
{
    alpha = new double*[numberOfObservedVars];
    for(int i = 0; i < numberOfObservedVars; i++)
        alpha[i] = new double[numberOfPossibleStatesZ];

    for(int i = 0; i < numberOfObservedVars; i++)
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
void Model::initializeGamma()
{
    gamma = new double *[numberOfObservedVars];
    for(int i = 0; i < numberOfObservedVars; i++)
        gamma[i] = new double [numberOfPossibleStatesZ];

    for(int i = 0; i < numberOfObservedVars; i++)
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            gamma[i][j] = 0;
}
void Model::initializeKsi()
{
    ksi = new double*[numberOfPossibleStatesZ];
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
            ksi[i] = new double[numberOfPossibleStatesZ];

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            ksi[i][j] = 0;
}
void Model::initializeMi()
{
    //cout << "usao u initializeMi" << endl;
    mi = new double*[numberOfObservedVars];

    for(int i = 0; i < numberOfObservedVars; i++)
        mi[i] = new double[numberOfPossibleStatesZ];

    for(int i = 0; i < numberOfObservedVars; i++)
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            mi[i][j] = 0;
}
void Model::computeAlpha()
{
    double sumC = 0;
    for (int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        alpha[0][i] = P[i] * emissionProbs[X[0]][i];
        sumC += alpha[0][i];
    }

    C[0] = 1.0 /sumC;

    for (int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        alpha[0][i] *= C[0];
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

        //for(int i = 0; i < numberOfPossibleStatesZ; i++)
            //alpha[k][i] *= C[k];
    }
}

void Model::computeBeta()
{
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        beta[numberOfObservedVars-1][i] = 1.0;

    for(int k = numberOfObservedVars-2; k >= 0; k--)
        for(int i = 0; i < numberOfPossibleStatesZ; i++)
            {
                beta[k][i] = 0;

                for(int j = 0; j < numberOfPossibleStatesZ; j++)
                    beta[k][i] += beta[k+1][j] * emissionProbs[X[k+1]][j] * transitionProbs[i][j];
                beta[k][i] *= C[k+1];
            }
}

void Model::computeGamma()
{
    for(int k = 0; k < numberOfObservedVars; k++)
    {
        double zbir = 0;

        for(int i = 0; i < numberOfPossibleStatesZ; i++)
        {
            zbir += alpha[k][i] * beta[k][i];
        }

        double normCoef = 1 / zbir;


        for(int i = 0; i < numberOfPossibleStatesZ; i++)
        {
            gamma[k][i] = alpha[k][i] * beta[k][i] * normCoef;
        }
    }
}
void Model::computeNextP()
{
    for(int i = 0; i < numberOfPossibleStatesZ; i++)

    P[i] = gamma[0][i];
}

void Model::computeCurrentT()
{
    double** newT = new double*[numberOfPossibleStatesZ];
    double** newE = new double*[numberOfPossibleStatesX];
    for (int i = 0; i<numberOfObservedVars; i++) newE[i] = new double[numberOfPossibleStatesZ];

        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            newT[j] = new double[numberOfPossibleStatesZ];

        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            for(int l = 0; l < numberOfPossibleStatesZ; l++)
                newT[j][l] = 0;

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        P[i] = alpha[0][i] * beta[0][i];

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
    {
        cout << "usao u  i petlju" << endl;
        double DD = 0.0;

        for(int j = 0; j < numberOfObservedVars-1; j++)
           DD += alpha[j][i] * beta[j][i];


        cout << "odradio DD prvi" << endl;
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
        {
            double NN = 0.0;
            for(int t = 0; t < numberOfObservedVars-1; t++)
            {
                cout << "For t = " << t << " " << emissionProbs[X[t+1]][j] << endl;
                NN += C[t+1] * alpha[t][i] * emissionProbs[X[t+1]][j] * beta[t + 1][j];
            }
             newT[i][j] = transitionProbs[i][j] *(NN / DD);
             cout << NN << endl;
        }

        cout << " odradio newT" << endl;



        for(int j = 0; j < numberOfPossibleStatesZ ; j++)
        {
            transitionProbs[i][j] = newT[i][j];
        }

        for(int j = 0; j < numberOfPossibleStatesX; j++)
            newE[j][i] = 0;

        DD += alpha[numberOfObservedVars - 1][i] * beta[numberOfObservedVars - 1][i];

        cout << "novo DD" << endl;

        for(int t = 0; t < numberOfObservedVars; t++)
        {
            newE[X[t]][i] += alpha[t][i] * beta[t][i];
        }

        cout << "emissionProbs[X[t]][i] += alpha[i][t] * beta[i][t];" << endl;

        for(int j = 0; j < numberOfPossibleStatesX; j++)
            newE[j][i] /= DD;
        cout << "kraj petlje" << endl;
        cout << DD << endl;
    }
    for (int i = 0; i<numberOfPossibleStatesX; i++) {
    	for (int j = 0; j<numberOfPossibleStatesZ; j++) emissionProbs[i][j] = newE[i][j];
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
void *Model::setArrayX(int *m_X){X = m_X;}


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
        cout <<"computeAlpha" <<endl;
        //printAlpha();

        //computeBeta();
        cout << "computeBeta()" << endl;
        //printBeta();

      //  computeCurrentT();
        cout << "computeCurrentT" << endl;

    //    printP();
        cout << "printP" << endl;
        cout << endl;

  //      printC();
        cout << " printC" << endl;
        cout << endl;
//
 //       printTrans();
        cout << "printTrans" << endl;
        cout << endl;

//        printEmission();

        cout << "printEmission" << endl;

    }
}
