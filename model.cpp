#include "model.h"
#include "unittest.h"
#include <algorithm>
#include <iostream>

using namespace std;




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

    for(int i = 0; i < numberOfObservedVars; i++)
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            ksi[i][j] = 0;
}
void Model::initializeMi()
{
    mi = new double*[numberOfObservedVars];

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        mi[i] = new double[numberOfPossibleStatesZ];

}
void Model::computeAlpha()
{
    for (int i = 0; i < numberOfPossibleStatesZ; i++)
        alpha[0][i] = P[i] * emissionProbs[X[0]][i];

    for(int k = 1; k < numberOfObservedVars; k++)
    {
        for(int i = 1; i < numberOfPossibleStatesZ; i++)
        {
            alpha[k][i] = 0;

                for(int j = 0; j < numberOfPossibleStatesZ; j++)
                    alpha[k][i] += alpha[k-1][j] * emissionProbs[X[k]][i] * transitionProbs[j][i];
        }
    }
}

void Model::computeBeta()
{
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        beta[numberOfObservedVars-1][i] = 1;

     for(int i = 0; i < numberOfObservedVars; i++)
            {
                 for(int j = 0; j < numberOfPossibleStatesZ; j++)
                    cout << transitionProbs[i][j] << " ";
                cout<<endl;
            }
    cout << endl;
    cout <<"========================"<<endl;
    cout << endl;
     for(int i = 0; i < numberOfObservedVars; i++)
        {
             for(int j = 0; j < numberOfPossibleStatesZ; j++)
                cout << emissionProbs[i][j] << " ";
            cout<<endl;
        }

    for(int k = numberOfObservedVars-2; k >= 0; k--)
        for(int i = 0; i < numberOfPossibleStatesZ; i++)
            {
                beta[k][i] = 0;

                for(int j = 0; j < numberOfPossibleStatesZ; j++)
                    beta[k][i] += beta[k+1][j] * emissionProbs[X[k+1]][j] * transitionProbs[i][j];
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
           // cout << "usao u drugi for: zbir += alpha[k][i] * beta[k][i];";
        }

        double normCoef = 1 / zbir;
        //cout << "nasao koeficijent" << endl;
        for(int i = 0; i < numberOfPossibleStatesZ; i++)
        {
            gamma[k][i] = alpha[k][i] * beta[k][i] * normCoef;
          //  cout << "racuna gama zapravo";
        }
    }
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
void Model::computeNextP()
{
    for(int i = 0; i < numberOfPossibleStatesZ; i++)

    P[i] = gamma[0][i];
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
void Model::computeMi()
{
    for(int i = 0; i < numberOfPossibleStatesX; i++)
        for(int l = 0; l < numberOfPossibleStatesZ; l++)
        {
            double br = 0;
            double im = 0;
	    for (int n=0; n<numberOfObservedVars; n++)
	    {
	    	if (X[n] == i) br += gamma[n][l];
		im += gamma[n][l];
	    }

            mi[i][l] = br / im;
        }
}
void Model::computeE()
{
    for(int a = 0; a < numberOfPossibleStatesX; a++)
        for(int b = 0; b < numberOfPossibleStatesZ; b++)
            {
                 emissionProbs[a][b] = 1;

                 helpX = new int[numberOfPossibleStatesX];
                 helpZ = new int[numberOfPossibleStatesZ];

                for(int i = 0; i < numberOfPossibleStatesX; i++)
                    {
                        if(i == a) helpX[i] = 1;
                        else helpX[i] = 0;
                    }
                for(int i = 0; i < numberOfPossibleStatesZ; i++)
                    {
                        if(i == b) helpZ[i] = 1;
                        else helpZ[i] = 0;
                    }

                for(int i = 0; i < numberOfPossibleStatesX; i++)
                    for(int j = 0; j < numberOfPossibleStatesZ; j++)
                        {
                            if( (helpX[i] * helpZ[j]) == 1) emissionProbs[a][b] *= mi[i][j];
                            else emissionProbs[a][b] *= 1;
                        }
                delete[] helpX;
                delete[] helpZ;
            }
}


//PUBLIC FUNCTIONS


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

void Model::printX()
{
    for(int i = 0; i < numberOfObservedVars; i++)
        cout << X[i] << endl;
}

void Model::train()
{
    int currentObservedState = 0;//OVO OVDE SE PROSLEDJUJE u train()
    /*initializeT();
    initializeE();
    initializeP();*/

    initializeAlpha();
    initializeBeta();
    initializeKsi();


    //for(int currentObservedState = 0; currentObservedState < numberOfPossibleStatesX; currentObservedState++)
      //  {
            computeAlpha();
            computeBeta();
            computeGamma();
            computeKsi(currentObservedState);
       // }


}
void Model::testPi()
{
    initializeAlpha();
    cout << "initializeAlpha();" << endl;
    initializeBeta();
    cout << " initializeBeta();" << endl;
    initializeGamma();
    cout << "initializeGamma();" << endl;

    computeAlpha();

    cout <<"computeAlpha" <<endl;
    computeBeta();
    cout << "computeBeta()" << endl;
    computeGamma();
    cout << "computeGamma" << endl;
    computeNextP();
    cout << "computeNextP" <<endl;

}
