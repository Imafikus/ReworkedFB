#include "unittest.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void UnitTest::getExpectedT(int StatesZ)
{
    ifstream inf("expectedT.txt");
    expectedT = new double*[StatesZ];

    for(int i = 0; i < StatesZ; i++)
        for(int j = 0; j < StatesZ; j++)
            {
                expectedT[i] = new double[StatesZ];
                inf >> expectedT[i][j];
            }
}

void UnitTest::getExpectedE(int StatesX, int StatesZ)
{
    ifstream inf("expectedE.txt");
    expectedE = new double*[StatesX];

    for(int i = 0; i < StatesX; i++)
        for(int j = 0; j < StatesZ; j++)
            {
                expectedE[i] = new double[StatesZ];
                inf >> expectedE[i][j];
            }
}

void UnitTest::getExpectedP(int StatesZ)
{
    ifstream inf("pi.txt");
    expectedP = new double[StatesZ];

    for(int i = 0; i < StatesZ; i++)
        inf >> expectedP[i];
}
void UnitTest::printExpectedP(int statesZ)
{
    getExpectedP(statesZ);
    for(int i = 0; i < statesZ; i++)
        cout << expectedP[i] << " ";
}
void UnitTest::testT(double **computedT, int numberOfPossibleStatesZ)
{
    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            if(abs(computedT[i][j] - expectedT[i][j]) > k)
            {
                cout << "Greska pri racnunanju T[ " << i << " " << j << " ]" << endl;
                return;
            }
    cout << "Uspesno zavrsen test!" <<endl;
}

void UnitTest::testE(double **computedE, int numberOfPossibleStatesX, int numberOfPossibleStatesZ)
{
    for(int i = 0; i < numberOfPossibleStatesX; i++)
        for(int j = 0; j < numberOfPossibleStatesZ; j++)
            if(abs(computedE[i][j] - expectedE[i][j]) > k)
            {
                cout << "Greska pri racnunanju T[ " << i << " " << j << " ]" << endl;
                return;
            }
    cout << "Uspesno zavrsen test!" <<endl;
}

void UnitTest::testP(double *computedP, int numberOfPossibleStatesZ)
{
    getExpectedP(numberOfPossibleStatesZ);

    for(int i = 0; i < numberOfPossibleStatesZ; i++)
        if(abs(computedP[i] - expectedP[i]) > k)
        {
            cout << "Greska pri izracunavanju P[ " << i << " ]" <<endl;
            return;
        }
    cout << "Test uspesno zavrsen!" << endl;
}
