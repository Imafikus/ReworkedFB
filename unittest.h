#ifndef UNITTEST_H
#define UNITTEST_H

class UnitTest
{
private:
    const double k = 1/1000000;


    double **expectedT;
    double **expectedE;
    double *expectedP;

    void getExpectedT(int StatesZ);
    void getExpectedE(int StatesX, int StatesZ);
    void getExpectedP(int StatesZ);

public:

   void testT(double **computedT, int numberOfPossibleStatesZ);
   void testE(double **computedE, int numberOfPossibleStatesX, int numberOfPossibleStatesZ);
   void testP(double *computedP, int numberOfPossibleStatesZ);

};


#endif // UNITTEST_H
