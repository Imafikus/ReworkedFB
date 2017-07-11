#ifndef MODEL_H
#define MODEL_H

#include <vector>

class Model
{
private:

    int numberOfObservedVars;
    int numberOfPossibleStatesZ;
    int numberOfPossibleStatesX;

    double **emissionProbs; // Emission probs
    double **transitionProbs; // Transitiom probs
    double *P; //probs for Z1 to take some of the possible states
    int *X; // array of inputs

    double **alpha; //forward part
    double **beta; // backward part

    //double probs[1000][1000];
    double **normalizedProbs;


    void initializeT();
    void initializeE();
    void initializeP();
    void initializeAlpha();
    void initializeBeta();
    void computeAlpha();
    void computeBeta();
    void computeNormalized();
        //fillArrayForGamma
        //computeKsi
        //computeCurrentT
        //computeCurrentP
        //computeCurrentE

public:

    Model(int numberOfObservedVars, int numberOfPossibleStatesX, int numberOfPossibleStatesZ):
        numberOfObservedVars(numberOfObservedVars),
        numberOfPossibleStatesX(numberOfPossibleStatesX),
        numberOfPossibleStatesZ(numberOfPossibleStatesZ) {}


    int getNumberOfObservedVars(){return numberOfObservedVars;}

    int getNumberOfPossibleStatesZ(){return numberOfPossibleStatesZ;}

    int getNumberOfPossibleStatesX(){return numberOfPossibleStatesX;}

    void setArrayX(std::vector<int> &values);

    void printNormalizedMatrix();
    void train();
};


#endif // MODEL_H

