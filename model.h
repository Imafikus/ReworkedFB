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

    double **ksi;
    double *gamma;

    //double probs[1000][1000];
    double **normalizedProbs;


    void initializeT();
    void initializeE();

    void initializeP();

    void initializeAlpha();

    void initializeBeta();

    void initializeGamma();

    void initializeKsi();

    void computeAlpha();

    void computeBeta();

    void computeGamma(int currentObservedState);

    void computeNextP(int currentObservedState);

    void computeKsi(int currentObservedState);

    void computeCurrentT(int currentObservedState);

    void computeNormalized();


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

    void printGamma();

    void printNormalizedMatrix();

    void train();
};


#endif // MODEL_H

