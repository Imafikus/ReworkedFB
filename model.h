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

    int *helpX;
    int *helpZ;

    double **alpha; //forward part
    double **beta; // backward part

    double **ksi;
    double **gamma;
    double **mi;


    double **normalizedProbs;


    void initializeT();
    void initializeE();

    void initializeP();

    void initializeAlpha();

    void initializeBeta();

    void initializeGamma();

    void initializeKsi();

    void initializeMi();

    void computeAlpha();

    void computeBeta();

    void computeGamma();

    void computeNextP();

    void computeKsi(int currentObservedState);

    void computeCurrentT(int currentObservedState);

    void computeMi();

    void computeE();



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

