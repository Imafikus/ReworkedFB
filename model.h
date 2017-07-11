class Model
{
private:

    double emissionProbs[1000][1000]; // Emission probs
    double transitionProbs[1000][1000]; // Transitiom probs
    double P[1000] //probs for Z1 to take some of the possible states
    int X[1000]; // array of inputs

    double alpha[1000][1000]; //forward part
    double beta[1000][1000]; // backward part

    //double probs[1000][1000];
    double normalizedProbs[1000][1000];

    double gama[1000];
    double ksi[1000];

    void initializeT()
        {
            for(int i = 1; i <= possibleStatesZ; i++)
                for(int j = 1; j <= possibleStatesZ; j++)
                    transitionProbs[i][j] = 1 / possibleStatesZ;
        }

    void initializeE()
        {
            for(int i = 1; i <= possibleStatesX; i++)
                for(int j = 1; j <= possibleStatesZ; j++)
                    emissionProbs][i]][j] = 1 / possibleStatesX;
        }

    void initializeP()
        {
            for(int i = 0; i <= possibleStatesZ; i++)
                P[i] = 1 / possibleStatesZ;
        }
    void initializeAlpha()
        {
            for(int i = 1; i <= possibleStatesZ; i++)
                for(int j = 1; j <= possibleStatesZ; i++)
                    alpha[i][j] = 0;
        }
    void initializeBeta()
        {
            for(int i = 1; i <= possibleStatesZ; i++)
                for(int j = 1; j <= possibleStatesZ; i++)
                    beta[i][j] = 0;
        }

    void computeAlpha()
        {
            for(int i = 1; i <= possibleStates; i++)
                alpfa[1][i] = P[i]*E[X[1]][i];

            for(int k = 2; k <= observedVars; k++)
                for(int i = 1; i <= possibleStates; i++)
                    for(int j = 1; j <= possibleStates; j++)
                        alpha[k][i] += alpha[k-1][i] * E[X[k]][j] * T[i][j];
        }
    void computeBeta()
        {
            for(int i = 1; i <= possibleStates; i++)
                beta[observedVars][i] = 1;

            for(int k = observedVars-1; k >= 1; k--)
                for(int i = 1; i <= possibleStates; i++)
                    for(int j = 1; j <= possibleStates; j++)
                        beta[k][i] += beta[k+1][i] * E[X[k+1]][j] * T[i][j];
        }
    void computeNormalized()
        {
            for (int k = 1; k <= observedVars; k++)
                {
                    float zbir = 0;

                    for(int i = 1; i <= possibleStates; i++)
                        zbir += alpha[k][i] * beta[k][i];

                    float x = 1 / zbir;

                    for (int i = 1; i <= possibleStates; i++)
                        normalized[k][i] = alpha[k][i] * beta[k][i] * x;
                }
        }
        //fillArrayForGamma
        //computeKsi
        //computeCurrentT
        //computeCurrentP
        //computeCurrentE

public:

    int numberOfObservedVars;
    int possibleStatesZ;
    int possibleStatesX;

    void getNumberOfObservedVars()
        {
            int a;
            cin << a;
            numberOfObservedVars = a;
        }
    void getNumberOfpossibleStatesZ()
        {
            int a;
            cin << a;
            possibleStatesZ = a;
        }
    void getNumberOfpossibleStatesX()
        {
            int a;
            cin << a;
            possibleStatesX = a;
        }
    void getArrayX()
        {
            for(int i = 1; i <= numberOfObservedVarsl; i++)
            {
                int a;
                cin >> a;
                X[i] = a;
            }
        }
    void printNormalizedMatrix()
        {
            for (int i = 1; i <= observedVars; i++)
                {
                    for (int j = 1; j <= possibleStates; j++)
                        cout << normalizedProbs[i][j];
                    cout << endl;
                }
        }
    void train()
        {
            getNumberOfpossibleStatesZ()
            getNumberOfObservedVars();
            getNumberOfpossibleStatesX();
            getArrayX();

            initializeT();
            initializeE();
            initializeP();
            initializeAlpha();
            initializeBeta();

            computeAlpha();
            computeBeta();
            computeNormalized();

            printNormalizedMatrix();

        }
};
