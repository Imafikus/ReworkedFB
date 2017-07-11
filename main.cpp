#include <iostream>
#include "model.h"
#include <vector>
#include <algorithm>

using namespace std;

void setValues(vector<int> &values, Model &model)
{
    for(int i = 1; i <= model.getNumberOfObservedVars(); i++)
    {
        int a;
        cin >> a;
        values.push_back(a);
    }
}
int main()
{
    Model model(4, 2, 4);
    //Model *model1= new Model(4, 2, 4);
    vector<int> values;
    setValues(values, model);
    model.setArrayX(values);
    model.train();
    return 0;
}
