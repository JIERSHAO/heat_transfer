#include "lagrangian.hpp"
#define N_k 30
#define N_e 30
int main()
{
    double* x = new double[4];
    double* y = new double[4];
    
    x[0]=0.0; x[1]=1.0; x[2]=1.0; x[3]=1.0;
    y[0]=0.0; y[1]=0.0; y[2]=4.0; y[3]=4.0;
   
    transform(x, y, N_k, N_e);

    delete [] x;
    delete [] y;
    return 0;
}