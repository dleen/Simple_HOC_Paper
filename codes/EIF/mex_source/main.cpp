//#include "mex.h"
#include "LIF_spike.h"

using namespace std;

int main()
{
	// To work with just a single LIF_spike unit:
	// Create the LIF_spike container

	double gamma, lambda, sigma;
	int N = 5;
	string neuron_model = "EIF";

	gamma = -60;
	lambda= 0;
	sigma = 6.23;

	LIF_spike Y(N);
	Y.create_XIF_data(gamma,lambda,sigma,neuron_model);

    
    vector<double> P_ret(N+1,0);
    double data1[4];
    
    Y.return_statistics(data1[0],data1[1],data1[2],data1[3],P_ret);
    
    
    
    for(int i=0; i<N+1; ++i)
    {
        cout<< P_ret.at(i)<<endl;
    }
}
