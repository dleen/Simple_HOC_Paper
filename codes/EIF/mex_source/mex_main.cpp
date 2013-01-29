#include "mex.h"
#include "LIF_spike.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// To work with just a single LIF_spike unit:
	// Create the LIF_spike container

	double gamma, lambda, sigma;
	int N = mxGetScalar(prhs[3]);
	string neuron_model = "EIF";

	gamma = mxGetScalar(prhs[0]);
	lambda= mxGetScalar(prhs[1]);
	sigma = mxGetScalar(prhs[2]);

	LIF_spike Y(N);
	Y.create_XIF_data(gamma,lambda,sigma,neuron_model);

	plhs[0] = mxCreateDoubleMatrix(4,1,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(N+1,1,mxREAL);
    
    double *data1 = mxGetPr(plhs[0]);
    double *data2 = mxGetPr(plhs[1]);
    
    vector<double> P_ret(N+1,0);
    
    Y.return_statistics(data1[0],data1[1],data1[2],data1[3],P_ret);
    
    for(int i=0; i<N+1; ++i)
    {
        data2[i] = P_ret.at(i);
    }
}
