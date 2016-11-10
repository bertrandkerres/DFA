#include "mex.h"
#include "matrix.h"
#include "DFA.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	/* Declare input variables in C++ data types */
    double *x;              // Signal data, double col vector, first input
    mwSize x_length;
    mwIndex *scale;         // Scale as number of samples, uint row vector, second input
    mwSize scale_length;
    double *q;              // Orders, double col vector, third input
    mwSize q_length;
    unsigned int m;         // Detrending order, uint scalar, fourth input
    
    /* Declare output variable (the structure function estimates) */
    void *F2;
    
    /* Check input arguments */
    if(nrhs != 4) {
        mexErrMsgIdAndTxt("Signals:DFA:nrhs", "Four inputs required.");
    }
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("Signals:DFA:nrhs", "One output required.");
    }
    // Check x vector
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("Signals:DFA:notDouble",
            "First parameter (x) must be of type double.");
    }
    if(mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("Signals:DFA:notColVector",
            "First parameter (x) must be a column vector.");
    }
    x = mxGetPr(prhs[0]);
    x_length = (mwSize) mxGetM(prhs[0]);
    
    // Check scale vector
    if( !mxIsUint32(prhs[1]) ) {
        mexErrMsgIdAndTxt("Signals:DFA:notUInt32",
            "Second parameter (scale) must be of type uint32");
    }
    if( mxGetM(prhs[1]) != 1 ) {
        mexErrMsgIdAndTxt("Signals:DFA:notColVector",
            "Second parameter (scale) must be a row vector.");
    }
    scale = (int *) mxGetData(prhs[1]);
    scale_length = (mwSize) mxGetN(prhs[1]);
    
    // Check q vector
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("Signals:DFA:notDouble",
            "Third parameter (q) must be of type double.");
    }
    if(mxGetN(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("Signals:DFA:notColVector",
            "Third parameter (q) must be a column vector.");
    }
    q = mxGetPr(prhs[2]);
    q_length = (mwSize) mxGetM(prhs[2]);
    
    // Check m scalar
    if( !mxIsUint32(prhs[3]) ) {
        mexErrMsgIdAndTxt("Signals:DFA:notUInt32",
            "Fourth parameter (m) must be of type uint32");
    }
    if( !(mxGetM(prhs[3])==1 && mxGetN(prhs[3])==1) ) {
        mexErrMsgIdAndTxt("Signals:DFA:notScalar",
            "Fourth parameter (m) must be a scalar array.");
    }
    m = (unsigned int) mxGetScalar(prhs[3]);

    try {
        F2 = mxMalloc(scale_length* q_length* sizeof(double));
        DFAm( x, x_length, (unsigned int*)scale, scale_length, q, q_length, m, F2 );
    }
    catch (std::exception &e) {
        mexPrintf (e.what());
        plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        mxFree (F2);
        return;
    }
       
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[0], (double*)F2);
    mxSetM(plhs[0], q_length);
    mxSetN(plhs[0], scale_length);
    
}