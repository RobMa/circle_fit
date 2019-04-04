/*
 * circle_fit.c - Combined circle fitting
 * First an algebraic fit is computed based on the GWAF method and the approximation by Taubin.
 * The result is used as initial value for the LM based least-squares fitting algorithm.
 *
 * Estimates a circle from sample points.
 * Inputs:
 *      x - Nx1 vector, where N >= 3
 *      y - Nx1 vector
 * 
 * Output:
 *      rxy - 3x1 vector where rxy(1), rxy(2) and rxy(3) are the radius, x-position and y-position of the circle, respectively.
 * 
 *
 * The calling syntax is:
 *
 *		rxy = circle_fit(x, y)
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include <circle_fit/c_interface.h>
 
/**
 * @brief The gateway function
 * 
 * @param nlhs Number of output (left-side) arguments, or the size of the plhs array.
 * @param plhs Array of output arguments.
 * @param nrhs Number of input (right-side) arguments, or the size of the prhs array.
 * @param prhs Array of input arguments.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
    {
        mexErrMsgIdAndTxt("circle_fit:nrhs", "Two inputs required.");
    }
    if (nlhs != 1)
    {
        mexErrMsgIdAndTxt("circle_fit:nlhs", "One output required.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    {
        mexErrMsgIdAndTxt("circle_fit:notDouble", "Input x must be type double.");
    }
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
    {
        mexErrMsgIdAndTxt("circle_fit:notDouble", "Input y must be type double.");
    }
    if(mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("circle_fit:notRowVector", "x must be column vector.");
    }
    if(mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("circle_fit:notRowVector", "y must be column vector.");
    }
    if(mxGetM(prhs[0]) < 3) {
        mexErrMsgIdAndTxt("circle_fit:invalidVector", "x must of length 3 or more.");
    }
    if(mxGetM(prhs[0]) != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("circle_fit:invalidVectors", "x and y must be of the same length.");
    }

    /* Read inputs */
    int N = mxGetM(prhs[0]);
    double *x = mxGetDoubles(prhs[0]);
    double *y = mxGetDoubles(prhs[1]);

    /* create the output vector */
    plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
    double *rxy = mxGetDoubles(plhs[0]);

    estimate_circle(x, y, N, rxy+1, rxy+2, rxy);
}

