#include "mex.h"
#include <math.h>
/* NOTE - THIS CODE WAS TRANSLATED FROM MATLAB BY CHATGPT */
/*
 * calculateRho1D11_mex.c
 *
 * This MEX function replicates the MATLAB function 'calculateRho1D11'.
 * Compile in MATLAB using: mex calculateRho1D11_mex.c
 *
 * Inputs:
 *   depths - [n x 1] vector
 *   rhos   - [n x 1] vector
 *   lambda - [11 x m] matrix
 *
 * Output:
 *   apparentResistivity - [1 x m] vector
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Input validation
    if (nrhs != 3) {
        mexErrMsgTxt("Three inputs required: depths, rhos, lambda.");
    }
    if (nlhs > 1) {
        mexErrMsgTxt("Only one output allowed.");
    }

    // Get input pointers
    double *depths = mxGetPr(prhs[0]);
    double *rhos = mxGetPr(prhs[1]);
    double *lambda = mxGetPr(prhs[2]);

    mwSize n = mxGetM(prhs[0]);  // number of layers (depths and rhos are column vectors)
    mwSize m = mxGetN(prhs[2]);  // number of electrode spacings

    // Check lambda is 11 x m
    if (mxGetM(prhs[2]) != 11) {
        mexErrMsgTxt("lambda must have 11 rows.");
    }

    // Compute layer thicknesses
    double *h = mxMalloc((n - 1) * sizeof(double));
    for (mwSize i = 0; i < n - 1; ++i) {
        h[i] = depths[i + 1] - depths[i];
    }

    // Filter coefficients (Guptasarma 1982)
    double filter[11] = {
        0.041873,
        -0.022258,
        0.387660,
        0.647103,
        1.84873,
        -2.96084,
        1.358412,
        -0.377590,
        0.097107,
        -0.024243,
        0.004046
    };

    // Create output
    plhs[0] = mxCreateDoubleMatrix(1, m, mxREAL);
    double *apparentResistivity = mxGetPr(plhs[0]);

    // Main loop over each electrode spacing
    for (mwSize i = 0; i < m; ++i) {
        double T[11];
        double lam[11];

        // Load lambda column
        for (int j = 0; j < 11; ++j) {
            lam[j] = lambda[j + 11 * i];
            T[j] = rhos[n - 1];  // Start with rho_n
        }

        // Loop over layers from n-2 downto 0
        for (int k = n - 2; k >= 0; --k) {
            for (int j = 0; j < 11; ++j) {
                double tanh_val = tanh(lam[j] * h[k]);
                double numerator = T[j] + rhos[k] * tanh_val;
                double denominator = 1 + (T[j] * tanh_val / rhos[k]);
                T[j] = numerator / denominator;
            }
        }

        // Compute apparent resistivity by summing filter * T
        double sum = 0.0;
        for (int j = 0; j < 11; ++j) {
            sum += T[j] * filter[j];
        }
        apparentResistivity[i] = sum;
    }

    mxFree(h);
}