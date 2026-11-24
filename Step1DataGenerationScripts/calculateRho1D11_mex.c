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
 *   depths - [n x 1] vector of interface depths. May be padded with NaN
 *            (e.g. [0; 1; 25; NaN; ... NaN]), as in the rjMCMC ensembles.
 *   rhos   - [n x 1] vector of layer resistivities, NaN-padded like depths
 *   lambda - [11 x m] matrix
 *
 * As in calculateRho1D11.m, NaN entries are removed from depths and rhos
 * before the calculation, so a padded model gives the same result as the
 * unpadded one.
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
    for (int a = 0; a < 3; ++a) {
        if (!mxIsDouble(prhs[a]) || mxIsComplex(prhs[a]) || mxIsSparse(prhs[a])) {
            mexErrMsgTxt("Inputs must be real, full double arrays.");
        }
    }

    // Get input pointers
    const double *depthsIn = mxGetPr(prhs[0]);
    const double *rhosIn = mxGetPr(prhs[1]);
    const double *lambda = mxGetPr(prhs[2]);

    mwSize nDepthsIn = mxGetNumberOfElements(prhs[0]);
    mwSize nRhosIn = mxGetNumberOfElements(prhs[1]);
    mwSize m = mxGetN(prhs[2]);  // number of electrode spacings

    // Check lambda is 11 x m
    if (mxGetM(prhs[2]) != 11) {
        mexErrMsgTxt("lambda must have 11 rows.");
    }

    // Drop NaN padding, as calculateRho1D11.m does with depths(~isnan(depths))
    double *depths = mxMalloc((nDepthsIn > 0 ? nDepthsIn : 1) * sizeof(double));
    double *rhos = mxMalloc((nRhosIn > 0 ? nRhosIn : 1) * sizeof(double));
    mwSize n = 0;
    for (mwSize i = 0; i < nDepthsIn; ++i) {
        if (!mxIsNaN(depthsIn[i])) depths[n++] = depthsIn[i];
    }
    mwSize nRhos = 0;
    for (mwSize i = 0; i < nRhosIn; ++i) {
        if (!mxIsNaN(rhosIn[i])) rhos[nRhos++] = rhosIn[i];
    }
    if (n == 0) {
        mxFree(depths); mxFree(rhos);
        mexErrMsgTxt("depths must contain at least one non-NaN value.");
    }
    if (nRhos < n) {
        mxFree(depths); mxFree(rhos);
        mexErrMsgTxt("rhos has fewer non-NaN values than depths.");
    }

    // Compute layer thicknesses
    double *h = mxMalloc((n > 1 ? n - 1 : 1) * sizeof(double));
    for (mwSize i = 0; i + 1 < n; ++i) {
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
        // (unsigned counter kk = k+1 avoids converting the unsigned n-2 to int)
        for (mwSize kk = n - 1; kk > 0; --kk) {
            mwSize k = kk - 1;
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
    mxFree(depths);
    mxFree(rhos);
}
