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

    // Compute layer thicknesses and 1/rho (one division per layer instead of per lambda)
    double *h = mxMalloc((n > 1 ? n - 1 : 1) * sizeof(double));
    double *invRho = mxMalloc(n * sizeof(double));
    for (mwSize i = 0; i + 1 < n; ++i) {
        h[i] = depths[i + 1] - depths[i];
    }
    for (mwSize i = 0; i < n; ++i) {
        invRho[i] = 1.0 / rhos[i];
    }

    // For x >= TANH_SAT, tanh(x) rounds to exactly 1.0 (1 - tanh(x) < 1e-17),
    // and then Constable A-4 gives T_k = rho_k whatever T_{k+1} is.
    const double TANH_SAT = 20.0;

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
        double sum = 0.0;

        for (int j = 0; j < 11; ++j) {
            const double lam = lambda[j + 11 * i];

            // Start at the shallowest saturated layer (T = rho_k there) if
            // there is one; otherwise at the half-space (T = rho_n).
            mwSize start = n - 1;
            for (mwSize k = 0; k + 1 < n; ++k) {
                if (lam * h[k] >= TANH_SAT) { start = k; break; }
            }
            double T = rhos[start];

            // Constable A-4 for layers start-1 downto 0
            // (unsigned counter kk = k+1 avoids converting the unsigned n-2 to int)
            for (mwSize kk = start; kk > 0; --kk) {
                mwSize k = kk - 1;
                double tanh_val = tanh(lam * h[k]);
                T = (T + rhos[k] * tanh_val) / (1.0 + T * tanh_val * invRho[k]);
            }

            // Apparent resistivity is the filter-weighted sum of T
            sum += T * filter[j];
        }
        apparentResistivity[i] = sum;
    }

    mxFree(h);
    mxFree(invRho);
    mxFree(depths);
    mxFree(rhos);
}
