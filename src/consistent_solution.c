#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "utils.h"
#include "consistent_solution.h"
#include <R_ext/RS.h> 
bool consistent_solution(
    const double p_data[],
    const int nconds,
    const int nrdata,
    const int k,
    const int tempk[],
    const unsigned int foundPI,
    const int p_implicants[],
    const int ck[],
    const int indx[],
    const int p_fsconds[],
    const double solcons,
    const double solcov) {
    int cindx, val;
    double *p_y = R_Calloc(1, double);
    double ymat[nrdata * k];
    double sumy = 0;
    for (int r = 0; r < nrdata; r++) {
        sumy += p_data[nconds * nrdata + r];
    }
    for (int i = 0; i < k; i++) { 
        int k2 = ck[tempk[i]];
        R_Free(p_y);
        p_y = R_Calloc(nrdata * k2, double);
        for (int c = 0; c < k2; c++) {
            cindx = indx[tempk[i] * nconds + c] - 1;
            val = p_implicants[tempk[i] * nconds + cindx] - 1;
            if (p_fsconds[cindx]) {
                bool negation = val == 0;
                for (int r = 0; r < nrdata; r++) {
                    p_y[c * nrdata + r] = negation ? (1 - p_data[cindx * nrdata + r]) : p_data[cindx * nrdata + r];
                }
            }
            else {
                for (int r = 0; r < nrdata; r++) {
                    p_y[c * nrdata + r] = (p_data[cindx * nrdata + r] == val) ? 1 : 0;
                }
            }
        }
        double pminx;
        for (int r = 0; r < nrdata; r++) {
            pminx = 1;
            for (int c = 0; c < k2; c++) {
                if (p_y[c * nrdata + r] < pminx) {
                    pminx = p_y[c * nrdata + r];
                }
            }
            ymat[i * nrdata + r] = pminx;
        }
    }
    double pmaxx;
    double sumx = 0, sumxy = 0;
    for (int r = 0; r < nrdata; r++) {
        pmaxx = 0;
        for (int c = 0; c < k; c++) {
            if (ymat[c * nrdata + r] > pmaxx) {
                pmaxx = ymat[c * nrdata + r];
            }
        }
        sumx += pmaxx;
        sumxy += ((pmaxx < p_data[nconds * nrdata + r]) ? pmaxx: p_data[nconds * nrdata + r]);
    }
    R_Free(p_y);
    return(agteb(sumxy / sumx, solcons) && agteb(sumxy / sumy, solcov));
}
