#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "find_consistent_models.h"
#include "utils.h"
#include "utils.h"
#include <R_ext/RS.h> 
void find_consistent_models(
    const int p_implicants[],
    const int indx[],
    const int ck[],
    const double p_data[],
    const int p_fuzzy[],
    const int nconds,
    const int nrdata,
    const int posrows,
    const double solcons,
    const double solcov,
    const bool allsol,
    const int soldepth,
    const unsigned int foundPI,
    const double maxcomb,
    int **solutions,
    int *nr,
    int *nc) {
    unsigned int estimsol = 1000;
    unsigned int maxk = posrows;
    if (foundPI < maxk) {
        maxk = foundPI;
    }
    if (soldepth < maxk && soldepth > 0) {
        maxk = soldepth;
    }
    int *p_sol = R_Calloc(maxk * estimsol, int);
    int *cksol = R_Calloc(estimsol, int);
    unsigned int solfound = 0;
    unsigned int prevfound = 0;
    bool keep_searching = true;
    int k = 1;
    double counter = 1;
    while (keep_searching && k <= maxk) {
            int tempk[k];
            for (int i = 0; i < k; i++) {
                tempk[i] = i; 
            }
            tempk[k - 1] -= 1; 
            int e = 0;
            int h = k;
            bool last = (foundPI == k);
            while (keep_searching && ((tempk[0] != foundPI - k) || last)) {
                increment(k, &e, &h, foundPI + last, tempk, 0);
                last = false;
                bool nonred = true;
                int i = 0;
                while (i < prevfound && nonred) {
                int sumeq = 0;
                int v = 0;
                while (sumeq == v && v < cksol[i]) {
                    for (int c = 0; c < k; c++) {
                        if (p_sol[i * maxk + v] == tempk[c]) {
                            sumeq++;
                        }
                    }
                    v++;
                }
                if (sumeq == v) { 
                    nonred = false; 
                }
                i++;
            }
            if (nonred) {
                if (consistent_solution(p_data, nconds, nrdata, k, tempk, foundPI, p_implicants, ck, indx, p_fuzzy, solcons, solcov)) {
                    for (int c = 0; c < k; c++) {
                        p_sol[solfound * maxk + c] = tempk[c];
                    }
                    cksol[solfound] = k;
                    solfound++;
                    if (solfound == estimsol) {
                        estimsol += 1000;
                        resize(&p_sol,  maxk,  estimsol, solfound);
                        resize(&cksol,   1,  estimsol, solfound);
                    }
                }
            }
            if (maxcomb > 0) {
                counter++;
                if ((counter / 1000000000) >= maxcomb) {
                    keep_searching = false;
                }
            }
        }
        prevfound = solfound;
        k += 1;
    }
    int *p_tempmat = R_Calloc(1, int);
    if (solfound > 0) {
        int finalrows = cksol[solfound - 1];
        R_Free(p_tempmat);
        p_tempmat = R_Calloc(finalrows * solfound, int);
        for (int c = 0; c < solfound; c++) {
            for (int r = 0; r < cksol[c]; r++) {
                p_tempmat[c * finalrows + r] = p_sol[c * maxk + r] + 1; 
            }
        }
        *nr = finalrows;
        *nc = solfound;
    }
    R_Free(p_sol);
    R_Free(cksol);
    R_Free(*solutions);
    *solutions = p_tempmat;
}
