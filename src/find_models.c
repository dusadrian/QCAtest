#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "find_models.h"
#include <R_ext/RS.h> 
void find_models(
    const int p_pichart[],
    const int pirows,
    const unsigned int picols,
    const bool allsol,
    const int k,
    const double maxcomb,
    const bool firstmin,
    int **solutions,
    int *nr,
    int *nc) {
    if (k == picols) {
        int* p_temp = R_Calloc(k, int);
        for (int i = 0; i < k; i++) {
            p_temp[i] = i + 1;
        }
        R_Free(*solutions);
        *solutions = p_temp;
        *nr = k;
        *nc = 1;
        return;
    }
    int *p_temp1 = R_Calloc(1, int);
    int *p_temp2 = R_Calloc(1, int);
    if (allsol) {
        int indmat[picols * pirows];
        int mintpis[pirows];
        for (int r = 0; r < pirows; r++) {
            mintpis[r] = 0;
            for (unsigned int c = 0; c < picols; c++) {
                if (p_pichart[c * pirows + r]) {
                    indmat[r * picols + mintpis[r]] = c;
                    mintpis[r]++;
                }
            }
        }
        R_Free(p_temp2);
        p_temp2 = R_Calloc(picols * mintpis[0], int);
        int *p_cols = R_Calloc(1, int); 
        for (int i = 0; i < mintpis[0]; i++) {
            p_temp2[i * picols + indmat[i]] = 1;
        }
        int tempcols = mintpis[0];
        for (int i = 1; i < pirows; i++) {
            R_Free(p_temp1);
            p_temp1 = R_Calloc(picols * tempcols * mintpis[i], int);
            for (int j = 0; j < mintpis[i]; j++) {
                Memcpy(&p_temp1[j * tempcols * picols], p_temp2, tempcols * picols);
                for (int tc = 0; tc < tempcols; tc++) {
                    p_temp1[(j * tempcols + tc) * picols + indmat[i * picols + j]] = 1;
                }
            }
            int temp2cols = tempcols * mintpis[i];
            R_Free(p_cols);
            p_cols = R_Calloc(temp2cols, int);
            for (int i = 0; i < temp2cols; i++) {
                p_cols[i] = true;
            }
            int survcols = temp2cols;
            super_rows(p_temp1, picols, &survcols, p_cols);
            R_Free(p_temp2);
            p_temp2 = R_Calloc(picols * survcols, int);
            Memcpy(p_temp2, p_temp1, picols * survcols);
            tempcols = survcols;
        }
        R_Free(p_temp1);
        p_temp1 = R_Calloc(picols * tempcols, int);
        R_Free(p_cols);
        p_cols = R_Calloc(tempcols, int);
        int maxr = 0;
        for (int c = 0; c < tempcols; c++) {
            for (int r = 0; r < picols; r++) {
                if (p_temp2[c * picols + r]) {
                    p_temp1[c * picols + p_cols[c]] = r + 1;
                    p_cols[c]++;
                }
                if (maxr < p_cols[c]) {
                    maxr = p_cols[c];
                }
            }
        }
        R_Free(p_temp2);
        p_temp2 = R_Calloc(maxr * tempcols, int);
        for (int c = 0; c < tempcols; c++) {
            for (int r = 0; r < maxr; r++) {
                p_temp2[c * maxr + r] = p_temp1[c * picols + r];
            }
        }
        int temp;
        int order[tempcols];
        for (int c = 0; c < tempcols; c++) {
            order[c] = c;
        }
        for (int r = maxr - 1; r >= 0; r--) {
            for (int c1 = 0; c1 < tempcols; c1++) {
                for (int c2 = c1 + 1; c2 < tempcols; c2++) {
                    if (p_temp2[order[c1] * maxr + r] > p_temp2[order[c2] * maxr + r]) {
                        temp = order[c2];
                        for (int i = c2; i > c1; i--) {
                            order[i] = order[i - 1];
                        }
                        order[c1] = temp;
                    }
                }
            }
        }
        for (int c1 = 0; c1 < tempcols; c1++) {
            for (int c2 = c1 + 1; c2 < tempcols; c2++) {
                if (p_cols[order[c1]] > p_cols[order[c2]]) {
                    temp = order[c2];
                    for (int i = c2; i > c1; i--) {
                        order[i] = order[i - 1];
                    }
                    order[c1] = temp;
                }
            }
        }
        R_Free(p_cols);
        R_Free(p_temp1);
        p_temp1 = R_Calloc(maxr * tempcols, int);
        for (int c = 0; c < tempcols; c++) {
            for (int r = 0; r < maxr; r++) {
                p_temp1[c * maxr + r] = p_temp2[order[c] * maxr + r];
            }
        }
        *nr = maxr;
        *nc = tempcols;
    }
    else {
        unsigned int solfound = 0;
        unsigned int estimsol = 100;
        R_Free(p_temp1);
        p_temp1 = R_Calloc(k * estimsol, int);
        int tempk[k];
        for (int i = 0; i < k; i++) {
            tempk[i] = i; 
        }
        tempk[k - 1] -= 1; 
        int e = 0;
        int h = k;
        bool keep_searching = true;
        bool last = (picols == k);
        double counter = 1;
        while (keep_searching && ((tempk[0] != picols - k) || last)) {
            increment(k, &e, &h, picols + last, tempk, 0);
            last = false;
            bool allrows = true;
            int r = 0;
            while (r < pirows && allrows) {
                bool covered = false;
                int c = 0;
                while (c < k && !covered) {
                    covered = p_pichart[tempk[c] * pirows + r];
                    c++;
                }
                allrows = covered;
                r++;
            }
            if (allrows) {
                for (int c = 0; c < k; c++) {
                    p_temp1[solfound * k + c] = tempk[c] + 1; 
                }
                solfound++;
                if (solfound == estimsol) {
                    estimsol += 100;
                    resize(&p_temp1, k, estimsol, solfound);
                }
            }
            if (firstmin && solfound > 0) {
                keep_searching = false;
            }
            if (maxcomb > 0) {
                counter++;
                if ((counter / 1000000000) >= maxcomb) {
                    keep_searching = false;
                }
            }
        }
        if (solfound > 0) {
            p_temp1 = Realloc(p_temp1, k * solfound, int);
        }
        else {
            R_Free(p_temp1);
            p_temp1 = R_Calloc(1, int);
        }
        *nr = k;
        *nc = solfound;
    }
    R_Free(p_temp2);
    R_Free(*solutions);
    *solutions = p_temp1;
}
