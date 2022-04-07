#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "row_dominance.h"
void row_dominance(int p_pichart[], int p_implicants[], int *p_ck, int pirows, unsigned int *foundPI, int nconds) {
    unsigned int picols = *foundPI;
    bool survcols[picols];
    int colsums[picols];
    int sortcol[picols];
    int temp;
    for (int c = 0; c < picols; c++) {
        colsums[c] = 0;
        for (int r = 0; r < pirows; r++) {
            colsums[c] += p_pichart[c * pirows + r];
        }
        sortcol[c] = c;
        survcols[c] = true;
    }
    for (int c1 = 0; c1 < picols; c1++) {
        for (int c2 = c1 + 1; c2 < picols; c2++) {
            if (colsums[sortcol[c1]] < colsums[sortcol[c2]]) {
                temp = sortcol[c1];
                sortcol[c1] = sortcol[c2];
                sortcol[c2] = temp;
            }
        }
    }
    for (int c1 = 0; c1 < picols; c1++) {
        if (survcols[sortcol[c1]]) {
            for (int c2 = c1 + 1; c2 < picols; c2++) {
                if (survcols[sortcol[c2]]) {
                    if (colsums[sortcol[c1]] > colsums[sortcol[c2]]) {
                        bool itcovers = true; 
                        int r = 0;
                        while (r < pirows && itcovers) {
                            if (p_pichart[sortcol[c2] * pirows + r]) {
                                itcovers = p_pichart[sortcol[c1] * pirows + r];
                            }
                            r++;
                        }
                        if (itcovers) {
                            survcols[sortcol[c2]] = false;
                            --(*foundPI);
                        }
                    }
                }
            }
        }
    }
    if (*foundPI < picols) {
        int s = 0;
        for (int c = 0; c < picols; c++) {
            if (survcols[c]) {
                for (int r = 0; r < pirows; r++) {
                    p_pichart[s * pirows + r] = p_pichart[c * pirows + r];
                }
                for (int r = 0; r < nconds; r++) {
                    p_implicants[s * nconds + r] = p_implicants[c * nconds + r];
                }
                p_ck[s] = p_ck[c];
                s++;
            }
        }
    }
    return;
}
