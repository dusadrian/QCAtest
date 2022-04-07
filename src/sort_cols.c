#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "sort_cols.h"
void sort_cols(int *p_matrix, int *sortcols, int *p_ck, const int nconds, unsigned int foundPI) {
    int temp;
    for (int i = nconds - 1; i >= 0; i--) {
        for (int c1 = 0; c1 < foundPI; c1++) {
            for (int c2 = c1 + 1; c2 < foundPI; c2++) {
                if (p_matrix[sortcols[c1] * nconds + i] < p_matrix[sortcols[c2] * nconds + i]) {
                    temp = sortcols[c2];
                    for (int c3 = c2; c3 > c1; c3--) {
                        sortcols[c3] = sortcols[c3 - 1];
                    }
                    sortcols[c1] = temp;
                }
            }
        }
        bool nonzero = true;
        int zeroidx = 0;
        while(zeroidx < foundPI && nonzero) {
            nonzero = p_matrix[sortcols[zeroidx] * nconds + i];
            zeroidx++;
        }
        zeroidx--;
        for (int c1 = 0; c1 < zeroidx; c1++) {
            for (int c2 = c1 + 1; c2 < zeroidx; c2++) {
                if (p_matrix[sortcols[c1] * nconds + i] > p_matrix[sortcols[c2] * nconds + i]) {
                    temp = sortcols[c2];
                    for (int c3 = c2; c3 > c1; c3--) {
                        sortcols[c3] = sortcols[c3 - 1];
                    }
                    sortcols[c1] = temp;
                }
            }
        }
    }
    for (int c1 = 0; c1 < foundPI; c1++) {
        for (int c2 = c1 + 1; c2 < foundPI; c2++) {
            if (p_ck[sortcols[c1]] > p_ck[sortcols[c2]]) {
                temp = sortcols[c2];
                for (int c3 = c2; c3 > c1; c3--) {
                    sortcols[c3] = sortcols[c3 - 1];
                }
                sortcols[c1] = temp;
            }
        }
    }
}
