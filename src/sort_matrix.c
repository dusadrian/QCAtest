#include "sort_matrix.h"
void sort_matrix(int *p_matrix, int *p_colindx, int *p_ck, int nconds, int foundPI) {
    for (int i = 0; i < foundPI; i++) {
        p_colindx[i] = i;
    }
    int temp;
    for (int i = nconds - 1; i >= 0; i--) {
        for (int c1 = 0; c1 < foundPI; c1++) {
            for (int c2 = c1 + 1; c2 < foundPI; c2++) {
                if (p_matrix[p_colindx[c1] * nconds + i] < p_matrix[p_colindx[c2] * nconds + i]) {
                    temp = p_colindx[c2];
                    for (int c3 = c2; c3 > c1; c3--) {
                        p_colindx[c3] = p_colindx[c3 - 1];
                    }
                    p_colindx[c1] = temp;
                }
            }
        }
        bool nonzero = true;
        int zeroidx = 0;
        while(zeroidx < foundPI && nonzero) {
            nonzero = p_matrix[p_colindx[zeroidx] * nconds + i];
            zeroidx++;
        }
        zeroidx--;
        for (int c1 = 0; c1 < zeroidx; c1++) {
            for (int c2 = c1 + 1; c2 < zeroidx; c2++) {
                if (p_matrix[p_colindx[c1] * nconds + i] > p_matrix[p_colindx[c2] * nconds + i]) {
                    temp = p_colindx[c2];
                    for (int c3 = c2; c3 > c1; c3--) {
                        p_colindx[c3] = p_colindx[c3 - 1];
                    }
                    p_colindx[c1] = temp;
                }
            }
        }
    }
    for (int c1 = 0; c1 < foundPI; c1++) {
        for (int c2 = c1 + 1; c2 < foundPI; c2++) {
            if (p_ck[p_colindx[c1]] > p_ck[p_colindx[c2]]) {
                temp = p_colindx[c2];
                for (int c3 = c2; c3 > c1; c3--) {
                    p_colindx[c3] = p_colindx[c3 - 1];
                }
                p_colindx[c1] = temp;
            }
        }
    }
}
