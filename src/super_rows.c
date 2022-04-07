#include "super_rows.h"
void super_rows(int *p_matrix, int rows, int *survcols, int *p_cols) {
    int cols = *survcols;
    int colsums[cols];
    int sortcol[cols];
    int temp;
    for (int c = 0; c < cols; c++) {
        colsums[c] = 0;
        for (int r = 0; r < rows; r++) {
            colsums[c] += p_matrix[c * rows + r];
        }
        sortcol[c] = c;
    }
    for (int c1 = 0; c1 < cols; c1++) {
        for (int c2 = c1 + 1; c2 < cols; c2++) {
            if (colsums[sortcol[c1]] > colsums[sortcol[c2]]) {
                temp = sortcol[c1];
                sortcol[c1] = sortcol[c2];
                sortcol[c2] = temp;
            }
        }
    }
    for (int c1 = 0; c1 < cols - 1; c1++) {
        if (p_cols[sortcol[c1]]) {
            for (int c2 = c1 + 1; c2 < cols; c2++) {
                if (p_cols[sortcol[c2]]) {
                    if (colsums[sortcol[c1]] <= colsums[sortcol[c2]]) {
                        bool cover = true; 
                        int r = 0;
                        while (r < rows && cover) {
                            if (p_matrix[sortcol[c1] * rows + r]) {
                                cover = p_matrix[sortcol[c2] * rows + r];
                            }
                            r++;
                        }
                        if (cover) {
                            p_cols[sortcol[c2]] = false;
                            --(*survcols);
                        }
                    }
                }
            }
        }
    }
    if (*survcols < cols) {
        int s = 0;
        for (int c = 0; c < cols; c++) {
            if (p_cols[c]) {
                for (int r = 0; r < rows; r++) {
                    p_matrix[s * rows + r] = p_matrix[c * rows + r];
                }
                s++;
            }
        }
    }
}
