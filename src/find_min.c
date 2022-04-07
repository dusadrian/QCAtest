#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "lp_min.h"
#include "find_min.h"
#include <R_ext/RS.h> 
void find_min(
    const int p_pichart[],
    const int pirows,
    const unsigned int picols,
    int *solmin,
    int p_indices[]
) {
    double *p_objective = R_Calloc(picols + 1, double);
    p_objective[0] = 0;
    for (int i = 0; i < picols; i++) {
        p_objective[i + 1] = 1;
    }
    double *p_constraints = R_Calloc(pirows * (picols + 2) + 1, double);
    p_constraints[0] = 0;
    int pos = 1;
    for (int r = 0; r < pirows; r++) { 
        for (int c = 0; c < picols; c++) {
            p_constraints[pos] = p_pichart[c * pirows + r] * 1;
            pos++;
        }
        p_constraints[pos] = 2;
        p_constraints[pos + 1] = 1;
        pos = pos + 2;
    }
    int *p_bin_vec = R_Calloc(picols, int);
    for (int i = 0; i < picols; i++) {
        p_bin_vec[i] = i + 1;
    }
    double *p_obj_val = R_Calloc(1, double);
    double *p_solution = R_Calloc(picols + 1, double);
    int *p_status = R_Calloc(1, int);
    lp_min(p_objective,
            pirows, 
            p_constraints,
            picols, 
            p_bin_vec,
            p_obj_val,
            p_solution,
            p_status);
    *solmin = 0;
    for (int c = 0; c < picols; c++) {
        if (p_solution[c] > 0) { 
            p_indices[*solmin] = c;
            *solmin += 1;
        }
    }
    R_Free(p_objective);
    R_Free(p_constraints);
    R_Free(p_bin_vec);
    R_Free(p_obj_val);
    R_Free(p_solution);
    R_Free(p_status);
    return;
}
