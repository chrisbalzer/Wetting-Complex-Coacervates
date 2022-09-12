#include "../SCFT.h"

void SCFT::calculateErrorRelax(double* F_ik, double* X_ik) {
    double temp = 0.0;

    // Set the error from the last iteration
    error_prev = error;

    // Calculate error (maximum relative deviation)
    error = 0.0;
    for (int i = 0; i < numComponents+1; i++) {
       for (int k = 0; k < M; k++) {
            temp = fabs(F_ik[i*M + k]);
            if (temp > error) error = temp;
        }
    }
}
