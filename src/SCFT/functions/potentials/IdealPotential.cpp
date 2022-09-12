#include "../../SCFT.h"

void SCFT::IdealPotential(double *phi_k, double* mu_k) {
    for(int i = 0; i < numComponents; i ++) {
        if (i < numPolymers) mu_k[i] = 1.0/N[i]*log(phi_k[i]/N[i]);
        else mu_k[i] = log(phi_k[i]);
    }
}

void SCFT::IdealPotential(double **phi_jk, double** mu_jk) {
    double* phi_k;
    phi_k = new double[numComponents];

    for(int k = 0; k < M; k++) {
        for(int j = 0; j < numComponents; j ++) {
            if (j < numPolymers) mu_jk[j][k] = 1.0/N[j]*log(phi_jk[j][k]/N[j]);
            else mu_jk[j][k] = log(phi_jk[j][k]);
        }
    }

    delete [] phi_k;
}
