#include "../../SCFT.h"

void SCFT::DHChainPotential(double *phi_k, double* mu_k) {
    double lnDelta, DlnDeltaDphi, kappa;

    // Calculate the inverse Debye length
    DHScreeningLength(phi_k,kappa);

    for (int j = 0; j < numComponents; j++) {
        // Initialize to zero
        mu_k[j] = 0.0;

        // Calculate logarithm of delta
        lnDelta = -Z[j]*Z[j]*BJ/D[j]/(1.0 + kappa*D[j]);

        for (int i = 0; i < numPolymers; i++) {
            // Calcluate derivative of ln(delta_i) with respect to volume fraction_j
            DlnDeltaDphi = 2.0*Pi*Z[i]*Z[i]*Z[j]*Z[j]*BJ*BJ/D[j]/D[j]/(1 + kappa*D[i])/(1 + kappa*D[i])/(kappa*D[j]);

            // Add contributions to chemical potential
            mu_k[j] += (1.0/N[i] - 1)*phi_k[i]*DlnDeltaDphi;
            if(i == j) mu_k[j] += (1.0/N[i] - 1)*lnDelta;
        }
    }
}

void SCFT::DHChainPotential(double **phi_jk, double** mu_jk) {
    double lnDelta, DlnDeltaDphi, kappa;
    double* phi_k;
    phi_k = new double[numComponents];

    for(int k = 0; k < M; k++) {
        for (int j = 0; j < numComponents; j++) {
            phi_k[j] = phi_jk[j][k];
        }

        // Calculate the inverse Debye length for position
        DHScreeningLength(phi_k,kappa);

        for (int j = 0; j < numComponents; j++) {
            // Initialize to zero
            mu_jk[j][k] = 0.0;

            // Calculate logarithm of delta
            lnDelta = -Z[j]*Z[j]*BJ/D[j]/(1.0 + kappa*D[j]);

            for (int i = 0; i < numPolymers; i++) {
                // Calcluate derivative of ln(delta_i) with respect to volume fraction_j
                DlnDeltaDphi = 2.0*Pi*Z[i]*Z[i]*Z[j]*Z[j]*BJ*BJ/D[j]/D[j]/(1 + kappa*D[i])/(1 + kappa*D[i])/(kappa*D[j]);

                // Add contributions to chemical potential
                mu_jk[j][k] += (1.0/N[i] - 1)*phi_k[i]*DlnDeltaDphi;
                if(i == j) mu_jk[j][k] += (1.0/N[i] - 1)*lnDelta;
            }
        }
    }

    delete [] phi_k;
}
