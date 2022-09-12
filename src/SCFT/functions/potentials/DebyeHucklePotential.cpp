#include "../../SCFT.h"

void SCFT::DHPotential(double *phi_k, double* mu_k) {
    double kappa, chi_i, chi_j, temp0, temp1, temp2, temp3,factor;

    // Calculate the inverse Debye length
    DHScreeningLength(phi_k,kappa);

    // Calculate parameters
    temp1 = 0.0;   // sum z_i^2 rho_i
    temp2 = 0.0;   // sum z_i^2 rho_i chi_i
    temp3 = 0.0;   // sum z_i^2 rho_i/(1 + kappa sigma_i)
    for (int i = 0; i < numComponents; i++) {
        temp0  = Z[i]*Z[i]*phi_k[i]/D[i]/D[i]/D[i];  // z_i^2 rho_i
        chi_i = 1/D[i]/D[i]/D[i] * (log(1 + kappa*D[i]) - kappa*D[i] + 0.5*kappa*kappa*D[i]*D[i]);

        temp1 += temp0;
        temp2 += temp0 * chi_i;
        temp3 += temp0 / (1 + kappa*D[i]);
    }

    // Calculate chemical potential (general form)
    for (int j = 0; j < numComponents; j++) {
        chi_j = 1/D[j]/D[j]/D[j] * (log(1 + kappa*D[j]) - kappa*D[j] + 0.5*kappa*kappa*D[j]*D[j]);
        factor = (chi_j + 2*Pi*kappa*BJ*temp3)/temp1  - temp2/(temp1*temp1);

        mu_k[j] = -Z[j]*Z[j]*factor/(4.0*Pi);
    }
}

void SCFT::DHPotential(double **phi_jk, double** mu_jk) {
    double kappa, chi_i, chi_j, temp0, temp1, temp2, temp3,factor;
    double* phi_k;
    phi_k = new double[numComponents];

    for(int k = 0; k < M; k++) {
        for (int j = 0; j < numComponents; j++) {
            phi_k[j] = phi_jk[j][k];
        }

        // Calculate the inverse Debye length
        DHScreeningLength(phi_k,kappa);

        // Calculate parameters
        temp1 = 0.0;   // sum z_i^2 rho_i
        temp2 = 0.0;   // sum z_i^2 rho_i chi_i
        temp3 = 0.0;   // sum z_i^2 rho_i/(1 + kappa sigma_i)
        for (int i = 0; i < numComponents; i++) {
            temp0  = Z[i]*Z[i]*phi_k[i]/D[i]/D[i]/D[i];  // z_i^2 rho_i
            chi_i = 1/D[i]/D[i]/D[i] * (log(1 + kappa*D[i]) - kappa*D[i] + 0.5*kappa*kappa*D[i]*D[i]);

            temp1 += temp0;
            temp2 += temp0 * chi_i;
            temp3 += temp0 / (1 + kappa*D[i]);
        }

        // Calculate chemical potential (general form)
        for (int j = 0; j < numComponents; j++) {
            chi_j = 1/D[j]/D[j]/D[j] * (log(1 + kappa*D[j]) - kappa*D[j] + 0.5*kappa*kappa*D[j]*D[j]);
            factor = (chi_j + 2*Pi*kappa*BJ*temp3)/temp1  - temp2/(temp1*temp1);

            mu_jk[j][k] = -Z[j]*Z[j]*factor/(4.0*Pi);
        }
    }

    delete [] phi_k;
}
