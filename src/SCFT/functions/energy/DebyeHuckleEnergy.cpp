#include "../../SCFT.h"

void SCFT::DHEnergy(double *phi_k, double& f_k) {
    double kappa, denom, numer, chi_i, temp;

    // Calculate the inverse Debye length
    DHScreeningLength(phi_k,kappa);

    // Calculate parameters
    denom = 0.0;
    numer = 0.0;
    for (int i = 0; i < numComponents; i++) {
        temp  = Z[i]*Z[i]*phi_k[i]/D[i]/D[i]/D[i];
        chi_i = 1/D[i]/D[i]/D[i] * (log(1 + kappa*D[i]) - kappa*D[i] + 0.5*kappa*kappa*D[i]*D[i]);

        denom += temp;
        numer += temp * chi_i;
    }

    // Calculate free energy
    f_k = -numer/(4*Pi*denom);
}

void SCFT::DHScreeningLength(double* phi_k,double& kappa) {
    kappa = 0.0;
    for (int i = 0; i < numComponents; i++) {
        kappa += Z[i]*Z[i]*phi_k[i]/D[i]/D[i]/D[i];
    }
    kappa = sqrt(4*Pi*BJ*kappa);
}
