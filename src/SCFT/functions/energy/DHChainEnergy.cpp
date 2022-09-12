#include "../../SCFT.h"

void SCFT::DHChainEnergy(double *phi_k, double& f_k) {
    double lnDelta, kappa;

    // Initialize energy to zero
    f_k = 0.0;

    // Calculate inverse Debye screening length
    DHScreeningLength(phi_k, kappa);

    // Loop through polymers
    for (int j = 0; j < numPolymers; j++) {
        // Calculate logarithm of delta
        lnDelta = -Z[j]*Z[j]*BJ/D[j]/(1.0 + kappa*D[j]);

        // Add contribution from polymer
        f_k += (1.0/N[j] - 1)*phi_k[j]*lnDelta;
    }
}
