#include "../../SCFT.h"

void SCFT::IdealEnergy(double *phi_k, double& f_k) {
    // Initialize free energy
    f_k = 0.0;

    // Add each component
    for(int i = 0; i < numComponents; i ++) {
        if (i < numPolymers) f_k += phi_k[i]/N[i]*(log(phi_k[i]/N[i]) - 1);
        else f_k +=  phi_k[i]*(log(phi_k[i]) - 1);
    }
}
