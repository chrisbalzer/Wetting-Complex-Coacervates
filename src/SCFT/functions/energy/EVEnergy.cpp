#include "../../SCFT.h"

void SCFT::EVEnergy(double *phi_k, double& f_k) {
    // Calculate total volume fraction and then energy
    double totalVolFrac = 0.0;
    for(int i = 0; i < numComponents; i++) {
        totalVolFrac += phi_k[i];
    }
    f_k = (1 - totalVolFrac)*(log(1 - totalVolFrac) - 1);
}
