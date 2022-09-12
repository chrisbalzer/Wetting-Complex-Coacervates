#include "../../SCFT.h"

void SCFT::EVPotential(double *phi_k, double* mu_k) {
  double totalVolFrac = 0.0;
  for(int i = 0; i < numComponents; i++) {
      totalVolFrac += phi_k[i];
  }
  for(int i = 0; i < numComponents; i++) {
      mu_k[i] = -log(1 - totalVolFrac);
  }
}


void SCFT::EVPotential(double **phi_jk, double** mu_jk) {
    double totalVolFrac;
    double* phi_k;
    phi_k = new double[numComponents];

    for(int k = 0; k < M; k++) {
        totalVolFrac = 0.0;
        for (int j = 0; j < numComponents; j++) {
            totalVolFrac += phi_jk[j][k];
        }

        for(int j = 0; j < numComponents; j++) {
            mu_jk[j][k] = -log(1 - totalVolFrac);
        }
    }

    delete [] phi_k;
}
