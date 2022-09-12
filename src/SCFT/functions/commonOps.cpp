#include "../SCFT.h"

void SCFT::refillX() {
    for (int i = 0; i < numComponents+1; i++) {
       for (int k = 0; k < M; k++) {
            if (i < numPolymers) {
              X[i*M + k] = sqrt(phi[i][k]);
            }
            else if (i < numComponents) {
                X[i*M + k] = phi[i][k];
            }
            else {
                X[i*M + k] = Psi[k];
            }
       }
    }
}

void SCFT::getX() {
    for (int i = 0; i < numComponents+1; i++) {
       for (int k = 0; k < M; k++) {
            if(i < numPolymers) {
                phi[i][k] = X[i*M + k]*X[i*M + k];
                varphi[i][k] = X[i*M + k];
            }
            else if (i < numComponents) {
                phi[i][k] = X[i*M + k];
            }
            else {
                Psi[k] = X[i*M + k];
            }
       }
    }
}


void SCFT::calcChargeDensity(double* rhoZ_k, double** phi_jk) {
    for(int k = 0; k < M; k++) {
      // Set to zero
      rhoZ_k[k] = 0.0;

      // Add contribution from each species
      for (int i = 0; i < numComponents; i++) {
          rhoZ_k[k] += Z[i]*phi_jk[i][k]/D[i]/D[i]/D[i];
      }
    }
}

void SCFT::densityNAN() {
    for (int i = 0; i < numComponents; i++) {
        for(int k = 0; k < M; k++) {
            if(isnan(phi[i][k])) {
                std::cout<<"SCFT - Density is NaN"<<std::endl;
                throw 1;
            }
        }
    }
}
