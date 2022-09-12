#include "../SCFT.h"

void SCFT::initializeDensity() {
      // Initialize to dilute bulk concentration
      for (int j = 0; j < numComponents; j++) {
          for (int k = 0; k < M; k++) {
              phi[j][k] = phiB[j];
          }
      }

      // Calculate the order parameter
      for (int j = 0; j < numPolymers; j++) {
          for (int k = 0; k < M; k++) {
              varphi[j][k] = sqrt(phiB[j]);
          }
      }

      // Fill array X with these values
      refillX();
}
