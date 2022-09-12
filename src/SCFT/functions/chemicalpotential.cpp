#include "../SCFT.h"

void SCFT::chemicalPotentialbulk() {
    // Define variables for contributions to chemical potential
    double *muID, *muDH, *muEV, *muDHChain;
    muID      = new double[numComponents];
    muDH      = new double[numComponents];
    muEV      = new double[numComponents];
    muDHChain = new double[numComponents];

    // Calculate contributions to the free energy
    IdealPotential(phiB,muID);
    DHPotential(phiB,muDH);
    EVPotential(phiB,muEV);
    DHChainPotential(phiB,muDHChain);

    // Sum parts of the chemical potential
    for (int i = 0; i < numComponents; i++) {
        mu[i]   = muID[i] + muDH[i] + muEV[i] + muDHChain[i];
        muEx[i] = muDH[i] + muEV[i] + muDHChain[i];
    }

    delete [] muID;
    delete [] muDH;
    delete [] muEV;
    delete [] muDHChain;
}


void SCFT::chemicalPotentiallocal(double** mu_local_jk, double** mu_localEx_jk, double**phi_jk) {
    // Define variables for contributions to chemical potential
    double **muID, **muDH, **muEV, **muDHChain;
    muID      = new double*[numComponents];
    muDH      = new double*[numComponents];
    muEV      = new double*[numComponents];
    muDHChain = new double*[numComponents];

    for(int i = 0; i < numComponents; i++) {
        muID[i]      = new double[M];
        muDH[i]      = new double[M];
        muEV[i]      = new double[M];
        muDHChain[i] = new double[M];
    }

    // Calculate contributions to the free energy
    IdealPotential(phi_jk,muID);
    DHPotential(phi_jk,muDH);
    EVPotential(phi_jk,muEV);
    DHChainPotential(phi_jk,muDHChain);

    // Sum parts of the chemical potential
    for (int i = 0; i < numComponents; i++) {
      for (int k = 0; k < M; k++) {
          mu_local_jk[i][k]   = muID[i][k] + muDH[i][k] + muEV[i][k] + muDHChain[i][k];
          mu_localEx_jk[i][k] = muDH[i][k] + muEV[i][k] + muDHChain[i][k];
      }
    }

    for(int i = 0; i < numComponents; i++) {
        delete [] muID[i];
        delete [] muDH[i];
        delete [] muEV[i];
        delete [] muDHChain[i];
    }
    delete [] muID;
    delete [] muDH;
    delete [] muEV;
    delete [] muDHChain;
}


void SCFT::chemicalPotentiallocalParallel(double** mu_local_jk, double** mu_localEx_jk, double**phi_jk) {
    #pragma omp parallel
    {
        double totalVolFrac;
        double kappa, chi_i, chi_j, temp0, temp1, temp2, temp3,factor;
        double lnDelta, DlnDeltaDphi;
        double phi_k[numComponents];
        double mu_temp;
        double muEx_temp;

        // Loop through each spatial position
        #pragma omp for
        for (int k = 0; k < M; k++) {
            // Get the density for this k
            for (int j = 0; j < numComponents; j++) {
                phi_k[j] = phi_jk[j][k];
            }

            // Calculate total volume fraction for excluded volume
            totalVolFrac = 0.0;
            for (int j = 0; j < numComponents; j++) {
                totalVolFrac += phi_k[j];
            }

            // Debye length
            DHScreeningLength(phi_k,kappa);

            // Calculate parameters for DH theory
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

            // Loop through components to calculate chemical potential
            for (int j = 0; j < numComponents; j++) {
                // Excluded volume
                mu_temp   = -log(1 - totalVolFrac);
                muEx_temp = -log(1 - totalVolFrac);

                // DH electrostatic correlation
                chi_j = 1/D[j]/D[j]/D[j] * (log(1 + kappa*D[j]) - kappa*D[j] + 0.5*kappa*kappa*D[j]*D[j]);
                factor = (chi_j + 2*Pi*kappa*BJ*temp3)/temp1  - temp2/(temp1*temp1);
                mu_temp   += -Z[j]*Z[j]*factor/(4.0*Pi);
                muEx_temp += -Z[j]*Z[j]*factor/(4.0*Pi);

                // Chain connectivity
                lnDelta = -Z[j]*Z[j]*BJ/D[j]/(1.0 + kappa*D[j]);
                for (int i = 0; i < numPolymers; i++) {
                    // Calcluate derivative of ln(delta_i) with respect to volume fraction_j
                    DlnDeltaDphi = 2.0*Pi*Z[i]*Z[i]*Z[j]*Z[j]*BJ*BJ/D[j]/D[j]/(1 + kappa*D[i])/(1 + kappa*D[i])/(kappa*D[j]);

                    // Add contributions to chemical potential
                    mu_temp += (1.0/N[i] - 1)*phi_k[i]*DlnDeltaDphi;
                    if(i == j) mu_temp += (1.0/N[i] - 1)*lnDelta;
                    muEx_temp += (1.0/N[i] - 1)*phi_k[i]*DlnDeltaDphi;
                    if(i == j) muEx_temp += (1.0/N[i] - 1)*lnDelta;
                }

                  // Ideal part of chemical potential
                  if (j < numPolymers) mu_temp += 1.0/N[j]*log(phi_k[j]/N[j]);
                  else mu_temp += log(phi_k[j]);

                  // Assign value in outside array
                  mu_local_jk[j][k]   = mu_temp;
                  mu_localEx_jk[j][k] = muEx_temp;
              }
          }
      }
}
