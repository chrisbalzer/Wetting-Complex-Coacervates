#include "../SCFT.h"

void SCFT::updateDensityRelax() {
    /*=============================
      Calculate based on current density profile
    =============================*/
    // Calculate local chemical potential based on current density
    if(numProcessors == 1) {
         chemicalPotentiallocal(mu_local,mu_localEx,phi);
    }
    else {
        chemicalPotentiallocalParallel(mu_local,mu_localEx,phi);
    }

    // Calculate the charge density
    calcChargeDensity(rhoZ, phi);

    /*=============================
     Calculate the error terms
    =============================*/
    // For convenience
    int indx;
    double dz = gridSize;

    // Calculate the difference (higher order)
    for (int i = 0; i < numComponents+1; i++) {
       for (int k = 0; k < M; k++) {
            indx = i*M + k;

            if(i < numPolymers) {
                if(k == 0) {
                    // Forward difference (4th order)
                    F[indx] = -uSurf[i] + 1.0/varphi[i][k]/6.0*(-0.25*varphi[i][k+4] + 4.0/3.0*varphi[i][k+3] - 3.0*varphi[i][k+2] + 4.0*varphi[i][k+1] - 25.0/12.0*varphi[i][k])/dz;
                }
                else if(k == 1 || k == M-2) {
                    // Centered difference (2nd order)
                    F[indx] =  mu[i] - mu_local[i][k] - Z[i]*Psi[k] + 1.0/varphi[i][k]/6.0*(varphi[i][k+1] - 2.0*varphi[i][k] + varphi[i][k-1])/dz/dz;
                }
                else if(k == M-1) {
                    F[indx] = mu[i] - mu_local[i][k];
                }
                else {
                    // Centered difference (4th order)
                    F[indx] =  mu[i] - mu_local[i][k] - Z[i]*Psi[k] + 1.0/varphi[i][k]/6.0*(-varphi[i][k+2] + 16.0*varphi[i][k+1] - 30.0*varphi[i][k] + 16.0*varphi[i][k-1] - varphi[i][k-2])/12.0/dz/dz;
                }
                F[indx] *= 2.0*varphi[i][k];
            }

            // Calculate salt ion error
            else if (i < numComponents) {

                if(k == M-1) {
                    F[indx] = mu[i] - mu_local[i][k];
                }
                else {
                    F[indx] = mu[i] - mu_local[i][k] - Z[i]*Psi[k];
                }
            }
            // Calculate potential error
            else {
                if(k == 0) {
                    if(BC_type == 0) {
                        F[indx] = surfBC - Psi[k];
                    }
                    else{
                        // Forward difference (4th order)
                        F[indx] = 4.0*Pi*BJ*surfBC + (-0.25*Psi[k+4] + 4.0/3.0*Psi[k+3] - 3.0*Psi[k+2] + 4.0*Psi[k+1] - 25.0/12.0*Psi[k])/dz;
                    }
                }
                else if(k == 1 || k == M-2) {
                    // Centered difference (2nd order)
                    F[indx] = (Psi[k+1] - 2.0*Psi[k] + Psi[k-1])/dz/dz + 4.0*Pi*BJ*rhoZ[k];
                }
                else if(k == M-1) {
                    F[indx] = 0.0;
                }
                else {
                    // Centered difference (4th order)
                    F[indx] = (-Psi[k+2]  + 16.0*Psi[k+1] - 30.0*Psi[k] + 16.0*Psi[k-1] - Psi[k-2])/12.0/dz/dz + 4.0*Pi*BJ*rhoZ[k];
                }
            }
       }
    }

    // Calculate error
    calculateErrorRelax(F,X);

    // Mix density
    mixDensityRelax();
}


void SCFT::calcFValues(double* F_temp, double* X_temp) {
    /*=============================
      Initialize temporary arrays
    =============================*/
    double*  Psi_k;
    double*  rhoZ_k;
    double**  phi_jk;
    double**  varphi_jk;
    double** mu_local_jk;
    double** mu_localEx_jk;

    Psi_k         = new double[M];
    rhoZ_k        = new double[M];
    phi_jk        = new double*[numComponents];
    mu_local_jk   = new double*[numComponents];
    mu_localEx_jk = new double*[numComponents];
    varphi_jk     = new double*[numPolymers];

    for (int i = 0; i < numComponents; i++) {
        phi_jk[i]        = new double[M];
        mu_local_jk[i]   = new double[M];
        mu_localEx_jk[i] = new double[M];
    }

    for (int i = 0; i < numPolymers; i++) {
        varphi_jk[i]    = new double[M];
    }

    /*=============================
      Assign phi and Psi from X
    =============================*/
    for (int i = 0; i < numComponents+1; i++) {
       for (int k = 0; k < M; k++) {
            if(i < numPolymers) {
                phi_jk[i][k] = X_temp[i*M + k]*X_temp[i*M + k];
                varphi_jk[i][k] = X_temp[i*M + k];
            }
            else if (i < numComponents) {
                phi_jk[i][k] = X_temp[i*M + k];
            }
            else {
                Psi_k[k] = X_temp[i*M + k];
            }
       }
    }

    /*=============================
      Calculate based on current density profile
    =============================*/
    // Calculate local chemical potential based on current density
    if(numProcessors == 1) {
         chemicalPotentiallocal(mu_local_jk,mu_localEx_jk,phi_jk);
    }
    else {
        chemicalPotentiallocalParallel(mu_local_jk,mu_localEx_jk,phi_jk);
    }

    // Calculate the charge density
    calcChargeDensity(rhoZ_k, phi_jk);

    /*=============================
     Calculate the error terms
    =============================*/
    // For convenience
    int indx;
    double dz = gridSize;

    // Calculate the difference
    for (int i = 0; i < numComponents+1; i++) {
         for (int k = 0; k < M; k++) {
              indx = i*M + k;

              if(i < numPolymers) {
                    // 2nd order accurate in space
                    if(k == 0) {
                        F_temp[indx] = -uSurf[i] + 1.0/varphi_jk[i][k]/6.0*(-varphi_jk[i][k+2] + 4.0*varphi_jk[i][k+1] - 3.0*varphi_jk[i][k])/dz/2.0;
                    }
                    else if(k == M-1) {
                        F_temp[indx] = mu[i] - mu_local_jk[i][k];
                    }
                    else {
                        F_temp[indx] =  mu[i] - mu_local_jk[i][k] - Z[i]*Psi_k[k] + 1.0/varphi_jk[i][k]/6.0*(varphi_jk[i][k+1] - 2.0*varphi_jk[i][k] + varphi_jk[i][k-1])/dz/dz;
                    }
                    F_temp[indx] *= 2*varphi_jk[i][k];
              }

              // Calculate salt ion error
              else if (i < numComponents) {
                  if(k == M-1) {
                      F_temp[indx] = mu[i] - mu_local_jk[i][k];
                  }
                  else {
                      F_temp[indx] = mu[i] - mu_local_jk[i][k] - Z[i]*Psi_k[k];
                  }
              }
              // Calculate potential error
              else {
                  if(k == 0) {
                      if(BC_type == 0) {
                          F_temp[indx] = surfBC - Psi_k[k];
                      }
                      else{
                          F_temp[indx] = 4.0*Pi*BJ*surfBC + (-Psi_k[k+2] + 4.0*Psi_k[k+1] - 3.0*Psi_k[k])/dz/2.0;
                      }
                  }
                  else if(k == M-1) {
                      F_temp[indx] = 0.0;
                  }
                  else {
                      F_temp[indx] = (Psi_k[k+1] - 2.0*Psi_k[k] + Psi_k[k-1])/dz/dz + 4.0*Pi*BJ*rhoZ_k[k];
                  }
              }
         }
    }

    /*=============================
      Delete temporary arrays
    =============================*/
    for (int i = 0; i < numComponents; i++) {
        delete [] phi_jk[i];
        delete [] mu_local_jk[i];
        delete [] mu_localEx_jk[i];
    }
    for (int i = 0; i < numPolymers; i++) {
        delete [] varphi_jk[i];
    }
    delete [] Psi_k;
    delete [] rhoZ_k;
    delete [] varphi_jk;
    delete [] mu_local_jk;
    delete [] mu_localEx_jk;
    delete [] phi_jk;
}


void SCFT::calcFValuesHO(double* F_temp, double* X_temp) {
    /*=============================
      Initialize temporary arrays
    =============================*/
    double*  Psi_k;
    double*  rhoZ_k;
    double**  phi_jk;
    double**  varphi_jk;
    double** mu_local_jk;
    double** mu_localEx_jk;

    Psi_k         = new double[M];
    rhoZ_k        = new double[M];
    phi_jk        = new double*[numComponents];
    mu_local_jk   = new double*[numComponents];
    mu_localEx_jk = new double*[numComponents];
    varphi_jk     = new double*[numPolymers];

    for (int i = 0; i < numComponents; i++) {
        phi_jk[i]        = new double[M];
        mu_local_jk[i]   = new double[M];
        mu_localEx_jk[i] = new double[M];
    }

    for (int i = 0; i < numPolymers; i++) {
        varphi_jk[i]    = new double[M];
    }

    /*=============================
      Assign phi and Psi from X
    =============================*/
    for (int i = 0; i < numComponents+1; i++) {
       for (int k = 0; k < M; k++) {
            if(i < numPolymers) {
                phi_jk[i][k] = X_temp[i*M + k]*X_temp[i*M + k];
                varphi_jk[i][k] = X_temp[i*M + k];
            }
            else if (i < numComponents) {
                phi_jk[i][k] = X_temp[i*M + k];
            }
            else {
                Psi_k[k] = X_temp[i*M + k];
            }
       }
    }

    /*=============================
      Calculate based on current density profile
    =============================*/
    // Calculate local chemical potential based on current density
    if(numProcessors == 1) {
         chemicalPotentiallocal(mu_local_jk,mu_localEx_jk,phi_jk);
    }
    else {
        chemicalPotentiallocalParallel(mu_local_jk,mu_localEx_jk,phi_jk);
    }

    // Calculate the charge density
    calcChargeDensity(rhoZ_k, phi_jk);

    /*=============================
     Calculate the error terms
    =============================*/
    // For convenience
    int indx;
    double dz = gridSize;

    // Calculate the difference
    for (int i = 0; i < numComponents+1; i++) {
         for (int k = 0; k < M; k++) {
              indx = i*M + k;

              if(i < numPolymers) {
                  if(k == 0) {
                      // 4th order accurate
                      F_temp[indx] = -uSurf[i] + 1.0/varphi_jk[i][k]/6.0*(-0.25*varphi_jk[i][k+4] + 4.0/3.0*varphi_jk[i][k+3] - 3.0*varphi_jk[i][k+2] + 4.0*varphi_jk[i][k+1] - 25.0/12.0*varphi_jk[i][k])/dz;
                  }
                  else if(k == 1 || k == M-2) {
                      // Centered difference (2nd order)
                      F_temp[indx] =  mu[i] - mu_local_jk[i][k] - Z[i]*Psi_k[k] + 1/varphi_jk[i][k]/6.0*(varphi_jk[i][k+1] - 2.0*varphi_jk[i][k] + varphi_jk[i][k-1])/dz/dz;
                  }
                  else if(k == M-1) {
                      F_temp[indx] = mu[i] - mu_local_jk[i][k];
                  }
                  else {
                      // Centered difference (4th order)
                      F_temp[indx] =  mu[i] - mu_local_jk[i][k] - Z[i]*Psi_k[k] + 1/varphi_jk[i][k]/6.0*(-varphi_jk[i][k+2] + 16.0*varphi_jk[i][k+1] - 30.0*varphi_jk[i][k] + 16.0*varphi_jk[i][k-1] - varphi_jk[i][k-2])/12.0/dz/dz;
                  }
                  F_temp[indx] *= 2.0*varphi_jk[i][k];
              }

              // Calculate salt ion error
              else if (i < numComponents) {
                  if(k == M-1) {
                      F_temp[indx] = mu[i] - mu_local_jk[i][k];
                  }
                  else {
                      F_temp[indx] = mu[i] - mu_local_jk[i][k] - Z[i]*Psi_k[k];
                  }
              }
              // Calculate potential error
              else {
                  if(k == 0) {
                      if(BC_type == 0) {
                          F_temp[indx] = surfBC - Psi_k[k];
                      }
                      else{
                          // Forward difference (4th order)
                          F_temp[indx] = 4.0*Pi*BJ*surfBC + (-0.25*Psi_k[k+4] + 4.0/3.0*Psi_k[k+3] - 3.0*Psi_k[k+2] + 4.0*Psi_k[k+1] - 25.0/12.0*Psi_k[k])/dz;
                      }
                  }
                  else if(k == 1 || k == M-2) {
                      // Centered difference (2nd order)
                      F_temp[indx] = (Psi_k[k+1] - 2.0*Psi_k[k] + Psi_k[k-1])/dz/dz + 4.0*Pi*BJ*rhoZ_k[k];
                  }
                  else if(k == M-1) {
                      F_temp[indx] = 0.0;
                  }
                  else {
                      // Centered difference (4th order)
                      F_temp[indx] = (-Psi_k[k+2]  + 16.0*Psi_k[k+1] - 30.0*Psi_k[k] + 16.0*Psi_k[k-1] - Psi_k[k-2])/12.0/dz/dz + 4.0*Pi*BJ*rhoZ_k[k];
                  }
              }
         }
    }

    /*=============================
      Delete temporary arrays
    =============================*/
    for (int i = 0; i < numComponents; i++) {
        delete [] phi_jk[i];
        delete [] mu_local_jk[i];
        delete [] mu_localEx_jk[i];
    }
    for (int i = 0; i < numPolymers; i++) {
        delete [] varphi_jk[i];
    }
    delete [] Psi_k;
    delete [] rhoZ_k;
    delete [] varphi_jk;
    delete [] mu_local_jk;
    delete [] mu_localEx_jk;
    delete [] phi_jk;
}
