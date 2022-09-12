#include "../SCFT.h"


void SCFT::calcProps() {
    // Open data file
    dataStream.open(dataFile.c_str(), std::fstream::out | std::ofstream::app);

    // Print surface charge density and surface potential
    double surfQ = -1/(4*Pi*BJ)*(-0.25*Psi[4] + 4.0/3.0*Psi[3] - 3.0*Psi[2] + 4.0*Psi[1] - 25.0/12.0*Psi[0])/gridSize;
    if(fabs(surfQ) < errorTol0) surfQ = 0;
    double Psi0;
    Psi0 = Psi[0];
    if(fabs(Psi0) < errorTol0) Psi0 = 0;
    dataStream << std::scientific << std::setprecision(5) << surfQ << " ";
    dataStream << std::scientific << std::setprecision(5) << Psi0 << " ";

    // Print the polymer boundary conditions (eta)
    for(int i = 0; i < numPolymers; i++) {
          dataStream << std::scientific << std::setprecision(5) << uSurf[i] << " ";
    }

    // Calculate the excess adsorption for each component
    calcExcessAdsorption();
    for(int i = 0; i < numComponents; i++) {
        // Add to outData
        dataStream << std::scientific << std::setprecision(5) << exAdsStore[i][stepIter%keepPhase] << " ";
    }

    // Calculate the surface tension
    double temp;
    double f_bulk;
    double P_bulk;
    double f_temp;
    double temp_grad;
    double chemPot_part;
    double uSurf_part;
    double integrand1,integrand2, integrand3;
    double GrandPotential;
    double Helmholtz;
    double gamma;

    double* phi_k = new double[numComponents];
    double* integrand1_z = new double[M];
    double* integrand2_z = new double[M];
    double* integrand3_z = new double[M];

    // Calculate bulk free energy density
    for (int j = 0; j < numComponents; j++) {
        phi_k[j] = phi[j][M-1];
    }
    calcEnergy(phi_k, f_bulk);

    // Calculate bulk pressure
    P_bulk = 0.0;
    for(int j = 0; j < numComponents; j++) {
        P_bulk += mu[j]*phi_k[j];
    }
    P_bulk -= f_bulk;

    // Calculate gradient of polymer parts
    uSurf_part = 0.0;
    for(int j = 0; j < numPolymers; j++) {
        uSurf_part += uSurf[j]*phi[j][0];
    }

    for(int k = 0; k < M; k++) {
        // Get density for point k
        for (int j = 0; j < numComponents; j++) {
            phi_k[j] = phi[j][k];
        }

        // Calculate free energy density at point k
        calcEnergy(phi_k, f_temp);

        // Calculate gradient of polymer parts
        temp_grad = 0.0;
        for(int j = 0; j < numPolymers; j++) {
            temp = 0.0;
            if(k == 0) {
                // 4th order forward difference
                temp = (-0.25*sqrt(phi[j][k+4]) + 4.0/3.0*sqrt(phi[j][k+3]) - 3.0*sqrt(phi[j][k+2]) + 4.0*sqrt(phi[j][k+1]) - 25.0/12.0*sqrt(phi[j][k]))/gridSize;
            }
            else if(k == 1 || k == M-2) {
                // 2nd order central difference
                temp = (sqrt(phi[j][k+1]) - sqrt(phi[j][k-1]))/2.0/gridSize;
            }
            else if(k == M-1) {
                // 2nd order backward difference
                temp = (sqrt(phi[j][k-2]) - 4.0*sqrt(phi[j][k-1]) + 3.0*sqrt(phi[j][k]))/2.0/gridSize;
            }
            else {
                // 4th order central difference
                temp = (-sqrt(phi[j][k+2]) + 8.0*sqrt(phi[j][k+1]) - 8.0*sqrt(phi[j][k-1]) + sqrt(phi[j][k-2]))/12.0/gridSize;
            }
            temp_grad += temp*temp/6.0;
        }

        // Chemical potential part
        chemPot_part = 0.0;
        for (int j = 0; j < numComponents; j++) {
            chemPot_part -= mu[j]*phi_k[j];
        }

        // Put it all together
        integrand1 = f_temp + temp_grad + 0.5*rhoZ[k]*Psi[k];
        integrand2 = f_temp + temp_grad + 0.5*rhoZ[k]*Psi[k] + chemPot_part;
        integrand3 = integrand2 + P_bulk;

        integrand1_z[k] = integrand1;
        integrand2_z[k] = integrand2;
        integrand3_z[k] = integrand3;
    }
    Helmholtz      = simps38(integrand1_z,0,M-1)*gridSize + uSurf_part + 0.5*surfQ*Psi[0];
    GrandPotential = simps38(integrand2_z,0,M-1)*gridSize + uSurf_part + 0.5*surfQ*Psi[0];
    gamma          = simps38(integrand3_z,0,M-1)*gridSize + uSurf_part - 0.5*surfQ*Psi[0];

    dataStream << std::scientific << std::setprecision(8) << Helmholtz << " ";
    dataStream << std::scientific << std::setprecision(8) << GrandPotential << " ";
    dataStream << std::scientific << std::setprecision(8) << P_bulk << " ";
    dataStream << std::scientific << std::setprecision(8) << gamma << " ";

    delete [] integrand1_z;
    delete [] integrand2_z;
    delete [] integrand3_z;
    delete [] phi_k;

    // Print bulk density to output file
    for(int i = 0; i < numComponents; i++) {
        // Add to outData
        dataStream << std::scientific << std::setprecision(8) << phiB[i] << " ";
    }

    // Print chemical potential to output file
    for(int i = 0; i < numComponents; i++) {
        // Add to outData
        dataStream << std::scientific << std::setprecision(8) << mu[i] << " ";
    }

    // Check if a phase transition has occured
    checkPhaseTransition();

    // Add line for s if doing continuation (and a counter for the number of turns)
    if(stepType == 2) {
      // Determine if a turn was made
      if(stepIter >= contKeep && pow(fabs(X_Cont/X_Store[(stepIter-1)%contKeep][M*(numComponents+1)]),pow(-1,turnCount+1)) <= 1.0  &&  ds0 > 0)  turnCount++;
      if(stepIter >= contKeep && pow(fabs(X_Cont/X_Store[(stepIter-1)%contKeep][M*(numComponents+1)]),pow(-1,turnCount  )) <= 1.0  &&  ds0 < 0)  turnCount++;

      // Track the surface tension
      if(turnCount%2 == 0 && X_Cont > X_ContTrack) {
          gammaTrack  = gamma;
          X_ContTrack = X_Cont;
          stopContinuation = false;
      }
      else {
          if(gamma/gammaTrack > 1.2) stopContinuation = true;
      }

      dataStream << std::scientific << std::setprecision(8) << s << " ";
      dataStream << std::scientific << std::setprecision(8) << turnCount << " ";
    }

    // End line and close
    dataStream << " \n";
    dataStream.close();
}

void SCFT::calcExcessAdsorption() {
    double exAds;
    for(int i = 0; i < numComponents; i++) {
        exAds = (phi[i][0] - phi[i][M-1])/2  + (phi[i][M-1] - phi[i][M-1])/2;
        for(int k = 1; k < M-1; k++) {
            exAds += (phi[i][k] - phi[i][M-1]);
        }
        exAds *= gridSize;

        // Store value to keep track of excess adsorption
        exAdsStore[i][stepIter%keepPhase] = exAds;
    }
}


void SCFT::checkPhaseTransition() {
    double threshJump = 5;
    double temp1, temp2, temp3,tempSlope1,tempSlope2;


    if(prePhaseTransition && stepIter > 2) {
        // Loop through the components to detect a jump in the
        for(int i = 0; i < numPolymers; i++) {
            // Temporarily store values
            temp1 = exAdsStore[i][stepIter%keepPhase];
            temp2 = exAdsStore[i][(stepIter-1)%keepPhase];
            temp3 = exAdsStore[i][(stepIter-2)%keepPhase];

            tempSlope1 = (temp1 - temp2)/temp2;
            tempSlope2 = (temp2 - temp3)/temp3;

            // Detect phase transition based on jump in excess adsorption (jump upward)
            if(fabs(tempSlope1/tempSlope2) > threshJump && fabs(tempSlope1) > threshJump/2 && temp1 > 1E-2 && stepType != 10 && stepType != 11 && stepType != 12 && stepType != 15 && stepType != 16) {
                // A phase transition has occured
                prePhaseTransition = false;
                writeLog("Phase transition (jump up) detected between steps " + std::to_string(stepIter - 1) + " and " + std::to_string(stepIter));

                // Set the step parameter so that the system ends where it started or stop entirely
                if(stepType < 5) {
                    if(stepIter + maxPhaseIter <= numSteps) numSteps = stepIter + maxPhaseIter;
                }
                else {
                    if(2*(stepIter + maxPhaseIter) <= numSteps) numSteps = 2*(stepIter + maxPhaseIter);
                }
                break;
            }
        }
    }

    for(int i = 0; i < numPolymers; i++) {
        // Temporarily store values
        temp1 = exAdsStore[i][stepIter%keepPhase];
        temp2 = exAdsStore[i][(stepIter-1)%keepPhase];
        temp3 = exAdsStore[i][(stepIter-2)%keepPhase];

        tempSlope1 = (temp1 - temp2)/temp2;
        tempSlope2 = (temp2 - temp3)/temp3;

        // Turn around if the excess adsorption is too large
        if(prePhaseTransition && temp1 > 1.0 && stepType != 2) {
            // A phase transition has occured
            prePhaseTransition = false;
            writeLog("Turning around or stopping at step " + std::to_string(stepIter) + " because of the excess adsorption is too high.");

            // Set the step parameter so that the system ends where it started
            if(stepType < 5 || stepType == 11) {
                numSteps = stepIter;
            }
            else {
                numSteps = 2*stepIter;
            }
            break;
        }
        else if(!prePhaseTransition && temp1 > 2.0 && stepType != 2) {
            // A phase transition has occured
            writeLog("Turning around or stopping at step " + std::to_string(stepIter) + " because of the excess adsorption is too high.");

            // Set the step parameter so that the system ends where it started
            if(stepType < 5 || stepType == 11) {
                numSteps = stepIter;
            }
            else {
                numSteps = 2*stepIter;
            }
            break;
        }

        // Continuation stopping condition
        if(stepType == 2 && ds0 > 0 && (fabs(temp1) > 1.5 || (temp1 > 1.0 && stopContinuation)) ) {
            writeLog("Stopping continuation at step " + std::to_string(stepIter) + " because of the excess adsorption exceeds the limit.");

            // Stop taking steps
            numSteps = stepIter;
            break;
        }
    }
}
