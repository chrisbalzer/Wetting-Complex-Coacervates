#include "../SCFT.h"

void SCFT::mixDensityRelax() {
    // Select the mixing type
    if(mixingType == 0) {
        mixDensityPicardRelax();
    }
    else if(mixingType == 1) {
        if(stepType == 2 && stepIter >= contKeep) mixDensityAndersonRelaxCont();
        else mixDensityAndersonRelax();
    }
    
    // Set values
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

    // Store mixing value
    prevMix[mixIter%mixKeep] = mixF;

    // Increase iterative parameter
    mixIter++;
}
