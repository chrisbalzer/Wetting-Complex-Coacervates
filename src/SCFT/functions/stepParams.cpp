#include "../SCFT.h"

void SCFT::stepParameters() {
    if(stepType == 0) {
        for (int i = 0; i < numPolymers; i++) {
            phiB[i] = pow(10, log10(phiB[i]) + stepSize);
        }

        // Calculate the chemical potential in bulk
        chemicalPotentialbulk();
    }
    if(stepType == 1) {
        if(stepIter <= numSteps/2) {
            for (int i = 0; i < numPolymers; i++) {
                phiB[i] = pow(10, log10(phiB[i]) + stepSize);
            }
        }
        else {
            for (int i = 0; i < numPolymers; i++) {
                phiB[i] = pow(10, log10(phiB[i]) - stepSize);
            }
        }
        
        // Calculate the chemical potential in bulk
        chemicalPotentialbulk();
    }
    if(stepType == 2) {
        // Index for the continuation parameter
        int lastIndx = M*(numComponents+1);

        // Other variables
        double a,b,c,ds1,ds2;
        double ratio = 1.0;
        double e1,e2;
        double maxAds0,maxAds1,maxAds2;

        // Store previous value and reset X0 to the last solution
        for (int i = 0; i < M*(numComponents+1); i++) {
            s_store[(stepIter - 1)%contKeep] = s;
            X_Store[(stepIter - 1)%contKeep][i] = X[i];
            X0[i] = X[i];
        }
        dYds_store[(stepIter - 1)%contKeep] = dYds[lastIndx];
        X_Store[(stepIter - 1)%contKeep][lastIndx] = X_Cont;
        X0[lastIndx] = X_Cont;

        // Take step
        if(stepIter < contKeep) {
            for (int i = 0; i < numPolymers; i++) {
                if(ds0 > 0) phiB[i] = pow(10, log10(phiB[i]) + 0.01);
                if(ds0 < 0) phiB[i] = pow(10, log10(phiB[i]) - 0.01);
            }
            X_Cont = log10(phiB[0]);
            X[lastIndx] = X_Cont;

            // Increment s
            s += 0.01;
        }
        else {
            // We only need the absolute value
            ds = fabs(ds);

            // Compute tangent
            // Nonuniform s spacing
            ds1 = s_store[(stepIter - 1)%contKeep] - s_store[(stepIter - 2)%contKeep];
            ds2 = s_store[(stepIter - 2)%contKeep] - s_store[(stepIter - 3)%contKeep];
            if(stepIter >= 2*contKeep) {
              c = (2*ds1 + ds2)/ds1/(ds1 + ds2);
              b = -(ds1 + ds2)/ds1/ds2;
              a = ds1/ds2/(ds1 + ds2);
            }
            else {
              c = 1.0/ds1;
              b = -1.0/ds1;
              a = 0.0;
            }

            double normdYds = 0.0;
            for (int i = 0; i < M*(numComponents+1) + 1; i++) {
                dYds[i]   = c*X_Store[(stepIter-1)%contKeep][i] + b*X_Store[(stepIter-2)%contKeep][i] + a*X_Store[(stepIter-3)%contKeep][i];
                normdYds +=  dYds[i]*dYds[i];
            }
            normdYds = sqrt(normdYds);

            // Normalize vector
            for (int i = 0; i < M*(numComponents+1) + 1; i++) {
                dYds[i] /= normdYds;
            }

            // Step size limits
            double maxDs = 0.05;
            double minDs = 0.01;

            // Calculate the last few adsorption
            maxAds0 = 0.0;
            maxAds1 = 0.0;
            maxAds2 = 0.0;
            for(int i = 0; i < numPolymers; i++) {
                if(fabs(exAdsStore[i][(stepIter-1)%keepPhase]) > maxAds0) maxAds0 = fabs(exAdsStore[i][(stepIter-1)%keepPhase]);
                if(fabs(exAdsStore[i][(stepIter-2)%keepPhase]) > maxAds1) maxAds1 = fabs(exAdsStore[i][(stepIter-2)%keepPhase]);
                if(fabs(exAdsStore[i][(stepIter-3)%keepPhase]) > maxAds2) maxAds2 = fabs(exAdsStore[i][(stepIter-3)%keepPhase]);
            }

            // Change "ds" based on the adsorption trends
            if(stepIter > 2*contKeep) {
              // Base on polymer adsorption
              e1  = (maxAds0 - maxAds1)/(X_Store[(stepIter - 1)%contKeep][lastIndx] - X_Store[(stepIter - 2)%contKeep][lastIndx]);
              e2  = (maxAds1 - maxAds2)/(X_Store[(stepIter - 2)%contKeep][lastIndx] - X_Store[(stepIter - 3)%contKeep][lastIndx]);
              ratio = e2/e1;
              if(ratio > 1.05)      ds *= 1.05;
              else if(ratio > 1.00) ds *= 1.02;
              else if(ratio < 0.85) ds *= 0.85;
            }
            if(ds < minDs) ds = minDs;
            if(ds > maxDs) ds = maxDs;
            if(maxAds0 < 1e-4)            ds = maxDs;
            if(maxAds0 > 0.3 && ds0 > 0)  ds = 2*maxDs;

            // Step parameters
            for (int i = 0; i < M*(numComponents + 1); i++) {
                X[i] = X0[i] + dYds[i]*ds;
            }
            X_Cont = X0[lastIndx] + dYds[lastIndx]*ds;

            // Store continuation parameter
            for(int i = 0; i < numPolymers; i++) {
                 phiB[i] = pow(10, X_Cont);
            }

            // Increment s
            s += ds;
        }

        // Output for testing
        // std::cout<<"----Continuation (Polymer Density)-----"<<std::endl;
        // std::cout<<" step      = "<<stepIter<<std::endl;
        // std::cout<<" s         = "<<s<<std::endl;
        // std::cout<<" ds        = "<<ds<<std::endl;
        // std::cout<<" dYds      = "<<dYds[lastIndx]<<std::endl;
        // if(stepIter > 2*contKeep) {
        //   std::cout<<" ratio     = "<<ratio<<std::endl;
        //   std::cout<<" e1        = "<<e1<<std::endl;
        //   std::cout<<" e2        = "<<e2<<std::endl;
        // }
        // std::cout<<" c_param   = "<<pow(10,X_Cont)<<std::endl;
        // std::cout<<" turnCount   "<<turnCount<<std::endl;
        // if(stepIter > 2) std::cout<<" gammaTrack  "<<gammaTrack<<std::endl;
    }
}
