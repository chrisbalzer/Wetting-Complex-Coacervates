#include "../SCFT.h"

/* Controller for solver */
 void SCFT::solve() {
    if(stepType == 2)      solveDensityContDensity();
    else                   solveDensity();
 }

 /* Main engine for SCFT */
 void SCFT::solveDensity() {
     // Initialize variables
     int    iter = 0;
     SysInfo::Timers timer;

     // Start timer
     START_CPU_TIMER(timer);

     // Open any necessary files for the given run
     openIter();

     // Initialize error to large value
     error = 100.0;
     mixIter = 0;
     anderIter = 0;

     // Set the error tolerance (possibly change based on the condition)
     errorTol = errorTol0;

     // Update density until convergence
     do {
         // Update new density, calculate error, and mix densities
         updateDensityRelax();

         // Output iteration info to iterative file
         if(iter < 1E5 && iter%100 == 0){
              writeIter(iter);
         }
         else if(iter < 5E5 && iter%1000 == 0){
              writeIter(iter);
         }
         else if(iter%10000 == 0) {
              writeIter(iter);
         }

         // Print output
         if(iter%100000 == 0 && printData == 1){
            writeDensPot();
         }

         // Increase iteration number
         iter++;
     }
     while(iter < maxIters && (error > errorTol || iter < minIters ) );

     // Stop timer
     STOP_CPU_TIMER(timer);

     // Process data and print results
     if(error < errorTol) {
        // Calculate properties from current profiles
        calcProps();

        // Print complete
        writeLog("Complete - Step " + std::to_string(stepIter) +  ":  duration (min) = " + std::to_string(timer.duration/60) + ", iters = " + std::to_string(iter) );

        // Write the density and potential
        if (printData == 1) writeDensPot();
     }
     else if(error < 1e5*errorTol || error < errorTol0) {
        // Calculate properties from current profiles
        calcProps();

        // Print complete
        writeLog("Nearly complete - Step " + std::to_string(stepIter) +  ":  duration (min) = " + std::to_string(timer.duration/60) + ", iters = " + std::to_string(iter) );

        // Write the density and potential
        if (printData == 1) writeDensPot();
     }
     else {
        // Write message that the step didn't converge
        writeLog("Incomplete - Step " + std::to_string(stepIter) +  ":  duration (min) = " + std::to_string(timer.duration/60) + ", iters = " + std::to_string(iter) );

        // Write the density and potential
        if (printData == 1) {
            writeDensPot();
            // Rename the density profiles and indicate they didn't converge
            redesignateFiles();
        }
     }

     // Increase iteration number
     stepIter++;
 }

 /* Continuation of density*/
 void SCFT::solveDensityContDensity() {
     // Initialize variables
     double temp;
     int    iter = 0;
     SysInfo::Timers timer;

     // Start timer
     START_CPU_TIMER(timer);

     // Open any necessary files for the given run
     openIter();

     // Initialize error to large value
     error = 100.0;
     mixIter = 0;
     anderIter = 0;

     // Set the error tolerance (possibly change based on the condition)
     errorTol = errorTol0;

     // Index for last entry
     int lastIndx = M*numComponents + M;

     // Store the initial profile
     getX();

     // Update density until convergence
     do {
         // Calculate bulk properties
         chemicalPotentialbulk();

         // Calculate residual based on this guess
         calcFValuesHO(F,X);

         // Calcualte error
         calculateErrorRelax(F,X);

         // Calculate last entry in F (new equation for continuation)
         X_Cont = log10(phiB[0]);
         if(stepIter >= contKeep) {
             temp = 0;
             for(int i = 0; i < M*numComponents+M; i++) {
                  temp += dYds[i]*(X[i] - X0[i]);
             }
             temp += dYds[lastIndx]*(X_Cont - X0[lastIndx]) - ds;
             F_Cont = -temp;
             if(fabs(temp) > error) error = fabs(temp);
         }

         // Update density
         mixDensityRelax();

         // NaN handling
         if(stepIter >= contKeep) {
             for (int i = 0; i < numComponents+1; i++) {
                 for(int k = 0; k < M; k++) {
                     if(isnan(X[i*M+k]) || error > 1E10) {
                         // Adjust the step size
                         s  -= fabs(ds);
                         ds *= 0.5;
                         s  += fabs(ds);

                         // Write the iteration number
                         writeIter(iter);

                         // Stop continuation if ds is too small
                         if(fabs(ds) < 1e-6) {
                           std::cout<<"ds is too small"<<std::endl;
                           throw 10;
                         }

                         // Take smaller step forward
                         for (int j = 0; j < M*(numComponents + 1); j++) {
                             X[j] = X0[j] + dYds[j]*ds;
                         }
                         X_Cont = X0[lastIndx] + dYds[lastIndx]*ds;

                         // Break loop
                         i = numComponents+1;
                         break;
                     }
                 }
             }
         }

         // Store new continuation parameter
         for(int i = 0; i < numPolymers; i++) {
              phiB[i] = pow(10, X_Cont);
         }

         // Output iteration info to iterative file
         if(iter < 1E5 && iter%100 == 0){
              writeIter(iter);
         }
         else if(iter < 5E5 && iter%1000 == 0){
              writeIter(iter);
         }
         else if(iter%10000 == 0) {
              writeIter(iter);
         }

         // Print output
         if(iter%100000 == 0 && printData == 1){
            writeDensPot();
         }

         // Increase iteration number
         iter++;
     }
     while(iter < maxIters && (error > errorTol || iter < minIters || (stepIter == 0 && (error > errorTol/1e3 || iter < 75*minIters)) ) );

     // Convert the X array to density and potential
     getX();

     // Calculate the charge density
     calcChargeDensity(rhoZ, phi);

     // Stop timer
     STOP_CPU_TIMER(timer);

     // Process data and print results
     if(error < errorTol) {
        // Calculate properties from current profiles
        calcProps();

        // Print complete
        writeLog("Complete - Step " + std::to_string(stepIter) +  ":  duration (min) = " + std::to_string(timer.duration/60) + ", iters = " + std::to_string(iter) );

        // Write the density and potential
        if (printData == 1) {
            writeDensPot();
        }
     }
     else if(error < 1e5*errorTol || error < errorTol0) {
        // Calculate properties from current profiles
        calcProps();

        // Print complete
        writeLog("Nearly complete - Step " + std::to_string(stepIter) +  ":  duration (min) = " + std::to_string(timer.duration/60) + ", iters = " + std::to_string(iter) );

        // Write the density and potential
        if (printData == 1) writeDensPot();
     }
     else {
        // Write message that the step didn't converge
        writeLog("Incomplete - Step " + std::to_string(stepIter) +  ":  duration (min) = " + std::to_string(timer.duration/60) + ", iters = " + std::to_string(iter) );

        // Write the density and potential
        if (printData == 1) writeDensPot();
     }

     // Increase iteration number
     stepIter++;
 }
