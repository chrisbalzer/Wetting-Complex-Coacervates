/* ===================================================================
Wetting of Complex Coacervates on Single Surface

This software is distributed under the MIT License.

Copyright (c) 2022 Chris Balzer
    Based on Balzer, Zhang, and Wang, Soft Matter, 18, 6326-6339 (2022) (https://doi.org/10.1039/D2SM00859A)
    See README.md in repository for more info and citation ---> https://github.com/chrisbalzer/Wetting-Complex-Coacervates
======================================================================*/

#include "global_external.h"
#include "global_internal.h"

using namespace SysInfo;
using namespace std;

int main(int argc, char* argv[]) {
    // Parse input arguments
    string inputPath;
    if(argc > 1) {
        inputPath = argv[1];
    }
    else {
        inputPath = "./";
    }

    // Initialize SCFT engine
    SCFT scft = SCFT(inputPath);

    // Initialize timers
    Timers MainTimer;
    Timers IterationTimer;

    // Start timer for main routine
    START_CPU_TIMER(MainTimer);

    // Loop over all conditions
    int step_count = 0;
    for (int i = 0; i <= scft.numSteps; i++) {
            // Start iteration timer
            START_CPU_TIMER(IterationTimer);

            // Solve for structure
            try {
                scft.solve();
            }
            catch (int e) {
                scft.writeLog("Error in solveDensity(): "+ std::to_string(e) + ", Step Number - " + std::to_string(step_count) + "\n");
                if(e == 10) break;
                scft.initializeDensity();
            }

            // Step the parameter (see input file)
            scft.stepParameters();

            // Stop iteration timer
            STOP_CPU_TIMER(IterationTimer);

            // Increase the total number of steps
            step_count++;
    }
    // Average iteration time
    IterationTimer.accum /= step_count;

    // Stop timer for main routine
    STOP_CPU_TIMER(MainTimer);

    // Write timing to output
    scft.writeLog("\nTotal Program Time         : " + std::to_string(MainTimer.duration/60) + " minutes");
    scft.writeLog("Average Time per Iteration : " + std::to_string(IterationTimer.accum/60) + " minutes");

    // Finalize I/O
    scft.finalizeIO();

    // Print the total time to the terminal
    std::cout << "\nTotal Program Time : " << MainTimer.duration/60 << " minutes" << std::endl;

    return 0;
}
