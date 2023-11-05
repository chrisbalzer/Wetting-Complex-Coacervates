/* ===================================================================
Copyright (c) 2022 Chris Balzer
    Based on Balzer, Zhang, and Wang, Soft Matter, 18, 6326-6339 (2022) (https://doi.org/10.1039/D2SM00859A)
    See README.md in repository for more info and citation ---> https://github.com/chrisbalzer/Wetting-Complex-Coacervates
======================================================================*/

#include "SCFT.h"

/* Default constructor */
SCFT::SCFT() = default;

/* Default destructor */
SCFT::~SCFT() {
    finalizeSystem();

    delete [] Z;
    delete [] D;
    delete [] phiB;
    delete [] mu;
    delete [] muEx;
    delete [] N;
    delete [] uSurf;
}

/* Constructor */
SCFT::SCFT(std::string inputPath) {
    // Number of components
    numPolymers = 2;
    numComponents = numPolymers + 2;

    // Initialize variables
    Z      = new short[numComponents];
    D      = new double[numComponents];
    phiB   = new double[numComponents];
    mu     = new double[numComponents];
    muEx   = new double[numComponents];
    N      = new double[numPolymers];
    uSurf  = new double[numPolymers];

    // Initialize output directory
    outputDir = "./";

    // Initialize bulk densities
    for (int i = 0; i < numComponents; i++) {
        phiB[i] = 0.0;
    }

    //  Read input file and open files
    readInput(inputPath);

    // Check if the system is charge neutral
    double chargeNeutral = 0.0;
    for (int i = 0; i < numComponents; i++) {
        chargeNeutral += Z[i]*phiB[i];
    }
    if (chargeNeutral < 0) {
      std::cerr<<"SCFT::SCFT - Negative net charge in bulk. Adjusting cation concentration."<<std::endl;
      logStream<<"SCFT::SCFT - Negative net charge in bulk. Adjusting cation concentration."<<"\n";

      phiB[numComponents-2] += fabs(chargeNeutral/Z[numComponents-2]);
    }
    if (chargeNeutral > 0) {
      std::cerr<<"SCFT::SCFT - Positive net charge in bulk. Adjusting anion concentration."<<std::endl;
      logStream<<"SCFT::SCFT - Positive net charge in bulk. Adjusting anion concentration."<<"\n";

      phiB[numComponents-1] += fabs(chargeNeutral/Z[numComponents-1]);
    }

    // Check if the total bulk volume fraction is more than 1
    double totalFrac = 0.0;
    for (int i = 0; i < numComponents; i++) {
        totalFrac = phiB[i];
    }
    if (totalFrac > 1.0) {
      std::cerr<<"SCFT::SCFT - Total volume fraction greater than 1 in bulk."<<std::endl;
      exit(0);
    }

    // Write initial message with input parameters
    writeInitial();

    // Initialize system variables
    initializeSystem();

    // Initialize volume fraction of components
    initializeDensity();
 }

void SCFT::initializeSystem() {
    // Calculate gridpoints and sizing
    M = ceil(systemSize/gridSize) + 1;
    gridSize = systemSize / ((double) M - 1);

    // Set the number of processors
    #ifdef _OPENMP
        omp_set_num_threads(numProcessors);
    #endif

    // Other variables
    stepIter = 0;

    // Allocate memory for spatially dependent arrays
    z          = new double[M];
    Psi        = new double[M];
    rhoZ       = new double[M];
    phi        = new double*[numComponents];
    phiNew     = new double*[numComponents];
    mu_local   = new double*[numComponents];
    mu_localEx = new double*[numComponents];
    varphi     = new double*[numPolymers];

    for (int i = 0; i < numComponents; i++) {
        phi[i]        = new double[M];
        phiNew[i]     = new double[M];
        mu_local[i]   = new double[M];
        mu_localEx[i] = new double[M];
    }

    for (int i = 0; i < numPolymers; i++) {
        varphi[i]    = new double[M];
    }

    // Initialize values of some arrays
    for (int k = 0; k < M; k++) {
        Psi[k]  = 0.0;
        rhoZ[k] = 0.0;
        z[k]    = k*gridSize;
    }

    /*======================
      Variables to track phase transition
    ======================*/
    prePhaseTransition = true;
    phaseIter          = 0;
    maxPhaseIter       = 5;
    keepPhase          = 5;
    exAdsStore = new double*[numComponents];
    for(int j = 0; j < numComponents; j++) {
        exAdsStore[j] = new double[keepPhase];
    }

    /*======================
      Mixing variables
    ======================*/
    mixKeep = 10;
    prevMix = new double[mixKeep];

    // Anderson mixing variable
    anderM = 5;

    // Minimum iterations
    minIters = 5e3;

    /*======================
      Mixing - Relaxation
    ======================*/
    // Simple Mixing (relaxation method)
    X = new double[M*numComponents + M];
    F = new double[M*numComponents + M];

    // Mixing (relaxation method)
    if(mixingType == 1) {
        // Anderson Mixing (relaxation method)
        G_And.resize(M*numComponents + M, anderM);
        X_And.resize(M*numComponents + M, anderM);
    }

    /*======================
     Continuation
    ======================*/
    turnCount   = 0;
    X_ContTrack = -10000;
    stopContinuation = false;
    contKeep = 5;
    s        = 0;
    ds0      = stepSize;
    ds       = ds0;
    X_Store = new double*[contKeep];
    for (int i = 0; i < contKeep; i++) {
        X_Store[i] = new double[M*numComponents + M + 1];
    }
    s_store    = new double[contKeep];
    dYds_store = new double[contKeep];
    X0         = new double[M*numComponents + M + 1];
    dYds       = new double[M*numComponents + M + 1];
    if(stepType == 2) {
        // Anderson Mixing (relaxation method)
        G_And_Cont.resize(M*numComponents + M+1, anderM);
        X_And_Cont.resize(M*numComponents + M+1, anderM);
    }
}

void SCFT::finalizeSystem() {
    /* Delete variables */
    for (int i = 0; i < numComponents; i++) {
        delete [] phi[i];
        delete [] phiNew[i];
        delete [] mu_local[i];
        delete [] mu_localEx[i];
        delete [] exAdsStore[i];
    }
    delete [] phi;
    delete [] phiNew;
    delete [] mu_local;
    delete [] mu_localEx;
    delete [] exAdsStore;

    for (int i = 0; i < numPolymers; i++) {
        delete [] varphi[i];
    }
    delete [] varphi;


    delete [] z;
    delete [] Psi;
    delete [] rhoZ;
    delete [] prevMix;
    
    delete [] X;
    delete [] F;

    for (int i = 0; i < contKeep; i++) {
        delete [] X_Store[i];
    }
    delete [] X_Store;
    delete [] s_store;
    delete [] dYds_store;
    delete [] X0;
    delete [] dYds;
}
