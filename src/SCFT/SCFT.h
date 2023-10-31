/* ===================================================================
Copyright (c) 2022 Chris Balzer
    Based on Balzer, Zhang, and Wang, Soft Matter, 18, 6326-6339 (2022) (https://doi.org/10.1039/D2SM00859A)
    See README.md in repository for more info and citation ---> https://github.com/chrisbalzer/Wetting-Complex-Coacervates
======================================================================*/

#ifndef SCFT_H
#define SCFT_H

/* Constants */
#define epsilon0 8.854187817E-12L
#define kB 1.3806504E-23L
#define e0 1.60217646E-19L
#define Na 6.02214179E+23L
#define Pi 3.1415926536

#include "../global_external.h"

/* Include utilities */
#include "../utilities/utilities.h"

/* Eigen dependencies */
#include <Eigen/Dense>

class SCFT {
private:
    /*=============================
    SCFT Variables
    =============================*/
    /* System variables */
    int    M;             // Total number of grid points
    int    numPolymers;   // Total number of polymer components
    int    numComponents; // Total number of components in the system
    double lunit;         // Length unit
    double systemSize;    // Starting separation (length units)
    double gridSize;      // Grid sizing (length units)
    double BJ;            // Bjerrum length (length units)

    /* Iterative variables */
    int    minIters;
    int    maxIters;
    double errorTol;
    double errorTol0;
    double error;
    double error_prev;

    /* Number of processors */
    short printData;

    /* Number of processors */
    int numProcessors;

    /* Tracking phase transition */
    int  phaseIter;
    int  maxPhaseIter;
    int  keepPhase;
    bool prePhaseTransition;
    double**  exAdsStore;

    /* Boundary variables */
    short  BC_type;       // Type of boundary condition
    double surfBC;        // Surface boundary condition

    /* Arrays for each component */
    short  *Z;            // Valency for each component
    double *D;            // Length scale for each component (length unit)
    double *N;            // Chain length for each component
    double *mu;           // Total chemcial potential for each component
    double *muEx;         // Excess chemical potential for each component
    double *phiB;         // Bulk volume fraction for each component
    double *uSurf;        // Interaction parameter with each surface (eta)

    /* Spatially defined arrays */
    double *z;
    double *Psi;
    double *rhoZ;
    double **phi;
    double **mu_local;
    double **mu_localEx;
    double **varphi;
    double **phiNew;

    /* Mixing parameters */
    int     mixIter;
    int     mixKeep;
    short   mixingType;
    double  mixF;
    double  minMix;
    double  maxMix;
    double* prevMix;

    /* Anderson-Acceleration parameters (solving ODEs) */
    int       anderM;
    int       anderIter;

    /* Mixing variables (relaxation) */
    double *F;
    double *X;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X_Relax;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> F_Relax;
    int      broydenIter;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X_Broyden;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> F_Broyden;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> J_Broyden;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> G_And;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X_And;

    /* Continuation variables */
    int    contKeep;
    int    turnCount;
    bool   stopContinuation;
    double s;
    double ds;
    double ds0;
    double X_Cont;
    double F_Cont;
    double gammaTrack;
    double X_ContTrack;
    double *dYds_store;
    double *s_store;
    double *X0;
    double *dYds;
    double **X_Store;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> G_And_Cont;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X_And_Cont;

    /*=============================
      Data I/O
    =============================*/
    /* Directories */
    std::string cwd;         // current working directory
    std::string outputDir;   // base level directory for data output (relative to cwd)

    /* File names */
    std::string logFilename   = "log.out";
    std::string iterFilename  = "iterations.out";
    std::string dataFilename  = "data.dat";


    /* Actual Files */
    std::string inputFile;
    std::string logFile;
    std::string iterFile;
    std::string outFile;
    std::string dataFile;

    /* std::ofstream to do the actual I/O */
    std::ifstream inputStream;
    std::ofstream logStream;
    std::ofstream iterStream;
    std::ofstream outStream;
    std::ofstream dataStream;

    void setOutputFileStructure();
    void closeFiles();
    void readInput(std::string inputPath);

    /*=============================
      Nuumerical tools
    =============================*/
    double trapz(double *f, int lower_bound, int upper_bound);
    double trapz(double* z, double *f, int lower_bound, int upper_bound);
    double gaussQuad(double* z, double *f, int lower_bound, int upper_bound);
    double simps1(double* f, int lower_bound, int upper_bound);
    double simps38(double* f, int lower_bound, int upper_bound);

public:
    /* Public variables*/
    /* Parameter stepping */
    int    numSteps;
    int    stepIter;
    short  stepType;
    double stepSize;

    /* Default constructor/destructor*/
    SCFT();
    ~SCFT();

    /* Constructor */
    SCFT(std::string inputPath);

    /*=============================
      Main Engine
    =============================*/
    void initializeSystem();
    void initializeDensity();
    void solve();
    void solveDensity();
    void solveDensityContDensity();
    void finalizeSystem();

    /*=============================
      Update density and mixing
    =============================*/
    // Relaxation type of update
    void updateDensityRelax();
    void mixDensityRelax();
    void mixDensityPicardRelax();
    void mixDensityAndersonRelax();
    void mixDensityAndersonRelaxCont();
    void calcFValues(double* F_temp, double* X_temp);
    void calcFValuesHO(double* F_temp, double* X_temp);

    /*=============================
      Energy
    =============================*/
    void calcEnergy(double *phi_k, double& f);
    void IdealEnergy(double *phi_k, double& f_k);
    void DHEnergy(double *phi_k, double& f_k);
    void EVEnergy(double *phi_k, double& f_k);
    void DHChainEnergy(double *phi_k, double& f_k);

    /*=============================
      Potentials
    =============================*/
    void chemicalPotentialbulk();
    void chemicalPotentiallocal(double** mu_local_jk, double** mu_localEx_jk, double**phi_jk);
    void chemicalPotentiallocalParallel(double** mu_local_jk, double** mu_localEx_jk, double**phi_jk);
    void IdealPotential(double *phi_k, double* mu_k);
    void IdealPotential(double **phi_jk, double** mu_jk);
    void DHPotential(double *phi_k, double* mu_k);
    void DHPotential(double **phi_jk, double** mu_jk);
    void EVPotential(double *phi_k, double* mu_k);
    void EVPotential(double **phi_jk, double** mu_jk);
    void DHChainPotential(double *phi_k, double* mu_k);
    void DHChainPotential(double **phi_jk, double** mu_jk);

    /*=============================
      Other Calculations
    =============================*/
    // Properties
    void calcChargeDensity(double* rhoZ_k, double** phi_jk);
    void DHScreeningLength(double* phi_k,double& kappa);
    void calcProps();
    void calcExcessAdsorption();
    void checkPhaseTransition();

    // Error
    void calculateErrorRelax(double* F_ik, double* X_ik);

    /*=============================
      Data operations
    =============================*/
    void densityNAN();
    void stepParameters();
    void refillX();
    void getX();

    /*=============================
      Data I/O
    =============================*/
    void writeInitial();
    void writeLog(std::string logText);
    void openIter();
    void writeIter(int iter);
    void writeDensPot();
    void redesignateFiles();
    void finalizeIO();
};

#endif //SCFT_H
