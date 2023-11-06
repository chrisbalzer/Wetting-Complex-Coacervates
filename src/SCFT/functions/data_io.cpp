#include "../SCFT.h"

void SCFT::readInput(std::string inputPath) {
    // Initialize some variables
    int item;
    std::string s;
    std::string delim = ";";
    std::string sub_s;
    std::size_t pos = 0;

    // Set the input file name
    inputFile = inputPath + "/input.dat";

    // Open the file
    inputStream.open(inputFile.c_str(), std::fstream::in);

    // Check if file has been opened
    if(inputStream) {
        // Read through input file
        while(std::getline(inputStream,s)) {

            if(s.compare("POLYMER_FRACTION:")==0) {
                if(std::getline(inputStream,s)) {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos) {
                        sub_s = s.substr(0, pos);

                        for (int j = 0; j < numPolymers; j++) {
                            if (item == j) phiB[j] = atof(sub_s.c_str());
                        }

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                }
                else {
                    std::cerr<<"SCFT::read_input - Could not read in POLYMER_FRACTION"<<std::endl;
                    std::exit(0);
                }
            }

            if(s.compare("SALT_FRACTION:")==0) {
                if(std::getline(inputStream,s)) {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos) {
                        sub_s = s.substr(0, pos);


                        for (int j = numPolymers; j < numComponents; j++) {
                            if (item == j - numPolymers) phiB[j] = atof(sub_s.c_str());
                        }

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                }
                else {
                    std::cerr<<"SCFT::read_input - Could not read in SALT_FRACTION"<<std::endl;
                    std::exit(0);
                }
            }

            if(s.compare("CHAIN_LENGTH:")==0) {
              if(std::getline(inputStream,s)) {
                  item = 0;
                  pos = 0;
                  while ((pos = s.find(delim)) != std::string::npos) {
                      sub_s = s.substr(0, pos);

                      for (int j = 0; j < numPolymers; j++) {
                          if (item == j) N[j] = atof(sub_s.c_str());
                      }

                      s.erase(0, pos + delim.length());
                      item++;
                  }
              }
              else {
                  std::cerr<<"SCFT::read_input - Could not read in CHAIN_LENGTH"<<std::endl;
                  std::exit(0);
              }
            }


            if(s.compare("SURFACE:")==0) {
                if(std::getline(inputStream,s)) {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos) {
                        sub_s = s.substr(0, pos);

                        if (item == 0) BC_type     = atof(sub_s.c_str());
                        if (item == 1) surfBC      = atof(sub_s.c_str());
                        for (int j = 0; j < numComponents; j++) {
                            if (item == 2+j) uSurf[j] = atof(sub_s.c_str());
                        }

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                }
                else {
                    std::cerr<<"SCFT::read_input - Could not read in SURFACE"<<std::endl;
                    std::exit(0);
                }
            }


            if(s.compare("VALENCY:")==0) {
                if(std::getline(inputStream,s)) {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos) {
                        sub_s = s.substr(0, pos);

                        for (int j = 0; j < numComponents; j++) {
                            if (item == j) Z[j] = atof(sub_s.c_str());
                        }

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                  }
                  else {
                      std::cerr<<"SCFT::read_input - Could not read in VALENCY"<<std::endl;
                      std::exit(0);
                  }
            }

            if(s.compare("DIAMETER:")==0) {
                if(std::getline(inputStream,s)) {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos) {
                        sub_s = s.substr(0, pos);

                        for (int j = 0; j < numComponents; j++) {
                            if (item == j) D[j] = atof(sub_s.c_str());
                        }

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                  }
                  else {
                      std::cerr<<"SCFT::read_input - Could not read in DIAMETER"<<std::endl;
                      std::exit(0);
                  }
            }

            if(s.compare("SYSTEM:")==0) {
                if(std::getline(inputStream,s)) {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos) {
                        sub_s = s.substr(0, pos);

                        if (item == 0) systemSize = atof(sub_s.c_str());
                        if (item == 1) gridSize = atof(sub_s.c_str());
                        if (item == 2) BJ = atof(sub_s.c_str());

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                  }
                  else {
                      std::cerr<<"SCFT::read_input - Could not read in SYSTEM"<<std::endl;
                      std::exit(0);
                  }
            }

            if(s.compare("ITERATIVE:")==0) {
                if(std::getline(inputStream,s)) {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos) {
                        sub_s = s.substr(0, pos);
                        if (item == 0) maxIters    = atof(sub_s.c_str());
                        if (item == 1) errorTol0   = atof(sub_s.c_str());
                        if (item == 2) mixingType  = atof(sub_s.c_str());
                        if (item == 3) mixF        = atof(sub_s.c_str());

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                  }
                  else {
                      std::cerr<<"SCFT::read_input - Could not read in ITERATIVE"<<std::endl;
                      std::exit(0);
                  }
            }

            if(s.compare("STEPS:")==0) {
                if(std::getline(inputStream,s)) {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos) {
                        sub_s = s.substr(0, pos);
                        if (item == 0) stepType    = atof(sub_s.c_str());
                        if (item == 1) numSteps    = atof(sub_s.c_str());
                        if (item == 2) stepSize    = atof(sub_s.c_str());


                        s.erase(0, pos + delim.length());
                        item++;
                    }
                  }
                  else {
                      std::cerr<<"SCFT::read_input - Could not read in STEPS"<<std::endl;
                      std::exit(0);
                  }
            }

            if(s.compare("CPUs:")==0) {
                if(std::getline(inputStream,s)) {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos) {
                        sub_s = s.substr(0, pos);
                        if (item == 0) numProcessors  = atof(sub_s.c_str());

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                  }
                  else {
                      std::cerr<<"SCFT::read_input - Could not read in CPUs"<<std::endl;
                      std::exit(0);
                  }
            }

            if(s.compare("OUT_DATA:")==0) {
                if(std::getline(inputStream,s)) {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos) {
                        sub_s = s.substr(0, pos);
                        if (item == 0) printData  = atof(sub_s.c_str());

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                  }
                  else {
                      std::cerr<<"SCFT::read_input - Could not read in OUT_DATA"<<std::endl;
                      std::exit(0);
                  }
            }


            if(s.compare("OUT_PATH:")==0) {
                if(std::getline(inputStream,s)) {
                    // Store the value
                    outputDir = outputDir  +  s + "/";
                }
                else {
                    std::cerr<<"SCFT::read_input - Could not read in OUT_PATH"<<std::endl;
                    std::exit(0);
                }
            }
      }

        // Close the input file
        inputStream.close();
    }
    else {
        std::cout<<"SCFT::read_input - No input file found. Using default."<<std::endl;
    }

    // Set up the file outputs
    setOutputFileStructure();
}

void SCFT::writeInitial() {
    logStream.open(logFile.c_str(), std::fstream::out | std::ofstream::app);
    // Write initial system
    logStream<<"=======================Wetting of Complex Coacervates========================="<<" \n";
    logStream<<"                  Copyright Owners: Chris Balzer (2022)"<<" \n";
    logStream<<"=============================================================================="<<" \n";

    logStream<<"\n==============="<<" \n";
    logStream<<"INITIALIZATION"<<" \n";
    logStream<<"==============="<<" \n";
    logStream<<"Data read from input file   : "<<inputFile<<" \n";
    logStream<<"Writing output to directory : "<<outputDir<<" \n";

    /* System parameters from input file */

    logStream<<"    Bulk Volume Fraction    : ";
    for (int j = 0; j < numComponents; j++) {
        logStream<<phiB[j]<<"; ";
    }
    logStream<<" \n";

    logStream<<"    Chain Length            : ";
    for (int j = 0; j < numPolymers; j++) {
        logStream<<N[j]<<"; ";
    }
    logStream<<" \n";

    logStream<<"    Valence                 : ";
    for (int j = 0; j < numComponents; j++) {
        logStream<<Z[j]<<"; ";
    }
    logStream<<" \n";

    logStream<<"    Diameter (Pos : Neg)    : ";
    for (int j = 0; j < numComponents; j++) {
        logStream<<D[j]<<"; ";
    }
    logStream<<" \n";

    logStream<<"    Surface Interaction     : ";
    for (int j = 0; j < numPolymers; j++) {
        logStream<<uSurf[j]<<"; ";
    }
    logStream<<" \n";


    if(BC_type == 0) {
        logStream<<"    Surface Type            : "<<"Fixed Potential"<<" \n";
        logStream<<"        Surface BC          : "<<surfBC<<" beta*e*V"<<" \n";
    }
    else {
        logStream<<"    Surface Type            : "<<"Fixed Surface Charge"<<" \n";
        logStream<<"        Surface BC          : "<<surfBC<<" e/sigma^2"<<" \n";
    }

    logStream<<"    System Size             : "<<systemSize<<" \n";
    logStream<<"    Grid Sizing             : "<<gridSize<<" \n";
    logStream<<"    Bjerrum Length          : "<<BJ<<" \n";
    logStream<<"    Error Tolerance         : "<<errorTol0<<" \n";

    if(stepType == 0) {
        logStream<<"    Step Type               : "<<"Polycation/polyanion bulk density (log)"<<" \n";
    }
    else if(stepType == 1) {
        logStream<<"    Step Type               : "<<"Polycation/polyanion bulk density (log) (Out and Back)"<<" \n";
        numSteps *= 2;
    }
    else if(stepType == 2) {
        logStream<<"    Step Type               : "<<"Continuation - Polymer fraction"<<" \n";

        // Switch to Anderson Acceleration
        mixingType = 1;
    }

    if(mixingType == 0) {
        logStream<<"    Mixing                  : "<<"Picard"<<" \n";
    }
    else if(mixingType == 1) {
        logStream<<"    Mixing                  : "<<"Anderson Acceleration"<<" \n";;
    }

    logStream<<"    Number of Steps         : "<<numSteps<<" \n";
    logStream<<"    Step Size               : "<<stepSize<<" \n";

    logStream<<"    Number of CPUs          : "<<numProcessors<<" \n";
    if (numProcessors == 1) {
          logStream<<"        Only 1 processor - not using OpenMP \n";
    }
    if (numProcessors > 1) {
        #ifdef _OPENMP
            logStream<<"        Using OpenMP \n";
        #else
            logStream<<"        Cannot find OpenMP (switching to 1 processor) \n";
            numProcessors = 1;
        #endif
    }

    logStream<<"\n==============="<<" \n";
    logStream<<"RUNNING SCFT"<<" \n";
    logStream<<"==============="<<" \n";

    // Close
    logStream.close();
}


void SCFT::writeLog(std::string logText) {
    logStream.open(logFile.c_str(), std::fstream::out | std::ofstream::app);
    logStream<<logText<<" \n";
    logStream.close();
}


void SCFT::openIter() {
    iterStream.open(iterFile.c_str(), std::fstream::out | std::ofstream::trunc);
    iterStream.close();
}

void SCFT::writeIter(int iter) {
    iterStream.open(iterFile.c_str(), std::fstream::out | std::ofstream::app);
    iterStream << "Iteration: " << iter << "  ";
    iterStream << "Error: " <<  std::scientific << std::setprecision(4) << error << "  ";
    iterStream << "Mix Coeff: " << std::scientific << std::setprecision(4) << mixF <<" \n";
    iterStream.close();
}


// Write the density profiles and electrostatic potential
void SCFT::writeDensPot() {
    // Get file name and open file
    outFile = outputDir + "density_potential_" + std::to_string(stepIter) + ".dat";
    outStream.open(outFile.c_str(), std::fstream::out | std::ofstream::trunc);

    // Print to file
    for ( int j = 0; j < M; j++) {
        outStream << std::fixed << std::setprecision(5) << z[j] << " ";
        for (int i = 0; i < numComponents; i++) {
            outStream << std::scientific << std::setprecision(8) << phi[i][j] << " ";
        }
        outStream << std::fixed << std::setprecision(8) << Psi[j] << " \n";
    }

    // Close the file
    outStream.close();
}

// Rename the density profiles and iteration files if the run fails
void SCFT::redesignateFiles() {
    std::string newFile;

    try {
        outFile = outputDir + "error_" + std::to_string(stepIter) + ".dat";
        newFile = outputDir + "UNCONVERGED_error_" + std::to_string(stepIter) + ".dat";
        std::rename(outFile.c_str(), newFile.c_str());

        outFile = outputDir + "mu_" + std::to_string(stepIter) + ".dat";
        newFile = outputDir + "UNCONVERGED_mu_" + std::to_string(stepIter) + ".dat";
        std::rename(outFile.c_str(), newFile.c_str());

        outFile = outputDir + "density_potential_" + std::to_string(stepIter) + ".dat";
        newFile = outputDir + "UNCONVERGED_density_potential_" + std::to_string(stepIter) + ".dat";
        std::rename(outFile.c_str(), newFile.c_str());

        // Keep the iteration error file
        newFile = outputDir + "iterations_" + std::to_string(stepIter) +  ".dat";
        std::rename(iterFile.c_str(), newFile.c_str());
    }
    catch(const char* e) {
        std::cerr<<e<<std::endl;
    }
}


void SCFT::finalizeIO() {
    // Write a message to end
    logStream<<"\n==============="<<" \n";
    logStream<<"END"<<" \n";
    logStream<<"==============="<<" \n";

    // Close all appending files
    closeFiles();
}

void SCFT::closeFiles() {
    logStream.close();
}

void SCFT::setOutputFileStructure() {
    // Set file streams
    logFile      = outputDir + logFilename;
    iterFile     = outputDir + iterFilename;
    dataFile     = outputDir + dataFilename;

    // Make the output directory if it isn't already made
    struct stat sb;
    if(stat(outputDir.c_str(), &sb) == -1) {
        std::string mkdirBase = "mkdir " + outputDir;
        int systemRet = std::system(mkdirBase.c_str());
        if(systemRet == -1)
        {
          std::cerr<<"SCFT : failed to create unique output directory"<<std::endl;
          exit(0);
        }
    }

    // Open files that will be open the whole time
    logStream.open(logFile.c_str(), std::fstream::out | std::ofstream::trunc);
    dataStream.open(dataFile.c_str(), std::fstream::out | std::ofstream::trunc);
    
    // Add data headers
    dataStream << "# (0) Surf. Charge, (1) Surf. Pot., (2-3) Eta Poly., (4-7) Excess Adsorption, (8) Helmholtz, (9) Grand Potential, (10) Osmotic pressure, (11) Surface Tension, (12-15) Bulk Density, (16-19) Chemical Potential";
    if(stepType == 2) dataStream << ", (20) s, (21) turnCount\n";
    else dataStream << "\n";
    
    logStream.close();
    dataStream.close();
}
