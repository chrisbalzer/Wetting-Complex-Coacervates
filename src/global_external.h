/* ===================================================================
Copyright (c) 2022 Chris Balzer
    Based on Balzer, Zhang, and Wang, Soft Matter, 18, 6326-6339 (2022) (https://doi.org/10.1039/D2SM00859A)
    See README.md in repository for more info and citation ---> https://github.com/chrisbalzer/Wetting-Complex-Coacervates
======================================================================*/

#ifndef GLOBALEXTERNAL_H
#define GLOBALEXTERNAL_H

/* IO dependencies */
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>

/* std::string dependencies */
#include <string>
#include <utility>
#include <iomanip>

/* Access system parameters */
#include <chrono>
#include <ctime>
#include <thread>

/* Math dedendencies */
#include <cmath>
#include <vector>
#include <random>
#include <functional>
#include <math.h>

/* Add OpenMP*/
#ifdef _OPENMP
    #include <omp.h>
#endif

#endif // GLOBALEXTERNAL_H
