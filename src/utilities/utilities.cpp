/* ===================================================================
Copyright (c) 2022 Chris Balzer
    Based on Balzer, Zhang, and Wang, Soft Matter, 18, 6326-6339 (2022) (https://doi.org/10.1039/D2SM00859A)
    See README.md in repository for more info and citation ---> https://github.com/chrisbalzer/Wetting-Complex-Coacervates
======================================================================*/

#include "utilities.h"

void GET_OMP_THREADS() {
  /* Set the number of OpenMP threads based on what is available */
  #ifdef _OPENMP
      int max_threads = omp_get_max_threads();
      if(max_threads > 12) omp_set_num_threads(6);
      else if(max_threads > 8) omp_set_num_threads(4);
      else if(max_threads > 4) omp_set_num_threads(2);
      else omp_set_num_threads(1);
  #endif
}

void START_CPU_TIMER(SysInfo::Timers &timers) {
  /* CPU timing */
  #ifdef _OPENMP
      timers.OMP = omp_get_wtime();
  #else
      timers.CPU = clock();
  #endif
}

void STOP_CPU_TIMER(SysInfo::Timers &timers)
{
  #ifdef _OPENMP
      timers.duration = omp_get_wtime() - timers.OMP;
      timers.accum += timers.duration;
  #else
      timers.CPU = clock() - timers.CPU;
      timers.duration = ((double) timers.CPU)/CLOCKS_PER_SEC;
      timers.accum += timers.duration;
  #endif
}
