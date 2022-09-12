/* ===================================================================
Copyright (c) 2022 Chris Balzer  (balzer@caltech.edu)
======================================================================*/

#ifndef UTILITES_H
#define UTILITES_H

#include "../global_external.h"


/* Namespace for structures*/
namespace SysInfo
{
    /* Structure for timing */
    struct Timers
    {
      clock_t CPU;
      double  OMP;
      double  duration;
      double  accum;

      Timers() {
          accum = 0.0;
      }
    };
}

/* Utility functions */
void GET_OMP_THREADS();
void START_CPU_TIMER(SysInfo::Timers &timers);
void STOP_CPU_TIMER(SysInfo::Timers &timers);

#endif //UTILITES_H
