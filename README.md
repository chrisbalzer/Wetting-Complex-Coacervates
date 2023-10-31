# Wetting of Complex Coacervates
This repository generally contains routines to study interfacial behavior for polyelectrolyte complex coacervates. This code calculates the density profiles of a polyanion, polyanion and salt (+/-) near a single surface based on [Wetting Behavior of Complex Coacervates](https://doi.org/10.1039/D2SM00859A).

**Author** - Christopher Balzer

## Compiling Code
Navigate to ``run/`` and compile using the Makefile with ``make``.  Note that the default compiler is ``g++-9``, but any ``g++`` compiler that has ``OpenMP`` support will work.

## Starting a calculation
Once you have compiled, edit the ``input.dat`` file as necessary. Then run the command

```
./WettingCC
```

The output will automatically be generated in the ``run/example/`` folder.

## Dependencies
This code requires the header library [Eigen](https://gitlab.com/libeigen/eigen) to run. For convenience, the necessary parts of Eigen are included in ``src/external``.
