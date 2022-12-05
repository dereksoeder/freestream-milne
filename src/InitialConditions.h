#pragma once

#include "Parameter.h"

void initializeZero(float *density, const parameters & params);
void initializeGauss(float *density, float b, const parameters & params); // b is the variance ('spherically' symmetric)
void initializeEllipticalGauss(float *density, float bx, float by, float beta, const parameters & params); // bx is the x variance etc...
void initializeMCGauss(float * density, float b, const parameters & params);
void initializeEllipticalMCGauss(float *density, float bx, float by, float beta, const parameters & params); // bx is the x variance etc...
void readEnergyDensitySuperMCBlock(float *density, const parameters & params);
void readEnergyDensityTRENTOBlock(float *density, const parameters & params);
void initialize2Gaussians(float *density, float bx, float by, float beta, const parameters & params); // bx is the x variance etc...
void readEnergyDensityTRENTO3DBlock(float *density, const parameters & params);
void readEnergyDensityCPUVH(float *density, const parameters & params);
