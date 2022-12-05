#pragma once

#include <stddef.h>
#include "Parameter.h"

float getX(size_t is, const parameters & params);
float getY(size_t is, const parameters & params);
float getEta(size_t is, const parameters & params);
void freeStream(float **density, float ***shiftedDensity, const parameters & params);
void convertInitialDensity(float *initialEnergyDensity, float **density, const parameters & params);
void convertInitialChargeDensity(float *initialChargeDensity, float **chargeDensity, const parameters & params);
float getEnergyDependentTau(float *initialEnergyDensity, const parameters & params);
