#pragma once

#include "Parameter.h"

void calculatePressure(float *energyDensity, float *baryonDensity, float *pressure, const parameters & params);
float temperatureFromEnergyDensity(float eps);
