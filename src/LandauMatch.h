#pragma once

#include "Parameter.h"

void calculateHypertrigTable(float ****hypertrigTable, const parameters & params);
void calculateStressTensor(float **stressTensor, float ***shiftedDensity, float ****hypertrigTable, const parameters & params);
void calculateBaryonCurrent(float **baryonCurrent, float ***shiftedChargeDensity, float ****hypertrigTable, const parameters & params);
void solveEigenSystem(float **stressTensor, float *energyDensity, float **flowVelocity, const parameters & params);
void calculateBulkPressure(float **stressTensor, float *energyDensity, float *pressure, float *bulkPressure, const parameters & params);
void calculateShearViscTensor(float **stressTensor, float *energyDensity, float **flowVelocity, float *pressure, float *bulkPressure, float **shearTensor, const parameters & params);
void calculateBaryonDensity(float *baryonDensity, float **baryonCurrent, float **flowVelocity, const parameters & params);
void calculateBaryonDiffusion(float **baryonDiffusion, float **baryonCurrent, float *baryonDensity, float **flowVelocity, const parameters & params);
void calculateThermalVorticityTensor(float *energyDensity, float **flowVelocity, float **thermalVelocityVector, float **thermalVorticityTensor, const parameters & params);
