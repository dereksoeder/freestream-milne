#pragma once

#include <stddef.h>
#include "Parameter.h"

void writeScalarToFile(float *var, char name[255], const parameters & params);
void writeVectorToFile(float **var, char name[255], size_t idx, const parameters & params);
void writeScalarToFileProjection(float *var, char name[255], const parameters & params);
void writeVectorToFileProjection(float **var, char name[255], size_t idx, const parameters & params);
void readDensityFile(float *density, char name[255], const parameters & params);
void readInParameters(struct parameters &params, const char * inputPath = nullptr);
//void readEoSTable();
