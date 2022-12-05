#pragma once

#include <stddef.h>  // for size_t

float ** calloc2dArrayf(float **array, size_t dim1, size_t dim2);
float *** calloc3dArrayf(float ***array, size_t dim1, size_t dim2, size_t dim3);
float **** calloc4dArrayf(float ****array, size_t dim1, size_t dim2, size_t dim3, size_t dim4);
float ***** calloc5dArrayf(float *****array, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5);
void free2dArrayf(float **array, size_t dim1);
void free3dArrayf(float ***array, size_t dim1, size_t dim2);
void free4dArrayf(float ****array, size_t dim1, size_t dim2, size_t dim3);
void free5dArrayf(float *****array, size_t dim1, size_t dim2, size_t dim3, size_t dim4);
float ** malloc2dArrayf(float **array, size_t dim1, size_t dim2);
float *** malloc3dArrayf(float ***array, size_t dim1, size_t dim2, size_t dim3);
float **** malloc4dArrayf(float ****array, size_t dim1, size_t dim2, size_t dim3, size_t dim4);
float ***** malloc5dArrayf(float *****array, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5);
