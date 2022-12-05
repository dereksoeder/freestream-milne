#pragma once

#include <stddef.h>  // for size_t

float linearInterp3D(float x0, float x1, float x2,
                      float a000, float a100, float a010, float a001,
                      float a110, float a101, float a011, float a111);

float linearInterp2D(float x0, float x1,
                      float a00, float a10, float a01, float a11);

size_t linearIndex(size_t ix, size_t iy, size_t iz, size_t dimx, size_t dimy);
