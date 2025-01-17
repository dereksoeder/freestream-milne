#include <stddef.h>
#include <math.h>
#include "FSUtilities.h"

float linearInterp3D(float x0, float x1, float x2,
                      float a000, float a100, float a010, float a001,
                      float a110, float a101, float a011, float a111)
{
  float result = 0.0;
  result = ((1-x0) * (1-x1) * (1-x2) * a000)
            + ((x0) * (1-x1) * (1-x2) * a100)
            + ((1-x0) * (x1) * (1-x2) * a010)
            + ((1-x0) * (1-x1) * (x2) * a001)
            + ((x0) * (x1) * (1-x2) * a110)
            + ((x0) * (1-x1) * (x2) * a101)
            + ((1-x0) * (x1) * (x2) * a011)
            + ((x0) * (x1) * (x2)  * a111);

  return result;
}

float linearInterp2D(float x0, float x1,
                      float a00, float a10, float a01, float a11)
{
  float result = 0.0;
  result = ((1-x0) * (1-x1) * a00)
            + ((x0) * (1-x1) * a10)
            + ((1-x0) * (x1) * a01)
            + ((x0) * (x1) * a11);

  return result;
}

size_t linearIndex(size_t ix, size_t iy, size_t iz, size_t dimx, size_t dimy)
{
  size_t is = ix + (dimx * iy) + (dimx * dimy * iz);
  return is;
}
