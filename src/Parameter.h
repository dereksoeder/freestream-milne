#pragma once

#include <stddef.h>  // for size_t

struct parameters
{
  int OUTPUTFORMAT;
  int BARYON;
  int IC_ENERGY;
  int IC_BARYON;
  float ETA_WIDTH;
  float ETA_FLAT;
  float SIGMA;
  float SIGMA_B;
  size_t DIM_X;
  size_t DIM_Y;
  size_t DIM_ETA;
  size_t DIM_RAP;
  size_t DIM_PHIP;
  float DX;
  float DY;
  float DETA;
  float DRAP;
  float DTAU;
  float TAU0;
  int EOS_TYPE;
  float E_FREEZE;
  int VISCOUS_MATCHING;
  int E_DEP_FS;
  float E_R;
  float TAU_R;
  float ALPHA;
  //these are computed based on the chosen parameters above; they are constrained
  size_t DIM;
  float TAU;
};
