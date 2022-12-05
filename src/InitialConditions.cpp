#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include "InitialConditions.h"

#ifndef THETA_FUNCTION
#define THETA_FUNCTION(X) ((float)X < (float)0 ? (float)0 : (float)1)
#endif

void initializeZero(float *density, const parameters & params)
{
  size_t DIM = params.DIM;
  for (size_t is = 0; is < DIM; is++)
  {
    density[is] = 0.0;
  }
}

void initializeGauss(float *density, float b, const parameters & params) // b is the variance ('spherically' symmetric)
{
  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (size_t is = 0; is < DIM; is++)
  {
    size_t ix = (is % DIM_X);
    size_t iy = ((is / DIM_X) % DIM_Y);
    size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    density[is] = e0 * exp(-(1.0 / b) * ((x * x) + (y * y) + (eta * eta)));
  }
}

void initializeEllipticalGauss(float *density, float bx, float by, float beta, const parameters & params) // bx is the x variance etc...
{
  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (size_t is = 0; is < DIM; is++)
  {
    size_t ix = (is % DIM_X);
    size_t iy = ((is / DIM_X) % DIM_Y);
    size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    density[is] = e0 * exp(-(1.0 / bx) * (x * x)) * exp(-(1.0 / by) * (y * y)) * exp(-(1.0 / beta) * (eta * eta));
  }
}

void initializeMCGauss(float * density, float b, const parameters & params)
{
  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (size_t is = 0; is < DIM; is++)
  {
    size_t ix = (is % DIM_X);
    size_t iy = ((is / DIM_X) % DIM_Y);
    size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    density[is] = e0 * ((float)rand() / RAND_MAX) * exp(-(1.0 / b) * ((x * x) + (y * y) + (eta * eta)));
  }
}

void initializeEllipticalMCGauss(float *density, float bx, float by, float beta, const parameters & params) // bx is the x variance etc...
{
  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (size_t is = 0; is < DIM; is++)
  {
    size_t ix = (is % DIM_X);
    size_t iy = ((is / DIM_X) % DIM_Y);
    size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    density[is] = e0 * ((float)rand() / RAND_MAX) * exp(-(1.0 / bx) * (x * x)) * exp(-(1.0 / by) * (y * y)) * exp(-(1.0 / beta) * (eta * eta));
  }
}

void readEnergyDensitySuperMCBlock(float *density, const parameters & params)
{
  float lower_tolerance = 1.0e-3;

  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float ETA_WIDTH = params.ETA_WIDTH;
  float ETA_FLAT = params.ETA_FLAT;
  float DETA = params.DETA;

  //first read in the transverse profile from superMC block data format
  float temp = 0.0;
  std::ifstream superMCFile;
  superMCFile.open("initial_profiles/ed.dat");
  if (superMCFile.is_open())
  {
    for (size_t ix = 0; ix < DIM_X; ix++)
    {
      for (size_t iy = 0; iy < DIM_Y; iy++)
      {
        superMCFile >> temp;
        for (size_t ieta = 0; ieta < DIM_ETA; ieta++) //copy the same value for all eta, then we will multiply by eta dependent function
        {
          size_t is = ix + (DIM_X * iy) + (DIM_X * DIM_Y * ieta); //the column packed index spanning x, y, z
          density[is] = temp;
        }
      }
    }
  }

  else
  {
    printf("Could not find initial profile in initial_superMC_ed!");
  }

  superMCFile.close();

  //now multiply by an eta-dependent profile; etaWidth is the width of the eta profile
  for (size_t is = 0; is < DIM; is++)
  {
    size_t ix = (is % DIM_X);
    size_t iy = ((is / DIM_X) % DIM_Y);
    size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;
    //here we use a the same profile as GPU-VH (see arXiv:1608.06577v1 p. 38)
    float arg = (-1.0) * (abs(eta) - ETA_FLAT) * (abs(eta) - ETA_FLAT) / (2.0 * ETA_WIDTH * ETA_WIDTH);
    arg = arg * THETA_FUNCTION(abs(eta) - ETA_FLAT);
    density[is] = density[is] * exp(arg) + lower_tolerance;
  }
}

void readEnergyDensityTRENTOBlock(float *density, const parameters & params)
{
  float lower_tolerance = 1.0e-3;

  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float ETA_WIDTH = params.ETA_WIDTH;
  float ETA_FLAT = params.ETA_FLAT;
  float DETA = params.DETA;

  //first read in the transverse profile from superMC block data format
  float temp = 0.0;
  std::ifstream superMCFile;
  superMCFile.open("initial_profiles/e.dat");
  if (superMCFile.is_open())
  {
    //skip the eight line (l) header
    std::string line;
    for (int l = 0; l < 12; l++) getline(superMCFile, line);
    for (size_t iy = 0; iy < DIM_Y; iy++)
    {
      for (size_t ix = 0; ix < DIM_X; ix++)
      {
        superMCFile >> temp;
        for (size_t ieta = 0; ieta < DIM_ETA; ieta++) //copy the same value for all eta, then we will multiply by eta dependent function
        {
          size_t is = ix + (DIM_X * iy) + (DIM_X * DIM_Y * ieta); //the column packed index spanning x, y, z
          density[is] = temp;
        }
      }
    }
  }

  else
  {
    printf("Could not find initial profile in initial_profiles!");
  }

  superMCFile.close();

  //now multiply by an eta-dependent profile; etaWidth is the width of the eta profile
  for (size_t is = 0; is < DIM; is++)
  {
    size_t ix = (is % DIM_X);
    size_t iy = ((is / DIM_X) % DIM_Y);
    size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;
    //here we use a the same profile as GPU-VH (see arXiv:1608.06577v1 p. 38)
    float arg = (-1.0) * (abs(eta) - ETA_FLAT) * (abs(eta) - ETA_FLAT) / (2.0 * ETA_WIDTH * ETA_WIDTH);
    arg = arg * THETA_FUNCTION(abs(eta) - ETA_FLAT);
    density[is] = density[is] * exp(arg) + lower_tolerance;
  }
}

void initialize2Gaussians(float *density, float bx, float by, float beta, const parameters & params) // bx is the x variance etc...
{
  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;

  float e0 = 0.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (size_t is = 0; is < DIM; is++)
  {
    size_t ix = (is % DIM_X);
    size_t iy = ((is / DIM_X) % DIM_Y);
    size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    float x1 = 3.0;
    float y1 = 0.0;
    float eta1 = 0.0;

    float x2 = -3.0;
    float y2 = 0.0;
    float eta2 = 0.0;
    density[is] = e0 * (exp(-(1.0 / bx) * ((x-x1) * (x-x1))) * exp(-(1.0 / by) * ((y-y1) * (y-y1))) * exp(-(1.0 / beta) * ((eta-eta1) * (eta-eta1)))
      + exp(-(1.0 / bx) * ((x-x2) * (x-x2))) * exp(-(1.0 / by) * ((y-y2) * (y-y2))) * exp(-(1.0 / beta) * ((eta-eta2) * (eta-eta2))) );
    density[is] *= ((float)rand() / RAND_MAX);
  }
}

void readEnergyDensityTRENTO3DBlock(float *density, const parameters & params)
{
  float lower_tolerance = 1.0e-3;

  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float ETA_WIDTH = params.ETA_WIDTH;
  float ETA_FLAT = params.ETA_FLAT;
  float DETA = params.DETA;

  //first read in the transverse profile from superMC block data format
  float temp = 0.0;
  std::ifstream superMCFile;
  superMCFile.open("initial_profiles/e.dat");
  if (superMCFile.is_open())
  {
    //skip the eight line (l) header
    //std::string line;
    //for (int l = 0; l < 8; l++) getline(superMCFile, line);
    for (size_t ix = 0; ix < DIM_X; ix++)
    {
      for (size_t iy = 0; iy < DIM_Y; iy++)
      {
        for (size_t ieta = 0; ieta < DIM_ETA; ieta++)
        {
          size_t is = ix + (DIM_X * iy) + (DIM_X * DIM_Y * ieta); //the column packed index spanning x, y, z
          superMCFile >> temp;
          density[is] = temp + lower_tolerance;
        }
      }
    }
  }

  else
  {
    printf("Could not find initial profile in initial_profiles!");
  }
  superMCFile.close();
}

void readEnergyDensityCPUVH(float *density, const parameters & params)
{
  float lower_tolerance = 1.0e-7;
  float scale_factor = 1.0;

  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float ETA_WIDTH = params.ETA_WIDTH;
  float ETA_FLAT = params.ETA_FLAT;
  float DETA = params.DETA;

  //first read in the transverse profile from superMC block data format
  float temp = 0.0;
  float dummy = 0.0;
  std::ifstream superMCFile;
  superMCFile.open("initial_profiles/e.dat");
  if (superMCFile.is_open())
  {
    for (size_t ieta = 0; ieta < DIM_ETA; ieta++)
    {
      for (size_t iy = 0; iy < DIM_Y; iy++)
      {
        for (size_t ix = 0; ix < DIM_X; ix++)
        {
          size_t is = ix + (DIM_X * iy) + (DIM_X * DIM_Y * ieta); //the column packed index spanning x, y, z
          superMCFile >> dummy >> dummy >> dummy >> temp;
          density[is] = temp * scale_factor + lower_tolerance;
        }
      }
    }
  }

  else
  {
    printf("Could not find initial profile in initial_profiles!");
  }
  superMCFile.close();
}
