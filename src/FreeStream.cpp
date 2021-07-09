#pragma once
#include <math.h>
#include "Parameter.h"
#include "FSConfig.h"
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>

#ifdef _OPENACC
#include <accelmath.h>
#endif

float getX(int is, parameters params)
{
  float xmin = (-1.0) * ((float)(params.DIM_X-1) / 2.0) * params.DX;
  return xmin + static_cast<float>(is % params.DIM_X) * params.DX;
}

float getY(int is, parameters params)
{
  float ymin = (-1.0) * ((float)(params.DIM_Y-1) / 2.0) * params.DY;
  return ymin + static_cast<float>((is / params.DIM_X) % params.DIM_Y) * params.DY;
}

float getEta(int is, parameters params)
{
  float etamin = (-1.0) * ((float)(params.DIM_ETA-1) / 2.0) * params.DETA;
  return etamin + static_cast<float>((is / params.DIM_X / params.DIM_Y) % params.DIM_ETA) * params.DETA;
}

float linearInterp3D(float x0, float x1, float x2,
                      float a000, float a100, float a010, float a001,
                      float a110, float a101, float a011, float a111)
{
  float result = 0.0;
  result = ( (1-x0) * (1-x1) * (1-x2) * a000 )
            + ( (x0) * (1-x1) * (1-x2) * a100 )
            + ( (1-x0) * (x1) * (1-x2) * a010 )
            + ( (1-x0) * (1-x1) * (x2) * a001 )
            + ( (x0) * (x1) * (1-x2) * a110 )
            + ( (x0) * (1-x1) * (x2) * a101 )
            + ( (1-x0) * (x1) * (x2) * a011 )
            + ( (x0) * (x1) * (x2)  * a111 );

  return result;
}

float linearInterp2D(float x0, float x1,
                      float a00, float a10, float a01, float a11)
{
  float result = 0.0;
  result = ( (1-x0) * (1-x1) * a00 )
            + ( (x0) * (1-x1) * a10 )
            + ( (1-x0) * (x1) * a01 )
            + ( (x0) * (x1) * a11 );

  return result;
}

//#pragma acc routine //build a copy of function to run on device
void freeStream(float **density, float ***shiftedDensity, parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  int DIM_RAP = params.DIM_RAP;
  int DIM_PHIP = params.DIM_PHIP;
  int DIM = params.DIM;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;
  //float DRAP = params.DRAP;
  float TAU0 = params.TAU0;
  float TAU = params.TAU;
  float DTAU = params.DTAU;

  float xmin = (-1.0) * ((float)(DIM_X-1) / 2.0) * DX;
  float ymin = (-1.0) * ((float)(DIM_Y-1) / 2.0) * DY;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;
  //float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;

  //this is to avoid GSL failure in very rare cases when freestreaming time is long for fluctuating events
  //gsl_error_handler_t * gsl_set_error_handler_off();

  //for case of 2+1D, set up the bicubic interpolating function of initial density
  std::cout << "Setting up bicubic splines ..."<< "\n";
  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  double x_vals[DIM_X];
  double y_vals[DIM_Y];
  double density_vals[DIM_X * DIM_Y];

  //fill arrays with values of coordinates and density values
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      //note that internal indexing in gsl is tranposed!
      int is = ix + (DIM_X * iy);
      x_vals[ix] = getX(is, params);
      y_vals[iy] = getY(is, params);
      int is_gsl = is;  // maybe it's no longer transposed relative to our coordinates, since now we're y-major too?
      density_vals[is_gsl] = density[is][0];
    }
  }

  double x_min_interp = x_vals[0];
  double y_min_interp = y_vals[0];

  double x_max_interp = x_vals[DIM_X-1];
  double y_max_interp = y_vals[DIM_Y-1];

  size_t nx = sizeof(x_vals) / sizeof(x_vals[0]);
  size_t ny = sizeof(y_vals) / sizeof(y_vals[0]);
  gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  gsl_spline2d_init(spline, x_vals, y_vals, density_vals, nx, ny);
  std::cout << "Bicubic splines initialized..."<< "\n";

  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    int ix = (is % DIM_X);
    int iy = ((is / DIM_X) % DIM_Y);
    int ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

    float x = (float)ix * DX  + xmin;
    float y = (float)iy * DY  + ymin;
    float eta = (float)ieta * DETA  + etamin;
    for (int irap = 0; irap < DIM_RAP; irap++)
    {
      //float rap = (float)irap * DRAP + rapmin;

      //try evaluating at values of rapidity y centered around y ~= eta
      //if (DIM_ETA > 1) rap = rap + eta;

      //w is an integration variable on the domain (-1,1) - careful not to include endpoints (nans)
      float w =  -.9975 + (float)irap * (1.995 / (float)(DIM_RAP - 1));
      float rap = eta + tan((M_PI / 2.0) * w );

      for (int iphip = 0; iphip < DIM_PHIP; iphip++)
      {
        float phip = float(iphip) * (2.0 * M_PI) / float(DIM_PHIP);

        //float eta_new = asinh( (TAU0 / TAU) * sinh(eta - rap) ) + rap;

        float eta_new, x_new, y_new; //the shifted coordinates
        if (DIM_ETA == 1)
        {
          eta_new = 0.0;
          x_new   = x - cos(phip) * DTAU;
          y_new   = y - sin(phip) * DTAU;
        }
        else if (DIM_ETA > 1)
        {
          eta_new = asinh( (TAU / TAU0) * sinh(eta - rap) ) + rap; //old formula works
          x_new   = x - cos(phip) * (TAU * cosh(rap - eta_new) - TAU0 * cosh(rap - eta));
          y_new   = y - sin(phip) * (TAU * cosh(rap - eta_new) - TAU0 * cosh(rap - eta));
        }

        //get fractions for linear interpolation routine
        float ix_new = (x_new - xmin) / DX;
        float iy_new = (y_new - ymin) / DY;
        float ieta_new;
        if (DIM_ETA > 1) ieta_new = (eta_new - etamin) / DETA;
        else if (DIM_ETA == 1) ieta_new = 0;

        float ix_new_f = floor(ix_new);
        float iy_new_f = floor(iy_new);
        float ieta_new_f = floor(ieta_new);

        float x_frac = ix_new - ix_new_f;
        float y_frac = iy_new - iy_new_f;
        float eta_frac = ieta_new - ieta_new_f;

        //prevent from going out of array bounds!
        if ( (ix_new_f >= 1) && (ix_new_f < DIM_X - 1) && (iy_new_f >= 1) && (iy_new_f < DIM_Y - 1) )
        {
          float interp = 0.0; //final result
          float interp_l = 0.0; // result of linear interpolation
          float interp_c = 0.0; //result of cubic spline interpolation

          //2+1D routine
          if (DIM_ETA == 1)
          {
            int is_new_11 = (int)ix_new_f + (DIM_X * (int)iy_new_f) + (DIM_X * DIM_Y * (int)ieta_new_f);
            int is_new_21 = ( (int)ix_new_f + 1) + (DIM_X * (int)iy_new_f) + (DIM_X * DIM_Y * (int)ieta_new_f);
            int is_new_12 = (int)ix_new_f + (DIM_X * ( (int)iy_new_f + 1) ) + (DIM_X * DIM_Y * (int)ieta_new_f);
            int is_new_22 = ( (int)ix_new_f + 1) + (DIM_X * ( (int)iy_new_f + 1) ) + (DIM_X * DIM_Y * (int)ieta_new_f);

            float a11 = density[is_new_11][irap];
            float a21 = density[is_new_21][irap];
            float a12 = density[is_new_12][irap];
            float a22 = density[is_new_22][irap];

            // biinear interpolation
            interp_l = linearInterp2D(x_frac, y_frac, a11, a21, a12, a22);
            // bicubic spline interpolation

            //first make sure x_new and y_new are inside bounds of interpolation
            if ( (x_min_interp < x_new) && (x_new < x_max_interp) && (y_min_interp < y_new) && (y_new < y_max_interp) )
            {
              //call bicubic interpolation
              interp_c = gsl_spline2d_eval(spline, x_new, y_new, xacc, yacc);
            }
            else interp_c = 0.0; //return zero will force us to use linear interpolation

            //if bicubic spline returns negative, use linear interpolation.
            if (interp_c > 0.0) interp = interp_c;
            else interp = interp_l;

          }

          else if ( (DIM_ETA > 1) && (ieta_new_f >= 1) && (ieta_new_f < DIM_ETA - 1) )
          {
            int is_new_000 = (int)ix_new_f + (DIM_X * (int)iy_new_f) + (DIM_X * DIM_Y * (int)ieta_new_f);
            int is_new_100 = ( (int)ix_new_f + 1) + (DIM_X * (int)iy_new_f) + (DIM_X * DIM_Y * (int)ieta_new_f);
            int is_new_010 = (int)ix_new_f + (DIM_X * ( (int)iy_new_f + 1) ) + (DIM_X * DIM_Y * (int)ieta_new_f);
            int is_new_110 = ( (int)ix_new_f + 1) + (DIM_X * ( (int)iy_new_f + 1) ) + (DIM_X * DIM_Y * (int)ieta_new_f);
            int is_new_001 = (int)ix_new_f + (DIM_X * (int)iy_new_f) + (DIM_X * DIM_Y * ( (int)ieta_new_f + 1) );
            int is_new_101 = ( (int)ix_new_f + 1) + (DIM_X * (int)iy_new_f) + (DIM_X * DIM_Y * ( (int)ieta_new_f + 1) );
            int is_new_011 = (int)ix_new_f + (DIM_X * ( (int)iy_new_f + 1) ) + (DIM_X * DIM_Y * ( (int)ieta_new_f + 1) );
            int is_new_111 = ( (int)ix_new_f + 1) + (DIM_X * ( (int)iy_new_f + 1) ) + (DIM_X * DIM_Y * ( (int)ieta_new_f + 1) );

            float a000 = density[is_new_000][irap];
            float a100 = density[is_new_100][irap];
            float a010 = density[is_new_010][irap];
            float a110 = density[is_new_110][irap];
            float a001 = density[is_new_001][irap];
            float a101 = density[is_new_101][irap];
            float a011 = density[is_new_011][irap];
            float a111 = density[is_new_111][irap];

            interp = linearInterp3D(x_frac, y_frac, eta_frac,
                                  a000, a100, a010, a001,
                                  a110, a101, a011, a111);
          }

          shiftedDensity[is][irap][iphip] = interp;
        }
      } //for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    } //for (int irap = 0; irap < DIM_RAP; irap++)
  } //for (int is = 0; is < DIM; is++)

  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
}

//this creates the initial G^(tau,tau) function, a function of spatial coordinates and rapidity
//in the special case if 2+1D, this calculates F^(tau,tau) in the notation of (PRC 91, 064906)
//rapidity dependence is determined by the assumption for rapidity dependence of the initial distribution function
void convertInitialDensity(float *initialEnergyDensity, float **density, parameters params)
{
  float SIGMA = params.SIGMA;
  int DIM = params.DIM;
  float TAU0 = params.TAU0;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;

  int DIM_RAP = params.DIM_RAP;
  //float DRAP = params.DRAP;
  float DETA = params.DETA;

  float n = sqrt(M_PI / 2.0) * SIGMA * (1.0 + exp(2.0 * SIGMA * SIGMA)); //the integral over cosh^2 * exp()
  float norm_factor = 1.0 / (2.0 * M_PI * n); //the normalization constant relating the intial energy density to the intial density profile G(tilde)^(tau,tau)

  //float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;

  if (DIM_ETA == 1) //catch the special case of 2+1D freestreaming; note DIM_RAP must also be 1 ! the normalization with SIGMA -> 0 does not generalize?
  {
    for (int is = 0; is < DIM; is++)
    {
      int ix = (is % DIM_X);
      int iy = ((is / DIM_X) % DIM_Y);

      for (int irap = 0; irap < DIM_RAP; irap++)
      {
        //density[is][irap] = initialEnergyDensity[is] * (TAU0 / (2.0 * M_PI)); //this is initial F^(tau,tau) in the notation of (PRC 91, 064906)
        density[is][irap] = initialEnergyDensity[is] / (2.0 * M_PI); //this is initial F^(tau,tau) in the notation of (PRC 91, 064906)

      }
    }
  }

  else
  {
    for (int is = 0; is < DIM; is++)
    {
      int ix = (is % DIM_X);
      int iy = ((is / DIM_X) % DIM_Y);
      int ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

      float eta = (float)ieta * DETA  + etamin;

      for (int irap = 0; irap < DIM_RAP; irap++)
      {
        //float rap = (float)irap * DRAP + rapmin;

        //try evaluating at values of rapidity y centered around y ~= eta
        //rap = rap + eta;

        //w is an integration variable on the domain (-1,1) - careful not to include endpoints (nans)
        float w =  -.9975 + (float)irap * (1.995 / (float)(DIM_RAP - 1));
        float rap = eta + tan((M_PI / 2.0) * w );

        float rap_factor = cosh(eta - rap) * cosh(eta - rap) * exp( (-1.0) * (eta - rap) * (eta - rap) / (2.0 * SIGMA * SIGMA) );
        density[is][irap] = initialEnergyDensity[is] * norm_factor * rap_factor; //this is initial G^(tau,tau)
      }
    }
  }
}
//this creates the initial J^(tau) function, a function of spatial coordinates and rapidity
//rapidity dependence is determined by the assumption for rapidity dependence of the initial baryon distribution function
void convertInitialChargeDensity(float *initialChargeDensity, float **chargeDensity, parameters params)
{
  float SIGMA_B = params.SIGMA_B;
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;

  int DIM_RAP = params.DIM_RAP;
  float DRAP = params.DRAP;
  float DETA = params.DETA;

  float n = sqrt(2.0 * M_PI) * SIGMA_B * exp(SIGMA_B * SIGMA_B / 2.0); //the integral over cosh * exp()
  float norm_factor = 1.0 / (2.0 * M_PI * n); //the normalization constant relating the intial baryon density to the intial charge density profile J(tilde)^(tau)

  float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;

  for (int is = 0; is < DIM; is++)
  {
    int ix = (is % DIM_X);
    int iy = ((is / DIM_X) % DIM_Y);
    int ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

    float eta = (float)ieta * DETA  + etamin;

    for (int irap = 0; irap < DIM_RAP; irap++)
    {
      float rap = (float)irap * DRAP + rapmin;
      float rap_factor = cosh(eta - rap) * exp((-1.0) * (eta - rap) * (eta - rap) / (SIGMA_B * SIGMA_B));
      chargeDensity[is][irap] = initialChargeDensity[is] * norm_factor * rap_factor; //this is initial J^(tau)
    }

  }
}

float getEnergyDependentTau(float *initialEnergyDensity, parameters params)
{
  float tau_R = params.TAU_R;
  float e_R = params.E_R;
  float alpha = params.ALPHA;
  float hbarc = 0.197326938;

  int DIM = params.DIM;
  float dx = params.DX;
  float dy = params.DY;

  //integrate over transverse plane to find transverse averaged energy density
  //defined by (int dx dy (eps^2) ) / (int dx dy eps)
  float numerator = 0.0;
  float denominator = 0.0;

  for (int is = 0; is < DIM; is++)
  {
    float eps = initialEnergyDensity[is]; // Multiply by hbarc for same units as e_R
    numerator += (eps*eps);
    denominator += eps;
  }

  float e_T = numerator / (denominator + 1.0e-7);
  e_T = e_T * hbarc; // Multiply by hbarc for same units as e_R

  float tau_fs = tau_R * pow( (e_T / e_R), alpha );
  return tau_fs;
}
