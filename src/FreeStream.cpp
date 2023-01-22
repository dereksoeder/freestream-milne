#include <stddef.h>
#include <math.h>
#include "FreeStream.h"
#include "FSConfig.h"
#include "FSUtilities.h"
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>

#ifndef _OPENMP  // from MUSIC/src/init.cxx
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#else
    #include <omp.h>
#endif

#ifdef _OPENACC
#include <accelmath.h>
#endif

float getX(size_t is, const parameters & params)
{
  float xmin = (-1.0) * ((float)(params.DIM_X-1) / 2.0) * params.DX;
  return xmin + static_cast<float>(is % params.DIM_X) * params.DX;
}

float getY(size_t is, const parameters & params)
{
  float ymin = (-1.0) * ((float)(params.DIM_Y-1) / 2.0) * params.DY;
  return ymin + static_cast<float>((is / params.DIM_X) % params.DIM_Y) * params.DY;
}

float getEta(size_t is, const parameters & params)
{
  float etamin = (-1.0) * ((float)(params.DIM_ETA-1) / 2.0) * params.DETA;
  return etamin + static_cast<float>((is / params.DIM_X / params.DIM_Y) % params.DIM_ETA) * params.DETA;
}

//#pragma acc routine //build a copy of function to run on device
void freeStream(float **density, float ***shiftedDensity, const parameters & params)
{
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  size_t DIM_RAP = params.DIM_RAP;
  size_t DIM_PHIP = params.DIM_PHIP;
  size_t DIM = params.DIM;
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
  for (size_t ix = 0; ix < DIM_X; ix++)
  {
    for (size_t iy = 0; iy < DIM_Y; iy++)
    {
      //note that internal indexing in gsl is tranposed!
      size_t is = ix + (DIM_X * iy);
      x_vals[ix] = getX(is, params);
      y_vals[iy] = getY(is, params);
      size_t is_gsl = is;  // maybe it's no longer transposed relative to our coordinates, since now we're y-major too?
      density_vals[is_gsl] = density[is][0];
    }
  }

  double cosphip[DIM_PHIP];  // precompute these to avoid DIM_ETA * DIM_Y * DIM_X * DIM_RAP redundant cos and sin calls
  double sinphip[DIM_PHIP];
  for (size_t iphip = 0; iphip < DIM_PHIP; iphip++)
  {
    float phip = float(iphip) * (2.0 * M_PI) / float(DIM_PHIP);
    cosphip[iphip] = cos(phip);
    sinphip[iphip] = sin(phip);
  }

  double x_min_interp = x_vals[0];
  double y_min_interp = y_vals[0];

  double x_max_interp = x_vals[DIM_X-1];
  double y_max_interp = y_vals[DIM_Y-1];

  size_t nx = sizeof(x_vals) / sizeof(x_vals[0]);
  size_t ny = sizeof(y_vals) / sizeof(y_vals[0]);
  gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);

  int naccs = omp_get_max_threads();  // allocate one interpolation accelerator per thread now, since it might not be safe to allocate/free them concurrently
  gsl_interp_accel * xaccs[naccs];    //   (accelerators are currently only used in the boost-invariant case, so since we're threading over eta, there won't be any concurrent use; still, let's leave it in case things change some day)
  gsl_interp_accel * yaccs[naccs];
  for (int iacc = 0; iacc < naccs; iacc++)
  {
    xaccs[iacc] = gsl_interp_accel_alloc();
    yaccs[iacc] = gsl_interp_accel_alloc();
  }

  gsl_spline2d_init(spline, x_vals, y_vals, density_vals, nx, ny);
  std::cout << "Bicubic splines initialized..."<< "\n";

 // separated the original `is` single-index loop into three loops, so that the (outermost) `ieta` loop can be parallelized,
 //  so threads do more work for the overhead (not optimal though if DIM_ETA << number of threads)

 #pragma omp parallel for
 for (size_t ieta = 0; ieta < DIM_ETA; ieta++)
 {
  size_t is = ieta * DIM_X * DIM_Y;

  gsl_interp_accel * xacc = nullptr;
  gsl_interp_accel * yacc = nullptr;

  int iacc = omp_get_thread_num();  // this should be unique among all threads currently running--we don't want two threads modifying the same accelerator at the same time
  if ((iacc >= 0) && (iacc < naccs))
  {
    xacc = xaccs[iacc];
    yacc = yaccs[iacc];
  }

  float eta = (float)ieta * DETA  + etamin;

  for (size_t iy = 0; iy < DIM_Y; iy++)
  {
   float y = (float)iy * DY  + ymin;

   for (size_t ix = 0; ix < DIM_X; is++, ix++)
   {
    float x = (float)ix * DX  + xmin;

    for (size_t irap = 0; irap < DIM_RAP; irap++)
    {
      //float rap = (float)irap * DRAP + rapmin;

      //try evaluating at values of rapidity y centered around y ~= eta
      //if (DIM_ETA > 1) rap = rap + eta;

      //w is an integration variable on the domain (-1,1) - careful not to include endpoints (nans)
      //float w =  -.9975 + (float)irap * (1.995 / (float)(DIM_RAP - 1));
      float w =  (DIM_RAP <= 1) ? 0. : -.975 + (float)irap * (1.95 / (float)(DIM_RAP - 1));  // still allows extreme rapidities (~25), but at least sinh and cosh won't return nan (float is limited to ~3.4E+38)
      float rap = eta + tan((M_PI / 2.0) * w );

      float eta_new;
      if (DIM_ETA == 1)
      {
        eta_new = 0.0;
      }
      else if (DIM_ETA > 1)
      {
        eta_new = asinh( (TAU / TAU0) * sinh(eta - rap) ) + rap; //old formula works
      }

      float coshrapeta    = cosh(rap - eta);
      float coshrapetanew = cosh(rap - eta_new);

      for (size_t iphip = 0; iphip < DIM_PHIP; iphip++)
      {
        float x_new, y_new; //the shifted coordinates
        if (DIM_ETA == 1)
        {
          x_new = x - cosphip[iphip] * DTAU;
          y_new = y - sinphip[iphip] * DTAU;
        }
        else if (DIM_ETA > 1)
        {
          x_new = x - cosphip[iphip] * (TAU * coshrapetanew - TAU0 * coshrapeta);
          y_new = y - sinphip[iphip] * (TAU * coshrapetanew - TAU0 * coshrapeta);
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
            size_t is_new_11 = (size_t)ix_new_f + (DIM_X * (size_t)iy_new_f) + (DIM_X * DIM_Y * (size_t)ieta_new_f);
            size_t is_new_21 = ( (size_t)ix_new_f + 1) + (DIM_X * (size_t)iy_new_f) + (DIM_X * DIM_Y * (size_t)ieta_new_f);
            size_t is_new_12 = (size_t)ix_new_f + (DIM_X * ( (size_t)iy_new_f + 1) ) + (DIM_X * DIM_Y * (size_t)ieta_new_f);
            size_t is_new_22 = ( (size_t)ix_new_f + 1) + (DIM_X * ( (size_t)iy_new_f + 1) ) + (DIM_X * DIM_Y * (size_t)ieta_new_f);

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
            size_t is_new_000 = (size_t)ix_new_f + (DIM_X * (size_t)iy_new_f) + (DIM_X * DIM_Y * (size_t)ieta_new_f);
            size_t is_new_100 = ( (size_t)ix_new_f + 1) + (DIM_X * (size_t)iy_new_f) + (DIM_X * DIM_Y * (size_t)ieta_new_f);
            size_t is_new_010 = (size_t)ix_new_f + (DIM_X * ( (size_t)iy_new_f + 1) ) + (DIM_X * DIM_Y * (size_t)ieta_new_f);
            size_t is_new_110 = ( (size_t)ix_new_f + 1) + (DIM_X * ( (size_t)iy_new_f + 1) ) + (DIM_X * DIM_Y * (size_t)ieta_new_f);
            size_t is_new_001 = (size_t)ix_new_f + (DIM_X * (size_t)iy_new_f) + (DIM_X * DIM_Y * ( (size_t)ieta_new_f + 1) );
            size_t is_new_101 = ( (size_t)ix_new_f + 1) + (DIM_X * (size_t)iy_new_f) + (DIM_X * DIM_Y * ( (size_t)ieta_new_f + 1) );
            size_t is_new_011 = (size_t)ix_new_f + (DIM_X * ( (size_t)iy_new_f + 1) ) + (DIM_X * DIM_Y * ( (size_t)ieta_new_f + 1) );
            size_t is_new_111 = ( (size_t)ix_new_f + 1) + (DIM_X * ( (size_t)iy_new_f + 1) ) + (DIM_X * DIM_Y * ( (size_t)ieta_new_f + 1) );

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
      } //for (size_t iphip = 0; iphip < DIM_PHIP; iphip++)
    } //for (size_t irap = 0; irap < DIM_RAP; irap++)
   } //for (size_t ix = 0; ix < DIM_X; is++, ix++)
  } //for (size_t iy = 0; iy < DIM_Y; iy++)
 } //for (size_t ieta = 0; ieta < DIM_ETA; ieta++)

  gsl_spline2d_free(spline);
  for (int iacc = 0; iacc < naccs; iacc++)
  {
    gsl_interp_accel_free(xaccs[iacc]);
    gsl_interp_accel_free(yaccs[iacc]);
  }
}

//this creates the initial G^(tau,tau) function, a function of spatial coordinates and rapidity
//in the special case if 2+1D, this calculates F^(tau,tau) in the notation of (PRC 91, 064906)
//rapidity dependence is determined by the assumption for rapidity dependence of the initial distribution function
void convertInitialDensity(float *initialEnergyDensity, float **density, const parameters & params)
{
  float SIGMA = params.SIGMA;
  size_t DIM = params.DIM;
  float TAU0 = params.TAU0;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;

  size_t DIM_RAP = params.DIM_RAP;
  //float DRAP = params.DRAP;
  float DETA = params.DETA;

  float n = (SIGMA <= 0.) ? 1. : sqrt(M_PI / 2.0) * SIGMA * (1.0 + exp(2.0 * SIGMA * SIGMA)); //the integral over cosh^2 * exp()
  float norm_factor = 1.0 / (2.0 * M_PI * n); //the normalization constant relating the intial energy density to the intial density profile G(tilde)^(tau,tau)

  //float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;

  if (DIM_ETA == 1) //catch the special case of 2+1D freestreaming; note DIM_RAP must also be 1 ! the normalization with SIGMA -> 0 does not generalize?
  {
    for (size_t is = 0; is < DIM; is++)
    {
      size_t ix = (is % DIM_X);
      size_t iy = ((is / DIM_X) % DIM_Y);

      for (size_t irap = 0; irap < DIM_RAP; irap++)
      {
        //density[is][irap] = initialEnergyDensity[is] * (TAU0 / (2.0 * M_PI)); //this is initial F^(tau,tau) in the notation of (PRC 91, 064906)
        density[is][irap] = initialEnergyDensity[is] / (2.0 * M_PI); //this is initial F^(tau,tau) in the notation of (PRC 91, 064906)

      }
    }
  }

  else
  {
    for (size_t is = 0; is < DIM; is++)
    {
      size_t ix = (is % DIM_X);
      size_t iy = ((is / DIM_X) % DIM_Y);
      size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

      float eta = (float)ieta * DETA  + etamin;

      for (size_t irap = 0; irap < DIM_RAP; irap++)
      {
        //float rap = (float)irap * DRAP + rapmin;

        //try evaluating at values of rapidity y centered around y ~= eta
        //rap = rap + eta;

        //w is an integration variable on the domain (-1,1) - careful not to include endpoints (nans)
        //float w =  -.9975 + (float)irap * (1.995 / (float)(DIM_RAP - 1));
        float w =  (DIM_RAP <= 1) ? 0. : -.975 + (float)irap * (1.95 / (float)(DIM_RAP - 1));  // still allows extreme rapidities (~25), but at least cosh won't return nan (float is limited to ~3.4E+38)
        float rap = eta + tan((M_PI / 2.0) * w );

        float rap_factor = (SIGMA <= 0.) ? 1. : cosh(eta - rap) * cosh(eta - rap) * exp( (-1.0) * (eta - rap) * (eta - rap) / (2.0 * SIGMA * SIGMA) );
        density[is][irap] = initialEnergyDensity[is] * norm_factor * rap_factor; //this is initial G^(tau,tau)
      }
    }
  }
}
//this creates the initial J^(tau) function, a function of spatial coordinates and rapidity
//rapidity dependence is determined by the assumption for rapidity dependence of the initial baryon distribution function
void convertInitialChargeDensity(float *initialChargeDensity, float **chargeDensity, const parameters & params)
{
  float SIGMA_B = params.SIGMA_B;
  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;

  size_t DIM_RAP = params.DIM_RAP;
  float DRAP = params.DRAP;
  float DETA = params.DETA;

  float n = (SIGMA_B <= 0.) ? 1. : sqrt(2.0 * M_PI) * SIGMA_B * exp(SIGMA_B * SIGMA_B / 2.0); //the integral over cosh * exp()
  float norm_factor = 1.0 / (2.0 * M_PI * n); //the normalization constant relating the intial baryon density to the intial charge density profile J(tilde)^(tau)

  float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;

  for (size_t is = 0; is < DIM; is++)
  {
    size_t ix = (is % DIM_X);
    size_t iy = ((is / DIM_X) % DIM_Y);
    size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

    float eta = (float)ieta * DETA  + etamin;

    for (size_t irap = 0; irap < DIM_RAP; irap++)
    {
      float rap = (float)irap * DRAP + rapmin;
      float rap_factor = (SIGMA_B <= 0.) ? 1. : cosh(eta - rap) * exp((-1.0) * (eta - rap) * (eta - rap) / (SIGMA_B * SIGMA_B));
      chargeDensity[is][irap] = initialChargeDensity[is] * norm_factor * rap_factor; //this is initial J^(tau)
    }

  }
}

float getEnergyDependentTau(float *initialEnergyDensity, const parameters & params)
{
  float tau_R = params.TAU_R;
  float e_R = params.E_R;
  float alpha = params.ALPHA;
  float hbarc = 0.197326938;

  size_t DIM = params.DIM;
  float dx = params.DX;
  float dy = params.DY;

  //integrate over transverse plane to find transverse averaged energy density
  //defined by (int dx dy (eps^2) ) / (int dx dy eps)
  float numerator = 0.0;
  float denominator = 0.0;

  for (size_t is = 0; is < DIM; is++)
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
