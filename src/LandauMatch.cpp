//trigTable is a table with 10 rows for ten combinations or p^(mu)p_(nu) normalized by the energy
#include <stddef.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
//#include <math.h>
#include "LandauMatch.h"
#include "EquationOfState.h"

void calculateHypertrigTable(float ****hypertrigTable, const parameters & params)
{
  size_t DIM_RAP = params.DIM_RAP;
  size_t DIM_PHIP = params.DIM_PHIP;
  size_t DIM_ETA = params.DIM_ETA;
  float DRAP = params.DRAP;
  float DETA = params.DETA;
  float TAU = params.TAU;

  float cosphip[DIM_PHIP];  // precompute these to avoid redundant cos and sin calls
  float sinphip[DIM_PHIP];
  for (size_t iphip = 0; iphip < DIM_PHIP; iphip++)
  {
    float phip = float(iphip) * (2.0 * M_PI) / float(DIM_PHIP);
    cosphip[iphip] = cos(phip);
    sinphip[iphip] = sin(phip);
  }

  float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;
  for (size_t irap = 0; irap < DIM_RAP; irap++)
  {
    //float rap = (float)irap * DRAP + rapmin;

    for (size_t iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      float cphip = cosphip[iphip];
      float sphip = sinphip[iphip];

      for (size_t ieta = 0; ieta < DIM_ETA; ieta++)
      {
        float eta = (float)ieta * DETA  + etamin;
        if (DIM_ETA == 1) eta = 0.0;

        //w is an integration variable on the domain (-1,1) - careful not to include endpoints (nans)
        //float w =  -.9975 + (float)irap * (1.995 / (float)(DIM_RAP - 1));
        float w =  (DIM_RAP <= 1) ? 0. : -.975 + (float)irap * (1.95 / (float)(DIM_RAP - 1));  // still allows extreme rapidities (~25), but at least cosh won't return nan (float is limited to ~3.4E+38)
        float rap = eta + tan((M_PI / 2.0) * w );
        if (DIM_ETA == 1) rap = 0.0;
        //try evaluating at values of rapidity y centered around y ~= eta
        //if (DIM_ETA > 1) rap = rap + eta;

        float coshrapeta = cosh(rap - eta);
        float tanhrapeta = tanh(rap - eta);

        hypertrigTable[0][irap][iphip][ieta] = 1.0; //p^tau, p^tau component
        hypertrigTable[1][irap][iphip][ieta] = cphip / coshrapeta; //p^tau, p^x
        hypertrigTable[2][irap][iphip][ieta] = sphip / coshrapeta; //p^tau, p^y
        hypertrigTable[3][irap][iphip][ieta] = (1.0 / TAU) * tanhrapeta; //p^tau, p^eta
        hypertrigTable[4][irap][iphip][ieta] = (cphip * cphip) / (coshrapeta * coshrapeta); //p^x, p^x
        hypertrigTable[5][irap][iphip][ieta] = (cphip * sphip) / (coshrapeta * coshrapeta); //p^x, p^y
        hypertrigTable[6][irap][iphip][ieta] = (1.0 / TAU) * (cphip * tanhrapeta) / coshrapeta; //p^x, p^eta
        hypertrigTable[7][irap][iphip][ieta] = (sphip * sphip) / (coshrapeta * coshrapeta); //p^y, p^y
        hypertrigTable[8][irap][iphip][ieta] = (1.0 / TAU) * (sphip * tanhrapeta) / coshrapeta; //p^y, p^eta
        hypertrigTable[9][irap][iphip][ieta] = (1.0 / (TAU * TAU)) * tanhrapeta * tanhrapeta; //p^eta, p^eta
      }
    }
  }
}

void calculateStressTensor(float **stressTensor, float ***shiftedDensity, float ****hypertrigTable, const parameters & params)
{
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  size_t DIM_RAP = params.DIM_RAP;
  size_t DIM_PHIP = params.DIM_PHIP;
  size_t DIM = params.DIM;
  //float DRAP = params.DRAP;
  float TAU = params.TAU;
  //float TAU = params.TAU;
  float weight_rap;

  float jacobian[DIM_RAP];  // precompute these to avoid 10 * DIM_ETA * DIM_Y * DIM_X redundant computations
  for (size_t irap = 0; irap < DIM_RAP; irap++)
  {
    //w is an integration variable on the domain (-1,1) - careful not to include endpoints (nans)
    //float w =  -.9975 + (float)irap * (1.995 / (float)(DIM_RAP - 1));
    //jacobian[irap] = (M_PI/2.0) / cos( (M_PI/2.0)*w ) / cos( (M_PI/2.0)*w ) * (1.995 / float(DIM_RAP - 1));
    float w =  (DIM_RAP <= 1) ? 0. : -.975 + (float)irap * (1.95 / (float)(DIM_RAP - 1));  // changed for consistency with hypertrigTable
    jacobian[irap] = (DIM_RAP <= 1) ? 1. : (M_PI/2.0) / cos( (M_PI/2.0)*w ) / cos( (M_PI/2.0)*w ) * (1.95 / float(DIM_RAP - 1));
  }

  float d_phip = (2.0 * M_PI) / float(DIM_PHIP);

  #pragma omp parallel for
  for (size_t ivar = 0; ivar < 10; ivar++)
  {
   float * stcompon = stressTensor[ivar];

   for (size_t is = 0, ieta = 0; ieta < DIM_ETA; ieta++)
   {
    for (size_t iy = 0; iy < DIM_Y; iy++)
    {
     for (size_t ix = 0; ix < DIM_X; is++, ix++)
     {
      float sum = 0.;
      for (size_t irap = 0; irap < DIM_RAP; irap++)
      {
        //try trapezoid rule for rapidity integral
        //if (irap == 0 || irap == DIM_RAP - 1) weight_rap = DRAP / 2.0;
        //else weight_rap = DRAP;

        for (size_t iphip = 0; iphip < DIM_PHIP; iphip++)
        {
          //check convergence!
          // T^(mu,nu) = int deta int dphip G^(mu,nu)
          //#pragma omp simd (+:stressTensor_tmp)
          if (DIM_ETA == 1) sum += shiftedDensity[is][irap][iphip] * hypertrigTable[ivar][irap][iphip][ieta];
          else sum += shiftedDensity[is][irap][iphip] * hypertrigTable[ivar][irap][iphip][ieta] * jacobian[irap];
        }
      }
      if (DIM_ETA == 1) stcompon[is] = sum * d_phip / TAU; //catch the special case of 2+1D FS (solution in PRC 91, 064906)
      else stcompon[is] = sum * d_phip; //multiply by common differential factor once
     } //for(ix)
    } //for(iy)
   } //for(ieta)
  } //for(ivar)
}

void calculateBaryonCurrent(float **baryonCurrent, float ***shiftedChargeDensity, float ****hypertrigTable, const parameters & params)
{
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  size_t DIM_RAP = params.DIM_RAP;
  size_t DIM_PHIP = params.DIM_PHIP;
  size_t DIM = params.DIM;
  float DRAP = params.DRAP;

  float d_phip = (2.0 * M_PI) / float(DIM_PHIP);

  for (size_t ivar = 0; ivar < 4; ivar++)
  {
    for (size_t is = 0; is < DIM; is++) //the column packed index for x, y and z
    {
      size_t ix = (is % DIM_X);
      size_t iy = ((is / DIM_X) % DIM_Y);
      size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

      for (size_t irap = 0; irap < DIM_RAP; irap++)
      {
        for (size_t iphip = 0; iphip < DIM_PHIP; iphip++)
        {
          //rather than gauss quadrature, just doing a elementary Riemann sum here; check convergence!
          // T^(mu,nu) = int deta int dphip G^(mu,nu)
          baryonCurrent[ivar][is] += shiftedChargeDensity[is][irap][iphip] * hypertrigTable[ivar][irap][iphip][ieta];
        }
      }
      baryonCurrent[ivar][is] = baryonCurrent[ivar][is] * DRAP * d_phip; //multiply by common differential factor once
    }
  }
}

void solveEigenSystem(float **stressTensor, float *energyDensity, float **flowVelocity, const parameters & params)
{
  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;
  float TAU = params.TAU;

  float tolerance = 1.0e-7;

  float gamma_max = 100.0; // the maximum allowed flow boost factor when searching for eigenvectors

  //#pragma omp parallel for  // would be good to parallelize this loop, but it's going to take some careful work (and note that GSL might already have some multithreading of its own)
  for (size_t is = 0; is < DIM; is++)
  {
    gsl_matrix *Tmunu; //T^(mu,nu) with two contravariant indices; we need to lower an index
    //using the metric to find the eigenvectors of T^(mu)_(nu) with one contravariant and one contravariant index
    Tmunu = gsl_matrix_alloc(4,4);
    gsl_matrix *gmunu;
    gmunu = gsl_matrix_alloc(4,4);
    gsl_matrix_complex *eigen_vectors;
    eigen_vectors = gsl_matrix_complex_alloc(4,4);
    gsl_vector_complex *eigen_values;
    eigen_values = gsl_vector_complex_alloc(4);

    //set the values of the energy momentum tensor

    //try adding a small value everywhere to T^\tau\tau make flow velocity -> 0 in dilute regions
    //gsl_matrix_set(Tmunu, 0, 0, stressTensor[0][is] + tolerance); //tau,tau
    gsl_matrix_set(Tmunu, 0, 0, stressTensor[0][is]); //tau,tau
    gsl_matrix_set(Tmunu, 0, 1, stressTensor[1][is]); //tau,x
    gsl_matrix_set(Tmunu, 0, 2, stressTensor[2][is]); //tau,y
    gsl_matrix_set(Tmunu, 0, 3, stressTensor[3][is]); //tau,eta
    gsl_matrix_set(Tmunu, 1, 1, stressTensor[4][is]); //x,x
    gsl_matrix_set(Tmunu, 1, 2, stressTensor[5][is]); //x,y
    gsl_matrix_set(Tmunu, 1, 3, stressTensor[6][is]); //x,eta
    gsl_matrix_set(Tmunu, 2, 2, stressTensor[7][is]); //y,y
    gsl_matrix_set(Tmunu, 2, 3, stressTensor[8][is]); //y,eta
    gsl_matrix_set(Tmunu, 3, 3, stressTensor[9][is]); //eta,eta
    gsl_matrix_set(Tmunu, 1, 0, stressTensor[1][is]); //x,tau
    gsl_matrix_set(Tmunu, 2, 0, stressTensor[2][is]); //y,tau
    gsl_matrix_set(Tmunu, 3, 0, stressTensor[3][is]); //eta,tau
    gsl_matrix_set(Tmunu, 2, 1, stressTensor[5][is]); //y,x
    gsl_matrix_set(Tmunu, 3, 1, stressTensor[6][is]); //eta,x
    gsl_matrix_set(Tmunu, 3, 2, stressTensor[8][is]); //eta,y

    //set the values of the "metric"; not really the metric, but the numerical constants
    //which are multiplied by the elements of T^(mu,nu) to get the values of T^(mu)_(nu)
    //note factors of TAU appropriate for milne coordinates g_(mu.nu) = diag(1,-1,-1,-TAU^2)
    gsl_matrix_set(gmunu, 0, 0, 1.0); //tau,tau
    gsl_matrix_set(gmunu, 0, 1, -1.0); //tau,x
    gsl_matrix_set(gmunu, 0, 2, -1.0); //tau,y
    gsl_matrix_set(gmunu, 0, 3, -1.0*TAU*TAU); //tau,eta
    gsl_matrix_set(gmunu, 1, 0, 1.0); //x,tau
    gsl_matrix_set(gmunu, 1, 1, -1.0); //x,x
    gsl_matrix_set(gmunu, 1, 2, -1.0); //x,y
    gsl_matrix_set(gmunu, 1, 3, -1.0*TAU*TAU); //x,eta
    gsl_matrix_set(gmunu, 2, 0, 1.0); //y,tau
    gsl_matrix_set(gmunu, 2, 1, -1.0); //y,x
    gsl_matrix_set(gmunu, 2, 2, -1.0); //y,y
    gsl_matrix_set(gmunu, 2, 3, -1.0*TAU*TAU); //y,eta
    gsl_matrix_set(gmunu, 3, 0, 1.0); //eta,tau
    gsl_matrix_set(gmunu, 3, 1, -1.0); //eta,x
    gsl_matrix_set(gmunu, 3, 2, -1.0); //eta,y
    gsl_matrix_set(gmunu, 3, 3, -1.0*TAU*TAU); //eta,eta
    //lower one index of the stress tensor; save it to the same matrix to save memory
    gsl_matrix_mul_elements(Tmunu, gmunu); //result stored in Tmunu !this multiplies element-wise, not ordinary matrix multiplication!
    gsl_eigen_nonsymmv_workspace *eigen_workspace;
    eigen_workspace = gsl_eigen_nonsymmv_alloc(4);
    gsl_eigen_nonsymmv(Tmunu, eigen_values, eigen_vectors, eigen_workspace);
    gsl_eigen_nonsymmv_free(eigen_workspace);

    //***does this have a solution for energy density and flow at every point?
    int eigenvalue_exists = 0;
    for (int i = 0; i < 4; i++)
    {
      gsl_complex eigenvalue = gsl_vector_complex_get(eigen_values, i);

      //if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
      //eigenvalue condition taken from JF's suggestion, test for robustness in dilute region
      if ( GSL_REAL(eigenvalue) > 0.0 && fabs( GSL_IMAG(eigenvalue) ) < ( fabs(GSL_REAL(eigenvalue)) * 1.0e-30) ) //choose eigenvalue
      {
        gsl_complex v0 = gsl_matrix_complex_get(eigen_vectors, 0 , i);
        gsl_complex v1 = gsl_matrix_complex_get(eigen_vectors, 1 , i);
        gsl_complex v2 = gsl_matrix_complex_get(eigen_vectors, 2 , i);
        gsl_complex v3 = gsl_matrix_complex_get(eigen_vectors, 3 , i);

        if (GSL_IMAG(v0) == 0 && (2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0 - (GSL_REAL(v3) * GSL_REAL(v3) * (TAU * TAU - 1.0) )) > 0) //choose timelike eigenvector
        {
          double minkowskiLength = GSL_REAL(v0)*GSL_REAL(v0) - (GSL_REAL(v1)*GSL_REAL(v1) + GSL_REAL(v2)*GSL_REAL(v2) + TAU*TAU*GSL_REAL(v3)*GSL_REAL(v3));
          double factor = 1.0 / sqrt(minkowskiLength);

          if (GSL_REAL(v0) < 0) factor=-factor;

          //ignore eigenvectors with gamma too large
          if ( (GSL_REAL(v0) * factor) < gamma_max)
          {
            eigenvalue_exists = 1;
            energyDensity[is] = GSL_REAL(eigenvalue);
            flowVelocity[0][is] = GSL_REAL(v0) * factor;
            flowVelocity[1][is] = GSL_REAL(v1) * factor;
            flowVelocity[2][is] = GSL_REAL(v2) * factor;
            flowVelocity[3][is] = GSL_REAL(v3) * factor;
          }
          /*
	        if ( (energyDensity[is] / tolerance < 1.0) && REGULATE)
	         {
	            energyDensity[is] = tolerance;
	            flowVelocity[0][is] = 1.0;
	            flowVelocity[1][is] = 0.0;
	            flowVelocity[2][is] = 0.0;
	            flowVelocity[3][is] = 0.0;
            }


            if (REGULATE)
            {
              //cut out unreasonable values of flow
              if (flowVelocity[1][is] > 10.0) {flowVelocity[1][is] = 10.0;}
              if (flowVelocity[2][is] > 10.0) {flowVelocity[2][is] = 10.0;}
              if (flowVelocity[1][is] < -10.0) {flowVelocity[1][is] = -10.0;}
              if (flowVelocity[2][is] < -10.0) {flowVelocity[2][is] = -10.0;}
              //enforce timelike condition on these cells

            }
            */

        } // if (GSL_IMAG(v0) == 0 && (2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0 - (GSL_REAL(v3) * GSL_REAL(v3) * (TAU * TAU - 1.0) )) > 0) //choose timelike eigenvector
      } // if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
    } //for (int i = 0; i < 4; ...)

    if (eigenvalue_exists == 0)
    {
      //in dilute regions where we can't find a timelike eigenvector, set e = 0, u^t = 1, u^x=u^y=u^n=0
      energyDensity[is] = 1.0e-3;
      flowVelocity[0][is] = 1.0;
      flowVelocity[1][is] = 0.0;
      flowVelocity[2][is] = 0.0;
      flowVelocity[3][is] = 0.0;
    }

  gsl_matrix_free(gmunu);
  gsl_matrix_free(Tmunu);
  gsl_matrix_complex_free(eigen_vectors);
  gsl_vector_complex_free(eigen_values);

  } // for (size_t is; is < DIM; ...)

  //try scaling the flow velocity by a smooth profile which goes to zero after some finite radius
  /*
  if (REGULATE)
  {
    printf("Regulating flow velocity profile in dilute regions \n");
    for (size_t is = 0; is < DIM; is++)
    {
      size_t ix = (is % DIM_X);
      size_t iy = ((is / DIM_X) % DIM_Y);
      size_t ieta = ((is / DIM_X / DIM_Y) % DIM_ETA);

      float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
      float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
      float r = sqrt(x*x + y*y);
      float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

      float R_WIDTH = 1.0;
      float R_FLAT = 7.5;
      float arg = (-1.0) * (r - R_FLAT) * (r - R_FLAT) / (2.0 * R_WIDTH * R_WIDTH);
      arg = arg * THETA_FUNCTION(r - R_FLAT);

      flowVelocity[1][is] = flowVelocity[1][is] * exp(arg);
      flowVelocity[2][is] = flowVelocity[2][is] * exp(arg);

      flowVelocity[0][is] = sqrt( 1 + flowVelocity[1][is]*flowVelocity[1][is] + flowVelocity[2][is]*flowVelocity[2][is] + TAU*TAU*flowVelocity[3][is]*flowVelocity[3][is]);
    }

  }
  */

} //solveEigenSystem()
void calculateBulkPressure(float **stressTensor, float *energyDensity, float *pressure, float *bulkPressure, const parameters & params)
{
  size_t DIM = params.DIM;
  float TAU = params.TAU;
  for (size_t is = 0; is < DIM; is++)
  {
    // PI = -1/3 * (T^(mu)_(mu) - epsilon) - p
    // T^(mu)_(mu) = T^(0,0) - T^(1,1) - T^(2,2) - (TAU^2)T^(3,3)
    float a =  stressTensor[0][is] - stressTensor[4][is] - stressTensor[7][is] - TAU*TAU*stressTensor[9][is];
    bulkPressure[is] = (-1.0/3.0) * (a - energyDensity[is]) - pressure[is];
  }
}
void calculateShearViscTensor(float **stressTensor, float *energyDensity, float **flowVelocity, float *pressure, float *bulkPressure, float **shearTensor, const parameters & params)
{
  size_t DIM = params.DIM;
  float TAU = params.TAU;
  for (size_t is = 0; is < DIM; is++)
  {
    // pi^(mu,nu) = T^(mu,nu) - epsilon * u^(mu)u^(nu) + (P + PI) * (g^(mu,nu) - u^(mu)u^(nu))
    //calculate ten components : upper triangular part
    float b = energyDensity[is] + pressure[is] + bulkPressure[is];
    float c = pressure[is] + bulkPressure[is];
    shearTensor[0][is] = stressTensor[0][is] - flowVelocity[0][is] * flowVelocity[0][is] * b + c; //pi^(tau,tau)
    shearTensor[1][is] = stressTensor[1][is] - flowVelocity[0][is] * flowVelocity[1][is] * b; //pi^(tau,x)
    shearTensor[2][is] = stressTensor[2][is] - flowVelocity[0][is] * flowVelocity[2][is] * b; //pi^(tau,y)
    shearTensor[3][is] = stressTensor[3][is] - flowVelocity[0][is] * flowVelocity[3][is] * b; //pi^(tau,eta)
    shearTensor[4][is] = stressTensor[4][is] - flowVelocity[1][is] * flowVelocity[1][is] * b - c; //pi^(x,x)
    shearTensor[5][is] = stressTensor[5][is] - flowVelocity[1][is] * flowVelocity[2][is] * b; //pi^(x,y)
    shearTensor[6][is] = stressTensor[6][is] - flowVelocity[1][is] * flowVelocity[3][is] * b; //pi^(x,eta)
    shearTensor[7][is] = stressTensor[7][is] - flowVelocity[2][is] * flowVelocity[2][is] * b - c; //pi^(y,y)
    shearTensor[8][is] = stressTensor[8][is] - flowVelocity[2][is] * flowVelocity[3][is] * b; //pi^(y,eta)
    shearTensor[9][is] = stressTensor[9][is] - flowVelocity[3][is] * flowVelocity[3][is] * b - c * (1.0/(TAU*TAU)); //pi^(eta,eta)
  }
}

// n_B = u^(mu)j_(mu)
void calculateBaryonDensity(float *baryonDensity, float **baryonCurrent, float **flowVelocity, const parameters & params)
{
  size_t DIM = params.DIM;
  float TAU = params.TAU;
  for (size_t is = 0; is < DIM; is++)
  {
    baryonDensity[is] = (flowVelocity[0][is] * baryonCurrent[0][is]) - (flowVelocity[1][is] * baryonCurrent[1][is]) - (flowVelocity[2][is] * baryonCurrent[2][is]) - (TAU * TAU * flowVelocity[3][is] * baryonCurrent[3][is]);
  }
}
// V^(mu) = j^(mu) - n_B * u^(mu)
void calculateBaryonDiffusion(float **baryonDiffusion, float **baryonCurrent, float *baryonDensity, float **flowVelocity, const parameters & params)
{
  size_t DIM = params.DIM;
  for (size_t ivar = 0; ivar < 4; ivar++)
  {
    for (size_t is = 0; is < DIM; is++)
    {
      baryonDiffusion[ivar][is] = baryonCurrent[ivar][is] - (baryonDensity[is] * flowVelocity[ivar][is]);
    }
  }
}


void calculateThermalVorticityTensor(float *energyDensity, float **flowVelocity, float **thermalVelocityVector, float **thermalVorticityTensor, const parameters & params)
{

  // omega_{\mu\nu} = 1/2 ( d_{\nu} \beta{\mu} - d{\mu} \beta_{\nu})
  size_t DIM = params.DIM;
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  float dx = params.DX;
  float dy = params.DY;
  float tau = params.TAU;
  for (size_t ix = 1; ix < DIM_X - 1; ix++)
  {
    for (size_t iy = 1; iy < DIM_Y - 1; iy++)
    {
      size_t is = ix + (DIM_X * iy);
      size_t is_px = (ix + 1) + (DIM_X * iy);
      size_t is_mx = (ix - 1) + (DIM_X * iy);
      size_t is_py = ix + (DIM_X * (iy + 1));
      size_t is_my = ix + (DIM_Y * (iy - 1));

      // omega_{\mu\nu} = 1/2 ( d_{\nu} \beta{\mu} - d{\mu} \beta_{\nu})
      //antisymmetric tensor, calculate 6 components (upper triangular)
      //calculate the covariant thermal flow vector \beta_{\mu}

      float u_t = flowVelocity[0][is];
      float u_t_px = flowVelocity[0][is_px];
      float u_t_mx = flowVelocity[0][is_mx];
      float u_t_py = flowVelocity[0][is_py];
      float u_t_my = flowVelocity[0][is_my];

      float u_x = -flowVelocity[1][is];
      float u_x_px = -flowVelocity[1][is_px];
      float u_x_mx = -flowVelocity[1][is_mx];
      float u_x_py = -flowVelocity[1][is_py];
      float u_x_my = -flowVelocity[1][is_my];

      float u_y = -flowVelocity[2][is];
      float u_y_px = -flowVelocity[2][is_px];
      float u_y_mx = -flowVelocity[2][is_mx];
      float u_y_py = -flowVelocity[2][is_py];
      float u_y_my = -flowVelocity[2][is_my];

      float u_n = 0.0;
      float u_n_px = 0.0;
      float u_n_mx = 0.0;
      float u_n_py = 0.0;
      float u_n_my = 0.0;

      float T = temperatureFromEnergyDensity(energyDensity[is]);
      float T_px = temperatureFromEnergyDensity(energyDensity[is_px]);
      float T_mx = temperatureFromEnergyDensity(energyDensity[is_mx]);
      float T_py = temperatureFromEnergyDensity(energyDensity[is_py]);
      float T_my = temperatureFromEnergyDensity(energyDensity[is_my]);

      float beta_t = u_t / T;
      float beta_x = u_x / T;
      float beta_y = u_y / T;
      float beta_n = u_n / T;

      thermalVelocityVector[0][is] = beta_t;
      thermalVelocityVector[1][is] = beta_x;
      thermalVelocityVector[2][is] = beta_y;
      thermalVelocityVector[3][is] = beta_n;

      float beta_t_px = u_t_px / T_px ;
      float beta_x_px = u_x_px / T_px ;
      float beta_y_px = u_y_px / T_px ;
      float beta_n_px = u_n_px / T_px ;

      float beta_t_mx = u_t_mx / T_mx;
      float beta_x_mx = u_x_mx / T_mx;
      float beta_y_mx = u_y_mx / T_mx;
      float beta_n_mx = u_n_mx / T_mx;

      float beta_t_py = u_t_py / T_py;
      float beta_x_py = u_x_py / T_py;
      float beta_y_py = u_y_py / T_py;
      float beta_n_py = u_n_py / T_py;

      float beta_t_my = u_t_my / T_my;
      float beta_x_my = u_x_my / T_my;
      float beta_y_my = u_y_my / T_my;
      float beta_n_my = u_n_my / T_my;

      //NOTE NEED TO FIX THESE, for the temporal derivatives we need to evolve for another time step
      thermalVorticityTensor[0][is] = 0.0; // w_tx
      thermalVorticityTensor[1][is] = 0.0; // w_ty
      thermalVorticityTensor[2][is] = 0.0; // w_tn
      thermalVorticityTensor[3][is] = 0.5 * ( (beta_x_py - beta_x_my) / 2.0 / dy - (beta_y_px - beta_y_mx) / 2.0 / dx ); // w_xy
      thermalVorticityTensor[4][is] = 0.0; // w_xn
      thermalVorticityTensor[5][is] = 0.0; // w_yn
    }
  }
}
