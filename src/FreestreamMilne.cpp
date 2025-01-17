//This file contains a wrapper class for freestream-milne

#include "Parameter.h"
#include "FreestreamMilne.h"
#include "EquationOfState.h"
#include "FileIO.h"
#include "FreeStream.h"
#include "InitialConditions.h"
#include "HydroValidity.h"
#include "LandauMatch.h"
#include "Memoryf.h"
#include "FSConfig.h"
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

FREESTREAMMILNE::FREESTREAMMILNE() {

}

FREESTREAMMILNE::~FREESTREAMMILNE() {
}

//use this function to initialize energy density within JETSCAPE
void FREESTREAMMILNE::initialize_from_vector(std::vector<float> energy_density_in) {
  init_energy_density = energy_density_in;
}

//use this function to return final hydro variables as vectors within JETSCAPE
void FREESTREAMMILNE::output_to_vectors(std::vector<double> &energy_density_out,
                                        std::vector<double> &pressure_out,
                                        std::vector<double> &ut_out,
                                        std::vector<double> &ux_out,
                                        std::vector<double> &uy_out,
                                        std::vector<double> &un_out,
                                        std::vector<double> &pitt_out,
                                        std::vector<double> &pitx_out,
                                        std::vector<double> &pity_out,
                                        std::vector<double> &pitn_out,
                                        std::vector<double> &pixx_out,
                                        std::vector<double> &pixy_out,
                                        std::vector<double> &pixn_out,
                                        std::vector<double> &piyy_out,
                                        std::vector<double> &piyn_out,
                                        std::vector<double> &pinn_out,
                                        std::vector<double> &Pi_out) {
  energy_density_out = final_energy_density;
  pressure_out = final_pressure;
  ut_out = final_ut;
  ux_out = final_ux;
  uy_out = final_uy;
  un_out = final_un;
  pitt_out = final_pitt;
  pitx_out = final_pitx;
  pity_out = final_pity;
  pitn_out = final_pitn;
  pixx_out = final_pixx;
  pixy_out = final_pixy;
  pixn_out = final_pixn;
  piyy_out = final_piyy;
  piyn_out = final_piyn;
  pinn_out = final_pinn;
  Pi_out = final_Pi;
}

//where the magic happens
int FREESTREAMMILNE::run_freestream_milne(const char * inputPath) {

float hbarc = 0.197326938;

if(PRINT_SCREEN) printf("Welcome to freestream-milne\n");

//declare parameter struct
struct parameters params;

//set default parameters in case of missing input file
params.OUTPUTFORMAT = 2;
params.BARYON = 0;
params.IC_ENERGY = 5;
params.IC_BARYON = 1;
params.ETA_WIDTH = 0.5;
params.ETA_FLAT = 0.5;
params.SIGMA = 1.0;
params.SIGMA_B = 1.0;
params.DIM_X = 101;
params.DIM_Y = 101;
params.DIM_ETA = 1;
params.DIM_RAP = 1;
params.DIM_PHIP = 51;
params.DX = 0.1;
params.DY = 0.1;
params.DETA = 0.1;
params.DRAP = 0.2;
params.DTAU = 0.5;
params.TAU0 = 0.1;
params.EOS_TYPE = 1;
params.E_FREEZE = 1.7;
params.VISCOUS_MATCHING = 1;
params.E_DEP_FS = 0;

//read in chosen parameters from input file if such a file exists
readInParameters(params, inputPath);

//define some useful combinations
params.DIM = params.DIM_X * params.DIM_Y * params.DIM_ETA;
params.TAU = params.TAU0 + params.DTAU;

size_t DIM_X = params.DIM_X;
size_t DIM_Y = params.DIM_Y;
size_t DIM_ETA = params.DIM_ETA;

if(PRINT_SCREEN)
  {
    printf("Parameters are ...\n");
    printf("(DIM_X, DIM_Y, DIM_ETA, DIM_PHIP, DIM_RAP) = (%zu, %zu, %zu, %zu, %zu)\n", params.DIM_X, params.DIM_Y, params.DIM_ETA, params.DIM_PHIP, params.DIM_RAP);
    printf("(DX, DY, DETA, DTAU) = (%.2f fm, %.2f fm, %.2f, %.2f fm/c)\n", params.DX, params.DY, params.DETA, params.DTAU);
    printf("TAU0 = %.2f fm/c\n", params.TAU0);
    printf("SIGMA = %.2f \n", params.SIGMA);
    printf("E_FREEZE = %.3f GeV / fm^3 \n", params.E_FREEZE);
    if (params.VISCOUS_MATCHING) printf("Will match to hydro including viscous part of Tmunu \n");
    else printf("Will match to hydro with ideal part of Tmunu \n");
    if (params.EOS_TYPE == 1) printf("Using EoS : Conformal \n");
    else if (params.EOS_TYPE == 2) printf("Using EoS : Wuppertal-Budhapest \n");
    else if (params.EOS_TYPE == 3) printf("Using EoS : Lattice QCD + HRG matched.\n");

    if (params.E_DEP_FS == 1)
    {
      printf("Freestreaming time dependent on transverse energy of profile \n");
      printf("tau_fs = tau_R * (e_T / e_R) ^ alpha \n");
      printf("e_R = %f GeV/fm \n", params.E_R);
      printf("tau_R = %f fm \n", params.TAU_R);
    }
  }
//allocate and initialize memory
if (PRINT_SCREEN) printf("Allocating memory\n");

//the initial energy density spatial profile
float *initialEnergyDensity;
initialEnergyDensity = (float *)calloc(params.DIM, sizeof(float));

//the initial baryon density spatial profile
float *initialChargeDensity = nullptr;
if(params.BARYON) initialChargeDensity = (float *)calloc(params.DIM, sizeof(float));

//the initial density G(tilde)^(tau,tau) at time tau_0
float **density = nullptr;
density = calloc2dArrayf(density, params.DIM, params.DIM_RAP); // function of x,y,eta and rapidity

//the initial density J(tilde)^(tau) at time tau_0
float **chargeDensity = nullptr;
if(params.BARYON) chargeDensity = calloc2dArrayf(chargeDensity, params.DIM, params.DIM_RAP); // function of x,y,eta and rapidity

//initialize energy density

//define a lower bound on energy density for all cells to regulate numerical noise in flow velocity in dilute regions
float lower_tolerance = 1.0e-7;

if (PRINT_SCREEN) printf("setting initial conditions on energy density : ");
if (params.IC_ENERGY == 1)
{
  initializeEllipticalGauss(initialEnergyDensity, 15.0, 15.0, 15.0, params);
  if(PRINT_SCREEN) printf("Smooth Oblate Gaussian \n");
}
else if (params.IC_ENERGY == 2)
{
  initializeEllipticalMCGauss(initialEnergyDensity, 15.0, 15.0, 15.0, params);
  if(PRINT_SCREEN) printf("Fluctuating Oblate Gaussian \n");
}
else if (params.IC_ENERGY == 3)
{
  readDensityFile(initialEnergyDensity, (char *)"initial_profiles/e", params);
  if(PRINT_SCREEN) printf("Reading from energy density file in initial_profiles/ \n");
}
else if (params.IC_ENERGY == 4)
{
  readEnergyDensityCPUVH(initialEnergyDensity, params);
  if(PRINT_SCREEN) printf("Reading from energy density file in initial_profiles/ \n");
}
else if (params.IC_ENERGY == 5)
{
  //read in initial energy density using the initiliaze_from_vector() function
  //note that this is not safe - if one passes an empty vector it will not throw an error
  //converting units of energy density from GeV / fm^3 to fm^(-4)
  if(PRINT_SCREEN) printf("Reading energy density from initial energy density vector\n");

  //check that the vector has the same dimensions as defined in input file!
  if ( params.DIM != init_energy_density.size() )
    {
      printf("Grid dimension of input vector does not match input file! \n");
      printf("Vector size : %zu , DIM = %zu \n", init_energy_density.size(), params.DIM);
      exit(-1);
    }
  //rescale initial distribution
  float rescale = 1.0;
  for (size_t i = 0; i < params.DIM; i++) initialEnergyDensity[i] = init_energy_density[i] * rescale / (float)hbarc + lower_tolerance;
  //for (size_t i = 0; i < params.DIM; i++) initialEnergyDensity[i] = init_energy_density[i] / (float)hbarc;
  //just doing this here for testing - try increasing normalization of initial distribution to improve stability
  //for (size_t i = 0; i < params.DIM; i++) initialEnergyDensity[i] = init_energy_density[i];
}
else if (params.IC_ENERGY == 6)
{
  printf("Initializing with 2 MC Gaussians \n");
  initialize2Gaussians(initialEnergyDensity, 1.0, 1.0, 1.0, params);
}
else
{
  printf("Not a valid initial Condition... Goodbye\n");
  return 0;
}

if (params.BARYON)
{
  //initialize baryon density
  if (PRINT_SCREEN) printf("setting initial conditions on baryon density : ");
  if (params.IC_BARYON == 1)
  {
    initializeEllipticalGauss(initialChargeDensity, 2.0, 3.0, 1.0, params);
    if (PRINT_SCREEN) printf("Smooth Oblate Gaussian \n");
  }
  else if (params.IC_BARYON == 2)
  {
    initializeEllipticalMCGauss(initialChargeDensity, 2.0, 3.0, 1.0, params);
    if (PRINT_SCREEN) printf("Fluctuating Oblate Gaussian \n");
  }
  else if (params.IC_BARYON == 3)
  {
    readDensityFile(initialChargeDensity, (char *)"initial_profiles/nB", params);
    if (PRINT_SCREEN) printf("Reading from baryon density file in initial_profiles/ \n");
  }
  else
  {
    printf("Not a valid initial Condition... Goodbye\n");
    return 0;
  }
}

//write initial energy density and baryon density to file
writeScalarToFile(initialEnergyDensity, (char *)"initial_e", params);
if (params.BARYON) writeScalarToFile(initialChargeDensity, (char *)"initial_nB", params);
writeScalarToFileProjection(initialEnergyDensity, (char *)"initial_e_projection", params);
if (params.BARYON) writeScalarToFileProjection(initialChargeDensity, (char *)"initial_nB_projection", params);

//calculate transverse energy dependent freestreaming time
if (params.E_DEP_FS == 1)
{
  printf("Calculating transverse energy dependent freestreaming time... \n");
  float tau_fs = getEnergyDependentTau(initialEnergyDensity, params);
  printf("Updating to DTAU = %f \n", tau_fs);
  params.DTAU = tau_fs;
  params.TAU = params.TAU0 + params.DTAU; //update Landau Matching Time
}

//set the value of the Landau matching time stored in class
tau_LandauMatch = params.TAU;

//now reset the number of points in phip based on freestreaming time
float min_dx_dy = min(params.DX, params.DY);
size_t nphip = size_t( ceil( (2.0 * M_PI * params.DTAU) / min_dx_dy ) );
if (nphip > params.DIM_PHIP) printf("Updating number of points in phi_p to %zu based on arc length \n", nphip);
params.DIM_PHIP = max(params.DIM_PHIP, nphip);


/////////////////////////////BEGIN TESTING FOR JETSCAPE//////////////////////////////
//make a toy plot of 1/tau * initial energy density to compare 2+1D freestreaming with only longitudinal (bjorken) dilution
float *scaledEnergyDensity = nullptr;
if (TEST_INTERPOL)
{
  printf("Calculating 1 / tau scaled profile for testing \n");
  scaledEnergyDensity = (float *)calloc(params.DIM, sizeof(float));
  for (size_t is = 0; is < params.DIM; is++) scaledEnergyDensity[is] = initialEnergyDensity[is] * (params.TAU0 / params.TAU);
  writeScalarToFile(scaledEnergyDensity, (char *)"scaled_e", params);
  writeScalarToFileProjection(scaledEnergyDensity, (char *)"scaled_e_projection", params);
}
/////////////////////////////END TESTING FOR JETSCAPE//////////////////////////////


//calculate total energy to check convergence
float totalEnergy = 0.0;
for (size_t is = 0; is < params.DIM; is++) totalEnergy += initialEnergyDensity[is];
if (params.DIM_ETA > 1) totalEnergy *= (params.TAU0 * params.DX * params.DY * params.DETA);
else totalEnergy *= (params.DX * params.DY);
printf("Total energy on grid before streaming : %f \n", totalEnergy);
//convert the energy density profile into the initial density profile to be streamed and free memory
convertInitialDensity(initialEnergyDensity, density, params);

if (!TEST_INTERPOL) free(initialEnergyDensity);
//convert the baryon density profile into the initial baryon density profile to be streamed and free memory
if (params.BARYON) convertInitialChargeDensity(initialChargeDensity, chargeDensity, params);
if (params.BARYON) free(initialChargeDensity);

//the shifted energy density profile G^(tau,tau) at time tau
float ***shiftedDensity = nullptr;
shiftedDensity = calloc3dArrayf(shiftedDensity, params.DIM, params.DIM_RAP, params.DIM_PHIP);

//the shifted baryon density profile J^(tau) at time tau
float ***shiftedChargeDensity = nullptr;
if(params.BARYON) shiftedChargeDensity = calloc3dArrayf(shiftedChargeDensity, params.DIM, params.DIM_RAP, params.DIM_PHIP);

//perform the free streaming time-update step and free up memory
//pretabulate trig and hypertrig functions before this step to save time?
if (PRINT_SCREEN) printf("performing the free streaming\n");

double sec = 0.0;
#ifdef _OPENMP
sec = omp_get_wtime();
#endif

///////////  BEGIN LOOP OVER TIME STEPS HERE ////////////////////////
////////// MOVE DECLARATIONS AND ALLOCATION OUTSIDE LOOP? //////////
freeStream(density, shiftedDensity, params);

//#pragma acc update host(shiftedDensity)
free2dArrayf(density, params.DIM);
if (params.BARYON) freeStream(chargeDensity, shiftedChargeDensity, params);
//#if(params.BARYON) pragma acc update host(shiftedChargeDensity)
if (params.BARYON) free2dArrayf(chargeDensity, params.DIM);

#ifdef _OPENMP
sec = omp_get_wtime() - sec;
#endif
if (PRINT_SCREEN) printf("Free streaming took %f seconds\n", sec);

//Landau matching to find the components of energy-momentum tensor
if (PRINT_SCREEN) printf("Landau matching to find hydrodynamic variables\n");

//the ten independent components of the stress tensor
float **stressTensor = nullptr;
stressTensor = calloc2dArrayf(stressTensor, 10, params.DIM);

//the four independent components of baryon current four-vector
float **baryonCurrent = nullptr;
if(params.BARYON) baryonCurrent = calloc2dArrayf(baryonCurrent, 4, params.DIM);

//a table containing 10 rows for 10 independent combinations of p_(mu)p_(nu)
//hypertrig table depends on TAU, so need to keep this inside loop
float ****hypertrigTable = nullptr;
hypertrigTable = calloc4dArrayf(hypertrigTable, 10, params.DIM_RAP, params.DIM_PHIP, params.DIM_ETA); //depends on eta because we have function of eta - y

if (PRINT_SCREEN) printf("calculating hypertrig table\n");
calculateHypertrigTable(hypertrigTable, params);

//calculate the ten independent components of the stress tensor by integrating over rapidity and phi_p
if (PRINT_SCREEN) printf("calculating independent components of stress tensor\n");
#ifdef _OPENMP
sec = omp_get_wtime();
#endif
calculateStressTensor(stressTensor, shiftedDensity, hypertrigTable, params);
free3dArrayf(shiftedDensity, params.DIM, params.DIM_RAP);
#ifdef _OPENMP
sec = omp_get_wtime() - sec;
#endif
if (PRINT_SCREEN) printf("calculating stress tensor took %f seconds\n", sec);

if (params.BARYON)
{
  //calculate the four independent components of the baryon current by integrating over rapidity and phi_p
  if (PRINT_SCREEN) printf("calculating independent components of baryon current\n");
  #ifdef _OPENMP
  sec = omp_get_wtime();
  #endif
  calculateBaryonCurrent(baryonCurrent, shiftedChargeDensity, hypertrigTable, params);
  free3dArrayf(shiftedChargeDensity, params.DIM, params.DIM_RAP);
  #ifdef _OPENMP
  sec = omp_get_wtime() - sec;
  #endif
  if (PRINT_SCREEN) printf("calculating baryon current took %f seconds\n", sec);
}

//done with hypertrig table as well
free4dArrayf(hypertrigTable, 10, params.DIM_RAP, params.DIM_PHIP);

//variables to store the hydrodynamic variables after the Landau matching is performed
//the energy density
float *energyDensity;
energyDensity = (float *)calloc(params.DIM, sizeof(float));

//the baryon density
float *baryonDensity = nullptr;
if(params.BARYON) baryonDensity = (float *)calloc(params.DIM, sizeof(float));

//the flow velocity
float **flowVelocity = nullptr;
flowVelocity = calloc2dArrayf(flowVelocity, 4, params.DIM);

//the pressure
float *pressure;
pressure = (float *)calloc(params.DIM, sizeof(float));

//the bulk pressure Pi
float *bulkPressure;
bulkPressure = (float *)calloc(params.DIM, sizeof(float));

//the shear stress tensor
float **shearTensor = nullptr;
shearTensor = calloc2dArrayf(shearTensor, 10, params.DIM); //calculate 10 components, can check tracelessness/orthogonality for accuracy

//the baryon diffusion current vector
float **baryonDiffusion = nullptr;
if(params.BARYON) baryonDiffusion = calloc2dArrayf(baryonDiffusion, 4, params.DIM);

//solve the eigenvalue problem for the energy density and flow velocity
if (PRINT_SCREEN) printf("solving eigenvalue problem for energy density and flow velocity\n");
#ifdef _OPENMP
sec = omp_get_wtime();
#endif
solveEigenSystem(stressTensor, energyDensity, flowVelocity, params);
#ifdef _OPENMP
sec = omp_get_wtime() - sec;
#endif
if (PRINT_SCREEN) printf("solving eigenvalue problem took %f seconds\n", sec);

if (params.BARYON)
{
  //calculate baryon density and diffusion current
  if (PRINT_SCREEN) printf("calculating baryon density and diffusion current \n");
  #ifdef _OPENMP
  sec = omp_get_wtime();
  #endif
  calculateBaryonDensity(baryonDensity, baryonCurrent, flowVelocity, params);
  calculateBaryonDiffusion(baryonDiffusion, baryonCurrent, baryonDensity, flowVelocity, params);
  #ifdef _OPENMP
  sec = omp_get_wtime() - sec;
  #endif
  if (PRINT_SCREEN) printf("calculating baryon density and diffusion current took %f seconds\n", sec);
}

calculatePressure(energyDensity, baryonDensity, pressure, params);

//viscous currents are nonzero if VISCOUS_MATCHING is on
if (params.VISCOUS_MATCHING)
{
  calculateBulkPressure(stressTensor, energyDensity, pressure, bulkPressure, params);
  calculateShearViscTensor(stressTensor, energyDensity, flowVelocity, pressure, bulkPressure, shearTensor, params);
}


/////////////////////////////BEGIN TESTING FOR JETSCAPE//////////////////////////////
if (TEST_INTERPOL)
{
  printf("approximating energy density profile at intermed. times by interpolating between initial and final profiles \n");
  float TAU = params.TAU;
  float TAU0 = params.TAU0;
  float DTAU = params.DTAU;
  float tau_i = TAU0 + (DTAU / 2.0); //some intermediate time
  float c_1 = (TAU0 / tau_i);
  float c_2 = (tau_i - TAU0) / DTAU / tau_i;
  for (size_t is = 0; is < params.DIM; is++) scaledEnergyDensity[is] = c_1 * initialEnergyDensity[is] + c_2 * ((TAU * energyDensity[is] - TAU0 * initialEnergyDensity[is]));
  writeScalarToFile(scaledEnergyDensity, (char *)"tau_interpolated_e", params);
  writeScalarToFileProjection(scaledEnergyDensity, (char *)"tau_interpolated_e_projection", params);
  free(scaledEnergyDensity);
}
/////////////////////////////END TESTING FOR JETSCAPE//////////////////////////////

//write the next entry in hdf5 spacetime evolution file
//WriteHistoryJETSCAPE();

///////////  END LOOP OVER TIME STEPS HERE ////////////////////////

float totalEnergyAfter = 0.0;
for (size_t is = 0; is < params.DIM; is++) totalEnergyAfter += stressTensor[0][is];
if (params.DIM_ETA > 1) totalEnergyAfter *= (params.TAU * params.DX * params.DY * params.DETA);
else totalEnergyAfter *= (params.TAU * params.DX * params.DY);
printf("Total energy on grid after streaming : %f \n", totalEnergyAfter);

//check which fraction of total energy lies within freezeout surface, which lies in 'corona'
float totalEnergyInsideHypersurf = 0.0;
for (size_t is = 0; is < params.DIM; is++)
{
  //if ( (energyDensity[is] * hbarc) > params.E_FREEZE) totalEnergyInsideHypersurf += energyDensity[is];
  if ( (energyDensity[is] * hbarc) > params.E_FREEZE) totalEnergyInsideHypersurf += stressTensor[0][is];
}
if (params.DIM_ETA > 1) totalEnergyInsideHypersurf *= (params.TAU * params.DX * params.DY * params.DETA);
else totalEnergyInsideHypersurf *= (params.TAU * params.DX * params.DY);
printf("Fraction of energy contained in Freezeout Hypersurface : %f \n", totalEnergyInsideHypersurf / totalEnergyAfter);

//write these to file
ofstream energy_file;
energy_file.open ("output/energy_inside_switch_surf.dat");
//energy_file << totalEnergyAfter * hbarc<< endl;
energy_file << totalEnergyInsideHypersurf * hbarc << endl;
energy_file.close();

//////////////////////////////////HYDRO VALIDITY//////////////////////////////////
//bulk inv reynolds #
float *R_Pi_Inv;
R_Pi_Inv = (float *)calloc(params.DIM, sizeof(float));
//shear inv reynolds #
float *R_pimunu_Inv;
R_pimunu_Inv = (float *)calloc(params.DIM, sizeof(float));
calculateBulkInvReynolds(pressure, bulkPressure, R_Pi_Inv, params);
calculateShearInvReynolds(energyDensity, pressure, shearTensor, R_pimunu_Inv, params);
writeScalarToFileProjection(R_Pi_Inv, (char *)"R_Pi_Inv_projection", params);
writeScalarToFileProjection(R_pimunu_Inv, (char *)"R_pimunu_Inv_projection", params);
size_t ctr = ((DIM_X - 1) / 2) + (DIM_X * ((DIM_Y - 1) / 2)) + (DIM_X * DIM_Y * ((DIM_ETA - 1) / 2));
printf("R_Pi_Inv at center : %f \n", R_Pi_Inv[ctr]);
printf("R_pimunu_Inv at center : %f \n", R_pimunu_Inv[ctr]);
free(R_Pi_Inv);
free(R_pimunu_Inv);


//check transversality and tracelesness
//components of pi^munu u_mu
float *pi_dot_u_tau;
pi_dot_u_tau = (float *)calloc(params.DIM, sizeof(float));
float *pi_dot_u_x;
pi_dot_u_x = (float *)calloc(params.DIM, sizeof(float));
float *pi_dot_u_y;
pi_dot_u_y = (float *)calloc(params.DIM, sizeof(float));
float *pi_dot_u_eta;
pi_dot_u_eta = (float *)calloc(params.DIM, sizeof(float));

calculate_pi_dot_u(flowVelocity, shearTensor, pi_dot_u_tau, pi_dot_u_x, pi_dot_u_y, pi_dot_u_eta, params);
writeScalarToFileProjection(pi_dot_u_tau, (char *)"pi_dot_u_tau_projection", params);
writeScalarToFileProjection(pi_dot_u_x, (char *)"pi_dot_u_x_projection", params);
writeScalarToFileProjection(pi_dot_u_y, (char *)"pi_dot_u_y_projection", params);
writeScalarToFileProjection(pi_dot_u_eta, (char *)"pi_dot_u_eta_projection", params);

free(pi_dot_u_tau);
free(pi_dot_u_x);
free(pi_dot_u_y);
free(pi_dot_u_eta);

//trace of pi^munu , pi^mu_mu
float *pi_mu_mu;
pi_mu_mu = (float *)calloc(params.DIM, sizeof(float));

calculate_pi_mu_mu(shearTensor, pi_mu_mu, params);
writeScalarToFileProjection(pi_mu_mu, (char *)"pi_mu_mu_projection", params);
free(pi_mu_mu);

//////////////////////////////////HYDRO VALIDITY//////////////////////////////////



////////////////////////////////////VORTICITY/////////////////////////////////////

float **thermalVorticityTensor = nullptr;
thermalVorticityTensor = calloc2dArrayf(thermalVorticityTensor, 6, params.DIM);
float **thermalVelocityVector = nullptr;
thermalVelocityVector = calloc2dArrayf(thermalVelocityVector, 4, params.DIM);
calculateThermalVorticityTensor(energyDensity, flowVelocity, thermalVelocityVector, thermalVorticityTensor, params);


if (PRINT_SCREEN) printf("writing hydro variables\n");

writeScalarToFile(energyDensity, (char *)"e", params);
writeScalarToFile(pressure, (char *)"p", params);
writeScalarToFile(bulkPressure, (char *)"bulk_PI", params);
writeScalarToFileProjection(energyDensity, (char *)"e_projection", params);
writeScalarToFileProjection(pressure, (char *)"p_projection", params);
writeScalarToFileProjection(bulkPressure, (char *)"bulk_PI_projection", params);

writeVectorToFile(flowVelocity, (char *)"ut", 0, params);
writeVectorToFile(flowVelocity, (char *)"ux", 1, params);
writeVectorToFile(flowVelocity, (char *)"uy", 2,params);
writeVectorToFile(flowVelocity, (char *)"un", 3,params);

writeVectorToFileProjection(flowVelocity, (char *)"ut_projection", 0,params);
writeVectorToFileProjection(flowVelocity, (char *)"ux_projection", 1,params);
writeVectorToFileProjection(flowVelocity, (char *)"uy_projection", 2,params);
writeVectorToFileProjection(flowVelocity, (char *)"un_projection", 3,params);


writeVectorToFile(shearTensor, (char *)"pitt", 0,params);
writeVectorToFile(shearTensor, (char *)"pitx", 1,params);
writeVectorToFile(shearTensor, (char *)"pity", 2,params);
writeVectorToFile(shearTensor, (char *)"pitn", 3,params);
writeVectorToFile(shearTensor, (char *)"pixx", 4,params);
writeVectorToFile(shearTensor, (char *)"pixy", 5,params);
writeVectorToFile(shearTensor, (char *)"pixn", 6,params);
writeVectorToFile(shearTensor, (char *)"piyy", 7,params);
writeVectorToFile(shearTensor, (char *)"piyn", 8,params);
writeVectorToFile(shearTensor, (char *)"pinn", 9,params);

writeVectorToFileProjection(shearTensor, (char *)"pitt_projection", 0,params);
writeVectorToFileProjection(shearTensor, (char *)"pitx_projection", 1,params);
writeVectorToFileProjection(shearTensor, (char *)"pity_projection", 2,params);
writeVectorToFileProjection(shearTensor, (char *)"pitn_projection", 3,params);
writeVectorToFileProjection(shearTensor, (char *)"pixx_projection", 4,params);
writeVectorToFileProjection(shearTensor, (char *)"pixy_projection", 5,params);
writeVectorToFileProjection(shearTensor, (char *)"pixn_projection", 6,params);
writeVectorToFileProjection(shearTensor, (char *)"piyy_projection", 7,params);
writeVectorToFileProjection(shearTensor, (char *)"piyn_projection", 8,params);
writeVectorToFileProjection(shearTensor, (char *)"pinn_projection", 9,params);

writeVectorToFileProjection(thermalVelocityVector, (char *)"beta_x_projection", 1, params);
writeVectorToFile(thermalVelocityVector, (char *)"beta_x", 1, params);
writeVectorToFileProjection(thermalVelocityVector, (char *)"beta_y_projection", 2, params);
writeVectorToFile(thermalVelocityVector, (char *)"beta_y", 2, params);

writeVectorToFileProjection(thermalVorticityTensor, (char *)"w_xy_projection", 3, params);
writeVectorToFile(thermalVorticityTensor, (char *)"w_xy", 3, params);

if (params.BARYON)
{
  writeScalarToFile(baryonDensity, (char *)"nB",params);
  writeScalarToFileProjection(baryonDensity, (char *)"nB_projection",params);
  writeVectorToFile(baryonDiffusion, (char *)"V_x", 1,params);
  writeVectorToFile(baryonDiffusion, (char *)"V_y", 2,params);
  writeVectorToFile(baryonDiffusion, (char *)"V_eta", 3,params);
}


//support for JETSCAPE - write hydro variables to vectors
final_energy_density.resize(params.DIM);
final_pressure.resize(params.DIM);
final_ut.resize(params.DIM);
final_ux.resize(params.DIM);
final_uy.resize(params.DIM);
final_un.resize(params.DIM);
final_pitt.resize(params.DIM);
final_pitx.resize(params.DIM);
final_pity.resize(params.DIM);
final_pitn.resize(params.DIM);
final_pixx.resize(params.DIM);
final_pixy.resize(params.DIM);
final_pixn.resize(params.DIM);
final_piyy.resize(params.DIM);
final_piyn.resize(params.DIM);
final_pinn.resize(params.DIM);
final_Pi.resize(params.DIM);

if ( (params.OUTPUTFORMAT == 2) || (params.OUTPUTFORMAT == 3) )
{
  for (size_t is = 0; is < params.DIM; is++)
  {
    //converting back to GeV / fm^3 for use in JETSCAPE
    final_energy_density[is] = (double)energyDensity[is] * hbarc;
    final_pressure[is] = (double)pressure[is] * hbarc;
    final_ut[is] = (double)flowVelocity[0][is];
    final_ux[is] = (double)flowVelocity[1][is];
    final_uy[is] = (double)flowVelocity[2][is];
    final_un[is] = (double)flowVelocity[3][is];
    final_pitt[is] = (double)shearTensor[0][is] * hbarc;
    final_pitx[is] = (double)shearTensor[1][is] * hbarc;
    final_pity[is] = (double)shearTensor[2][is] * hbarc;
    final_pitn[is] = (double)shearTensor[3][is] * hbarc;
    final_pixx[is] = (double)shearTensor[4][is] * hbarc;
    final_pixy[is] = (double)shearTensor[5][is] * hbarc;
    final_pixn[is] = (double)shearTensor[6][is] * hbarc;
    final_piyy[is] = (double)shearTensor[7][is] * hbarc;
    final_piyn[is] = (double)shearTensor[8][is] * hbarc;
    final_pinn[is] = (double)shearTensor[9][is] * hbarc;
    final_Pi[is] = (double)bulkPressure[is] * hbarc;
  }
}

//free the memory
free2dArrayf(stressTensor, 10);
free2dArrayf(thermalVorticityTensor, 6);
free2dArrayf(thermalVelocityVector, 4);
free(energyDensity);
free2dArrayf(flowVelocity, 4);
free(pressure);
free(bulkPressure);
free2dArrayf(shearTensor, 10);

if (params.BARYON)
{
  free2dArrayf(baryonCurrent, 4);
  free(baryonDensity);
  free2dArrayf(baryonDiffusion, 4);
}

if (PRINT_SCREEN) printf("Done... Goodbye!\n");

//change this to return a different int status if something goes wrong?
return 0;
}
