#pragma once

#include "Parameter.h"

void calculateBulkInvReynolds(float *pressure, float *bulkPressure, float *R_Pi_Inv, const parameters & params);
void calculateShearInvReynolds(float *energyDensity, float *pressure, float **shearTensor, float *R_pimunu_Inv, const parameters & params);
void calculate_pi_dot_u(float **flowVelocity, float **shearTensor, float *pi_dot_u_tau, float *pi_dot_u_x, float *pi_dot_u_y, float *pi_dot_u_eta, const parameters & params);
void calculate_pi_mu_mu(float **shearTensor, float *pi_mu_mu, const parameters & params);
