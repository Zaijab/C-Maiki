/* FILE: model.c
 * -------------
 * Doubleended to implement the model found at
 * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1691475/
 * 
 * In the paper there is a set of difference equations
 * The difference equations are implemented in model_updater
 * model_updater is then used in model_solver to generate a sequence of the number of cases
 * The sequence is compared to real data and a least squares error value can be computed
 * The error value solely depends on the parameters for the model
 * Thus GSL will be called to minimize the error function by finding optimal parameters
 * The optimal parameters can then be used in a simulation, be plotted, and documented.
 *
 * The code adapts the model equations by using array indicies to function as subscripts in the model.
 * This means that E_m(t) in math gets translated to E[1][t] in C.
 * However arrays are indexed by natural numbers, thus the following key associates groups to natural numbers.
 * h -> 0
 * m -> 1
 * c -> 2
 * c (the community group) will be 2 for the time being, but the plans are for the model to expand.
 * So there can be more community indicies later. 
 * Thus, the code is designed to allow for c to be any value greater than 2.
 */

#include <math.h> // Using exp()
#include <stdio.h> // Using printf()
#include "model.h"
 
/* FUNCTION: main
 * --------------
 * Entry podouble of the program.
 * To doublestantiate initial values of the parameters then run model_optimizer to mutate them. 
 *
 * PARAMETERS:
 * -----------
 * argc : double : the number of commandline arguments when running the program.
 * argv : char** : array of character arrays, each index holds one commandline argument.
 *
 * RETURN:
 * -------
 * 0 if the function runs without error.
 * 1 if the function runs with error.
 */
double main(double argc, char* argv[]) {
  FILE *data_file = fopen("new_cases.txt", "w+");
  

  double seven_day_avg[SIMULATION_TIME - 7] = {0};

  for(int i = 0; i < SIMULATION_TIME - 7; i++) {
    seven_day_avg[i] = (result[i] + result[i+1] +result[i+2] +result[i+3] +result[i+4] +result[i+5] +result[i+6]) / 7;
  }
    
  for(int i = 0; i < SIMULATION_TIME - 7; i++) {
    fprintf(data_file, "%d %f\n", i, seven_day_avg[i]);
  }

  fclose(data_file);
  
  return 0;
}

/* FUNCTION: model_optimizer
 * -------------------------
 * To be called in main
 * Takes in podoubleers to the parameters of the model, then mutates them to be values that minimize the least squares residual.
 * Calls on model_solver over and over again.
 *
 * PARAMETERS:
 * -----------
 * l -
 * p - 
 * q - 
 * h - 
 * r - 
 */
void model_optimizer(double S[3][SIMULATION_TIME], double E[3][14][SIMULATION_TIME], double I[3][5][SIMULATION_TIME], double R[3][SIMULATION_TIME],
		     double alpha, double beta[3], double epsilon, double rho, double gamma, double kappa, double eta, double delta, double p[14], double q[14], double h[3][5], double lambda[3][SIMULATION_TIME], double r,
		     double result[SIMULATION_TIME]) {
  /* Do litterally nothing relating to optimization atm */
  model_solver(S, E, I, R, alpha, beta, epsilon, rho, gamma, kappa, eta, delta, p, q, h, lambda, r, result);

  for(int i = 0; i < SIMULATION_TIME; i++) {
    //prdoublef("Day %d: %d\n", i, result[i]);
  }
}

/* FUNCTION: model_solver
 * ----------------------
 * Initializes state variables S, E, I, R, and t and takes in the parameters for the model.
 * Will call on model_updater (t-1) times to get (t) many days of simulated data.
 * Returns an array of doubleegers representing the daily counts (the thing we need to match the data).
 * model_optimizer will use this returned array to generate an error value.
 *
 * PARAMETERS:
 * -----------
 * l -
 * p - 
 * q - 
 * h - 
 * r - 
 */

void model_solver(double S[3][SIMULATION_TIME], double E[3][14][SIMULATION_TIME], double I[3][5][SIMULATION_TIME], double R[3][SIMULATION_TIME],
		  double alpha, double beta[3], double epsilon, double rho, double gamma, double kappa, double eta, double delta, double p[14], double q[14], double h[3][5], double lambda[3][SIMULATION_TIME], double r,
		  double result[SIMULATION_TIME]) {

  for(int i = 0; i < SIMULATION_TIME; i++) {
    result[i] = model_updater(S, E, I, R, i, alpha, beta, epsilon, rho, gamma, kappa, eta, delta, p, q, h, lambda, r, result);
  }
}

/* FUNCTION: model_updater
 * -----------------------
 * Implementation of the difference equations. 
 * The function takes in state variables S, E, I, R, and t, along with parameters, l, p, q, h, and r.
 * Using the given state variables and parameters, the model can advance one day.
 *
 * PARAMETERS:
 * -----------
 * S[Group][Day] - Number of succeptible people in a given Group on a given Day
 * E[Group][Stage][Day] - Number of exposed people in a given Group in a given stage on a given Day
 * I[Group][Stage][Day] - Number of infected people in a given Group in a given stage on a given Day
 * R[Group][Day] - Number of recoverd people in a given Group on a given Day
 * t - The given day
 * l -
 * p - 
 * q - 
 * h - 
 * r - 
 *
 * RETURNS:
 * --------
 * The new infected count on that day.
 * The new count is calculated by looking at the positive flux doubleo the managed pool.
 */
double model_updater(double S[3][SIMULATION_TIME], double E[3][14][SIMULATION_TIME], double I[3][5][SIMULATION_TIME], double R[3][SIMULATION_TIME], int t,
		  double alpha, double beta[3], double epsilon, double rho, double gamma, double kappa, double eta, double delta, double p[14], double q[14], double h[3][5], double lambda[3][SIMULATION_TIME], double r,
		  double result[SIMULATION_TIME]) {

  /* Sum certain variables to be used in formulae */
  double exposed_sum_h = 0;
  double infected_sum_h = 0;
  double exposed_sum_m = 0;
  double infected_sum_m = 0;
  double exposed_sum_c = 0;
  double infected_sum_c = 0;
  
  for(int i = 0; i < 14; i++) {
    exposed_sum_h += E[0][i][t];
  }
  
  for(int i = 0; i < 5; i++) {
    infected_sum_h += I[0][i][t];
  }
  
  for(int i = 0; i < 14; i++) {
    exposed_sum_m += E[1][i][t];
  }
  
  for(int i = 0; i < 5; i++) {
    infected_sum_m += I[1][i][t];
  }
  
  for(int i = 0; i < 14; i++) {
    exposed_sum_c += E[2][i][t];
  }

  for(int i = 0; i < 5; i++) {
    infected_sum_c += I[2][i][t];
  }

  /* Travelers
  double Er = 1435;
  double Et = 15984;
  */
  
  /* Set beta depending on the day t */
  double beta_t = 0;

  if(t < 25) {
    beta_t = beta[0];
  } else if(t < 80) {
    beta_t = beta[1];
  } else {
    beta_t = beta[2];
  }
    

  double N_c = S[2][t] + exposed_sum_c + infected_sum_c + R[2][t] + rho * (S[0][t] + exposed_sum_h + infected_sum_h + R[0][t]);
  double N_h = S[0][t] + exposed_sum_h + infected_sum_h + R[0][t] + exposed_sum_m + infected_sum_m;

  
  /* lambda[2][t] = (beta_t * (infected_sum_c + epsilon * exposed_sum_c) +
		  rho * beta_t * (infected_sum_h + epsilon * exposed_sum_h) +
		  gamma * beta_t * epsilon * exposed_sum_m +
		  beta_t * delta * (Er + alpha * Et)) / N_c; */
  
  lambda[2][t] = (beta_t * (infected_sum_c + epsilon * exposed_sum_c) +
		  rho * beta_t * (infected_sum_h + epsilon * exposed_sum_h) +
		  gamma * beta_t * epsilon * exposed_sum_m) / N_c;
  
  lambda[0][t] = rho * lambda[2][t] + eta * beta_t * (infected_sum_h + epsilon * exposed_sum_h + kappa * infected_sum_m) / N_h;


  /* COMMUNITY EQNS */
  S[2][t+1] = exp(-lambda[2][t]) * S[2][t];
  E[2][0][t+1] = (1 - exp(-lambda[2][t])) * S[2][t];
  
  for(int i = 1; i < 14; i++) {
    E[2][i][t+1] = (1 - p[i-1]) * (1 - q[i-1]) * E[2][i-1][t];
  }
  
  for(int i = 0; i < 13; i++) {
    I[2][0][t+1] += p[i]*(1-q[i])*E[2][i][t];
  }
  I[2][1][t+1] = (1-h[2][0]) * I[2][0][t];
  I[2][2][t+1] = (1-h[2][1]) * I[2][1][t] + (1 - r) * (1 - h[2][2]) * I[2][2][t];
  for(int j = 3; j < 5; j++) {
    I[2][j][t] = r * (1-h[2][j-1])*I[2][j-1][t] + (1-r)*(1-h[2][j])*I[2][j][t];
  }
  R[2][t+1] = R[2][t] + r*I[2][4][t]+r*I[1][4][t] + E[2][13][t] + E[1][13][t];
  /* END OF COMMUNITY */

  /* HCW EQNS */
  S[0][t+1] = exp(-lambda[0][t]) * S[0][t];
  E[0][0][t+1] = (1 - exp(-lambda[0][t])) * S[0][t];
  
  for(int i = 1; i < 14; i++) {
    E[0][i][t+1] = (1 - p[i-1]) * E[0][i-1][t];
  }
  
  for(int i = 0; i < 13; i++) {
    I[0][0][t+1] += p[i]*(1-q[i])*E[0][i][t];
  }
  I[0][1][t+1] = (1-h[0][0]) * I[0][0][t];
  I[0][2][t+1] = I[0][1][t] + (1 - r) * (1 - h[0][2]) * I[0][2][t];
  for(int j = 3; j < 5; j++) {
    I[0][j][t] = r * (1-h[0][j-1])*I[0][j-1][t] + (1-r)*(1-h[0][j])*I[0][j][t];
  }
  R[0][t+1] = R[0][t] + r*I[0][4][t]+r*I[0][4][t] + E[0][13][t];
  /* END OF HCW */
  
  /* MANAGED EGNS */ 
  for(int i = 1; i < 14; i++) {
    E[1][i][t+1] = (1 - p[i - 1])*(q[i - 1] * E[2][i - 1][t] + E[1][i - 1][t]);
  }

  for(int i = 0; i < 13; i++) {
    I[1][0][t+1] += p[i]*(q[i]*E[2][i][t] + E[1][i][t]);
  }			  
  I[1][1][t+1] = h[2][0] * I[2][0][t] + h[0][0] * I[0][0][t] + I[1][0][t];
  I[1][2][t + 1] = (h[2][1] * I[2][1][t] + h[0][1] * I[0][1][t] + I[1][1][t])
    + (1 - r)*(h[2][2] * I[2][2][t] + h[0][2] * I[0][2][t] + I[1][2][t]);
  
  for(int j = 3; j < 5; j++) {
    I[1][j][t + 1] = r * (h[2][j - 1] * I[2][j - 1][t] + h[0][j - 1] * I[0][j - 1][t] + I[1][j - 1][t])
      + (1 - r)*(h[2][j] * I[2][j][t] + h[0][j] * I[0][j][t] + I[1][j][t]);
  }
  /* END OF MANAGED */
  /* MANAGED EGNS */ 
  for(int i = 1; i < 14; i++) {
    E[1][i][t+1] = (1 - p[i - 1])*(q[i - 1] * E[2][i - 1][t] + E[1][i - 1][t]);
  }

  for(int i = 0; i < 13; i++) {
    I[1][0][t+1] += p[i]*(q[i]*E[2][i][t] + E[1][i][t]);
  }			  
  I[1][1][t+1] = h[2][0] * I[2][0][t] + h[0][0] * I[0][0][t] + I[1][0][t];
  I[1][2][t + 1] = (h[2][1] * I[2][1][t] + h[0][1] * I[0][1][t] + I[1][1][t])
    + (1 - r)*(h[2][2] * I[2][2][t] + h[0][2] * I[0][2][t] + I[1][2][t]);
  
  for(int j = 3; j < 5; j++) {
    I[1][j][t + 1] = r * (h[2][j - 1] * I[2][j - 1][t] + h[0][j - 1] * I[0][j - 1][t] + I[1][j - 1][t])
      + (1 - r)*(h[2][j] * I[2][j][t] + h[0][j] * I[0][j][t] + I[1][j][t]);
  }
  /* END OF MANAGED */
  double return_value = h[2][0] * I[2][0][t] + h[2][1] * I[2][1][t] +
    h[2][2] * I[2][2][t] + h[2][3] * I[2][3][t] + (1 - r) * h[2][4] * I[2][4][t] +
    h[0][0] * I[0][0][t] + h[0][1] * I[0][1][t] + h[0][2] * I[0][2][t] +
    h[0][3] * I[0][3][t] + (1-r) * h[0][4] * I[0][4][t];

  if(DEBUG) {
    printf("\n\n---------- DAY %d ----------\n", t);
    printf("\tS_c[t] : %f\n", S[2][t]);
    printf("\tE_c[t] : %f\n", exposed_sum_c);
    printf("\tI_c[t] : %f\n", infected_sum_c);
    printf("\tR_c[t] : %f\n", R[2][t]);
    printf("\tE_m[t] : %f\n", exposed_sum_m);
    printf("\tI_m[t] : %f\n", infected_sum_m);
    printf("\tI_c,1[t] : %f\n", I[2][0][t]);
    printf("\tS_h[t] : %f\n", S[0][t]);
    printf("\tE_h[t] : %f\n", exposed_sum_h);
    printf("\tI_h[t] : %f\n", infected_sum_h);
    printf("\tR_h[t] : %f\n", R[0][t]);
    printf("\tl_c[t] : %.12f\n", lambda[2][t]);
    printf("\tl_h[t] : %f\n", lambda[0][t]);
    printf("\tNewly Infected : %f\n", return_value);
    printf("----------------------------\n\n");
  }

  return return_value;
}
