/* FILE: model.c
 * -------------
 * Intended to implement the model found at
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

/* I define a time to run the simulation here. Change the variable then recompile to get different lengths of simulation.
 * The reason is to make the array indicies known at compile time rather than using a variable on the stack / heap.
 * C gets mad if I don't do this.
 */
#define SIMULATION_TIME 100
#define DEBUG true

void model_optimizer(int S[3][SIMULATION_TIME], int E[3][14][SIMULATION_TIME], int I[3][5][SIMULATION_TIME], int R[3][SIMULATION_TIME],
		     double alpha, double beta[3], double epsilon, double rho, double gamma, double kappa, double eta, double delta, double p[14], double q[14], double h[3][5], double lambda[3][SIMULATION_TIME], double r,
		     int result[SIMULATION_TIME]);

void model_solver(int S[3][SIMULATION_TIME], int E[3][14][SIMULATION_TIME], int I[3][5][SIMULATION_TIME], int R[3][SIMULATION_TIME],
		  double alpha, double beta[3], double epsilon, double rho, double gamma, double kappa, double eta, double delta, double p[14], double q[14], double h[3][5], double lambda[3][SIMULATION_TIME], double r,
		  int result[SIMULATION_TIME]);

int model_updater(int S[3][SIMULATION_TIME], int E[3][14][SIMULATION_TIME], int I[3][5][SIMULATION_TIME], int R[3][SIMULATION_TIME], int t,
		  double alpha, double beta[3], double epsilon, double rho, double gamma, double kappa, double eta, double delta, double p[14], double q[14], double h[3][5], double lambda[3][SIMULATION_TIME], double r,
		  int result[SIMULATION_TIME]);
 
/* FUNCTION: main
 * --------------
 * Entry point of the program.
 * To intstantiate initial values of the parameters then run model_optimizer to mutate them. 
 *
 * PARAMETERS:
 * -----------
 * argc : int : the number of commandline arguments when running the program.
 * argv : char** : array of character arrays, each index holds one commandline argument.
 *
 * RETURN:
 * -------
 * 0 if the function runs without error.
 * 1 if the function runs with error.
 */
int main(int argc, char* argv[]) {
  /* Initial Conditions for State Variables */
  int S[3][SIMULATION_TIME] = {0};
  int E[3][14][SIMULATION_TIME] = {0};
  int I[3][5][SIMULATION_TIME] = {0};
  int R[3][SIMULATION_TIME] = {0};
  S[0][0] = 15000;
  S[2][0] = 937711;
  I[2][0][0] = 1;

  /* Initial Conditions for Parameters */
  double alpha = 0.5;
  double beta[3] = {0.4549, 0.1042, 0.2054};
  double epsilon = 0.525;
  double rho = 0.8;
  double gamma = 0.2;
  double kappa = 0.5;
  double eta = 0.5;
  double delta = 1;
  double p[14] = {0.002*3, 0.005*3, 0.01*3, 0.04*3, 0.08*3, 0.04*3, 0.01*3, 0.005*3, 0.002*3, 0.002*3, 0.002*3, 0.001*3, 0.001*3, 0};
  double q[14] = {0};
  double h[3][5] = {{0.2, 0.5, 0.9, 0.98, 0.99}, {0,0,0,0,0}, {0.2, 0.5, 0.9, 0.98, 0.99}};
  double lambda[3][SIMULATION_TIME] = {0};
  double r = 0.48;
  
  /* Array to hold information to match to data */
  int result[SIMULATION_TIME] = {0};

  /* Give state variables, parameters, and a data-simulation array to model optimizer */
  model_optimizer(S, E, I, R, alpha, beta, epsilon, rho, gamma, kappa, eta, delta, p, q, h, lambda, r, result);
  
  return 0;
}

/* FUNCTION: model_optimizer
 * -------------------------
 * To be called in main
 * Takes in pointers to the parameters of the model, then mutates them to be values that minimize the least squares residual.
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
void model_optimizer(int S[3][SIMULATION_TIME], int E[3][14][SIMULATION_TIME], int I[3][5][SIMULATION_TIME], int R[3][SIMULATION_TIME],
		     double alpha, double beta[3], double epsilon, double rho, double gamma, double kappa, double eta, double delta, double p[14], double q[14], double h[3][5], double lambda[3][SIMULATION_TIME], double r,
		     int result[SIMULATION_TIME]) {
  /* Do litterally nothing relating to optimization atm */
  model_solver(S, E, I, R, alpha, beta, epsilon, rho, gamma, kappa, eta, delta, p, q, h, lambda, r, result);

  for(int i = 0; i < SIMULATION_TIME; i++) {
    //printf("Day %d: %d\n", i, result[i]);
  }
}

/* FUNCTION: model_solver
 * ----------------------
 * Initializes state variables S, E, I, R, and t and takes in the parameters for the model.
 * Will call on model_updater (t-1) times to get (t) many days of simulated data.
 * Returns an array of integers representing the daily counts (the thing we need to match the data).
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

void model_solver(int S[3][SIMULATION_TIME], int E[3][14][SIMULATION_TIME], int I[3][5][SIMULATION_TIME], int R[3][SIMULATION_TIME],
		  double alpha, double beta[3], double epsilon, double rho, double gamma, double kappa, double eta, double delta, double p[14], double q[14], double h[3][5], double lambda[3][SIMULATION_TIME], double r,
		  int result[SIMULATION_TIME]) {

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
 * The new count is calculated by looking at the positive flux into the managed pool.
 */
int model_updater(int S[3][SIMULATION_TIME], int E[3][14][SIMULATION_TIME], int I[3][5][SIMULATION_TIME], int R[3][SIMULATION_TIME], int t,
		  double alpha, double beta[3], double epsilon, double rho, double gamma, double kappa, double eta, double delta, double p[14], double q[14], double h[3][5], double lambda[3][SIMULATION_TIME], double r,
		  int result[SIMULATION_TIME]) {
 
  int exposed_sum_h = 0;
  for(int i = 0; i < 14; i++) {
    exposed_sum_h += E[0][i][t];
  }
  int infected_sum_h = 0;
  for(int i = 0; i < 5; i++) {
    infected_sum_h += I[0][i][t];
  }
  int exposed_sum_c = 0;
  for(int i = 0; i < 14; i++) {
    exposed_sum_c += E[2][i][t];
  }
  int infected_sum_c = 0;
  for(int i = 0; i < 5; i++) {
    infected_sum_c += I[2][i][t];
  }
  int exposed_sum_m = 0;
  for(int i = 0; i < 14; i++) {
    exposed_sum_m += E[1][i][t];
  }
  int infected_sum_m = 0;
  for(int i = 0; i < 5; i++) {
    infected_sum_m += I[1][i][t];
  }

  int Er = 1435;
  int Et = 15984;
  double beta_t = 0;

  if(t < 25) {
    beta_t = beta[0];
  }
  else if(t <= 54) {
    beta_t = beta[1];
  }
  else {
    beta_t = beta[2];
  }
    

  int N_c = S[2][t] + exposed_sum_c + infected_sum_c + R[2][t] + rho * (S[0][t] + exposed_sum_h + infected_sum_h + R[0][t]);
  int N_h = S[0][t] + exposed_sum_h + infected_sum_h + R[0][t] + exposed_sum_m + infected_sum_m;

  
  /* lambda[2][t] = (beta_t * (infected_sum_c + epsilon * exposed_sum_c) +
		  rho * beta_t * (infected_sum_h + epsilon * exposed_sum_h) +
		  gamma * beta_t * epsilon * exposed_sum_m +
		  beta_t * delta * (Er + alpha * Et)) / N_c; */
  
  lambda[2][t] = (beta_t * (infected_sum_c + epsilon * exposed_sum_c) +
		  rho * beta_t * (infected_sum_h + epsilon * exposed_sum_h) +
		  gamma * beta_t * epsilon * exposed_sum_m) / N_c;
  
  lambda[0][t] = rho * lambda[2][t] + eta * beta_t * (infected_sum_h + epsilon * exposed_sum_h + kappa * infected_sum_m) / N_h;
  
  S[2][t+1] = exp(-lambda[2][t]) * S[2][t];
  E[2][1][t+1] = (1 - exp(-lambda[2][t])) * S[2][t];
  
  for(int i = 2; i < 14; i++) {
    E[2][i][t+1] = (1 - p[i-1]) * E[2][i-1][t];
  }
  
  for(int i = 0; i < 14; i++) {
    I[2][1][t+1] += p[i]*(1-q[i])*E[2][i][t];
  }
  I[2][2][t+1] = (1-h[2][1])*I[2][1][t];
  I[2][3][t+1] = (1-h[2][2])*I[2][2][t] + (1-r)*(1-h[2][3])*I[2][3][t];
  for(int j = 4; j < 5; j++) {
    I[2][j][t] = (r) * (1-h[2][j-1])*I[2][j-1][t] + (1-r)*(1-h[2][j])*I[2][j][t];
  }
  R[2][t+1] = R[2][t] + r*I[2][5][t]+r*I[1][5][t];
  
  S[0][t+1] = exp(-lambda[0][t]) * S[0][t];
  E[0][1][t+1] = (1 - exp(-lambda[0][t])) * S[0][t];
  
  for(int i = 2; i < 14; i++) {
    E[0][i][t+1] = (1 - p[i-1]) * E[0][i-1][t];
  }
  
  for(int i = 0; i < 14; i++) {
    I[0][1][t+1] += p[i]*(1-q[i])*E[0][i][t];
  }
  I[0][2][t+1] = (1-h[0][1])*I[0][1][t];
  I[0][3][t+1] = (1-h[0][2])*I[0][2][t] + (1-r)*(1-h[0][3])*I[0][3][t];
  for(int j = 4; j <= 5; j++) {
    I[0][j][t] = r*(1-h[0][j-1])*I[0][j-1][t] + (1-r)*(1-h[0][j])*I[0][j][t];
  }
  R[0][t+1] = R[0][t] + r*I[0][5][t]+r*I[1][5][t];
  
  for(int i = 2; i < 14; i++) {
    E[1][i][t+1] = (1 - p[i - 1])*(q[i] * E[2][i - 1][t] + E[1][i - 1][t]);
  }

  I[1][1][t+1] = 0;
  for(int i = 1; i < 14; i++) {
    I[1][1][t+1] += p[i]*(q[i]*E[2][i][t] + E[1][i][t]);
  }			  
  I[1][2][t+1] = h[2][1] * I[2][1][t] + h[0][1] * I[0][1][t] + I[1][1][t];
  I[1][3][t + 1] = (h[2][2] * I[2][2][t] + h[0][2] * I[0][2][t] + I[1][2][t])
    + (1 - r)*(h[2][3] * I[2][3][t] + h[0][3] * I[0][3][t] + I[1][1][t]);
  
  for(int j = 4; j < 5; j++) {
    I[1][j][t + 1] = r * (h[2][j - 1] * I[2][j - 1][t] + h[0][j - 1] * I[0][j - 1][t] + I[1][j - 1][t])
      + (1 - r)*(h[2][j] * I[2][j][t] + h[0][j] * I[0][j][t] + I[1][j][t]);
  }

  double return_value = 0;
  return_value = h[2][0] * I[2][0][t] + h[2][1] * I[2][1][t] +
    h[2][2] * I[2][2][t] + h[2][3] * I[2][3][t] + (1 - r) * h[2][4] * I[2][4][t] +
    h[0][0] * I[0][0][t] + h[0][1] * I[0][1][t] + h[0][2] * I[0][2][t] +
    h[0][3] * I[0][3][t] + (1-r) * h[0][4] * I[0][4][t];

  /* 
  for(int i = 2; i < 14; i++) {
    return_value += (1 - p[i - 1])*(q[i] * E[2][i - 1][t]); // Will always be +0 due to q = {0}
  }

  for(int i = 1; i < 14; i++) {
    return_value += p[i]*(q[i]*E[2][i][t]); // Will always be +0 due to q = {0}
  }			  
  return_value += h[2][1] * I[2][1][t] + h[0][1] * I[0][1][t];
  return_value += (h[2][2] * I[2][2][t] + h[0][2] * I[0][2][t])
    + (1 - r)*(h[2][3] * I[2][3][t] + h[0][3] * I[0][3][t] + I[1][1][t]);
  
  for(int j = 4; j < 5; j++) {
    return_value += r * (h[2][j - 1] * I[2][j - 1][t] + h[0][j - 1] * I[0][j - 1][t])
      + (1 - r)*(h[2][j] * I[2][j][t] + h[0][j] * I[0][j][t]);
  }
  */
  printf("Day %d: %f\n", t, return_value);
  return (int) return_value;
}
