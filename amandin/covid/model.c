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

void model_optimizer(double l[3][SIMULATION_TIME], double p[14], double q[14], double h[3][5], double r);

void model_solver(double l[3][SIMULATION_TIME], double p[14], double q[14], double h[3][5], double r, int result[SIMULATION_TIME]);

int model_updater(double S[3][SIMULATION_TIME], double E[3][14][SIMULATION_TIME], double I[3][5][SIMULATION_TIME], double R[3][SIMULATION_TIME], int t,
		  double l[3][SIMULATION_TIME], double p[14], double q[14], double h[3][5], double r);
 
/* FUNCTION: main
 * --------------
 * Entry point of the program.
 * To intstantiate initial values of the parameters then run model_optimizer to mutate them. 
 *
 * PARAMETERS:
 * -----------
 * argc : int : the number of commandline arguments when running the program
 * argv : char** : array of character arrays, each index holds one commandline argument
 */
int main(int argc, char** argv) {
  double l[3][SIMULATION_TIME]; // [Group, Last day of simuation]
  double p[14]; // [Day of exposure]
  double q[14]; // [Day of exposure]
  double h[3][5]; // [Group, Day of infection]
  double r; // Chance to recover
  
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < SIMULATION_TIME; j++) {
      l[i][j] = 0;
    }
  }
  
  for(int i = 0; i < 14; i++) {
    p[i] = 0;
  }
  
  for(int i = 0; i < 14; i++) {
    q[i] = 0;
  }
  
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < SIMULATION_TIME; j++) {
      l[i][j] = 0;
    }
  }
  
  model_optimizer(l, p, q, h, r);
  
  return 0;
}

/* FUNCTION: model_updater
 * -----------------------
 * Implementation of the difference equations. 
 * The function takes in state variables S, E, I, R, and t, along with parameters, l, p, q, h, and r.
 * Using the given state variables and parameters, the model can advance one day.
 *
 * PARAMETERS:
 * -----------
 * S 
 * E
 * I
 * R
 * t
 * l
 * p
 * q
 * h
 * r
 *
 * RETURNS:
 * --------
 * The new infected count on that day.
 * The new count is calculated by looking at the positive flux into the managed pool
 */
int model_updater(double S[3][SIMULATION_TIME], double E[3][14][SIMULATION_TIME], double I[3][5][SIMULATION_TIME], double R[3][SIMULATION_TIME], int t,
		  double l[3][SIMULATION_TIME], double p[14], double q[14], double h[3][5], double r) {

  S[2][t+1] = exp(-l[2][t]) * S[2][t];
  E[2][1][t+1] = (1 - exp(-l[2][t])) * S[2][t];
  
  for(int i = 2; i < 14; i++) {
    E[2][i][t+1] = (1 - p[i-1]) * E[2][i-1][t];
  }
  
  for(int i = 1; i < 14; i++) {
    I[2][1][t+1] += p[i]*(1-q[i])*E[2][i][t];
  }
  I[2][2][t+1] = (1-h[2][1])*I[2][1][t];
  I[2][3][t+1] = (1-h[2][2])*I[2][2][t] + (1-r)*(1-h[2][3])*I[2][3][t];
  for(int j = 4; j < 5; j++) {
    I[2][j][t] = (r) * (1-h[2][j-1])*I[2][j-1][t] + (1-r)*(1-h[2][j])*I[2][j][t];
  }
  R[2][t+1] = R[2][t] + r*I[2][5][t]+r*I[1][5][t];
  
  S[0][t+1] = exp(-l[0][t]) * S[0][t];
  E[0][1][t+1] = (1 - exp(-l[0][t])) * S[0][t];
  
  for(int i = 2; i < 14; i++) {
    E[0][i][t+1] = (1 - p[i-1]) * E[0][i-1][t];
  }
  
  for(int i = 1; i < 14; i++) {
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

  int return_value = 0;
  for(int i = 2; i < 14; i++) {
    return_value += (1 - p[i - 1])*(q[i] * E[2][i - 1][t]);
  }

  for(int i = 1; i < 14; i++) {
    return_value += p[i]*(q[i]*E[2][i][t]);
  }			  
  return_value += h[2][1] * I[2][1][t] + h[0][1] * I[0][1][t];
  return_value += (h[2][2] * I[2][2][t] + h[0][2] * I[0][2][t])
    + (1 - r)*(h[2][3] * I[2][3][t] + h[0][3] * I[0][3][t] + I[1][1][t]);
  
  for(int j = 4; j < 5; j++) {
    return_value += r * (h[2][j - 1] * I[2][j - 1][t] + h[0][j - 1] * I[0][j - 1][t])
      + (1 - r)*(h[2][j] * I[2][j][t] + h[0][j] * I[0][j][t]);
  }

  return return_value;
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
 */
void model_solver(double l[3][SIMULATION_TIME], double p[14], double q[14], double h[3][5], double r, int result[SIMULATION_TIME]) {
  double S[3][SIMULATION_TIME];
  double E[3][14][SIMULATION_TIME];
  double I[3][5][SIMULATION_TIME];
  double R[3][SIMULATION_TIME];

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < SIMULATION_TIME; j++) {
      S[i][j] = 0;
    }
  }
  
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 14; j++) {
      for(int k = 0; k < SIMULATION_TIME; k++) {
	E[i][j][k] = 0;
      }
    }
  }
  
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 5; j++) {
      for(int k = 0; k < SIMULATION_TIME; k++) {
	I[i][j][k] = 0;
      }
    }
  }
  
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < SIMULATION_TIME; j++) {
      R[i][j] = 0;
    }
  }
  
  for(int i = 0; i < SIMULATION_TIME; i++) {
    result[i] = model_updater(S, E, I, R, i, l, p, q, h, r); // Using the state variables and parameters, update the model each day
  }
}

/* FUNCTION: model_optimizer
 * -------------------------
 * To be called in main
 * Takes in pointers to the parameters of the model, then mutates them to be values that minimize the least squares residual.
 * 
 *
 * PARAMETERS:
 * -----------
 */
void model_optimizer(double l[3][SIMULATION_TIME], double p[14], double q[14], double h[3][5], double r) {
  
  int result[SIMULATION_TIME];
  
  for(int i = 0; i < SIMULATION_TIME; i++) {
    result[i] = 0;
  }
  
  model_solver(l, p, q, h, r, result);

  for(int i = 0; i < SIMULATION_TIME; i++) {
    printf("Day %d: %d\n", i, result[i]);
  }
}
