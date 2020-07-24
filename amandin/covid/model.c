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
 * c -> >=2
 * c (the community group) will be 2 for the time being, but the plans are for the model to expand.
 * So there can be more community indicies later.
 */

#include <math.h>
#include <stdio.h>

int model_updater(double** S, double*** E, double*** I, double** R, int t,
		  double** l, double* p, double* q, double** h, double r);

void model_solver(int t, double** l, double* p, double* q, double** h, double r, int* result);

void model_optimizer(int t, double*** l, double** p, double** q, double*** h, double* r);
 
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
int main(int argc, char** argv){
  int t = 100;
  double l[3][t]; // [Group, Last day of simuation]
  double p[14]; // [Day of exposure]
  double q[14]; // [Day of exposure]
  double h[3][5]; // [Group, Day of infection]
  double r; // Chance to recover

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < t; j++) {
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
    for(int j = 0; j < 5; j++) {
      h[i][j] = 0;
    }
  }

  r = 0;
  
  model_optimizer(t, &l, &p, &q, &h, &r);
  
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
int model_updater(double** S, double*** E, double*** I, double** R, int t,
		   double** l, double* p, double* q, double** h, double r) {
  S[2][t+1] = exp(-l[2][t]) * S[2][t];
  E[2][1][t+1] = (1 - exp(-l[2][t])) * S[2][t];
  
  for(int i = 2; i <= 14; i++) {
    E[2][i][t+1] = (1 - p[i-1]) * E[2][i-1][t];
  }
  
  for(int i = 1; i <= 14; i++) {
    I[2][1][t+1] += p[i]*(1-q[i])*E[2][i][t];
  }
  I[2][2][t+1] = (1-h[2][1])*I[2][1][t];
  I[2][3][t+1] = (1-h[2][2])*I[2][2][t] + (1-r)*(1-h[2][3])*I[2][3][t];
  for(int j = 4; j <= 5; j++) {
    I[2][j][t] = (r) * (1-h[2][j-1])*I[2][j-1][t] + (1-r)*(1-h[2][j])*I[2][j][t];
  }
  R[2][t+1] = R[2][t] + r*I[2][5][t]+r*I[1][5][t];
  
  S[0][t+1] = exp(-l[0][t]) * S[0][t];
  E[0][1][t+1] = (1 - exp(-l[0][t])) * S[0][t];
  
  for(int i = 2; i <= 14; i++) {
    E[0][i][t+1] = (1 - p[i-1]) * E[0][i-1][t];
  }
  
  for(int i = 1; i <= 14; i++) {
    I[0][1][t+1] += p[i]*(1-q[i])*E[0][i][t];
  }
  I[0][2][t+1] = (1-h[0][1])*I[0][1][t];
  I[0][3][t+1] = (1-h[0][2])*I[0][2][t] + (1-r)*(1-h[0][3])*I[0][3][t];
  for(int j = 4; j <= 5; j++) {
    I[0][j][t] = r*(1-h[0][j-1])*I[0][j-1][t] + (1-r)*(1-h[0][j])*I[0][j][t];
  }
  R[0][t+1] = R[0][t] + r*I[0][5][t]+r*I[1][5][t];
  
  for(int i = 2; i <= 14; i++) {
    E[1][i][t+1] = (1 - p[i - 1])*(q[i] * E[2][i - 1][t] + E[1][i - 1][t]);
  }

  I[1][1][t+1] = 0;
  for(int i = 1; i <= 14; i++) {
    I[1][1][t+1] += p[i]*(q[i]*E[2][i][t] + E[1][i][t]);
  }			  
  I[1][2][t+1] = h[2][1] * I[2][1][t] + h[0][1] * I[0][1][t] + I[1][1][t];
  I[1][3][t + 1] = r * (h[2][2] * I[2][2][t] + h[0][2] * I[0][2][t] + I[1][2][t])
    + (1 - r)*(h[2][3] * I[2][3][t] + h[0][3] * I[0][3][t] + I[1][1][t]);
  
  for(int j = 4; j <= 5; j++) {
    I[1][j][t + 1] = r * (h[2][j - 1] * I[2][j - 1][t] + h[0][j - 1] * I[0][j - 1][t] + I[1][j - 1][t])
      + (1 - r)*(h[2][j] * I[2][j][t] + h[0][j] * I[0][j][t] + I[1][j][t]);
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
 */
void model_solver(int t, double*** l, double* p, double* q, double** h, double r, int* result) {
  double S[3][t];
  double E[3][14][t];
  double I[3][5][t];
  double R[3][t];

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < t; j++) {
	S[i][j] = 0;
    }
  }
  
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 14; j++) {
      for(int k = 0; k < t; k++) {
	E[i][j][k] = 0;
      }
    }
  }
  
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 5; j++) {
      for(int k = 0; k < t; k++) {
	I[i][j][k] = 0;
      }
    }
  }
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < t; j++) {
	R[i][j] = 0;
    }
  }

  for(int i = 0; i < t; i++) {
    result[i] = model_updater(S, E, I, R, i, l, p, q, h, r) // Using the state variables and parameters, update the model each day
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
void model_optimizer(int t, double**** l, double** p, double** q, double*** h, double* r) {

  int result[t];

  for(int i = 0; i < t; i++) {
    result[i] = 0;
  }

  model_updater(t, l*, p*, q*, h*, r*, result);

  for(i = 0; i < t; i++) {
    printf(result[i]);
  }
}
