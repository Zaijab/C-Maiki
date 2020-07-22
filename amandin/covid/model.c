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

void model_updater(float S[][], float E[][][], float I[][][], float R[][], float l[][], float p[], float q[], float h[][], float r, int t);

int* model_solver(float l[][], float p[], float q[], float h[][], float r, int t);

void model_optimizer(float* l[][], float* p[], float* q[], float* h[][], float* r, int t);

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
 */
void model_updater(float S[][], float E[][][], float I[][][], float R[][][], int t,
		   float l[][], float p[], float q[], float h[][], float r) {
  S[2][t+1] = exp(-l[2][t]) * S[2][t];
  E[2][1][t+1] = (1 - exp(-l[2][t])) * S[2][t];
  
  for(int i = 2; i <= 14; i++) {
    E[2][i][t+1] = (1 - p[i-1]) * E[2][i-1][t];
  }
  
  for(int i = 1; i <= 14; i++) {
    I[2][1][t+1] += p[i]*(1-q[i])*E[2][i][t];
  }
  I[2][2][t+1] = (1-h[2][1])*I[2][1][t];
  I[2][3][t+1] = (1-h[2][2])*I[2][2][t] + (1-r)(1-h[2][3])*I[2][3][t];
  for(int j = 4; j <= 5; j++) {
    I[2][j] = r*(1-h[2][j-1])*I[2][j-1][t] + (1-r)(1-h[2][j])*I[2][j][t];
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
  I[0][3][t+1] = (1-h[0][2])*I[0][2][t] + (1-r)(1-h[0][3])*I[0][3][t];
  for(int j = 4; j <= 5; j++) {
    I[0][j] = r*(1-h[0][j-1])*I[0][j-1][t] + (1-r)(1-h[0][j])*I[0][j][t];
  }
  R[0][t+1] = R[0][t] + r*I[0][5][t]+r*I[1][5][t];
  
  for(int i = 2; i <= 14; i++) {
    E[1][t+1] = (1 - p[i - 1])*(q[i] * E[2][i - 1][t] + E[1][i - 1][t]);
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
int* model_solver(int t, float l[][][], float p[], float q[], float h[][], float r) {
  float S[3][t];
  float E[3][14][t];
  float I[3][5][t];
  float R[3][t];
  
  return {1};
}

/* FUNCTION: model_optimizer
 * -------------------------
 * To be called in main
 * Takes in pointers to the parameters of the model, then mutates them to be values that minimize the least squares residual.
 * 
 * PARAMETERS:
 * -----------
 */
void model_optimizer(float* R[][][], float* l[][][], float* p[], float* q[], float* h[][], float* r, int t) {

}
