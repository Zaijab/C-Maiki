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
 */

#include <math.h>

void model_updater(float S[][], float E[][][], float I[][][], float R[][], float l[][], float p[], float q[], float h[][], float r, int t);

int* model_solver(float l[][], float p[], float q[], float h[][], float r, int t);

void model_optimizer(float* l[][], float* p[], float* q[], float* h[][], float* r, int t);

/* FUNCTION: main
 * --------------
 * Entry point of the program
 * To intstantiate initial values of the parameters then run model_optimizer to mutate them
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
 * 
 *
 * PARAMETERS:
 * -----------
 */
void model_updater(float S[][], float E[][][], float I[][][], float R[][][], float l[][], float p[], float q[], float h[][], float r, int t) {
  S[h][t+1] = exp(-l[h][t]) * S[h][t];
  E[h][1][t+1] = (1 - exp(-l[h][t])) * S[h][t];

  for(int i = 2; i <= 10; i++) {
    E[h][i][t+1] = (1 - p[i-1]) * E[h][i-1][t];
  }
  
  for(int i = 1; i <= 10; i++) {
    I[h][1][t+1] += p[i]*(1-q[i])*E[h][i][t]
  }
  I[h][2][t+1] = (1-h[h][1])*I[h][1][t];
  I[h][3][t+1] = (1-h[h][2])*I[h][2][t] + (1-r)(1-h[h][3])*I[h][3][t];
  for(int j = 4; j <= 5; j++) {
    I[h][j] = r*(1-h[h][j-1])*I[h][j-1][t] + (1-r)(1-h[h][j])*I[h][j][t];
  }
  R[h][t+1] = R[h][t] + r*I[h][5][t]+r*I[m][5][t]
}

/* FUNCTION: model_solver
 * ----------------------
 *
 *
 * PARAMETERS:
 * -----------
 */
int* model_solver(float l[][][], float p[], float q[], float h[][], float r, int t) {
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
