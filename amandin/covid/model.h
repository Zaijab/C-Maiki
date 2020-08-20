#ifndef MODEL_H
#define MODEL_H
/* Define a number of array sizes relevent to the model
 * Useful to make the code more readable and modular
 * Necessary due to the inability to have variable length arrays
 */

/* Length of time to run the simulation */
#define SIMULATION_TIME 150

/* Partitioning the groups of the model */
#define HEALTHCARE 0
#define MANAGED 1
#define COMMUNITY 2
#define NUM_GROUPS 3

/* Stages of incubation and infection */
#define INCUBATION_PERIOD 14
#define INFECTION_STAGES 5

/* The number of Betas that will be used in the simulation */
#define NUM_BETAS 3

/* Change and recompile to have debug print statements */
#define DEBUG 0

/* Allocate memory for the models state variables and parameters */
struct model {
  double S[NUM_GROUPS][SIMULATION_TIME];
  double E[NUM_GROUPS][INCUBATION_PERIOD][SIMULATION_TIME];
  double I[NUM_GROUPS][INFECTION_STAGES][SIMULATION_TIME];
  double R[NUM_GROUPS][SIMULATION_TIME];
  double beta[NUM_BETAS];
  double epsilon;
  double rho;
  double gamma;
  double kappa;
  double eta;
  double delta;
  double p[INCUBATION_PERIOD];
  double q[INCUBATION_PERIOD];
  double h[NUM_GROUPS][INFECTION_STAGES];
  double lambda[NUM_GROUPS][SIMULATION_TIME];
  double r;
  int t;
} const struct inital_model = {S = {0}, E = {0}, I = {0}, R = {0},
  beta = {0.4549, 0.1042, 0.2054}, epsilon = 0.525, rho = 0.8, gamma = 0.2,
  kappa = 0.5, eta = 0.5, delta = 1,
  p = {0.006, 0.015, 0.12, 0.32, 0.44, 0.38, 0.3, 0.16, 0.12, 0.1, 0.06, 0.03, 0.03, 0}, q = {0},
  h = {{0.2, 0.5, 0.9, 0.98, 0.99}, {0,0,0,0,0}, {0.2, 0.5, 0.9, 0.98, 0.99}}, lambda = {0}, r = 0.48, t = 0};
#endif
