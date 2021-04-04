#ifndef MODEL_H
#define MODEL_H
/* Define a number of array sizes relevent to the model
 * Useful to make the code more readable and modular
 * Necessary due to the inability to have variable length arrays
 */

/* Length of time to run the simulation */
#define SIMULATION_TIME 20

/* Partitioning the groups of the model */
#define HEALTHCARE_UNMANAGED 0
#define HEALTHCARE_MANAGED 1
#define COMMUNITY_UNMANAGED 2
#define COMMUNITY_MANAGED 3
#define NUM_GROUPS 4

/* Stages of incubation and infection */
#define INCUBATION_PERIOD 14
#define INFECTION_STAGES 5

/* The number of Betas that will be used in the simulation */
#define NUM_BETAS 3

/* Change and recompile to have debug print statements */
#define DEBUG 1

/* Allocate memory for the models state variables and parameters */
struct model {
  double S[NUM_GROUPS][SIMULATION_TIME];
  double E[NUM_GROUPS][INCUBATION_PERIOD][SIMULATION_TIME];
  double I[NUM_GROUPS][INFECTION_STAGES][SIMULATION_TIME];
  double R[NUM_GROUPS][SIMULATION_TIME];
  double E_sum[NUM_GROUPS][SIMULATION_TIME];
  double I_sum[NUM_GROUPS][SIMULATION_TIME];
  double betas[NUM_BETAS];
  double beta;
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
  double nu;
  double newly_infected[SIMULATION_TIME];
  double N[NUM_GROUPS][SIMULATION_TIME];
  int t;
} const initial_model = {
    .S = {[HEALTHCARE_UNMANAGED] = {[0] = 15000}, [COMMUNITY_UNMANAGED] = {[0] = 937711}},
    .E = {0},
    .I = {[COMMUNITY_UNMANAGED] = {[0] = {[0] = 1}}},
    .R = {0},
    .E_sum = {0},
    .I_sum = {0},
    .betas = {0.4549, 0.1042, 0.2054},
    .beta = 0,
    .epsilon = 0.525,
    .rho = 0.8,
    .gamma = 0.2,
    .kappa = 0.5,
    .eta = 0.5,
    .delta = 1,
    .p = {0.006, 0.015, 0.12, 0.32, 0.44, 0.38, 0.3, 0.16, 0.12, 0.1, 0.06,
          0.03, 0.03, 0},
    .q = {0},
    .h = {{0.2, 0.5, 0.9, 0.98, 0.99},
          {0, 0, 0, 0, 0},
          {0.2, 0.5, 0.9, 0.98, 0.99}},
    .lambda = {0},
    .r = 0.48,
    .nu = 0.08,
    .newly_infected = {0},
    .N = {0},
    .t = 0};

/* Implement the difference equations in this function */
void model_updater(struct model *m);
#endif
