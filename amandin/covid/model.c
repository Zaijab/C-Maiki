/* FILE: model.c
 * -------------
 * To implement the model found at
 * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1691475/
 */

#include "model.h"
#include <math.h>  // Using exp()
#include <stdio.h> // Using printf()

/* FUNCTION: main
 * --------------
 * Entry point of the program.
 * Instantiates the model, runs model_updator, then writes the newly infected
 * count to a text file.
 *
 * PARAMETERS:
 * -----------
 * argc : double : the number of commandline arguments when running the program.
 * argv : char** : array of character arrays, each index holds one commandline
 * argument.
 *
 * RETURN:
 * -------
 * 0 if the function runs without error.
 * 1 if the function runs with error.
 */
int main(int argc, char *argv[]) {
  /* Initialize the model */
  struct model model_values = initial_model;

  /* Update the simulation SIMULATION_TIME many times */
  while (model_values.t < SIMULATION_TIME) {
    model_updater(&model_values);
  }

  /* Write values of simulation to file for plotting purposes */
  FILE *data_file = fopen("new_cases.txt", "w+");
  for (int i = 0; i < SIMULATION_TIME; i++) {
    fprintf(data_file, "%d %f\n", i, model_values.newly_infected[i]);
  }
  fclose(data_file);

  return 0;
}

void model_updater(struct model *m) {
  /* UPDATE t VARIABLES */
  if (m->t < 25) {
    m->beta = m->betas[0];
  } else if (m->t < 80) {
    m->beta = m->betas[1];
  } else {
    m->beta = m->betas[2];
  }

  for (int i = 0; i < INFECTION_STAGES; i++) {
    m->I_sum[HEALTHCARE_UNMANAGED][m->t] += m->I[HEALTHCARE_UNMANAGED][i][m->t];
    m->I_sum[HEALTHCARE_MANAGED][m->t] += m->I[HEALTHCARE_MANAGED][i][m->t];
    m->I_sum[COMMUNITY_UNMANAGED][m->t] += m->I[COMMUNITY_UNMANAGED][i][m->t];
    m->I_sum[COMMUNITY_MANAGED][m->t] += m->I[COMMUNITY_MANAGED][i][m->t];
  }

  for (int i = 0; i < INCUBATION_PERIOD; i++) {
    m->E_sum[HEALTHCARE_UNMANAGED][m->t] += m->E[HEALTHCARE_UNMANAGED][i][m->t];
    m->E_sum[HEALTHCARE_MANAGED][m->t] += m->E[HEALTHCARE_MANAGED][i][m->t];
    m->E_sum[COMMUNITY_UNMANAGED][m->t] += m->E[COMMUNITY_UNMANAGED][i][m->t];
    m->E_sum[COMMUNITY_MANAGED][m->t] += m->E[COMMUNITY_MANAGED][i][m->t];
  }

  m->N[COMMUNITY_UNMANAGED][m->t] =
      m->S[COMMUNITY_UNMANAGED][m->t] + m->E_sum[COMMUNITY_UNMANAGED][m->t] +
      m->I_sum[COMMUNITY_UNMANAGED][m->t] + m->R[COMMUNITY_UNMANAGED][m->t] +
      m->rho * (m->S[HEALTHCARE_UNMANAGED][m->t] +
                m->E_sum[HEALTHCARE_UNMANAGED][m->t] +
                m->I_sum[HEALTHCARE_UNMANAGED][m->t] +
                m->R[HEALTHCARE_UNMANAGED][m->t]);

  m->N[HEALTHCARE_UNMANAGED][m->t] =
      m->S[HEALTHCARE_UNMANAGED][m->t] + m->E_sum[HEALTHCARE_UNMANAGED][m->t] +
      m->I_sum[HEALTHCARE_UNMANAGED][m->t] + m->R[HEALTHCARE_UNMANAGED][m->t] +
      m->nu * (m->I_sum[COMMUNITY_MANAGED][m->t] +
               m->I_sum[HEALTHCARE_MANAGED][m->t]);

  m->lambda[COMMUNITY_UNMANAGED][m->t] =
      m->beta *
      (m->I_sum[COMMUNITY_UNMANAGED][m->t] +
       m->epsilon * m->E_sum[COMMUNITY_UNMANAGED][m->t] +
       m->rho *
           (m->I_sum[HEALTHCARE_UNMANAGED][m->t] +
            m->epsilon * m->E_sum[HEALTHCARE_UNMANAGED][m->t] +
            m->gamma * ((1 - m->nu) * m->I_sum[HEALTHCARE_UNMANAGED][m->t] +
                        m->epsilon * m->E_sum[HEALTHCARE_UNMANAGED][m->t])) +
       m->gamma * ((1 - m->nu) * m->I_sum[COMMUNITY_UNMANAGED][m->t] +
                   m->epsilon * m->E_sum[COMMUNITY_UNMANAGED][m->t])) /
      m->N[COMMUNITY_UNMANAGED][m->t];

  m->lambda[HEALTHCARE_UNMANAGED][m->t] =
      m->rho * m->lambda[COMMUNITY_UNMANAGED][m->t] +
      m->eta * m->beta *
          (m->I_sum[HEALTHCARE_UNMANAGED][m->t] +
           m->epsilon * m->E_sum[HEALTHCARE_UNMANAGED][m->t] +
           m->kappa * m->nu *
               (m->I_sum[HEALTHCARE_UNMANAGED][m->t] +
                m->I_sum[COMMUNITY_UNMANAGED][m->t])) /
          m->N[HEALTHCARE_UNMANAGED][m->t];
  /* END UPDATING t level variables */

  if (DEBUG) {
    printf("\n---------- DAY %d ----------\n", m->t);
    printf("\tNewly Managed: %f\n", m->newly_infected[m->t]);
    printf("\tlambda_c(t): %f\n", m->lambda[COMMUNITY_UNMANAGED][m->t]);
  }

  /* UPDATE STATE VARIABLES TO t+1
   * Update procedure is the same for Healthcare & Community
   * Iterate over the different groups.
   */
  for (int g = 0; g < NUM_GROUPS; g += 2) {
    m->S[g][m->t + 1] = exp(-m->lambda[g][m->t]) * m->S[g][m->t];
    m->E[g][0][m->t + 1] = (1 - exp(-m->lambda[g][m->t])) * m->S[g][m->t];

    for (int i = 1; i < INCUBATION_PERIOD; i++) {
      m->E[g][i][m->t + 1] =
          (1 - m->p[i - 1]) * (1 - m->q[i - 1]) * m->E[g][i - 1][m->t];
    }

    for (int i = 0; i < INCUBATION_PERIOD; i++) {
      m->I[g][0][m->t + 1] += m->p[i] * (1 - m->q[i]) * m->E[g][i][m->t];
    }

    m->I[g][1][m->t + 1] = (1 - m->h[g][0]) * m->I[g][0][m->t];
    m->I[g][2][m->t + 1] = (1 - m->h[g][1]) * m->I[g][1][m->t] +
                           (1 - m->r) * (1 - m->h[g][2]) * m->I[g][2][m->t];
    for (int i = 3; i < INFECTION_STAGES; i++) {
      m->I[g][i][m->t + 1] =
          m->r * (1 - m->h[g][i - 1]) * m->I[g][i - 1][m->t] +
          (1 - m->r) * (1 - m->h[g][i]) * m->I[g][i][m->t];
    }

    for (int i = 1; i < INCUBATION_PERIOD; i++) {
      m->E[g + 1][i][m->t + 1] =
          (1 - m->p[i - 1]) *
          (m->q[i - 1] * m->E[g][i - 1][m->t] + m->E[g + 1][i - 1][m->t]);
    }

    for (int i = 0; i < INCUBATION_PERIOD; i++) {
      m->I[g + 1][0][m->t + 1] +=
          m->p[i] * (m->q[i] * m->E[g][i][m->t] + m->E[g + 1][i][m->t]);
    }

    m->I[g + 1][1][m->t + 1] =
        m->h[g][0] * m->I[g][0][m->t] + m->I[g + 1][0][m->t];

    m->I[g + 1][2][m->t + 1] =
        m->h[g][1] * m->I[g][1][m->t] + m->I[g + 1][1][m->t] +
        (1 - m->r) * (m->h[g][2] * m->I[g][2][m->t] + m->I[g + 1][2][m->t]);

    for (int i = 3; i < INCUBATION_PERIOD; i++) {
      m->I[g][i][m->t + 1] =
          m->r * (m->h[g][i - 1] * m->I[g][i - 1][m->t] +
                  m->I[g + 1][i - 1][m->t]) +
          (1 - m->r) * (m->h[g][i] * m->I[g][i][m->t] + m->I[g + 1][i][m->t]);
    }

    m->R[g][m->t + 1] = m->R[g][m->t] +
                        m->r * (m->I[g][INFECTION_STAGES - 1][m->t] +
                                m->I[g + 1][INFECTION_STAGES - 1][m->t]) +
                        (1 - m->p[INCUBATION_PERIOD - 1]) *
                            m->E[g][INCUBATION_PERIOD - 1][m->t] +
                        (1 - m->p[INCUBATION_PERIOD - 1]) *
                            m->E[g][INCUBATION_PERIOD - 1][m->t];

    for (int i = 0; i < INFECTION_STAGES; i++) {
      m->newly_infected[m->t + 1] += m->h[g][i] * m->I[g + 1][i][m->t];
    }
    m->newly_infected[m->t + 1] += (1 - m->r) * m->h[g][INFECTION_STAGES - 1] *
                                   m->I[g + 1][INFECTION_STAGES - 1][m->t];
  }

  m->t += 1;
}
