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
double main(double argc, char *argv[]) {
  /* Run the simulation */
  struct model model_values = initial_model;

  while (model_values.t < SIMULATION_TIME) {
    model_updater(&model_values);
  }

  /* Plot the simulation */

  return 0;
}

void model_updater(struct model *m) {

  if (m->t < 25) {
    m->beta = m->betas[0];
  } else if (m->t < 80) {
    m->beta = m->betas[1];
  } else {
    m->beta = m->betas[2];
  }

  sums(m);

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
    m->nu*(m->I_sum[COMMUNITY_MANAGED][m->t] + m->I_sum[HEALTHCARE_MANAGED][m->t]);

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

  if (DEBUG) {
    printf("\n---------- DAY %d ----------\n", m->t);
    printf("\tS[C][t]: %f\n", m->S[COMMUNITY_UNMANAGED][m->t]);
  }

  m->S[COMMUNITY_UNMANAGED][m->t + 1] =
      exp(-m->lambda[COMMUNITY_UNMANAGED][m->t]) * m->S[COMMUNITY_UNMANAGED][m->t];
  m->E[COMMUNITY_UNMANAGED][0][m->t + 1] =
      (1 - exp(-m->lambda[COMMUNITY_UNMANAGED][m->t])) * m->S[COMMUNITY_UNMANAGED][m->t];

  for (int i = 1; i < INCUBATION_PERIOD; i++) {
    m->E[COMMUNITY_UNMANAGED][i][m->t + 1] =
        (1 - m->p[i - 1]) * (1 - m->q[i - 1]) *
        m->E[COMMUNITY_UNMANAGED][i - 1][m->t];
  }

  for (int i = 0; i < INCUBATION_PERIOD; i++) {
    m->I[COMMUNITY_UNMANAGED][0][m->t + 1] +=
        m->p[i] * (1 - m->q[i]) * m->E[COMMUNITY_UNMANAGED][i][m->t];
  }

  m->I[COMMUNITY_UNMANAGED][1][m->t + 1] =
      (1 - m->h[COMMUNITY_UNMANAGED][0]) * m->I[COMMUNITY_UNMANAGED][0][m->t];
  m->I[COMMUNITY_UNMANAGED][2][m->t + 1] =
      (1 - m->h[COMMUNITY_UNMANAGED][1]) * m->I[COMMUNITY_UNMANAGED][1][m->t] +
      (1 - m->r) * (1 - m->h[COMMUNITY_UNMANAGED][2]) *
          m->I[COMMUNITY_UNMANAGED][2][m->t];
  for (int i = 3; i < INFECTION_STAGES; i++) {
    m->I[COMMUNITY_UNMANAGED][i][m->t + 1] =
        m->r * (1 - m->h[COMMUNITY_UNMANAGED][i - 1]) *
            m->I[COMMUNITY_UNMANAGED][i - 1][m->t] +
        (1 - m->r) * (1 - m->h[COMMUNITY_UNMANAGED][i]) *
            m->I[COMMUNITY_UNMANAGED][i][m->t];
  }

  for (int i = 1; i < INCUBATION_PERIOD; i++) {
    m->E[COMMUNITY_MANAGED][i][m->t + 1] =
        (1 - m->p[i - 1]) *
        (m->q[i - 1] * m->E[COMMUNITY_UNMANAGED][i - 1][m->t] +
         m->E[COMMUNITY_MANAGED][i - 1][m->t]);
  }

  for (int i = 0; i < INCUBATION_PERIOD; i++) {
    m->I[COMMUNITY_UNMANAGED][0][m->t + 1] +=
        m->p[i] * (m->q[i] * m->E[COMMUNITY_UNMANAGED][i][m->t] +
                   m->E[COMMUNITY_UNMANAGED][i][m->t]);
  }

  m->I[COMMUNITY_UNMANAGED][1][m->t + 1] = m->h[COMMUNITY_UNMANAGED][0]
                                           m->t += 1;
}

void sums(struct model *m) {
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
}
