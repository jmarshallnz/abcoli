// exp_samp.cpp : Defines the exported functions for the DLL application.
//

#include <cstdlib>
#include <math.h>
#include <stdio.h>

#include <R.h>

// TODO: Better rand() function would be nice - this one is icky to say the least
//       eg we could call the R functions directly here
inline double runif(double high)
{
  return (high * rand()) / RAND_MAX;
}

inline double rexp(double rate)
{
  double u = ((double)rand()) / RAND_MAX;
  while (!u)
    u = ((double)rand()) / RAND_MAX;
  return -log(u) / rate;
}

inline int rdiscrete(const unsigned int n, const double *const pdf)
{
  double scale = 0;
  for (unsigned int j = 0; j < n; j++)
    scale += pdf[j];

  const double x = unif_rand() * scale;
  double cdf = 0;
  for (unsigned int j = 0; j < n; j++)
  {
    cdf += pdf[j];
    if (x < cdf)
      return j;
  }
  return n-1;
}

extern "C"
{
  /*! \brief Get the next exponential event, given the appropriate parameters
   \param num_groups Number of groups in the SI model
   \param S susceptible animals [1..num_groups]
   \param I infected animals [1..num_groups]
   \param E local environments [..num_groups]
   \param G global environment
   \param params parameters ([..num_groups*7]
   \param event [out] returned event details [1..5]
   */
  //__declspec(dllexport) 
  void get_exp_event(int *num_groups, int *num_sub_groups, double *S, double *I, double *E, double *G, double *params, double *event)
  {
    const unsigned int n = (unsigned int)*num_groups-1;
    const unsigned int m = (unsigned int)*num_sub_groups;
    double *rates = new double[11 * n];
    double *r = rates;
    double *s_j = S;
    double *i_j = I;
    for (unsigned int j = 0; j < n; j++)
    {
      double s = 0;
      double i = 0;
      double *p = params + j*7;
      for (unsigned int k = 0; k < m; k++)
      {
        s += *s_j++;
        i += *i_j++;
      }
      *r++ = p[0] * s;        // death of susceptible
      *r++ = p[0] * i;        // death of infected
      *r++ = p[1] * s;        // culling of susceptible
      *r++ = p[1] * i;        // culling of infected
      *r++ = p[2] * s * i;    // direct infection
      *r++ = p[3] * s * E[j]; // indirect infection via consumption from local env
      *r++ = p[3] * i * E[j]; // consumption from local env by infected
      *r++ = p[4] * s * (*G); // indirect infection via consumption from global env
      *r++ = p[4] * i * (*G); // consumption from global env by infected
      *r++ = p[5] * i;        // recovery by infected
      *r++ = p[6] * (i+s);    // birthrate throughout year (non-seasonal)
    }
    double total_rate = 0;
    for (unsigned int j = 0; j < 11*n; j++)
      total_rate += rates[j];

    // generate a random exponential for time to next event, followed by a uniform for event type
    GetRNGstate();
    event[0] = exp_rand() / total_rate;
    int value = rdiscrete(11*n, rates);
    int group = value / 11;
    event[1] = (value % 11) + 1;              // type of event
    event[2] = group + 1;                     // group that it affects
    event[3] = rdiscrete(m, S + group*m) + 1; // subgroup of S
    event[4] = rdiscrete(m, I + group*m) + 1; // subgroup of I
    PutRNGstate();

    delete[] rates;
  }
}


/*
rdiscrete <- function(values)
{
  u <- runif(1, 0, sum(values))

  cumm <- 0;
  n <-length(values);
  for(j in 1:n)
  {
    cumm <- cumm + values[j];
    if (u < cumm)
      return(j);
  }
  return(length(values))
}

get_next_event_R <- function(S, I, E, G, params)
{
  num_groups = length(S[,1]);

  s <- rowSums(S);
  i <- rowSums(I);

  rates <- matrix(0, num_groups, 10)
  rates[, 1] <- params$death      * s;     # death of susceptible
  rates[, 2] <- params$death      * i;     # death of infected
  rates[, 3] <- params$culling    * s;     # culling of susceptible
  rates[, 4] <- params$culling    * i;     # culling of infected
  rates[, 5] <- params$inf_direct * s * i; # direct infection
  rates[, 6] <- params$inf_local  * s * E; # indirect infection via consumption from local env
  rates[, 7] <- params$inf_local  * i * E; # consumption from local env by infected
  rates[, 8] <- params$inf_global * s * G; # indirect infection via consumption from global env
  rates[, 9] <- params$inf_global * i * G; # consumption from global env by infected
  rates[,10] <- params$recovery   * i;     # recovery by infected

  total_rate = sum(rates);
#  cat("tr=", total_rate, "rates=", rates, "\n")

  # generate a random exponential for time to next event, followed by a uniform for event type

  event <- rexp(1, total_rate);
  event_num <- rdiscrete(rates)
  event[2] <- ((event_num-1) %/% num_groups) + 1;
  event[3] <- ((event_num-1) %% num_groups) + 1;
  event[4] <- rdiscrete(S[event[3],])
  event[5] <- rdiscrete(I[event[3],])
  return(event)
}
*/
