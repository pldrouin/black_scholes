#ifndef _BLACK_SCHOLES_
#define _BLACK_SCHOLES_

#include <math.h>
#include <string.h>

#include "root_finder.h"

#ifndef M_SQRT2
#define M_SQRT2        1.41421356237309504880
#endif

typedef struct
{
  double IV; //Annualised implied (relative) volatility
  double OP; //Option price at time t

  double delta;
  double gamma;
  double dgammadS;
  double vega;
  double theta;
  double rho;
} bs_partial_data;

typedef struct
{
  //Common values between call and put options
  double T; //Time of option expiration
  double K; //strike price
  double t; //Evaluation time
  double St; //Asset Price at time t
  double r; //Annualised risk-free rate

  bs_partial_data call;
  bs_partial_data put;
} bs_data;

float bs_logf(float x);
float bs_erfinv(float z);

inline static double bs_N(const double x) {return 0.5*(erf(x/M_SQRT2)+1);}
inline static float bs_Ninv(const double z) {return bs_erfinv(2*z - 1)*M_SQRT2;}
inline static double bs_dNdx(const double x) {return 0.398942280401432702863*exp(-0.5*x*x);}

inline static void bs_compute_call_put_prices(bs_data* bsd, const double sigma)
{
  if(bsd->t>=bsd->T) {

    if(bsd->St>bsd->K) {
      bsd->call.OP=bsd->St-bsd->K;
      bsd->put.OP=0;

    } else {
      bsd->call.OP=0;
      bsd->put.OP=bsd->K-bsd->St;
    }
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double Kemrtau = bsd->K * exp(-bsd->r * tau);
  const double sqrttau = sqrt(tau);
  bsd->call.IV=bsd->put.IV=sigma;
  const double sigsqrttau = sigma * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * sigma * sigma) * tau) / sigsqrttau;
  bsd->call.OP = bs_N(d1) * bsd->St - bs_N(d1-sigsqrttau) * Kemrtau;
  bsd->put.OP = bsd->call.OP + Kemrtau - bsd->St;
}

inline static void bs_compute_call_put_prices_and_greeks(bs_data* bsd, const double sigma)
{
  bs_partial_data* cd=&bsd->call;
  bs_partial_data* pd=&bsd->put;

  if(bsd->t>=bsd->T) {
    memset(cd,0,sizeof(bs_partial_data));
    memset(pd,0,sizeof(bs_partial_data));

    if(bsd->St>bsd->K)
      cd->OP=bsd->St-bsd->K;

    else
      pd->OP=bsd->K-bsd->St;

    return;
  }
  const double tau = bsd->T-bsd->t;
  const double Kemrtau = bsd->K * exp(-bsd->r * tau);
  const double sqrttau = sqrt(tau);
  cd->IV=pd->IV=sigma;
  const double sigsqrttau = sigma * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * sigma * sigma) * tau) / sigsqrttau;
  const double Nd1 = bs_N(d1);
  const double Nd2 = bs_N(d1 - sigsqrttau);
  cd->OP = Nd1 * bsd->St - Nd2 * Kemrtau;
  pd->OP = cd->OP + Kemrtau - bsd->St;

  const double Npd1 = bs_dNdx(d1);
  cd->delta = Nd1;
  pd->delta = cd->delta - 1;
  cd->gamma = pd->gamma = Npd1 / (bsd->St * cd->IV * sqrttau);
  cd->vega = pd->vega = bsd->St * Npd1 * sqrttau;
  const double theta_first_term = -0.5 * bsd->St * Npd1 * cd->IV / sqrttau;
  cd->theta = theta_first_term - bsd->r * Kemrtau * Nd2;
  pd->theta = theta_first_term + bsd->r * Kemrtau * (1 - Nd2);
  cd->rho = Kemrtau * tau * Nd2;
  pd->rho = Kemrtau * tau * (Nd2-1);
}

inline static void bs_compute_call_price(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->call.OP = (bsd->St>bsd->K?bsd->St-bsd->K:0);
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double sigsqrttau = bsd->call.IV * sqrt(tau);
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->call.IV * bsd->call.IV) * tau) / sigsqrttau;
  bsd->call.OP = bs_N(d1) * bsd->St - bs_N(d1 - sigsqrttau) * bsd->K * exp(-bsd->r * tau);
}

inline static void bs_compute_call_delta(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->call.delta = 0;
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double sigsqrttau = bsd->call.IV * sqrt(tau);
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->call.IV * bsd->call.IV) * tau) / sigsqrttau;
  bsd->call.delta = bs_N(d1);
}

inline static void bs_compute_call_gamma(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->call.gamma = 0;
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double sqrttau = sqrt(tau);
  const double sigsqrttau = bsd->call.IV * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->call.IV * bsd->call.IV) * tau) / sigsqrttau;
  bsd->call.gamma = bs_dNdx(d1) / (bsd->St * bsd->call.IV * sqrttau);
}

inline static void bs_compute_call_delta_gamma(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->call.delta = 0;
    bsd->call.gamma = 0;
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double sqrttau = sqrt(tau);
  const double sigsqrttau = bsd->call.IV * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->call.IV * bsd->call.IV) * tau) / sigsqrttau;
  bsd->call.delta = bs_N(d1);
  bsd->call.gamma = bs_dNdx(d1) / (bsd->St * bsd->call.IV * sqrttau);
}

inline static void bs_compute_call_gamma_dgammadS(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->call.gamma = 0;
    bsd->call.dgammadS = 0;
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double sqrttau = sqrt(tau);
  const double sigsqrttau = bsd->call.IV * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->call.IV * bsd->call.IV) * tau) / sigsqrttau;
  bsd->call.gamma = bs_dNdx(d1) / (bsd->St * bsd->call.IV * sqrttau);
  bsd->call.dgammadS = -bsd->call.gamma/bsd->St * (1 + d1 / sigsqrttau);
}

inline static void bs_compute_call_price_and_greeks(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->call.OP = (bsd->St>bsd->K?bsd->St-bsd->K:0);
    bsd->call.delta = 0;
    bsd->call.gamma = 0;
    bsd->call.vega = 0;
    bsd->call.theta = 0;
    bsd->call.rho = 0;
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double Kemrtau = bsd->K * exp(-bsd->r * tau);
  const double sqrttau = sqrt(tau);
  const double sigsqrttau = bsd->call.IV * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->call.IV * bsd->call.IV) * tau) / sigsqrttau;
  const double Nd1 = bs_N(d1);
  const double Nd2 = bs_N(d1 - sigsqrttau);
  const double Npd1 = bs_dNdx(d1);
  bsd->call.OP = Nd1 * bsd->St - bs_N(d1 - sigsqrttau) * bsd->K * exp(-bsd->r * tau);
  bsd->call.delta = Nd1;
  bsd->call.gamma = Npd1 / (bsd->St * bsd->call.IV * sqrttau);
  bsd->call.vega = bsd->St * Npd1 * sqrttau;
  bsd->call.theta = -0.5 * bsd->St * Npd1 * bsd->call.IV / sqrttau - bsd->r * Kemrtau * Nd2;
  bsd->call.rho = Kemrtau * tau * Nd2;
}

inline static void bs_compute_put_price(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->put.OP = (bsd->St<bsd->K?bsd->K-bsd->St:0);
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double sigsqrttau = bsd->put.IV * sqrt(tau);
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->put.IV * bsd->put.IV) * tau) / sigsqrttau;
  bsd->put.OP = (bs_N(d1) - 1) * bsd->St + (1 - bs_N(d1 - sigsqrttau)) * bsd->K * exp(-bsd->r * tau);
}

inline static void bs_compute_put_delta(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->put.delta = 0;
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double sigsqrttau = bsd->put.IV * sqrt(tau);
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->put.IV * bsd->put.IV) * tau) / sigsqrttau;
  bsd->put.delta =  bs_N(d1) - 1;
}

inline static void bs_compute_put_gamma(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->put.gamma = 0;
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double sqrttau = sqrt(tau);
  const double sigsqrttau = bsd->put.IV * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->put.IV * bsd->put.IV) * tau) / sigsqrttau;
  bsd->put.gamma = bs_dNdx(d1) / (bsd->St * bsd->put.IV * sqrttau);
}

inline static void bs_compute_put_delta_gamma(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->put.delta = 0;
    bsd->put.gamma = 0;
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double sqrttau = sqrt(tau);
  const double sigsqrttau = bsd->put.IV * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->put.IV * bsd->put.IV) * tau) / sigsqrttau;
  bsd->put.delta = bs_N(d1) - 1;
  bsd->put.gamma = bs_dNdx(d1) / (bsd->St * bsd->put.IV * sqrttau);
}

inline static void bs_compute_put_gamma_dgammadS(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->put.gamma = 0;
    bsd->put.dgammadS = 0;
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double sqrttau = sqrt(tau);
  const double sigsqrttau = bsd->put.IV * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->put.IV * bsd->put.IV) * tau) / sigsqrttau;
  bsd->put.gamma = bs_dNdx(d1) / (bsd->St * bsd->put.IV * sqrttau);
  bsd->put.dgammadS = -bsd->put.gamma/bsd->St * (1 + d1 / sigsqrttau);
}

inline static void bs_compute_put_price_and_greeks(bs_data* bsd)
{
  if(bsd->t>=bsd->T) {
    bsd->put.OP = (bsd->St<bsd->K?bsd->K-bsd->St:0);
    bsd->put.delta = 0;
    bsd->put.gamma = 0;
    bsd->put.vega = 0;
    bsd->put.theta = 0;
    bsd->put.rho = 0;
    return;
  }
  const double tau = bsd->T-bsd->t;
  const double Kemrtau = bsd->K * exp(-bsd->r * tau);
  const double sqrttau = sqrt(tau);
  const double sigsqrttau = bsd->put.IV * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * bsd->put.IV * bsd->put.IV) * tau) / sigsqrttau;
  const double Nd1 = bs_N(d1);
  const double Nd2 = bs_N(d1 - sigsqrttau);
  const double Npd1 = bs_dNdx(d1);
  bsd->put.OP = (Nd1 - 1) * bsd->St + (1 - bs_N(d1 - sigsqrttau)) * bsd->K * exp(-bsd->r * tau);
  bsd->put.delta = Nd1 - 1;
  bsd->put.gamma = Npd1 / (bsd->St * bsd->put.IV * sqrttau);
  bsd->put.vega = bsd->St * Npd1 * sqrttau;
  bsd->put.theta = -0.5 * bsd->St * Npd1 * bsd->put.IV / sqrttau + bsd->r * Kemrtau * (1 - Nd2);
  bsd->put.rho = Kemrtau * tau * (Nd2-1);
}

inline static void bs_compute_strike_from_call_delta(bs_data* bsd)
{
  const double tau = bsd->T-bsd->t;
  bsd->K = bsd->St * exp((bsd->r + 0.5 * bsd->call.IV*bsd->call.IV ) * tau - bsd->call.IV * sqrt(tau) * bs_Ninv(bsd->call.delta));
  //printf("St=%f, tau=%f, r=%f, IV=%f, delta=%f\n, K=%f\n",bsd->St,tau,bsd->r,bsd->call.IV,bsd->call.delta,bsd->K);
}

inline static void bs_compute_strike_from_put_delta(bs_data* bsd)
{
  const double tau = bsd->T-bsd->t;
  bsd->K = bsd->St * exp((bsd->r + 0.5 * bsd->put.IV*bsd->put.IV ) * tau - bsd->put.IV * sqrt(tau) * bs_Ninv(bsd->put.delta + 1));
  //printf("St=%f, tau=%f, r=%f, IV=%f, delta=%f\n, K=%f\n",bsd->St,tau,bsd->r,bsd->put.IV,bsd->put.delta,bsd->K);
}

/**
 * @brief Computes the probability of the call option to be ITM based on its
 * delta value, IV and T-t
 */
inline static double bs_compute_call_PITM_from_call_delta(bs_data* bsd)
{
  return bs_N(bs_Ninv(bsd->call.delta) - bsd->call.IV * sqrt(bsd->T-bsd->t));
}

/**
 * @brief Computes the probability of the put option to be ITM based on its
 * delta value, IV and T-t
 */
inline static double bs_compute_put_PITM_from_put_delta(bs_data* bsd)
{
  return 1 - bs_N(bs_Ninv(bsd->put.delta + 1) - bsd->put.IV * sqrt(bsd->T-bsd->t));
}

/**
 * @brief Computes the delta of a call option that is necessary to obtained the
 * requested probability for this option to be ITM, based on its IV and T-t
 * values
 */
inline static void bs_compute_call_delta_from_call_PITM(bs_data* bsd, const double pitm)
{
  bsd->call.delta =  bs_N(bs_Ninv(pitm) + bsd->call.IV * sqrt(bsd->T-bsd->t));
}

/**
 * @brief Computes the delta of a put option that is necessary to obtained the
 * requested probability for this option to be ITM, based on its IV and T-t
 * values
 */
inline static void bs_compute_put_delta_from_put_PITM(bs_data* bsd, const double pitm)
{
  bsd->put.delta = bs_N(bs_Ninv(1 - pitm) + bsd->put.IV * sqrt(bsd->T-bsd->t)) - 1;
}

/**
 * @brief IV value estimate using the Newton's method.
 *
 * @param iv: Current iv value (input/output).
 * @param diff: Current value for the option price discrepancy. (output)
 * @param params: Pointer to bsd handle. (input)
 */
inline static double _bs_compute_call_iv_estimate(const double iv, double* diff, void* params)
{
  bs_data* bsd=(bs_data*)params;
  //printf("IV: %f -> ",*iv);

  const double tau = bsd->T - bsd->t;
  const double sqrttau = sqrt(tau);
  const double sigsqrttau = iv * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * iv * iv) * tau) / sigsqrttau;
  const double Kemrtau = bsd->K * exp(-bsd->r * tau);
  const double op = bs_N(d1) * bsd->St - bs_N(d1 - sigsqrttau) * Kemrtau;
  *diff = op - bsd->call.OP;
  const double d1oiv = d1/iv;
  const double dCdIV = bs_dNdx(d1) * (-d1oiv + sqrttau) * bsd->St + bs_dNdx(d1 - sigsqrttau) * d1oiv * Kemrtau; 
  //printf("%f, price %f, diff %f\n",iv-*diff / dCdIV,op,*diff);
  return -*diff / dCdIV;
}

/**
 * @brief IV value estimate using the Newton's method.
 *
 * @param iv: Current iv value (input/output).
 * @param diff: Current value for the option price discrepancy. (output)
 * @param params: Pointer to bsd handle. (input)
 */
inline static double _bs_compute_put_iv_estimate(const double iv, double* diff, void* params)
{
  bs_data* bsd=(bs_data*)params;
  //printf("IV: %f -> ",iv);

  const double tau = bsd->T - bsd->t;
  const double sqrttau = sqrt(tau);
  const double sigsqrttau = iv * sqrttau;
  const double d1 = (log(bsd->St / bsd->K) + (bsd->r + 0.5 * iv * iv) * tau) / sigsqrttau;
  const double Kemrtau = bsd->K * exp(-bsd->r * tau);
  const double op = (bs_N(d1) - 1) * bsd->St + (1 - bs_N(d1 - sigsqrttau)) * Kemrtau;
  *diff = op - bsd->put.OP;
  const double d1oiv = d1/iv;
  const double dPdIV = bs_dNdx(d1) * (-d1oiv + sqrttau) * bsd->St + bs_dNdx(d1 - sigsqrttau) * d1oiv * Kemrtau; 
  //printf("%f, price %f, diff %f\n",iv-*diff / dPdIV,op,*diff);
  return -*diff / dPdIV;
}

int bs_compute_call_iv(bs_data* bsd);
int bs_compute_put_iv(bs_data* bsd);

#endif
