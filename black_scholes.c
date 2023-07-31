#include "black_scholes.h"

/*
https://stackoverflow.com/questions/27229371/inverse-error-function-in-c .
M. Giles, "Approximating the erfinv function." In GPU Computing Gems Jade Edition, pp. 109-116. 2011 (https://people.maths.ox.ac.uk/gilesm/codes/erfinv/gems.pdf)

The solution below generally follows the approach by Giles, but simplifies it in not requiring the square root for the tail portion, i.e. it uses two approximations of the type x · P(w). The code takes maximum advantage of the fused multiply-add operation FMA, which is exposed via the standard math functions fma() and fmaf() in C. Many common compute platforms, such as IBM Power, Arm64, x86-64, and GPUs offer this operation in hardware. Where no hardware support exists, the use of fma{f}() will likely make the code below unacceptably slow as the operation needs to be emulated by the standard math library. Also, functionally incorrect emulations of FMA are known to exist.

The accuracy of the standard math library's logarithm function logf() will have some impact on the accuracy of my_erfinvf() below. As long as the library provides a faithfully-rounded implementation with error < 1 ulp, the stated error bound should hold and it did for the few libraries I tried. For improved reproducability, I have included my own portable faithfully-rounded implementation, bs_logf().
*/

float bs_logf(float x)
{
  float i, m, r, s, t;
  int e;

  m = frexpf (x, &e);
  if (m < 0.666666667f) { // 0x1.555556p-1
    m = m + m;
    e = e - 1;
  }
  i = (float)e;
  /* m in [2/3, 4/3] */
  m = m - 1.0f;
  s = m * m;
  /* Compute log1p(m) for m in [-1/3, 1/3] */
  r =             -0.130310059f;  // -0x1.0ae000p-3
  t =              0.140869141f;  //  0x1.208000p-3
  r = fmaf (r, s, -0.121484190f); // -0x1.f19968p-4
  t = fmaf (t, s,  0.139814854f); //  0x1.1e5740p-3
  r = fmaf (r, s, -0.166846052f); // -0x1.55b362p-3
  t = fmaf (t, s,  0.200120345f); //  0x1.99d8b2p-3
  r = fmaf (r, s, -0.249996200f); // -0x1.fffe02p-3
  r = fmaf (t, m, r);
  r = fmaf (r, m,  0.333331972f); //  0x1.5554fap-2
  r = fmaf (r, m, -0.500000000f); // -0x1.000000p-1
  r = fmaf (r, s, m);
  r = fmaf (i,  0.693147182f, r); //  0x1.62e430p-1 // log(2)
  if (!((x > 0.0f) && (x <= 3.40282346e+38f))) { // 0x1.fffffep+127
    r = x + x;  // silence NaNs if necessary
    if (x  < 0.0f) r = ( 0.0f / 0.0f); //  NaN
    if (x == 0.0f) r = (-1.0f / 0.0f); // -Inf
  }
  return r;
}

float bs_erfinv(float z)
{
  float p, r, t;
  t = fmaf (z, 0.0f - z, 1.0f);
  t = bs_logf (t);
  if (fabsf(t) > 6.125f) { // maximum ulp error = 2.35793
    p =              3.03697567e-10f; //  0x1.4deb44p-32
    p = fmaf (p, t,  2.93243101e-8f); //  0x1.f7c9aep-26
    p = fmaf (p, t,  1.22150334e-6f); //  0x1.47e512p-20
    p = fmaf (p, t,  2.84108955e-5f); //  0x1.dca7dep-16
    p = fmaf (p, t,  3.93552968e-4f); //  0x1.9cab92p-12
    p = fmaf (p, t,  3.02698812e-3f); //  0x1.8cc0dep-9
    p = fmaf (p, t,  4.83185798e-3f); //  0x1.3ca920p-8
    p = fmaf (p, t, -2.64646143e-1f); // -0x1.0eff66p-2
    p = fmaf (p, t,  8.40016484e-1f); //  0x1.ae16a4p-1
  } else { // maximum ulp error = 2.35002
    p =              5.43877832e-9f;  //  0x1.75c000p-28
    p = fmaf (p, t,  1.43285448e-7f); //  0x1.33b402p-23
    p = fmaf (p, t,  1.22774793e-6f); //  0x1.499232p-20
    p = fmaf (p, t,  1.12963626e-7f); //  0x1.e52cd2p-24
    p = fmaf (p, t, -5.61530760e-5f); // -0x1.d70bd0p-15
    p = fmaf (p, t, -1.47697632e-4f); // -0x1.35be90p-13
    p = fmaf (p, t,  2.31468678e-3f); //  0x1.2f6400p-9
    p = fmaf (p, t,  1.15392581e-2f); //  0x1.7a1e50p-7
    p = fmaf (p, t, -2.32015476e-1f); // -0x1.db2aeep-3
    p = fmaf (p, t,  8.86226892e-1f); //  0x1.c5bf88p-1
  }
  r = z * p;
  return r;
}

int bs_compute_call_iv(bs_data* bsd)
{

  if(!isfinite(bsd->call.IV) || bsd->call.IV <= 0) {
      fprintf(stderr,"%s: Error: Initial IV value must be finite and strictly positive!\n",__func__);
      return -1;
  }
  const double eps=1e-10;
  root_finder* rf=root_finder_init(_bs_compute_call_iv_estimate, bsd);
  double diff;

  //bsd->call.IV = 2.5;
  int ret=root_finder_find(rf, eps, 100, 0, 5, &bsd->call.IV, &diff);

  root_finder_free(rf);

  if(ret) {

    if(ret==-3) fprintf(stderr,"%s: Warning: Convergence seems to have been reached, but the root discrepancy (%22.15e) is larger than required (%22.15e)!\n",__func__,diff,eps);

    else if(ret==-2) {
      fprintf(stderr,"%s: Error: Root could not be found!\n",__func__);
      return -2;
    }
  }
  return 0;
}

int bs_compute_put_iv(bs_data* bsd)
{

  if(!isfinite(bsd->put.IV) || bsd->put.IV <= 0) {
      fprintf(stderr,"%s: Error: Initial IV value must be finite and strictly positive!\n",__func__);
      return -1;
  }
  const double eps=1e-10;
  root_finder* rf=root_finder_init(_bs_compute_put_iv_estimate, bsd);
  double diff;

  //bsd->put.IV = 2.5;
  int ret=root_finder_find(rf, eps, 100, 0, 5, &bsd->put.IV, &diff);

  root_finder_free(rf);

  if(ret) {

    if(ret==-3) fprintf(stderr,"%s: Warning: Convergence seems to have been reached, but the root discrepancy (%22.15e) is larger than required (%22.15e)!\n",__func__,diff,eps);

    else if(ret==-2) {
      fprintf(stderr,"%s: Error: Root could not be found!\n",__func__);
      return -2;
    }
  }
  return 0;
}
