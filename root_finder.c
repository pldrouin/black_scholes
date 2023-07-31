/**
 * @file root_finder.c
 * @brief Finds the root of a function depending on a single parameter.
 * @author <pldrouin@gmail.com>
 */

#include "root_finder.h"

static inline bool root_finder_autosolve(root_finder* rf, double* x, double* diff, double left_bound, bool left_overshoot, double right_bound, bool right_overshoot, uint32_t* rem_count)
{
  double dx;
 
  while(*rem_count) {
    --(*rem_count);
    *x = 0.5 * (left_bound + right_bound);
    dx = rf ->func(*x, diff, rf->params);

    *x += dx;
    printf("%s: %f, %i, %f, %i, %i -> %f\n",__func__,left_bound,left_overshoot,right_bound,right_overshoot,*rem_count,*x);

    if(dx > 0) {

      if(*x < right_bound) return true;

      if(right_overshoot == false) {
	left_bound = *x - dx;
	left_overshoot = false;

      } else {
	right_bound = *x - dx;
	right_overshoot = true;
      }

    } else {

      if(*x > left_bound) return true;

      if(left_overshoot == false) {
	right_bound = *x - dx;
	right_overshoot = false;

      } else {
	left_bound = *x - dx;
	left_overshoot = true;
      }
    }
  }
  return false;
}

int root_finder_find(root_finder* rf, const double eps, const uint32_t maxiter, const double xmin, const double xmax, double* x, double* calcdiff)
{
  printf("%s(%e, %i, %f, %f, %f)\n",__func__,eps,maxiter,xmin,xmax,*x);
  double diff=0;
  uint32_t rem_iter=maxiter;
  uint32_t maxxiter=0;
  uint32_t minxiter=0;
  bool samex=false;
  const double griddx=(xmax-xmin)/maxiter;

  if(!(*x >= xmin)) *x=xmin;
  //printf("x=%22.15e (xmin=%22.15e, xmax=%22.15e)\n",*x,xmin,xmax);
  double dx=rf->func(*x, &diff, rf->params);

  if(dx==0) samex=true;

  else samex=false;

  bool success=true;

  while(!(fabs(diff) < eps) && --rem_iter > 0) {
    *x += dx;
    //printf("x=%22.15e, diff=%22.15e\n",*x,diff);

    if(*x >= xmax) {
      printf("Max hit\n");

      if(success) {

	if(rf->func(xmax, &diff, rf->params) < 0) success=root_finder_autosolve(rf, x, &diff, *x - dx, false, xmax, false, &rem_iter);

	else success=root_finder_autosolve(rf, x, &diff, xmin, (rf->func(xmin, &diff, rf->params) < 0), *x - dx, true, &rem_iter);

	rem_iter=maxiter;
      }

      if(!success) {
	*x = xmax - (++maxxiter)*griddx;

	while(*x > xmin) {
	  dx=rf->func(*x, &diff, rf->params);
	  printf("Returned x is %e\n",*x+dx);

	  if(!(dx <= 0) || !(*x+dx >= xmin)) *x = xmax - (++maxxiter)*griddx;

	  else {
	    *x += dx;
	    break;
	  }
	}
      }
    }

    if(*x <= xmin) {
      printf("Min hit\n");

      if(success) {

	if(rf->func(xmin, &diff, rf->params) > 0) success=root_finder_autosolve(rf, x, &diff, xmin, false, *x - dx, false, &rem_iter);

	else success=root_finder_autosolve(rf, x, &diff, *x - dx, true, xmax, (rf->func(xmax, &diff, rf->params) > 0), &rem_iter);

	rem_iter=maxiter;
      }

      if(!success) {
	*x = xmin + (++minxiter)*griddx;

	while(*x < xmax) {
	  dx=rf->func(*x, &diff, rf->params);
	  printf("Returned x is %e\n",*x+dx);

	  if(!(dx >= 0) || !(*x+dx <= xmax)) *x = xmin + (++minxiter)*griddx;

	  else {
	    *x += dx;
	    break;
	  }
	}
      }
    }

    //printf(" (corrected x=%22.15e)\n",*x);
    dx=rf->func(*x, &diff, rf->params);

    if(dx==0) samex=true;

    else samex=false;
  }

  printf("x is %e, diff is %e\n",*x,diff);

  if(calcdiff) *calcdiff=diff;
  return -samex-2*(!rem_iter);
}
