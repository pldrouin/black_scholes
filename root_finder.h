/**
 * @file root_finder.h
 * @brief Finds the root of a function depending on a single parameter.
 * @author <pldrouin@gmail.com>
 */

#ifndef _ROOT_FINDER_
#define _ROOT_FINDER_

#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * Root finding algorithm data structure
 */
typedef struct
{
  double (*func)(double const x, double* diff, void* params);	//!< Pointer to the user-provided function.
  void* params;						//!< Pointer to the data to be passed to the user-function.
} root_finder;

/**
 * @brief Initialises the root finding algorithm.
 *
 * This function initialises the root finding algorithm and returns a handle to
 * the allocated memory.
 *
 * @param func: User function used to perform parameter estimate.
 * @param params: Pointer to data sent to the user function.
 * @return pointer to the allocated memory.
 */
inline static root_finder* root_finder_init(double (*func)(double const x, double* diff, void* params), void* params){
  root_finder* ret=(root_finder*)malloc(sizeof(root_finder));

  if(ret) {
    ret->func=func;
    ret->params=params;
  }
  return ret;
}

/**
 * @brief Frees the dynamic memory used by the root finding algorithm.
 *
 * @param rf: Pointer to memory allocated by the root finding algorithm.
 */
inline static void root_finder_free(root_finder* rf){free(rf);}

/**
 * @brief Finds the root of a function depending on a single variable.
 *
 * This function searches for the root of the user-provided function.
 *
 * @param rf: handle for the root finding algorithm.
 * @param eps: Stop point for the algorithm. The algorithm stops once the
 * absolute value for the second argument from the user-defined function is
 * smaller than eps.
 * @param maxiter: Maximum number of iterations for the algorithm.
 * @param xmin: Minimum value for the location of the root.
 * @param xmax: Maximum value for the location of the root.
 * @param x: Initial estimate for the root location, and location of the
 * computed root.
 * @return 0 if successful and a non-zero value otherwise.
 */
int root_finder_find(root_finder* rf, const double eps, const uint32_t maxiter, const double xmin, const double xmax, double* x, double* calcdiff);

#endif
