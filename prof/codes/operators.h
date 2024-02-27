#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>
#include <stdlib.h>
#include "matrix.h"

#ifndef _OPERATORS_
#define _OPERATORS_

int init_operator_angular_z(gsl_complex ***sz, int **ipair, int print);
double coeff_ladder(double j, double m, int sign);
int init_pauli(gsl_complex ***pauli, int print);
int init_operator_L_jj(gsl_complex **Lplus, gsl_complex **Lminus, gsl_complex **Lz, gsl_complex **transform, int print);
int init_operator_J_jj(gsl_complex **Lplus, gsl_complex **Lminus, gsl_complex **Lz, int print);
int init_operator_S_jj(gsl_complex **Lplus, gsl_complex **Lminus, gsl_complex **Lz, gsl_complex **transform, int print);
int test_angular_operators(gsl_complex **transform);

#endif
