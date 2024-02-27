#include "lattice.h"

#ifndef __MyARP__
#define __MyARP__

int check_map(int *istart, int *column, int block);
double arpack_areig(int *istart, double *Hdia, gsl_complex *Hamil, int *column, int block, gsl_complex **vector_out, double *overlap, int *lancount, double *values, double *res, double Tolerance);
double arpack_areig_Hc(int *istart, double *Hdia, gsl_complex *Hamil, int *column, int block, gsl_complex **vector_out, double *overlap, int *lancount, double *values, double *res, double Tolerance, int *nconv);

#endif
