#include <stdio.h>
#ifndef __Hamil__
#define __Hamil__

#include "ls_basis.h"
#include "lattice.h"
#include "matrix.h"

extern double Uinter, SOC, J;
extern int myrank, nctot, nb, Njj, Ndd;
extern int Powns, Powns2, Spinmax, Blocks, Nonzero_SOC;
extern gsl_complex zero, **Ecluster;
extern long ranseed;

float ran2(long *idum);

int construct_hamil_direct(int *istart, double *Hdia, gsl_complex *Hamil, int *column, gsl_complex **matrix, int block);
double LiuDavidson(int *istart, double *Hdia, gsl_complex *Hamil, int *column, int block, gsl_complex **vector_in0, gsl_complex **vector_out, double *overlap, int *lancount, double *values, double *res, double Tolerance);
double Davidson   (int *istart, double *Hdia, gsl_complex *Hamil, int *column, int block, gsl_complex **vector_in0, gsl_complex **vector_out, double *overlap, int *lancount, double *values, double *res, double Tolerance);
void H_vector_map_complex(int *istart, double *Hdia, gsl_complex *Hamil, int *column, gsl_complex *vectorin, gsl_complex *vectorout, int block);
int print_map(char *save_directory, int count, int n, typebasis **basis, int *istart, sectorinfo ginfo[], int gindex, PNSpara egv);
int read_map(int n, typebasis **bitbasis, int *istart, gsl_complex *Hamil, int *column, sectorinfo ginfo[], int gindex, PNSpara egv);
int compute_Hdia(int n, double *Hdia, sectorinfo ginfo[], int gindex, typebasis **basis, PNSpara egv);

double MLanczos_Nvector(int *istart, double *Hdia, gsl_complex *Hamil, int *column, int NN, gsl_complex *vector_in0, gsl_complex *vector_out0, double *overlap, int *lancount, double *values, double *res);
gsl_complex **lanczos_vector_complex(	int *istart, double *Hdia, gsl_complex *Hamil, int *column, int *dimension, gsl_complex *p0, double ai[], double bi[], int grblock, int stdlaniter, int Vsave);

int translate(Coulomb *DD, int Ndd);
int init_dd(double U, double J, Coulomb *DD);
int init_others(double J, Coulomb *JJ);

int build_mapinfo(gsl_complex **Ecluster, typebasis **bitbasis, mapele *mapinfo, int *istart, sectorinfo *ginfo, int gindex);
int construct_mapnumber(gsl_complex **Ecluster, typebasis **bitbasis, sectorinfo *ginfo, int init, int final, int *Ndd, int *Njj);

Coulomb *mkvectorC(int NN);
mapele *mkvectormapele(int mapnumber);

#endif
