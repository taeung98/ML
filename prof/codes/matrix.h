#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "f2c.h"
#include "clapack.h"

#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_eigen.h"
#include "lattice.h"

#ifndef _MYMATRIX_
#define _MYMATRIX_

//inline double gsl_complex_mul_r(gsl_complex *z1, gsl_complex *z2);
//inline double gsl_complex_mul_i(gsl_complex *z1, gsl_complex *z2);

extern inline double gsl_complex_mul_r(gsl_complex *z1, gsl_complex *z2)	{return z1->dat[0]*z2->dat[0]-z1->dat[1]*z2->dat[1];}
extern inline double gsl_complex_mul_i(gsl_complex *z1, gsl_complex *z2)	{return z1->dat[0]*z2->dat[1]+z1->dat[1]*z2->dat[0];}

int cross_product(double a1[3], double a2[3], double out[3]);
double normalize_gsc(gsl_complex *vec, int block);
int Gram_Schmidt_gsl(gsl_complex **V, int NN, int block, int all);
int test_GS(gsl_complex **V, int NN, int block);

char **mkmatrixc(int row, int column);
int **mkmatrixi(int row, int column);
int ***mktritensori(int row, int column, int third);
typebasis **mkmatrixb(int row, int column);
double **mkmatrixd(int row, int column);
void fprint_gscmatrixd(FILE *fp, char *comment, gsl_complex **matrix, int Nx, int Ny);
void fprint_gscmatrixd_math(FILE *fp, char *comment, gsl_complex **matrix, int Nx, int Ny);
void print_gscmatrixd_math(char *comment, gsl_complex **matrix, int Nx, int Ny);
gsl_complex **mkgscmatrixd(int row, int column);
double ***mktritensord(int row, int column, int third);
double ****mktetratensord(int row, int column, int third, int fourth);
gsl_complex  ***mkgsctritensord(int row, int column, int third);
gsl_complex ****mkgsctetratensord(int row, int column, int third, int fourth);
double *****mkpentatensord(int row, int column, int third, int fourth, int fifth);
gsl_complex *****mkgscpentatensord(int row, int column, int third, int fourth, int fifth);
gsl_complex  ******mkgschexatensord(int row, int column, int third, int fourth, int fifth, int sixth);
char *mkvectorc(int NN);
unsigned long *mkvectorul(unsigned long NN);
typebasis *mkvectorb(int NN);
int *mkvectori(int NN);
double *mkvectord(int NN);
gsl_complex *mkgscvectord(int NN);
void initmatrixi(int **matrix, int row, int column);
void initvectori(int *vector, int NN);
void initvectorui(unsigned int *vector, int NN);
void initvectorb(typebasis *vector, int NN);
void initmatrixd(double **matrix, int row, int column);
void initvectord(double *vector, int NN);
void initgscmatrixd(gsl_complex **matrix, int row, int column);
void initgscvectord(gsl_complex *vector, int NN);
void freematrixd(double **matrix, int row);
void freetritensori(int ***tensor, int row, int col);
void freetetratensori(int ****tensor, int row, int col, int third);
void freetritensord(double ***tensor, int row, int col);
void freetetratensord(double ****tensor, int row, int col, int third);
void freegscmatrixd(gsl_complex **matrix, int row);
void freegsctritensord(gsl_complex ***tensor, int row, int col);
void freegsctetratensord(gsl_complex ****tensor, int row, int col, int third);
void freepentatensord(double *****tensor, int row, int col, int third, int fourth);
void freegscpentatensord(gsl_complex *****tensor, int row, int col, int third, int fourth);
void freegschexatensord(gsl_complex ******tensor, int row, int col, int third, int fourth, int sixth);
void freematrixc(char **matrix, int row);
void freematrixi(int **matrix, int row);
void freematrixui(unsigned int **matrix, int row);
void freematrixb(typebasis **matrix, int row);
void print_gscvectord(char *comment, gsl_complex *vector, int Nx);
void print_gscmatrixd_compact(char *comment, gsl_complex **matrix, int Nx, int Ny);
void print_gscmatrixd(char *comment, gsl_complex **matrix, int Nx, int Ny);
int print_gsctetratensord_nonzero(char *comment, gsl_complex ****matrix, int row, int column, int third, int fourth);
void print_vectorc(char *comment, char *vector, int Nx);
void print_vectorui(char *comment, unsigned int *vector, int Nx);
void print_vectori(char *comment, int *vector, int Nx);
void print_vectord(char *comment, double *vector, int Nx);
void print_matrixi(char *comment, int **matrix, int Nx, int Ny);
void print_matrixd(char *comment, double **matrix, int Nx, int Ny);
int print_tetratensord_nonzero(char *comment, double ****matrix, int row, int column, int third, int fourth);
void inverse_matrixd_lapack(double **input, double **output, int NN);
void inverse_complex_matrix_lapack(gsl_complex **input, gsl_complex **output, int NN);
void inverse_complex_matrix_all_omp(gsl_complex ***input, gsl_complex ***output, int mystart, int myend, int NN);
int check_unitary(gsl_complex **transform, int NN);
int iunitary(gsl_complex **original, gsl_complex **target, gsl_complex **transform, int NN);
int similarity(gsl_complex **original, gsl_complex **target, double **transform, int NN);
int unitary(gsl_complex **original, gsl_complex **target, gsl_complex **transform, int NN);
int gsc_mul(gsl_complex **left, gsl_complex **right, gsl_complex **output, int NN);
int eigen_lapack(gsl_complex **matrix, double *eval, gsl_complex **evec, int NN, int ASC);
int eigen_gsl(gsl_complex **matrix, double *eval, gsl_complex **evec, int NN);

#endif
