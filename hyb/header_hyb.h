#include "header/main.h"
#include "header/gsl.h"

#ifndef GLOBAL_H
#define GLOBAL_H

#define Nmax	8192
#define B		128
#define ETA		0.15
#define PI		M_PI
#define Lmax	128
#define Amax	128
#define Wn(n, B)	((2*(n)+1)*M_PI/(B))
#define Wr(n, N) (-8. + (n)*16./(N-1)) // e_min=-8, e_max=8
#define N_EIGEN 1024
#define N_DOS 1024
#define lim		128
#define nquad	128

int nearest, Nb;
double  ti[Lmax], wi[Lmax];	//Gaussian node, weight
char FILE_h5[512], Hamilton[512];
double mu, *wn, *w;

//n,m,h: 2,3,4 hopping term 
struct E_params { 
	double n; 
	double m; 
	double h; 
	int near;
};

struct G_params { 
	int n; 
	double complex *W;
	struct E_params e_params;
};

struct complex_func { 
	double complex* Giwn;
	double complex* Diwn;
	double complex* Gw;
	double complex* Dw;
	double complex* Gt;
	double complex* Dt;
	double complex* G_coeff;
	double complex* D_coeff;
};

struct Malloc {
    double* Chi_rank;
    double* Bath_rank;
    double* SortChi_rank;
    double* SortBath_rank;
    long long int* Iter_rank;
    int* Index_rank;
    long long int* SortIter_rank;
    int extrema_index[32];
    int num_extrema;
    int peak_num;
    int* half_index;
    int SortNum_rank;
};

double bisection(double desired, double* arr_x, double* arr_y, int size);

double NNH(double x, void* parmas); //Next Neighbor Hopping 

void pathname(void* params);

double Gre(double x, void *p);

double Gim(double x, void *p);

void calculate_GD(double *W, double complex *G, double complex *D, void*p, char A);

void fermiNdos(void *params);

void hyb(double complex *D, double *bath, char domain);

void L_coeff(double complex* func, double complex* A, double *ti, double *wi);

void LT( double *ti, double complex *A, double complex *Y, int num);

void FT(double complex* b, double* node1, double complex* af, double* node2, double cycle);

/*
void AN(double complex *func, double* x, double complex* A, double *ti, double *wi, char* fname);

void LT( double *ti, double complex *A, double complex *Y, int num);

void FT(double complex* b, double* node1, double complex* af, double* node2, double cycle);

void inverse( double* ti, double complex* A, double complex* y, char* folder);
*/

double integral_imag(double *x, double complex *y, int size);

double abs_sign(double x);

void Gaussian(double a,double b, int len, double *ni, double *wi);

void sort_bath(const gsl_vector *v, double *ep, double *V);

double obj_f(const gsl_vector *v, void *params);

void obj_df(const gsl_vector* v, void* params, gsl_vector *df);

void obj_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);

double derivative(double* x, double complex* y, int index, int size);

int find_extrema(double *x, double complex *y, int* extrema_index);

double randnum(double min,double max);

void generate_V(double* bath, double err_range, double complex* Dw);

void randSplit(int N, int n, int* splits);

void generate_ep(double* bath, int* ex_index, int* half_index, int peak_num);

void init_bath( double* bath, int* ex_index, int* half_index, int peak_num , double complex* Dw);

void FWHM(double* x, double complex* y, int* peak_index, int peak_num, int peak, int* half_index); //Find Half Maximum

void CGM(double *bath_init, double* bath, double* chi, double complex* Diwn, long long int* iter_save);

void print_bath(double *bath, double *chi, int num);

void printT(char* fname);

void swap_int(int* a, int* b);

void swap_double(double* a, double* b);

// Partition the array into two sub-arrays and return the pivot index
int partition(double* arr, int* indices, int low, int high);

// Main quick sort function
void quickSort(double* arr, int* indices, int low, int high);

void sort_test(double *v, double *ep, double *V);

double test_f(double *v, void *params);

void save_freq(hid_t file_id, double* W, double* WN);

void save_data(hid_t file_id, hid_t group1_id, double complex* data, char* dataset_name, int len);

void load_chi(hid_t path_id, int chi_num, double *chi_val);

void process_datasets(hid_t path_id);

void save_bath(double* sdata, double* chi, int snum, long long int* Iter, long long int TotIter, int seed);

int compare_chi(const void *a, const void *b);

#endif
