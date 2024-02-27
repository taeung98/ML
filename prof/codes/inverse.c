#include <gsl/gsl_complex.h>
#include <omp.h>
#include "matrix.h"

extern int nctot;

void inverse_complex_matrix_22(gsl_complex input[2][2], gsl_complex output[2][2]){
	gsl_complex deter;
	deter = gsl_complex_sub(
			gsl_complex_mul(
				input[0][0],
				input[1][1]
			),
			gsl_complex_mul(
				input[0][1],
				input[1][0]
			)
		);
	output[0][0] = gsl_complex_div( input[1][1], deter );
	output[1][1] = gsl_complex_div( input[0][0], deter );
	output[1][0]
		= gsl_complex_div(
			gsl_complex_mul_real( input[0][1], -1 ),
			deter
		);
	output[0][1]
		= gsl_complex_div(
			gsl_complex_mul_real( input[1][0], -1 ),
			deter
		);
}

void inverse_complex_matrix_lapack(gsl_complex **input, gsl_complex **output, int NN){
        int mu, nu;

	integer N = NN;
	integer LWORK = N*N+N;

	doublecomplex *temp = (doublecomplex*) malloc ( N*N*sizeof(doublecomplex) );
	doublecomplex *WORK = (doublecomplex*) malloc ( LWORK*sizeof(doublecomplex) );
	integer *perm = (integer*) malloc ( 2*N*sizeof(integer) );

	integer INFO;
	for(mu=0; mu<N; mu++) for(nu=0; nu<N; nu++){
		temp[mu*N+nu].r = GSL_REAL(input[nu][mu]);
		temp[mu*N+nu].i = GSL_IMAG(input[nu][mu]);
	}
	zgetrf_(&N, &N, temp, &N, perm, &INFO );
	if( (int)INFO != 0 ){
		printf("Error in zgetrf: INFO = %d (inverse_complex_matrix_lapack)\n", (int)INFO);
		print_gscmatrixd("input", input, NN, NN);
		exit(1);
	}
	zgetri_(&N, temp, &N, perm, WORK, &LWORK, &INFO);
	if( (int)INFO != 0 ){
		printf("Error in zgetri: INFO = %d (inverse_complex_matrix_lapack)\n", (int)INFO);
		print_gscmatrixd("input", input, NN, NN);
		exit(1);
	}
	for(mu=0; mu<N; mu++) for(nu=0; nu<N; nu++){
		GSL_SET_COMPLEX(&output[nu][mu], temp[mu*N+nu].r, temp[mu*N+nu].i );
	}
	free(temp); free(WORK); free(perm);
}

void inverse_complex_matrix_all_omp(gsl_complex ***input, gsl_complex ***output, int mystart, int myend, int NN){
        int i, k, mu, nu;
	int Omp_num = omp_get_max_threads();

	integer N = NN;
	integer LWORK = N*N+N;

	doublecomplex **temp = (doublecomplex**) malloc( Omp_num* sizeof(doublecomplex*));
	doublecomplex **WORK = (doublecomplex**) malloc( Omp_num* sizeof(doublecomplex*));
	integer **perm = (integer**) malloc( Omp_num* sizeof(integer*));
	for(k=0; k<Omp_num; k++){
		temp[k] = (doublecomplex*) malloc ( N*N*sizeof(doublecomplex) );
		WORK[k] = (doublecomplex*) malloc ( LWORK*sizeof(doublecomplex) );
		perm[k] = (integer*) malloc ( 2*N*sizeof(integer) );
	}

#pragma omp parallel for default(none)	\
	private(i,mu,nu) shared(mystart, myend, input, N, output, temp, perm, WORK, LWORK, NN)
	for(i=mystart; i<myend; i++){
		integer INFO;
		int mythread = omp_get_thread_num();
		for(mu=0; mu<N; mu++) for(nu=0; nu<N; nu++){
			temp[mythread][mu*N+nu].r = GSL_REAL(input[i][nu][mu]);
			temp[mythread][mu*N+nu].i = GSL_IMAG(input[i][nu][mu]);
		}
		zgetrf_(&N, &N, temp[mythread], &N, perm[mythread], &INFO );
		if( (int)INFO != 0 ){
			printf("Error in zgetrf: INFO = %d (inverse_complex_matrix_all_omp)\n", (int)INFO);
			print_gscmatrixd("input[i]", input[i], NN, NN);
			exit(1);
		}
		zgetri_(&N, temp[mythread], &N, perm[mythread], WORK[mythread], &LWORK, &INFO);
		if( (int)INFO != 0 ){
			printf("Error in zgetri: INFO = %d (inverse_complex_matrix_all_omp)\n", (int)INFO);
			print_gscmatrixd("input[i]", input[i], NN, NN);
			exit(1);
		}
		for(mu=0; mu<N; mu++) for(nu=0; nu<N; nu++){
			GSL_SET_COMPLEX(&output[i][nu][mu], temp[mythread][mu*N+nu].r, temp[mythread][mu*N+nu].i );
		}
	}
	for(k=0; k<Omp_num; k++){
		free(temp[k]);
		free(WORK[k]);
		free(perm[k]);
	}
	free(temp); free(WORK); free(perm);
}
