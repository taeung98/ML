#include <cmath>
#include "arp.h"
#include "gsl/gsl_complex.h"
#include "arscomp.h"
#include "arcomp.h"
#include "areig.h"
#include "lapackc.h"
#include "hc.cpp"

int check_map(int *istart, int *column, int block){
	int i, j, mapnumber=0;
	for(i=0; i<block; i++){
		for(j=istart[i]; j<istart[i+1]; j++){
			if( column[j]<0 || column[j]>block-1 )	printf("wrong column! i=%d, mapnumber=%d, column=%d\n", i, j, column[j]);
			if( j>istart[i] && column[j]<column[j-1] )	printf("decreasing column! i=%d, mapnumber=%d, column=%d\n", i, j, column[j]);
		}
		for(j=istart[i]; column[j]<i && j<istart[i+1]; j++){
			if( j>istart[i] && (column[j] == column[j-1]) ){
				mapnumber--;
				printf("1same column: i=%d, starting from %d, mapnumber=%d, column=%d, previous %d next %d, block=%d\n", i, istart[i], mapnumber, column[j], column[j-1], column[j+1], block);	fflush(stdout);
			}
			mapnumber++;
		}
		if( i == column[j] ){
			printf("diagonal elements are not supposed to be here!\n");
			printf("i=%d, mapnumber=%d\n", i, mapnumber);	fflush(stdout);
		}
		mapnumber++;

		for(; j<istart[i+1]; j++){
			if( j>istart[i] && (column[j] == column[j-1]) ){
				mapnumber--;
				printf("2same column: i=%d, starting from %d, mapnumber=%d, column=%d, previous %d next %d\n", i, istart[i], mapnumber, column[j], column[j-1], column[j+1]);	fflush(stdout);
			}
			mapnumber++;
		}
	}

	return mapnumber;
}

double arpack_areig_Hc(int *istart, double *Hdia, gsl_complex *Hamil, int *column, int block, gsl_complex **vector_out, double *overlap, int *lancount, double *values, double *res, double Tolerance, int *nconv){
	int i;
	double *EigVal;
	arcomplex<double> *EigVec;
	MyHamil<std::complex<double> > A(block, istart, Hdia, Hamil, column);

	EigVal = new double[Nvec];
	EigVec = new arcomplex<double>[block];

 	ARCompStdEig < double, MyHamil<std::complex<double> > >
	dprob(block, Nvec, &A, &MyHamil<std::complex<double> >::H_vector_map_complex, "SR", 2*Nvec+1, Tolerance, 1000*Nvec);

	*nconv = dprob.FindEigenvectors();

	if( *nconv < Nvec ){
		printf("Ar %d/%d:\t", *nconv, Nvec);
		for(i=0; i<*nconv; i++)
			printf("%.18lf\t", real(dprob.Eigenvalue(i)) );
		printf("\n");	fflush(stdout);
	}
	for(i=0; i<*nconv; i++){
		values[i] = real(dprob.Eigenvalue(i));
		for(int j=0; j<block; j++){
			GSL_REAL( vector_out[i][j] ) = real(dprob.RawEigenvector(i)[j]);
			GSL_IMAG( vector_out[i][j] ) = imag(dprob.RawEigenvector(i)[j]);
		}
	}

	for(i=0; i<*nconv; i++){
		res[i] = 0;
		A.H_vector_map_complex(dprob.RawEigenvector(i), EigVec);
		for(int j=0; j<block; j++)
			res[i] += norm(values[i] * dprob.RawEigenvector(i)[j] - EigVec[j]);
	}

	delete [] EigVal;
	delete [] EigVec;
	A.free_hamil();

	return values[0];
}

double arpack_areig(int *istart, double *Hdia, gsl_complex *Hamil, int *column, int block, gsl_complex **vector_out, double *overlap, int *lancount, double *values, double *res, double Tolerance){
	int i, j, mapnumber=0, *arcol, *aristart, nnz = istart[block] + block;
	arcomplex<double> *Hamil_arpack, *EigVal, *EigVec;

	//nnz = check_map(istart, column, block);
	//printf("duplicates = %d, original = %d, nnz = %d\n", istart[block] + block - nnz, istart[block] + block, nnz );	fflush(stdout);
	//nnz = istart[block] + block;

	Hamil_arpack = new arcomplex<double>[nnz];
	aristart = new int[block+1];
	arcol = new int[nnz];
	EigVal = new arcomplex<double>[Nvec];
	EigVec = new arcomplex<double>[Nvec*block];

	for(i=0; i<nnz; i++)	Hamil_arpack[i] = arcomplex<double>(0, 0);


	aristart[0] = 0;
	for(i=0; i<block; i++){
		//printf("i=%d\n", i);	fflush(stdout);
		for(j=istart[i]; column[j]<i && j<istart[i+1]; j++){
			//printf("j=%d, istart[%d] = %d, column[%d] = %d, mapnumber=%d\n", j, i, istart[i], j, column[j], mapnumber);	fflush(stdout);
			if( j>istart[i] && (column[j] == column[j-1]) ){
				nnz--;
				mapnumber--;
			}
			Hamil_arpack[mapnumber] += arcomplex<double>(GSL_REAL(Hamil[j]), -GSL_IMAG(Hamil[j]));
			arcol[mapnumber] = column[j];
			mapnumber++;
		}
		//printf("diag start\n");	fflush(stdout);
		Hamil_arpack[mapnumber] = arcomplex<double>(Hdia[i], 0);
		arcol[mapnumber] = i;
		mapnumber++;
		//printf("diag complet\n");	fflush(stdout);

		for(; j<istart[i+1]; j++){
			if( j>istart[i] && (column[j] == column[j-1]) ){
				nnz--;
				mapnumber--;
			}
			Hamil_arpack[mapnumber] += arcomplex<double>(GSL_REAL(Hamil[j]), -GSL_IMAG(Hamil[j]));
			arcol[mapnumber] = column[j];
			mapnumber++;
		}
		aristart[i+1] = mapnumber;
	}
	if( mapnumber != nnz )	printf("arpack_areig:: inconsistent matrix!\n");
	else			printf("mapnumber = %d\n", mapnumber);				fflush(stdout);

	int nconv = AREig(EigVal, EigVec, block, nnz, Hamil_arpack, arcol, aristart, Nvec, "SR", 2*Nvec+1, Tolerance, 300*Nvec);

	if( nconv < Nvec ){
		printf("Ar %d/%d:\t", nconv, Nvec);
		for(i=0; i<nconv; i++)
			printf("%.18lf\t", real(EigVal[i]));
		printf("\n");	fflush(stdout);
	}
	for(i=0; i<nconv; i++){
		values[i] = real(EigVal[i]);
		for(j=0; j<block; j++){
			GSL_REAL( vector_out[i][j] ) = real(EigVec[block*i+j]);
			GSL_IMAG( vector_out[i][j] ) = imag(EigVec[block*i+j]);
		}
	}
	double Eg = real(EigVal[0]);

	delete [] aristart;
	delete [] arcol;
	delete [] Hamil_arpack;
	delete [] EigVal;
	delete [] EigVec;

	return Eg;


}
