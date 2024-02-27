#include "lattice.h"
#include "gsl/gsl_complex.h"
#include <cmath>

template <class T>
class MyHamil{
private:
	int block;
	T *Hamil;
	double *Hdia;
	int *istart;
	int *column;
public:
	void H_vector_map_complex(T* vectorin, T* vectorout);
	void free_hamil();

	MyHamil(int input_block, int *input_istart, double *input_Hdia, gsl_complex *input_Hamil, int *input_column){
		block = input_block;
		istart = input_istart;
		Hdia = input_Hdia;
		column = input_column;
		Hamil = new std::complex<double>[istart[block]];
		for(int i=0; i<istart[block]; i++)	Hamil[i] = std::complex<double>(GSL_REAL(input_Hamil[i]), GSL_IMAG(input_Hamil[i]));
	}

}; //MyHamil

template<>
void MyHamil<std::complex<double> >::free_hamil(){
	delete [] Hamil;
}

template<>
void MyHamil<std::complex<double> >::H_vector_map_complex(std::complex<double> *vectorin, std::complex<double> *vectorout){
	int i, j;
	std::complex<double> tempG;

#pragma omp parallel for default(none)	\
	private(i,j,tempG) shared(vectorin, vectorout)
	for(j=0; j<block; j++){
		tempG = vectorin[j] * Hdia[j];

		for(i=istart[j]; i<istart[j+1]; i++)
			tempG += vectorin[column[i]] * Hamil[i];
		vectorout[j] = tempG;
	}
	return;
	//print_gscvectord("vectorout", vectorout, block);
}//end of H_vector_map_complex


