#include <omp.h>
#include "matrix.h"

//inline double gsl_complex_mul_r(gsl_complex *z1, gsl_complex *z2)	{return z1->dat[0]*z2->dat[0]-z1->dat[1]*z2->dat[1];}
//inline double gsl_complex_mul_i(gsl_complex *z1, gsl_complex *z2)	{return z1->dat[0]*z2->dat[1]+z1->dat[1]*z2->dat[0];}

int cross_product(double a1[3], double a2[3], double out[3]){
	out[0] = a1[1]*a2[2] - a1[2]*a2[1];
	out[1] = a1[2]*a2[0] - a1[0]*a2[2];
	out[2] = a1[0]*a2[1] - a1[1]*a2[0];
	return 0;
}

double normalize_gsc(gsl_complex *vec, int block){
	int i;
	double sum=0;

#pragma omp parallel for default(none) private(i) shared(block, vec) reduction(+:sum)
	for(i=0; i<block; i++)
		sum += gsl_complex_abs2( vec[i] );
	sum = sqrt(sum);
#pragma omp parallel for default(none) private(i) shared(block, vec, sum)
	for(i=0; i<block; i++)
		vec[i] = gsl_complex_div_real( vec[i], sum );

	return sum;
}

int Gram_Schmidt_gsl(gsl_complex **V, int NN, int block, int all){
	int kmin, k, i, j;
	double inn_r, inn_i, norm;
	gsl_complex inn;
	if( NN < 1 ){
		printf("Gram_Schmidt: NN should be larger than 0\n");
		exit(1);
	}
	else if( NN == 1 ){
		normalize_gsc(V[0], block);

		return 0;
	}
	if( all ){
		kmin = 1;
		for(k=0; k<NN; k++){
			norm = 0;
#pragma omp parallel for default(none) private(i) shared(k, block, V) reduction(+:norm)
			for(i=0; i<block; i++){
				norm += gsl_complex_inner( V[k][i], V[k][i] );
			}
			norm = sqrt(norm);
#pragma omp parallel for default(none) private(i) shared(k, block, V, norm)
			for(i=0; i<block; i++)
				V[k][i] = gsl_complex_div_real( V[k][i], norm);
		}
	}
	else		kmin = (NN>1) ? NN-1 : 0;

	for(k=kmin; k<NN; k++){
		for(j=0; j<k; j++) {
			inn_r = 0; inn_i = 0;
#pragma omp parallel for default(none) private(i) shared(k, j, block, V) reduction(+:inn_r, inn_i)
			for(i=0; i<block; i++){
				inn_r += gsl_complex_in_r( V[j][i], V[k][i] );
				inn_i += gsl_complex_in_i( V[j][i], V[k][i] );
			}
			inn = gsl_complex_rect( inn_r, inn_i );
#pragma omp parallel for default(none) private(i) shared(inn, k, j, block, V)
			for(i=0; i<block; i++)
				V[k][i] = gsl_complex_sub( V[k][i], gsl_complex_mul( V[j][i], inn ));
		}
		normalize_gsc(V[k], block);
		//if( normalize_gsc(V[k], block) < 1e-12 ) break;
	}
	return 0;
}

int test_GS(gsl_complex **V, int NN, int block){
	int j,k,i;
	double norm, **over;
	over = mkmatrixd(NN, NN);
	for(k=0; k<NN; k++){
		norm = 0;
#pragma omp parallel for default(none) private(i) shared(V, block, k) reduction(+:norm)
		for(i=0; i<block; i++)
			norm += gsl_complex_abs2(V[k][i]);
		if( fabs(norm - 1) > 1e-6 )
			printf("k=%d, norm =%.10lf\n", k, norm);
	}
	for(j=0; j<NN; j++) for(k=0; k<NN; k++){
		norm = 0;
#pragma omp parallel for default(none) private(i) shared(V, block, j, k) reduction(+:norm)
		for(i=0; i<block; i++)
			norm += gsl_complex_inner(V[j][i], V[k][i]);
		over[j][k] = norm;
	}
	for(j=0; j<NN; j++) for(k=j+1; k<NN; k++)
		if( fabs(over[j][k]) > 1e-6  ) {
			printf("test_GS: Large overlap!\n");
			print_matrixd("over", over, NN, NN);
			return 1;
		}
	freematrixd(over, NN);
	return 0;
}


int my_alloc_error(char *name){
	printf("Allocation error in %s\n", name);
	exit(1);
}

char **mkmatrixc(int row, int column){
	int i;
	char **matrix;
	matrix = (char **) malloc ( row * sizeof(char*) );
	for(i=0; i<row; i++)	matrix[i] = (char *) malloc ( column * sizeof(char) );
	if( matrix == NULL )	my_alloc_error("mkmatrixi");
	return matrix;
}

int **mkmatrixi(int row, int column){
	int i, **matrix;
	matrix = (int **) malloc ( row * sizeof(int*) );
	for(i=0; i<row; i++)	matrix[i] = (int *) malloc ( column * sizeof(int) );
	if( matrix == NULL )	my_alloc_error("mkmatrixi");
	return matrix;
}

int ***mktritensori(int row, int column, int third){
	int i, j;
	int ***tensor;
	tensor = (int ***) malloc ( row * sizeof(int**) );
	for(i=0; i<row; i++){
		tensor[i] = (int **) malloc ( column * sizeof(int*) );
		for(j=0; j<column; j++)	tensor[i][j] = (int *) malloc (  third * sizeof(int) );
	}
	if( tensor == NULL )	my_alloc_error("mktritensori");
	return tensor;
}
typebasis **mkmatrixb(int row, int column){
	int i, j;
	typebasis **matrix;
	matrix = (typebasis **) malloc ( row * sizeof(typebasis*) );
	for(i=0; i<row; i++)	matrix[i] = (typebasis *) malloc ( column * sizeof(typebasis) );
	for(i=0; i<row; i++) for(j=0; j<column; j++)	matrix[i][j] = Basis_One;
	if( matrix == NULL )	my_alloc_error("mkmatrixb");
	return matrix;
}

double **mkmatrixd(int row, int column){
	int i;
	double **matrix;
	matrix = (double **) malloc ( row * sizeof(double*) );
	for(i=0; i<row; i++)	matrix[i] = (double *) malloc ( column * sizeof(double) );
	if( matrix == NULL )	my_alloc_error("mkmatrixd");
	return matrix;
}

void fprint_gscmatrixd(FILE *fp, char *comment, gsl_complex **matrix, int Nx, int Ny){
	int i, j;
	fprintf(fp, "<:: %s ::> ( %d x %d )--------------\n", comment, Nx, Ny);
	for(i=0; i<Nx; i++){
		for(j=0; j<Ny; j++)
			fprintf(fp, "%13.10lf, %13.10lf    ", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
		fprintf(fp, "\n");
	}
	fprintf(fp, "------------------------------\n\n");
}
void fprint_gscmatrixd_math(FILE *fp, char *comment, gsl_complex **matrix, int Nx, int Ny){
	int i, j;
	fprintf(fp, "<:: %s ::> ( %d x %d )--------------\n", comment, Nx, Ny);
	fprintf(fp, "{\n");
	for(i=0; i<Nx-1; i++){
		fprintf(fp, "{ ");
		for(j=0; j<Ny-1; j++)
			if(gsl_complex_abs2(matrix[i][j])>1e-16)
				if( fabs(GSL_IMAG(matrix[i][j]))<1e-16 )
					fprintf(fp, "%21.18lf , ", GSL_REAL(matrix[i][j]));
				else
					fprintf(fp, "%21.18lf + I * %21.18lf, ", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
			else    fprintf(fp, "0, ");
		if(gsl_complex_abs2(matrix[i][j])>1e-16)
			if( fabs(GSL_IMAG(matrix[i][j]))<1e-16 )
				fprintf(fp, "%21.18lf },\n", GSL_REAL(matrix[i][j]));
			else
				fprintf(fp, "%21.18lf + I * %21.18lf },\n", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
		else    fprintf(fp, "0 },\n");
	}
	fprintf(fp, "{ ");
	for(j=0; j<Ny-1; j++)
		if(gsl_complex_abs2(matrix[i][j])>1e-16)
			if( fabs(GSL_IMAG(matrix[i][j]))<1e-16 )
				fprintf(fp, "%21.18lf , ", GSL_REAL(matrix[i][j]));
			else
				fprintf(fp, "%21.18lf + I * %21.18lf, ", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
		else    fprintf(fp, "0, ");
	if(gsl_complex_abs2(matrix[i][j])>1e-16)
		if( fabs(GSL_IMAG(matrix[i][j]))<1e-16 )
			fprintf(fp, "%21.18lf }\n}\n", GSL_REAL(matrix[i][j]));
		else
			fprintf(fp, "%21.18lf + I * %21.18lf }\n}\n", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
	else    fprintf(fp, "0 }\n}\n");
	fprintf(fp, "------------------------------\n\n");
}

void print_gscmatrixd_math(char *comment, gsl_complex **matrix, int Nx, int Ny){
	int i, j;
	printf("<:: %s ::>--------------\n", comment);
	printf("{\n");
	for(i=0; i<Nx-1; i++){
		printf("{ ");
		for(j=0; j<Ny-1; j++)
			if(gsl_complex_abs2(matrix[i][j])>1e-16)
				if( fabs(GSL_IMAG(matrix[i][j]))<1e-16 )
					printf("%21.18lf , ", GSL_REAL(matrix[i][j]));
				else
					printf("%21.18lf + I * %21.18lf, ", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
			else    printf("0, ");
		if(gsl_complex_abs2(matrix[i][j])>1e-16)
			if( fabs(GSL_IMAG(matrix[i][j]))<1e-16 )
				printf("%21.18lf },\n", GSL_REAL(matrix[i][j]));
			else
				printf("%21.18lf + I * %21.18lf },\n", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
		else    printf("0 },\n");
	}
	printf("{ ");
	for(j=0; j<Ny-1; j++)
		if(gsl_complex_abs2(matrix[i][j])>1e-16)
			if( fabs(GSL_IMAG(matrix[i][j]))<1e-16 )
				printf("%21.18lf , ", GSL_REAL(matrix[i][j]));
			else
				printf("%21.18lf + I * %21.18lf, ", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
		else    printf("0, ");
	if(gsl_complex_abs2(matrix[i][j])>1e-16)
		if( fabs(GSL_IMAG(matrix[i][j]))<1e-16 )
			printf("%21.18lf }\n}\n", GSL_REAL(matrix[i][j]));
		else
			printf("%21.18lf + I * %21.18lf }\n}\n", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
	else    printf("0 }\n}\n");
	printf("------------------------------\n\n");
}

gsl_complex **mkgscmatrixd(int row, int column){
	int i;
	gsl_complex **matrix;
	matrix = (gsl_complex **) malloc ( row * sizeof(gsl_complex*) );
	for(i=0; i<row; i++)	matrix[i] = (gsl_complex *) malloc ( column
			* sizeof(gsl_complex) );
	if( matrix == NULL )	my_alloc_error("mkgscmatrixd");
	return matrix;
}

double ***mktritensord(int row, int column, int third){
	int i, j;
	double ***tensor;
	tensor = (double ***) malloc ( row * sizeof(double**) );
	for(i=0; i<row; i++){
		tensor[i] = (double **) malloc ( column * sizeof(double*) );
		for(j=0; j<column; j++)	tensor[i][j] = (double *) malloc (  third * sizeof(double) );
	}
	if( tensor == NULL )	my_alloc_error("mktritensord");
	return tensor;
}
double ****mktetratensord(int row, int column, int third, int fourth){
	int i, j, k;
	double ****tensor;
	tensor = (double ****) malloc ( row * sizeof(double***) );
	for(i=0; i<row; i++){
		tensor[i] = (double ***) malloc ( column * sizeof(double**) );
		for(j=0; j<column; j++){
			tensor[i][j] = (double **) malloc (  third * sizeof(double*) );
			for(k=0; k<third; k++)
				tensor[i][j][k] = (double *) malloc ( fourth * sizeof(double) );
		}
	}
	if( tensor == NULL )	my_alloc_error("mktetratensord");
	return tensor;
}

gsl_complex  ***mkgsctritensord(int row, int column, int third){
	int i, j;
	gsl_complex  ***tensor;
	tensor = (gsl_complex ***) malloc ( row * sizeof(gsl_complex**) );
	for(i=0; i<row; i++){
		tensor[i] = (gsl_complex **) malloc ( column * sizeof(gsl_complex*) );
		for(j=0; j<column; j++)	tensor[i][j] = (gsl_complex *) malloc (  third * sizeof(gsl_complex) );
	}
	if( tensor == NULL )	my_alloc_error("mkgsctritensord");
	return tensor;
}

gsl_complex ****mkgsctetratensord(int row, int column, int third, int fourth){
	int i, j, k;
	gsl_complex  ****tensor;
	tensor = (gsl_complex ****) malloc ( row * sizeof(gsl_complex***) );
	for(i=0; i<row; i++){
		tensor[i] = (gsl_complex ***) malloc ( column * sizeof(gsl_complex**) );
		for(j=0; j<column; j++){
			tensor[i][j] = (gsl_complex **) malloc (  third * sizeof(gsl_complex*) );
			for(k=0; k<third; k++)	tensor[i][j][k] = (gsl_complex *) malloc ( fourth * sizeof(gsl_complex) );
		}
	}
	if( tensor == NULL )	my_alloc_error("mkgsctetratensord");
	return tensor;
}
double *****mkpentatensord(int row, int column, int third, int fourth, int fifth){
	int i, j, k, l;
	double  *****tensor;
	tensor = (double *****) malloc ( row * sizeof(double****) );
	for(i=0; i<row; i++){
		tensor[i] = (double ****) malloc ( column * sizeof(double***) );
		for(j=0; j<column; j++){
			tensor[i][j] = (double ***) malloc (  third * sizeof(double**) );
			for(k=0; k<third; k++){
				tensor[i][j][k] = (double **) malloc ( fourth * sizeof(double*) );
				for(l=0; l<fourth; l++)
					tensor[i][j][k][l] = (double *) malloc ( fifth * sizeof(double) );
			}
		}
	}
	if( tensor == NULL )	my_alloc_error("mkpentatensord");
	return tensor;
}
gsl_complex *****mkgscpentatensord(int row, int column, int third, int fourth, int fifth){
	int i, j, k, l;
	gsl_complex  *****tensor;
	tensor = (gsl_complex *****) malloc ( row * sizeof(gsl_complex****) );
	for(i=0; i<row; i++){
		tensor[i] = (gsl_complex ****) malloc ( column * sizeof(gsl_complex***) );
		for(j=0; j<column; j++){
			tensor[i][j] = (gsl_complex ***) malloc (  third * sizeof(gsl_complex**) );
			for(k=0; k<third; k++){
				tensor[i][j][k] = (gsl_complex **) malloc ( fourth * sizeof(gsl_complex*) );
				for(l=0; l<fourth; l++)
					tensor[i][j][k][l] = (gsl_complex *) malloc ( fifth * sizeof(gsl_complex) );
			}
		}
	}
	if( tensor == NULL )	my_alloc_error("mkgscpentatensord");
	return tensor;
}

gsl_complex  ******mkgschexatensord(int row, int column, int third, int fourth, int fifth, int sixth){
	int i, j, k, l, m;
	gsl_complex  ******tensor;
	tensor = (gsl_complex ******) malloc ( row * sizeof(gsl_complex*****) );
	for(i=0; i<row; i++){
		tensor[i] = (gsl_complex *****) malloc ( column * sizeof(gsl_complex****) );
		for(j=0; j<column; j++){
			tensor[i][j] = (gsl_complex ****) malloc (  third * sizeof(gsl_complex***) );
			for(k=0; k<third; k++){
				tensor[i][j][k] = (gsl_complex ***) malloc ( fourth * sizeof(gsl_complex**) );
				for(l=0; l<fourth; l++){
					tensor[i][j][k][l] = (gsl_complex **) malloc ( fifth * sizeof(gsl_complex*) );
					for(m=0; m<fifth; m++)
						tensor[i][j][k][l][m] = (gsl_complex *) malloc ( sixth * sizeof(gsl_complex) );
				}
			}
		}
	}
	if( tensor == NULL )	my_alloc_error("mkgschexatensord");
	return tensor;
}

char *mkvectorc(int NN){
	char *vector;
	vector = (char *) malloc ( NN * sizeof(char) );
	if( vector == NULL )	my_alloc_error("mkvectorc");
	return vector;
}

unsigned long *mkvectorul(unsigned long NN){
	unsigned long *vector;
	vector = (unsigned long *) malloc ( NN * sizeof(unsigned long) );
	if( vector == NULL )	my_alloc_error("mkvectorul");
	return vector;
}
typebasis *mkvectorb(int NN){
	typebasis *vector;
	vector = (typebasis *) malloc ( NN * sizeof(typebasis) );
	if( vector == NULL )	my_alloc_error("mkvectorb");
	return vector;
}

int *mkvectori(int NN){
	int *vector;
	vector = (int *) malloc ( NN * sizeof(int) );
	if( vector == NULL )	my_alloc_error("mkvectori");
	return vector;
}

double *mkvectord(int NN){
	double *vector;
	vector = (double *) malloc ( NN * sizeof(double) );
	if( vector == NULL )	my_alloc_error("mkvectord");
	return vector;
}

gsl_complex *mkgscvectord(int NN){
	gsl_complex *vector;
	vector = (gsl_complex *) malloc ( NN * sizeof(gsl_complex) );
	if( vector == NULL )	my_alloc_error("mkgscvectord");
	return vector;
}

void initmatrixi(int **matrix, int row, int column){
	int i,j;
	for(i=0; i<row; i++) for(j=0; j<column; j++)	matrix[i][j] = 0;
}

void initvectori(int *vector, int NN){
	int i;
	for(i=0; i<NN; i++)	vector[i] = 0;
}

void initvectorui(unsigned int *vector, int NN){
	int i;
	for(i=0; i<NN; i++)	vector[i] = 0;
}

void initvectorb(typebasis *vector, int NN){
	int i;
	for(i=0; i<NN; i++)	vector[i] = 0ULL;
}
void initmatrixd(double **matrix, int row, int column){
	int i,j;
	for(i=0; i<row; i++) for(j=0; j<column; j++)	matrix[i][j] = 0;
}

void initvectord(double *vector, int NN){
	int i;
	for(i=0; i<NN; i++)	vector[i] = 0;
}

void initgscmatrixd(gsl_complex **matrix, int row, int column){
	int i,j;
	for(i=0; i<row; i++) for(j=0; j<column; j++)	GSL_SET_COMPLEX(&matrix[i][j],0,0) ;
}
void initgscvectord(gsl_complex *vector, int NN){
	int i;
	for(i=0; i<NN; i++)	GSL_SET_COMPLEX(&vector[i],0,0) ;
}

void freematrixd(double **matrix, int row){
	int i;
	for(i=0; i<row; i++)	free(matrix[i]);
	free(matrix);
}

void freetritensori(int ***tensor, int row, int col){
	int i, j;
	for(i=0; i<row; i++)	{
		for (j=0;j<col;j++) free(tensor[i][j]);
		free(tensor[i]);
	}
	free(tensor);
}

void freetetratensori(int ****tensor, int row, int col, int third){
	int i, j, k;
	for(i=0; i<row; i++)	{
		for (j=0;j<col;j++){
			for (k=0;k<third;k++)	free( tensor[i][j][k] );
			free(tensor[i][j]);
		}
		free(tensor[i]);
	}
	free(tensor);
}
void freetritensord(double ***tensor, int row, int col){
	int i, j;
	for(i=0; i<row; i++)	{
		for (j=0;j<col;j++) free(tensor[i][j]);
		free(tensor[i]);
	}
	free(tensor);
}

void freetetratensord(double ****tensor, int row, int col, int third){
	int i, j, k;
	for(i=0; i<row; i++)	{
		for (j=0;j<col;j++){
			for (k=0;k<third;k++)	free( tensor[i][j][k] );
			free(tensor[i][j]);
		}
		free(tensor[i]);
	}
	free(tensor);
}


void freegscmatrixd(gsl_complex **matrix, int row){
	int i;
	for(i=0; i<row; i++)	free(matrix[i]);
	free(matrix);
}

void freegsctritensord(gsl_complex ***tensor, int row, int col){
	int i, j;
	for(i=0; i<row; i++)	{
		for (j=0;j<col;j++) free(tensor[i][j]);
		free(tensor[i]);
	}
	free(tensor);
}
void freegsctetratensord(gsl_complex ****tensor, int row, int col, int third){
	int i, j, k;
	for(i=0; i<row; i++)	{
		for (j=0;j<col;j++){
			for (k=0;k<third;k++)	free( tensor[i][j][k] );
			free(tensor[i][j]);
		}
		free(tensor[i]);
	}
	free(tensor);
}

void freepentatensord(double *****tensor, int row, int col, int third, int fourth){
	int i, j, k, l;
	for(i=0; i<row; i++)	{
		for (j=0;j<col;j++){
			for (k=0;k<third;k++){
				for(l=0; l<fourth; l++)	free( tensor[i][j][k][l] );
				free( tensor[i][j][k] );
			}
			free(tensor[i][j]);
		}
		free(tensor[i]);
	}
	free(tensor);
}

void freegscpentatensord(gsl_complex *****tensor, int row, int col, int third, int fourth){
	int i, j, k, l;
	for(i=0; i<row; i++)	{
		for (j=0;j<col;j++){
			for (k=0;k<third;k++){
				for(l=0; l<fourth; l++)	free( tensor[i][j][k][l] );
				free( tensor[i][j][k] );
			}
			free(tensor[i][j]);
		}
		free(tensor[i]);
	}
	free(tensor);
}

void freegschexatensord(gsl_complex ******tensor, int row, int col, int third, int fourth, int sixth){
	int i, j, k, l, m;
	for(i=0; i<row; i++){
		for (j=0;j<col;j++){
			for (k=0;k<third;k++){
				for(l=0; l<fourth; l++){
					for(m=0; m<sixth; m++)	free( tensor[i][j][k][l][m] );
					free( tensor[i][j][k][l] );
				}
				free( tensor[i][j][k] );
			}
			free(tensor[i][j]);
		}
		free(tensor[i]);
	}
	free(tensor);
}

void freematrixc(char **matrix, int row){
	int i;
	for(i=0; i<row; i++)	free(matrix[i]);
	free(matrix);
}
void freematrixi(int **matrix, int row){
	int i;
	for(i=0; i<row; i++)	free(matrix[i]);
	free(matrix);
}
void freematrixui(unsigned int **matrix, int row){
	int i;
	for(i=0; i<row; i++)	free(matrix[i]);
	free(matrix);
}
void freematrixb(typebasis **matrix, int row){
	int i;
	for(i=0; i<row; i++)	free(matrix[i]);
	free(matrix);
}

void print_vectord(char *comment, double *vector, int Nx){
	int i;
	printf("<:: %s ::> ( %d )--------------\n", comment, Nx);
	for(i=0; i<Nx; i++){
		printf("%10.7f  ", vector[i] );
	}
	printf("\n------------------------------\n\n");
}

void print_vectorui(char *comment, unsigned int *vector, int Nx){
	int i;
	printf("<:: %s ::> ( %d )--------------\n", comment, Nx);
	for(i=0; i<Nx; i++){
		printf("%10u  ", vector[i] );
	}
	printf("\n------------------------------\n\n");
}

void print_vectorc(char *comment, char *vector, int Nx){
	int i;
	printf("<:: %s ::> ( %d )--------------\n", comment, Nx);
	for(i=0; i<Nx; i++){
		printf("%c  ", vector[i] );
	}
	printf("\n------------------------------\n\n");
}

void print_vectori(char *comment, int *vector, int Nx){
	int i;
	printf("<:: %s ::> ( %d )--------------\n", comment, Nx);
	for(i=0; i<Nx; i++){
		printf("%10d  ", vector[i] );
	}
	printf("\n------------------------------\n\n");
}

void print_matrixi(char *comment, int **matrix, int Nx, int Ny){
	int i, j;
	printf("<:: %s ::>--------------\n", comment);
	for(i=0; i<Nx; i++){
		printf("%3d:", i);
		for(j=0; j<Ny; j++)
			printf("%10d  ", matrix[i][j] );
		printf("\n");
	}
	printf("------------------------------\n\n");
}

void print_gscvectord(char *comment, gsl_complex *vector, int Nx){
	int i;
	printf("<:: %s ::> ( %d )--------------\n", comment, Nx);
	for(i=0; i<Nx; i++){
		printf("%13.9lf, %13.9lf    ", GSL_REAL(vector[i]), GSL_IMAG(vector[i]));
	}
	printf("\n");
	printf("------------------------------\n\n");
}

void print_gscmatrixd_compact(char *comment, gsl_complex **matrix, int Nx, int Ny){
	int i, j;
	printf("<:: %s ::> ( %d x %d )--------------\n", comment, Nx, Ny);
	for(i=0; i<Nx; i++){
		for(j=0; j<Ny; j++)
			printf("%11.8lf, %11.8lf    ", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
		printf("\n");
	}
	printf("------------------------------\n\n");
}

void print_gscmatrixd(char *comment, gsl_complex **matrix, int Nx, int Ny){
	int i, j;
	printf("<:: %s ::> ( %d x %d )--------------\n", comment, Nx, Ny);
	for(i=0; i<Nx; i++){
		for(j=0; j<Ny; j++)
			printf("%13.9lf, %13.9lf    ", GSL_REAL(matrix[i][j]), GSL_IMAG(matrix[i][j]));
		printf("\n");
	}
	printf("------------------------------\n\n");
}

int print_gsctetratensord_nonzero(char *comment, gsl_complex ****matrix, int row, int column, int third, int fourth){
	int i, j, k, l, shown=0;
	printf("<:: %s ::> ( %d x %d x %d x %d )--------------\n", comment, row, column, third, fourth);
	for(i=0; i<row; i++) for(j=0; j<column; j++) for(k=0; k<third; k++) for(l=0; l<fourth; l++){
		if( gsl_complex_abs2(matrix[i][j][k][l]) > 1e-6 ){
			printf("[%2d,%2d,%2d,%2d]\t%10.6lf  %10.6lf\n", i, j, k, l, GSL_REAL(matrix[i][j][k][l]), GSL_IMAG(matrix[i][j][k][l]) );
			shown++;
		}
	}
	printf("------------------------------ total %d\n\n", shown);
	return 0;
}

void print_matrixd(char *comment, double **matrix, int Nx, int Ny){
	int i, j;
	printf("<:: %s ::>--------------\n", comment);
	for(i=0; i<Nx; i++){
		for(j=0; j<Ny; j++)
			printf("%9.6lf\t", matrix[i][j] );
		printf("\n");
	}
	printf("------------------------------\n\n");
}

int print_tetratensord_nonzero(char *comment, double ****matrix, int row, int column, int third, int fourth){
	int i, j, k, l, shown=0;
	printf("<:: %s ::> ( %d x %d x %d x %d )--------------\n", comment, row, column, third, fourth);
	for(i=0; i<row; i++) for(j=0; j<column; j++) for(k=0; k<third; k++) for(l=0; l<fourth; l++){
		if( fabs(matrix[i][j][k][l]) > 1e-6 ){
			printf("[%2d,%2d,%2d,%2d]\t%10.6lf\n", i, j, k, l, matrix[i][j][k][l]);
			shown++;
		}
	}
	printf("------------------------------ total %d\n\n", shown);
	return 0;
}

void inverse_matrixd_lapack(double **input, double **output, int NN){
        int mu, nu;

	integer N = NN;
	integer LWORK = N*N+N;

	double *temp = (double*) malloc ( N*N*sizeof(double) );
	double *WORK = (double*) malloc ( LWORK*sizeof(double) );
	integer *perm = (integer*) malloc ( 2*N*sizeof(integer) );

	integer INFO;
	for(mu=0; mu<N; mu++) for(nu=0; nu<N; nu++){
		temp[mu*N+nu] = input[nu][mu];
	}
	dgetrf_(&N, &N, temp, &N, perm, &INFO );
	if( (int)INFO != 0 ){
		printf("Error in zgetrf: INFO = %d (inverse_matrixd_lapack)\n", (int)INFO);
		print_matrixd("input", input, NN, NN);
		exit(1);
	}
	dgetri_(&N, temp, &N, perm, WORK, &LWORK, &INFO);
	if( (int)INFO != 0 ){
		printf("Error in zgetri: INFO = %d (inverse_matrixd_lapack)\n", (int)INFO);
		print_matrixd("input", input, NN, NN);
		exit(1);
	}
	for(mu=0; mu<N; mu++) for(nu=0; nu<N; nu++){
		output[nu][mu] = temp[mu*N+nu];
	}
	free(temp); free(WORK); free(perm);
}

int similarity(gsl_complex **original, gsl_complex **target, double **transform, int NN){
	int mu, nu, k;
	gsl_complex **backup, **right, zero = gsl_complex_rect(0,0);
	double **inv;

	right = mkgscmatrixd(NN, NN);
	backup = mkgscmatrixd(NN, NN);
	inv = mkmatrixd(NN, NN);

	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
		backup[mu][nu] = original[mu][nu];
		target[mu][nu] = zero;
		right[mu][nu] = zero;
		inv[mu][nu] = transform[nu][mu];
	}
	//inverse_complex_matrix_lapack( transform, inv, NN );

	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++){
		right[mu][nu]
			= gsl_complex_add(
				right[mu][nu],
				gsl_complex_mul_real(backup[mu][k], transform[k][nu])
			);
	}
	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
		target[mu][nu] =
			gsl_complex_add(
				target[mu][nu],
				gsl_complex_mul_real( right[k][nu], inv[mu][k] )
			);
		
	//print_gscmatrixd("original", original, NN, NN);
	//print_gscmatrixd("target", target, NN, NN);

	freegscmatrixd(right, NN);
	freegscmatrixd(backup, NN);
	freematrixd(inv, NN);

	return 0;
}

int unitary(gsl_complex **original, gsl_complex **target, gsl_complex **transform, int NN){//inverse comes left
	int mu, nu, k;
	gsl_complex **backup, **right, zero = gsl_complex_rect(0,0), **inv;

	right = mkgscmatrixd(NN, NN);
	backup = mkgscmatrixd(NN, NN);
	inv = mkgscmatrixd(NN, NN);

	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
		backup[mu][nu] = original[mu][nu];
		target[mu][nu] = zero;
		right[mu][nu] = zero;
		inv[mu][nu] = gsl_complex_conjugate( transform[nu][mu] );
	}
	//inverse_complex_matrix_lapack( transform, inv, NN );

	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++){
		target[mu][nu]
			= gsl_complex_add(
				target[mu][nu],
				gsl_complex_mul(inv[mu][k], transform[k][nu])
			);
		right[mu][nu]
			= gsl_complex_add(
				right[mu][nu],
				gsl_complex_mul(backup[mu][k], transform[k][nu])
			);
	}
	for(mu=0; mu<NN; mu++){
		if( fabs(1-gsl_complex_abs2(target[mu][mu]))> 1e-6 ){
			printf("wrong diagonal part\n");
			exit(1);
		}
		for(nu=mu+1; nu<NN; nu++){
			if( gsl_complex_abs2(target[mu][nu]) > 1e-6 ){
				printf("wrong off-diagonal part\n");
				exit(1);
			}
		}
		for(nu=0; nu<NN; nu++)	target[mu][nu] = zero;
	}
	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
		target[mu][nu] =
			gsl_complex_add(
				target[mu][nu],
				gsl_complex_mul( inv[mu][k], right[k][nu] )
			);
		
	//print_gscmatrixd("original", original, NN, NN);
	//print_gscmatrixd("target", target, NN, NN);

	freegscmatrixd(right, NN);
	freegscmatrixd(backup, NN);
	freegscmatrixd(inv, NN);

	return 0;
}

int check_unitary(gsl_complex **transform, int NN){
	int mu, nu, k;
	gsl_complex **target, zero = gsl_complex_rect(0,0), **inv;

	inv = mkgscmatrixd(NN, NN);
	target = mkgscmatrixd(NN, NN);

	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
		target[mu][nu] = zero;
		inv[mu][nu] = gsl_complex_conjugate( transform[nu][mu] );
	}
	//inverse_complex_matrix_lapack( transform, inv, NN );

	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++){
		target[mu][nu]
			= gsl_complex_add(
				target[mu][nu],
				gsl_complex_mul(inv[mu][k], transform[k][nu])
			);
	}
	for(mu=0; mu<NN; mu++){
		if( fabs(1-gsl_complex_abs2(target[mu][mu]))> 1e-6 ){
			printf("wrong diagonal part\n");
			exit(1);
		}
		for(nu=mu+1; nu<NN; nu++){
			if( gsl_complex_abs2(target[mu][nu]) > 1e-6 ){
				printf("wrong off-diagonal part\n");
				exit(1);
			}
		}
	}
	print_gscmatrixd("check_unitary: target", target, NN, NN);
	freegscmatrixd(inv, NN);
	freegscmatrixd(target, NN);
	return 0;
}

int iunitary(gsl_complex **original, gsl_complex **target, gsl_complex **transform, int NN){
	int mu, nu, k;
	gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv;

	backup = mkgscmatrixd(NN, NN);
	right = mkgscmatrixd(NN, NN);
	pose = mkgscmatrixd(NN, NN);
	inv = mkgscmatrixd(NN, NN);

	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
		pose[mu][nu] = gsl_complex_conjugate(transform[nu][mu]);
		backup[mu][nu] = original[mu][nu];
		target[mu][nu] = zero;
		right[mu][nu] = zero;
		inv[mu][nu] = transform[mu][nu];
	}

	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
		right[mu][nu]
			= gsl_complex_add(
				right[mu][nu],
				gsl_complex_mul(backup[mu][k], pose[k][nu])
			);
	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
		target[mu][nu] =
			gsl_complex_add(
				target[mu][nu],
				gsl_complex_mul( inv[mu][k], right[k][nu] )
			);
		
	//print_gscmatrixd("original", backup, NN, NN);
	//print_gscmatrixd("target", target, NN, NN);

	freegscmatrixd(right, NN);
	freegscmatrixd(pose, NN);
	freegscmatrixd(inv, NN);
	freegscmatrixd(backup, NN);

	return 0;
}

int gsc_mul(gsl_complex **left, gsl_complex **right, gsl_complex **output, int NN){
	int mu, nu, k;
	gsl_complex zero = gsl_complex_rect(0,0);
	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
		output[mu][nu] = zero;
	}
	for(mu=0; mu<NN; mu++) for(k=0; k<NN; k++) for(nu=0; nu<NN; nu++){
		output[mu][nu]
			= gsl_complex_add(
				output[mu][nu],
				gsl_complex_mul( left[mu][k], right[k][nu] )
			);
	}

	return 0;
}

int eigen_gsl(gsl_complex **matrix, double *eval, gsl_complex **evec, int NN){
	int mu, nu;
	gsl_eigen_hermv_workspace *work = gsl_eigen_hermv_alloc( NN );
	gsl_matrix_complex *A = gsl_matrix_complex_alloc(NN, NN);
	gsl_matrix_complex *eve = gsl_matrix_complex_alloc(NN, NN);
	gsl_vector *ev = gsl_vector_alloc(NN);
	for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
		gsl_matrix_complex_set(A, mu, nu, matrix[mu][nu]);
	}

	gsl_eigen_hermv(A, ev, eve, work);
	gsl_eigen_hermv_sort(ev, eve, GSL_EIGEN_SORT_VAL_DESC);

	for(mu=0; mu<NN; mu++){
		eval[mu] = gsl_vector_get(ev, mu);
		for(nu=0; nu<NN; nu++)
			evec[mu][nu] = gsl_matrix_complex_get(eve, nu, mu);	//column is the eigenvector in gsl, so transpose
	}

	gsl_matrix_complex_free(A);
	gsl_matrix_complex_free(eve);
	gsl_vector_free(ev);

	gsl_eigen_hermv_free(work);
	return 0;
}

int eigen_lapack(gsl_complex **matrix, double *eval, gsl_complex **evec, int NN, int ASC){
	int mu, nu;
	char JOBZ = 'V', UPLO='U';
	integer N=NN, LDA = NN, LWORK = -1, LRWORK=-1, INFO;
	doublecomplex *WORK, wkopt;
	doublereal *RWORK, W[NN*NN];

	doublecomplex subspace[NN*NN];
	for(mu=0; mu<NN*NN; mu++){
		subspace[mu].r = 0;
		subspace[mu].i = 0;
	}
	for(mu=0; mu<NN; mu++) for(nu=mu; nu<NN; nu++){
		subspace[nu*NN+mu].r =  GSL_REAL(matrix[mu][nu]);
		subspace[nu*NN+mu].i =  GSL_IMAG(matrix[mu][nu]);
	}
	LRWORK = 3*NN-2;
	RWORK = (doublereal *) malloc ( LRWORK * sizeof(doublereal) );
	zheev_(&JOBZ, &UPLO, &N, subspace, &LDA, W, &wkopt, &LWORK, RWORK, &INFO);
	LWORK = wkopt.r;
	WORK = (doublecomplex *) malloc ( LWORK * sizeof(doublecomplex) );

	zheev_(&JOBZ, &UPLO, &N, subspace, &LDA, W, WORK, &LWORK, RWORK, &INFO);
	if( ASC ) for(mu=0; mu<NN; mu++){
		eval[mu] = W[mu];
		for(nu=0; nu<NN; nu++)
			evec[mu][nu] = gsl_complex_rect(subspace[mu*NN+nu].r, subspace[mu*NN+nu].i);	//evec[0][] = lowest ev
	}
	else for(mu=0; mu<NN; mu++){
		eval[mu] = W[NN-mu-1];
		for(nu=0; nu<NN; nu++)
			GSL_SET_COMPLEX( &evec[mu][nu],
					subspace[(NN-mu-1)*NN+nu].r, subspace[(NN-mu-1)*NN+nu].i);	//evec[0][] = lowest ev
	}

	free(RWORK);
	free(WORK);

	return (int)INFO;
}

