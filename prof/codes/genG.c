#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"

#include "f2c.h"
#include "ran2.c"
#include "clapack.h"

#include "matrix.h"
#include "inverse.c"
#include "basic.h"
#include "hamil.h"
#include "lattice.h"
#include "ls_basis.h"
#include "response.c"
#include "hop.c"

#define Epsilon		0.02

#define Naxis_default	7

#define max_x		1.
#define Plot_range	20.
#define PERIOD		1
#define Nkunit		100

char get_freq_code(gsl_complex *freq, int Nfreq);
double compute_ldos(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory, gsl_complex *freq, int Nfreq);
int generate_green(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory);
int continued_fraction_general(gsl_complex ***Green, int mu, int nu, int direction, double *ai, double *bi, double gm, double p0inner, int dimension, gsl_complex *freq, int Nfreq);
int build_bath_general(int n, gsl_complex ***Green_bath, PNSpara egv, gsl_complex *freq, int Nfreq);
int print_green_lattice_block(char *para, gsl_complex ****Green, gsl_complex *freq, int Nfreq);

void print_complex(gsl_complex cmnumber);

long ranseed = -1;
double Hz, JTD, SOC, J, MU, DELTA, Ts, Uinter, Efermi;
int beta, Nmax, nctot, nb, tnb, myrank=0, size, beta_real, Nonzero_SOC, ROTATE, Blocks, Powns, Powns2, Ndd, Njj, AF=0;
gsl_complex zero, ***Greennew, **Ecluster;
MPI_Datatype MPI_GSL_COMPLEX;

int main(int argc, char *argv[]){
	if(argc != 12){
		printf("Usage :: %s <input.mak> <Hz> <JTD> <SOC> <J> <U_init> <U_final> <DELTA> <save_directory> <beta_real> <AF>\n", argv[0]);
		exit(1);
	}
	int i, j, k, count = 0;
	int degeneracy[Ni], tol;

	PNSpara egv[Ni];
	double U_init, U_final, re, im;
	double ***ginformation;

	FILE *fd; 
	char pathinform[1024], save_directory[1024];

	MPI_Datatype itype[2] = {MPI_DOUBLE, MPI_DOUBLE};
	MPI_Aint idisp[2] = {0, sizeof(double)};
	int iblock[2] = {1, 1};//real and imaginary part, respectively

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Type_create_struct(2, iblock, idisp, itype, &MPI_GSL_COMPLEX);
	MPI_Type_commit( &MPI_GSL_COMPLEX );

	nctot = NC;
	nb = NB;
	tnb = 2*NB;
	Powns	= (int)pow(2,Ns);
	Powns2	= (int)pow(4,Ns);
	zero = gsl_complex_rect(0,0);
	Hz	= atof( argv[2] );
	JTD	= atof( argv[3] );
	SOC	= atof( argv[4] );
	J	= atof( argv[5] );
	U_init	= atof( argv[6] );
	U_final	= atof( argv[7] );
	DELTA	= atof( argv[8] );
	beta_real	= atoi( argv[10] );
	AF	= atoi( argv[11] );

	if( fabs(SOC) < 1e-6 && !SpinFlip ){
		ROTATE = 0;
		Nonzero_SOC = 0;
		Blocks = (Ns+1)*(Ns+1);
	}
	else if( fabs(SOC) < 1e-6 && SpinFlip ) {
		ROTATE = 0;
		Nonzero_SOC = 0;
	}
	else{
		ROTATE = 1;
		Nonzero_SOC = 1;
		Blocks	= 2*Ns+1;
	}

	sprintf(save_directory, "%s", argv[9]);

	fd = fopen(argv[1], "r");
	fscanf(fd, "%d", &count);
	fscanf(fd, "%lf", &Uinter);
	fscanf(fd, "%lf", &MU);
	fscanf(fd, "%d", &tol);
	fscanf(fd, "%d", &beta);
	fscanf(fd, "%d", &Nmax);

	ginformation = (double ***) malloc( Ni*sizeof(double**) );

	int n;
	for(n=0; n<Ni; n++){
		egv[n].egbath = mkvectord(tnb);	
		egv[n].hybrid = mkgscmatrixd(tNC, tnb);
		fscanf(fd, "%d", &degeneracy[n]);

		ginformation[n] = mkmatrixd(degeneracy[n], 3);
		for(i=0; i<degeneracy[n]; i++) for(j=0; j<3; j++)	fscanf(fd, "%lf", &ginformation[n][i][j]);

		for(k=0; k<tnb; k++){
			fscanf(fd, "%lf", &re);
			egv[n].egbath[k] = re;
		}
		for(j=0; j<tNC; j++) for(k=0; k<tnb; k++){
			fscanf(fd, "%lf %lf", &re, &im);
			egv[n].hybrid[j][k] = gsl_complex_rect( re, im );
		}
	}
	fclose(fd);

	printf("U = %.2lf, MU = %.2lf (of UF%.2lf UT%.2lf) from '%s'\n", Uinter, MU, U_init, U_final, argv[1]);

	Ecluster = mkgscmatrixd( NU, NU);
	char paraEcluster[1024];
	sprintf(paraEcluster, "%s/Ecluster.dat", save_directory);
	load_Ecluster(paraEcluster, Ecluster);
	update_chem(Ecluster);

	sprintf(pathinform, "%s/lattice/matsu", save_directory);
	mkdirectory(pathinform);
	generate_green(ginformation, egv, degeneracy, count, save_directory);

	for(n=0; n<Ni; n++)	freematrixd(ginformation[n], degeneracy[n]);
	free(ginformation);
	freegscmatrixd(Ecluster, NU);

	MPI_Finalize();
	return 0;
}

int load_Green_spin(char *para, gsl_complex ******Green, int *degeneracy, double ***ginformation, gsl_complex *freq, int Nfreq){
	FILE *fp;
	double *ai, *bi, p0inner, gm;
	int n, k, mu, nu, direction, find, dimension, degen, spin;
	fp = fopen(para, "r");	nofile(fp, para);

	do{
		if( fscanf(fp, "%d %d %d %d %d %d %d", &n, &degen, &spin, &mu, &nu, &direction, &find) == EOF )	break;
		printf("%d:: [ %d / %d ] spin%d (%d, %d) d_%d (%d)\t", n, degen, degeneracy[n], spin, mu, nu, direction, find);

		gm = ginformation[n][degen][0];

		if( find == 0 ){
			printf("\n");
			continue;
		}
		else{
			fscanf(fp, "%d %lf", &dimension, &p0inner);
			printf("p0inner_complex = %.10lf\n", p0inner);

			ai = mkvectord(dimension);	bi = mkvectord(dimension);

			for(k=0; k<dimension; k++)	fscanf(fp, "%lf", &ai[k]); 
			for(k=0; k<dimension; k++)	fscanf(fp, "%lf", &bi[k]); 

			continued_fraction_general(Green[n][degen][spin], mu, nu, direction, ai, bi, gm, p0inner, dimension, freq, Nfreq);

			free(ai);	free(bi);
		}
	}while(1);
	fclose(fp);

	return 0;
}
int load_Green_nospin(char *para, gsl_complex *****Green, int *degeneracy, double ***ginformation, gsl_complex *freq, int Nfreq){
	FILE *fp;
	double *ai, *bi, p0inner, gm;
	int n, k, mu, nu, direction, find, dimension, degen;
	fp = fopen(para, "r");	nofile(fp, para);

	do{
		if( fscanf(fp, "%d %d %d %d %d %d", &n, &degen, &mu, &nu, &direction, &find) == EOF )	break;
		printf("%d:: [ %d / %d ] (%d, %d) d_%d (%d)\t", n, degen, degeneracy[n], mu, nu, direction, find);

		gm = ginformation[n][degen][0];

		if( find == 0 ){
			printf("\n");
			continue;
		}
		else{
			fscanf(fp, "%d %lf", &dimension, &p0inner);
			printf("p0inner_complex = %.10lf\n", p0inner);

			ai = mkvectord(dimension);	bi = mkvectord(dimension);

			for(k=0; k<dimension; k++)	fscanf(fp, "%lf", &ai[k]); 
			for(k=0; k<dimension; k++)	fscanf(fp, "%lf", &bi[k]); 

			continued_fraction_general(Green[n][degen], mu, nu, direction, ai, bi, gm, p0inner, dimension, freq, Nfreq);

			free(ai);	free(bi);
		}
	}while(1);
	fclose(fp);

	return 0;
}

int compute_self_general(double ***ginformation, PNSpara *egv, char *save_directory, int *degeneracy, int count, gsl_complex ****Self, gsl_complex *freq, int Nfreq){
	int n, mu, nu, i, sp, degen;
	double myweight, partition = 0, myenergy, gm;
	gsl_complex ****Green_bath, ****Ginverse, ****Green_total;
	gsl_complex *****tildekai, *****kai;
	gsl_complex ******Giplus, ******Green_on, **Grplus;

	Giplus		= (gsl_complex ******) malloc ( Ni * sizeof(gsl_complex*****) );
	Green_on	= (gsl_complex ******) malloc ( Ni * sizeof(gsl_complex*****) );
	tildekai	= (gsl_complex *****) malloc ( Ni * sizeof(gsl_complex****) );
	kai		= (gsl_complex *****) malloc ( Ni * sizeof(gsl_complex****) );
	for(n=0; n<Ni; n++){
		Giplus[n]	= mkgscpentatensord(degeneracy[n], Spmax, nctot, nctot, Nfreq);
		Green_on[n]	= mkgscpentatensord(degeneracy[n], Spmax, nctot, nctot, Nfreq);
		tildekai[n]	= mkgsctetratensord(degeneracy[n], nctot, nctot, Nfreq);
		kai[n]		= mkgsctetratensord(degeneracy[n], nctot, nctot, Nfreq);
	}

	Ginverse	= mkgsctetratensord(Ni, Nfreq, tNC, tNC);
	Green_total	= mkgsctetratensord(Ni, Nfreq, tNC, tNC);

	Green_bath	= mkgsctetratensord(Ni, Nfreq, tNC, tNC);
	Grplus		= mkgscmatrixd(Spmax, Nfreq);

	char paragreenon[1024], paragiplus[1024], paratildekai[1024], parakai[1024];

	sprintf(paragreenon, "%s/datafile/greenon_u%.2lf_%dth.ab", save_directory, Uinter, count);
	sprintf(paragiplus, "%s/datafile/giplus_u%.2lf_%dth.ab", save_directory, Uinter, count);
	sprintf(paratildekai, "%s/datafile/tildekai_u%.2lf_%dth.ab", save_directory, Uinter, count);
	sprintf(parakai, "%s/datafile/kai_u%.2lf_%dth.ab", save_directory, Uinter, count);

	for(n=0; n<Ni; n++){
		for(degen=0; degen<degeneracy[n]; degen++) for(sp=0; sp<Spmax; sp++) for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nfreq; i++){
			Green_on[n][degen][sp][mu][nu][i] = zero;
			Giplus[n][degen][sp][mu][nu][i] = zero;
		}
		for(degen=0; degen<degeneracy[n]; degen++) for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nfreq; i++){
			tildekai[n][degen][mu][nu][i] = zero;
			kai[n][degen][mu][nu][i] = zero;
		}
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) for(i=0; i<Nfreq; i++){
			Green_total[n][i][mu][nu] = zero;
		}
	}
	if( nctot > 1 ) load_Green_spin(paragiplus, Giplus, degeneracy, ginformation, freq, Nfreq);
	load_Green_spin(paragreenon, Green_on, degeneracy, ginformation, freq, Nfreq);
	if( Nonzero_SOC ){
		load_Green_nospin(paratildekai, tildekai, degeneracy, ginformation, freq, Nfreq);
		load_Green_nospin(parakai, kai, degeneracy, ginformation, freq, Nfreq);
	}

	gsl_complex	onemi, onepi;
	GSL_SET_COMPLEX(&onemi, 1, -1);
	GSL_SET_COMPLEX(&onepi, 1, 1);

	for(n=0; n<Ni; n++){
		gm = ginformation[n][0][0];
		partition = 0;
		for(degen=0; degen<degeneracy[n]; degen++){
			myenergy = ginformation[n][degen][0];
			if( beta_real )	myweight	= exp( - ( myenergy - gm )*beta_real );
			else		myweight	= 1;
			partition += myweight;
			for(sp=0; sp<Spmax; sp++) {
				for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nfreq; i++){
					Giplus[n][degen][sp][mu][nu][i] = gsl_complex_mul( Giplus[n][degen][sp][mu][nu][i], gsl_complex_rect(0, -1 ) );
				}
				for(mu=0; mu<nctot; mu++) for(nu=mu+1; nu<nctot; nu++) for(i=0; i<Nfreq; i++){
					Grplus[sp][i] = Green_on[n][degen][sp][mu][nu][i];
					Green_on[n][degen][sp][mu][nu][i]
						= gsl_complex_sub(
								Green_on[n][degen][sp][mu][nu][i],
								gsl_complex_add(
									gsl_complex_mul( Green_on[n][degen][sp][mu][mu][i], onemi ),
									gsl_complex_mul( Green_on[n][degen][sp][nu][nu][i], onemi )
									)
								);
					Green_on[n][degen][sp][mu][nu][i] = gsl_complex_add( Green_on[n][degen][sp][mu][nu][i], Giplus[n][degen][sp][mu][nu][i] );
					Green_on[n][degen][sp][mu][nu][i] = gsl_complex_div_real( Green_on[n][degen][sp][mu][nu][i], 2 );

					Green_on[n][degen][sp][nu][mu][i]
						= gsl_complex_sub(
								Grplus[sp][i],
								gsl_complex_add(
									gsl_complex_add( Green_on[n][degen][sp][mu][mu][i], Green_on[n][degen][sp][nu][nu][i] ),
									Green_on[n][degen][sp][mu][nu][i]
									)
								);
				}
				for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nfreq; i++){
					Green_total[n][i][2*mu+sp][2*nu+sp] = gsl_complex_add( Green_total[n][i][2*mu+sp][2*nu+sp], gsl_complex_mul_real( Green_on[n][degen][sp][mu][nu][i], myweight) ) ;
				}
			}
			if( Nonzero_SOC ) for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nfreq; i++){
				Green_total[n][i][2*mu][2*nu+1]
					= gsl_complex_add(
							Green_total[n][i][2*mu][2*nu+1],
							gsl_complex_mul_real(
								gsl_complex_sub(
									gsl_complex_add(
										kai[n][degen][mu][nu][i],
										gsl_complex_mul(
											tildekai[n][degen][mu][nu][i],
											gsl_complex_rect(0,-1)
											)
										),
									gsl_complex_mul(
										gsl_complex_add(
											Green_on[n][degen][0][mu][mu][i],
											Green_on[n][degen][1][nu][nu][i]
											),
										onemi
										)
									),
								myweight/2.
								)
								);
				Green_total[n][i][2*nu+1][2*mu]
					= gsl_complex_add(
							Green_total[n][i][2*nu+1][2*mu],
							gsl_complex_mul_real(
								gsl_complex_sub(
									gsl_complex_add(
										kai[n][degen][mu][nu][i],
										gsl_complex_mul(
											tildekai[n][degen][mu][nu][i],
											gsl_complex_rect(0,1)
											)
										),
									gsl_complex_mul(
										gsl_complex_add(
											Green_on[n][degen][0][mu][mu][i],
											Green_on[n][degen][1][nu][nu][i]
											),
										onepi
										)
									),
								myweight/2.
								)
								);
			}
		}
		for(i=0; i<Nfreq; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
			Green_total[n][i][mu][nu] = gsl_complex_div_real( Green_total[n][i][mu][nu], partition);
		}
		build_bath_general(n, Green_bath[n], egv[n], freq, Nfreq);
		inverse_complex_matrix_all_omp( Green_total[n], Ginverse[n], 0, Nfreq, tNC);
		for(i=0; i<Nfreq; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
			Self[n][i][mu][nu] = gsl_complex_sub( Green_bath[n][i][mu][nu], Ginverse[n][i][mu][nu] );
		}
	}

	char code = get_freq_code(freq, Nfreq);
	char paratemp[1024], parari[1024];
	if( code == 'r' )	sprintf(parari, "real");
	else if( code == 'i' )	sprintf(parari, "imag");
	else			my_error("code error");
	sprintf(paratemp, "%s/lattice/matsu/green%s_Nplot%d_ep%.3lf_%dth.dat", save_directory, parari, Nfreq, Epsilon, count);
	print_green_lattice_block(paratemp, Green_total, freq, Nfreq);

	sprintf(paratemp, "%s/lattice/matsu/self%s_Nplot%d_ep%.3lf_%dth.dat", save_directory, parari, Nfreq, Epsilon, count);
	print_green_lattice_block(paratemp, Self, freq, Nfreq);
	gsl_complex ****Self_t2g = mkgsctetratensord(Ni, Nfreq, tNC, tNC);


	gsl_complex **transform = mkgscmatrixd(NU, NU);
	if( ROTATE ){
		init_transform(transform);
		for(n=0; n<Ni; n++) for(i=0; i<Nfreq; i++)
			unitary(Self[n][i], Self_t2g[n][i], transform, tNC);
	}
	else{
		for(n=0; n<Ni; n++) for(i=0; i<Nfreq; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)
			Self_t2g[n][i][mu][nu] = Self[n][i][mu][nu];
	}
	sprintf(paratemp, "%s/lattice/matsu/tself%s_Nplot%d_ep%.3lf_%dth.dat", save_directory, parari, Nfreq, Epsilon, count);
	print_green_lattice_block(paratemp, Self_t2g, freq, Nfreq);
	freegscmatrixd(transform, NU);

	freegsctetratensord(Self_t2g, Ni, Nfreq, tNC);

	if( code == 'r' ){
		sprintf(paratemp, "%s/lattice/matsu/dos_Nplot%d_ep%.3lf_%dth.dat", save_directory, Nfreq, Epsilon, count);
		FILE *ftemp = fopen(paratemp, "w");
		double local[Ni], sum, accu=0, dw = GSL_REAL(freq[1]) - GSL_REAL(freq[0]);
		for(i=0; i<Nfreq; i++){
			sum = 0;
			for(n=0; n<Ni; n++){
				local[n] = 0;
				for(mu=0; mu<tNC; mu++){
					local[n] += -GSL_IMAG(Green_total[n][i][mu][mu])/PI;
				}
				sum += local[n];
			}
			accu += sum;
			fprintf(ftemp, "%lf\t%lf\t\t", GSL_REAL(freq[i]), sum);
			for(n=0; n<Ni; n++)
				fprintf(ftemp, "%lf\t", local[n]); fprintf(ftemp, "\t");
			fprintf(ftemp, "%lf\t\t", accu*dw);
			for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) fprintf(ftemp, "%lf\t", -GSL_IMAG(Green_total[n][i][mu][mu])/PI);
			fprintf(ftemp, "\n");
		}
		fclose(ftemp);
	}

	for(n=0; n<Ni; n++){
		freegscpentatensord(Giplus[n],	degeneracy[n], Spmax, nctot, nctot);
		freegscpentatensord(Green_on[n],	degeneracy[n], Spmax, nctot, nctot);
		freegsctetratensord(tildekai[n],	degeneracy[n], nctot, nctot);
		freegsctetratensord(kai[n],		degeneracy[n], nctot, nctot);
	}
	free(Giplus); free(Green_on); free(tildekai); free(kai);
	freegsctetratensord(Green_bath, Ni, Nfreq, tNC);
	freegsctetratensord(Green_total, Ni, Nfreq, tNC);
	freegsctetratensord(Ginverse, Ni, Nfreq, tNC);

	freegscmatrixd(Grplus, Spmax);

	return 0;
}

int compute_Glocalk_general(gsl_complex **Glocalk, gsl_complex **Glocalkinverse, gsl_complex ****Self, gsl_complex **hopmatrix, double *diag, int i, gsl_complex omega){
	int n, mu, nu, same;
	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
		same = (int)(mu==nu);
		Glocalkinverse[mu][nu]
			= gsl_complex_sub(
					gsl_complex_mul_real( gsl_complex_add_real( omega,  -diag[mu] ), same ),
					hopmatrix[mu][nu]
			);
	}
	for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) 
		Glocalkinverse[tNC*n+mu][tNC*n+nu]
			= gsl_complex_rect(
				 GSL_REAL( Glocalkinverse[tNC*n+mu][tNC*n+nu] ) - GSL_REAL(Self[n][i][mu][nu]) ,
				 GSL_IMAG( Glocalkinverse[tNC*n+mu][tNC*n+nu] ) - GSL_IMAG(Self[n][i][mu][nu]) 
				);
	inverse_complex_matrix_lapack( Glocalkinverse, Glocalk, NU);
	return 0;
}

int read_ground(int n, gsl_complex *ground, int offset, int gndblock, int count, char *save_directory){
	FILE *fgnd;
	char paragnd[1024];
	sprintf(paragnd, "%s/%c%.2lf_%dth_n%d.gnd", save_directory, control_char, control_para, count, n);
	fgnd = fopen(paragnd, "r");
	fseek(fgnd, offset*sizeof(gsl_complex), SEEK_SET );
	fread(ground, sizeof(gsl_complex), gndblock, fgnd);

	fclose(fgnd);

	return 0;
}

double **read_pos(){
	double **pos = mkmatrixd(Ni, 3);
	FILE *fp = fopen("inputs/AVECTORS", "r");
	int n, i, j;
	double dummy;

	for(i=0; i<3; i++) for(j=0; j<3; j++)	fscanf(fp, "%lf", &dummy);

	fscanf(fp, "%d", &n);
	if( n != Ni ){
		printf("In 'AVECTORS', n=%d should be the same with Ni=%d!\n", n, Ni);
		exit(1);
	}

	for(n=0; n<Ni; n++) for(j=0; j<3; j++){
		fscanf(fp, "%lf", &pos[n][j]);
	}
	print_matrixd("pos", pos, Ni, 3);

	fclose(fp);
	return pos;
}

int periodize_green(gsl_complex **Glocal, gsl_complex **Gpd, double kx, double ky, double kz, double **pos){
	int i, j, mu, nu;
	gsl_complex fac;

	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
		Gpd[mu][nu] = zero;
		for(i=0; i<Ni; i++) for(j=0; j<Ni; j++){
			fac = gsl_complex_rect(
				cos( ( pos[j][0] - pos[i][0] ) * kx + ( pos[j][1] - pos[i][1] ) * ky + ( pos[j][2] - pos[i][2] ) * kz),
				sin( ( pos[j][0] - pos[i][0] ) * kx + ( pos[j][1] - pos[i][1] ) * ky + ( pos[j][2] - pos[i][2] ) * kz)
			);
			Gpd[mu][nu] = gsl_complex_add(
					Gpd[mu][nu], 
					gsl_complex_mul( Glocal[tNC*i+mu][tNC*j+nu], fac )
			);
		}
	}
	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) Gpd[mu][nu] = gsl_complex_div_real( Gpd[mu][nu], Ni );

	return 0;
}

int save_green(char *paraname, gsl_complex *freq, gsl_complex ***Green, int Ng, char type, int Npmax, int saveall){
	int i, mu;
	gsl_complex dxG;
	double omega;
	FILE *fk = fopen(paraname, "w");
	for(i=0; i<Npmax; i++){
		dxG = zero;
		for(mu=0; mu<Ng; mu++) 	dxG = gsl_complex_add(dxG, Green[i][mu][mu]);
		omega = type == 'r' ? GSL_REAL(freq[i]) : GSL_IMAG(freq[i]);
		if( saveall < 0 )	fprintf(fk, "%.15lf\t%.15lf\t%.15lf\t\t", omega, GSL_IMAG(dxG)/Ng, GSL_REAL(dxG)/Ng );
		else			fprintf(fk, "%.15lf\t%.15lf\t%.15lf\t\t", omega, GSL_REAL(dxG)/Ng, GSL_IMAG(dxG)/Ng );

		if( abs(saveall) > 1 ) for(mu=0; mu<tNC; mu++)	fprintf(fk, "%.15lf\t%.15lf\t", GSL_REAL(Green[i][mu][mu]), GSL_IMAG(Green[i][mu][mu]) );
		fprintf(fk, "\n");
	}
	fclose(fk);
	return 0;
}

int generate_green(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory){
	double **pos = read_pos();
	double stemp[Naxis_default][3] = {
		{0.5000,  0.5000,  0.500},
		{0.0000,  0.0000,  0.000},
		{0.5000,  0.0000,  0.500},
		{0.5000,  0.2500,  0.750},
		{0.5000,  0.5000,  0.500},
		{0.3750,  0.7500,  0.375},
		{0.0000,  0.0000,  0.000},
	};
	double bvec_default[3][3] = {
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1}
	};
	int Naxis = Naxis_default;
	int n, mu, nu, i, k, axis;
	double kx, ky, kz;

	char indices_default[Naxis_default] = {'L', 'G', 'X', 'W', 'L', 'K','G'};
	gsl_complex ***Glocalkinverse, ***Glocalk, ***Gt2g, ****Self;

	double **sympoint = mkmatrixd(Naxis, 3), **bvec = mkmatrixd( 3, 3 );
	char **indices = mkmatrixc( Naxis, 128 );

	FILE *fb = fopen("inputs/HISYM", "r");
	if( fb != NULL ){
		fscanf(fb, "%d", &Naxis);
		freematrixc(indices, Naxis);
		freematrixd(sympoint, Naxis_default);
		sympoint = mkmatrixd(Naxis, 3);
		indices = mkmatrixc(Naxis, 128);
		for(i=0; i<Naxis; i++){
			for(k=0; k<3; k++){
				fscanf(fb, "%lf", &sympoint[i][k]);
				if( PERIOD )	sympoint[i][k] *= 4*PI;
				else		sympoint[i][k] *= 2*PI;
			}
			fscanf(fb, "%s", indices[i]);
		}
		fclose(fb);
	}
	else{
		for(i=0; i<Naxis; i++) for(k=0; k<3; k++)	sympoint[i][k] = stemp[i][k]*2*PI;
		for(i=0; i<Naxis; i++) indices[i][0] = indices_default[i];
	}

	fb = fopen("inputs/BVECTORS", "r");
	if( fb != NULL ){
		for(i=0; i<3; i++) for(k=0; k<3; k++)	fscanf(fb, "%lf", &bvec[i][k]);
		fclose(fb);
	}
	else{
		for(i=0; i<3; i++) for(k=0; k<3; k++)	bvec[i][k] = bvec_default[i][k];
	}
	printf("Naxis = %d\n", Naxis);
	print_matrixd("sympoint", sympoint, Naxis, 3);
	print_matrixd("bvec", bvec, 3, 3);

	int Nkdx[Naxis-1];
	double dis=0, displace[3] = {0}, dis_tot=0;
	for(i=0; i<Naxis-1; i++){
		kx = sympoint[i+1][0] - sympoint[i][0];
		ky = sympoint[i+1][1] - sympoint[i][1];
		kz = sympoint[i+1][2] - sympoint[i][2];
		for(k=0; k<3; k++)	displace[k] = kx*bvec[0][k] + ky*bvec[1][k] + kz*bvec[2][k];
		print_vectord("displace", displace, 3);
		dis_tot += sqrt(displace[0]*displace[0] + displace[1]*displace[1] + displace[2]*displace[2]);
		printf("%lf\n", dis_tot);	fflush(stdout);
	}
	for(i=0; i<Naxis-1; i++){
		kx = sympoint[i+1][0] - sympoint[i][0];
		ky = sympoint[i+1][1] - sympoint[i][1];
		kz = sympoint[i+1][2] - sympoint[i][2];
		for(k=0; k<3; k++)	displace[k] = kx*bvec[0][k] + ky*bvec[1][k] + kz*bvec[2][k];
		dis = sqrt(displace[0]*displace[0] + displace[1]*displace[1] + displace[2]*displace[2]);
		Nkdx[i] = (int)(dis/dis_tot * (Naxis-1) * Nkunit);
	}
	print_vectori("Nkdx", Nkdx, Naxis-1);

	int Nk_total=0, Nkrefer[Naxis];
	double dk[3];
	char parak[1024];

	Nkrefer[0] = 0;
	for(axis=0; axis<Naxis-1; axis++){
		Nk_total += Nkdx[axis];
		Nkrefer[axis+1] = Nk_total;
	}
	gsl_complex ***hopmatrix, *oneline;
	hopmatrix = mkgsctritensord( Nk_total+1, NU, NU );
	oneline = mkgscvectord( (Nk_total+1)*NU*NU );

	int line;
	Latt *data;
	Quad quadrature;
	data = read_lattice(WANNIER, &line, &Efermi);
	init_quadrature(&quadrature);
	char parahop[1024];
	FILE *fhop;
	if( PERIOD )
		sprintf(parahop, "%s/mdisp/post_pd_naxis%d_Nk%d_", Path_hop, Naxis, Nk_total+1);
	else
		sprintf(parahop, "%s/mdisp/post_naxis%d_Nk%d_", Path_hop, Naxis, Nk_total+1);
	for(i=0; i<Naxis; i++)	sprintf(parahop, "%s%s", parahop, indices[i]);
	sprintf(parahop, "%s.hop", parahop);
	fhop = fopen(parahop, "rb");
	if( fhop == NULL ){
		i=0;
		for(axis=0; axis<Naxis-1; axis++) {
			for(k=0; k<3; k++)	dk[k] = sympoint[axis+1][k]-sympoint[axis][k];
			for(k=0; k<Nkdx[axis]; k++){
				kx = sympoint[axis][0] + k*dk[0]/Nkdx[axis];
				ky = sympoint[axis][1] + k*dk[1]/Nkdx[axis];
				kz = sympoint[axis][2] + k*dk[2]/Nkdx[axis];
				compute_hopmatrix(data, hopmatrix[i], kx, ky, kz, line);
				i++;
			}
		}
		axis--;
		kx = sympoint[axis][0] + k*dk[0]/Nkdx[axis];
		ky = sympoint[axis][1] + k*dk[1]/Nkdx[axis];
		kz = sympoint[axis][2] + k*dk[2]/Nkdx[axis];
		compute_hopmatrix(data, hopmatrix[i], kx, ky, kz, line);

		i=0;
		for(k=0; k<Nk_total+1; k++){
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				oneline[i] = hopmatrix[k][mu][nu];
				i++;
			}
		}
		fhop = fopen(parahop, "wb");
		fwrite(oneline, sizeof(gsl_complex), (Nk_total+1)*NU*NU, fhop);
		fclose(fhop);
	}
	else{
		fhop = fopen(parahop, "rb");
		fread(oneline, sizeof(gsl_complex), (Nk_total+1)*NU*NU, fhop);
		fclose(fhop);

		i=0;
		for(k=0; k<Nk_total+1; k++){
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				hopmatrix[k][mu][nu] = oneline[i];
				i++;
			}
		}
	}
	free(oneline);

	double distortion[Ni][tNC];
	init_distortion(distortion);
	for(k=0; k<Nk_total+1; k++) for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++)
		GSL_REAL(hopmatrix[k][tNC*n+mu][tNC*n+mu]) = GSL_REAL(hopmatrix[k][tNC*n+mu][tNC*n+mu]) + distortion[n][mu];

	FILE *fp;
	char fname[1024];
	sprintf(fname, "%s/LinG.txt", save_directory);
	fp = fopen(fname, "r");

	gsl_complex **rot;
	if(fp != NULL){
		rot = mkgscmatrixd(NU, NU);
		build_local_to_global(rot, fname);
		for(i=0; i<Nk_total; i++)	iunitary(hopmatrix[i], hopmatrix[i], rot, NU);

		freegscmatrixd(rot, NU);
		fclose(fp);
	}
	gsl_complex **transform = mkgscmatrixd(NU, NU);
	if( ROTATE ){
		printf("performing ROTATION\n");
		init_transform(transform);
		for(i=0; i<Nk_total+1; i++)	iunitary(hopmatrix[i], hopmatrix[i], transform, NU);
	}

	int Nplot = 512;

	Glocalkinverse	= mkgsctritensord( Nplot, NU, NU );
	Glocalk		= mkgsctritensord( Nplot, NU, NU );
	Gt2g		= mkgsctritensord( Nplot, NU, NU );

	double lower=-Plot_range/2., upper=Plot_range/2.;
	gsl_complex **freq = mkgscmatrixd( 2, Nplot );
	for(i=0; i<Nplot; i++){
		freq[0][i] = gsl_complex_rect( 0, (2*i+1)*PI/beta );
		freq[1][i] = gsl_complex_rect( lower +  i*(upper-lower)/(Nplot-1), Epsilon );
	}

	Self		= mkgsctetratensord(Ni, Nplot, tNC, tNC);

	double diag[NU];
	update_diag(diag);
	print_vectord("diag", diag, NU);
	gsl_complex ***Gpd = mkgsctritensord(Nplot, tNC, tNC), ***Gsave;
	char type;
	int myNp, myNg;

	for(int code=0; code<2; code++){
		type = get_freq_code(freq[code], Nplot);
		compute_self_general(ginformation, egv, save_directory, degeneracy, count, Self, freq[code], Nplot);
		for(axis=0; axis<Naxis; axis++){
#pragma omp parallel for default(none) private(n, mu,nu, i) shared(hopmatrix,Glocalk,Glocalkinverse,Self,Nkrefer,axis, diag, ROTATE, transform, Gt2g, Nplot, code, freq)
			for(i=0; i<Nplot; i++){
				compute_Glocalk_general(Glocalk[i], Glocalkinverse[i], Self, hopmatrix[Nkrefer[axis]], diag, i, freq[code][i]);
				if( ROTATE ) unitary(Glocalk[i], Gt2g[i], transform, NU);
			}
			kx = sympoint[axis][0];
			ky = sympoint[axis][1];
			kz = sympoint[axis][2];

			if( PERIOD ){
				for(i=0; i<Nplot; i++)
					periodize_green(Glocalk[i], Gpd[i], kx, ky, kz, pos);
				Gsave = Gpd;
				myNg = tNC;
			}
			else{
				Gsave = Glocalk;
				myNg = NU;
			}
			for(int dev=1; dev<7; dev++){
				myNp = dev*Nplot/6;
				sprintf(parak, "%s/lattice/matsu/%c_AF%d_k%s%d_u%.2lf_ts%.2lf_d%.3lf_ep%.2lf_%dth.dat", save_directory, type, AF, indices[axis], dev, Uinter, Ts, DELTA, Epsilon, count);
				save_green(parak, freq[code], Gsave, myNg, type, myNp, -code - code*dev*(int)(dev==6) );
			}
			if( ROTATE ){
				if( PERIOD ){
					for(i=0; i<Nplot; i++)
						periodize_green(Gt2g[i], Gpd[i], kx, ky, kz, pos);
				}
				else	Gsave = Gt2g;

				for(int dev=1; dev<7; dev++){
					myNp = dev*Nplot/6;
					sprintf(parak, "%s/lattice/matsu/%c_AF%d_tk%s%d_u%.2lf_ts%.2lf_d%.3lf_ep%.2lf_%dth.dat", save_directory, type, AF, indices[axis], dev, Uinter, Ts, DELTA, Epsilon, count);
					save_green(parak, freq[code], Gsave, myNg, type, myNp, -code - code*dev* (int)(dev==6) );
				}
			}
		}
	}
	freegsctritensord(Gpd, Nplot, tNC);

	gsl_complex *myhop, **khop = mkgscmatrixd(NU, NU);
	myhop = mkgscvectord(Nintx*Ninty*Nintz*NU*NU);
	compute_hopmatrix_all(myhop, &quadrature);
	gsl_complex ***Glocal = mkgsctritensord( Nplot, NU, NU );

	double myfactor;
	for(int code=0; code<2; code++){
		for(i=0; i<Nplot; i++) for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++) Glocal[i][mu][nu] = zero;
		type = get_freq_code(freq[code], Nplot);
		compute_self_general(ginformation, egv, save_directory, degeneracy, count, Self, freq[code], Nplot);
		for(int kx=0; kx<Nintx; kx++) for(int ky=0; ky<Ninty; ky++) for(int kz=0; kz<Nintz; kz++) {
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++) khop[mu][nu] = Evaluate(myhop, kx, ky, kz, mu, nu);
			myfactor = (quadrature.weight[0])[kx]*(quadrature.weight[1])[ky]*(quadrature.weight[2])[kz]/Area;
#pragma omp parallel for default(none) private(n, mu,nu, i) shared(Glocalk,Glocalkinverse,Self, diag, Nplot, code, freq, khop, Glocal, myfactor)
			for(i=0; i<Nplot; i++){
				compute_Glocalk_general(Glocalk[i], Glocalkinverse[i], Self, khop, diag, i, freq[code][i]);
				for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
					Glocal[i][mu][nu] = gsl_complex_add( Glocal[i][mu][nu], gsl_complex_mul_real( Glocalk[i][mu][mu], myfactor ) );
			}
		}

		for(int dev=1; dev<7; dev++){
			myNp = dev*Nplot/6;
			sprintf(parak, "%s/lattice/matsu/%c_AF%d_LDOS%d_u%.2lf_ts%.2lf_d%.3lf_ep%.2lf_%dth.dat", save_directory, type, AF, dev, Uinter, Ts, DELTA, Epsilon, count);
			save_green(parak, freq[code], Glocal, NU, type, myNp, -code - code*dev*(int)(dev==6) );
		}
	}

	free(data);
	freematrixd(sympoint, Naxis);
	freematrixd(bvec, 3);
	freematrixc(indices, Naxis);

	freegscmatrixd(freq, 2);
	freegsctritensord(Glocalkinverse, Nplot, NU);
	freegsctritensord(Glocalk, Nplot, NU);
	freegsctritensord(Glocal, Nplot, NU);

	freegscmatrixd(khop, NU);
	freegsctetratensord(Self, Ni, Nplot, tNC);
	freegsctritensord(hopmatrix, Nk_total+1, NU);
	freematrixd(pos, Ni);

	return 0;
}

void print_complex(gsl_complex cmnumber){
	printf("%lf + %lf i\n", GSL_REAL(cmnumber), GSL_IMAG(cmnumber) );
}

int continued_fraction_general(gsl_complex ***Green, int mu, int nu, int direction, double *ai, double *bi, double gm, double p0inner, int dimension, gsl_complex *freq, int Nfreq){
	int i, k;
	gsl_complex tempg;
	for(i=0; i<Nfreq; i++) {
		tempg = zero;
		for(k=dimension-1; k>0; k--){
			tempg = gsl_complex_div(
					gsl_complex_rect(bi[k], 0),
					gsl_complex_sub(
						gsl_complex_add_real( freq[i], -1*ai[k] * direction + gm * direction),
						tempg
					)
				);
		}
		tempg = gsl_complex_sub(
				gsl_complex_add_real( freq[i], -1*ai[0] * direction + gm * direction),
				tempg
			);
		Green[mu][nu][i]
			= gsl_complex_add(
					Green[mu][nu][i],
					gsl_complex_mul_real(
						gsl_complex_inverse(tempg),
						p0inner 
					)
			);
	}
	return 0;
}//end of continued_fraction_general

int build_bath_general(int n, gsl_complex ***Green_bath, PNSpara egv, gsl_complex *freq, int Nfreq){
	gsl_complex tempGreen;
	int mu, nu, k, i;
	for(i=0; i<Nfreq; i++) {
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
			tempGreen = zero;
			for(k=0; k<tnb; k++){
				tempGreen = gsl_complex_sub(
						tempGreen,
						gsl_complex_div(
							gsl_complex_mul(
								gsl_complex_conjugate( egv.hybrid[nu][k] ),
								egv.hybrid[mu][k] 
							),
							gsl_complex_sub_real( freq[i], egv.egbath[k] )
						)
					);
					//(-) sign is included.
			}
			Green_bath[i][mu][nu]
				= gsl_complex_add(
						gsl_complex_sub(
							gsl_complex_mul_real(freq[i], (int)(mu==nu) ),
							Ecluster[tNC*n+mu][tNC*n+nu]
						),
						tempGreen
					);
		}
	}
	return 0;
}//end of build_bath_general

char get_freq_code(gsl_complex *freq, int Nfreq){
	char code;
	gsl_complex dw;
	if( Nfreq > 1 )	dw = gsl_complex_sub( freq[1], freq[0] );
	else		dw = zero;
	if( fabs(GSL_REAL( dw ) ) < 1e-10 )	code = 'i';
	else					code = 'r';
	return code;
}

int print_green_lattice_block(char *para, gsl_complex ****Green, gsl_complex *freq, int Nfreq){
	FILE *fp = fopen(para, "w");
	char code = get_freq_code(freq, Nfreq);
	int mu, nu, i, n;

	double *ff = mkvectord(Nfreq);
	if( code == 'r' )	for(i=0; i<Nfreq; i++)	ff[i] = GSL_REAL(freq[i]);
	else if( code == 'i' )	for(i=0; i<Nfreq; i++)	ff[i] = GSL_IMAG(freq[i]);
	else{
		my_error("code error in 'print_self_general'");
	}
	fprintf(fp, "#wn\t\t\t");
	for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	fprintf(fp, "%d_%d%d\t\t\t\t\t\t\t", n, mu, nu );
	fprintf(fp, "\n");
	for(i=0; i<Nfreq; i++){
		fprintf(fp, "%.18lf\t", ff[i]);
		for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	fprintf(fp, "%20.16lf\t%20.16lf\t\t", GSL_REAL( Green[n][i][mu][nu] ), GSL_IMAG( Green[n][i][mu][nu] ) );
		fprintf(fp, "\n");
	}
	free(ff);
	fclose(fp);
	return 0;
}
