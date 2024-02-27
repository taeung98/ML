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

#define Epsilon		0.03

#define Nplot_dense	2048
#define Nplot_sparse	1024

#define Nplot		1625
#define max_x		1.
#define Plot_range	16.

#define delta_plot	(Plot_range/(Nplot-1))
#define wnp		GSL_REAL(freq[i])

char get_freq_code(gsl_complex *freq, int Nfreq);
double compute_ldos(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory, gsl_complex *freq, int Nfreq);
int compute_energy(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory);
int continued_fraction_general(gsl_complex ***Green, int mu, int nu, int direction, double *ai, double *bi, double gm, double p0inner, int dimension, gsl_complex *freq, int Nfreq);
int build_bath_general(int n, gsl_complex ***Green_bath, PNSpara egv, gsl_complex *freq, int Nfreq);
int print_green_lattice_block(char *para, gsl_complex ****Green, gsl_complex *freq, int Nfreq);

void print_complex(gsl_complex cmnumber);

long ranseed = -1;
double Hz, JTD, SOC, J, MU, DELTA, Ts, Uinter;
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

	double lower=-Plot_range/2, upper=Plot_range/2;
	gsl_complex *fre = mkgscvectord( Nplot );
	for(i=0; i<Nplot; i++){
		fre[i] = gsl_complex_rect( lower +  i*(upper-lower)/(Nplot-1), Epsilon );
	}

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

	sprintf(pathinform, "%s/local", save_directory);
	mkdirectory(pathinform);
	sprintf(pathinform, "%s/lattice/%c%.3lf", save_directory, control_char, control_para);
	mkdirectory(pathinform);
	double integrated = compute_ldos(ginformation, egv, degeneracy, count, save_directory, fre, Nplot);
	compute_energy(ginformation, egv, degeneracy, count, save_directory);

	FILE *fp;
	sprintf(pathinform, "%s/local/lp0_Nw%d_Nintx%d_Ninty%d_Nintz%d.dat", save_directory, Nplot, Nintx, Ninty, Nintz);
	fp = fopen(pathinform, "a");
	fprintf(fp, "%.3lf\t%.16lf\n", control_para, integrated);
	fclose(fp);

	free(fre);
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
	sprintf(paratemp, "%s/lattice/%c%.3lf/green%s_Nplot%d_ep%.3lf_%dth.dat", save_directory, control_char, control_para, parari, Nfreq, Epsilon, count);
	print_green_lattice_block(paratemp, Green_total, freq, Nfreq);

	sprintf(paratemp, "%s/lattice/%c%.3lf/self%s_Nplot%d_ep%.3lf_%dth.dat", save_directory, control_char, control_para, parari, Nfreq, Epsilon, count);
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
	sprintf(paratemp, "%s/lattice/%c%.3lf/tself%s_Nplot%d_ep%.3lf_%dth.dat", save_directory, control_char, control_para, parari, Nfreq, Epsilon, count);
	print_green_lattice_block(paratemp, Self_t2g, freq, Nfreq);
	freegscmatrixd(transform, NU);

	freegsctetratensord(Self_t2g, Ni, Nfreq, tNC);

	if( code == 'r' ){
		sprintf(paratemp, "%s/lattice/%c%.3lf/dos_Nplot%d_ep%.3lf_%dth.dat", save_directory, control_char, control_para, Nfreq, Epsilon, count);
		FILE *ftemp = fopen(paratemp, "w");
		double local[Ni], sum, accu=0;
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
			fprintf(ftemp, "%lf\t%lf\t\t", wnp, sum);
			for(n=0; n<Ni; n++)
				fprintf(ftemp, "%lf\t", local[n]); fprintf(ftemp, "\t");
			fprintf(ftemp, "%lf\t\t", accu*delta_plot);
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

int compute_expectations(int n, double **ginformation, int degeneracy, int count, char *save_directory, gsl_complex **first, gsl_complex ****second){
	typebasis **basis;
	basis = mkmatrixb(Powns2, 2);
	sectorinfo *ginfo = (sectorinfo *) malloc( Blocks*sizeof(sectorinfo) );
	int gnd, gndblock, mu, nu, p, k, l, gindex;
	gsl_complex *ground, *p0;

	int offset;
	build_bitbasis(Ns, basis, ginfo, Nonzero_SOC);
	offset = 0;
	//print_ginfo(Ns, ginfo, Blocks);
	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) first[mu][nu] = zero;
	for(int degen=0; degen<degeneracy; degen++){
		gnd = (int)ginformation[degen][1];
		gndblock = (int)ginformation[degen][2];
		ground = mkgscvectord(gndblock);
		p0 = mkgscvectord(gndblock);
		gindex = find_gindex(ginfo, gnd, gndblock, Blocks);

		read_ground(n, ground, offset, gndblock, count, save_directory);
		offset += gndblock;

		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
			if( (mu%2==nu%2) || Nonzero_SOC ){
				for(p=0; p<gndblock; p++)	p0[p] = zero;
				apply_hop(basis, &ginfo[gindex], mu, nu, p0, gsl_complex_rect(1,0), ground);
				for(p=0; p<gndblock; p++)	first[mu][nu] = gsl_complex_add( first[mu][nu], gsl_complex_mul(gsl_complex_conjugate(p0[p]), ground[p] ) );
			}
		}

		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) for(k=0; k<tNC; k++) for(l=0; l<tNC; l++)	second[mu][nu][k][l] = zero;
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) for(k=0; k<tNC; k++) for(l=0; l<tNC; l++){
			if( (mu%2+nu%2 == k%2+l%2) || Nonzero_SOC ){
				for(p=0; p<gndblock; p++)	p0[p] = zero;
				apply_four(basis, &ginfo[gindex], mu, nu, k, l, p0, gsl_complex_rect(1,0), ground);
				for(p=0; p<gndblock; p++)	second[mu][nu][k][l] = gsl_complex_add( second[mu][nu][k][l], gsl_complex_mul(gsl_complex_conjugate(p0[p]), ground[p] ) );
			}
		}
		free(ground);
		free(p0);

	}
	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) first[mu][nu] = gsl_complex_div_real( first[mu][nu], degeneracy );
	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) for(k=0; k<tNC; k++) for(l=0; l<tNC; l++)	second[mu][nu][k][l] = gsl_complex_div_real( second[mu][nu][k][l], degeneracy );

	freematrixb(basis, Powns2);
	return 0;
}

int compute_coefficients(gsl_complex **S0, gsl_complex **S1, gsl_complex **first, gsl_complex ****second){
	double ****U = mktetratensord(tNC, tNC, tNC, tNC);
	double Ktemp, Ltemp;
	int mu, nu, p, q, k, l, i, j;
	int Ndd, Njj;
	Coulomb JJ[tNC*tNC*tNC*tNC];	//In fact, tNC! is enough
	Coulomb DD[tNC*(tNC-1)/2];
	Ndd = init_dd(Uinter, J, DD);
	Njj = init_others(J, JJ);
	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) for(k=0; k<tNC; k++) for(l=0; l<tNC; l++)	U[mu][nu][k][l] = 0;
	for(p=0; p<Ndd; p++){
		printf("%d %d %d %d\t%10.6lf\n", DD[p].i, DD[p].j, DD[p].k, DD[p].l, DD[p].V);
		U[DD[p].i][DD[p].j][DD[p].k][DD[p].l] = DD[p].V/ 4.;
		U[DD[p].j][DD[p].i][DD[p].k][DD[p].l] = DD[p].V/-4.;
		U[DD[p].i][DD[p].j][DD[p].l][DD[p].k] = DD[p].V/-4.;
		U[DD[p].j][DD[p].i][DD[p].l][DD[p].k] = DD[p].V/ 4.;
	}
	for(p=0; p<Njj; p++){
		printf("%d %d %d %d\t%10.6lf\n", JJ[p].i, JJ[p].j, JJ[p].k, JJ[p].l, JJ[p].V);
		U[JJ[p].i][JJ[p].j][JJ[p].k][JJ[p].l] = JJ[p].V/ 4.;
		U[JJ[p].j][JJ[p].i][JJ[p].k][JJ[p].l] = JJ[p].V/-4.;
		U[JJ[p].i][JJ[p].j][JJ[p].l][JJ[p].k] = JJ[p].V/-4.;
		U[JJ[p].j][JJ[p].i][JJ[p].l][JJ[p].k] = JJ[p].V/ 4.;
	}
	print_tetratensord_nonzero("U", U, tNC, tNC, tNC, tNC);
	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
		S0[mu][nu] = zero;
		for(i=0; i<tNC; i++) for(j=0; j<tNC; j++) S0[mu][nu] = gsl_complex_add( S0[mu][nu], gsl_complex_mul_real(first[i][j], 4*U[mu][i][j][nu]) );
	}

	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
		S1[mu][nu] = zero;
		for(i=0; i<tNC; i++) for(j=0; j<tNC; j++) for(k=0; k<tNC; k++) for(l=0; l<tNC; l++){
			Ktemp = 0;
			for(p=0; p<tNC; p++)	Ktemp += ( 4*U[mu][p][k][l] * U[i][j][p][nu]  - 16*U[mu][i][k][p] * U[p][j][l][nu] );
			S1[mu][nu] = gsl_complex_add(S1[mu][nu], gsl_complex_mul_real( second[i][j][k][l], Ktemp ) );
		}
		for(i=0; i<tNC; i++) for(j=0; j<tNC; j++) {
			Ltemp = 0;
			for(p=0; p<tNC; p++) for(q=0; q<tNC; q++)	Ltemp -= 8*U[mu][i][p][q] * U[p][q][j][nu];
			S1[mu][nu] = gsl_complex_add(S1[mu][nu], gsl_complex_mul_real( first[i][j], Ltemp ) );
		}
		for(i=0; i<tNC; i++) S1[mu][nu] = gsl_complex_sub( S1[mu][nu], gsl_complex_mul( S0[mu][i], S0[i][nu] ) );
	}
	freetetratensord(U, tNC, tNC, tNC);

	return 0;
}

int test_self_high_frequencies(double ***ginformation, int *degeneracy, int count, char *save_directory, gsl_complex ****Self, gsl_complex *freq, int Nfreq){
	char parafp[1024];
	int i, mu, nu;
	gsl_complex **first = mkgscmatrixd(tNC, tNC), ****second = mkgsctetratensord(tNC, tNC, tNC, tNC);
	gsl_complex **S0 = mkgscmatrixd(tNC, tNC), **S1 = mkgscmatrixd(tNC, tNC);
	gsl_complex *****Self_high	= mkgscpentatensord(2, Ni, Nfreq, tNC, tNC);

	int n;
	for(n=0; n<Ni; n++){
		compute_expectations(n, ginformation[n], degeneracy[n], count, save_directory, first, second);
		compute_coefficients(S0, S1, first, second);

		sprintf(parafp, "first_%d", n);
		print_gscmatrixd(parafp, first, tNC, tNC);
		sprintf(parafp, "second_%d", n);
		print_gsctetratensord_nonzero(parafp, second, tNC, tNC, tNC, tNC);
		sprintf(parafp, "S0_%d", n);
		print_gscmatrixd(parafp, S0, tNC, tNC);
		sprintf(parafp, "S1_%d", n);
		print_gscmatrixd(parafp, S1, tNC, tNC);

		for(i=0; i<Nfreq; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) Self_high[0][n][i][mu][nu] = S0[mu][nu];
		for(i=0; i<Nfreq; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) Self_high[1][n][i][mu][nu] = gsl_complex_div( S1[mu][nu], freq[i] );

	}
	sprintf(parafp, "%s/lattice/%c%.3lf/highself_Nplot%d_%dth_o0.dat", save_directory, control_char, control_para, Nfreq, count);
	print_green_lattice_block(parafp, Self_high[0], freq, Nfreq);

	sprintf(parafp, "%s/lattice/%c%.3lf/highself_Nplot%d_%dth_o1.dat", save_directory, control_char, control_para, Nfreq, count);
	print_green_lattice_block(parafp, Self_high[1], freq, Nfreq);

	for(n=0; n<Ni; n++) for(i=0; i<Nfreq; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) 
		Self_high[0][n][i][mu][nu] = gsl_complex_sub( gsl_complex_add(Self_high[0][n][i][mu][nu], Self_high[1][n][i][mu][nu]), Self[n][i][mu][nu] );
	sprintf(parafp, "%s/lattice/%c%.3lf/errorself_Nplot%d_%dth_o1.dat", save_directory, control_char, control_para, Nfreq, count);
	print_green_lattice_block(parafp, Self_high[0], freq, Nfreq);

	sprintf(parafp, "%s/lattice/%c%.3lf/selfimag_Nplot%d_%dth_full.dat", save_directory, control_char, control_para, Nfreq, count);
	print_green_lattice_block(parafp, Self, freq, Nfreq);

	for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
		if( gsl_complex_abs( Self_high[0][n][Nfreq-1][mu][nu] ) > 1 ){
			print_gscmatrixd("Self", Self[n][Nfreq-1], tNC, tNC);
			print_gscmatrixd("Self_high", Self_high[1][n][Nfreq-1], tNC, tNC);
			my_error("invalid tail expansion");
		}
		if( gsl_complex_abs( Self_high[0][n][Nfreq-1][mu][nu] ) > 1e-6 ){
			print_gscmatrixd("Self", Self[n][Nfreq-1], tNC, tNC);
			print_gscmatrixd("Self_high", Self_high[1][n][Nfreq-1], tNC, tNC);
			my_error("Perhaps your w_high is not large enough");
		}
	}

	freegscmatrixd(S0, tNC);
	freegscmatrixd(S1, tNC);
	freegscmatrixd(first, tNC);
	freegsctetratensord(second, tNC, tNC, tNC);
	freegscpentatensord(Self_high, 2, Ni, Nfreq, tNC);
	return 0;
}

int compute_energy(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory){
	gsl_complex ***first = mkgsctritensord(Ni, tNC, tNC), *****second = mkgscpentatensord(Ni, tNC, tNC, tNC, tNC);
	gsl_complex ***S0 = mkgsctritensord(Ni, tNC, tNC), ***S1 = mkgsctritensord(Ni, tNC, tNC), **S0_tilde = mkgscmatrixd(NU, NU), **S1_tilde = mkgscmatrixd(NU, NU);

	int n, mu, nu, i, p;

	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++) {
		S0_tilde[mu][nu] = zero;
		S1_tilde[mu][nu] = zero;
	}
	double filling_imp = 0;
	for(n=0; n<Ni; n++){
		compute_expectations(n, ginformation[n], degeneracy[n], count, save_directory, first[n], second[n]);
		compute_coefficients(S0[n], S1[n], first[n], second[n]);
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
			S0_tilde[n*tNC+mu][n*tNC+nu] = S0[n][mu][nu];
			S1_tilde[n*tNC+mu][n*tNC+nu] = S1[n][mu][nu];
		}
		for(mu=0; mu<tNC; mu++)	filling_imp += GSL_REAL(first[n][mu][mu]);
	}

	int kx, ky, kz, Nfreq = Nplot_dense + Nplot_sparse;
	gsl_complex ***Glocalkinverse, ***Glocalk, ****Self, *freq = mkgscvectord( Nfreq );
	double *position, *weight, w_high = 10000;
	position = mkvectord( Nfreq );
	weight = mkvectord( Nfreq );
	gauleg( 0, 20, position-1, weight-1, Nplot_dense);
	gauleg( 20, w_high, position+Nplot_dense-1, weight+Nplot_dense-1, Nplot_sparse);

	for(i=0; i<Nfreq; i++)	freq[i] = gsl_complex_rect(0, position[i]);
	Self		= mkgsctetratensord(Ni, Nfreq, tNC, tNC);

	compute_self_general(ginformation, egv, save_directory, degeneracy, count, Self, freq, Nfreq);
	test_self_high_frequencies(ginformation, degeneracy, count, save_directory, Self, freq, Nfreq);

	Glocalkinverse	= mkgsctritensord( Nfreq, NU, NU );
	Glocalk		= mkgsctritensord( Nfreq, NU, NU );

	int same;

	gsl_complex *hopmatrix, *g2, *g3;
	Quad quadrature;
	hopmatrix = mkgscvectord(Nintx*Ninty*Nintz*NU*NU);
	g2 = mkgscvectord(Nintx*Ninty*Nintz*NU*NU);	initgscvectord(g2, Nintx*Ninty*Nintz*NU*NU);
	g3 = mkgscvectord(Nintx*Ninty*Nintz*NU*NU);	initgscvectord(g3, Nintx*Ninty*Nintz*NU*NU);
	compute_hopmatrix_all(hopmatrix, &quadrature);
	double *g2_ksum, *g3_ksum, **Integrand;
	double H0g1_ksum = 0, H0g2_ksum = 0, H0g3_ksum = 0;
	double Sg1_ksum = 0, Sg2_ksum = 0, Sg3_ksum = 0;
	double *H0Gk = mkvectord( Nfreq );
	double *SGk = mkvectord( Nfreq );

	Integrand	= mkmatrixd( Nfreq, NU );
	g2_ksum		= mkvectord( NU );
	g3_ksum		= mkvectord( NU );

	double diag[NU];
	update_diag(diag);

	for(i=0; i<Nfreq; i++){
		H0Gk[i] = 0;
		SGk[i] = 0;
		for(mu=0; mu<NU; mu++){
			Integrand[i][mu] = 0;
		}
	}
	for(mu=0; mu<NU; mu++){
		H0g1_ksum += GSL_REAL(gsl_complex_sub_real(Ecluster[mu][nu], diag[mu])) /2.;
		Sg1_ksum += GSL_REAL( S0_tilde[mu][mu] )/2.;
		g2_ksum[mu] = 0;
		g3_ksum[mu] = 0;
	}
	double myfactor;
	for(kx=0; kx<Nintx; kx++) for(ky=0; ky<Ninty; ky++) for(kz=0; kz<Nintz; kz++) {
		for(mu=0; mu<NU; mu++) hopmatrix(kx, ky, kz, mu, mu) = gsl_complex_add_real( hopmatrix(kx, ky, kz, mu, mu), diag[mu] );
		for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++) {
			Evaluate(g2, kx, ky, kz, mu, nu) = gsl_complex_add(hopmatrix(kx, ky, kz, mu, nu), S0_tilde[mu][nu] );
			Evaluate(g3, kx, ky, kz, mu, nu) = S1_tilde[mu][nu];
		}

		for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++) for(p=0; p<NU; p++) {
			Evaluate(g3, kx, ky, kz, mu, nu)
				= gsl_complex_add( Evaluate(g3, kx, ky, kz, mu, nu),
					gsl_complex_mul( Evaluate(g2, kx, ky, kz, mu, p), Evaluate(g2, kx, ky, kz, p, nu) )
				);
		}

		myfactor = (quadrature.weight[0])[kx]*(quadrature.weight[1])[ky]*(quadrature.weight[2])[kz]/Area;
		for(mu=0; mu<NU; mu++){
			g2_ksum[mu] -= GSL_REAL(Evaluate(g2, kx, ky, kz, mu, mu)) * myfactor /w_high / PI;			//multiplied by 2/2pi= 1/pi
			g3_ksum[mu] -= GSL_IMAG(Evaluate(g3, kx, ky, kz, mu, mu)) * myfactor /w_high/w_high/ PI / 2.;
		}

		for(mu=0; mu<NU; mu++){
			Sg2_ksum -= GSL_REAL( S1_tilde[mu][mu] ) * myfactor / w_high / PI;
			for(p=0; p<NU; p++) {
				H0g2_ksum -= GSL_REAL( gsl_complex_mul( hopmatrix(kx, ky, kz, mu, p), Evaluate(g2, kx, ky, kz, p, mu) ) ) * myfactor / w_high / PI;
				H0g3_ksum -= GSL_IMAG( gsl_complex_mul( hopmatrix(kx, ky, kz, mu, p), Evaluate(g3, kx, ky, kz, p, mu) ) ) * myfactor / w_high/w_high / PI / 2.;

				Sg2_ksum -= GSL_REAL( gsl_complex_mul( S0_tilde[mu][p], Evaluate(g2, kx, ky, kz, p, mu) ) ) * myfactor / w_high / PI;
				Sg3_ksum -= (
						  GSL_IMAG( gsl_complex_mul( S0_tilde[mu][p], Evaluate(g3, kx, ky, kz, p, mu) ) )
						+ GSL_IMAG( gsl_complex_mul( S1_tilde[mu][p], Evaluate(g2, kx, ky, kz, p, mu) ) )
					) * myfactor / w_high/w_high / PI / 2.;
			}
		}
#pragma omp parallel for default(none) private(n, mu,nu,same, p, i) shared(kx,ky,kz,Integrand,H0Gk,SGk, hopmatrix,diag,Glocalk,Glocalkinverse,Self,freq,Nfreq, quadrature, myfactor, weight)
		for(i=0; i<Nfreq; i++){
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				same = (int)(mu==nu);
				Glocalkinverse[i][mu][nu]
					= gsl_complex_sub(
							gsl_complex_mul_real( freq[i], same ), 
							hopmatrix(kx,ky,kz,mu,nu)
							);
			}
			for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) 
				Glocalkinverse[i][tNC*n+mu][tNC*n+nu]
					= gsl_complex_rect(
							GSL_REAL( Glocalkinverse[i][tNC*n+mu][tNC*n+nu] ) - GSL_REAL(Self[n][i][mu][nu]) ,
							GSL_IMAG( Glocalkinverse[i][tNC*n+mu][tNC*n+nu] ) - GSL_IMAG(Self[n][i][mu][nu]) 
							);
			inverse_complex_matrix_lapack( Glocalkinverse[i], Glocalk[i], NU);

			for(mu=0; mu<NU; mu++){
				Integrand[i][mu] += GSL_REAL(Glocalk[i][mu][mu]) * weight[i] * myfactor /PI;	//multiplied by 2 to cover wn<0
				for(p=0; p<NU; p++){
					H0Gk[i] += GSL_REAL( gsl_complex_mul( hopmatrix(kx, ky, kz, mu, p), Glocalk[i][p][mu] ) ) * weight[i] * myfactor /PI;
				}
			}
			for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(p=0; p<NU; p++) {
				SGk[i] += GSL_REAL( gsl_complex_mul( Self[n][i][mu][p%tNC], Glocalk[i][p][n*tNC+mu] ) ) * weight[i] * myfactor /PI /2.;
			}
		}
	}

	FILE *fp;
	char parafp[1024];

	sprintf(parafp, "%s/lattice/%c%.3lf/Integrand_Nplot%d_ep%.3lf_%dth.dat", save_directory, control_char, control_para, Nfreq, Epsilon, count);
	fp = fopen(parafp, "w");
	fprintf(fp, "#wn\t\t\t");
	for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++)	fprintf(fp, "%d_%d\t\t\t", n, mu);
	fprintf(fp, "\n");
	for(i=0; i<Nfreq; i++){
		fprintf(fp, "%20.16lf\t", GSL_IMAG(freq[i]));
		for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) fprintf(fp, "%20.16lf\t", Integrand[i][n*tNC+mu]);
		fprintf(fp, "\n");
	}
	fclose(fp);


	double sum, occupancy=0, kinetic=0, correlation=0, sum_H0Gk=0, sum_SGk=0;

	for(i=0; i<Nfreq; i++){
		sum_H0Gk += H0Gk[i];
		sum_SGk += SGk[i];
	}
	kinetic		= sum_H0Gk + H0g1_ksum + H0g2_ksum + H0g3_ksum;
	correlation	=  sum_SGk +  Sg1_ksum +  Sg2_ksum +  Sg3_ksum;


	sprintf(parafp, "%s/lattice/filling_Nplot%d_wh%g_%dth.dat", save_directory, Nfreq, w_high, count);
	fp = fopen(parafp, "w");
	fprintf(fp, "#mu\tsum\t\t\tglatt\t\t\tg2\t\t\tg3\n");
	for(mu=0; mu<NU; mu++){
		sum = 0;
		for(i=0; i<Nfreq; i++) sum += Integrand[i][mu];
		fprintf(fp, "%d\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\n", mu, 0.5 + sum + g2_ksum[mu] + g3_ksum[mu], sum, g2_ksum[mu], g3_ksum[mu] );
		occupancy += 0.5 + sum + g2_ksum[mu] + g3_ksum[mu];
	}
	fclose(fp);
	sprintf(parafp, "%s/lattice/occupancy_Nplot%d_wh%g_%dth.dat", save_directory, Nfreq, w_high, count);
	fp = fopen(parafp, "w");
	fprintf(fp, "#lattice\t\t\timpurity\n");
	fprintf(fp, "%20.16lf\t%20.16lf\n", occupancy, filling_imp);
	fclose(fp);


	sprintf(parafp, "%s/lattice/energy_Nplot%d_wh%g_%dth.dat", save_directory, Nfreq, w_high, count);
	fp = fopen(parafp, "w");
	fprintf(fp, "#total\t\t\tchemical\t\t\tkinetic\t\t\tH0Gk\t\t\tH0g1\t\t\tH0g2\t\t\tH0g3\t\t\t");
	fprintf(fp, "correlation\t\t\tSGk\t\t\tSg1\t\t\tSg2\t\t\tSg3\n");

	fprintf(fp, "%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\n",
			kinetic + correlation - MU * occupancy, MU*occupancy, kinetic, sum_H0Gk, H0g1_ksum, H0g2_ksum, H0g3_ksum,
			correlation, sum_SGk, Sg1_ksum, Sg2_ksum, Sg3_ksum );

	fclose(fp);

	print_vectord("g2_ksum", g2_ksum, NU);
	print_vectord("g3_ksum", g3_ksum, NU);
	printf("occupancy = %lf\n", occupancy);

	freematrixd(Integrand, Nfreq);
	free(g2_ksum);
	free(g3_ksum);
	freegsctritensord(Glocalkinverse, Nfreq, NU);
	freegsctritensord(Glocalk, Nfreq, NU);

	free(hopmatrix);
	free(g2); free(g3);
	free(position); free(weight);
	freegsctetratensord(Self, Ni, Nfreq, tNC);

	freegscmatrixd(S0_tilde, NU);
	freegscmatrixd(S1_tilde, NU);
	free(H0Gk);
	free(SGk);
	freegsctritensord(S0, Ni, tNC);
	freegsctritensord(S1, Ni, tNC);
	freegsctritensord(first, Ni, tNC);
	freegscpentatensord(second, Ni, tNC, tNC, tNC);
	free(freq);
	return 0;
}



double compute_ldos(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory, gsl_complex *freq, int Nfreq){
	int n, mu, nu, i;
	int kx, ky, kz;

	gsl_complex ***Glatt, ***Glocalkinverse, ***Glocalk, ****Self;

	Self		= mkgsctetratensord(Ni, Nfreq, tNC, tNC);
	compute_self_general(ginformation, egv, save_directory, degeneracy, count, Self, freq, Nfreq);

	Glocalkinverse	= mkgsctritensord( Nfreq, NU, NU );
	Glocalk		= mkgsctritensord( Nfreq, NU, NU );
	Glatt		= mkgsctritensord( Nfreq, NU, NU );

	int same;

	gsl_complex *hopmatrix;
	Quad quadrature;
	hopmatrix = mkgscvectord(Nintx*Ninty*Nintz*NU*NU);
	compute_hopmatrix_all(hopmatrix, &quadrature);
	double **Ldos = mkmatrixd( Nfreq, NU );

	FILE *fp;
	char fname[1024];
	sprintf(fname, "%s/LinG.txt", save_directory);
	fp = fopen(fname, "r");

	gsl_complex **rot;
	if(fp != NULL){
		rot = mkgscmatrixd(NU, NU);
		build_local_to_global(rot, fname);
		rotate_hopmatrix(hopmatrix, rot);

		freegscmatrixd(rot, NU);
		fclose(fp);
	}
	gsl_complex **transform = mkgscmatrixd(NU, NU);
	if( ROTATE ){
		init_transform(transform);
		rotate_hopmatrix(hopmatrix, transform);
	}

	double diag[NU];
	update_diag(diag);

	for(i=0; i<Nfreq; i++){
		for(mu=0; mu<NU; mu++){
			Ldos[i][mu] = 0;
			for(nu=0; nu<NU; nu++)	Glatt[i][mu][nu] = zero;
		}
	}
	for(kx=0; kx<Nintx; kx++) for(ky=0; ky<Ninty; ky++) for(kz=0; kz<Nintz; kz++) {
#pragma omp parallel for default(none) private(n, mu,nu,same, i) shared(kx,ky,kz,Ldos,hopmatrix,diag,Glocalk,Glocalkinverse,Self,freq,Nfreq, quadrature, Glatt)
		for(i=0; i<Nfreq; i++){
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				same = (int)(mu==nu);
				Glocalkinverse[i][mu][nu]
					= gsl_complex_sub(
							gsl_complex_mul_real( gsl_complex_sub_real( freq[i], diag[mu] ), same ), 
							hopmatrix(kx,ky,kz,mu,nu)
							);
			}
			for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) 
				Glocalkinverse[i][tNC*n+mu][tNC*n+nu]
					= gsl_complex_rect(
							GSL_REAL( Glocalkinverse[i][tNC*n+mu][tNC*n+nu] ) - GSL_REAL(Self[n][i][mu][nu]) ,
							GSL_IMAG( Glocalkinverse[i][tNC*n+mu][tNC*n+nu] ) - GSL_IMAG(Self[n][i][mu][nu]) 
							);
			inverse_complex_matrix_lapack( Glocalkinverse[i], Glocalk[i], NU);
			for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++)
				Ldos[i][n*tNC+mu] -= GSL_IMAG(Glocalk[i][n*tNC+mu][n*tNC+mu])/PI/Area*(quadrature.weight[0])[kx]*(quadrature.weight[1])[ky]*(quadrature.weight[2])[kz];
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
				Glatt[i][mu][nu]
					= gsl_complex_add(
						Glatt[i][mu][nu],
						gsl_complex_mul_real( Glocalk[i][mu][nu], (quadrature.weight[0])[kx]*(quadrature.weight[1])[ky]*(quadrature.weight[2])[kz]/Area)
					);

		}
	}
#pragma omp parallel for default(none) private(n, mu,nu,same, i) shared(Glocalk,Glocalkinverse,Self,freq,Nfreq, Glatt, Ecluster,diag)
	for(i=0; i<Nfreq; i++){
		inverse_complex_matrix_lapack( Glatt[i], Glocalkinverse[i], NU);
		for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)
			Glocalkinverse[i][tNC*n+mu][tNC*n+nu] = gsl_complex_add( Glocalkinverse[i][tNC*n+mu][tNC*n+nu], Self[n][i][mu][nu] );
		for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
			same = (int)(mu==nu);
			Glocalkinverse[i][mu][nu]
				= gsl_complex_sub(
						Glocalkinverse[i][mu][nu],
						gsl_complex_mul_real( gsl_complex_sub_real( freq[i], diag[mu] ), same )
					);
		}
	}

	char parafp[1024];
	sprintf(parafp, "%s/lattice/ldos_%c%.2lf_Nplot%d.dat", save_directory, control_char, control_para, Nfreq);
	fp = fopen(parafp, "w");
	double sum = 0;
	for(i=0; i<Nfreq; i++){
		fprintf(fp, "%lf\t", wnp);
		for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++){
			fprintf(fp, "%.16lf\t", Ldos[i][n*tNC+mu]);
			sum += Ldos[i][n*tNC+mu];
		}
		fprintf(fp, "%lf\n", sum*delta_plot);
	}
	fclose(fp);

	sprintf(parafp, "%s/lattice/delta_%c%.2lf_Nplot%d.dat", save_directory, control_char, control_para, Nfreq);
	fp = fopen(parafp, "w");
	for(i=0; i<Nfreq; i++){
		fprintf(fp, "%lf\t", wnp);
		for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
			fprintf(fp, "%.16lf\t%.16lf\t\t", GSL_REAL(Glocalkinverse[i][mu][nu]), GSL_IMAG(Glocalkinverse[i][mu][nu]) );
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	sprintf(parafp, "%s/lattice/%c%.3lf/ldos_Nplot%d_ep%.3lf_%dth.dat", save_directory, control_char, control_para, Nfreq, Epsilon, count);
	fp = fopen(parafp, "w");
	sum = 0;
	for(i=0; i<Nfreq; i++){
		fprintf(fp, "%lf\t", wnp);
		for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++){
			fprintf(fp, "%.16lf\t", Ldos[i][n*tNC+mu]);
			sum += Ldos[i][n*tNC+mu];
		}
		fprintf(fp, "%lf\n", sum*delta_plot);
	}
	fclose(fp);

	sprintf(parafp, "%s/lattice/%c%.3lf/delta_Nplot%d_ep%.3lf_%dth.dat", save_directory, control_char, control_para, Nfreq, Epsilon, count);
	fp = fopen(parafp, "w");
	for(i=0; i<Nfreq; i++){
		fprintf(fp, "%lf\t", wnp);
		for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
			fprintf(fp, "%.16lf\t%.16lf\t\t", GSL_REAL(Glocalkinverse[i][mu][nu]), GSL_IMAG(Glocalkinverse[i][mu][nu]) );
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	sum = 0;
	for(i=0; i<Nfreq/2; i++){
		for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++)
			sum += Ldos[i][n*tNC+mu];
	}
	for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++)
		sum -= (Ldos[0][n*tNC+mu] + Ldos[Nfreq-1][n*tNC+mu])/2.;


	if(ROTATE){
		gsl_complex **tdiag = mkgscmatrixd(NU, NU);
		gsl_complex ***rself = mkgsctritensord(Nfreq, NU, NU);
		for(i=0; i<Nfreq; i++) for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)	rself[i][mu][nu] = zero;
		for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)	tdiag[mu][nu] = zero;
		for(mu=0; mu<NU; mu++){
			GSL_REAL(tdiag[mu][mu]) = diag[mu];
		}
		unitary(tdiag, tdiag, transform, NU);

#pragma omp parallel for default(none) private(n, mu,nu,same, i) shared(Glocalk,Glocalkinverse,Self,freq,Nfreq, Glatt, Ecluster,diag, transform, tdiag, rself)
		for(i=0; i<Nfreq; i++){
			unitary(Glatt[i], Glatt[i], transform, NU);
			for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	rself[i][tNC*n+mu][tNC*n+nu] = Self[n][i][mu][nu];
			unitary(rself[i], rself[i], transform, NU);
			inverse_complex_matrix_lapack( Glatt[i], Glocalkinverse[i], NU);
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				same = (int)(mu==nu);
				Glocalkinverse[i][mu][nu]
					= gsl_complex_sub(
							gsl_complex_add(Glocalkinverse[i][mu][nu], rself[i][mu][nu]),
							gsl_complex_mul_real( gsl_complex_sub( freq[i], tdiag[mu][mu] ), same )
						);
			}
		}
		freegscmatrixd(tdiag, NU);

		for(i=0; i<Nfreq; i++) for(mu=0; mu<NU; mu++){
			Ldos[i][mu] = -GSL_IMAG(Glatt[i][mu][mu])/PI;
		}

		FILE *fp;
		char parafp[1024];
		sprintf(parafp, "%s/lattice/tldos_%c%.2lf_Nplot%d.dat", save_directory, control_char, control_para, Nfreq);
		fp = fopen(parafp, "w");
		double sum = 0;
		for(i=0; i<Nfreq; i++){
			fprintf(fp, "%lf\t", wnp);
			for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++){
				fprintf(fp, "%.16lf\t", Ldos[i][n*tNC+mu]);
				sum += Ldos[i][n*tNC+mu];
			}
			fprintf(fp, "%lf\n", sum*delta_plot);
		}
		fclose(fp);

		sprintf(parafp, "%s/lattice/tselfreal_%c%.2lf_Nplot%d.dat", save_directory, control_char, control_para, Nfreq);
		fp = fopen(parafp, "w");
		for(i=0; i<Nfreq; i++){
			fprintf(fp, "%lf\t", wnp);
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				fprintf(fp, "%.16lf\t%.16lf\t\t", GSL_REAL(rself[i][mu][nu]), GSL_IMAG(rself[i][mu][nu]) );
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

		sprintf(parafp, "%s/lattice/tdelta_%c%.2lf_Nplot%d.dat", save_directory, control_char, control_para, Nfreq);
		fp = fopen(parafp, "w");
		for(i=0; i<Nfreq; i++){
			fprintf(fp, "%lf\t", wnp);
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				fprintf(fp, "%.16lf\t%.16lf\t\t", GSL_REAL(Glocalkinverse[i][mu][nu]), GSL_IMAG(Glocalkinverse[i][mu][nu]) );
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

		sprintf(parafp, "%s/lattice/%c%.3lf/tldos_Nplot%d_ep%.3lf_%dth.dat", save_directory, control_char, control_para, Nfreq, Epsilon, count);
		fp = fopen(parafp, "w");
		sum = 0;
		for(i=0; i<Nfreq; i++){
			fprintf(fp, "%lf\t", wnp);
			for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++){
				fprintf(fp, "%.16lf\t", Ldos[i][n*tNC+mu]);
				sum += Ldos[i][n*tNC+mu];
			}
			fprintf(fp, "%lf\n", sum*delta_plot);
		}
		fclose(fp);

		sprintf(parafp, "%s/lattice/%c%.3lf/tdelta_Nplot%d_ep%.3lf_%dth.dat", save_directory, control_char, control_para, Nfreq, Epsilon, count);
		fp = fopen(parafp, "w");
		for(i=0; i<Nfreq; i++){
			fprintf(fp, "%lf\t", wnp);
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				fprintf(fp, "%.16lf\t%.16lf\t\t", GSL_REAL(Glocalkinverse[i][mu][nu]), GSL_IMAG(Glocalkinverse[i][mu][nu]) );
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

		freegsctritensord(rself, Nfreq, NU);
	}
	freegsctetratensord(Self, Ni, Nfreq, tNC);
	freegscmatrixd(transform, NU);
	freematrixd(Ldos, Nfreq);
	freegsctritensord(Glocalkinverse, Nfreq, NU);
	freegsctritensord(Glocalk, Nfreq, NU);
	freegsctritensord(Glatt, Nfreq, NU);
	free(hopmatrix);

	return sum*delta_plot;
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
