#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"

#include "f2c.h"
#include "clapack.h"

#include "matrix.h"
#include "inverse.c"
#include "basic.h"
#include "lattice.h"
#include "ls.c"

#define Epsilon_default	0.02

#define Naxis_default	7
#define Nplot		1024
#define max_x		1.
#define Plot_range	16.
#define PLOT_BAND	1
#define Nkunit		100

#define delta_plot	(Plot_range/(Nplot-1))
#define wnp		fre[i]
#define PERIOD		1

void mklattice(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory);
void mkfermi(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory);
int continued_fraction_real(gsl_complex ***Green, int mu, int nu, int direction, double *ai, double *bi, double gm, double p0inner, int dimension);
int build_bath_real(int n, gsl_complex ***Green_bath, PNSpara egv);

int print_tail(FILE *fdx, int x, int y);
int print_head(FILE *fdx, int x, int y, char *paradata);
void print_complex(gsl_complex cmnumber);

long seed = -11;
//long seed = -134551;
double Hz, JTD, SOC, J, MU, DELTA, Ts, Uinter, *fre, Efermi, Epsilon;
int beta, Nmax, nctot, nb, tnb, myrank=0, beta_real, Nonzero_SOC, ROTATE, AF=0;
gsl_complex zero, ***Greennew, **Ecluster;

int main(int argc, char *argv[]){
	if(argc != 12 && argc != 13){
		printf("Usage :: %s <input.mak> <Hz> <JTD> <SOC> <J> <U_init> <U_final> <DELTA> <save_directory> <beta_real> <AF> <Epsilon: optional>\n", argv[0]);
		exit(1);
	}
	int i, j, k, count = 0;
	int degeneracy[Ni], tol;

	PNSpara egv[Ni];
	double U_init, U_final, re, im;
	double ***ginformation;

	FILE *fd; 
	char pathinform[1024], save_directory[1024];

	fre = mkvectord( Nplot );
	for(i=0; i<Nplot; i++)	fre[i] = (double)(i-Nplot/2) * delta_plot;

	nctot = NC;
	nb = NB;
	tnb = 2*NB;
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
	if( argc == 13 )	
		Epsilon = atof( argv[12] );
	else	Epsilon = Epsilon_default;


	if( fabs(SOC) < 1e-6 ){
		ROTATE = 0;
		Nonzero_SOC = 0;
	}
	else{
		ROTATE = 1;
		Nonzero_SOC = 1;
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
	print_gscmatrixd("Ecluster", Ecluster, NU, NU);

	gsl_complex ***locals = mkgsctritensord(Ni, tNC, tNC);
	gsl_complex **change = mkgscmatrixd(tNC, tNC);
	int n1, n2, mu, nu;
	for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)
		locals[n][mu][nu] = Ecluster[n*tNC+mu][n*tNC+nu];
	for(n1=0; n1<Ni; n1++) for(n2=n1+1; n2<Ni; n2++){
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)
			change[mu][nu] = gsl_complex_sub( locals[n1][mu][nu], locals[n2][mu][nu] );
		print_gscmatrixd("diff", change, tNC, tNC);
	} 

	freegsctritensord(locals, Ni, tNC);
	freegscmatrixd(change, tNC);

	sprintf(pathinform, "%s/green", save_directory);
	mkdirectory(pathinform);
	sprintf(pathinform, "%s/lattice/fermi", save_directory);
	mkdirectory(pathinform);
	sprintf(pathinform, "%s/lattice/vdx", save_directory);
	mkdirectory(pathinform);
	sprintf(pathinform, "%s/lattice/%c%.3lf", save_directory, control_char, control_para);
	mkdirectory(pathinform);

	mklattice(ginformation, egv, degeneracy, count, save_directory);
	//mkfermi(ginformation, egv, degeneracy, count, save_directory);

	free(fre);
	for(n=0; n<Ni; n++)	freematrixd(ginformation[n], degeneracy[n]);
	free(ginformation);
	freegscmatrixd(Ecluster, tNC);
	return 0;
}

int load_Green_spin(char *para, gsl_complex ******Green, int *degeneracy, double ***ginformation){
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

			continued_fraction_real(Green[n][degen][spin], mu, nu, direction, ai, bi, gm, p0inner, dimension);

			free(ai);	free(bi);
		}
	}while(1);
	fclose(fp);

	return 0;
}
int load_Green_nospin(char *para, gsl_complex *****Green, int *degeneracy, double ***ginformation){
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

			continued_fraction_real(Green[n][degen], mu, nu, direction, ai, bi, gm, p0inner, dimension);

			free(ai);	free(bi);
		}
	}while(1);
	fclose(fp);

	return 0;
}

int compute_self(double ***ginformation, PNSpara *egv, char *save_directory, int *degeneracy, int count, gsl_complex ****Self){
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
		Giplus[n]	= mkgscpentatensord(degeneracy[n], Spmax, nctot, nctot, Nplot);
		Green_on[n]	= mkgscpentatensord(degeneracy[n], Spmax, nctot, nctot, Nplot);
		tildekai[n]	= mkgsctetratensord(degeneracy[n], nctot, nctot, Nplot);
		kai[n]		= mkgsctetratensord(degeneracy[n], nctot, nctot, Nplot);
	}

	Ginverse	= mkgsctetratensord(Ni, Nplot, tNC, tNC);
	Green_total	= mkgsctetratensord(Ni, Nplot, tNC, tNC);

	Green_bath	= mkgsctetratensord(Ni, Nplot, tNC, tNC);
	Grplus		= mkgscmatrixd(Spmax, Nplot);

	char paragreenon[400], paragiplus[400], paratildekai[400], parakai[400];

	sprintf(paragreenon, "%s/datafile/greenon_u%.2lf_%dth.ab", save_directory, Uinter, count);
	sprintf(paragiplus, "%s/datafile/giplus_u%.2lf_%dth.ab", save_directory, Uinter, count);
	sprintf(paratildekai, "%s/datafile/tildekai_u%.2lf_%dth.ab", save_directory, Uinter, count);
	sprintf(parakai, "%s/datafile/kai_u%.2lf_%dth.ab", save_directory, Uinter, count);

	for(n=0; n<Ni; n++){
		for(degen=0; degen<degeneracy[n]; degen++) for(sp=0; sp<Spmax; sp++) for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nplot; i++){
			Green_on[n][degen][sp][mu][nu][i] = zero;
			Giplus[n][degen][sp][mu][nu][i] = zero;
		}
		for(degen=0; degen<degeneracy[n]; degen++) for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nplot; i++){
			tildekai[n][degen][mu][nu][i] = zero;
			kai[n][degen][mu][nu][i] = zero;
		}
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) for(i=0; i<Nplot; i++){
			Green_total[n][i][mu][nu] = zero;
		}
	}
	if( nctot > 1 ) load_Green_spin(paragiplus, Giplus, degeneracy, ginformation);
	load_Green_spin(paragreenon, Green_on, degeneracy, ginformation);
	if( Nonzero_SOC ){
		load_Green_nospin(paratildekai, tildekai, degeneracy, ginformation);
		load_Green_nospin(parakai, kai, degeneracy, ginformation);
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
				for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nplot; i++){
					Giplus[n][degen][sp][mu][nu][i] = gsl_complex_mul( Giplus[n][degen][sp][mu][nu][i], gsl_complex_rect(0, -1 ) );
				}
				for(mu=0; mu<nctot; mu++) for(nu=mu+1; nu<nctot; nu++) for(i=0; i<Nplot; i++){
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
				for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nplot; i++){
					Green_total[n][i][2*mu+sp][2*nu+sp] = gsl_complex_add( Green_total[n][i][2*mu+sp][2*nu+sp], gsl_complex_mul_real( Green_on[n][degen][sp][mu][nu][i], myweight) ) ;
				}
			}
			if( Nonzero_SOC ) for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nplot; i++){
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
		for(i=0; i<Nplot; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
			Green_total[n][i][mu][nu] = gsl_complex_div_real( Green_total[n][i][mu][nu], partition);
		}
		build_bath_real(n, Green_bath[n], egv[n]);
		inverse_complex_matrix_all_omp( Green_total[n], Ginverse[n], 0, Nplot, tNC);
		for(i=0; i<Nplot; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
			Self[n][i][mu][nu] = gsl_complex_sub( Green_bath[n][i][mu][nu], Ginverse[n][i][mu][nu] );
		}
	}

	FILE *ftemp;
	char paratemp[400];
	sprintf(paratemp, "%s/lattice/%c%.3lf/greenreal_ep%.3lf_%dth.dat", save_directory, control_char, control_para, Epsilon, count);
	ftemp = fopen(paratemp, "w");
	for(i=0; i<Nplot; i++){
		fprintf(ftemp, "%lf\t", wnp);
		for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
			fprintf(ftemp, "%lf\t%lf\t", GSL_REAL(Green_total[n][i][mu][nu]), GSL_IMAG(Green_total[n][i][mu][nu]));
		}
		fprintf(ftemp, "\n");
	}
	fclose(ftemp);

	sprintf(paratemp, "%s/lattice/%c%.3lf/selfreal_ep%.3lf_%dth.dat", save_directory, control_char, control_para, Epsilon, count);
	ftemp = fopen(paratemp, "w");
	for(i=0; i<Nplot; i++){
		fprintf(ftemp, "%lf\t", wnp);
		for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
			fprintf(ftemp, "%lf\t%lf\t", GSL_REAL(Self[n][i][mu][nu]), GSL_IMAG(Self[n][i][mu][nu]));
		}
		fprintf(ftemp, "\n");
	}
	fclose(ftemp);

	gsl_complex ****Self_t2g = mkgsctetratensord(Ni, Nplot, tNC, tNC);


	gsl_complex **transform = mkgscmatrixd(NU, NU);
	if( ROTATE ){
		init_transform(transform);
		for(n=0; n<Ni; n++) for(i=0; i<Nplot; i++)
			unitary(Self[n][i], Self_t2g[n][i], transform, tNC);
	}
	else{
		for(n=0; n<Ni; n++) for(i=0; i<Nplot; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)
			Self_t2g[n][i][mu][nu] = Self[n][i][mu][nu];
	}
	sprintf(paratemp, "%s/lattice/%c%.3lf/tselfreal_ep%.3lf_%dth.dat", save_directory, control_char, control_para, Epsilon, count);
	ftemp = fopen(paratemp, "w");
	for(i=0; i<Nplot; i++){
		fprintf(ftemp, "%lf\t", wnp);
		for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
			fprintf(ftemp, "%lf\t%lf\t", GSL_REAL(Self[n][i][mu][nu]), GSL_IMAG(Self[n][i][mu][nu]));
		}
		fprintf(ftemp, "\n");
	}
	fclose(ftemp);
	freegscmatrixd(transform, NU);

	freegsctetratensord(Self_t2g, Ni, Nplot, tNC);

	sprintf(paratemp, "%s/lattice/%c%.3lf/dos_ep%.3lf_%dth.dat", save_directory, control_char, control_para, Epsilon, count);
	ftemp = fopen(paratemp, "w");
	double local[Ni], sum, accu=0;
	for(i=0; i<Nplot; i++){
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

	for(n=0; n<Ni; n++){
		freegscpentatensord(Giplus[n],	degeneracy[n], Spmax, nctot, nctot);
		freegscpentatensord(Green_on[n],	degeneracy[n], Spmax, nctot, nctot);
		freegsctetratensord(tildekai[n],	degeneracy[n], nctot, nctot);
		freegsctetratensord(kai[n],		degeneracy[n], nctot, nctot);
	}
	free(Giplus); free(Green_on); free(tildekai); free(kai);
	freegsctetratensord(Green_bath, Ni, Nplot, tNC);
	freegsctetratensord(Green_total, Ni, Nplot, tNC);
	freegsctetratensord(Ginverse, Ni, Nplot, tNC);

	freegscmatrixd(Grplus, Spmax);

	return 0;
}

void mkfermi(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory){
	int n, mu, nu, i, k, axis;
	double dxG, kx, ky, kz;

	FILE *fdx, *fb;
	char paradx[1024];
	gsl_complex **Glocalkinverse, **Glocalk, **Gt2g, ****Self;

	Self		= mkgsctetratensord(Ni, Nplot, tNC, tNC);
	compute_self(ginformation, egv, save_directory, degeneracy, count, Self);

	double **kvec = mkmatrixd( 3, 3 );
	char knames[3];
	fb = fopen("inputs/FS", "r");
	nofile(fb, "FS");
	for(i=0; i<3; i++){
		for(k=0; k<3; k++)	fscanf(fb, "%lf", &kvec[i][k]);
		fscanf(fb, "%s", &knames[i]);
	}
	for(i=0; i<3; i++) for(k=0; k<3; k++)	kvec[i][k] *= 2*PI;
	fclose(fb);

	print_matrixd("kvec", kvec, 3, 3);
	print_vectorc("knames", knames, 3);


	double **bvec = mkmatrixd( 3, 3 );
	double bvec_default[3][3] = {
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1}
	};
	fb = fopen("inputs/BVECTORS", "r");
	if( fb != NULL ){
		for(i=0; i<3; i++) for(k=0; k<3; k++)	fscanf(fb, "%lf", &bvec[i][k]);
		fclose(fb);
	}
	else{
		for(i=0; i<3; i++) for(k=0; k<3; k++)	bvec[i][k] = bvec_default[i][k];
	}

	int same, Nkx, Nky, Nkdx[2];
	double displace[3], dis_tot, dis[2];
	for(k=0; k<2; k++){
		kx = kvec[k+1][0] - kvec[0][0];
		ky = kvec[k+1][1] - kvec[0][1];
		kz = kvec[k+1][2] - kvec[0][2];
		for(i=0; i<3; i++)	displace[i] = kx*bvec[0][i] + ky*bvec[1][i] + kz*bvec[2][i];
		dis[k] = sqrt(displace[0]*displace[0] + displace[1]*displace[1] + displace[2]*displace[2]);
		dis_tot += dis[k];
	}
	for(k=0; k<2; k++){
		Nkdx[k] = (int)(dis[k]/dis_tot * 300);
	}
	Nkx = Nkdx[0];
	Nky = Nkdx[1];
	int Nk_total=Nkx*Nky;
	double dk[2][3];
	Glocalkinverse	= mkgscmatrixd( NU, NU );
	Glocalk		= mkgscmatrixd( NU, NU );
	Gt2g		= mkgscmatrixd( NU, NU );

	gsl_complex ***hopmatrix = mkgsctritensord( Nk_total, NU, NU );
	gsl_complex *oneline = mkgscvectord( (Nk_total)*NU*NU );

	int line;
	Latt *data;
	Quad quadrature;
	data = read_lattice(WANNIER, &line, &Efermi);
	init_quadrature(&quadrature);
	if( fabs(Uinter)<1e-6 )	find_Ef(data, &quadrature, line, 5./6.);
	char parahop[1024];
	FILE *fhop;
	sprintf(parahop, "%s/mdisp/fermi_Nkx%d_Nky%d", Path_hop, Nkx, Nky);
	for(i=0; i<3; i++){
		sprintf(parahop, "%s%c", parahop, knames[i]);
	}
	sprintf(parahop, "%s.hop", parahop);

	fhop = fopen(parahop, "rb");
	if( fhop == NULL ){
		for(k=0; k<3; k++){
			dk[0][k] = (kvec[1][k] - kvec[0][k])/(Nkx-1);
			dk[1][k] = (kvec[2][k] - kvec[0][k])/(Nky-1);
		}
		for(i=0; i<Nkx; i++) for(k=0; k<Nky; k++){
			kx = kvec[0][0] + i*dk[0][0] + k*dk[1][0];
			ky = kvec[0][1] + i*dk[0][1] + k*dk[1][1];
			kz = kvec[0][2] + i*dk[0][2] + k*dk[1][2];
			compute_hopmatrix(data, hopmatrix[i*Nky+k], kx, ky, kz, line);
		}
		i=0;
		for(k=0; k<Nk_total; k++){
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				oneline[i] = hopmatrix[k][mu][nu];
				i++;
			}
		}
		fhop = fopen(parahop, "wb");
		fwrite(oneline, sizeof(gsl_complex), (Nk_total)*NU*NU, fhop);
		fclose(fhop);
	}
	else{
		fhop = fopen(parahop, "rb");
		fread(oneline, sizeof(gsl_complex), (Nk_total)*NU*NU, fhop);
		fclose(fhop);

		i=0;
		for(k=0; k<Nk_total; k++){
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				hopmatrix[k][mu][nu] = oneline[i];
				i++;
			}
		}
	}
	free(oneline);
	double distortion[Ni][tNC];
	init_distortion(distortion);
	for(k=0; k<Nk_total; k++) for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++)
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
		for(i=0; i<Nk_total; i++)	iunitary(hopmatrix[i], hopmatrix[i], transform, NU);
	}

	FILE *fd;
	char paradata[1024], pararoot[1024];

	sprintf(paradx, "%s/lattice/fermi/Nkx%d_Nky%d_ep%.2lf_%dth.dx", save_directory, Nkx, Nky, Epsilon, count);
	fdx = fopen(paradx, "w");

	sprintf(pararoot, "u%.2lf_ts%.2lf_d%.3lf_ep%.2lf_%dth.txt", Uinter, Ts, DELTA, Epsilon, count);
	fprintf(fdx, "Dim\t%d\n", 2);

	fprintf(fdx, "KNAMES\t");
	for(i=0; i<3; i++)	fprintf(fdx, "%c ", knames[i]);
	fprintf(fdx, "\n");
	fprintf(fdx, "Nx\t%d\n", Nkx);
	fprintf(fdx, "Ny\t%d\n", Nky);
	fprintf(fdx, "datafile\t%s\n", pararoot);
	fprintf(fdx, "\n");
	fprintf(fdx, "xrange\t%d\t%d\n", 0, Nkx);
	fprintf(fdx, "yrange\t%d\t%d\n", 0, Nky);

	sprintf(paradata, "%s/lattice/fermi/%s", save_directory, pararoot);
	fd = fopen(paradata, "w");

	char parasingle[1024];
	sprintf(parasingle, "%s/lattice/fermi/band_Nk%d.dat", save_directory, Nk_total);
	FILE *fs = fopen(parasingle, "w");
	double *eval;
	gsl_complex **evec;
	if( PLOT_BAND ){
		fprintf(fdx, "BAND\tlattice/fermi/band_Nk%d.dat\n", Nk_total);
		eval = mkvectord(NU);
		evec = mkgscmatrixd(NU, NU);
	}
	fclose(fdx);
	double ***decomposed = mktritensord(2, NU, Nk_total);

	double diag[NU];
	update_diag(diag);
	i=Nplot/2;
	for(axis=0; axis<Nkx; axis++) {
		for(k=0; k<Nky; k++) {
			dxG = 0;
			//fprintf(ftemp, "%lf\t", wnp);
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
				same = (int)(mu==nu);
				Glocalkinverse[mu][nu]
					= gsl_complex_sub(
							gsl_complex_rect( (wnp-diag[mu]) *same , Epsilon*same),
							hopmatrix[axis*Nky + k][mu][nu]
							);
			}

			for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) 
				Glocalkinverse[tNC*n+mu][tNC*n+nu]
					= gsl_complex_rect(
							GSL_REAL( Glocalkinverse[tNC*n+mu][tNC*n+nu] ) - GSL_REAL(Self[n][i][mu][nu]) ,
							GSL_IMAG( Glocalkinverse[tNC*n+mu][tNC*n+nu] ) - GSL_IMAG(Self[n][i][mu][nu]) 
							);
			inverse_complex_matrix_lapack( Glocalkinverse, Glocalk, NU);

			for(mu=0; mu<NU; mu++){
				dxG -= GSL_IMAG(Glocalk[mu][mu]) ;
				decomposed[0][mu][axis*Nky+k] = -GSL_IMAG(Glocalk[mu][mu])/PI ;
			}
			fprintf(fd, "%.10lf\t", dxG/PI/NU);
			if( ROTATE ){
				unitary(Glocalk, Gt2g, transform, NU);
				for(mu=0; mu<NU; mu++){
					decomposed[1][mu][axis*Nky+k] = -GSL_IMAG(Gt2g[mu][mu])/PI ;
				}
			}

			if( PLOT_BAND ){
				eigen_lapack(hopmatrix[axis*Nky+k], eval, evec, NU, 1);
				for(mu=0; mu<NU; mu++) if( fabs(eval[mu]-Efermi) < 1e-3 )
					fprintf(fs, "%lf\t%lf\n", (double)axis/Nkx, (double)k/Nky);
			}
		}
		fprintf(fd, "\n");
	}
	fclose(fd);
	if( PLOT_BAND)	fclose(fs);

	if(1){
		sprintf(pararoot, "decomp_u%.2lf_ts%.2lf_d%.3lf_ep%.2lf_%dth.txt", Uinter, Ts, DELTA, Epsilon, count);
		sprintf(paradata, "%s/lattice/fermi/%s", save_directory, pararoot);
		fd = fopen(paradata, "w");
		for(mu=0; mu<NU; mu++){
			for(k=0; k<Nk_total; k++)
				fprintf(fd, "%lf\t", decomposed[0][mu][k]);
			fprintf(fd, "\n");
		}
		fprintf(fd, "\n");

		if( ROTATE ) for(mu=0; mu<NU; mu++){
			for(k=0; k<Nk_total; k++)
				fprintf(fd, "%lf\t", decomposed[1][mu][k]);
			fprintf(fd, "\n");
		}
		fclose(fd);

		sprintf(paradx, "%s/lattice/fermi/decomp_Nkx%d_Nky%d_ep%.2lf_%dth.dx", save_directory, Nkx, Nky, Epsilon, count);
		fs = fopen(paradx, "w");

		fprintf(fs, "Ni\t%d\n", Ni);
		for(n=0; n<Ni; n++)
			fprintf(fs, "No%d\t%d\n", n, tNC);
		fprintf(fs, "KNAMES\t");
		for(i=0; i<3; i++)	fprintf(fs, "%c ", knames[i]);
		fprintf(fs, "\n");
		fprintf(fs, "Nx\t%d\n", Nkx);
		fprintf(fs, "Ny\t%d\n", Nky);
		fprintf(fs, "datafile\t%s\n", pararoot);
		fprintf(fs, "\n");
		fprintf(fs, "xrange\t%.18lf\t%.18lf\n", 0., 1.);
		fprintf(fs, "yrange\t%.18lf\t%.18lf\n", 0., 1.);
		fclose(fs);
		freegscmatrixd(transform, NU);
	}
	freetritensord(decomposed, 2, NU);
	freegscmatrixd(evec, NU);
	free(eval);
	free(data);
	freematrixd(kvec, 3);

	freegscmatrixd(Glocalkinverse, NU);
	freegscmatrixd(Glocalk, NU);

	freegsctetratensord(Self, Ni, Nplot, tNC);
	freegsctritensord(hopmatrix, Nk_total, NU);


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

int compute_Glocalk(gsl_complex **Glocalk, gsl_complex **Glocalkinverse, gsl_complex ****Self, gsl_complex **hopmatrix, double *diag, int i){
	int n, mu, nu, same;
	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
		same = (int)(mu==nu);
		Glocalkinverse[mu][nu]
			= gsl_complex_sub(
					gsl_complex_rect( (wnp-diag[mu]) *same , Epsilon*same),
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

void mklattice(double ***ginformation, PNSpara *egv, int *degeneracy, int count, char *save_directory){
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
	double dxG, kx, ky, kz;

	FILE *fdx;
	char paradx[1024], indices_default[Naxis_default] = {'L', 'G', 'X', 'W', 'L', 'K','G'};
	gsl_complex ***Glocalkinverse, ***Glocalk, ***Gt2g, ****Self;

	Self		= mkgsctetratensord(Ni, Nplot, tNC, tNC);
	compute_self(ginformation, egv, save_directory, degeneracy, count, Self);

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
	Glocalkinverse	= mkgsctritensord( Nplot, NU, NU );
	Glocalk		= mkgsctritensord( Nplot, NU, NU );
	Gt2g		= mkgsctritensord( Nplot, NU, NU );

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
	FILE *fk, *fhop;
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

	double diag[NU];
	update_diag(diag);
	print_vectord("diag", diag, NU);
	gsl_complex **Gpd = mkgscmatrixd(tNC, tNC);

	for(axis=0; axis<Naxis; axis++){
		sprintf(parak, "%s/lattice/vdx/AF%d_k%s_u%.2lf_ts%.2lf_d%.3lf_ep%.2lf_%dth.dat", save_directory, AF, indices[axis], Uinter, Ts, DELTA, Epsilon, count);
#pragma omp parallel for default(none) private(n, mu,nu, i) shared(hopmatrix,Glocalk,Glocalkinverse,Self,fre,Nkrefer,axis, diag, ROTATE, transform, Gt2g)
		for(i=0; i<Nplot; i++){
			compute_Glocalk(Glocalk[i], Glocalkinverse[i], Self, hopmatrix[Nkrefer[axis]], diag, i);
			if( ROTATE ) unitary(Glocalk[i], Gt2g[i], transform, NU);
		}
		kx = sympoint[axis][0];
		ky = sympoint[axis][1];
		kz = sympoint[axis][2];

		fk = fopen(parak, "w");
		for(i=0; i<Nplot; i++){
			dxG = 0;
			if( PERIOD ){
				periodize_green(Glocalk[i], Gpd, kx, ky, kz, pos);
				for(mu=0; mu<tNC; mu++) dxG -= GSL_IMAG(Gpd[mu][mu])/tNC ;
				fprintf(fk, "%.10lf\t%.10lf\t\t", wnp, dxG/PI);
				for(mu=0; mu<tNC; mu++)	fprintf(fk, "%.10lf\t", -GSL_IMAG(Gpd[mu][mu])/PI);
			}
			else{
				for(mu=0; mu<NU; mu++) dxG -= GSL_IMAG(Glocalk[i][mu][mu])/NU ;
				fprintf(fk, "%.10lf\t%.10lf\t\t", wnp, dxG/PI);
				for(mu=0; mu<NU; mu++)	fprintf(fk, "%.10lf\t", -GSL_IMAG(Glocalk[i][mu][mu])/PI);
			}

			fprintf(fk, "\n");
		}
		fclose(fk);
		if( ROTATE ){
			sprintf(parak, "%s/lattice/vdx/tk%s_u%.2lf_ts%.2lf_d%.3lf_ep%.2lf_%dth.dat", save_directory, indices[axis], Uinter, Ts, DELTA, Epsilon, count);
			fk = fopen(parak, "w");
			for(i=0; i<Nplot; i++){
				dxG = 0;
				if( PERIOD ){
					periodize_green(Gt2g[i], Gpd, kx, ky, kz, pos);
					for(mu=0; mu<tNC; mu++) dxG -= GSL_IMAG(Gpd[mu][mu])/tNC ;
					fprintf(fk, "%.10lf\t%.10lf\t\t", wnp, dxG/PI);
					for(mu=0; mu<tNC; mu++)	fprintf(fk, "%.10lf\t", -GSL_IMAG(Gpd[mu][mu])/PI);
				}
				else{
					for(mu=0; mu<NU; mu++) dxG -= GSL_IMAG(Gt2g[i][mu][mu])/NU ;
					fprintf(fk, "%.10lf\t%.10lf\t\t", wnp, dxG/PI);
					for(mu=0; mu<NU; mu++)	fprintf(fk, "%.10lf\t", -GSL_IMAG(Gt2g[i][mu][mu])/PI);
				}
				fprintf(fk, "\n");
			}
			fclose(fk);
		}
	}
	freegscmatrixd(Gpd, tNC);

	FILE *fd;
	char paradata[1024], pararoot[1024];

	int sum=0;
	for(axis=0; axis<Naxis-1; axis++){
		printf("%s: %lf\n", indices[axis], max_x/Nk_total*sum);
		sum += Nkdx[axis];
	}

	sprintf(paradx, "%s/lattice/vdx/naxis%d_Nplot%d_Nk%d_ep%.2lf_%dth.dx", save_directory, Naxis, Nplot, Nk_total+1, Epsilon, count);
	fdx = fopen(paradx, "w");

	sprintf(pararoot, "u%.2lf_ts%.2lf_d%.3lf_ep%.2lf_%dth.txt", Uinter, Ts, DELTA, Epsilon, count);
	fprintf(fdx, "KNAMES\t");
	for(axis=0; axis<Naxis; axis++) fprintf(fdx, "%s\t", indices[axis]);
	fprintf(fdx, "\nKPOINTS\t");
	for(axis=0; axis<Naxis-1; axis++) fprintf(fdx, "%d\t", Nkdx[axis]);
	fprintf(fdx, "\n");

	char parasingle[1024];
	sprintf(parasingle, "%s/lattice/vdx/band_naxis%d_Nk%d.dat", save_directory, Naxis, Nk_total+1);
	FILE *fs = fopen(parasingle, "w");
	double *eval;
	gsl_complex **evec;
	if( PLOT_BAND ){
		fprintf(fdx, "BAND\tlattice/vdx/band_naxis%d_Nk%d.dat\n", Naxis, Nk_total+1);
		eval = mkvectord(NU);
		evec = mkgscmatrixd(NU, NU);
	}
	print_head(fdx, Nplot, Nk_total+1, pararoot);
	print_tail(fdx,	Nplot, Nk_total+1);
	double ****decomposed	= mktetratensord(2, NU, Nplot, Nk_total+1);
	double ****decom_pd	= mktetratensord(2, tNC, Nplot, Nk_total+1);


	int onemore = 1;
#pragma omp parallel for default(none) private(Gpd, i, axis, k, onemore, kx, ky, kz, mu, nu, dk) shared(Naxis, Efermi, pos, sympoint, hopmatrix,Glocalk,Glocalkinverse,Self,fre,Nkrefer,Nkdx, diag, ROTATE, transform, Gt2g, fs, decomposed, decom_pd, eval, evec)
	for(i=0; i<Nplot; i++) {
		Gpd = mkgscmatrixd(tNC, tNC);
		for(axis=0; axis<Naxis-1; axis++) {
			for(k=Nkrefer[axis]; k<Nkrefer[axis]+Nkdx[axis]; k++){
				if( !onemore && axis == Naxis-2 && k == Nkrefer[axis]+Nkdx[axis]-1 )	k++;
				compute_Glocalk(Glocalk[i], Glocalkinverse[i], Self, hopmatrix[k], diag, i);

				for(mu=0; mu<NU; mu++) decomposed[0][mu][i][k] = -GSL_IMAG(Glocalk[i][mu][mu])/PI ;

				for(int j=0; j<3; j++)	dk[j] = sympoint[axis+1][j]-sympoint[axis][j];
				kx = sympoint[axis][0] + (k-Nkrefer[axis])*dk[0]/Nkdx[axis];
				ky = sympoint[axis][1] + (k-Nkrefer[axis])*dk[1]/Nkdx[axis];
				kz = sympoint[axis][2] + (k-Nkrefer[axis])*dk[2]/Nkdx[axis];

				if( PERIOD ){
					periodize_green(Glocalk[i], Gpd, kx, ky, kz, pos);
					for(mu=0; mu<tNC; mu++) decom_pd[0][mu][i][k] = -GSL_IMAG( Gpd[mu][mu] )/PI;
				}
				if( ROTATE ){
					unitary(Glocalk[i], Gt2g[i], transform, NU);
					for(mu=0; mu<NU; mu++) decomposed[1][mu][i][k] = -GSL_IMAG(Gt2g[i][mu][mu])/PI ;
					if( PERIOD ){
						periodize_green(Gt2g[i], Gpd, kx, ky, kz, pos);
						for(mu=0; mu<tNC; mu++) decom_pd[1][mu][i][k] = -GSL_IMAG(Gpd[mu][mu])/PI ;
					}
				}
				if( PLOT_BAND && i==0 ){
					eigen_lapack(hopmatrix[k], eval, evec, NU, 1);
					for(mu=0; mu<NU; mu++)	fprintf(fs, "%lf\t", eval[mu]-Efermi);
					fprintf(fs, "\n");
				}
			}
			if( axis == Naxis-2 && k == Nkrefer[axis]+Nkdx[axis]-1 && onemore ){
				k--;
				onemore = 0;
			}
		}
		freegscmatrixd(Gpd, tNC);
	}
	if( PLOT_BAND)	fclose(fs);

	sprintf(paradata, "%s/lattice/vdx/%s", save_directory, pararoot);
	fd = fopen(paradata, "w");
	for(i=0; i<Nplot; i++) for(k=0; k<Nk_total+1; k++){
		dxG = 0;
		if( PERIOD )	for(mu=0; mu<tNC; mu++) dxG += decom_pd[0][mu][i][k];
		else		for(mu=0; mu<NU; mu++) dxG += decomposed[0][mu][i][k];
		fprintf(fd, "%.10lf\t", dxG);
	}
	fclose(fd);

	if(1){
		sprintf(pararoot, "decomp_u%.2lf_ts%.2lf_d%.3lf_ep%.2lf_%dth.txt", Uinter, Ts, DELTA, Epsilon, count);
		sprintf(paradata, "%s/lattice/vdx/%s", save_directory, pararoot);
		fd = fopen(paradata, "w");
		for(mu=0; mu<NU; mu++){
			for(i=0; i<Nplot; i++) for(k=0; k<Nk_total+1; k++)
				fprintf(fd, "%lf\t", decomposed[0][mu][i][k]);
			fprintf(fd, "\n");
		}
		fprintf(fd, "\n");

		if( ROTATE ) for(mu=0; mu<NU; mu++){
			for(i=0; i<Nplot; i++) for(k=0; k<Nk_total+1; k++)
				fprintf(fd, "%lf\t", decomposed[1][mu][i][k]);
			fprintf(fd, "\n");
		}
		fclose(fd);

		sprintf(paradx, "%s/lattice/vdx/decomp_naxis%d_Nplot%d_Nk%d_ep%.2lf_%dth.dx", save_directory, Naxis, Nplot, Nk_total+1, Epsilon, count);
		fs = fopen(paradx, "w");

		fprintf(fs, "Ni\t%d\n", Ni);
		for(n=0; n<Ni; n++)
			fprintf(fs, "No%d\t%d\n", n, tNC);
		fprintf(fs, "Nx\t%d\n", Nk_total+1);
		fprintf(fs, "Ny\t%d\n", Nplot);
		fprintf(fs, "datafile\t%s\n", pararoot);
		fprintf(fs, "KNAMES\t");
		for(axis=0; axis<Naxis; axis++) fprintf(fs, "%s\t", indices[axis]);
		fprintf(fs, "\nKPOINTS\t");
		for(axis=0; axis<Naxis-1; axis++) fprintf(fs, "%d\t", Nkdx[axis]);
		fprintf(fs, "\n");
		fprintf(fs, "xrange\t%.18lf\t%.18lf\n", 0., 1.);
		fprintf(fs, "yrange\t%.18lf\t%.18lf\n", fre[0], fre[Nplot-1]);
		fclose(fs);
		freegscmatrixd(transform, NU);
	}
	freetetratensord(decomposed, 2, NU, Nplot);
	freegscmatrixd(evec, NU);
	free(eval);
	free(data);
	freematrixd(sympoint, Naxis);
	freematrixd(bvec, 3);
	freematrixc(indices, Naxis);

	freegsctritensord(Glocalkinverse, Nplot, NU);
	freegsctritensord(Glocalk, Nplot, NU);

	freegsctetratensord(Self, Ni, Nplot, tNC);
	freegsctritensord(hopmatrix, Nk_total+1, NU);
	freematrixd(pos, Ni);
	freetetratensord(decom_pd, 2, tNC, Nplot);
}

void print_complex(gsl_complex cmnumber){
	printf("%lf + %lf i\n", GSL_REAL(cmnumber), GSL_IMAG(cmnumber) );
}

int print_head(FILE *fdx, int x, int y, char *paradata){
	fprintf(fdx, "object 1 class array items\t%d data file %s\n", x*y, paradata);
	return 0;
}
int print_tail(FILE *fdx, int x, int y){
	int i;

	fprintf(fdx, "\n");
	fprintf(fdx, "attribute \"dep\" string \"positions\"\n");
	fprintf(fdx, "\nobject 2 class gridpositions counts  %d\t%d\n", x, y);
	i=0;
	fprintf(fdx, "origin   %.10lf   %.10lf\n", 0., wnp);
	fprintf(fdx, "delta    0\t%.10lf\n", (double)delta_plot);
	fprintf(fdx, "delta    %.10lf\t0\n", (double)max_x/(y-1) );
	fprintf(fdx, "\nobject 3 class gridconnections counts  %d  %d\n", x, y);
	fprintf(fdx, "attribute \"elements type\" string \"cubes\"\n");
	fprintf(fdx, "attribute \"ref\" string \"positions\"\n");
	fprintf(fdx, "\nobject \"density\" class field\n");
	fprintf(fdx, "component \"data\" 1\n");
	fprintf(fdx, "component \"positions\" 2\n");
	fprintf(fdx, "component \"connections\" 3\n");

	fclose(fdx);
	return 0;
}

int continued_fraction_real(gsl_complex ***Green, int mu, int nu, int direction, double *ai, double *bi, double gm, double p0inner, int dimension){
	int i, k;
	gsl_complex tempg;
	for(i=0; i<Nplot; i++) {
		tempg = zero;
		for(k=dimension-1; k>0; k--){
			tempg = gsl_complex_div(
					gsl_complex_rect(bi[k], 0),
					gsl_complex_sub(
						gsl_complex_rect(-1*ai[k] * direction + gm * direction + wnp, Epsilon ),
						tempg
					)
				);
		}
		tempg = gsl_complex_sub(
				gsl_complex_rect(-1*ai[0] * direction + gm * direction + wnp, Epsilon ),
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
}//end of continued_fraction_real

int build_bath_real(int n, gsl_complex ***Green_bath, PNSpara egv){
	gsl_complex tempGreen;
	int mu, nu, k, i;
	for(i=0; i<Nplot; i++) {
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
							gsl_complex_rect( wnp - egv.egbath[k], Epsilon )
						)
					);
					//(-) sign is included.
			}
			Green_bath[i][mu][nu]
				= gsl_complex_add(
						gsl_complex_sub(
							gsl_complex_rect( wnp*(mu==nu), Epsilon*(mu==nu) ),
							Ecluster[tNC*n+mu][tNC*n+nu]
						),
						tempGreen
					);
		}
	}
	return 0;
}//end of build_bath_real

