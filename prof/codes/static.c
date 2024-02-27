#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "matrix.h"
#include "basic.h"
#include "lattice.h"
#include "ls_basis.c"
#include "ls.c"

int compute_static(typebasis **basis, int count, char *save_directory);

//long seed = -134551;
double J, MU, DELTA, Uinter, Hz, JTD, SOC;
int beta, Nmax, nctot, nb, tnb, Spinmax, Blocks,Powns, Powns2, myrank=0, beta_real, Nonzero_SOC, ROTATE=0, AF=0;
gsl_complex zero, ***Greennew, **Ecluster;

int main(int argc, char *argv[]){
	if(argc != 12){
		printf("Usage :: %s <input.mak> <Hz> <JTD> <SOC> <J> <U_init> <U_final> <DELTA> <save_directory> <beta_real> <AF>\n", argv[0]);
		exit(1);
	}
	sectorinfo **ginfo;
	int i, n, count = 0;
	double U_init, U_final;

	FILE *fd; 
	char pathinform[1024], save_directory[1024];

	nctot = NC;
	nb = NB;
	tnb = 2*NB;
	Powns	= (int)pow(2,Ns);
	Powns2	= (int)pow(2,2*Ns);
	Blocks	= 2*Ns+1;
	Spinmax	= 2*Ns+1;
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

	if( fabs(SOC) < 1e-6 ){
		ROTATE = 0;
		Nonzero_SOC = 0;
		Blocks = (Ns+1)*(Ns+1);
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
	fclose(fd);

	typebasis **basis;
	basis = mkmatrixb(Powns2, 2);
	ginfo = (sectorinfo **) malloc( Ni * sizeof(sectorinfo*));
	for(n=0; n<Ni; n++){
		ginfo[n] = (sectorinfo *) malloc( Blocks*sizeof(sectorinfo) );
		build_bitbasis(Ns, basis, ginfo[n], Nonzero_SOC);
		print_ginfo(Ns, ginfo[n], Blocks);
	}
	printf("U = %.2lf, MU = %.2lf (of UF%.2lf UT%.2lf) from '%s'\n", Uinter, MU, U_init, U_final, argv[1]);

	/*
	Ecluster = mkgscmatrixd( NU, NU);
	char paraEcluster[1024];
	sprintf(paraEcluster, "%s/Ecluster.dat", save_directory);
	load_Ecluster(paraEcluster, Ecluster);
	update_chem(Ecluster);
	*/

	sprintf(pathinform, "%s/result/%c%.3lf", save_directory, control_char, control_para);
	mkdirectory(pathinform);

	compute_static(basis, count, save_directory);

	for(n=0; n<Powns2; n++)	free(basis[n]);		free(basis);
	for(n=0; n<Ni; n++){
		for(i=0; i<Blocks; i++)
			freegscmatrixd(ginfo[n][i].ground, Nvec);
		free( ginfo[n] );
	}
	free(ginfo);

	return 0;
}

int read_mak(PNSpara *egv, int *degeneracy, double ***ginformation, char *save_directory, int count){
	int tol, n, i, j, k;
	char para[1024];
	double re, im;
	FILE *fd;
	sprintf(para, "%s/%c%.2lf_%dth.mak", save_directory, control_char, control_para, count);
	//printf("trying to read %s\n", para);

	fd = fopen(para, "r");
	if( fd == NULL ){
		printf("%s does not exist!\n", para);
		return 1;
	}
	fscanf(fd, "%d", &count);
	fscanf(fd, "%lf", &Uinter);
	fscanf(fd, "%lf", &MU);
	fscanf(fd, "%d", &tol);
	fscanf(fd, "%d", &beta);
	fscanf(fd, "%d", &Nmax);

	for(n=0; n<Ni; n++){
		fscanf(fd, "%d", &degeneracy[n]);

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

	return 0;
}

int compute_static(typebasis **basis, int count, char *save_directory){
	PNSpara egv[Ni];
	typebasis bin[2];
	int n, degen, gndblock, ic, degeneracy[Ni], mu, nu, i, gnd, phase, *table, index, final_nonzero=0, kramers[tNC];
	char paragnd[1024], paramag[1024], paradouble[1024], parafilling[1024], paraext[1024], parat2g[1024], paraalpha[1024];
	FILE *fmag, *ffil, *ffil_t2g, *fdoc, *fext, *fground, *ft2g, *falpha;
	gsl_complex *ground;
	gsl_complex ***excitonic = mkgsctritensord(Ni, tNC, tNC);
	double lsum;
	double myweight, partition = 0, myenergy, gm;
	gsl_complex ***excitonic_t2g = mkgsctritensord(Ni, tNC, tNC);
	double doccu[Ni][NC], filling[Ni][NC][2], mag[Ni][NC], sum[3], magfactor[tNC];
	double ***ginformation;
	gsl_complex **transform = mkgscmatrixd(NU, NU);
	if( ROTATE ){
		init_transform(transform);
		double jeff_pair[6] = {0, 1, 1, 0, 2, 2};
		double jeff[6] = {1.5, 0.5, -0.5, -1.5, 0.5, -0.5};
		for(mu=0; mu<tNC; mu++){
			kramers[mu] = jeff_pair[mu];
			magfactor[mu] = jeff[mu];
		}
	}
	else{
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	transform[mu][nu] = zero;
		for(mu=0; mu<tNC; mu++)	transform[mu][mu] = gsl_complex_rect(1, 0);
		for(mu=0; mu<tNC; mu++){
			kramers[mu] = mu/2;
			magfactor[mu] = 0.5 * pow(-1, mu);
		}
	}
	ginformation = mktritensord(Ni, 100, 3);

	sprintf(paramag, "%s/result/%c%.3lf/mag.dat", save_directory, control_char, control_para);
	sprintf(paraalpha, "%s/result/%c%.3lf/alpha.dat", save_directory, control_char, control_para);
	sprintf(paradouble, "%s/result/%c%.3lf/double.dat", save_directory, control_char, control_para);
	sprintf(parafilling, "%s/result/%c%.3lf/filling.dat", save_directory, control_char, control_para);
	sprintf(paraext, "%s/result/%c%.3lf/excitonic.dat", save_directory, control_char, control_para);
	sprintf(parat2g, "%s/result/%c%.3lf/t2g_excitonic.dat", save_directory, control_char, control_para);
	fmag = fopen(paramag, "w");
	falpha = fopen(paraalpha, "w");
	fdoc = fopen(paradouble, "w");
	ffil = fopen(parafilling, "w");
	fext = fopen(paraext, "w");
	ft2g = fopen(parat2g, "w");
	sprintf(parafilling, "%s/result/%c%.3lf/t2g_filling.dat", save_directory, control_char, control_para);
	ffil_t2g = fopen(parafilling, "w");


	table = mkvectori( Powns2 );
	for(i=0; i<Powns2; i++){
		index = (basis[i][0]<<Ns) + basis[i][1];
		table[index] = i;
	}

	for(n=0; n<Ni; n++){
		egv[n].egbath = mkvectord(tnb);	
		egv[n].hybrid = mkgscmatrixd(tNC, tnb);
	}
	for(ic=1; ic<count+1; ic++){
		if( read_mak(egv, degeneracy, ginformation, save_directory, ic) ){
			if( ic == count ){
				printf("no mak file exists!\n");
				exit(1);
			}
			continue;
		}

		fprintf(falpha, "%d\t", ic);
		fprintf(fmag, "%d\t", ic);
		fprintf(fdoc, "%d\t", ic);
		fprintf(ffil, "%d\t", ic);
		fprintf(ffil_t2g, "%d\t", ic);
		fprintf(fext, "%d\t", ic);
		fprintf(ft2g, "%d\t", ic);
		sum[0] = 0;
		sum[1] = 0;
		sum[2] = 0;
		for(n=0; n<Ni; n++){
			sprintf(paragnd, "%s/u%.2lf_%dth_n%d.gnd", save_directory, Uinter, ic, n);
			fground = fopen(paragnd, "rb");
			nofile(fground, paragnd);

			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
				excitonic[n][mu][nu] = zero;
				excitonic_t2g[n][mu][nu] = zero;
			}
			for(mu=0; mu<nctot; mu++){
				doccu[n][mu] = 0;
				filling[n][mu][0] = 0;
				filling[n][mu][1] = 0;
				mag[n][mu] = 0;
			}

			gm = ginformation[n][0][0];
			partition = 0;
			for(degen=0; degen<degeneracy[n]; degen++){
				myenergy = ginformation[n][degen][0];
				if( beta_real )	myweight	= exp( - ( myenergy - gm )*beta_real );
				else		myweight	= 1;
				partition += myweight;
				gnd = (int)ginformation[n][degen][1];
				gndblock = (int)ginformation[n][degen][2];

				ground = mkgscvectord( gndblock );
				fread(ground, sizeof(gsl_complex), gndblock, fground);
				printf("ic=%d, n=%d, gnd=%d, gndblock=%d, degen=%d\n", ic, n, gnd, gndblock, degen);

				for(mu=0; mu<nctot; mu++){
					for(i = gnd; i<gnd+gndblock; i++){	
						doccu[n][mu] += myweight * gsl_complex_abs2(ground[i-gnd])* One(basis[i][0], mu) * One(basis[i][1], mu);
						mag[n][kramers[2*mu]]	+= magfactor[2*mu] * myweight * gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][1], mu));
						mag[n][kramers[2*mu+1]]	+= magfactor[2*mu+1] * myweight * gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][0], mu));

						filling[n][mu][0] += myweight * gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][0], mu) );
						filling[n][mu][1] += myweight * gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][1], mu) );
						excitonic[n][2*mu][2*mu]	= gsl_complex_add_real( excitonic[n][2*mu][2*mu],	myweight * gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][0], mu) ) );
						excitonic[n][2*mu+1][2*mu+1]	= gsl_complex_add_real( excitonic[n][2*mu+1][2*mu+1],	myweight * gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][1], mu) ) );
					}
				}
				for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
					for(i = gnd; i<gnd+gndblock; i++){	
						if( Zero( basis[i][mu%2], mu/2) && One( basis[i][nu%2], nu/2 ) ){
							bin[0] = basis[i][0];
							bin[1] = basis[i][1];
							phase  = permu( bin, nu/2, nu%2);	bin[nu%2] = Turnoff( bin[nu%2], nu/2 );
							phase *= permu( bin, mu/2, mu%2);	bin[mu%2] = Turnon( bin[mu%2], mu/2 );
							index = (bin[0]<<Ns) + bin[1];
							if( table[index]-gnd >= 0 && table[index]-gnd < gndblock )
								excitonic[n][mu][nu]
									= gsl_complex_add(
										excitonic[n][mu][nu],
										gsl_complex_mul( gsl_complex_conjugate(ground[table[index]-gnd]), gsl_complex_mul_real( ground[i-gnd], myweight * phase ) )
								);
						}
					}
				}
				free(ground);
			}
			fclose(fground);
			
			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	excitonic[n][mu][nu] = gsl_complex_div_real( excitonic[n][mu][nu], partition );
			unitary(excitonic[n], excitonic_t2g[n], transform, tNC);
			for(mu=0; mu<nctot; mu++){
				mag[n][mu] /= partition;
				doccu[n][mu] /= partition;
				filling[n][mu][0] /= partition;
				filling[n][mu][1] /= partition;
				sum[0] += mag[n][mu];
				sum[1] += doccu[n][mu];
				sum[2] += filling[n][mu][0] + filling[n][mu][1];
			}
		}
		fprintf(fmag, "%9.6lf\t\t", sum[0]);
		fprintf(fdoc, "%9.6lf\t\t", sum[1]);
		fprintf(ffil, "%9.6lf\t\t", sum[2]);
		fprintf(ffil_t2g, "%9.6lf\t\t", sum[2]);
		for(n=0; n<Ni; n++){
			lsum = 0;
			sum[0] = 0;
			sum[1] = 0;
			sum[2] = 0;
			for(mu=0; mu<nctot; mu++){
				sum[0] += mag[n][mu];
				sum[1] += doccu[n][mu];
				sum[2] += filling[n][mu][0] + filling[n][mu][1];
				lsum += GSL_REAL(excitonic_t2g[n][2*mu][2*mu]) + GSL_REAL(excitonic_t2g[n][2*mu+1][2*mu+1]);
			}
			fprintf(falpha, "%9.6lf\t\t", (2-GSL_REAL(excitonic_t2g[n][0][0])-GSL_REAL(excitonic_t2g[n][1][1]))/(tNC-lsum) );
			fprintf(fmag, "\t%9.6lf\t", sum[0]);
			fprintf(fdoc, "\t%9.6lf\t", sum[1]);
			fprintf(ffil, "\t%9.6lf\t", sum[2]);
			fprintf(ffil_t2g, "\t%9.6lf\t", sum[2]);
			for(mu=0; mu<nctot; mu++) fprintf(fmag, "%9.6lf\t", mag[n][mu]);
			for(mu=0; mu<nctot; mu++) fprintf(fdoc, "%9.6lf\t", doccu[n][mu]);
			for(mu=0; mu<nctot; mu++) fprintf(ffil, "%9.6lf\t%9.6lf\t\t", filling[n][mu][0], filling[n][mu][1]);
			for(mu=0; mu<nctot; mu++) fprintf(ffil_t2g, "%9.6lf\t%9.6lf\t\t", GSL_REAL(excitonic_t2g[n][mu][0]), GSL_REAL(excitonic_t2g[n][mu][1]));
			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
				fprintf(fext, "%9.6lf\t%9.6lf\t\t", GSL_REAL(excitonic[n][mu][nu]), GSL_IMAG(excitonic[n][mu][nu]));
				fprintf(ft2g, "%9.6lf\t%9.6lf\t\t", GSL_REAL(excitonic_t2g[n][mu][nu]), GSL_IMAG(excitonic_t2g[n][mu][nu]));
				if( ic==count && gsl_complex_abs2( excitonic[n][mu][nu] ) > 1e-6 )	final_nonzero = 1;
			}
		}
		fprintf(falpha, "\n");
		fprintf(fmag, "\n");
		fprintf(fdoc, "\n");
		fprintf(ffil, "\n");
		fprintf(ffil_t2g, "\n");
		fprintf(fext, "\n");
		fprintf(ft2g, "\n");
	}
	if(final_nonzero){
		for(n=0; n<Ni; n++)	print_gscmatrixd("excitonic", excitonic[n], tNC, tNC);
		for(n=0; n<Ni; n++)	print_gscmatrixd("excitonic_t2g", excitonic_t2g[n], tNC, tNC);
	}
	freetritensord(ginformation, Ni, SDmax);
	fclose(fmag);
	fclose(fdoc);
	fclose(ffil);
	fclose(ffil_t2g);
	fclose(fext);
	freegsctritensord(excitonic, Ni, tNC);
	freegsctritensord(excitonic_t2g, Ni, tNC);
	for(n=0; n<Ni; n++){
		free( egv[n].egbath );
		freegscmatrixd( egv[n].hybrid, tNC );
	}
	free(table);
	freegscmatrixd(transform, NU);

	return 0;
}

