#include <string.h>
#include <gsl/gsl_sort.h>
#include "lattice.h"
#include "gauleg.c"
extern long seed;
extern int nctot, nb, Nmax, beta, tnb, myrank, size, Nonzero_SOC, AF;
extern double DELTA, MU, Hz, JTD, SOC;
extern gsl_complex zero;

int shift_bath(PNSpara *egv, double delta){
	for(int i=0; i<tnb; i++)	egv->egbath[i] = egv->egbath[i] + delta;
	return 0;
}

int update_diag(double *diag){
	int n, mu;
	for(n=0; n<Ni; n++){
		for(mu=0; mu<tNC; mu++)
			diag[n*tNC+mu] = -MU;
		if( Nonzero_SOC && tNC == 6 )	for(mu=4; mu<6; mu++)
			diag[n*tNC+mu] = -MU+1.5*SOC;
	}
	return 0;
}

double set_chem_pot(int Ne, double U, double J){
	double chem[6] = {
		0,
		9/2.* U - 11.	* J,
		7/2.* U - 9.	* J,
		5/2.* U - 15/2. * J,
		3/2.* U - 7.	* J,
		1/2.* U - 4.	* J,
	};
	if( Ne<1 || Ne>5 ){
		printf("Trivial solution! Ne=%d\n", Ne);
		exit(1);
	}

	return chem[Ne];
}


int init_distortion(double distortion[Ni][tNC]){
	int n, mu;
	double dist_eg[tNC] = { JTD/2, JTD/2, -JTD/2, -JTD/2};
	double moments_eg[tNC] = { 0.5, -0.5, 0.5, -0.5};
	double dist_t2g[tNC] = { JTD, JTD, 0, 0, 0, 0};
	double moments_jeff[tNC] = { 1.5, 0.5, -0.5, -1.5, 0.5, -0.5};
	double moments_t2g[tNC] = { 0.5, -0.5, 0.5, -0.5, 0.5, -0.5 };
	double *mymoment;
	double coeff[4][Ni] = {
		{1, 1, 1, 1, 1, 1, 1, 1},
		{1, -1, 1, -1, 1, -1, 1, -1},
		{1, -1, -1, 1, 1, -1, -1, 1},
		{1, -1, -1, 1, -1, 1, 1, -1},
	};

	if( Nonzero_SOC )	mymoment = moments_jeff;
	else			mymoment = moments_t2g;

	if( tNC == 6 ) for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++){
		distortion[n][mu] = dist_t2g[mu] + 0.5 * Hz * coeff[AF][n] * mymoment[mu];
	}
	else if( tNC == 4 ) for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++){
		distortion[n][mu] = dist_eg[mu] + 0.5 * Hz * moments_eg[mu] * (2 *(n%2) - 1);
	}
	return 0;
}

int init_quadrature_nint(Quad *quadrature, int Nint){
	quadrature->Nint[0] = Nint;
	quadrature->Nint[1] = Nint;
	quadrature->Nint[2] = Nint;


	for(int k=0; k<3; k++){
		quadrature->position[k] = mkvectord( quadrature->Nint[k] );
		quadrature->weight[k] = mkvectord( quadrature->Nint[k] );
		gauleg( -PI, PI, (quadrature->position[k])-1, (quadrature->weight[k])-1, quadrature->Nint[k]);
		printf("axis%d: %d points generated\n", k, quadrature->Nint[k]);
	}
	return 0;
}

int init_quadrature(Quad *quadrature){
	quadrature->Nint[0] = Nintx;
	quadrature->Nint[1] = Ninty;
	quadrature->Nint[2] = Nintz;


	for(int k=0; k<3; k++){
		quadrature->position[k] = mkvectord( quadrature->Nint[k] );
		quadrature->weight[k] = mkvectord( quadrature->Nint[k] );
		gauleg( -PI, PI, (quadrature->position[k])-1, (quadrature->weight[k])-1, quadrature->Nint[k]);
		printf("axis%d: %d points generated\n", k, quadrature->Nint[k]);
	}
	return 0;
}

int free_quadrature(Quad *quadrature){
	for(int i=0; i<3; i++){
		free(quadrature->position[i]);
		free(quadrature->weight[i]);
	}
	return 0;
}

int compute_hopmatrix(Latt *data, gsl_complex **hopmatrix, double k1, double k2, double k3, int line){
	int i, mu, nu;

	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
		hopmatrix[mu][nu] = zero;

	for(i=0; i<line; i++){
		hopmatrix[data[i].mu][data[i].nu]
			= gsl_complex_add(
					hopmatrix[data[i].mu][data[i].nu],
					gsl_complex_mul(
						data[i].tmunu,
						gsl_complex_rect(
							cos( k1*data[i].da1+ k2*data[i].da2 + k3*data[i].da3 ),
							sin( k1*data[i].da1+ k2*data[i].da2 + k3*data[i].da3 )
					)
				)
			);
	}

	return 0;
}//end of compute_hopmatrix




int init_pair_t2g(int **pair_t2g){
	pair_t2g[0][0] = 0; pair_t2g[0][1] = 1;	pair_t2g[0][2] = 1;
	pair_t2g[1][0] = 2; pair_t2g[1][1] = 3;	pair_t2g[1][2] = 1;
	pair_t2g[2][0] = 4; pair_t2g[2][1] = 5;	pair_t2g[2][2] = 1;
	return 0;
}

int init_pair_jj(int **pair_jj){
	int i, j;
	int pair_t2g[tNC/2][3] = {
		{ 0, 3,		3 },
		{ 1, 2, 	1 },
		{ 4, 5, 	1 },
	};
	int pair_eg[tNC/2][3] = {
		{ 0, 1,		1 },
		{ 2, 3, 	-1 },
	};
	if( tNC == 6 ) for(i=0; i<nctot; i++) for(j=0; j<3; j++)
		pair_jj[i][j] = pair_t2g[i][j];
	else if( tNC == 4 ) for(i=0; i<nctot; i++) for(j=0; j<3; j++)
		pair_jj[i][j] = pair_eg[i][j];
	return 0;
}

double find_Ef(Latt *data, Quad *quadrature, int lines, double occupancy){
	double distortion[Ni][tNC];
	double diag[NU];
	init_distortion(distortion);
	update_diag(diag);

	int Nenergy = (1<<18), Nk_brill=1,n;
	int *histo = mkvectori(Nenergy);
	int Omp_num = omp_get_max_threads();
	gsl_complex ***hamil = mkgsctritensord(Omp_num, NU, NU), ***evec = mkgsctritensord(Omp_num, NU, NU);
	double emin = 1e6, emax = -1e6;
	int kx, ky, kz, mu, i;
	for(i=0; i<3; i++)	Nk_brill *= quadrature->Nint[i];
	double **eval = mkmatrixd(Nk_brill, NU);
	for(i=0; i<Nenergy; i++)	histo[i] = 0;

#pragma omp parallel default(none) private(kx,ky,kz,mu,i,n) shared(data,lines,quadrature,eval,evec,hamil, zero, emin, emax, distortion, diag)
	{
		int mythread = omp_get_thread_num(), myindex;
		double mymax=-1e6, mymin=1e6;
#pragma omp for
		for(kx=0; kx<Nintx; kx++) for(ky=0; ky<Ninty; ky++) for(kz=0; kz<Nintz; kz++){
			myindex = kx*quadrature->Nint[1]*quadrature->Nint[2] + ky*quadrature->Nint[2] + kz;
			compute_hopmatrix(data, hamil[mythread], quadrature->position[0][kx], quadrature->position[1][ky], quadrature->position[2][kz], lines);

			for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) 
				GSL_REAL(hamil[mythread][tNC*n+mu][tNC*n+mu]) = GSL_REAL(hamil[mythread][tNC*n+mu][tNC*n+mu]) + distortion[n][mu] + diag[n*tNC+mu];

			eigen_lapack(hamil[mythread], eval[myindex], evec[mythread], NU, 1);
			if( eval[myindex][0]	< mymin )	mymin = eval[myindex][0];
			if( eval[myindex][NU-1]	> mymax )	mymax = eval[myindex][NU-1];
		}

#pragma omp critical
		{
			if( mymin < emin )	emin = mymin;
			if( mymax > emax )	emax = mymax;
		}
	}
	emin -= 1e-6;
	emax += 1e-6;
	double dE = (emax - emin) / Nenergy;
	int ie;
	printf("<find_fermi_emergy>::-----------------------------\nEnergy in (%lf:%lf), dE = %lf with %d grids\n", emin, emax, dE, Nenergy);
	for(i=0; i<Nk_brill; i++) for(mu=0; mu<NU; mu++) {
		ie = (int)(( eval[i][mu] - emin )/dE);
		histo[ie] += 1;
	}

	int accumulated = 0;
	for(i=0; i<Nenergy; i++) {
		accumulated += histo[i];
		if( accumulated > occupancy*Nk_brill*NU )	break;
	}
	double Efermi = emin + i * dE - dE/2.;

	printf("Fermi energy = %lf\n--------------------------------------------------\n", Efermi);

	free(histo);
	freegsctritensord(hamil, Omp_num, NU);
	freegsctritensord(evec, Omp_num, NU);
	freematrixd(eval, Nk_brill);
	return Efermi;
}

Latt* read_lattice(char *name, int *line, double *Efermi){
	int mu, nu, n, m, l;
	double da1, da2, da3;
	double re, im;
	Latt *data;
	FILE *fla;
	if( !test_file(WANNIER) ){
		printf("No Wannier file is found!\n");
		exit(1);
	}
	fla = fopen(name , "r");
	int num = 0;

	fscanf(fla, "%lf", Efermi);

	while( fscanf(fla, "%d %d %d %d %d %lf %lf", &n, &m, &l, &mu, &nu, &re, &im) != EOF ){
		num++;
	}
	fclose(fla);
	*line = num;

	data = (Latt*) malloc( sizeof(Latt)*num );
	printf("%d line\n", num);	fflush(stdout);

	fla = fopen(name , "r");


	fscanf(fla, "%lf", Efermi);
	*Efermi = *Efermi * Hop;
	num = 0;

	while( fscanf(fla, "%d %d %d %d %d %lf %lf", &n, &m, &l, &mu, &nu, &re, &im) != EOF ){
		//dx = n * R[0][0] + m * R[1][0] + l * R[2][0];
		//dy = n * R[0][1] + m * R[1][1] + l * R[2][1];
		//dz = n * R[0][2] + m * R[1][2] + l * R[2][2];
		da1 = n;
		da2 = m;
		da3 = l;
		//dy = n+m;
		//dz = m+l;
		mu--; nu--;

		data[num].da1 = da1;
		data[num].da2 = da2;
		data[num].da3 = da3;
		data[num].mu = mu;
		data[num].nu = nu;
		data[num].tmunu = gsl_complex_rect(Hop*re, Hop*im);
		//print_complex(data[num].tmunu);

		num ++;
	}
	fclose(fla);

	return data;
}

void print_p(char *name, double p[]){
	int i;
	printf("\n----------------------- < %s > -----------------------------------------------\n", name);
	for(i=0; i<Np; i++)	printf("p[%2d] = %lf\t", i, p[i]);
	printf("\n");
	printf("--------------------------------------------------------------------------------\n\n");
}

int init_nametag_cluster(char nametag[tNC][16], char basis){
	if( basis == 'j'){
		if( tNC == 6 ){
			strcpy(nametag[0], " 3/2");
			strcpy(nametag[1], " 1/2");
			strcpy(nametag[2], "-1/2");
			strcpy(nametag[3], "-3/2");
			strcpy(nametag[4], " 1/2");
			strcpy(nametag[5], "-1/2");
		}
		else if( tNC == 4 ){
			strcpy(nametag[0], " 1/2");
			strcpy(nametag[1], "-1/2");
			strcpy(nametag[2], " 1/2");
			strcpy(nametag[3], "-1/2");
		}
	}
	else if( basis == 't'){
		if( tNC == 6 ){
			strcpy(nametag[0], "xy_u");
			strcpy(nametag[1], "xy_d");
			strcpy(nametag[2], "yz_u");
			strcpy(nametag[3], "yz_d");
			strcpy(nametag[4], "zx_u");
			strcpy(nametag[5], "zx_d");

		}
		else if ( tNC == 4 ){
			strcpy(nametag[0], "x2-y2 ");
			strcpy(nametag[1], "x2-y2 ");
			strcpy(nametag[2], "3z2-r2");
			strcpy(nametag[3], "3z2-r2");
		}
	}
	else{
		printf("undefined basis name '%c'\n", basis);
		exit(1);
	}
	return 0;
}

void print_egv(char *name, PNSpara egv){
	int sp2, mu, k;
	char nametag_cluster[tNC][16], nametag_bath_spin[Spmax][16];
	if( Nonzero_SOC )	init_nametag_cluster(nametag_cluster, 'j');
	else			init_nametag_cluster(nametag_cluster, 't');
	strcpy(nametag_bath_spin[0], " s");
	strcpy(nametag_bath_spin[1], "-s");
	printf("\n----------------------- < %s > -----------------------------------------------\n", name);
	for(sp2=0; sp2<Spmax; sp2++){
		for(k=0; k<nb; k++)	printf("egbath[%d][%s]\t\t= %lf\t\t\t", k, nametag_bath_spin[sp2], egv.egbath[2*k+sp2]);
		printf("\n");
	}
	for(sp2=0; sp2<Spmax; sp2++){
		for(mu=0; mu<tNC; mu++){
			for(k=0; k<nb; k++)	printf("hybrid[%s][%d,%s]\t= (%9.6lf, %9.6lf) \t", nametag_cluster[mu], k, nametag_bath_spin[sp2], GSL_REAL(egv.hybrid[mu][2*k+sp2]), GSL_IMAG(egv.hybrid[mu][2*k+sp2]));
			printf("\n");
		}
		printf("\n");
	}
	printf("--------------------------------------------------------------------------------\n\n");
}
void fprint_egv(FILE *file, char *name, PNSpara egv){
	int sp1, sp2, i, j;
	fprintf(file, "\n----------------------- < %s > -----------------------------------------------\n", name);
	for(sp1=0; sp1<Spmax; sp1++){
		for(i=0; i<nb; i++)	fprintf(file, "egbath[%d][%i]\t= %lf\t\t\t", sp1, i, egv.egbath[2*i+sp1]);
		fprintf(file, "\n");
	}
	for(sp1=0; sp1<Spmax; sp1++) for(sp2=0; sp2<Spmax; sp2++){
		for(i=0; i<nctot; i++){
			for(j=0; j<nb; j++)	fprintf(file, "hybrid[%d,%d][%d,%d]\t= (%9.6lf, %9.6lf) \t", i, sp1, j, sp2, GSL_REAL(egv.hybrid[2*i+sp1][2*j+sp2]), GSL_IMAG(egv.hybrid[2*i+sp1][2*j+sp2]));
			fprintf(file, "\n");
		}
		fprintf(file, "\n");
	}
	fprintf(file, "--------------------------------------------------------------------------------\n\n");
}

int sort_bath_zeroSOC(PNSpara *egv, int sp){
	size_t p[2*NB];
	int i, j;
	double *egbath = mkvectord(nb);
	gsl_complex **hybrid = mkgscmatrixd(nctot, nb);
	print_egv("before sorting zeroSOC", *egv);
	for(i=0; i<nb; i++){
		egbath[i] = egv->egbath[2*i+sp];
	}

	gsl_sort_index(p, egbath, 1, nb);
	for(i=0; i<nctot; i++) for(j=0; j<nb; j++)	hybrid[i][j] = egv->hybrid[2*i+sp][2*j+sp];

	for(i=0; i<nb; i++){
		egv->egbath[2*i+sp] = egbath[p[i]];
	}
	for(i=0; i<nctot; i++) for(j=0; j<nb; j++)	egv->hybrid[2*i+sp][2*j+sp] = hybrid[i][p[j]];
	print_egv("after sorting zeroSOC", *egv);

	free(egbath);
	freegscmatrixd(hybrid, nctot);

	return 0;
}

int sort_bath(PNSpara *egv){
	size_t p[2*NB];
	int i, j;
	double *egbath = mkvectord(tnb);
	gsl_complex **hybrid = mkgscmatrixd(tNC, tnb);
	print_egv("before sorting", *egv);
	for(i=0; i<tnb; i++){
		egbath[i] = egv->egbath[i];
	}

	gsl_sort_index(p, egbath, 1, tnb);
	for(i=0; i<tNC; i++) for(j=0; j<tnb; j++)	hybrid[i][j] = egv->hybrid[i][j];

	for(i=0; i<tnb; i++){
		egv->egbath[i] = egbath[p[i]];
	}
	for(i=0; i<tNC; i++) for(j=0; j<tnb; j++)	egv->hybrid[i][j] = hybrid[i][p[j]];
	print_egv("after sorting", *egv);

	free(egbath);
	freegscmatrixd(hybrid, tNC);

	return 0;
}

int set_inter_spin_zero(PNSpara *egv){
	for(int mu=0; mu<nctot; mu++) for(int l=0; l<nb; l++){
		egv->hybrid[2*mu][2*l+1] = zero;
		egv->hybrid[2*mu+1][2*l] = zero;
	}
	return 0;
}

int symmetrize_egv(PNSpara *egv, gsl_complex **transform){
	gsl_complex **V_t2g = mkgscmatrixd(tNC, tnb);
	gsl_complex **Tinv = mkgscmatrixd(tNC, tNC);
	int mu, nu, l, k;

	PNSpara tegv;
	tegv.egbath = mkvectord(tnb);
	tegv.hybrid = mkgscmatrixd(tNC, tnb);

	sort_bath(egv);

	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	Tinv[mu][nu] = gsl_complex_conjugate( transform[nu][mu] );

	for(mu=0; mu<tNC; mu++) for(l=0; l<tnb; l++){
		V_t2g[mu][l] = zero;
		for(k=0; k<tNC; k++)
			V_t2g[mu][l] = gsl_complex_add( V_t2g[mu][l], gsl_complex_mul( Tinv[mu][k], egv->hybrid[k][l] ) );
	}
	print_gscmatrixd_compact("V_t2g", V_t2g, tNC, tnb);
	//symmetrize
	for(l=0; l<nb; l++){
		tegv.egbath[2*l] = (egv->egbath[2*l] + egv->egbath[2*l+1])/2.;
		tegv.egbath[2*l+1] = tegv.egbath[2*l];
	}
	for(mu=0; mu<nctot; mu++) for(l=0; l<nb; l++){
		tegv.hybrid[2*mu][2*l]
			= gsl_complex_div_real(
				gsl_complex_add(
					V_t2g[2*mu][2*l],
					gsl_complex_conjugate( V_t2g[2*mu+1][2*l+1] )
				),
				2.0
			);
		tegv.hybrid[2*mu+1][2*l+1] = gsl_complex_conjugate( tegv.hybrid[2*mu][2*l] );

		tegv.hybrid[2*mu+1][2*l]
			= gsl_complex_div_real(
					gsl_complex_sub(
						V_t2g[2*mu+1][2*l],
						gsl_complex_conjugate( V_t2g[2*mu][2*l+1] )
					),
					2.0
				);
		tegv.hybrid[2*mu][2*l+1] = gsl_complex_mul_real( gsl_complex_conjugate( tegv.hybrid[2*mu+1][2*l] ), -1 );
	}
	if( !Nonzero_SOC ) for(mu=0; mu<nctot; mu++) for(l=0; l<nb; l++){
		tegv.hybrid[2*mu][2*l+1] = zero;
		tegv.hybrid[2*mu+1][2*l] = zero;
	}

	//permutation
	int Sunit=3, num_units=NB/Sunit, sp;
	double local_sum_eg;
	gsl_complex local_hyb[tNC];
	for(l=0; l<num_units; l++){
		for(sp=0; sp<2; sp++){
			local_sum_eg = 0;
			for(mu=0; mu<tNC; mu++) local_hyb[mu] = tegv.hybrid[mu][2*(0 + l*Sunit) + sp];	//showing 0 to emphasize the first element in Sunit
			for(k=0; k<Sunit; k++){
				local_sum_eg += tegv.egbath[2*(k + l*Sunit) + sp];
			}
			for(k=0; k<Sunit; k++){
				tegv.egbath[2*(k + l*Sunit) + sp] = local_sum_eg / Sunit;
				for(mu=0; mu<tNC; mu++){
					tegv.hybrid[(2*k+mu)%tNC][2*(k + l*Sunit) + sp] = local_hyb[mu];
				}
			}
		}
	}
	for(l=0; l<tnb; l++) egv->egbath[l] = tegv.egbath[l];
	for(mu=0; mu<tNC; mu++) for(l=0; l<tnb; l++)	V_t2g[mu][l] = tegv.hybrid[mu][l];
	print_gscmatrixd_compact("symmetrized V_t2g", V_t2g, tNC, tnb);

	for(mu=0; mu<tNC; mu++) for(l=0; l<tnb; l++){
		egv->hybrid[mu][l] = zero;
		for(k=0; k<tNC; k++)
			egv->hybrid[mu][l] = gsl_complex_add( egv->hybrid[mu][l], gsl_complex_mul( transform[mu][k], V_t2g[k][l] ) );
	}
	print_egv("after symmetrization", *egv);

	for(l=0; l<tnb; l++)	egv->egbath[l] = tegv.egbath[l];

	free(tegv.egbath);
	freegscmatrixd(tegv.hybrid, tNC);
	freegscmatrixd(V_t2g, tNC);
	freegscmatrixd(Tinv, tNC);
	return 0;
}

int convert_egv_to_p_zeroSOC(PNSpara egv, double p[], int sp){
	int i, j;

	for(i=0; i<nb; i++)				p[i] = egv.egbath[2*i+sp];
	for(i=0; i<nctot; i++) for(j=0; j<nb; j++){
		p[nb + nb*i + j] = GSL_REAL( egv.hybrid[2*i+sp][2*j+sp] );
		p[nb + nb*nctot + nb*i +j]  = GSL_IMAG( egv.hybrid[2*i+sp][2*j+sp] );
	}

	return 0;
}

int convert_p_to_egv_zeroSOC(PNSpara egv, double p[], int sp){
	int i, j;
	for(i=0; i<nb; i++)				egv.egbath[2*i+sp] = p[i];
	for(i=0; i<nctot; i++) for(j=0; j<nb; j++){
		egv.hybrid[2*i+sp][2*j+sp] =
			gsl_complex_rect(
					p[nb + nb*i + j],
					p[nb + nb*nctot + nb*i +j] 
					);
	}
	return 0;
}

int convert_egv_to_p(PNSpara egv, double p[]){
	int i, j;

	for(i=0; i<tnb; i++)		p[i] = egv.egbath[i];
	for(i=0; i<tNC; i++) for(j=0; j<tnb; j++){
		p[tnb + 2*tnb*i +2*j] = GSL_REAL( egv.hybrid[i][j] );
		p[tnb + 2*tnb*i +2*j+1] = GSL_IMAG( egv.hybrid[i][j] );
	}

	//print_egv("convert_egv_to_p", egv);
	//print_p("convert_egv_to_p", p);
	return 0;
}

int convert_p_to_egv(PNSpara egv, double p[]){
	int i, j;
	for(i=0; i<tnb; i++)				egv.egbath[i] = p[i];
	for(i=0; i<tNC; i++) for(j=0; j<tnb; j++){
		egv.hybrid[i][j] =
			gsl_complex_rect(
					p[tnb + 2*tnb*i +2*j ],
					p[tnb + 2*tnb*i +2*j + 1 ]
					);
	}
	//print_egv("convert_p_to_egv", egv);
	return 0;
}

/*
int convert_egv_to_p(PNSpara egv, double p[]){
	int i, j;

	for(i=0; i<tnb; i++)		p[i] = egv.egbath[i];
	for(i=0; i<tNC; i++) for(j=0; j<tnb; j++){
		p[tnb + tnb*i + j] = GSL_REAL( egv.hybrid[i][j] );
		p[tnb + tnb*tNC + tnb*i +j] = GSL_IMAG( egv.hybrid[i][j] );
	}

	//print_egv("convert_egv_to_p", egv);
	//print_p("convert_egv_to_p", p);
	return 0;
}

int convert_p_to_egv(PNSpara egv, double p[]){
	int i, j;
	for(i=0; i<tnb; i++)				egv.egbath[i] = p[i];
	for(i=0; i<tNC; i++) for(j=0; j<tnb; j++){
		egv.hybrid[i][j] =
			gsl_complex_rect(
					p[tnb + tnb*i + j],
					p[tnb + tnb*tNC + tnb*i +j] 
					);
	}
	//print_egv("convert_p_to_egv", egv);
	return 0;
}



int convert_egv_to_p(PNSpara egv, double p[]){
	int j, Neg;
	Neg = tnb/2;

	for(j=0; j<Neg; j++)		p[j] = egv.egbath[2*j];
	for(j=0; j<tnb; j++){
		p[Neg +2*j]	= GSL_REAL( egv.hybrid[0][j] );
		p[Neg +2*j+1]	= GSL_IMAG( egv.hybrid[0][j] );
	}

	//print_egv("convert_egv_to_p", egv);
	//print_p("convert_egv_to_p", p);
	return 0;
}

int convert_p_to_egv(PNSpara egv, double p[]){
	int j, Neg;
	int fa[2] = {1,-1};
	Neg = tnb/2;
	for(j=0; j<Neg; j++){
		egv.egbath[2*j]		= p[j];
		egv.egbath[2*j+1]	= p[j];
	}
	for(j=0; j<tnb; j++){
		egv.hybrid[0][j] =
			gsl_complex_rect(
					p[Neg + 2*j ],
					p[Neg + 2*j + 1 ]
					);

		egv.hybrid[1][j +1-2*(j%2)] =
			gsl_complex_rect(
					 p[Neg + 2*j ]		* fa[j%2],
					-p[Neg + 2*j + 1]	* fa[j%2]
					);
	}
	//print_egv("convert_p_to_egv", egv);
	return 0;
}
*/

int init_sublattice(int *sublattice){
	int mu;
	for(mu=0; mu<nctot; mu++)	sublattice[mu] = mu%2;
	return 0;
}

int init_del(double **delx){
	int mu, nu;
	for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++)
		delx[mu][nu] = nu-mu;
	return 0;
}


int update_chem(gsl_complex **matrix){
	int mu;
	double diag[NU];
	update_diag(diag);
	for(mu=0; mu<NU; mu++)	matrix[mu][mu] = gsl_complex_rect( GSL_REAL(matrix[mu][mu]) + diag[mu], 0 );
	return 0;
}
int confirm_hopmatrix(gsl_complex **matrix, gsl_complex *hopmatrix, Quad *quadrature){
	int mu, nu, kx, ky, kz;

	gsl_complex **sum;

	if( myrank == 0 )	printf("<:: computing hopping matrix ::>-------------------------\n");
	sum		= mkgscmatrixd( NU, NU );
	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
		sum[mu][nu] = zero;
	}
	for(kx=0; kx<Nintx; kx++) for(ky=0; ky<Ninty; ky++) for(kz=0; kz<Nintz; kz++) {
		for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
			sum[mu][nu]
				= gsl_complex_add(
						sum[mu][nu], 
						gsl_complex_mul_real(hopmatrix(kx,ky,kz,mu,nu), (quadrature->weight[0])[kx]*(quadrature->weight[1])[ky]*(quadrature->weight[2])[kz])
				);
		}
	}
	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
		sum[mu][nu] = gsl_complex_div_real( sum[mu][nu], Area );
		//sum[mu][nu] = gsl_complex_div_real( gsl_complex_add(sum[mu][nu], gsl_complex_conjugate(sum[nu][mu])), 2.*Factor );
	}
	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
		matrix[mu][nu] = sum[mu][nu];
	}
	int nonherm=0;
	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
		if( gsl_complex_abs2( gsl_complex_sub( matrix[mu][nu], gsl_complex_conjugate(matrix[nu][mu]) ) ) > 1e-6 ){
			printf("confirm_hopmatrix:: non-Hermitian (%d,%d)\n", mu, nu);
			nonherm++;
		}
	}
	if( myrank == 0 ){
		printf("-----------------------------------------------------\n");
		print_gscmatrixd("Ecluster: confirm_Ecluster", matrix, NU, NU);
		printf("-----------------------------------------------------\n");
		if( nonherm )	exit(1);
	}
	freegscmatrixd(sum, NU);
	//free(data);
	return 0;
}


int init_Ecluster(gsl_complex **matrix, gsl_complex *hopmatrix, Quad *quadrature){
	int mu, nu;
	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
		matrix[mu][nu] = zero;

	confirm_hopmatrix(matrix, hopmatrix, quadrature);

	return 0;
}

int save_Ecluster(char *para, gsl_complex **matrix){
	int mu, nu;
	FILE *fe;

	fe = fopen(para, "w");
	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
		fprintf(fe, "%d\t%d\t%.18lf\t%.18lf\n", mu, nu, GSL_REAL(matrix[mu][nu]), GSL_IMAG(matrix[mu][nu]));
	fclose(fe);

	return 0;
}

int load_Ecluster(char *para, gsl_complex **matrix){
	int mu, nu;
	double re, im;
	FILE *fe;

	fe = fopen(para, "r");
	nofile(fe, para);
	while( fscanf(fe, "%d %d %lf %lf", &mu, &nu, &re, &im) != EOF ){
		matrix[mu][nu] = gsl_complex_rect(re, im);
	}
	fclose(fe);
	return 0;
}

int init_transform(gsl_complex **transform){
	gsl_complex trans[6][6] = {
		{ zero,					zero,				gsl_complex_rect(-1/sqrt(2), 0),	zero,					gsl_complex_rect(0, -1/sqrt(2)),	zero	},
		{ gsl_complex_rect(sqrt(2./3.), 0),	zero,				zero,					gsl_complex_rect(-1/sqrt(6), 0),	zero,					gsl_complex_rect(0, -1/sqrt(6)) },
		{ zero,					gsl_complex_rect(sqrt(2./3.),0),gsl_complex_rect(1/sqrt(6), 0),		zero,					gsl_complex_rect(0,-1/sqrt(6)),		zero	},
		{ zero,					zero,				zero,					gsl_complex_rect(1/sqrt(2), 0),		zero,					gsl_complex_rect(0, -1/sqrt(2))	},
		{ gsl_complex_rect(-1/sqrt(3), 0),	zero,				zero,					gsl_complex_rect(-1/sqrt(3), 0),	zero,					gsl_complex_rect(0, -1/sqrt(3))	},
		{ zero,					gsl_complex_rect(1/sqrt(3),0),	gsl_complex_rect(-1/sqrt(3), 0),	zero,					gsl_complex_rect(0, 1/sqrt(3)),		zero				},
	};
	int mu, nu, n;
	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
		transform[mu][nu] = zero;
	}
	if( tNC == 6 ) for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)
		transform[tNC*n+mu][tNC*n+nu] = gsl_complex_conjugate( trans[mu][nu] );

	else if( tNC == 4 ) for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) 
		transform[tNC*n+mu][tNC*n+mu] = gsl_complex_rect(1, 0);
	if( myrank == 0 )	print_gscmatrixd("init_transform", transform, NU, NU);
	return 0;
}

gsl_complex ***init_pauli(){
	gsl_complex pi[3][2][2] = {
		{
			{ zero,				gsl_complex_rect(1, 0) },
			{ gsl_complex_rect(1, 0),	zero			},
		},
		{
			{ zero,				gsl_complex_rect(0, -1) },
			{ gsl_complex_rect(0, 1),	zero			},
		},
		{
			{ gsl_complex_rect(1, 0),	zero			},
			{ zero,				gsl_complex_rect(-1, 0) },
		}
	};
	gsl_complex ***pauli = mkgsctritensord(3, 2, 2);
	for(int n=0; n<3; n++){
		for(int i=0; i<2; i++) for(int j=0; j<2; j++)	pauli[n][i][j] = pi[n][i][j];
		print_gscmatrixd("pauli", pauli[n], 2, 2);
	}
	return pauli;
}

int build_local_to_global(gsl_complex **rot, char *fname){
	FILE *fp;
	double rot_33[3][3] = {{0}};
	double u[3], sum, theta;
	int shift, i, j, id, sp1, sp2, n;
	gsl_complex rot_22[2][2];

	gsl_complex ***pauli = init_pauli();
	
	fp = fopen(fname, "r");

	//cyclic shift 0;1;2; = x,y; y,z; z,x;
	fscanf(fp, "%d", &shift);
	for(i=0; i<2; i++) for(j=0; j<3; j++)
		fscanf(fp, "%lf", &rot_33[(i+shift)%3][j]);
	fclose(fp);

	for(i=0; i<2; i++){
		id = (i+shift)%3;
		sum = 0;
		for(j=0; j<3; j++)	sum += rot_33[id][j] * rot_33[id][j];
		for(j=0; j<3; j++)	rot_33[id][j] /= sqrt(sum);
	}
	cross_product(rot_33[shift], rot_33[(shift+1)%3], rot_33[(shift+2)%3]);

	for(i=0; i<3; i++)	u[i] = rot_33[(i+2)%3][(i+1)%3] - rot_33[(i+1)%3][(i+2)%3];
	print_vectord("rotation vector", u, 3);
	sum = 0;
	for(i=0; i<3; i++)	sum += u[i] * u[i];
	sum = sqrt(sum);
	for(i=0; i<3; i++)	u[i] /= sum;
	theta = asin(sum/2.);
	printf("|u| = %lf, theta = %lf\n", sum, theta);

	for(i=0; i<2; i++) for(j=0; j<2; j++)	rot_22[i][j] = zero;
	for(i=0; i<2; i++)			rot_22[i][i] = gsl_complex_rect( cos(theta/2), 0 );
	for(i=0; i<2; i++) for(j=0; j<2; j++) for(n=0; n<3; n++){
		rot_22[i][j]
			= gsl_complex_add(
				rot_22[i][j],
				gsl_complex_mul(
					gsl_complex_rect( 0, sin(theta/2) ),
					gsl_complex_mul_real( pauli[n][i][j], u[n] )
				)
			);
	}
	printf("<:: rot_22 ::> ( 2 x 2 )------------");
	for(i=0; i<2; i++){
		for(j=0; j<2; j++){
			printf("%lf, %lf\t", GSL_REAL( rot_22[i][j] ), GSL_REAL( rot_22[i][j] ) );
		}
		printf("\n");
	}


	for(i=0; i<NU; i++) for(j=0; j<NU; j++)	rot[i][j] = zero;
	for(n=0; n<Ni; n++) for(i=0; i<3; i++) for(j=0; j<3; j++) for(sp1=0; sp1<2; sp1++) for(sp2=0; sp2<2; sp2++){//transposed (inversed), i=row index (supposed to be a column index, not row)
		rot[n*tNC + 2*i+sp1][n*tNC + 2*j+sp2] = gsl_complex_mul_real( rot_22[sp2][sp1], rot_33[i][j] );
	}
	print_gscmatrixd("rot_all", rot, NU, NU);

	check_unitary(rot, NU);
	freegsctritensord(pauli, 3, 2);
	return 0;
}


