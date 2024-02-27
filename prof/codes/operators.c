#include "operators.h"

int init_operator_angular_z(gsl_complex ***sz, int **ipair, int print){
	int mu, k, l;
	gsl_complex zero = gsl_complex_rect(0, 0);
	for(mu=0; mu<3; mu++) for(k=0; k<6; k++) for(l=0; l<6; l++)	sz[mu][k][l] = zero;
	for(mu=0; mu<3; mu++){
		sz[mu][ipair[mu][0]][ipair[mu][0]] = gsl_complex_rect( ipair[mu][2] * SPFACTOR, 0 );
		sz[mu][ipair[mu][1]][ipair[mu][1]] = gsl_complex_rect(-ipair[mu][2] * SPFACTOR, 0 );
	}
	if( print ){
		char para[1024];
		for(mu=0; mu<3; mu++){
			sprintf(para, "sz, mu = %d", mu);
			print_gscmatrixd(para, sz[mu], 6, 6);
		}
	}
	return 0;
}

int test_angular_operators(gsl_complex **transform){
	int error = 0;
	gsl_complex **Lplus = mkgscmatrixd(NU, NU), **Lminus = mkgscmatrixd(NU, NU), **Lz = mkgscmatrixd(NU, NU);
	gsl_complex **Jplus = mkgscmatrixd(NU, NU), **Jminus = mkgscmatrixd(NU, NU), **Jz = mkgscmatrixd(NU, NU);
	gsl_complex **Splus = mkgscmatrixd(NU, NU), **Sminus = mkgscmatrixd(NU, NU), **Sz = mkgscmatrixd(NU, NU);

	init_operator_S_jj(Splus, Sminus, Sz, transform, 1);
	init_operator_L_jj(Lplus, Lminus, Lz, transform, 1);
	init_operator_J_jj(Jplus, Jminus, Jz, 1);

	int i, j;
	gsl_complex difference;
	for(i=0; i<6; i++) for(j=0; j<6; j++){
		difference = gsl_complex_sub( Jplus[i][j], gsl_complex_add(Lplus[i][j], Splus[i][j]) );
		difference = gsl_complex_sub( Jminus[i][j], gsl_complex_add(Lminus[i][j], Sminus[i][j]) );
		difference = gsl_complex_sub( Jz[i][j], gsl_complex_add(Lz[i][j], Sz[i][j]) );
		if( gsl_complex_abs2(difference) > 1e-6 )	error++;
	}

	if( error )	printf("<test_angular_operators>:: error = %d\n", error);

	freegscmatrixd(Lplus, NU);
	freegscmatrixd(Lminus, NU);
	freegscmatrixd(Lz, NU);

	freegscmatrixd(Jplus, NU);
	freegscmatrixd(Jminus, NU);
	freegscmatrixd(Jz, NU);

	freegscmatrixd(Splus, NU);
	freegscmatrixd(Sminus, NU);
	freegscmatrixd(Sz, NU);

	return error;
}

double coeff_ladder(double j, double m, int sign){
	if( (j - sign*m) * (j + sign*m + 1 ) > 0 )
		return sqrt( (j - sign*m) * (j + sign*m + 1 ) );
	else
		return 0;
}

int init_pauli(gsl_complex ***pauli, int print){
	gsl_complex zero = gsl_complex_rect(0, 0);
	gsl_complex cone = gsl_complex_rect(1, 0);
	gsl_complex sx[2][2] = { {zero, cone}, {cone, zero} };
	gsl_complex sy[2][2] = { {zero, gsl_complex_rect(0, -1)}, {gsl_complex_rect(0, 1), zero} };
	gsl_complex sz[2][2] = { {cone, zero}, {zero, gsl_complex_rect(-1, 0)} };

	for(int i=0; i<2; i++) for(int j=0; j<2; j++){
		pauli[0][i][j] = sx[i][j];
		pauli[1][i][j] = sy[i][j];
		pauli[2][i][j] = sz[i][j];
	}
	char para[1024];
	if( print ) for(int i=0; i<3; i++){
		sprintf(para, "pauli %d", i);
		print_gscmatrixd(para, pauli[i], 2, 2);
	}

	return 0;
}

int init_operator_S_jj(gsl_complex **Splus, gsl_complex **Sminus, gsl_complex **Sz, gsl_complex **transform, int print){
	gsl_complex zero = gsl_complex_rect(0, 0);
	gsl_complex ***pauli = mkgsctritensord(3, 2, 2), cone = gsl_complex_rect(1, 0);
	init_pauli(pauli, print);
	gsl_complex unit33[3][3] = {
		{ cone, zero, zero },
		{ zero, cone, zero},
		{ zero, zero, cone},
	};

	int i, j, s1, s2;

	for(i=0; i<3; i++) for(j=0; j<3; j++) for(s1=0; s1<2; s1++) for(s2=0; s2<2; s2++){
		Splus[2*i+s1][2*j+s2]	= gsl_complex_mul( unit33[i][j] , gsl_complex_div_real(gsl_complex_add(pauli[0][s1][s2], gsl_complex_mul_imag(pauli[1][s1][s2], 1)), 2.) );
		Sminus[2*i+s1][2*j+s2]	= gsl_complex_mul( unit33[i][j] , gsl_complex_div_real(gsl_complex_sub(pauli[0][s1][s2], gsl_complex_mul_imag(pauli[1][s1][s2], 1)), 2.) );
		Sz[2*i+s1][2*j+s2]	= gsl_complex_mul( unit33[i][j] , gsl_complex_div_real(pauli[2][s1][s2], 2.) );
	}

	if( print ){
		print_gscmatrixd("Splus", Splus, 6, 6);
		print_gscmatrixd("Sminus", Sminus, 6, 6);
		print_gscmatrixd("Sz", Sz, 6, 6);
	}

	iunitary(Splus, Splus, transform, 6);
	iunitary(Sminus, Sminus, transform, 6);
	iunitary(Sz, Sz, transform, 6);

	if( print ){
		print_gscmatrixd("Splus_jj", Splus, 6, 6);
		print_gscmatrixd("Sminus_jj", Sminus, 6, 6);
		print_gscmatrixd("Sz_jj", Sz, 6, 6);
	}

	return 0;
}

int init_operator_J_jj(gsl_complex **Jplus, gsl_complex **Jminus, gsl_complex **Jz, int print){
	double jz_diag[6] = { 1.5, 0.5, -0.5, -1.5, 0.5, -0.5 };
	double j_mag[6] = { 1.5, 1.5, 1.5, 1.5, 0.5, 0.5 };

	int i, j;

	for(i=0; i<6; i++) for(j=0; j<6; j++) {
		Jplus[i][j]	= gsl_complex_rect( (int)( fabs( j_mag[i] - j_mag[j] )<1e-6 ) * (int)( fabs(jz_diag[i]-jz_diag[j]-1)<1e-6 ) * coeff_ladder(j_mag[j], jz_diag[j], 1), 0 );
		Jminus[i][j]	= gsl_complex_rect( (int)( fabs( j_mag[i] - j_mag[j] )<1e-6 ) * (int)( fabs(jz_diag[i]-jz_diag[j]+1)<1e-6 ) * coeff_ladder(j_mag[j], jz_diag[j], -1), 0 );
		Jz[i][j] = gsl_complex_rect( (int)(i==j) * jz_diag[i], 0 );
	}

	if( print ){
		print_gscmatrixd("Jplus", Jplus, 6, 6);
		print_gscmatrixd("Jminus", Jminus, 6, 6);
		print_gscmatrixd("Jz", Jz, 6, 6);
	}

	return 0;
}

int init_operator_L_jj(gsl_complex **Lplus, gsl_complex **Lminus, gsl_complex **Lz, gsl_complex **transform, int print){
	gsl_complex zero = gsl_complex_rect(0, 0);
	double tbtunit[2][2] = {{-1, 0}, {0, -1}};
	gsl_complex lz_perm[3][3] = {
		{ zero, zero, gsl_complex_rect(0, -1) },
		{ zero, zero, zero},
		{ gsl_complex_rect(0,1), zero, zero }
	};

	gsl_complex lp_perm[3][3] = {
		{ zero, gsl_complex_rect(0, 1), zero },
		{ gsl_complex_rect(0, -1), zero, gsl_complex_rect(-1,0)},
		{ zero, gsl_complex_rect(1,0), zero }
	};
	gsl_complex lz[3][3], lp[3][3];
	int permutation[3] = {1, 2, 0};
	int i, j, s1, s2;

	for(i=0; i<3; i++) for(j=0; j<3; j++) {
		lz[i][j] = lz_perm[permutation[i]][permutation[j]];
		lp[i][j] = lp_perm[permutation[i]][permutation[j]];
	}
	for(i=0; i<3; i++) for(j=0; j<3; j++) for(s1=0; s1<2; s1++) for(s2=0; s2<2; s2++)
		Lz[2*i+s1][2*j+s2] = gsl_complex_mul_real( lz[i][j] , tbtunit[s1][s2] );

	for(i=0; i<3; i++) for(j=0; j<3; j++) for(s1=0; s1<2; s1++) for(s2=0; s2<2; s2++)
		Lplus[2*i+s1][2*j+s2] = gsl_complex_mul_real( lp[i][j] , tbtunit[s1][s2] );

	for(i=0; i<3; i++) for(j=0; j<3; j++) for(s1=0; s1<2; s1++) for(s2=0; s2<2; s2++)
		Lminus[2*i+s1][2*j+s2] = gsl_complex_mul_real( gsl_complex_conjugate(lp[j][i]) , tbtunit[s1][s2] );

	if( print ){
		print_gscmatrixd("Lplus", Lplus, 6, 6);
		print_gscmatrixd("Lminus", Lminus, 6, 6);
		print_gscmatrixd("Lz", Lz, 6, 6);
	}

	iunitary(Lplus, Lplus, transform, 6);
	iunitary(Lminus, Lminus, transform, 6);
	iunitary(Lz, Lz, transform, 6);

	if( print ){
		print_gscmatrixd("Lplus_jj", Lplus, 6, 6);
		print_gscmatrixd("Lminus_jj", Lminus, 6, 6);
		print_gscmatrixd("Lz_jj", Lz, 6, 6);
	}

	return 0;
}


