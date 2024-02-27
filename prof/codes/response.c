#include "response.h"
int free_table(Glist_res *table, int tmax){
	for(int i=0; i<tmax; i++){
		freegsctritensord(table[i].mat, table[i].num_op, table[i].dim_op);
		free(table[i].aibi);
	}

	return 0;
}


int init_table_response(Glist_res *table, int tmax, gsl_complex ***matrix_operator, int num_op){
	int i, j, direc, p, q, k, l, dim_op = tNC;

	i=0;
	if( num_op == 2 ){
		for(direc=-1; direc<2; direc+=2){
			for(p=0; p<nctot; p++) for(q=p; q<nctot; q++){//orbital
				table[i].mat = mkgsctritensord(num_op, dim_op, dim_op);

				table[i].direction = direc;

				table[i].ipair[0] = p;
				for(k=0; k<6; k++) for(l=0; l<6; l++) table[i].mat[0][k][l] = matrix_operator[p][k][l];

				table[i].ipair[1] = q;
				for(k=0; k<6; k++) for(l=0; l<6; l++) table[i].mat[1][k][l] = matrix_operator[q][k][l];

				table[i].type = 'k';
				i++;
			}

			for(p=0; p<nctot; p++) for(q=p+1; q<nctot; q++){
				table[i].mat = mkgsctritensord(num_op, dim_op, dim_op);

				table[i].direction = direc;

				table[i].ipair[0] = p;
				for(k=0; k<6; k++) for(l=0; l<6; l++) table[i].mat[0][k][l] = matrix_operator[p][k][l];

				table[i].ipair[1] = q;
				for(k=0; k<6; k++) for(l=0; l<6; l++) table[i].mat[1][k][l] = gsl_complex_mul_imag(matrix_operator[q][k][l], -direc);

				table[i].type = 't';
				i++;
			}
		}
		if( i != tmax ){
			printf("Error in init_table_response!\n");	exit(1);
		}
		for(j=0; j<tmax; j++){
			table[j].dim_op = dim_op;
			table[j].num_op = num_op;
			table[j].aibi = mkvectord(SDmax*2*Upper);
			for(k=0; k<SDmax*2*Upper; k++)	table[j].aibi[k] = 0;
		}
	}
	else	my_error("No rule has been defined");

	return 0;
}

int print_table_response(Glist_res *table, int tmax){
	int i, p, num_op = table[0].num_op, k, l;
	printf("i :  direc ip1 ip2      ");
	for(p=0; p<num_op; p++)
		printf("O%d: mu%d nu%d sp%d  coeff   ", p, p, p, p);
	printf("type\n");
	for(i=0; i<tmax; i++){
		printf("i%2d: %c %2d < %d, %d >:\t",  i, table[i].type, table[i].direction, table[i].ipair[0], table[i].ipair[1]);
		for(p=0; p<num_op; p++){
			if( p )	printf("\t\t\t");
			for(k=0; k<6; k++) for(l=0; l<6; l++)	if(gsl_complex_abs2(table[i].mat[p][k][l]) > 1e-6){
				printf("[%d]%2d,%2d (%5.2lf, %5.2lf)    ",
						p, k, l, GSL_REAL(table[i].mat[p][k][l]), GSL_IMAG(table[i].mat[p][k][l]));
			}
			printf("\n");
		}
	}
	return 0;
}

int free_table_response(Glist_res *table, int tmax){
	int i;
	for(i=0; i<tmax; i++)
		free(table[i].aibi);
	return 0;
}

int apply_four(typebasis **basis, sectorinfo *ginfo, int mu, int nu, int k, int l, gsl_complex *p0, gsl_complex coeff_TT, gsl_complex *ground){
	int i, phase, opern, refer, block, index, *table;
	int o1 = mu/2, sp1 = mu%2, o2 = nu/2, sp2 = nu%2, o3=k/2, sp3=k%2, o4=l/2, sp4=l%2;
	typebasis grstate[2];

	if( mu==nu || k==l )	return 0;

	table = mkvectori( Powns2 );
	for(i=0; i<Powns2; i++){
		index = (basis[i][0]<<Ns) + basis[i][1];
		table[index] = i;
	}

	refer = ginfo->refer;
	block = ginfo->block;

	gsl_complex *plocal = mkgscvectord( block );
	initgscvectord(plocal, block);

	for(i=refer; i<refer+block; i++){
		if( One(basis[i][sp3], o3) && One(basis[i][sp4], o4) ){
			grstate[0] = basis[i][0];	grstate[1] = basis[i][1];
			phase  = permu(grstate, o4, sp4);	grstate[sp4] = Turnoff(grstate[sp4], o4);
			phase *= permu(grstate, o3, sp3);	grstate[sp3] = Turnoff(grstate[sp3], o3);
			if( Zero(grstate[sp1], o1) && Zero(grstate[sp2], o2) ){
				phase *= permu(grstate, o2, sp2);	grstate[sp2] = Turnon(grstate[sp2], o2);
				phase *= permu(grstate, o1, sp1);	grstate[sp1] = Turnon(grstate[sp1], o1);

				index = (grstate[0]<<Ns) + grstate[1];
				opern = (int)table[index] - refer;
				if( opern < 0 || opern >= block ){
					printf("opern %d is not in [%d:%d)\n", opern, 0, block);
					exit(1);
				}
				plocal[opern] = gsl_complex_add( plocal[opern], gsl_complex_mul_real( ground[i-refer], phase) );
			}
		}
	}

	for(i=0; i<block; i++)	p0[i] = gsl_complex_add(p0[i], gsl_complex_mul(plocal[i], coeff_TT) );

	/*
	double sum=0;
	for(i=0; i<block; i++)	sum += gsl_complex_abs2(plocal[i]);
	printf("apply_hop, %d %d , o1 %d, sp1 %d, o2 %d, sp2 %d, %lf, coeff = (%lf, %lf)\n", mu, nu, o1, sp1, o2, sp2, sum, GSL_REAL(coeff_TT), GSL_IMAG(coeff_TT));
	*/

	free(table);
	free(plocal);

	return 0;
}

int apply_hop(typebasis **basis, sectorinfo *ginfo, int mu, int nu, gsl_complex *p0, gsl_complex coeff_TT, gsl_complex *ground){
	int i, phase, opern, refer, block, index, *table;
	int o1 = mu/2, sp1 = mu%2, o2 = nu/2, sp2 = nu%2;
	typebasis grstate[2];

	table = mkvectori( Powns2 );
	for(i=0; i<Powns2; i++){
		index = (basis[i][0]<<Ns) + basis[i][1];
		table[index] = i;
	}

	refer = ginfo->refer;
	block = ginfo->block;

	gsl_complex *plocal = mkgscvectord( block );
	initgscvectord(plocal, block);

	for(i=refer; i<refer+block; i++){
		if( One(basis[i][sp2], o2) ){
			grstate[0] = basis[i][0];	grstate[1] = basis[i][1];
			phase = permu(grstate, o2, sp2);	grstate[sp2] = Turnoff(grstate[sp2], o2);
			if( Zero(grstate[sp1], o1) ){
				phase *= permu(grstate, o1, sp1);	grstate[sp1] = Turnon(grstate[sp1], o1);

				index = (grstate[0]<<Ns) + grstate[1];
				opern = (int)table[index] - refer;
				if( opern < 0 || opern >= block ){
					printf("opern %d is not in [%d:%d)\n", opern, 0, block);
					exit(1);
				}
				plocal[opern] = gsl_complex_add( plocal[opern], gsl_complex_mul_real( ground[i-refer], phase) );
			}
		}
	}

	for(i=0; i<block; i++)	p0[i] = gsl_complex_add(p0[i], gsl_complex_mul(plocal[i], coeff_TT) );

	/*
	double sum=0;
	for(i=0; i<block; i++)	sum += gsl_complex_abs2(plocal[i]);
	printf("apply_hop, %d %d , o1 %d, sp1 %d, o2 %d, sp2 %d, %lf, coeff = (%lf, %lf)\n", mu, nu, o1, sp1, o2, sp2, sum, GSL_REAL(coeff_TT), GSL_IMAG(coeff_TT));
	*/

	free(table);
	free(plocal);

	return 0;
}

int apply_operator(typebasis **basis, sectorinfo ginfo[], int gindex, int direction, gsl_complex *p0, gsl_complex *ground, gsl_complex **operator_matrix){
	int alpha, beta;
	gsl_complex coeff;

	for(alpha=0; alpha<tNC; alpha++) for(beta=0; beta<tNC; beta++){
		coeff = operator_matrix[alpha][beta];
		if( gsl_complex_abs2( coeff ) > Th_hop )
			apply_hop(basis, &ginfo[gindex], alpha, beta, p0, coeff, ground);
	}

	return 0;
}//end of apply_operator

int print_p0inner(int direction, Glist_res *table, double p0inner){
	printf("%2d: %2d\t < %d, %d >:\t ", direction, table->direction, table->ipair[0], table->ipair[1]);
	printf("%c\t%lf\n", table->type, p0inner);
	return 0;
}

int build_p0_operator(typebasis **basis, sectorinfo ginfo[], int gindex, gsl_complex *ground, gsl_complex *p0, double *p0inner, Glist_res *table){
	int i, k, grblock, direction = table->direction, num_op = table->num_op;

	grblock = ginfo[gindex].block;

	initgscvectord(p0, grblock);

	for(i=0; i<num_op; i++){
		apply_operator(basis, ginfo, gindex, direction, p0, ground, table->mat[i]);
	}

	*p0inner = 0.;
	for(k=0; k<grblock; k++)	*p0inner += gsl_complex_abs2( p0[k] );
	if( fabs(*p0inner) > 1e-6 )	for(k=0; k<grblock; k++)	p0[k] = gsl_complex_div_real( p0[k], sqrt(*p0inner) );
	if( (table->ipair[0] == table->ipair[1]) && table->num_op == 2 )	*p0inner /= 4.;

	print_p0inner(direction, table, *p0inner);
	//if(direction == -1)	printf("%c: Op1 = (n%d-n%d) *%d*%.1lf, Op2 = (n%d-n%d) *%d*%.1lf\tp0inner\t= %.10lf\n", table->type, table->mu[0], table->nu[0], table->sp[0], SPFACTOR, table->mu[1], table->nu[1], table->sp[1], SPFACTOR, *p0inner);
	return 0;
}//end of build_p0_operator

int build_spin_correlation(typebasis **basis, sectorinfo ginfo[], int gindex, gsl_complex *ground, PNSpara egv, Glist_res *table, int degen, int n, gsl_complex **transform){
	int p0index, grblock, mapnumber, *column, dimension, *istart, k;
	double *ai, *bi, *Hdia, p0inner, *aibi = (table->aibi) + degen*2*Upper;
	gsl_complex *Hamil, *p0_complex;

	p0index = gindex;
	grblock = ginfo[p0index].block;

	istart	= mkvectori(grblock+1);
	mapnumber = ginfo[p0index].mapnumber;
	Hdia	= mkvectord(grblock);
	Hamil	= mkgscvectord(mapnumber);
	column	= mkvectori(mapnumber);
	read_map(n, basis, istart, Hamil, column, ginfo, p0index, egv);
	compute_Hdia(n, Hdia, ginfo, p0index, basis, egv);

	ai = mkvectord(grblock+1);	bi = mkvectord(grblock+1);

	p0_complex = mkgscvectord( grblock );
	dimension = ( grblock < Upper ? grblock : Upper );

	build_p0_operator(basis, ginfo, gindex, ground, p0_complex, &p0inner, table);
	lanczos_vector_complex(istart, Hdia, Hamil, column, &dimension, p0_complex, ai, bi, grblock, dimension, 0);
	for(k=0; k<dimension; k++){
		aibi[k] = ai[k];
		aibi[k+Upper] = bi[k];
	}
	aibi[Upper] = p0inner;
	table->dimension = dimension;

	free(istart);	free(Hdia);	free(Hamil);	free(column);
	free(ai);	free(bi);
	free(p0_complex);

	return 0;
}//end of build_spin_correlation



