#define Repeat_rank 30

#include "header/main.h"
#include "header/gsl.h"
#include "header/malloc.h"
#include "header_hyb.h"
//#include "header/legendre_poly.h"

//#define DEBUG

void near4(struct G_params* p, struct Malloc* mem , struct complex_func* cfunc);
void near123(struct G_params* p, struct Malloc* mem, struct complex_func* cfunc);
void prepare(struct G_params* p, struct Malloc* mem, struct complex_func* cfunc);
void fitting(struct Malloc* mem, struct complex_func* cfunc, int seed);
void bundle(struct G_params* p, struct Malloc* mem , struct complex_func* cfunc);
void save_cfunc(struct complex_func* cfunc);
int CheckError();
int find_seed();

int list_idx[3];

int main(int argc, char *argv[]){
	Nb = atoi(argv[1]);
	nearest = atoi(argv[2]); 
	for(int i = 0; i < 3; i++) 
		list_idx[i] = atoi(argv[i+3]);

	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	struct G_params* p = (struct G_params*)malloc(sizeof(struct G_params));		p->W = mdcomplex(Nmax);
	struct Malloc* mem = (struct Malloc*)malloc(sizeof(struct Malloc));
	struct complex_func* cfunc = (struct complex_func*)malloc(sizeof(struct complex_func));
	
	w = mdouble(Nmax);
	wn = mdouble(Nmax); 

	cfunc->Diwn = mdcomplex(Nmax);	
	cfunc->Giwn = mdcomplex(Nmax); 	
	cfunc->Gw = mdcomplex(Nmax);		
	cfunc->Dw = mdcomplex(Nmax);
	cfunc->Gt = mdcomplex(Lmax);		
	cfunc->Dt = mdcomplex(Lmax);
	cfunc->G_coeff = mdcomplex(nquad);		
	cfunc->D_coeff = mdcomplex(nquad);

	p->e_params.near = nearest;
	// malloc for each rank /////////////////////////////////////////////
	mem->SortNum_rank = (Repeat_rank > 10) ? 10 : Repeat_rank;
	mem->Chi_rank = mdouble( Repeat_rank * sizeof(double) );
	mem->Iter_rank = mI64( Repeat_rank * sizeof(long long int) );
	mem->Bath_rank = mdouble( Repeat_rank * 4*Nb * sizeof(double) );
	mem->Index_rank = mint( Repeat_rank * sizeof(int) );
	// malloc to be sorted for each rank ///////////////////////////////
	mem->SortChi_rank = mdouble( mem->SortNum_rank * sizeof(double) );
	mem->SortIter_rank = mI64( mem->SortNum_rank * sizeof(long long int) );
	mem->SortBath_rank = mdouble( mem->SortNum_rank * 4*Nb *  sizeof(double) );

	Gaussian(0, B, Lmax, ti, wi);

	if ( nearest == 4 )
		near4(p, mem, cfunc);
	else 	
		near123(p, mem, cfunc);

	free(mem->Bath_rank);		free(mem->Chi_rank); free(mem->Index_rank); free(mem->Iter_rank);
	free(mem->SortChi_rank);	free(mem->SortBath_rank);	free(mem->SortIter_rank);

	free(w);					free(wn);
	free(p->W);					free(p);

	free(cfunc->Giwn);			free(cfunc->Diwn);
	free(cfunc->Gw);			free(cfunc->Dw);			
	free(cfunc->Gt);			free(cfunc->Dt);
	free(cfunc->G_coeff);		free(cfunc->D_coeff);
	return 0; 
}

void prepare(struct G_params* p, struct Malloc* mem, struct complex_func* cfunc)
{
	fermiNdos(p);	// calculate fermi value & DOS	
	// if 3rd value is 'r',calculate G in real domain, else in imag domain 	
	calculate_GD(w, cfunc->Gw, cfunc->Dw, p, 'r');			
	calculate_GD(wn, cfunc->Giwn, cfunc->Diwn, p, 'i'); 
	FT(cfunc->Giwn, wn, cfunc->Gt, ti, B);
	FT(cfunc->Diwn, wn, cfunc->Dt, ti, B);
	L_coeff(cfunc->Gt, cfunc->G_coeff, ti, wi);
	L_coeff(cfunc->Dt, cfunc->D_coeff, ti, wi);

	mem->peak_num = find_extrema(w, cfunc->Dw, mem->extrema_index);
	mem->half_index = (int*)malloc( 2 * mem->peak_num * sizeof(int) );
	mem->num_extrema = 2*mem->peak_num -1; 
	for(int i=0; i < mem->peak_num; i++)
		FWHM(w, cfunc->Dw, mem->extrema_index, mem->num_extrema, 2*i, mem->half_index);

	pathname(p);
	printf("%s\n", Hamilton);
	
#ifndef DEBUG
	if( nearest == 4)
		sprintf(FILE_h5, "database_test/Nb%d/n%d_%d%d%d.h5", Nb, nearest, list_idx[0], list_idx[1], list_idx[2]);
	else
		sprintf(FILE_h5, "database_test/Nb%d/n%d_%d%d.h5", Nb, nearest, list_idx[0], list_idx[1]);
#else
	sprintf(FILE_h5, "test.h5");
#endif

	save_cfunc( cfunc );
}

void save_cfunc(struct complex_func* cfunc){
	hid_t F_id, G1_id;
	char G1[1024], G2[1024]; 
	sprintf(G1, "/%s", Hamilton);	sprintf(G2, "bath%d", Nb);

	F_id = H5Fopen(FILE_h5, H5F_ACC_RDWR, H5P_DEFAULT); 
	if(F_id < 0)
		F_id = H5Fcreate(FILE_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	G1_id = H5Gopen2(F_id, G1, H5P_DEFAULT);
	if( G1_id < 0)
		G1_id = H5Gcreate2(F_id, G1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	save_data(F_id, G1_id, cfunc->Gw, "Gw", Nmax);
	save_data(F_id, G1_id, cfunc->Dw, "Dw", Nmax);
	save_data(F_id, G1_id, cfunc->Giwn, "Giwn", Nmax);
	save_data(F_id, G1_id, cfunc->Diwn, "Diwn", Nmax);
	save_data(F_id, G1_id, cfunc->G_coeff, "Gcff", Lmax);
	save_data(F_id, G1_id, cfunc->D_coeff, "Dcff", Lmax);
	save_freq(F_id, w, wn);

	H5Gclose(G1_id);
	H5Fclose(F_id);
}

void near123(struct G_params* p, struct Malloc* mem, struct complex_func* cfunc)
{
	int s1 = list_idx[0]-1;	//0~17
	int s2 = list_idx[1]-1;	//0~3
	int s3[] = {2, 7, 12, 16, 20};

	p->e_params.n = 0;
	p->e_params.m = 0;
	p->e_params.h = 0;
	p->e_params.near = 1;

	if( s1 == 0 )
		bundle(p, mem, cfunc);

	for(int n = 2 + s1 ; n < 3 + s1 ; n++) 
	{
		p->e_params.near = 2;
		p->e_params.n = 0.05 * n;

		bundle(p, mem, cfunc);
			
		for(int m = s3[s2]; m < s3[s2+1]; m++) 
		{
			p->e_params.near = 3;
			p->e_params.m = 0.05 * n * 0.05 * m;

			bundle(p, mem, cfunc);
		}
	}
}

void near4(struct G_params* p, struct Malloc* mem, struct complex_func* cfunc)
{
	p->e_params.near = 4;
	int s1 = list_idx[0] - 1;
	int s2 = 3 * (list_idx[1] - 1);
	int s3 = 3 * (list_idx[2] - 1);

	for (int i = 0; i < 3; i++) {
		int n = 2 + s1 + 6*i;
		for (int m = 2 + s2; m < 5 + s2; m++) {
			for (int h = 2 + s3; h < 5 + s3; h++) 
			{
				p->e_params.n = 0.05 * n;
				p->e_params.m = 0.05 * n * 0.05 * m;
				p->e_params.h = 0.05 * n * 0.05 * m * 0.05 * h;

				bundle(p, mem, cfunc);
			}
		}
	}
}

int find_seed()
{
	int seed;
	hid_t file_id, group_id, seed_id;

	file_id = H5Fopen(FILE_h5, H5F_ACC_RDONLY, H5P_DEFAULT);
	group_id = H5Gopen2(file_id, Hamilton, H5P_DEFAULT);
	seed_id = H5Aopen(group_id, "seed", H5P_DEFAULT);
	if( seed_id < 0)
	{
		H5Aclose(seed_id);
		H5Gclose(group_id);
		H5Fclose(file_id);
#ifdef DEBUG
		printf("not yet start\n");
#endif
		return 0;
	}

	H5Aread(seed_id, H5T_NATIVE_INT, &seed);

	H5Aclose(seed_id);
	H5Gclose(group_id);
	H5Fclose(file_id);
#ifdef DEBUG
		printf("latest seed = %d\n", seed);
#endif
	return seed+1;
}

int CheckError()
{
	hid_t file_id, group_id, group2_id, P0_id, P9_id, P0chi_id, P9chi_id, iter_id;
	double chi_m, chi_M;	
	int iter, result;

	char bathname[128]; sprintf(bathname, "bath%d", Nb);
	
	file_id = H5Fopen(FILE_h5, H5F_ACC_RDONLY, H5P_DEFAULT);
	
	group_id = H5Gopen2(file_id, Hamilton, H5P_DEFAULT);

	group2_id = H5Gopen2(group_id, bathname , H5P_DEFAULT);
	if( group2_id <  0 ){
		printf("bath is not exist yet\n");
		H5Gclose(group2_id);
		H5Gclose(group_id);
		H5Fclose(file_id);
		return 0;
	}
	iter_id = H5Aopen(group2_id, "Tot_Iter", H5P_DEFAULT);
	H5Aread(iter_id, H5T_STD_I64LE, &iter);
	
	P0_id = H5Dopen2(group2_id, "P0", H5P_DEFAULT);	
	P0chi_id = H5Aopen(P0_id, "chi", H5P_DEFAULT);
	H5Aread(P0chi_id, H5T_NATIVE_DOUBLE, &chi_m);
	
	P9_id = H5Dopen2(group2_id, "P9", H5P_DEFAULT);	
	P9chi_id = H5Aopen(P9_id, "chi", H5P_DEFAULT);
	H5Aread(P9chi_id, H5T_NATIVE_DOUBLE, &chi_M);
	
	H5Aclose(P0chi_id);	H5Aclose(P9chi_id);	H5Aclose(iter_id);
	H5Dclose(P0_id);	H5Dclose(P9_id);	
	H5Gclose(group_id);	H5Gclose(group2_id);	
	H5Fclose(file_id);	

	double relativeError = fabs(chi_m - chi_M) / fabs(chi_m);
	double allowedRange = 0.05;
//	printf("tot iter:%d  rel: %e all: %e\n", iter, relativeError, allowedRange);

	result = (relativeError <= allowedRange) ? 1 : 0;
	printf("chck: %d\n", result);
	return result; 
}

void fitting(struct Malloc* mem, struct complex_func* cfunc, int seed)
{
#ifdef DEBUG
	printf("seed: %d\n", seed);
#endif

	long long int TotIter_rank = 0;
	// bath optimization
	for(int i = 0; i < Repeat_rank; i++)	mem->Index_rank[i] = i;
	for(int i = 0; i < Repeat_rank; i++)
	{
#ifdef DEBUG
		time_t stime = time(NULL);
#endif
		init_bath( &mem->Bath_rank[i*4*Nb], mem->extrema_index, mem->half_index, mem->peak_num , cfunc->Dw );
		CGM( &mem->Bath_rank[i*4*Nb], &mem->Bath_rank[i*4*Nb + 2*Nb], &mem->Chi_rank[i], cfunc->Diwn, &mem->Iter_rank[i]);
		TotIter_rank += mem->Iter_rank[i];
#ifdef DEBUG
		printf("chi: %e, CG time: %e\n", mem->Chi_rank[i], difftime(time(NULL), stime));
#endif
	}

	quickSort(mem->Chi_rank, mem->Index_rank, 0, Repeat_rank-1);	
	for(int i = 0; i < mem->SortNum_rank; i++) // Sorted by rank and save only 10 data per rank
	{
		mem->SortChi_rank[i] = mem->Chi_rank[i];
		mem->SortIter_rank[i] = mem->Iter_rank[ mem->Index_rank[i] ];
		for(int j = 0; j < 4*Nb; j++)
			mem->SortBath_rank[ i*4*Nb + j ] = mem->Bath_rank[ mem->Index_rank[i] * 4*Nb + j ];
	}
	save_bath( mem->SortBath_rank, mem->SortChi_rank, mem->SortNum_rank, mem->SortIter_rank, TotIter_rank, seed);
}

void bundle(struct G_params* p, struct Malloc* mem , struct complex_func* cfunc)
{
	prepare(p, mem, cfunc);
	int seed = find_seed();
	while( CheckError() == 0 )
	{
		srand(seed);
		fitting(mem, cfunc, seed);
		seed++;
	}
}
