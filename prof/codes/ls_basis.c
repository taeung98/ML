#include "ls_basis.h"

int find_refer(typebasis **basis, int refer, int block, typebasis *bin){
	int nbasis;
	for(nbasis=refer; nbasis<refer+block; nbasis++)
		if( basis[nbasis][0] == bin[0] && basis[nbasis][1] == bin[1] )	return nbasis;
	
	puts("Error in find_refer! No matching basis!");
	exit(1);
}
int permu_print(typebasis *bin, int site, int spin){
	int nonzero;
	nonzero = count_bit(site, bin[0]) + count_bit(site, bin[1]);
	printf("(%d,%d)[%d]", site, spin, nonzero);
	nonzero += spin * ( (bin[0]>>site) & 0x01 );
	printf("[%d] ", nonzero);
	if( site>Ns || site<0 || spin<0 || spin>1){
		printf("check the variables, site and spin\n");
	}
	return (int) pow(-1,nonzero);
}
int permu(typebasis *bin, int site, int spin){
	int nonzero;
	nonzero = count_bit(site, bin[0]) + count_bit(site, bin[1]);
	//printf("[%d]", nonzero);
	nonzero += spin * ( (bin[0]>>site) & 0x01 );
	//printf("[%d]", nonzero);
	return (int) pow(-1,nonzero);
}

void build_bitsubbasis(int ntotal, typebasis *subset){
	int site, nsub=1, begin=0;
	//subset[nbasis][site] is equal to subset[nbasis] >> site & 1
	while( begin != nsub ){
		site = 0;
		while( site < ntotal && (~subset[begin]>>site&0x01)  ){
			subset[nsub] = subset[begin];
			subset[nsub] = subset[nsub] | 0x01<<site ;
			nsub++; site++;
		}
		begin++;
	}
}

void build_bitbasis(int ntotal, typebasis **basis, sectorinfo *ginfo, int Nonzero_SOC){
	int nsub0,nsub1,i, ptlold;
	int nbasis;
	int *subptl, *ptl;
	int Spinmax = 2*ntotal+1, index;
	typebasis *subbasis;

	subbasis = mkvectorb((0x01<<ntotal));
	initvectorb( subbasis, (0x01<<ntotal));
	build_bitsubbasis(ntotal, subbasis);
	//print_bitsubbasis(ntotal, bitsubbasis);

	ptl	= mkvectori((0x01<<(2*ntotal)));
	subptl	= mkvectori((0x01<<ntotal));

	for(nsub1=0; nsub1<(0x01<<ntotal); nsub1++) subptl[nsub1] = count_bit(ntotal, subbasis[nsub1] );


	int Blocks;

	if( !Nonzero_SOC )	Blocks = (ntotal+1)*(ntotal+1);
	else			Blocks	= 2*ntotal+1;
	for(i=0; i<Blocks; i++)	ginfo[i].block=0;
	nbasis=0;

	if( !Nonzero_SOC ){
		for(ptlold=0; ptlold<Spinmax; ptlold++) for(i=0; i<ptlold+1; i++) {
			for(nsub0=0; nsub0<(0x01<<ntotal); nsub0++) for(nsub1=0; nsub1<(0x01<<ntotal); nsub1++) {
				if ( subptl[nsub0]==i && subptl[nsub1]==ptlold-i ){
					basis[nbasis][0] = subbasis[nsub0];
					basis[nbasis][1] = subbasis[nsub1];
					ptl[nbasis] = subptl[nsub0] + subptl[nsub1];
					index = find_index(ntotal, ptlold, subptl[nsub1]-subptl[nsub0], Blocks);
					nbasis++;
					ginfo[index].block =  ginfo[index].block+1;
				}
			}
		}
	}
	else{
		for(ptlold=0; ptlold<Spinmax; ptlold++) for(i=0; i<ptlold+1; i++) {
			for(nsub0=0; nsub0<(0x01<<ntotal); nsub0++) for(nsub1=0; nsub1<(0x01<<ntotal); nsub1++) {
				if ( subptl[nsub0]==i && subptl[nsub1]==ptlold-i ){
					basis[nbasis][0] = subbasis[nsub0];
					basis[nbasis][1] = subbasis[nsub1];
					ptl[nbasis] = subptl[nsub0] + subptl[nsub1];
					//index = find_index(ntotal, ptlold, subptl[nsub1]-subptl[nsub0]);
					//printf("n0 %d n1 %d, ptl %d, sp%d, index %d\n", nsub0, nsub1, ptlold, subptl[nsub1]-subptl[nsub0], index);
					nbasis++;
					ginfo[ptlold].block =  ginfo[ptlold].block+1;
				}
			}
		}
	}
	for(i=0; i<Blocks; i++) ginfo[i].ground = mkgscmatrixd( Nvec, ginfo[i].block );
	//end of basis construction

	//store useful information to speed up
	for(i=0; i<Blocks; i++){
		ginfo[i].ndeg = 0;
		if( i == 0 )	ginfo[i].refer = 0;
		else		ginfo[i].refer	= ginfo[i-1].refer + ginfo[i-1].block;
		ginfo[i].ptl = count_bit(ntotal, basis[ginfo[i].refer][1]) + count_bit(ntotal, basis[ginfo[i].refer][0]);
		if( !Nonzero_SOC ){
			ginfo[i].spin = count_bit(ntotal, basis[ginfo[i].refer][1]) - count_bit(ntotal, basis[ginfo[i].refer][0]);
			index = find_index(ntotal, ginfo[i].ptl, ginfo[i].spin, Blocks);
			if( i != index ){
				printf("build_bitbasis:: find_index fails in indexing %d != %d\n", i, index);
				exit(1);
			}
		}
		else			ginfo[i].spin = 0;
	}
	free(ptl);
	free(subptl);
	free(subbasis);
}
int print_bitbasis_coeff(int ntotal, typebasis **basis, double *coeff){
	int i;
	printf("==========================================================\n");
	for(i=0; i<(0x01<<(2*ntotal)); i++){
		printf("%4dth [%3u,%3u]  [ptl %2d]  [sp %2d]\t",
				i,
				basis[i][0],
				basis[i][1],
				count_bit(ntotal, basis[i][1])+count_bit(ntotal, basis[i][0]),
				count_bit(ntotal, basis[i][1])-count_bit(ntotal, basis[i][0])
		);
		print_bit(ntotal, basis[i][1] );
		print_bit(ntotal, basis[i][0] );	printf("\t%6.4lf\n", coeff[i]);
	}
	printf("==========================================================\n");
	return 0;
}
int print_bitbasis(int ntotal, typebasis **basis){
	int i;
	for(i=0; i<(0x01<<(2*ntotal)); i++){
		printf("%4dth [%3u,%3u]  [ptl %2d]  [sp %2d]\t",
				i,
				basis[i][0],
				basis[i][1],
				count_bit(ntotal, basis[i][1])+count_bit(ntotal, basis[i][0]),
				count_bit(ntotal, basis[i][1])-count_bit(ntotal, basis[i][0])
		);
		print_bit(ntotal, basis[i][1] );
		print_bit(ntotal, basis[i][0] );	printf("\n");
	}
	return 0;
}

int print_bitsubbasis(int ntotal, typebasis *subset){
	int i;
	for(i=0; i<(0x01<<ntotal); i++){
		printf("%d: %d(%d)\t::\t", i, subset[i], count_bit(ntotal, subset[i]) );
		print_bit(ntotal, subset[i] );	printf("\n");
	}
	return 0;
}

int print_bit(int ntotal, typebasis number){
	int i;
	//static unsigned int size=8*sizeof(typebasis);
	int sz=ntotal;

	char *bin;
	bin = mkvectorc(8*sizeof(typebasis)+1);

	for(i=sz-1; i>=0; i--, number >>= 1 )
		bin[i] = (0x01&number)+'0';
	bin[sz] = '\0';
	printf("%s", bin);
	free(bin);
	return 0;
}
int count_bit(int ntotal, typebasis bit){
	int i, nonzero=0;
	//static int size=8*sizeof(typebasis);
	for(i=0; i<ntotal; i++, bit>>=1)	nonzero += 0x01 & bit;
	return nonzero;
}

int find_gindex(sectorinfo *ginfo, int refer, int block, int Nblocks){
	int i;
	for(i=0; i<Nblocks; i++){
		if( ginfo[i].refer == refer && ginfo[i].block == block )	return i;
	}
	puts("Error in find_gindex!");
	exit(1);

	return 0;
}

int find_index(int ntotal, int ptlgr, int netspingr, int Blocks){
	if(ptlgr < ntotal + 1)	return ptlgr*(ptlgr+1)/2 + ( -netspingr + ptlgr )/2;
	else			return Blocks - (2*ntotal-ptlgr + 1)*(2*ntotal-ptlgr + 2)/2 + (-netspingr + 2*ntotal-ptlgr)/2;
}

int print_ginfo(int ntotal, sectorinfo *ginfo, int Blocks){
	int gindex;
	for(gindex=0; gindex<Blocks; gindex++){
		printf("%03d\t%10d\t%8d\tptl%d\tspin%d\n", gindex, ginfo[gindex].refer, ginfo[gindex].block, ginfo[gindex].ptl, ginfo[gindex].spin);
	}
	return 0;
}

int init_trans_cubic_to_spherical(gsl_complex **YtoC){//	txy = YtoC Ylm
	gsl_complex zero = gsl_complex_rect(0,0);
	gsl_complex temp[NC][NC] = {
		{zero,}
		//{ zero, gsl_complex_rect(1,0), zero},
		//{ gsl_complex_rect( -1/sqrt(2.), 0 ), zero, gsl_complex_rect( 1/sqrt(2.), 0 ) },
		//{ gsl_complex_rect( 0, 1/sqrt(2.) ), zero, gsl_complex_rect( 0, 1/sqrt(2.) ) },
	};
	double unit[2][2] = { {1,0}, {0,1} };
	int i, j;
	for(i=0; i<tNC; i++) for(j=0; j<tNC; j++)
		YtoC[i][j] = gsl_complex_mul_real( temp[i/2][j/2] , unit[i%2][j%2] );

	return 0;
}




