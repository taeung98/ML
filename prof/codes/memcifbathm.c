#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lattice.h"
#include "matrix.h"
#include "hamil.h"

#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
#define NR_END 1
#define FREE_ARG char*

int sort_map_column(mapele *mapinfo, int *istart, sectorinfo ginfo[], int gindex);
int init_all(int init);

int count_bit(int ntotal, typebasis bit);
extern int Powns, Powns2, nctot, Spinmax, Blocks, nb;

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

void free_ivector(int *v, long nl)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void indexx(unsigned long n, typebasis arr[], unsigned long indx[])
{
	unsigned long i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	typebasis a;

	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=1;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[l]]) {
				SWAP(indx[l+1],indx[l])
			}
			i=l+1;
			j=ir;
			indxt=indx[l];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI


int sort_map_column(mapele *mapinfo, int *istart, sectorinfo ginfo[], int gindex){
	mapele *tempmap;
	int nn;
	unsigned long *index;
	int i, row, block;
	typebasis *column;

	block = ginfo[gindex].block;

	for(row=0; row<block; row++){
		nn = istart[row+1] - istart[row];
		if( nn < 2 )	continue;
		column = mkvectorb( nn+1 );
		index = mkvectorul( nn+1 );

		for(i=0; i<nn; i++)	column[i+1] = mapinfo[istart[row]+i].column;

		indexx(nn, column, index);

		tempmap = mkvectormapele(nn);
		for(i=0; i<nn; i++)	tempmap[i] = mapinfo[istart[row]+i];
		for(i=0; i<nn; i++)	mapinfo[istart[row]+i] = tempmap[index[i+1]-1];
		free(tempmap);
		free(index);
		free(column);
	}
	return 0;
}



