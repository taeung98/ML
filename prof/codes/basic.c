#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "matrix.h"


#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void my_error(char *error_text) {
	fprintf(stderr,"<-------------------------- run-time error -------------------------->\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"<-------------------------- now exiting to system ------------------->\n");
	exit(0);
}

void indexxd(int n, double arr[], int indx[])
{
	int i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	double a;

	istack = mkvectori(NSTACK+2);
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
			if (jstack > NSTACK) my_error("NSTACK too small in indexx.");
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
	for(i=0; i<n; i++)	indx[i+1] -= 1;
	free(istack);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI

int nofile(FILE *fp, char *filename){
	if(fp == NULL){
		printf("There is no file '%s'!\n", filename);
		exit(1);
	}
	return 0;
}

int test_file(char *dir){
	struct stat buf;
	int found;
	int err = stat(dir, &buf);
	if( err == -1 )	found = 0;
	else		found = 1;
	return found;
}

int check_parent(char *dir){
	int i, j, len=0;
	char parent[1024];
	while( dir[len] !='\0' )	len++;
	for(i=len; i>0; i--)
		if( dir[i] == '/' ){
			for(j=0; j<i; j++)	parent[j] = dir[j];
			parent[j] = '\0';
			if( mkdir(parent, S_IRUSR|S_IWUSR|S_IXUSR) != 0 )
				check_parent(dir);
			else{
				mkdir(dir, S_IRUSR|S_IWUSR|S_IXUSR);
				return 1;
			}
		}

	return 0;
}
int mkdirectory(char *dir){
	if( test_file(dir)==0 ){
		if( mkdir(dir, S_IRUSR|S_IWUSR|S_IXUSR) != 0 ){
			if( check_parent(dir) != 0 )
				return 0;
			else
				printf("Fail to generate the folder '%s'\n", dir);
			exit(1);
		}
		else
			printf("A new directory '%s' has been generated\n", dir);
	}
	return 0;
}

int mk_savedir(char *pre){
	int tail=0;
	char savedir[400], order[600];
	sprintf(savedir, "%s", pre);
	do{
		sprintf(order, "test -d %s/", savedir);
		if( system(order) != 0 ){
			sprintf(order, "%s/datafile", savedir);
			mkdirectory(order);
			sprintf(order, "%s/iter", savedir);
			mkdirectory(order);
			sprintf(order, "%s/result", savedir);
			mkdirectory(order);
			sprintf(order, "%s/ongoing", savedir);
			mkdirectory(order);
			break;
		}
		else{
			tail++;
			sprintf(savedir, "%s_%d", pre, tail);
		}

	}while(1);
	sprintf(pre, "%s", savedir);

	return tail;
}

