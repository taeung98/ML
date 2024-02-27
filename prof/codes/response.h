#ifndef _Response_
#define _Response_

typedef struct{
	int num_op;
	int dim_op;
	int direction;
	int dimension;
	gsl_complex ***mat;
	int ipair[2];
	char type;
	double *aibi;
} Glist_res;


#endif
