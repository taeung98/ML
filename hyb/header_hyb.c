#include "header/main.h"
#include "header/gsl.h"
#include "header_hyb.h"
#include "header/legendre_poly.h"

//#define DEBUG

double bisection(double desired, double* arr_x, double* arr_y, int size)
{
	int left = 0;
	int right = size - 1;
	while( left<=right ){
		int mid = (left+right)/2;
		double mid_value = (arr_y[left]+arr_y[right])/2.0;
		if( fabs(mid_value-desired) < 1e-6 ){
			return arr_x[mid];
		} else if ( mid_value-desired > 1e-6){
			right = mid - 1;
		} else {
			left = mid + 1;
		}
	}
	double closest_left = arr_x[right];
	double closest_right = arr_x[left];
	return (closest_left + closest_right)/2.0;
}

double NNH(double x, void* parmas) //Next Neighbor Hopping 
{
	struct E_params* p = (struct E_params*)parmas;
	double n,m,h;
	n = p->n; m = p->m; h = p->h; 		

	//normalization of bandwidth
	return -.5 * ( cos(x) + n*cos(2*x) + m*cos(3*x) + h*cos(4*x) ) / ( 1 + n + m + h);
}

void pathname(void* params){
	struct G_params* p = (struct G_params*)params;
	int N = p->e_params.near;
	int n = p->e_params.n * 20;
	int m = p->e_params.m * 20 * 20 / n; 
	int h = p->e_params.h * 20 * 20 * 20 / n / m;
	switch(N){
		case 1: 
			sprintf(Hamilton, "N%1d", N);								break;
		case 2:
			sprintf(Hamilton, "N%1d_n%02d", N, n);						break;
		case 3:
			sprintf(Hamilton, "N%1d_n%02d_m%02d", N, n, m);				break;
		case 4:
			sprintf(Hamilton, "N%1d_n%02d_m%02d_h%02d", N, n, m, h);	break;
		default:
			printf("N-Nearest Neighbor is larger than 4\n");	exit(1);
	}
}

double Gre(double x, void *p){
	struct G_params* params = (struct G_params*)p;
	struct E_params* e_params  = &(params->e_params);
	int n = (params->n);
	double complex W = *(params->W + n);
	double (*E)(double, void*) = NNH;
		
	return creal(  1. / ( W +  mu - E(x,e_params) ) );
}

double Gim(double x, void *p){
	struct G_params* params = (struct G_params*)p;
	struct E_params* e_params  = &(params->e_params);
	int n = (params->n);
	double complex W = *(params->W + n);
	double (*E)(double, void*) = NNH;
	
	return cimag( 1. / ( W +  mu - E(x,e_params)) );
}

void calculate_GD(double *W, double complex *G, double complex *D, void*p, char A){ 
	struct G_params* params = (struct G_params*)p;
//	struct E_params* e_params = &(params->e_params);
	double gre, gim, abserr; abserr=0;
	double epsabs = 1e-10; double epsrel = 1e-10;	
//	e_params->N = 1;
	
	if(A == 'r'){
		for(int n=0; n<Nmax; n++){
			W[n] = Wr(n,Nmax);
			params->W[n] = W[n]+(1.*I)*ETA;
		}
	}else{
		for(int n=0; n<Nmax; n++){
			W[n] = Wn(n,B);
			params->W[n] = (1.*I)*W[n];	
		}
	}
	
	gsl_integration_workspace *S = gsl_integration_workspace_alloc(N_EIGEN);
	gsl_function F1,F2;
	F1.function = &Gre;
	F1.params = params;
	F2.function = &Gim;
	F2.params = params;
	
	for(int n=0; n<Nmax; n++){
		gre = 0; gim =0;
		params->n = n;
					
		gsl_integration_qag(&F1,-PI,PI,epsabs,epsrel,N_EIGEN,6,S,&gre,&abserr);
		gsl_integration_qag(&F2,-PI,PI,epsabs,epsrel,N_EIGEN,6,S,&gim,&abserr);
		
		G[n] = (gre + (1.*I)*gim)/(2*PI);
		D[n] = params->W[n] + mu - 1./G[n];
	}
	
	gsl_integration_workspace_free(S);
}

void fermiNdos(void *params){
	mu = 0;
	struct G_params* p = (struct G_params*)params;
//	struct E_params* e_params = &(p->e_params);

	double epsabs = 1e-10; double epsrel = 1e-10;
	double DOS[N_DOS];	  memset(DOS, 0, sizeof(DOS));
	double energy[N_EIGEN];
	double sum_DOS[N_DOS-1]; memset(sum_DOS, 0, sizeof(sum_DOS)); // # of states(electorns)
	double dos, abserr;
	
	for(int n=0; n<N_EIGEN; n++)
		energy[n] = Wr(n, N_EIGEN);
			
	gsl_integration_workspace *S = gsl_integration_workspace_alloc(N_DOS);
	gsl_function F;
	F.function = &Gim;
	F.params = p; 
	for(int n=0;n<N_DOS;n++){
		dos = 0; abserr=0;
		p->n = n;
		p->W[n] = energy[n] + (1.*I)*ETA;
		gsl_integration_qag(&F, -PI, PI, epsabs, epsrel, N_EIGEN, 6, S, &dos, &abserr);
		DOS[n] = -dos/(2*PI)/PI; // (1/2PI): because k sum & (-1/PI): kk relation	
	}
	gsl_integration_workspace_free(S);
	
	double dE = (energy[N_EIGEN-1] - energy[0])/(N_EIGEN-1);
	sum_DOS[0] = (DOS[1]+DOS[0])/2.*dE;	
	for(int i=0;i<N_DOS-1-1;i++)
		sum_DOS[i+1] = sum_DOS[i] + (DOS[i+2]+DOS[i+1])/2.*dE;
	
	mu = bisection(0.5, energy, sum_DOS, N_DOS-1);
}

void hyb(double complex *D, double *bath, char domain){
	if( domain == 'r' ){
		for(int n = 0; n < Nmax; n++){
			D[n] = 0.;
			for(int i = 0 ; i<Nb; i++){
				D[n] += bath[i]*bath[i]/(w[n] + (1.*I)*ETA - bath[i+Nb]);
			}
		}
	}else{
		for(int n = 0; n < Nmax; n++){
			D[n] = 0.;
			for(int i = 0; i < Nb; i++){
				D[n] += bath[i]*bath[i]/((1.*I)*wn[n] - bath[i+Nb]);
			}
		}
	}
}

void L_coeff(double complex* func, double complex* A, double *ti, double *wi) 
{// calculate Legendre coeff using Gaussian quad
    double (*p)(int, double);
    p = P;

    for(int i=0; i<nquad; i++){
        double complex c=0;
        for(int j=0; j<Lmax; j++){
            c += func[j]*p(i, 2*ti[j]/B-1)*wi[j];
        }
        A[i] = c*(2*i+1)/B;
    }
}

void LT( double *ti, double complex *A, double complex *Y, int num)
{// Legendre transformation
    double (*p)(int, double);
    p = P;
	double complex f;
    for(int i=0;i<Lmax; i++){
        f = 0;
        for(int j=0; j<num; j++){
            f += A[j]*p(j, 2*ti[i]/B-1);
        }
        Y[i] = f;
    }
}

void FT(double complex* b, double* node1, double complex* af, double* node2, double cycle)
{// Fourier transformation using Gaussian quad
	int b_point_num = Nmax;
    int af_point_num = Lmax;
    double complex f =0;
    for(int i=0; i<af_point_num; i++){  //tau index
        f = 0.;
        for(int j=0; j<b_point_num; j++){   // Wn index
            f += (cos(node1[j]*node2[i]) - (1.*I) * sin(node1[j]*node2[i])) * (b[j]-1 / ( (1.*I) * node1[j] ));
        }
        af[i] = f/cycle-0.5;
    }
}

/*
void L_coeff(double complex *func, double* x, double complex* A, double *ti, double *wi, char* fname){
// calculate Legendre coeff using Gaussian quad
	double (*p)(int, double);
	p = P;

	for(int i=0; i<Amax; i++){
		double complex c=0;
		x[i] = (double)i;
		for(int j=0; j<L; j++){
			c += func[j]*p(i, 2*ti[j]/B-1)*wi[j];
		}
		A[i] = c*(2*i+1)/B;
	}
}

void LT( double *ti, double complex *A, double complex *Y, int num){
// Legendre transformation
	double (*p)(int, double);
	p = P;
	double complex f;
	for(int i=0;i<L; i++){
		f = 0;
		for(int j=0; j<num; j++){
			f += A[j]*p(j, 2*ti[i]/B-1);
		}
		Y[i] = f;
	}
}

void FT(double complex* b, double* node1, double complex* af, double* node2, double cycle){
// Fourier transformation using Gaussian quad
	int b_point_num = Nmax;
	int af_point_num = L;
	double complex f =0;
	for(int i=0; i<af_point_num; i++){  //tau index
		f = 0.;
		for(int j=0; j<b_point_num; j++){   // Wn index
			f += (cos(node1[j]*node2[i]) - (1.*I) * sin(node1[j]*node2[i])) * (b[j]-1 / ( (1.*I) * node1[j] ));
		}
		af[i] = f/cycle-0.5;
	}
}

void inverse( double* ti, double complex* A, double complex* y, char* folder){
	int num[5]; for(int i=0;i<3;i++){ num[i]=3+2*i; }; num[3]=30; num[4] = Amax;
	char path[50]="";	snprintf(path, sizeof(path), "%s/%s", Hamilton, folder);
	for(int i=0;i<5;i++){
		char fname[1024];	sprintf(fname, "%d", num[i]);
		LT(ti, A, y, num[i]);
		//save_re(ti, y, L, path, fname);
	}
}
*/

double integral_imag(double *x, double complex *y, int size){
	double result = 0.;
	double dx = (x[size-1] - x[0]) / (size-1);
	for(int i = 1; i < size; i++){
		result += ( cimag(y[i-1]) + cimag(y[i]) ) * dx /2.;
	}
	return result;
}

double abs_sign(double x){
	if(x<0){
		return -1.;
	}
	else{
		return 1.;
	}
}

void Gaussian(double a,double b, int len, double *ni, double *wi){
	//Gaussian Quadrature node & weights
	gsl_integration_glfixed_table *w;
	w = gsl_integration_glfixed_table_alloc(len);
	for(int i=0;i<len;i++){
		gsl_integration_glfixed_point(a,b,i,&ni[i],&wi[i],w);
	} 
	gsl_integration_glfixed_table_free(w);
}

void sort_bath(const gsl_vector *v, double *ep, double *V)
{
	int index[Nb];
	for(int i=0;i<Nb;i++){
		ep[i] = gsl_vector_get(v,Nb+i); // Nb ~ 2*Nb-1
		index[i] = i;
	}
	quickSort(ep, index, 0, Nb-1);
	for(int i = 0; i < Nb; i++)
		V[i] = gsl_vector_get(v, index[i]);	 
}

double obj_f(const gsl_vector *v, void *params){
	double b[2*Nb], result = 0.;
	double complex D, iwn;
	double complex* Diwn = (double complex*)params;
	int n , j;

	for(j=0; j < 2*Nb; j++)
		b[j] = gsl_vector_get(v, j);// 0~Nb-1: V, Nb~2Nb-1: epslion
	
//	sort_bath(v, ep, V);
	for(n=0; n < lim; n++){
		D = 0.;
		iwn = (1.*I)*wn[n];
		for(j=0; j < Nb; j++){
			D += b[j] * b[j] / ( iwn - b[Nb+j] );
		}
		//result += creal( Diwn[n] - D ) * creal( Diwn[n] - D ) + cimag( Diwn[n] - D ) * cimag( Diwn[n] - D ); // hyb diff
		result +=  cabs(Diwn[n] - D) * cabs( Diwn[n] - D);
	}
	return result/lim;
}

void obj_df(const gsl_vector* v, void* params, gsl_vector *df)
{
	double b[2*Nb], dFdV, dFde;
	double complex D, diwn, iwn;
	double complex* Diwn = (double complex*)params;
	double V, ep, omg, ep2, V2, omg2, real, imag;
	int i, j, n;

	for(j=0; j < 2*Nb; j++)
		b[j] = gsl_vector_get(v, j);// 0~Nb-1: V, Nb~2Nb-1: epslion

//	sort_bath(v, ep, V);
	for( i=0; i<Nb; i++){
		dFdV = 0.;	dFde = 0.;
		V = b[i]; ep = b[Nb+i];
		V2 = V*V; ep2 = ep*ep;
		for( n=0; n<lim; n++){
			D=0.;
			iwn = (1.*I)*wn[n]; diwn = Diwn[n]; omg = wn[n]; omg2 = omg*omg;
			for( j = 0; j < Nb; j++)
				D += b[j]*b[j]/(iwn - b[Nb+j]);

			real = creal(diwn - D); imag = cimag(diwn - D);
			//dFdV += ( 4*cabs(V[i])*sign*(ep[i]*creal(Diwn[n]-D) + wn[n]*cimag(Diwn[n]-D)) / (wn[n]*wn[n]+ep[i]*ep[i]) );// wn[n];
			dFdV += ( 4*V*(ep*real + omg*imag) / (omg2+ep2) );// wn[n];
			dFde += ( -2*V2*( ( 2*ep2/(ep2+omg2) - 1.0 )*real + 2*omg*ep*imag/(ep2+omg2) )/(ep2+omg2) );// wn[n];
		}
		gsl_vector_set(df,i,dFdV/lim);
		gsl_vector_set(df,Nb+i,dFde/lim);
	}
}

void obj_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df){
	*f = obj_f(x, params);
	obj_df(x,params,df);
}

double derivative(double* x, double complex* y, int index, int size)
{
	double h = (x[size-1] - x[0])/(size-1);

	if(index>0 && index<size-1){
		return cimag(y[index+1] - y[index-1])/(2.*h);
	}else if (index==0){
		return cimag(y[1] -y[0])/(2.*h);
	}else if (index==size-1){
		return cimag(y[size-1] - y[size-2])/h;
	}else {
		return NAN;
	}
}

int find_extrema(double *x, double complex *y, int* extrema_index)
{
	int peak_num = 0;
	int iter = 0;
	int c[100];	memset(c,0,sizeof(c));	//save epsilon
	double prev_dy = derivative(x,y,0,Nmax);
	
	for(int i=1;i<Nmax;i++){
		double dy = derivative(x,y,i,Nmax);
		if( (prev_dy < 0 &&  dy > 0) || (prev_dy > 0 && dy <0) ){
//			printf("w[%d]=%f: prev_dy = %f, dy = %f\n",i,w[i],prev_dy,dy);
			c[iter] = ( (i-1) + i ) / 2;
			iter++;
			if(prev_dy<0 && dy>0){
				peak_num++;
			}
		}
		prev_dy = dy;
	}
	for(int i=0;i<iter;i++)
		extrema_index[i] = c[i];

	return peak_num;
}	

double randnum(double min,double max){
	return min + ( (double)rand() / RAND_MAX ) * (max - min);
}

void generate_V(double* bath, double err_range, double complex* Dw){
	double v[Nb];
	double sumOfSquares = 0.0;
	double const_sum = -integral_imag(w, Dw, Nmax)/PI/Nb;
	double MIN_RANGE_PERCENT = 1. - err_range;
	double MAX_RANGE_PERCENT = 1. + err_range;
	// Generate random values for v_0 to v_(Nb-2) within the specified range
	for (int i = 0; i < Nb - 1; i++) {
		double minValue = (MIN_RANGE_PERCENT * sqrt(const_sum));
		double maxValue = (MAX_RANGE_PERCENT * sqrt(const_sum));
		v[i] = randnum(minValue,maxValue);
		sumOfSquares += v[i] * v[i];
	}
	// Calculate the square of v_10 to maintain the desired sum of squares
	double lastValueSquared = fabs(const_sum*Nb - sumOfSquares);
	
	v[Nb- 1] = sqrt(lastValueSquared);
	for (int i = 0; i < Nb; i++) 
		bath[i] = v[i];
}

void randSplit(int N, int n, int* splits){
	for(int i=0;i<n;i++)
		splits[i] = 0;
	
	for(int i=0;i<n-1;i++){
		splits[i] = rand() % (N-1);
		N -= splits[i];
	}
	splits[n-1] = N;
}

void generate_ep(double* bath, int* ex_index, int* half_index, int peak_num )
{
	int splits[peak_num], idx[Nb], N; N=Nb;

	randSplit(Nb, peak_num, splits);
	for(int i = 0; i < peak_num; i++){
		for(int j = 0; j < splits[i]; j++)
			bath[N+j] = randnum(w[half_index[2*i]], w[half_index[2*i+1]]);
		
		N += splits[i];
	}
	for(int i = 0; i < Nb; i++)	idx[i] =  i;
	
	quickSort(&bath[Nb], idx, 0, Nb-1); 
}

void init_bath( double* bath, int* ex_index, int* half_index, int peak_num , double complex* Dw){
	generate_V(bath, 0., Dw);
	generate_ep(bath, ex_index, half_index, peak_num);
}

//peak_index: ¿¿¿¿¿¿¿ index, peak_num: ¿¿ ¿¿, peak: ¿¿¿¿¿ ¿¿ ¿¿¿ ¿¿ ¿¿¿¿ index
void FWHM(double* x, double complex* y, int* peak_index, int peak_num, int peak, int* half_index){ //Find Half Maximum
	double half;
	if(peak_num != 1 ){
		if(peak != peak_num-1){
			half = cimag(y[peak_index[peak]]) - ( cimag(y[peak_index[peak]]) - cimag(y[peak_index[peak+1]]) ) * 0.65;
		}else if(peak == peak_num-1){
			half = cimag(y[peak_index[peak]]) - ( cimag(y[peak_index[peak]]) - cimag(y[peak_index[peak-1]]) ) * 0.65;
		}
	}else{
		half = cimag(y[peak_index[peak]]) * 0.5;
	}
	
	int left_index = peak_index[peak]-2;
	while(left_index > 0 && cimag(y[left_index]) < half){
		left_index--;
	}
	half_index[peak] = left_index;
	
	int right_index = peak_index[peak]+2;
	while(right_index <Nmax-1 && cimag(y[right_index]) < half ){
		right_index++;
	}	
	half_index[peak+1] = right_index;
}

void CGM(double *bath_init, double* bath, double* chi, double complex* Diwn, long long int* iter_save){
	long long int iter = 0;
	int status;	
	double step_size, tol, epsabs;
//	time_t stime = time(NULL);
	
	step_size = .5;
	tol = 0.58;
	epsabs = 1e-12;
	
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	gsl_vector *x,*y;
	gsl_multimin_function_fdf obj_func;
	
	obj_func.n = 2*Nb;
	obj_func.f = obj_f;
	obj_func.df = obj_df;
	obj_func.fdf = obj_fdf;
	obj_func.params = Diwn;
	
	// sets the value of the i-th element of a vector v to x.
	x = gsl_vector_alloc(2*Nb);
	for(int i=0;i<2*Nb;i++) gsl_vector_set(x, i, bath_init[i]);		
#ifdef DEBUG
	for(int i=0;i<2*Nb;i++) gsl_vector_set(x, i, randnum(w[0], w[Nmax]));		
#endif
	T = gsl_multimin_fdfminimizer_conjugate_pr;
	s =  gsl_multimin_fdfminimizer_alloc (T, 2*Nb);
	gsl_multimin_fdfminimizer_set(s, &obj_func, x, step_size, tol);	
	
//	printf("init f = %e/n", s->f);
//	conjugate_state_t* vstate = (conjugate_state_t*)s->state; 
	do{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);
	
		if(status == GSL_ENOPROG) /* iteration is not making progress towards solution (return 27)*/
			break;
		
		status = gsl_multimin_test_gradient(s->gradient, epsabs);
	}while(status == GSL_CONTINUE); 

//	gsl_multimin_fdfminimizer_restart(s);
//	printf("status: %d, iter = %d, tol = %e, step_size = %e, epsabs = %e, chi = %e\n", status, iter, vstate->tol, vstate->step, epsabs, s->f);
//	printf("chi: %.2e\n",  s->f);

	y = gsl_vector_alloc(2*Nb);
	for(int i = 0; i < 2*Nb; i++) gsl_vector_set(y, i, gsl_vector_get(s->x, i));
	sort_bath(y, &bath[Nb], &bath[0]);

	*chi = obj_f(y, Diwn);
	*iter_save = iter;

	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);
}

void print_bath(double *bath, double *chi, int num)
{
	printf("chi = %e, test_f = %e\n", chi[num], test_f(&bath[num*4*Nb+2*Nb],NULL));
/*	for(int i = 0; i < 4*Nb; i++){
		if(i == 0) printf("init V:\n");
		if(i == Nb) printf("\ninit e:\n ");
		if(i == 2*Nb) printf("\nfinl V:\n ");
		if(i == 3*Nb) printf("\nfinl e:\n ");
		printf("bath[%d] = %e;\n ", i%(2*Nb), bath[num*4*Nb + i]);
	}printf("\n\n");*/
}

void printT(char* fname) 
{
	time_t current_time;
	time(&current_time);
	char format_time[100];
	struct tm* time_info;
	time_info = localtime(&current_time);
	strftime(format_time, sizeof(format_time), "%Y-%m-%d %H:%M:%S", time_info);
	printf("%s time: %s\n", fname, format_time);
}

void swap_int(int* a, int* b) {
	int temp = *a;
	*a = *b;
	*b = temp;
}

void swap_double(double* a, double* b){
	double temp = *a;
	*a = *b;
	*b = temp;
}

// Partition the array into two sub-arrays and return the pivot index
int partition(double* arr, int* indices, int low, int high){
	double pivot = arr[high];
	int i = (low - 1);

	for (int j = low; j <= high - 1; j++) {
		if (arr[j] < pivot) {
			i++;
			swap_double(&arr[i], &arr[j]);
			swap_int(&indices[i], &indices[j]);
		}
	}

	swap_double(&arr[i + 1], &arr[high]);
	swap_int(&indices[i + 1], &indices[high]);
	return (i + 1);
}

// Main quick sort function
void quickSort(double* arr, int* indices, int low, int high) {
	if (low < high) {
		int pi = partition(arr, indices, low, high);
		quickSort(arr, indices,  low,  pi - 1);
		quickSort(arr, indices, pi + 1, high);
	}
}

void sort_test(double *v, double *ep, double *V){
	int index[Nb];
	for(int i=0; i < Nb; i++){
		ep[i] = v[Nb+i]; // Nb ~ 2*Nb-1
		index[i] = i;
	}
	quickSort(ep, index, 0, Nb-1);
	for(int i = 0; i < Nb; i++)
		V[i] = v[index[i]];  
} 

double test_f(double *v, void *params){
	double V[Nb], ep[Nb], result = 0.; 
	double complex D;
	double complex* Diwn = (double complex*)params;

	sort_test(v, ep, V); 

	for(int n=0;n<lim;n++){
		D = 0.; 
		for(int j=0;j<Nb;j++){
			D += cabs(V[j])*cabs(V[j])/((1.*I)*wn[n]-ep[j]);
		}   
		result += (creal( Diwn[n] - D ) * creal( Diwn[n] - D ) + cimag( Diwn[n] - D ) * cimag( Diwn[n] - D )); // hyb diff
	}   
	return result/lim;
}

void save_freq(hid_t file_id, double* W, double* WN){
	hid_t dataset_id, dataspace_id, attribute_id, attr_dataspace_id;
	hsize_t dims[2]; 
	hsize_t attr_dims[1];
	
	dims[0] = 2; dims[1] = Nmax;
	double freq[2][Nmax];	// 0: real_domain, 1: iamg_domain 
	for(int i = 0; i<Nmax; i++){
		freq[0][i] = W[i];
		freq[1][i] = WN[i];
	}
/*	
	file_id = H5Fopen(FILE_h5, H5F_ACC_RDWR, H5P_DEFAULT); 
	if(file_id < 0)
		file_id	= H5Fcreate(FILE_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
*/	
	char name[16]; sprintf(name, "/omega");	
	dataspace_id = H5Screate_simple(2, dims, NULL);
	dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
	if(dataset_id < 0){
		dataset_id = H5Dcreate2(file_id, name, H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}else{
		H5Sclose(dataspace_id);
//		H5Fclose(file_id);	
		return;
	}
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, freq);	
		
	int attr_data = Nmax; 
	attr_dims[0] = 1;	
	attr_dataspace_id = H5Screate_simple(1, attr_dims, NULL);	
	attribute_id = H5Acreate2(dataset_id, "Nmax", H5T_STD_I32BE, attr_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data);
	
	H5Aclose(attribute_id);
	H5Sclose(attr_dataspace_id);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);	
//	H5Fclose(file_id);

	return;
}

void save_data(hid_t file_id, hid_t group1_id, double complex* data, char* dataset_name, int len)
{
	hid_t /*file_id, group1_id,*/ dataset_id, dataspace_id;
	hsize_t func_dims[2]; func_dims[0] = 2; func_dims[1] = len; 
	double func[2][len];//0: real_part, 1: imag_part

	for(int i = 0; i < len; i++)
	{
		func[0][i] = creal(data[i]);
		func[1][i] = cimag(data[i]);
	}
/*
	file_id = H5Fopen(FILE_h5, H5F_ACC_RDWR, H5P_DEFAULT); 
	if(file_id < 0)
		file_id = H5Fcreate(FILE_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	char Group1[1024]; 	sprintf(Group1, "/%s", Hamilton); 
	group1_id = H5Gopen2(file_id, Group1, H5P_DEFAULT);
	if( group1_id < 0)
		group1_id = H5Gcreate2(file_id, Group1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
*/
	dataspace_id = H5Screate_simple(2, func_dims, NULL);
	char Dname[16]; sprintf(Dname, "%s", dataset_name);

	dataset_id = H5Dopen2(group1_id, Dname, H5P_DEFAULT);
	if(dataset_id < 0){
		dataset_id = H5Dcreate2(group1_id, Dname, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}else{	
		H5Sclose(dataspace_id);	
		H5Dclose(dataset_id);
//		H5Gclose(group1_id);
//		H5Fclose(file_id);
		return;
	}
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, func);	
		
	H5Sclose(dataspace_id);	
	H5Dclose(dataset_id);
//	H5Gclose(group1_id);
//	H5Fclose(file_id);
	
	return;
}

void load_chi(hid_t path_id, int chi_num, double *chi_val)
{
	hid_t dataset_id, attr_id;

	for (int i = 0; i < chi_num; i++) {
		char dataset_name[20];
		sprintf(dataset_name, "P%d", i);
		// Open dataset
		dataset_id = H5Dopen2(path_id, dataset_name, H5P_DEFAULT);
#ifdef DEBUG
		double dataset_data[2][2*Nb];
		hid_t dataspace_id;
		hsize_t dims[2];
		if (dataset_id < 0) {
			H5Dclose(dataset_id);
			fprintf(stderr, "Unable to open dataset: %s\n", dataset_name);
			return;
		}
		// Open dataspace
		dataspace_id = H5Dget_space(dataset_id);
		H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
		// Read the dataset into the array
		if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_data) < 0) {
			fprintf(stderr, "Error reading dataset\n");
			H5Sclose(dataspace_id);
			H5Dclose(dataset_id);
			return;
		}
		printf("Entire dataset:\n");
		for (int j = 0; j < dims[1]; j++) 
			printf("%e\t%e\n ", dataset_data[0][j], dataset_data[1][j]);
		printf("\n");
#endif
		// Read 'chi' attribute
		attr_id = H5Aopen(dataset_id, "chi", H5P_DEFAULT);
#ifdef DEBUG
		if (attr_id < 0) {
			fprintf(stderr, "Unable to open attribute 'chi' for dataset: %s\n", dataset_name);
			H5Aclose(attr_id);
			H5Dclose(dataset_id);
			return;
		}
		H5Sclose(dataspace_id);
#endif		

		H5Aread(attr_id, H5T_NATIVE_DOUBLE, &(chi_val[i]));
//		printf("chi = %e, test_f = %e\n", chi_val[i], test_f(dataset_data[1], NULL));
		// Close dataset and dataspace
		H5Aclose(attr_id);
		H5Dclose(dataset_id);
	}
}

void process_datasets(hid_t path_id) // path_id -> bath{Nb}
{
	hid_t dataset_id;
//	herr_t status;
	hsize_t num_objs;	H5Gget_num_objs(path_id, &num_objs);
	int dataset_count = (int)num_objs;	
	double chi_values[dataset_count];

	if(dataset_count == 10) 
		return;
#ifdef DEBUG
	printf("\n accumulate data\n");
	printf("\n dataset_count = %d\n", dataset_count);
#endif
	load_chi(path_id, dataset_count, chi_values);

	int chi_index[dataset_count];
	for(int i = 0; i < dataset_count; i++)
	{
		 chi_index[i] = i;
	//	 printf("load: chi[%d] = %e\n", i ,chi_values[i]);
	}

	// Sort chi_values in descending order
	quickSort(chi_values, chi_index, 0, dataset_count-1);	

	// Open datasets in descending order of chi_values
	for (int i = 0; i < dataset_count; i++)
	{
		char dataset_name[20];
		sprintf(dataset_name, "P%d", chi_index[i]);
		dataset_id = H5Dopen2(path_id, dataset_name, H5P_DEFAULT);
		// Remove datasets with higher chi_values
		if (i < 10) {
			// Rename datasets with lower chi_values
			char new_dataset_name[20];
			sprintf(new_dataset_name, "temp%d", i);
/*status=*/	H5Lmove(path_id, dataset_name, path_id, new_dataset_name, H5P_DEFAULT, H5P_DEFAULT);
/*
			if (status < 0) 
				fprintf(stderr, "Error renaming dataset: %s to %s\n", dataset_name, new_dataset_name);
*/
		}else {
/*status*/	H5Gunlink(path_id, dataset_name);
			H5Gget_num_objs(path_id, &num_objs);
/*
			if (status < 0) 
				fprintf(stderr, "Error removing dataset: %s\n", dataset_name);
*/
		}
		H5Dclose(dataset_id);
	}
	for (int i = 0; i < 10; i++)
	{
		char dataset_name[20];	sprintf(dataset_name, "temp%d", i);
		dataset_id = H5Dopen2(path_id, dataset_name, H5P_DEFAULT);
		char new_dataset_name[20]; sprintf(new_dataset_name, "P%d", i);
/*status*/ H5Lmove(path_id, dataset_name, path_id, new_dataset_name, H5P_DEFAULT, H5P_DEFAULT);
		H5Dclose(dataset_id);
	}	
#ifdef DEBUG	
	printf("\n after sort: \n");	
	load_chi(path_id, dataset_count, chi_values);
#endif
}

void save_bath(double* sdata, double* chi, int snum, long long int* Iter, long long int TotIter, int seed)
{
	hid_t file_id, dataset_id, dataspace_id, group1_id, group2_id, chi_id, chi_dataspace_id, Iter_id, Iterspace_id, TotIter_id, TotIterspace_id, seed_id, seedspace_id;
	hsize_t dims[2], D1_dim[1]; D1_dim[0] = 1;	
	long long int Iter_data;
	double chi_data, bath[2][2*Nb];

	file_id = H5Fopen(FILE_h5, H5F_ACC_RDWR, H5P_DEFAULT);
	if(file_id < 0)
		file_id = H5Fcreate(FILE_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	dims[0] = 2; dims[1] = 2*Nb;
	
	char Group1[1024], Group2[1024], Data_num[128]; 
	sprintf(Group1, "/%s", Hamilton); sprintf(Group2, "%s%d", "bath", Nb);

	group1_id = H5Gopen2(file_id, Group1, H5P_DEFAULT);
	if( group1_id < 0)
		group1_id = H5Gcreate2(file_id, Group1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	seedspace_id = H5Screate_simple(1, D1_dim, NULL);	
	seed_id = H5Acreate2(group1_id, "seed", H5T_NATIVE_INT, seedspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(seed_id, H5T_NATIVE_INT, &seed);
	H5Sclose(seedspace_id);		H5Aclose(seed_id);
		
	group2_id = H5Gopen2(group1_id, Group2, H5P_DEFAULT);	
	if( group2_id < 0)
		group2_id = H5Gcreate2(group1_id, Group2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 

	TotIter_id = H5Aopen(group2_id, "Tot_Iter", H5P_DEFAULT);
	if ( TotIter_id < 0 ) {
		TotIterspace_id = H5Screate_simple(1, D1_dim, NULL);
		TotIter_id = H5Acreate2(group2_id, "Tot_Iter", H5T_STD_I64LE, TotIterspace_id, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(TotIter_id, H5T_STD_I64LE, &TotIter);
		H5Sclose(TotIterspace_id);
	}else {
		long long int TotIter_up;
		H5Aread(TotIter_id, H5T_STD_I64LE, &TotIter_up);
		TotIter_up += TotIter;
		H5Awrite(TotIter_id, H5T_STD_I64LE, &TotIter_up);
	}
	H5Aclose(TotIter_id);

	hsize_t num_objs; H5Gget_num_objs( group2_id, &num_objs);

	int snum1 = (int)num_objs; 
//	printf("snum1 = %d\n", snum1);
	int snum2 = snum1 + snum;
	
	for(int i = snum1; i < snum2; i++)
	{
		int k = i - snum1;
//		printf("k = %d\n", k);
		for(int j = 0; j < 2*Nb; j++)
		{
			bath[0][j] = sdata[k*4*Nb + j];
			bath[1][j] = sdata[k*4*Nb + 2*Nb + j];
//			printf("%e\t%e\n", bath[0][j], bath[1][j]);
		}

//		printf("chi = %e, test_f = %e\n\n", chi[k], test_f(&sdata[k*4*Nb + 2*Nb], NULL) );

		sprintf(Data_num,"%c%d",'P',i);
		dataspace_id = H5Screate_simple(2, dims, NULL);
		dataset_id = H5Dcreate2(group2_id, Data_num, H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, bath);	
		
		chi_data = chi[k];	D1_dim[0] = 1;	
		Iter_data = Iter[k];
		//printf("%d\n", Iter_data);
		chi_dataspace_id = H5Screate_simple(1, D1_dim, NULL);	
		Iterspace_id = H5Screate_simple(1, D1_dim, NULL);	
	
		chi_id = H5Acreate2(dataset_id, "chi", H5T_IEEE_F64BE, chi_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		Iter_id = H5Acreate2(dataset_id, "Iter", H5T_STD_I64LE, Iterspace_id, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(chi_id, H5T_NATIVE_DOUBLE, &chi_data);
		H5Awrite(Iter_id, H5T_STD_I64LE, &Iter_data);

		H5Sclose(Iterspace_id);		H5Aclose(Iter_id);		
		H5Sclose(chi_dataspace_id);	H5Aclose(chi_id);	
		H5Sclose(dataspace_id);		H5Dclose(dataset_id);
	}
	process_datasets(group2_id);

	H5Gclose(group2_id);
	H5Gclose(group1_id);
	H5Fclose(file_id);
}

int compare_chi(const void *a, const void *b) {
	double chi_a = *((double *)a);
	double chi_b = *((double *)b);

	if (chi_a < chi_b) return -1;
	if (chi_a > chi_b) return 1;
	return 0;
}
