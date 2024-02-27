#include "hamil.h"
#include "basic.h"
#include "memcifbathm.c"

int construct_hamil_direct(int *istart, double *Hdia, gsl_complex *Hamil, int *column, gsl_complex **matrix, int block){
	int mu, nu;
	gsl_complex *vectorin = mkgscvectord(block);
	gsl_complex *vectorout = mkgscvectord(block);
	for(mu=0; mu<block; mu++){
		for(nu=0; nu<block; nu++)
			matrix[mu][nu] = zero;
		//matrix[mu][mu] = gsl_complex_rect(Hdia[mu], 0);
	}

	for(mu=0; mu<block; mu++){
		initgscvectord(vectorin, block);
		vectorin[mu] = gsl_complex_rect(1, 0);
		H_vector_map_complex(istart, Hdia, Hamil, column, vectorin, vectorout, block);
		for(nu=0; nu<block; nu++)
			matrix[nu][mu] = vectorout[nu];

	}
	//print_gscmatrixd("hamil", matrix, block, block);
	free(vectorin);
	free(vectorout);
	return 0;
}

double LiuDavidson(int *istart, double *Hdia, gsl_complex *Hamil, int *column, int block, gsl_complex **vector_in0, gsl_complex **vector_out, double *overlap, int *lancount, double *values, double *res, double Tolerance){
	int i, nk, numk, mu, nu, Nend, Nstart, k, subdim = (Nvec<block? Nvec:block), ngs;
	double inn_r, inn_i, sum;
	char JOBZ = 'V', UPLO='U';
	integer N, LDA, LWORK = (JDmax*SDmax)*(JDmax*SDmax)+(JDmax*SDmax), INFO;
	int nwork_temp = LWORK;

	doublecomplex *WORK = (doublecomplex*) malloc( nwork_temp * sizeof(doublecomplex) );
	doublereal *RWORK = (doublereal*) malloc( nwork_temp * sizeof(doublereal) );
	doublereal *W = (doublereal*) malloc( nwork_temp * sizeof(doublereal) );


	if( block < 100 ){
		gsl_complex **ha = mkgscmatrixd( block, block );
		gsl_complex **ev = mkgscmatrixd( block, block );
		double *energy = mkvectord( block );
		construct_hamil_direct(istart, Hdia, Hamil, column, ha, block);

		eigen_lapack( ha, energy, ev, block, 1);
		for(nk=0; nk<subdim; nk++) for(i=0; i<block; i++) vector_out[nk][i] = ev[nk][i];
		for(nk=0; nk<(SDmax<block?SDmax:block); nk++){
			res[nk] = 0;
			values[nk] = energy[nk];
		}
		freegscmatrixd(ha, block);
		freegscmatrixd(ev, block);
		free(energy);
		return values[0];
	}
	else{
		gsl_complex **V, **HV, **r, **u, **t;
		doublecomplex subspace[(JDmax*SDmax)*(JDmax*SDmax)];
		for(mu=0; mu<(JDmax*SDmax)*(JDmax*SDmax); mu++){
			subspace[mu].r = 0;
			subspace[mu].i = 0;
		}

		Nend =  block < JDmax*SDmax ? block/SDmax : JDmax;
		//Nstart = (*lancount) < 1 ? 0 : 1;
		Nstart = (*lancount) < 1 ? 0 : 1;

		u = mkgscmatrixd(Nend*SDmax, block);
		t = mkgscmatrixd(SDmax, block);
		r = mkgscmatrixd(SDmax, block);
		HV = mkgscmatrixd(Nend*SDmax, block);
		V = mkgscmatrixd((Nend+1)*SDmax, block);

		nk = Nstart;
		numk = (Nstart+1)*SDmax;
		if( nk )	Gram_Schmidt_gsl(vector_in0, numk, block, 1);
		for(k=0; k<numk; k++){
			for(i=0; i<block; i++)
				V[k][i] = vector_in0[k][i];
			H_vector_map_complex(istart, Hdia, Hamil, column, V[k], HV[k], block);
		}
		N=k;

		do{
			//printf("nk=%d, numk=%d\n", nk, numk);	fflush(stdout);
			//print_gscmatrixd("V", V, numk, block);
			//test_GS(V, numk, block);

			//compute HV for new vectors in the subspace: not required at 1st iteration, where numk==N
			for(nu=0; nu<numk-N; nu++)
				H_vector_map_complex(istart, Hdia, Hamil, column, V[N+nu], HV[N+nu], block);
			N = numk;		//size of the subspace is now updated
			LDA = N;
			for(mu=0; mu<numk; mu++) for(nu=mu; nu<numk; nu++){
				inn_r = 0; inn_i=0;
#pragma omp parallel for default(none) private(i) shared(HV, V, block, mu, nu) reduction(+:inn_r, inn_i)
				for(i=0; i<block; i++){
					inn_r += gsl_complex_in_r( V[mu][i], HV[nu][i] );
					inn_i += gsl_complex_in_i( V[mu][i], HV[nu][i] );
				}
				subspace[mu*numk+nu].r =  0;
				subspace[mu*numk+nu].i =  0;

				subspace[nu*numk+mu].r =  inn_r;
				subspace[nu*numk+mu].i =  inn_i;
			}
		   	//print_gscmatrixd("subspace", matrix_temp, numk, numk);

			zheev_(&JOBZ, &UPLO, &N, subspace, &LDA, W, WORK, &LWORK, RWORK, &INFO);
			if( (int)INFO != 0 ){
				printf("LiuDavidson: Failed in zheev_, INFO = %d, Nk = %d, numk=%d, lancount=%d\n", (int)INFO, nk, numk, *lancount);	fflush(stdout);
				print_gscmatrixd("vectorin", vector_in0, nk*SDmax, block);
				print_vectord("eval", W, numk);
				print_vectord("res", res, SDmax);
				gsl_complex **matrix = mkgscmatrixd(SDmax*JDmax, SDmax*JDmax);

				for(mu=0; mu<numk; mu++) for(nu=mu; nu<numk; nu++){
					inn_r = 0; inn_i=0;
#pragma omp parallel for default(none) private(i) shared(HV, V, block, mu, nu) reduction(+:inn_r, inn_i)
					for(i=0; i<block; i++){
						inn_r += gsl_complex_in_r( HV[mu][i], V[nu][i] );
						inn_i += gsl_complex_in_i( HV[mu][i], V[nu][i] );
					}
					subspace[mu*numk+nu].r =  0;
					subspace[mu*numk+nu].i =  0;

					subspace[nu*numk+mu].r =  inn_r;
					subspace[nu*numk+mu].i =  inn_i;

					//matrix[nu][mu] = gsl_complex_rect(inn_r, inn_i);
					//matrix[mu][nu] = gsl_complex_rect(inn_r, -inn_i);
				}
				print_gscmatrixd("V", V, numk, block);
				print_gscmatrixd("HV", HV, numk, block);
				for(mu=0; mu<numk; mu++) for(nu=0; nu<numk; nu++){
					matrix[mu][nu] = gsl_complex_rect( subspace[nu*numk+mu].r, subspace[nu*numk+mu].i );
				}
				print_gscmatrixd("error zheev_", matrix, numk, numk);

				freegscmatrixd(matrix, SDmax*JDmax);
				exit(1);
			}
			for(i=0; i<SDmax; i++)
				values[i] = W[i];

			for(nu=0; nu<numk; nu++){
#pragma omp parallel for default(none) private(i, mu) shared(u, V, block, nu, zero, numk, subspace)
				for(i=0; i<block; i++){
					u[nu][i] = zero;
					for(mu=0; mu<numk; mu++) {
						u[nu][i] = gsl_complex_add(
								u[nu][i],
								gsl_complex_mul(
									gsl_complex_rect(subspace[nu*numk+mu].r, subspace[nu*numk+mu].i),
									//gsl_matrix_complex_get(evec, mu, nu),
									V[mu][i]
									)
								);
					}
				}
			}
			for(nu=0; nu<SDmax; nu++){
#pragma omp parallel for default(none) private(i, mu) shared(nu,numk, u, r, HV, block, W, zero, subspace)
				for(i=0; i<block; i++){
					r[nu][i] = zero;
					for(mu=0; mu<numk; mu++) {
						r[nu][i] = gsl_complex_add(
								r[nu][i],
								gsl_complex_mul(
									gsl_complex_rect(subspace[mu+numk*nu].r, subspace[mu+numk*nu].i),
									//gsl_matrix_complex_get(evec, mu, 0),
									HV[mu][i]
									)
								);
					}
					r[nu][i] = gsl_complex_sub(
							r[nu][i],
							gsl_complex_mul_real(u[nu][i], W[nu])
							);
				}
			}

			for(nu=0; nu<SDmax; nu++){
				res[nu] = 0;
#pragma omp parallel for default(none) private(i) shared(nu, t, r, block, W, Hdia)
				for(i=0; i<block; i++)
					t[nu][i] = gsl_complex_div_real( r[nu][i], W[nu] - Hdia[i] );
				inn_r=0;
#pragma omp parallel for default(none) private(i) shared(nu, t, r, block, W, Hdia) reduction(+:inn_r)
				for(i=0; i<block; i++)
					inn_r += gsl_complex_abs2(r[nu][i]);
				res[nu] = inn_r;
			}

				//print_gscmatrixd("t", t, SDmax, block);

			for(nu=0; nu<SDmax; nu++){
				normalize_gsc(t[nu], block);
				for(i=0; i<block; i++)	V[numk][i] = t[nu][i];
				sum = normalize_gsc(V[numk], block);
				if( sum == NAN ){
					printf("NAN: numk=%d\n", numk);	fflush(stdout);
					continue;
				}
				numk++;
				Gram_Schmidt_gsl(V, numk, block, 0);
				Gram_Schmidt_gsl(V, numk, block, 0);
			}
			Gram_Schmidt_gsl(V, numk, block, 1);

			ngs = 0;
			while( test_GS(V, numk, block) != 0 ){
				Gram_Schmidt_gsl(V, numk, block, 1);
				ngs++;
				if( ngs > 5 ){
					printf("LiuDavidson: overlap is too large\n");
					exit(1);
				}
			}
			if( (int)N == numk ){
				for(nu=0; nu<SDmax; nu++){
					for(i=0; i<block; i++)	V[numk][i] = gsl_complex_rect( ran2(&ranseed)-0.5, ran2(&ranseed)-0.5 );
					normalize_gsc(V[numk], block);
					numk++;
					Gram_Schmidt_gsl(V, numk, block, 0);
					Gram_Schmidt_gsl(V, numk, block, 0);
				}
			}

			//print_vectord("res", res, SDmax);
			//printf("in the end: nk=%d, numk=%d\n", nk, numk);

			nk++;
		}while( nk<Nend && numk<block && res[0]>Tolerance && res[1]>Tolerance );
		//freegscmatrixd(matrix_temp, SDmax*JDmax);

		for(nk=0; nk<SDmax; nk++){
			for(i=0; i<block; i++){
				vector_out[nk][i] = u[nk][i];
			}
			normalize_gsc(vector_out[nk], block);
		}
		if( nk == 1 ){
			for(nk=0; nk<SDmax; nk++){
				for(i=0; i<block; i++){
					vector_out[nk+SDmax][i] = gsl_complex_rect( ran2(&ranseed)-0.5, ran2(&ranseed)-0.5 );
				}
				normalize_gsc(vector_out[nk+SDmax], block);
			}
			Gram_Schmidt_gsl(vector_out, SDmax, block, 1);
		}
		else{
			for(nk=0; nk<SDmax; nk++){
				for(i=0; i<block; i++){
					vector_out[nk+SDmax][i] = V[nk+SDmax][i];
				}
				normalize_gsc(vector_out[nk], block);
			}
		}


		freegscmatrixd(HV, Nend*SDmax);
		freegscmatrixd(V, (Nend+1)*SDmax);
		freegscmatrixd(u, Nend*SDmax);
		freegscmatrixd(t, SDmax);
		freegscmatrixd(r, SDmax);

	}

	//print_gscmatrixd("vector_out", vector_out, Nvec, block);

	inn_r = 0; inn_i = 0;
#pragma omp parallel for default(none) private(i) shared(block, vector_in0, vector_out) reduction(+:inn_r, inn_i)
	for(i=0; i<block; i++){
		inn_r += gsl_complex_in_r( vector_out[0][i], vector_in0[0][i] );
		inn_i += gsl_complex_in_i( vector_out[0][i], vector_in0[0][i] );
	}
	(*overlap) = sqrt(inn_r*inn_r + inn_i*inn_i);
	//freegscmatrixd(matrix, JDmax);

	int ndeg=0;
	for(nk=1; nk<SDmax; nk++)
		if( values[nk] - values[0] < 10e-14)	ndeg++;

	if(ndeg){
		for(nk=1; nk<ndeg+1; nk++){
			inn_r=0; inn_i=0;
#pragma omp parallel for default(none) private(i) shared(block, vector_in0, vector_out, nk) reduction(+:inn_r, inn_i)
			for(i=0; i<block; i++){
				inn_r += gsl_complex_in_r( vector_out[nk][i], vector_in0[0][i] );
				inn_i += gsl_complex_in_i( vector_out[nk][i], vector_in0[0][i] );
			}
			(*overlap) = sqrt((*overlap)*(*overlap) + inn_r*inn_r + inn_i*inn_i);
		}
	}
	(*lancount) = (*lancount) + 1;

	//printf("ndeg=%d, W[0] = %.18lf\toverlap = %.18lf\n", ndeg, W[0], *overlap);
	//for(nk=0; nk<SDmax; nk++) printf("%lf\t", values[nk]); printf("\n");
	free(WORK); free(RWORK); free(W);
	return values[0];
}

double Davidson(int *istart, double *Hdia, gsl_complex *Hamil, int *column, int block, gsl_complex **vector_in0, gsl_complex **vector_out, double *overlap, int *lancount, double *values, double *res, double Tolerance){
	int i, nk, numk, mu, nu, Nend, Nstart, k;
	double theta, inn_r, inn_i;
	char JOBZ = 'V', UPLO='U';
	integer N, LDA, LWORK = (JDmax)*(JDmax)+(JDmax), INFO;
	doublecomplex WORK[(JDmax)*(JDmax)+(JDmax)];
	doublereal RWORK[(JDmax)*(JDmax)], W[(JDmax)*(JDmax)];


	//gsl_complex **matrix;
	//matrix = mkgscmatrixd((JDmax), (JDmax));

	if( block < 100 ){
		gsl_complex **ha = mkgscmatrixd( block, block );
		gsl_complex **ev = mkgscmatrixd( block, block );
		double *energy = mkvectord( block );
		int subdim = (Nvec<block? Nvec:block);
		construct_hamil_direct(istart, Hdia, Hamil, column, ha, block);

		eigen_lapack( ha, energy, ev, block, 1);
		for(nk=0; nk<subdim; nk++) for(i=0; i<block; i++) vector_out[nk][i] = ev[nk][i];
		for(nk=0; nk<(Nvec<block?Nvec:block); nk++){
			res[nk] = 0;
			values[nk] = energy[nk];
		}
		freegscmatrixd(ha, block);
		freegscmatrixd(ev, block);
		free(energy);
		return values[0];
	}
	else{
		gsl_complex **V, **HV, *r, **u, *t;
		doublecomplex subspace[(JDmax)*(JDmax)];
		for(mu=0; mu<(JDmax)*(JDmax); mu++){
			subspace[mu].r = 0;
			subspace[mu].i = 0;
		}

		Nend =  block < (JDmax) ? block : (JDmax);
		Nstart = (*lancount) < 1 ? 0 : Nvec-1;

		u = mkgscmatrixd(Nend, block);
		t = mkgscvectord(block);
		r = mkgscvectord(block);
		HV = mkgscmatrixd(Nend, block);
		V = mkgscmatrixd(Nend+1, block);

		nk = Nstart;
		numk = Nstart+1;
		Gram_Schmidt_gsl(vector_in0, numk, block, 1);
		for(k=0; k<numk; k++){
			for(i=0; i<block; i++)
				V[k][i] = vector_in0[k][i];
		}
		for(k=0; k<numk; k++)
			H_vector_map_complex(istart, Hdia, Hamil, column, V[k], HV[k], block);
		//print_gscvectord("input", vector_in0, block);
		do{
			//print_gscmatrixd_math("V", V, numk, block);
			//test_GS(V, numk, block);
			Gram_Schmidt_gsl(V, numk, block, 0);
			Gram_Schmidt_gsl(V, numk, block, 0);
			//test_GS(V, numk, block);
			H_vector_map_complex(istart, Hdia, Hamil, column, V[nk], HV[nk], block);
			N = numk;
			LDA = N;
			for(mu=0; mu<numk; mu++) for(nu=mu; nu<numk; nu++){
				inn_r = 0; inn_i=0;
#pragma omp parallel for default(none) private(i) shared(HV, V, block, mu, nu) reduction(+:inn_r, inn_i)
				for(i=0; i<block; i++){
					inn_r += gsl_complex_in_r( HV[mu][i], V[nu][i] );
					inn_i += gsl_complex_in_i( HV[mu][i], V[nu][i] );
				}
				subspace[mu*numk+nu].r =  0;
				subspace[mu*numk+nu].i =  0;

				subspace[nu*numk+mu].r =  inn_r;
				subspace[nu*numk+mu].i =  inn_i;

				//matrix[nu][mu] = gsl_complex_rect(inn_r, inn_i);
				//matrix[mu][nu] = gsl_complex_rect(inn_r, -inn_i);
			}

			/////////////////////////////////
			//print_gscmatrixd_math("matrix", matrix, numk, numk);
			/////////////////////////////////

			zheev_(&JOBZ, &UPLO, &N, subspace, &LDA, W, WORK, &LWORK, RWORK, &INFO);

			if( (int)INFO != 0 ){
				printf("Davidson: Failed in zheev_, INFO = %d, Nk = %d, numk=%d, lancount=%d\n", (int)INFO, nk, numk, *lancount);	fflush(stdout);
				print_gscmatrixd("vectorin", vector_in0, nk*SDmax, block);
				print_vectord("eval", W, numk);
				print_vectord("res", res, SDmax);
				gsl_complex **matrix = mkgscmatrixd(SDmax*JDmax, SDmax*JDmax);

				for(mu=0; mu<numk; mu++) for(nu=mu; nu<numk; nu++){
					inn_r = 0; inn_i=0;
#pragma omp parallel for default(none) private(i) shared(HV, V, block, mu, nu) reduction(+:inn_r, inn_i)
					for(i=0; i<block; i++){
						inn_r += gsl_complex_in_r( HV[mu][i], V[nu][i] );
						inn_i += gsl_complex_in_i( HV[mu][i], V[nu][i] );
					}
					subspace[mu*numk+nu].r =  0;
					subspace[mu*numk+nu].i =  0;

					subspace[nu*numk+mu].r =  inn_r;
					subspace[nu*numk+mu].i =  inn_i;

					//matrix[nu][mu] = gsl_complex_rect(inn_r, inn_i);
					//matrix[mu][nu] = gsl_complex_rect(inn_r, -inn_i);
				}
				print_gscmatrixd("V", V, numk, block);
				print_gscmatrixd("HV", HV, numk, block);
				for(mu=0; mu<numk; mu++) for(nu=0; nu<numk; nu++){
					matrix[mu][nu] = gsl_complex_rect( subspace[nu*numk+mu].r, subspace[nu*numk+mu].i );
				}
				print_gscmatrixd("error zheev_", matrix, numk, numk);

				freegscmatrixd(matrix, SDmax*JDmax);
				exit(1);
			}

			/////////////////////////////////
			/*
			   for(mu=0; mu<numk; mu++){
			   for(nu=0; nu<numk; nu++)
			   printf("(%lf, %lf)\t", subspace[mu*numk+nu].r, subspace[mu*numk+nu].i );
			   printf("\n");
			   }
			   printf("\n");
			   for(mu=0; mu<numk; mu++){
			   for(nu=0; nu<numk; nu++)
			   matrix[mu][nu] = gsl_complex_rect( subspace[mu*numk+nu].r, subspace[mu*numk+nu].i );
			   }
			   print_gscmatrixd_math("eigenvectors", matrix, numk, numk);
			   */
			/////////////////////////////////

			//print_gsl_matrix_complex("subspace", subspace, numk, numk);

			if(nk==Nend-1) for(i=0; i<Nend; i++)
				values[i] = W[i];

			theta = W[0];
			for(nu=0; nu<numk; nu++){
#pragma omp parallel for default(none) private(i, mu) shared(u, V, block, nu, zero, numk, subspace)
				for(i=0; i<block; i++){
					u[nu][i] = zero;
					for(mu=0; mu<numk; mu++) {
						u[nu][i] = gsl_complex_add(
								u[nu][i],
								gsl_complex_mul(
									gsl_complex_rect(subspace[nu*numk+mu].r, subspace[nu*numk+mu].i),
									//gsl_matrix_complex_get(evec, mu, nu),
									V[mu][i]
									)
								);
					}
				}
			}
#pragma omp parallel default(none) private(i, mu) shared(nk, numk, u, r, HV, block, theta, zero, res, subspace, inn_r)
			{
#pragma omp for 
				for(i=0; i<block; i++){
					r[i] = zero;
					for(mu=0; mu<numk; mu++) {
						r[i] = gsl_complex_add(
								r[i],
								gsl_complex_mul(
									gsl_complex_rect(subspace[mu].r, subspace[mu].i),
									//gsl_matrix_complex_get(evec, mu, 0),
									HV[mu][i]
									)
								);
					}
					r[i] = gsl_complex_sub(
							r[i],
							gsl_complex_mul_real(u[0][i], theta)
							);
				}

			}
				inn_r = 0;
#pragma omp parallel for default(none) private(i) shared(block, r) reduction(+:inn_r)
				for(i=0; i<block; i++)
					(inn_r) += gsl_complex_abs2(r[i]);
				*res = inn_r;
			//print_gscmatrixd("u", u, numk, block);

#pragma omp parallel for default(none) private(i) shared(block, r, t, theta, Hdia)
			for(i=0; i<block; i++)
				t[i] = gsl_complex_div_real( r[i], theta - Hdia[i] );
			nk++;
			numk++;
			//printf("nk=%d, numk=%d\n", nk, numk);
			for(i=0; i<block; i++)	V[nk][i] = t[i];
		}while( nk<Nend && *res > Tolerance );

		for(nk=0; nk<Nvec; nk++){
			for(i=0; i<block; i++){
				vector_out[nk][i] = u[nk][i];
			}
			normalize_gsc(vector_out[nk], block);
		}

		//print_gscmatrixd("vector_out", vector_out, Nvec, block);

		inn_r = 0; inn_i = 0;
#pragma omp parallel for default(none) private(i) shared(block, vector_in0, vector_out) reduction(+:inn_r, inn_i)
		for(i=0; i<block; i++){
			inn_r += gsl_complex_in_r( vector_out[0][i], vector_in0[0][i] );
			inn_i += gsl_complex_in_i( vector_out[0][i], vector_in0[0][i] );
		}
		(*overlap) = sqrt(inn_r*inn_r + inn_i*inn_i);
		//freegscmatrixd(matrix, (JDmax));

		int ndeg=0;
		for(nk=1; nk<Nvec; nk++)
			if( values[nk] - values[0] < 10e-14)	ndeg++;

		if(ndeg){
			for(nk=1; nk<ndeg+1; nk++){
				inn_r=0; inn_i=0;
#pragma omp parallel for default(none) private(i) shared(block, vector_in0, vector_out, nk) reduction(+:inn_r, inn_i)
				for(i=0; i<block; i++){
					inn_r += gsl_complex_in_r( vector_out[nk][i], vector_in0[0][i] );
					inn_i += gsl_complex_in_i( vector_out[nk][i], vector_in0[0][i] );
				}
				(*overlap) = sqrt((*overlap)*(*overlap) + inn_r*inn_r + inn_i*inn_i);
			}
		}
		(*lancount) = (*lancount) + 1;

		freegscmatrixd(HV, Nend);
		freegscmatrixd(V, Nend+1);
		freegscmatrixd(u, Nend);
		free(t); free(r);
	}

	/*
	printf("ndeg=%d, theta = %.18lf\toverlap = %.18lf\n", ndeg, theta, *overlap);
	for(nk=0; nk<Nvec; nk++)
		printf("%lf\t", values[nk]);
	printf("\n");
	*/
	return theta;
}

void H_vector_map_complex(int *istart, double *Hdia, gsl_complex *Hamil, int *column, gsl_complex *vectorin, gsl_complex *vectorout, int block){
	int i, j;
	gsl_complex tempG;

#pragma omp parallel for default(none)	\
	private(i,j,tempG) shared(Hdia, istart, column, Hamil, vectorin, vectorout, block)
	for(j=0; j<block; j++){
		tempG = gsl_complex_mul_real( vectorin[j], Hdia[j] );

		for(i=istart[j]; i<istart[j+1]; i++){
			tempG = gsl_complex_add(
					tempG,
					gsl_complex_rect(
						gsl_complex_mul_r( &vectorin[column[i]], &Hamil[i] ),
						gsl_complex_mul_i( &vectorin[column[i]], &Hamil[i] )
					)
					//gsl_complex_mul( vectorin[column[i]], Hamil[i] )
				);
			//printf("Hamil[%d] = %lf + i %lf\n", i, GSL_REAL(Hamil[i]), GSL_IMAG(Hamil[i]));
			//printf("vectorout[%d] += vectorin[%d] (%lf + i %lf)*Hamil[i]\n", j, column[i], GSL_REAL(vectorin[column[i]]), GSL_IMAG(vectorin[column[i]]) );
		}
		vectorout[j] = tempG;
	}
	//print_gscvectord("vectorout", vectorout, block);
}//end of H_vector_map_complex

int print_Hdia(FILE *fp, int n, sectorinfo ginfo[], int gindex, typebasis **basis, PNSpara egv){
	int i, j, sp, k, refer, block, x;
	double element=0;
	Coulomb *DD;

	DD = mkvectorC(Ndd);
	init_dd(Uinter, J, DD);

	refer = ginfo[gindex].refer;
	block = ginfo[gindex].block;

	for(i=0; i<block; i++){
		j = refer + i;
		element = 0.;
		for(x=0; x<Ndd; x++){
			element += DD[x].V * (int)One(basis[j][DD[x].spa], DD[x].a) * (int)One(basis[j][DD[x].spb], DD[x].b);
			if( (int)One(basis[j][DD[x].spa], DD[x].a) && (int)One(basis[j][DD[x].spb], DD[x].b)){
				fprintf(fp, "Hamil[[%2d]][[%2d]]\t= Hamil[[%2d]][[%2d]]\t+ DD[[%2d]];\t\t(* %9.6lf *)\n", i+1, i+1, i+1, i+1, x+1, DD[x].V * (int)One(basis[j][DD[x].spa], DD[x].a) * (int)One(basis[j][DD[x].spb], DD[x].b));
			}
		}
		for(sp=0; sp<2; sp++) {
			for(k=0; k<nctot; k++){
				element += GSL_REAL(Ecluster[n*tNC+2*k+sp][n*tNC+2*k+sp]) * One(basis[j][sp], k);
				if(One(basis[j][sp], k))	fprintf(fp, "Hamil[[%2d]][[%2d]]\t= Hamil[[%2d]][[%2d]]\t+ Ecluster[[%2d]][[%2d]];\t\t(* %9.6lf *)\n", i+1, i+1, i+1, i+1, 2*k+sp+1, 2*k+sp+1, GSL_REAL(Ecluster[n*tNC+2*k+sp][n*tNC+2*k+sp]) * One(basis[j][sp], k));
			}
			for(k=nctot; k<nctot+nb; k++) element += egv.egbath[2*(k-nctot)+sp] * One(basis[j][sp], k); 
		}
	}
	free(DD);
	return 0;
}//end of compute_Hdia



int print_map(char *save_directory, int count, int n, typebasis **basis, int *istart, sectorinfo ginfo[], int gindex, PNSpara egv){
	int i, j, block, mapnumber;
	Coulomb *JJ;
	JJ = mkvectorC( Njj );
	init_others(J, JJ);


	//for(i=0; i<Njj; i++) printf("%d: %lf\n", i, JJ[i].V);

	gsl_complex element=zero;
	mapele *mapinfo;
	block = ginfo[gindex].block;

	if( block > 100 ){
		printf("print_map:: Too large dimension for gindex = %d (block = %d)\n", gindex, block);
		return 0;
	}

	mapinfo = mkvectormapele(ginfo[gindex].mapnumber);
	build_mapinfo(Ecluster, basis, mapinfo, istart, ginfo, gindex);
	mapnumber = ginfo[gindex].mapnumber;

	char para[1024];
	FILE *fp;
	sprintf(para, "%s/debug", save_directory);
	mkdirectory(para);
	sprintf(para, "%s/debug/H%03d_block%d_%dth.txt", save_directory, gindex, block, count);
	fp = fopen(para, "w");
	print_Hdia(fp, n, ginfo, gindex, basis, egv);

	for(j=0; j<block; j++) for(i=istart[j]; i<istart[j+1]; i++){
		if( mapinfo[i].code == 'c' ){
			element = gsl_complex_mul_real( Ecluster[n*tNC+ (int)(mapinfo[i].mu)][n*tNC+ (int)(mapinfo[i].p)], mapinfo[i].phase );
			fprintf(fp, "Hamil[[%2d]][[%2d]]\t= Hamil[[%2d]][[%2d]]\t+ Ecluster[[%2d]][[%2d]] * %2d;\t(* %9.6lf, %9.6lf *)\n", j+1, mapinfo[i].column+1, j+1, mapinfo[i].column+1, mapinfo[i].mu+1, mapinfo[i].p+1, mapinfo[i].phase, GSL_REAL(element), GSL_IMAG(element));
		}
		else if( mapinfo[i].code == 's' ){
			element = gsl_complex_mul_real( egv.hybrid[(int)(mapinfo[i].mu)][(int)(mapinfo[i].p)], mapinfo[i].phase );
			fprintf(fp, "Hamil[[%2d]][[%2d]]\t= Hamil[[%2d]][[%2d]]\t+ hybrid[[%2d]][[%2d]] * %2d;\t(* %9.6lf, %9.6lf *)\n", j+1, mapinfo[i].column+1, j+1, mapinfo[i].column+1, mapinfo[i].mu+1, mapinfo[i].p+1, mapinfo[i].phase, GSL_REAL(element), GSL_IMAG(element));
		}
		else if( mapinfo[i].code == 'h' ){
			element = gsl_complex_mul_real(
					gsl_complex_conjugate( egv.hybrid[(int)(mapinfo[i].mu)][(int)(mapinfo[i].p)] ),
					mapinfo[i].phase
					);
			fprintf(fp, "Hamil[[%2d]][[%2d]]\t= Hamil[[%2d]][[%2d]]\t+ chybrid[[%2d]][[%2d]] * %2d;\t(* %9.6lf, %9.6lf *)\n", j+1, mapinfo[i].column+1, j+1, mapinfo[i].column+1, mapinfo[i].mu+1, mapinfo[i].p+1, mapinfo[i].phase, GSL_REAL(element), GSL_IMAG(element));
		}
		else if( mapinfo[i].code == 'j' ){
			element = gsl_complex_rect( JJ[mapinfo[i].mu].V * mapinfo[i].phase, 0 );
			fprintf(fp, "Hamil[[%2d]][[%2d]]\t= Hamil[[%2d]][[%2d]]\t+ JJ[[%2d]] * %2d;\t(* %9.6lf, %9.6lf *)\n", j+1, mapinfo[i].column+1, j+1, mapinfo[i].column+1, mapinfo[i].mu +1, mapinfo[i].phase, GSL_REAL(element), GSL_IMAG(element));
		}
		else{
			fprintf(fp, "Error in read_map\n");
			fprintf(fp, "myrank=%d, gindex=%d, j=%d, i=%d, code=%c, mapnumber=%d, phase = %d\n", myrank, gindex, j, i, mapinfo[i].code, mapnumber, (int)mapinfo[i].phase);
			exit(1);
		}
	}
	fclose(fp);

	free(mapinfo);
	free(JJ);

	return 0;
}//end of print_map


int read_map(int n, typebasis **basis, int *istart, gsl_complex *Hamil, int *column, sectorinfo ginfo[], int gindex, PNSpara egv){
	int i, j, block, mapnumber;
	Coulomb *JJ;
	JJ = mkvectorC( Njj );
	init_others(J, JJ);

	//for(i=0; i<Njj; i++) printf("%d: %lf\n", i, JJ[i].V);

	gsl_complex element=zero;
	mapele *mapinfo;
	block = ginfo[gindex].block;

	mapinfo = mkvectormapele(ginfo[gindex].mapnumber);
	build_mapinfo(Ecluster, basis, mapinfo, istart, ginfo, gindex);
	mapnumber = ginfo[gindex].mapnumber;

#pragma omp parallel for default(none)	\
	private(i,element,j) shared(mapnumber,istart,Ecluster,nctot,Hamil,column,mapinfo,egv,block, myrank, gindex, JJ, n)
	for(j=0; j<block; j++) for(i=istart[j]; i<istart[j+1]; i++){
		if( mapinfo[i].code == 'c' )	element = gsl_complex_mul_real( Ecluster[n*tNC+ (int)(mapinfo[i].mu)][n*tNC+ (int)(mapinfo[i].p)], mapinfo[i].phase );
		else if( mapinfo[i].code == 's' ){
			element = gsl_complex_mul_real( egv.hybrid[(int)(mapinfo[i].mu)][(int)(mapinfo[i].p)], mapinfo[i].phase );
		}
		else if( mapinfo[i].code == 'h' ){
			element = gsl_complex_mul_real(
					gsl_complex_conjugate( egv.hybrid[(int)(mapinfo[i].mu)][(int)(mapinfo[i].p)] ),
					mapinfo[i].phase
				);
		}
		else if( mapinfo[i].code == 'j' ){
			element = gsl_complex_rect( JJ[mapinfo[i].mu].V * mapinfo[i].phase, 0 );
			//printf("i%d, %d, %lf+ I %lf\n", i, mapinfo[i].mu, GSL_REAL(element), GSL_IMAG(element));
		}
		else{
			printf("Error in read_map\n");
			printf("myrank=%d, gindex=%d, j=%d, i=%d, code=%c, mapnumber=%d, phase = %d\n", myrank, gindex, j, i, mapinfo[i].code, mapnumber, (int)mapinfo[i].phase);
			exit(1);
		}
		Hamil[i] = element;
		column[i] = mapinfo[i].column;
	}
	free(mapinfo);
	free(JJ);

	return 0;
}//end of read_map

int compute_Hdia(int n, double *Hdia, sectorinfo ginfo[], int gindex, typebasis **basis, PNSpara egv){
	int i, j, sp, k, refer, block, x;
	double element=0;
	Coulomb *DD;

	DD = mkvectorC(Ndd);
	init_dd(Uinter, J, DD);

	//for(i=0; i<Ndd; i++) printf("%d: %lf\n", i, DD[i].V);

	refer = ginfo[gindex].refer;
	block = ginfo[gindex].block;

#pragma omp parallel for default(none)	\
	private(i,sp,k,element,j,x) shared(block,refer,Hdia,Uinter,basis,egv,Ecluster,nctot,nb,DD,Ndd,n)
	for(i=0; i<block; i++){
		j = refer + i;
		element = 0.;
		for(x=0; x<Ndd; x++){
			element += DD[x].V * (int)One(basis[j][DD[x].spa], DD[x].a) * (int)One(basis[j][DD[x].spb], DD[x].b);
		}
		for(sp=0; sp<2; sp++) {
			for(k=0; k<nctot; k++)
				element += GSL_REAL(Ecluster[n*tNC+2*k+sp][n*tNC+2*k+sp]) * One(basis[j][sp], k);
			for(k=nctot; k<nctot+nb; k++) element += egv.egbath[2*(k-nctot)+sp] * One(basis[j][sp], k); 
		}
		Hdia[i] = element;
	}
	free(DD);
	return 0;
}//end of compute_Hdia


double MLanczos_Nvector(int *istart, double *Hdia, gsl_complex *Hamil, int *column, int NN, gsl_complex *vector_in0, gsl_complex *vector_out0, double *overlap, int *lancount, double *values, double *res){
	int i, j;
	gsl_complex **ortho, *input;
	double WORK[2*Nvec-1];
	double sum = 0;

	char JOB='V';
	integer resi, N=Nvec;
	doublereal ai[Nvec+1], bi[Nvec+1];
	doublereal subspace[Nvec*Nvec];

	input = mkgscvectord( NN );
	for(i=0; i<NN; i++){
		input[i] = vector_in0[i];
		vector_out0[i] = zero;
	}

	ortho = lanczos_vector_complex(istart, Hdia, Hamil, column, &i, input, ai, bi, NN, Nvec, 1); 
	for(i=0; i<Nvec-1; i++)	bi[i] = sqrt( bi[i+1] ); 

	//if( myrank == 10 )	for(i=0; i<Nvec; i++) printf("ai[%d] = %lf, bi[%d] = %lf\n", i, ai[i], i, bi[i]);

	dstev_(&JOB, &N, ai, bi, subspace, &N, WORK, &resi);


	for(i=0; i<Nvec; i++){
		values[i] = ai[i];
		for(j=0; j<NN; j++)
			vector_out0[j]
				= gsl_complex_add(
					vector_out0[j],
					gsl_complex_mul_real( ortho[i][j], subspace[i])
					);
	}

	sum = 0;
	for(i=0; i<NN; i++)	sum += gsl_complex_abs2( vector_out0[i] );
	sum = sqrt(sum);
	*overlap = 0;

	for(i=0; i<NN; i++){
		vector_out0[i] = gsl_complex_div_real( vector_out0[i], sum );
		*overlap += gsl_complex_inner( vector_in0[i], vector_out0[i] );
	}

	freegscmatrixd(ortho, Nvec);
	free( input );
	*lancount = *lancount + 1;

	//compute res
	gsl_complex *Hx = mkgscvectord(NN);
	H_vector_map_complex(istart, Hdia, Hamil, column, vector_out0, Hx, NN);
	sum = 0;
	for(i=0; i<NN; i++){
		Hx[i] = gsl_complex_sub( Hx[i], gsl_complex_mul_real( vector_out0[i], values[0] ) );
		sum += gsl_complex_abs2( Hx[i] );
	}
	*res = sum;

	free(Hx);

	return values[0];
}

gsl_complex ** lanczos_vector_complex(int *istart, double *Hdia, gsl_complex *Hamil, int *column, int *dimension, gsl_complex *p1, double ai[], double bi[], int grblock, int stdlaniter, int Vsave){
	int i, psub = 0;
	gsl_complex *Hp0, *Hp1, *p0, *p2;
	gsl_complex **ortho=NULL, **reortho=NULL;
	double p2inner = 0, temp;

	Hp0	= mkgscvectord( grblock );
	Hp1	= mkgscvectord( grblock );
	p0	= mkgscvectord( grblock );
	p2	= mkgscvectord( grblock );
	initvectord(ai,stdlaniter);
	initvectord(bi,stdlaniter);
	for(i=0; i<grblock; i++){
		Hp0[i]	= zero;
		Hp1[i]	= zero;
		p0[i]	= zero;
		p2[i]	= zero;
	}
	

	if( Vsave )
		ortho = mkgscmatrixd(Nvec, grblock);
	if( Ngram && !Vsave )
		reortho = mkgscmatrixd(Ngram, grblock);

	do{
		int Ngram_div = Ngram;
		if( Vsave ) for(i=0; i<grblock; i++)
			ortho[psub][i] = p1[i]; 
		if( Ngram && !Vsave ){
			for(i=0; i<grblock; i++)
				reortho[psub % Ngram_div][i] = p1[i]; 
			if( psub % Ngram_div == Ngram-1 ){
				Gram_Schmidt_gsl( reortho, Ngram, grblock, 0 );
				for(i=0; i<grblock; i++)
					p1[i] = reortho[psub % Ngram_div][i];
			}
		}

		H_vector_map_complex(istart, Hdia, Hamil, column, p1, Hp1, grblock);

		temp = sqrt( bi[psub] );
		p2inner = 0;
#pragma omp parallel for default(none) shared(grblock, Hp1, p0, p1, p2, temp) reduction(+:p2inner)
		for(i=0; i<grblock; i++){
			//p2[i] = gsl_complex_sub( Hp1[i], gsl_complex_mul_real( p0[i], temp ));
			p2[i] = gsl_complex_rect(
					GSL_REAL(Hp1[i]) - temp*GSL_REAL(p0[i]),
					GSL_IMAG(Hp1[i]) - temp*GSL_IMAG(p0[i]) );
			p2inner += gsl_complex_inner( p1[i], p2[i] );
		}

		ai[psub] = p2inner;

		p2inner = 0;
		temp = ai[psub];
#pragma omp parallel for default(none) shared(grblock, p1, p2, temp) reduction(+:p2inner)
		for(i=0; i<grblock; i++){
			//p2[i] = gsl_complex_sub( p2[i], gsl_complex_mul_real( p1[i], temp ));
			p2[i] = gsl_complex_rect(
					GSL_REAL(p2[i]) - temp*GSL_REAL(p1[i]),
					GSL_IMAG(p2[i]) - temp*GSL_IMAG(p1[i]) );
			p2inner += gsl_complex_abs2( p2[i] );
		}
		bi[psub+1] = p2inner;
		p2inner = sqrt( p2inner );
#pragma omp parallel for default(none) shared(grblock, p0, p1, p2, p2inner)
		for(i=0; i<grblock; i++){
			p0[i] = p1[i];
			p1[i] = gsl_complex_div_real( p2[i], p2inner ) ;
		}
		psub++;	
	}while( psub < stdlaniter && bi[psub] > 1e-14); 	
	*dimension = psub;
	if( grblock == 1 )	*dimension = 1;
	free(Hp0); free(Hp1); free(p0); free(p2);
	if( Ngram && !Vsave )	freegscmatrixd( reortho, Ngram );

	return ortho;
}//lanczos_vector_complex

int translate(Coulomb *DD, int Ndd){
	int x;
	int dic[tNC][3] = {
		{0}
	};
	for(x=0; x<tNC; x++){
		dic[x][0] = x;
		dic[x][1] = x/2;
		dic[x][2] = x%2;
	}
	for(x=0; x<Ndd; x++){
		DD[x].a = dic[DD[x].i][1]; DD[x].spa = dic[DD[x].i][2];
		DD[x].b = dic[DD[x].j][1]; DD[x].spb = dic[DD[x].j][2];
		DD[x].c = dic[DD[x].k][1]; DD[x].spc = dic[DD[x].k][2];
		DD[x].d = dic[DD[x].l][1]; DD[x].spd = dic[DD[x].l][2];
	}
	return 0;
}

int init_dd(double U, double J, Coulomb *DD){
	int Ndd;
	int x;

	double input_SOC[15][3] = {
		{1, 2, U - 7*J/3},
		{1, 3, U - 7*J/3},
		{1, 4, U - J    },
		{1, 5, U - 5*J/3},
		{1, 6, U - 8*J/3},
		{2, 3, U - J    },
		{2, 4, U - 7*J/3},
		{2, 5, U - 2*J  },
		{2, 6, U - 7*J/3},
		{3, 4, U - 7*J/3},
		{3, 5, U - 7*J/3},
		{3, 6, U - 2*J  },
		{4, 5, U - 8*J/3},
		{4, 6, U - 5*J/3},
		{5, 6, U - 4*J/3}
	};

	double input_eg[6][3] = {
		{1, 2, U},
		{1, 3, U - 3*J},
		{1, 4, U - 2*J},
		{2, 3, U - 2*J},
		{2, 4, U - 3*J},
		{3, 4, U},
	};

	double input_t2g[15][3] = {
		{1, 2, U},
		{1, 3, U - 3*J},
		{1, 4, U - 2*J},
		{1, 5, U-3*J},
		{1, 6, U-2*J},
		{2, 3, U - 2*J},
		{2, 4, U - 3*J},
		{2, 5, U-2*J},
		{2, 6, U-3*J},
		{3, 4, U},
		{3, 5, U-3*J},
		{3, 6, U-2*J},
		{4, 5, U -2*J},
		{4, 6, U -3*J},
		{5, 6, U},
	};
	if( Nonzero_SOC ){
		Ndd = 15;
		for(x=0; x<Ndd; x++){
			DD[x].i = (int) input_SOC[x][0]-1;
			DD[x].j = (int) input_SOC[x][1]-1;
			DD[x].k = (int) input_SOC[x][1]-1;
			DD[x].l = (int) input_SOC[x][0]-1;
			DD[x].V = input_SOC[x][2];
		}
	}
	else{
		if( NC == 2 ){
			Ndd = 6;
			for(x=0; x<Ndd; x++){
				DD[x].i = (int) input_eg[x][0]-1;
				DD[x].j = (int) input_eg[x][1]-1;
				DD[x].k = (int) input_eg[x][1]-1;
				DD[x].l = (int) input_eg[x][0]-1;
				DD[x].V = input_eg[x][2];
			}
		}
		else if( NC == 3 ){
			Ndd = 15;
			for(x=0; x<Ndd; x++){
				DD[x].i = (int) input_t2g[x][0]-1;
				DD[x].j = (int) input_t2g[x][1]-1;
				DD[x].k = (int) input_t2g[x][1]-1;
				DD[x].l = (int) input_t2g[x][0]-1;
				DD[x].V = input_t2g[x][2];
			}
		}
	}

	translate(DD, Ndd);
	return Ndd;
}

int init_others(double J, Coulomb *JJ){
	int Njj, x;
	
	double input_SOC[32][5] = {
		{1, 2, 1, 5, -(2*sqrt(2.)* J)/3},
		{1, 3, 1, 6, -(sqrt(2.)* J)/3},
		{1, 3, 2, 5, -sqrt(2./3.)* J},
		{1, 4, 2, 3, 4*J/3},
		{1, 4, 2, 6, -sqrt(2.)*J/3},
		{1, 4, 3, 5, -sqrt(2.)*J/3},
		{1, 4, 5, 6, -5*J/3},
		{1, 5, 1, 2, -2*sqrt(2.)* J/3},
		{1, 6, 1, 3, -sqrt(2.)* J/ 3},
		{1, 6, 2, 5, -J/sqrt(3.)},
		{2, 3, 1, 4, 4*J/3},
		{2, 3, 2, 6, -sqrt(2.)* J/3},
		{2, 3, 3, 5, -sqrt(2.)* J/3},
		{2, 3, 5, 6, 5*J/ 3},
		{2, 4, 3, 6, -sqrt(2./3.)* J},
		{2, 4, 4, 5, -sqrt(2.)* J/3},
		{2, 5, 1, 3, -sqrt(2./3.)* J},
		{2, 5, 1, 6, -J/sqrt(3.)},
		{2, 6, 1, 4, -sqrt(2.)* J/3},
		{2, 6, 2, 3, -sqrt(2.)* J/3},
		{2, 6, 3, 5, -2*J/3},
		{3, 4, 4, 6, -2*sqrt(2.)* J/3},
		{3, 5, 1, 4, -sqrt(2.)* J/3},
		{3, 5, 2, 3, -sqrt(2.)* J/3},
		{3, 5, 2, 6, -2*J/3},
		{3, 6, 2, 4, -sqrt(2./3.)* J},
		{3, 6, 4, 5, -J/sqrt(3.)},
		{4, 5, 2, 4, -sqrt(2.)* J/3},
		{4, 5, 3, 6, -J/sqrt(3.)},
		{4, 6, 3, 4, -2*sqrt(2.)* J/3},
		{5, 6, 1, 4, -5*J/3.},
		{5, 6, 2, 3, 5*J/3.}
	};

	double input_eg[4][5] = {
		{1, 2, 4, 3, J},
		{1, 4, 2, 3,  J},

		{2, 3, 4, 1, -J},
		{3, 4, 2, 1, J},
	};

	double input_t2g[12][5] = {
		{1, 2, 4, 3, J},
		{1, 2, 6, 5, J},
		{1, 4, 2, 3,  J},
		{1, 6, 5, 2, -J},
		{3, 4, 6, 5, J},
		{3, 6, 5, 4, -J},

		{2, 3, 4, 1, -J},
		{2, 5, 6, 1, -J},
		{3, 4, 2, 1, J},
		{4, 5, 6, 3, -J},
		{5, 6, 2, 1, J},
		{5, 6, 4, 3, J},
	};

	if( Nonzero_SOC ){
		Njj = 32;
		for(x=0; x<Njj; x++){
			JJ[x].i = (int) input_SOC[x][0]-1;
			JJ[x].j = (int) input_SOC[x][1]-1;
			JJ[x].k = (int) input_SOC[x][2]-1;
			JJ[x].l = (int) input_SOC[x][3]-1;
			JJ[x].V = input_SOC[x][4];
		}
	}
	else{
		if( NC == 2 ){
			Njj = 4;
			for(x=0; x<Njj; x++){
				JJ[x].i = (int) input_eg[x][0]-1;
				JJ[x].j = (int) input_eg[x][1]-1;
				JJ[x].k = (int) input_eg[x][2]-1;
				JJ[x].l = (int) input_eg[x][3]-1;
				JJ[x].V = input_eg[x][4];
			}
		}

		else if( NC == 3 ){
			Njj = 12;
			for(x=0; x<Njj; x++){
				JJ[x].i = (int) input_t2g[x][0]-1;
				JJ[x].j = (int) input_t2g[x][1]-1;
				JJ[x].k = (int) input_t2g[x][2]-1;
				JJ[x].l = (int) input_t2g[x][3]-1;
				JJ[x].V = input_t2g[x][4];
			}
		}
	}

	translate(JJ, Njj);
	return Njj;
}

int construct_nonzero_cluster(gsl_complex **Ecluster, int **nonzero){
	int mu, nu;
	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
		if( Bilinear_all || gsl_complex_abs2(Ecluster[mu][nu]) > Th_hop ){
			nonzero[mu][nu] = 1;
		}
		else	nonzero[mu][nu] = 0;
	}
	return 0;
}

int count_mapnumber(gsl_complex **Ecluster, typebasis **bitbasis, sectorinfo *ginfo, int gindex, int *Ndd, int *Njj){
	int row, mu, nu, sp1, sp2, refer, block, mapnumber=0, **nonzero, x;
	typebasis bin[2];
	Coulomb JJ[tNC*tNC*tNC*tNC];	//In fact, tNC! is enough
	Coulomb DD[tNC*(tNC-1)/2];

	*Ndd = init_dd(0, 0, DD);
	*Njj = init_others(0, JJ);

	nonzero = mkmatrixi(tNC, tNC);
	construct_nonzero_cluster(Ecluster, nonzero);

	refer = ginfo[gindex].refer;
	block = ginfo[gindex].block;
	//printf("gindex = %d, refer = %d, block = %d\n", gindex, refer, block);
	//print_gbitbasis(Ns, bitbasis, ginfo, gindex);

	//printf("<------------------ %dth, ptl = %d ---------------------->\n", gindex, ginfo[gindex].ptl);
	for(row=0; row<block; row++){
		for(x=0; x<*Njj; x++){
			if( One(bitbasis[refer+row][JJ[x].spd], JJ[x].d ) && One(bitbasis[refer+row][JJ[x].spc], JJ[x].c ) ){
				bin[0] = bitbasis[refer+row][0]; bin[1] = bitbasis[refer+row][1];

				bin[JJ[x].spd] = Turnoff(bin[JJ[x].spd], JJ[x].d );	//phi_d
				bin[JJ[x].spc] = Turnoff(bin[JJ[x].spc], JJ[x].c );	//phi_c
				if( Zero(bin[JJ[x].spb], JJ[x].b) && Zero(bin[JJ[x].spa], JJ[x].a) ){
					mapnumber++;
				}
			}
		}
		for(sp1=0; sp1<2; sp1++){
			int spin_init, spin_final;
			if( Nonzero_SOC ){
				spin_init = 0;
				spin_final = 2;
			}
			else{
				spin_init = sp1;
				spin_final = sp1+1;
			}
			for(sp2=spin_init; sp2<spin_final; sp2++) {
				/******************************************** Diagonal ******************************************************/
				// diagonal part is always relevant, so it can be omit.
				/************************************************************************************************************/
				for(nu=0; nu<nctot; nu++) for(mu=0; mu<nctot; mu++){//intracluster hopping
					if( Zero( bitbasis[refer+row][sp1], mu ) &&  One( bitbasis[refer+row][sp2], nu ) ){
						if( nonzero[2*mu+sp1][2*nu+sp2] ){
							mapnumber++;
							//printf("F %d T%d, %d of %d -> %d of %d\n", row, element.column, sp2, nu, sp1, mu);
						}
					}
				}
				for(mu=nctot; mu<Ns; mu++) for(nu=0; nu<nctot; nu++){//hybridization
					if( Zero( bitbasis[refer+row][sp1], mu ) &&  One( bitbasis[refer+row][sp2], nu ) ){
						mapnumber++;
						//printf("F %d T%d, %d of %d -> %d of %d\n", row, element.column, sp2, nu, sp1, mu);
					}
					else if( One( bitbasis[refer+row][sp1], mu ) &&  Zero( bitbasis[refer+row][sp2], nu ) ){
						mapnumber++;
						//printf("F %d T%d, %d of %d -> %d of %d\n", row, element.column, sp2, nu, sp1, mu);
					}
				}
			}
		}
	}
	ginfo[gindex].mapnumber = mapnumber;
	freematrixi(nonzero, 2*NC);

	return 0;
}

int build_mapinfo(gsl_complex **Ecluster, typebasis **bitbasis, mapele *mapinfo, int *istart, sectorinfo *ginfo, int gindex){
	int i, *table, row, mu, nu, phase, sp1, sp2, refer, block, mapnumber=0, index, **nonzero, x, Njj;
	typebasis bin[2];
	mapele element;
	Coulomb JJ[tNC*tNC*tNC*tNC];	//In fact, tNC! is enough

	Njj = init_others(0, JJ);

	table = mkvectori( Powns2 );
	for(i=0; i<Powns2; i++){
		index = (bitbasis[i][0]<<Ns) + bitbasis[i][1];
		table[index] = i;
	}

	nonzero = mkmatrixi(tNC, tNC);
	construct_nonzero_cluster(Ecluster, nonzero);

	refer = ginfo[gindex].refer;
	block = ginfo[gindex].block;
	//printf("gindex = %d, refer = %d, block = %d\n", gindex, refer, block);
	//print_gbitbasis(Ns, bitbasis, ginfo, gindex);

	//printf("<------------------ %dth, ptl = %d ---------------------->\n", gindex, ginfo[gindex].ptl);
	for(row=0; row<block; row++){
		istart[row] = mapnumber;
		for(x=0; x<Njj; x++){
			if( One(bitbasis[refer+row][JJ[x].spd], JJ[x].d ) && One(bitbasis[refer+row][JJ[x].spc], JJ[x].c ) ){
				bin[0] = bitbasis[refer+row][0]; bin[1] = bitbasis[refer+row][1];

				phase  = permu( bin, JJ[x].d, JJ[x].spd);	bin[JJ[x].spd] = Turnoff(bin[JJ[x].spd], JJ[x].d );	//phi_d
				phase *= permu( bin, JJ[x].c, JJ[x].spc);	bin[JJ[x].spc] = Turnoff(bin[JJ[x].spc], JJ[x].c );	//phi_c
				if( Zero(bin[JJ[x].spb], JJ[x].b) && Zero(bin[JJ[x].spa], JJ[x].a) ){
					phase *= permu( bin, JJ[x].b, JJ[x].spb);	bin[JJ[x].spb] = Turnon(bin[JJ[x].spb], JJ[x].b );	//phi_b
					phase *= permu( bin, JJ[x].a, JJ[x].spa);	bin[JJ[x].spa] = Turnon(bin[JJ[x].spa], JJ[x].a );	//phi_a
					element.code	= 'j';
					index = (bin[0]<<Ns) + bin[1];
					element.column	= table[index]-refer;
					element.mu	= x;
					element.p	= -1;	//dummy: shouldn't be used
					element.phase	= phase;
					mapinfo[mapnumber] = element;
					/*
					printf("(%d) %c: %dth phase %d (%d,%d)(%d,%d)(%d,%d)(%d,%d)\n", mapnumber, element.code, x, phase, JJ[x].a, JJ[x].spa, JJ[x].b, JJ[x].spb, JJ[x].c, JJ[x].spc, JJ[x].d, JJ[x].spd);
					print_bit(Ns, bitbasis[refer+row][0]);
					print_bit(Ns, bitbasis[refer+row][1]);
					printf("\n");
					print_bit(Ns, bitbasis[refer+element.column][0]);
					print_bit(Ns, bitbasis[refer+element.column][1]);
					printf("\n");
					*/
					mapnumber++;
				}
			}
		}
		for(sp1=0; sp1<2; sp1++){ 
			int spin_init, spin_final;
			if( Nonzero_SOC ){
				spin_init = 0;
				spin_final = 2;
			}
			else{
				spin_init = sp1;
				spin_final = sp1+1;
			}
			for(sp2=spin_init; sp2<spin_final; sp2++) {
				/******************************************** Diagonal ******************************************************/
				// diagonal part is always relevant, so it can be omit.
				/************************************************************************************************************/

				//print_bit(Ns, bitbasis[refer+row][sp]);
				//print_bit(Ns, bitbasis[refer+column][sp]);
				//printf("(%d,%d)\t%d\t%d\n", row, column, bitbasis[refer+row][sp], bitbasis[refer+column][sp]);
				for(nu=0; nu<nctot; nu++) for(mu=0; mu<nctot; mu++){//intracluster hopping
					if( Zero( bitbasis[refer+row][sp1], mu ) &&  One( bitbasis[refer+row][sp2], nu ) ){
						if( nonzero[2*mu+sp1][2*nu+sp2] ){
							bin[0] = bitbasis[refer+row][0];
							bin[1] = bitbasis[refer+row][1];
							phase  = permu( bin, nu, sp2);	bin[sp2] = Turnoff( bin[sp2], nu );
							phase *= permu( bin, mu, sp1);	bin[sp1] = Turnon( bin[sp1], mu );
							index = (bin[0]<<Ns) + bin[1];
							element.code	= 'c';
							element.column	= table[index]-refer;
							element.mu	= 2*mu+sp1;
							element.p	= 2*nu+sp2;
							element.phase	= phase;
							mapinfo[mapnumber] = element;
							mapnumber++;
							//printf("F %d T%d, %d of %d -> %d of %d\n", row, element.column, sp2, nu, sp1, mu);
						}
					}
				}
				for(mu=nctot; mu<Ns; mu++) for(nu=0; nu<nctot; nu++){//hybridization
					if( Zero( bitbasis[refer+row][sp1], mu ) &&  One( bitbasis[refer+row][sp2], nu ) ){
						bin[0] = bitbasis[refer+row][0];
						bin[1] = bitbasis[refer+row][1];
						phase  = permu( bin, nu, sp2);	bin[sp2] = Turnoff( bin[sp2], nu );
						phase *= permu( bin, mu, sp1);	bin[sp1] = Turnon( bin[sp1], mu );
						index = (bin[0]<<Ns) + bin[1];
						element.code	= 'h';
						element.column	= table[index]-refer;
						element.mu	= 2*nu+sp2;
						element.p	= 2*(mu-nctot)+sp1;
						element.phase	= phase;
						mapinfo[mapnumber] = element;
						mapnumber++;
						//printf("F %d T%d, %d of %d -> %d of %d\n", row, element.column, sp2, nu, sp1, mu);
					}
					else if( One( bitbasis[refer+row][sp1], mu ) &&  Zero( bitbasis[refer+row][sp2], nu ) ){
						bin[0] = bitbasis[refer+row][0];
						bin[1] = bitbasis[refer+row][1];
						phase  = permu( bin, mu, sp1);	bin[sp1] = Turnoff( bin[sp1], mu );
						phase *= permu( bin, nu, sp2);	bin[sp2] = Turnon( bin[sp2], nu );
						index = (bin[0]<<Ns) + bin[1];
						element.code	= 's';
						element.column	= table[index]-refer;
						element.mu	= 2*nu+sp2;
						element.p	= 2*(mu-nctot)+sp1;
						element.phase	= phase;
						mapinfo[mapnumber] = element;
						mapnumber++;
						//printf("F %d T%d, %d of %d -> %d of %d\n", row, element.column, sp2, nu, sp1, mu);
					}
				}
			}
		}
	}
	istart[block] = mapnumber;
	//printf("mapnumber = %d\n", mapnumber);

	sort_map_column(mapinfo, istart, ginfo, gindex);

	freematrixi(nonzero, tNC);
	free(table);

	return 0;
}

int construct_mapnumber(gsl_complex **Ecluster, typebasis **bitbasis, sectorinfo *ginfo, int init, int final, int *Ndd, int *Njj){
	int gindex;

	if(final > Blocks){
		printf("final must be less than Blocks!\n");
		exit(1);
	}
	//perform_old(Ns);
	//print_bitbasis(Ns, bitbasis);
	//print_ginfo(ginfo);

	for(gindex=init; gindex<final; gindex++){
		count_mapnumber(Ecluster, bitbasis, ginfo, gindex, Ndd, Njj);
	}
	return 0;
}



Coulomb *mkvectorC(int NN){
	Coulomb *vector;
	vector = (Coulomb *) malloc( NN * sizeof(Coulomb) );
	return vector;
}

mapele *mkvectormapele(int mapnumber){
	mapele *mapvector;
	mapvector = (mapele *) malloc ( mapnumber * sizeof( mapele ) );
	return mapvector;
}
