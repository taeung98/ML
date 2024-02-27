#include <mpi.h>
#include "matrix.h"
#include "lattice.h"
#include "ls.c"

extern int myrank, size;
extern double DELTA, MU, JTD;
extern gsl_complex zero;
extern MPI_Datatype MPI_GSL_COMPLEX;



double compute_hopmatrix_all(gsl_complex *hopmatrix, Quad *quadrature){
	int i, j, kx, ky, kz, mu, nu, line=0, num=NU*NU*Nintx*Ninty*Nintz;
	double Efermi;

	char para[1024];
	sprintf(para, "%s/%s", Path_hop, WANNIER);
	Latt *data;
	data = read_lattice(para, &line, &Efermi);
	init_quadrature(quadrature);

	FILE *fhop;
	sprintf(para, "%s/mdisp", Path_hop);
	if( myrank == 0 )	mkdirectory(para);
	char parahop[1024];
	int *number, *offset;
	sprintf(parahop, "%s/line%d_Nintx%d_Ninty%d_Nintz%d.hop", para, line, Nintx, Ninty, Nintz);
	fhop = fopen(parahop, "rb");

	if( fhop == NULL ){
		int mystart, myend;
		gsl_complex ex;
		mystart = myrank * (Nintx/size);
		myend = (myrank+1) * (Nintx/size);
		if( myrank == size )	myend = Nintx;
		if( Nintx%size != 0 ){
			printf("Nintx %d should be multiple of size %d!\n", Nintx, size);
			exit(1);
		}
		number = mkvectori( size );
		offset = mkvectori( size );

		for(kx=mystart; kx<myend; kx++) for(ky=0; ky<Ninty; ky++) for(kz=0; kz<Nintz; kz++) for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++) 
			hopmatrix(kx,ky,kz,mu,nu) = zero;

		for(i=0; i<line; i++) for(kx=mystart; kx<myend; kx++){
#pragma omp parallel for default(none)	\
		private(ky,kz,ex) shared(i,kx,data,hopmatrix,quadrature)
			for(ky=0; ky<Ninty; ky++) for(kz=0; kz<Nintz; kz++){
				ex = gsl_complex_rect(
						cos( (quadrature->position[0])[kx]*data[i].da1+ (quadrature->position[1])[ky]*data[i].da2 + (quadrature->position[2])[kz]*data[i].da3 ),
						sin( (quadrature->position[0])[kx]*data[i].da1+ (quadrature->position[1])[ky]*data[i].da2 + (quadrature->position[2])[kz]*data[i].da3 )
						);
				hopmatrix(kx,ky,kz,data[i].mu,data[i].nu)
					= gsl_complex_rect(
							GSL_REAL(hopmatrix(kx,ky,kz,data[i].mu,data[i].nu)) + gsl_complex_mul_r(&(data[i].tmunu), &ex ),
							GSL_IMAG(hopmatrix(kx,ky,kz,data[i].mu,data[i].nu)) + gsl_complex_mul_i(&(data[i].tmunu), &ex )
							);
			}
		}
		j=0;
		for(kx=mystart; kx<myend; kx++) for(ky=0; ky<Ninty; ky++) for(kz=0; kz<Nintz; kz++) for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++) {
			j++;
		}
		MPI_Allgather(&j, 1, MPI_INT, number, 1, MPI_INT, MPI_COMM_WORLD);
		offset[0] = 0;
		for(i=1; i<size; i++){
			offset[i] = offset[i-1] + number[i-1];
		}
		MPI_Allgatherv(MPI_IN_PLACE, number[myrank], MPI_GSL_COMPLEX, hopmatrix, number, offset, MPI_GSL_COMPLEX, MPI_COMM_WORLD);
		
		if( myrank == 0 ){
			fhop = fopen(parahop, "wb");
			fwrite(hopmatrix, sizeof(gsl_complex), num, fhop);
			fclose(fhop);
			for(i=1; i<size; i++)
				printf("i=%d, %d   ", i, number[i]);
			printf("\n");
		}
		free(number);
		free(offset);
	}
	else{
		if( myrank == 0 ){
			fread(hopmatrix, sizeof(gsl_complex), num, fhop);
		}
		fclose(fhop);
		MPI_Bcast(hopmatrix, num, MPI_GSL_COMPLEX, 0, MPI_COMM_WORLD);
	}
	free(data);

	int n;
	double distortion[Ni][tNC];
	init_distortion(distortion);
	for(kx=0; kx<Nintx; kx++) for(ky=0; ky<Ninty; ky++) for(kz=0; kz<Nintz; kz++) for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) 
		GSL_REAL(hopmatrix(kx,ky,kz,(n*tNC+mu),(n*tNC+mu))) = GSL_REAL(hopmatrix(kx,ky,kz,(n*tNC+mu),(n*tNC+mu))) + distortion[n][mu];

	return Efermi;
}

int rotate_hopmatrix_nint(gsl_complex *hopmatrix, gsl_complex **transform, int Nint){
	int kx, ky, kz, mu, nu, Omp_num = omp_get_max_threads();
	gsl_complex ***original = mkgsctritensord(Omp_num, NU, NU);
	gsl_complex ***rotated = mkgsctritensord(Omp_num, NU, NU);

	for(kx=0; kx<Nint; kx++){
#pragma omp parallel default(none)	\
		private(ky,kz,mu,nu) shared(kx,hopmatrix,transform, original, rotated, Nint)
		{
			int mythread = omp_get_thread_num();
		#pragma omp for
			for(ky=0; ky<Nint; ky++) for(kz=0; kz<Nint; kz++){
				for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
					original[mythread][mu][nu] = Evaluate_nint(hopmatrix, kx, ky, kz, mu, nu);
				iunitary(original[mythread], rotated[mythread], transform, NU);
				for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
					Evaluate_nint(hopmatrix, kx, ky, kz, mu, nu) = rotated[mythread][mu][nu];
			}
		}
	}

	freegsctritensord(original, Omp_num, NU);
	freegsctritensord(rotated, Omp_num, NU);
	return 0;
}

int rotate_hopmatrix(gsl_complex *hopmatrix, gsl_complex **transform){
	int kx, ky, kz, mu, nu, Omp_num = omp_get_max_threads();
	gsl_complex ***original = mkgsctritensord(Omp_num, NU, NU);
	gsl_complex ***rotated = mkgsctritensord(Omp_num, NU, NU);

	for(kx=0; kx<Nintx; kx++){
#pragma omp parallel default(none)	\
		private(ky,kz,mu,nu) shared(kx,hopmatrix,transform, original, rotated)
		{
			int mythread = omp_get_thread_num();
		#pragma omp for
			for(ky=0; ky<Ninty; ky++) for(kz=0; kz<Nintz; kz++){
				for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
					original[mythread][mu][nu] = hopmatrix(kx, ky, kz, mu, nu);
				iunitary(original[mythread], rotated[mythread], transform, NU);
				for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
					hopmatrix(kx, ky, kz, mu, nu) = rotated[mythread][mu][nu];
			}
		}
	}

	freegsctritensord(original, Omp_num, NU);
	freegsctritensord(rotated, Omp_num, NU);
	return 0;
}

int compute_lattice_vectors(double bvec[3][3], double rvec[3][3]){
	int i, j;
	double bvect[3][3] = {
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1}
	};
	for(i=0; i<3; i++) for(j=0; j<3; j++)	bvec[i][j] = bvect[i][j];
	FILE *fb = fopen("BVECTORS", "r");
	if( fb != NULL ){
		for(i=0; i<3; i++) for(j=0; j<3; j++)	fscanf(fb, "%lf", &bvec[i][j]);
		fclose(fb);
	}
	else{
		printf("BVECTORS are not given: Assuming simple cubic lattice...");
	}
	double volume=0, cross[3];

	cross_product(bvec[0], bvec[1], cross);
	for(i=0; i<3; i++)	volume += cross[i] * bvec[2][i];

	double normfactor = volume/2/PI;
	printf("volume=%lf, normfactor = %lf\n", volume, normfactor);
	for(i=0; i<3; i++){
		cross_product(bvec[(i+1)%3], bvec[(i+2)%3], cross);
		for(j=0; j<3; j++)	rvec[i][j] = cross[j]/normfactor;
	}

	for(i=0; i<3; i++){
		printf("r[%d]: ", i);
		for(j=0; j<3; j++){
			printf("%lf\t", rvec[i][j]);
		}
		printf("\n");
	}

	return 0;
}

int compute_dh_all(gsl_complex *hop, gsl_complex *dh, gsl_complex *kin, Quad *quadrature, int Nint){
	int i, j, kx, ky, kz, mu, nu, line, num=NU*NU*Nint*Nint*Nint;
	double Efermi;

	char para[1024];
	sprintf(para, "%s/%s", Path_hop, WANNIER);
	Latt *data;
	data = read_lattice(para, &line, &Efermi);
	init_quadrature_nint(quadrature, Nint);

	double rvec[3][3], bvec[3][3];
	compute_lattice_vectors(bvec, rvec);

	FILE *fhop, *fdh, *fkin;
	if( myrank == 0 )	mkdirectory("./mdisp");
	char parahop[1024], paradh[1024], parakin[1024];
	int *number, *offset;
	gsl_complex mi = gsl_complex_rect(0, -1);

	sprintf(parahop, "./mdisp/line%d_Nint%d.hop", line, Nint); fhop = fopen(parahop, "rb");
	sprintf(paradh , "./mdisp/line%d_Nint%d.dh" , line, Nint); fdh  = fopen(paradh , "rb");
	sprintf(parakin, "./mdisp/line%d_Nint%d.kin", line, Nint); fkin = fopen(parakin, "rb");

	if( fhop == NULL || fdh == NULL || fkin == NULL ){
		int mystart, myend, myindex;
		gsl_complex ex, ex1, ex2;
		mystart = myrank * (Nint/size);
		myend = (myrank+1) * (Nint/size);
		if( myrank == size )	myend = Nint;
		if( Nint%size != 0 ){
			printf("Nint %d should be multiple of size %d!\n", Nint, size);
			exit(1);
		}
		number = mkvectori( size );
		offset = mkvectori( size );

		for(i=0; i<num; i++){
			hop[i] = zero;
			dh[i] = zero;
			kin[i] = zero;
		}

		for(i=0; i<line; i++) for(kx=mystart; kx<myend; kx++){
#pragma omp parallel for default(none)	\
		private(ky,kz,ex,ex1,ex2,myindex) shared(i,kx,data,hop,dh,kin,line,quadrature,rvec,mi,Nint)
			for(ky=0; ky<Nint; ky++) for(kz=0; kz<Nint; kz++){
				ex = gsl_complex_rect(
						cos( (quadrature->position[0])[kx]*data[i].da1+ (quadrature->position[1])[ky]*data[i].da2 + (quadrature->position[2])[kz]*data[i].da3 ),
						sin( (quadrature->position[0])[kx]*data[i].da1+ (quadrature->position[1])[ky]*data[i].da2 + (quadrature->position[2])[kz]*data[i].da3 )
						);
				ex1 = gsl_complex_mul_real(
						gsl_complex_mul( ex, mi ), data[i].da1*rvec[0][0] + data[i].da2*rvec[1][0] + data[i].da3*rvec[2][0] 
					);
				ex2 = gsl_complex_mul_real( ex,	- pow(     data[i].da1*rvec[0][0] + data[i].da2*rvec[1][0] + data[i].da3*rvec[2][0], 2 ));

				myindex = kx*Nint*Nint*NU*NU + ky*Nint*NU*NU + kz*NU*NU + data[i].mu*NU + data[i].nu;
				hop[myindex]
					= gsl_complex_rect(
							GSL_REAL(hop[myindex]) + gsl_complex_mul_r(&(data[i].tmunu), &ex ),
							GSL_IMAG(hop[myindex]) + gsl_complex_mul_i(&(data[i].tmunu), &ex )
							);

				dh[myindex]
					= gsl_complex_rect(
							GSL_REAL(dh[myindex]) + gsl_complex_mul_r(&(data[i].tmunu), &ex1 ),
							GSL_IMAG(dh[myindex]) + gsl_complex_mul_i(&(data[i].tmunu), &ex1 )
							);

				kin[myindex]
					= gsl_complex_rect(
							GSL_REAL(kin[myindex]) + gsl_complex_mul_r(&(data[i].tmunu), &ex2 ),
							GSL_IMAG(kin[myindex]) + gsl_complex_mul_i(&(data[i].tmunu), &ex2 )
							);
			}
		}
		j=0;
		for(kx=mystart; kx<myend; kx++) for(ky=0; ky<Nint; ky++) for(kz=0; kz<Nint; kz++) for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++) {
			j++;
		}
		MPI_Allgather(&j, 1, MPI_INT, number, 1, MPI_INT, MPI_COMM_WORLD);
		offset[0] = 0;
		for(i=1; i<size; i++){
			offset[i] = offset[i-1] + number[i-1];
		}
		if( myrank==0 ) MPI_Allgatherv(MPI_IN_PLACE, number[myrank], MPI_GSL_COMPLEX, hop, number, offset, MPI_GSL_COMPLEX, MPI_COMM_WORLD);
		else		MPI_Allgatherv(hop, number[myrank], MPI_GSL_COMPLEX, hop, number, offset, MPI_GSL_COMPLEX, MPI_COMM_WORLD);

		if( myrank==0 ) MPI_Allgatherv(MPI_IN_PLACE, number[myrank], MPI_GSL_COMPLEX, dh, number, offset, MPI_GSL_COMPLEX, MPI_COMM_WORLD);
		else		MPI_Allgatherv(dh, number[myrank], MPI_GSL_COMPLEX, dh, number, offset, MPI_GSL_COMPLEX, MPI_COMM_WORLD);

		if( myrank==0 ) MPI_Allgatherv(MPI_IN_PLACE, number[myrank], MPI_GSL_COMPLEX, kin, number, offset, MPI_GSL_COMPLEX, MPI_COMM_WORLD);
		else		MPI_Allgatherv(kin, number[myrank], MPI_GSL_COMPLEX, kin, number, offset, MPI_GSL_COMPLEX, MPI_COMM_WORLD);
		
		if( myrank == 0 ){
			fhop = fopen(parahop, "wb");
			fdh = fopen(paradh, "wb");
			fkin = fopen(parakin, "wb");
			fwrite(hop, sizeof(gsl_complex), num, fhop);
			fwrite(dh, sizeof(gsl_complex), num, fdh);
			fwrite(kin, sizeof(gsl_complex), num, fkin);
			fclose(fhop);
			fclose(fdh);
			fclose(fkin);
			for(i=1; i<size; i++)
				printf("i=%d, %d   ", i, number[i]);
			printf("\n");
			fflush(stdout);
		}
		free(number);
		free(offset);
	}
	else{
		if( myrank == 0 ){
			fread(hop, sizeof(gsl_complex), num, fhop);
			fread(dh, sizeof(gsl_complex), num, fdh);
			fread(kin, sizeof(gsl_complex), num, fkin);
		}
		fclose(fhop);
		fclose(fdh);
		fclose(fkin);
		MPI_Bcast(hop, num, MPI_GSL_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Bcast(dh, num, MPI_GSL_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Bcast(kin, num, MPI_GSL_COMPLEX, 0, MPI_COMM_WORLD);
	}
	free(data);

	int n;
	double distortion[Ni][tNC];
	init_distortion(distortion);
	for(kx=0; kx<Nint; kx++) for(ky=0; ky<Nint; ky++) for(kz=0; kz<Nint; kz++) for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) {
		GSL_REAL(Evaluate_nint(hop, kx,ky,kz,(n*tNC+mu),(n*tNC+mu))) = GSL_REAL(Evaluate_nint(hop, kx,ky,kz,(n*tNC+mu),(n*tNC+mu))) + distortion[n][mu];
	}

	return 0;
}


