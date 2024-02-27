#ifndef _Lattice_
#define _Lattice_

#define Hop	1.

#define Path_hop	"./"
#define Path_save	"./"
#define Bilinear_all	1
#define Post_process	1
#define Debug		0
#define Save_davidson	1

#define Ns	12		//number of orbitals: bath + impurity (NB+NC)	12
#define NC	3		//number of correlated (impurity) orbitals
#define tNC	6		//2*NC
#define NB	9		//number of bath orbitals			9
#define Np	234		//number of independent bath parameters (onsite energies, hybridization strengths)	(2*nc)*(2*nb)*2+2*nb, 156 for nb=6 234 for nb=9
#define Np_zeroSOC	63	//number of independent bath parameters (no SOC)	42 for nb=6 , 63 for nb=9   (2*nc+1)*nb
#define Ni	8		//number of impurities (# of Ru)
#define NU	48		//tNC * Ni



#define Spmax	2

#define Tolgreen	14		//ground state accuracy 1e-14
#define N_pre		20		//Lanczos steps prior to Davidson
#define Upper		600		//# of continued fraction coefficients
#define Ngram		0
#define SPFACTOR	0.5
#define PARAMAGNETIC	0
#define SpinFlip	0

#define control_para	Uinter
#define control_char	'u'

#define JDmax		12
#define SDmax		64
#define Dmax		SDmax
#define Nvec		8
#define ind(n)		n

#define Nintx		64
#define Ninty		64
#define Nintz		64
#define Maxmin		2000
#define WANNIER		"inputs/lattice.txt"

#define Th_hop		1e-10

#define Area		(8*PI*PI*PI)
#define PI		3.14159265358979323846264338
#define delta_w		(double)(PI/beta)
//#define wn		((double)( (2*i+1) * PI/beta ))
#define wn		Matsu[i]
#define diff		1e-8						//dx for numerical derivation.

#define gsl_complex_inner(z1, z2)	(GSL_REAL(z1)*GSL_REAL(z2)+GSL_IMAG(z1)*GSL_IMAG(z2))
#define gsl_complex_in_r(z1, z2)	(GSL_REAL(z1)*GSL_REAL(z2)+GSL_IMAG(z1)*GSL_IMAG(z2))
#define gsl_complex_in_i(z1, z2)	(GSL_REAL(z1)*GSL_IMAG(z2)-GSL_IMAG(z1)*GSL_REAL(z2))

#define hopmatrix(kx,ky,kz,mu,nu)	hopmatrix[(kx)*Ninty*Nintz*NU*NU + (ky)*Nintz*NU*NU + (kz)*NU*NU + (mu)*NU + nu]
#define Evaluate(hop, kx,ky,kz,mu,nu)	hop[(kx)*Ninty*Nintz*NU*NU + (ky)*Nintz*NU*NU + (kz)*NU*NU + (mu)*NU + nu]
#define Evaluate_nint(hop, kx,ky,kz,mu,nu)	hop[(kx)*Nint*Nint*NU*NU + (ky)*Nint*NU*NU + (kz)*NU*NU + (mu)*NU + nu]

#define hyb(i, j)			gsl_complex_rect(p[tnb + 2*tnb*i +2*j], p[tnb + 2*tnb*i +2*j+1])

#define Turnoff(a, n)	((a)^(0x01<<(n)))
#define Turnon(a, n)	((a)|(0x01<<(n)))
#define Zero(a,n)	((~(a)>>(n))&0x01)
#define One(a,n)	(((a)>>(n))&0x01)
#define Basis_One	1u
//#define prefix		"./work"

#include <gsl/gsl_complex.h>

typedef unsigned int	typebasis;
typedef struct{
	int i;
	int j;
	int k;
	int l;

	int a;
	int b;
	int c;
	int d;

	int spa;
	int spb;
	int spc;
	int spd;

	double V;
} Coulomb;

typedef struct{
	double *egbath;
	gsl_complex **hybrid;
} PNSpara;

typedef struct{
	double *position[3];
	double *weight[3];
	int Nint[3];
} Quad;

typedef struct{
	int da1;
	int da2;
	int da3;
	int mu;
	int nu;
	gsl_complex tmunu;
} Latt;

typedef struct{
	int direction;
	int dimension;
	int mu[2];
	int sp[2];
	gsl_complex co[2];
	char type;
	double *aibi;
} Glist;

typedef struct{
	int column;
	char code;
	unsigned char mu;
	unsigned char p;
	signed char phase;
} mapele;

typedef struct{
	//double energy;
	double energy[Dmax];
	int refer;
	int block;
	int ptl;
	int spin;
	int ndeg;
	int mapnumber;
	int *istart;
	mapele *mapinfo;
	gsl_complex **ground;
} sectorinfo; 

#endif


