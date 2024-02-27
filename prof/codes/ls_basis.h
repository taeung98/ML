#include "lattice.h"
#include "matrix.h"
void build_bitsubbasis(int ntotal, typebasis *subset);
void build_bitbasis(int ntotal, typebasis **bitbasis, sectorinfo *ginfo, int Nonzero_SOC);
int print_bitbasis(int ntotal, typebasis **basis);
int print_bitsubbasis(int ntotal, typebasis *subset);
int print_bit(int ntotal, typebasis number);
int count_bit(int ntotal, typebasis bit);
int find_refer(typebasis **basis, int refer, int block, typebasis *bin);
int permu(typebasis *bin, int site, int spin);
int permu_test(typebasis *bin, int site, int spin);
int find_gindex(sectorinfo *ginfo, int refer, int block, int Nblocks);
int find_index(int ntotal, int ptlgr, int netspingr, int Blocks);
int print_ginfo(int ntotal, sectorinfo *ginfo, int Blocks);
int print_bitbasis_coeff(int ntotal, typebasis **basis, double *coeff);

extern int nctot;

