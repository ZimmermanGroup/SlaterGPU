#ifndef REDUCEH
#define REDUCEH


//reduce over single precision
void reduce_2c1(int s1, int s2, int gs, float** val1, float** val3, int iN, int N, float* An);
void reduce_2c2(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float** val4, int iN, int N, float* An);
void reduce_2c3(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int iN, int N, float* An);
void reduce_3c1(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float* valt1, int N, int Naux, int imaxN, int imaxNa, float* C);
void reduce_3c1(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float* valt1, int N, int Naux, int imaxN, int imaxNa, float* C, int tid);
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, int N, int Naux, int imaxN, int imaxNa, float* C);
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, int N, int Naux, int imaxN, int imaxNa, float* C, int tid);
void reduce_3c2(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float* valt1, float* valt2, int N, int Naux, int imaxN, int imaxNa, float* C);
void reduce_3c2b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int N, int Naux, int imaxN, int imaxNa, float* C);
void reduce_3c3(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, float* valt1, float* valt2, float* valt3, int N, int Naux, int imaxN, int imaxNa, float* C);
void reduce_3c3b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, int N, int Naux, int imaxN, int imaxNa, float* C);

//6d integrals
void reduce_4c_1d(int s1, int s2, int s3, int s4, int gs, int gsp, int M, float* grid1, float* grid2, float** val1a, float** val2a, float** val3a, float** val4a, float* wt1, float* wt2, float* gt);
void reduce_4c_1c(int s1, int s2, int s3, int s4, int gs, int gsp, int M, float* grid1, float* grid2, float** val1a, float** val2a, float** val3a, float** val4a, float* wt1, float* wt2, float* gt);
void reduce_4c_1b(int s1, int s2, int s3, int s4, int gs, int gsp, int M, float* grid1, float* grid2, float** val1a, float** val2a, float** val3a, float** val4a, float* wt1, float* wt2, float* gt);
void reduce_4c_1(int s1, int s2, int s3, int s4, int gs, int gsp, int M, float* grid1, float* grid2, float** val1a, float** val2a, float** val3a, float** val4a, float* wt1, float* wt2, float* gt);

//gradients

void reduce_2c1d(int m1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val3, float** val1x, float** val3x, int iN, int N, int natoms, double scalar, double* xyz_grad);
void reduce_2c1d(int m1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val3, float** val1x, float** val3x, int iN, int N, int natoms, double scalar, float* xyz_g);

void reduce_2c2d(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val1x, float** val2x, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, double* xyz_grad);
void reduce_2c2d(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val1x, float** val2x, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, float* xyz_grad);

void reduce_2c1ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, double* xyz_grad);
void reduce_2c1ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, float* xyz_g);

void reduce_2c2ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val1x, float** val2x, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, double* xyz_grad);
void reduce_2c2ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val1x, float** val2x, float** val3x, float** val4x, int iN, int N, int natoms, double scalar, float* xyz_grad);

void reduce_2c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, int iN, int N, int natoms, double scalar, double* xyz_grad);
void reduce_2c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x,int iN, int N, int natoms, double scalar, float* xyz_grad);

void reduce_2c3de(int m1, int s1, int s2, int s3, int s4, int gs, double* norms, double* Pao, float** val1, float** val2, float** val3, float** val1x, float** val2x, float** val3x, int N, int imaxN, int natoms, double scalar, double* xyz_grad);
void reduce_2c3de(int m1, int s1, int s2, int s3, int s4, int gs, double* norms, double* Pao, float** val1, float** val2, float** val3, float** val1x, float** val2x, float** val3x, int N, int imaxN, float scalar, float* xyz_grad);

void reduce_3c2d112(int m1, int n1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, int N, int Naux, int imaxN, int imaxNa, int natoms, double scalar, double* xyz_grad);
void reduce_3c2d112(int m1, int n1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, int N, int Naux, int imaxN, int imaxNa, float scalar, float* xyz_grad);

void reduce_3c2d122(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, int N, int Naux, int imaxN, int imaxNa, int natoms, double scalar, double* xyz_grad);
void reduce_3c2d122(int m1, int n1, int s1, int s2, int s3, int s4, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, int N, int Naux, int imaxN, int imaxNa, float scalar, float* xyz_grad);

void reduce_3c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, float** val7x, float** val8x, float** val9x, int N, int Naux, int imaxN, int imaxNa, int natoms, double scalar, double* xyz_grad);
void reduce_3c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, double* norms1, double* norms2, double* dC, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** valx, float** val8, float** val9, float** val1x, float** val2x, float** val3x, float** val4x, float** val5x, float** val6x, float** val7x, float** val8x, float** val9x, int N, int Naux, int imaxN, int imaxNa, float scalar, float* xyz_grad);


//reduce over double precision
void reduce_2c1(int s1, int s2, int gs, float** val1, float** val3, int iN, int N, double* An);
void reduce_2c2(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float** val4, int iN, int N, double* An);
void reduce_2c3(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int iN, int N, double* An);
void reduce_3c1(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float* valt1, int N, int Naux, int imaxN, int imaxNa, double* C);
void reduce_3c1(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, float* valt1, int N, int Naux, int imaxN, int imaxNa, double* C, int tid);
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, int N, int Naux, int imaxN, int imaxNa, double* C);
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, float** val1, float** val2, float** val3, int N, int Naux, int imaxN, int imaxNa, double* C, int tid);
void reduce_3c2(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float* valt1, float* valt2, int N, int Naux, int imaxN, int imaxNa, double* C);
void reduce_3c2b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int N, int Naux, int imaxN, int imaxNa, double* C);
void reduce_3c3(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, float* valt1, float* valt2, float* valt3, int N, int Naux, int imaxN, int imaxNa, double* C);
void reduce_3c3b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, float** val7, float** val8, float** val9, int N, int Naux, int imaxN, int imaxNa, double* C);



void reduce_2c2v(int p1, int s1, int s2, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, int iN, int N, int nc, double scalar, double* V);
void reduce_2c2v(int p1, int s1, int s2, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, int iN, int N, int nc, double scalar, float* V);
void reduce_2c2vd(int p1, int s1, int s2, int gs, double* norms, double* dpq, float** val1x, float** val2x, float** val3, float** val4, int iN, int N, int nc, double scalar, double* dV);
void reduce_2c2vd(int p1, int s1, int s2, int gs, double* norms, double* dpq, float** val1x, float** val2x, float** val3, float** val4, int iN, int N, int nc, double scalar, float* dV);


void reduce_2c3v(int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int iN, int N, int nc, double scalar, double* V);
void reduce_2c3v(int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1, float** val2, float** val3, float** val4, float** val5, float** val6, int iN, int N, int nc, double scalar, float* V);
void reduce_2c3vd(int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1x, float** val2x, float** val3x, float** val4, float** val5, float** val6, int iN, int N, int nc, double scalar, double* dV);
void reduce_2c3vd(int p1, int s1, int s2, int s3, int s4, int gs, double* norms, double* dpq, float** val1x, float** val2x, float** val3x, float** val4, float** val5, float** val6, int iN, int N, int nc, double scalar, float* dV);

//reduce fully double precision
void reduce_2c1(int s1, int s2, int gs, double** val1, double** val3, int iN, int N, double* An);
void reduce_2c1(int s1, int s2, int s3, int s4, int gs, double** val1, double** val3, int iN, int N, double* An);
void reduce_2c1(int s1, int s2, int s3, int s4, int gs, double** val1, double** val3, int iN, int N, double* An, int tid);
void reduce_2c2(int s1, int s2, int s3, int s4, int gs, double** val1, double** val2, double** val3, double** val4, int iN, int N, double* An);
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, double** val1, double** val2, double** val3, int N, int Naux, int imaxN, int imaxNa, double* C);
void reduce_3c1b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, double** val1, double** val2, double** val3, int N, int Naux, int imaxN, int imaxNa, double* C);
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, double** val1, double** val2, double** val3, int N, int Naux, int imaxN, int imaxNa, double* C, int tid); 
void reduce_3c1b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, double** val1, double** val2, double** val3, int N, int Naux, int imaxN, int imaxNa, double* C, int tid);

#endif
