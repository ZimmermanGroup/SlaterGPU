#ifndef REDUCEH
#define REDUCEH
#include "fp_def.h"

#include "fp_def.h"
//reduce over single precision
void reduce_2c1(int s1, int s2, int gs, FP1** val1, FP1** val3, int iN, int N, FP1* An);
void reduce_2c2(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, int iN, int N, FP1* An);
void reduce_2c3(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int iN, int N, FP1* An);
void reduce_3c1(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1* valt1, int N, int Naux, int imaxN, int imaxNa, FP1* C);
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, int N, int Naux, int imaxN, int imaxNa, FP1* C);
void reduce_3c2(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1* valt1, FP1* valt2, int N, int Naux, int imaxN, int imaxNa, FP1* C);
void reduce_3c2b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int N, int Naux, int imaxN, int imaxNa, FP1* C);
void reduce_3c3(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, FP1* valt1, FP1* valt2, FP1* valt3, int N, int Naux, int imaxN, int imaxNa, FP1* C);
void reduce_3c3b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, int N, int Naux, int imaxN, int imaxNa, FP1* C);


//gradients

void reduce_2c1d(int m1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val3, FP1** val1x, FP1** val3x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad);
void reduce_2c1d(int m1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val3, FP1** val1x, FP1** val3x, int iN, int N, int natoms, FP2 scalar, FP1* xyz_g);

void reduce_2c2d(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad);
void reduce_2c2d(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP1* xyz_grad);

void reduce_2c1ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad);
void reduce_2c1ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP1* xyz_g);

void reduce_2c2ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad);
void reduce_2c2ds(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, int iN, int N, int natoms, FP2 scalar, FP1* xyz_grad);

void reduce_2c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, int iN, int N, int natoms, FP2 scalar, FP2* xyz_grad);
void reduce_2c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x,int iN, int N, int natoms, FP2 scalar, FP1* xyz_grad);

void reduce_2c3de(int m1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* Pao, FP1** val1, FP1** val2, FP1** val3, FP1** val1x, FP1** val2x, FP1** val3x, int N, int imaxN, int natoms, FP2 scalar, FP2* xyz_grad);
void reduce_2c3de(int m1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* Pao, FP1** val1, FP1** val2, FP1** val3, FP1** val1x, FP1** val2x, FP1** val3x, int N, int imaxN, FP1 scalar, FP1* xyz_grad);

void reduce_3c2d112(int m1, int n1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, int N, int Naux, int imaxN, int imaxNa, int natoms, FP2 scalar, FP2* xyz_grad);
void reduce_3c2d112(int m1, int n1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, int N, int Naux, int imaxN, int imaxNa, FP1 scalar, FP1* xyz_grad);

void reduce_3c2d122(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, int N, int Naux, int imaxN, int imaxNa, int natoms, FP2 scalar, FP2* xyz_grad);
void reduce_3c2d122(int m1, int n1, int s1, int s2, int s3, int s4, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, int N, int Naux, int imaxN, int imaxNa, FP1 scalar, FP1* xyz_grad);

void reduce_3c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, FP1** val7x, FP1** val8x, FP1** val9x, int N, int Naux, int imaxN, int imaxNa, int natoms, FP2 scalar, FP2* xyz_grad);
void reduce_3c3d(int m1, int n1, int p1, int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP2* norms1, FP2* norms2, FP2* dC, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** valx, FP1** val8, FP1** val9, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4x, FP1** val5x, FP1** val6x, FP1** val7x, FP1** val8x, FP1** val9x, int N, int Naux, int imaxN, int imaxNa, FP1 scalar, FP1* xyz_grad);


//reduce over FP2 precision
void reduce_2c1(int s1, int s2, int gs, FP1** val1, FP1** val3, int iN, int N, FP2* An);
void reduce_2c2(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, int iN, int N, FP2* An);
void reduce_2c3(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int iN, int N, FP2* An);
void reduce_3c1(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, FP1* valt1, int N, int Naux, int imaxN, int imaxNa, FP2* C);
void reduce_3c1b(int s1, int s2, int s3, int s4, int gs, FP1** val1, FP1** val2, FP1** val3, int N, int Naux, int imaxN, int imaxNa, FP2* C);
void reduce_3c2(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1* valt1, FP1* valt2, int N, int Naux, int imaxN, int imaxNa, FP2* C);
void reduce_3c2b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int N, int Naux, int imaxN, int imaxNa, FP2* C);
void reduce_3c3(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, FP1* valt1, FP1* valt2, FP1* valt3, int N, int Naux, int imaxN, int imaxNa, FP2* C);
void reduce_3c3b(int s1, int s2, int s3, int s4, int s5, int s6, int gs, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, FP1** val7, FP1** val8, FP1** val9, int N, int Naux, int imaxN, int imaxNa, FP2* C);



void reduce_2c2v(int p1, int s1, int s2, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, int iN, int N, int nc, FP2 scalar, FP2* V);
void reduce_2c2v(int p1, int s1, int s2, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, int iN, int N, int nc, FP2 scalar, FP1* V);
void reduce_2c2vd(int p1, int s1, int s2, int gs, FP2* norms, FP2* dpq, FP1** val1x, FP1** val2x, FP1** val3, FP1** val4, int iN, int N, int nc, FP2 scalar, FP2* dV);
void reduce_2c2vd(int p1, int s1, int s2, int gs, FP2* norms, FP2* dpq, FP1** val1x, FP1** val2x, FP1** val3, FP1** val4, int iN, int N, int nc, FP2 scalar, FP1* dV);


void reduce_2c3v(int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int iN, int N, int nc, FP2 scalar, FP2* V);
void reduce_2c3v(int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1, FP1** val2, FP1** val3, FP1** val4, FP1** val5, FP1** val6, int iN, int N, int nc, FP2 scalar, FP1* V);
void reduce_2c3vd(int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4, FP1** val5, FP1** val6, int iN, int N, int nc, FP2 scalar, FP2* dV);
void reduce_2c3vd(int p1, int s1, int s2, int s3, int s4, int gs, FP2* norms, FP2* dpq, FP1** val1x, FP1** val2x, FP1** val3x, FP1** val4, FP1** val5, FP1** val6, int iN, int N, int nc, FP2 scalar, FP1* dV);


#endif
