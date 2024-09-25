#ifndef _GAMMA_H_
#define _GAMMA_H_

#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>

// #include "accel.h"
// Single precision values
#define MAXLOGF 88.72283905206835
#define MACHEPF 5.9604644775390625E-8
#define BIGF 16777216.0

// Double precision values
#define MAXLOG 7.08396418532264106224E2
#define MACHEP 1.11022302462515654042E-16
#define BIG 4.503599627370496e15


//Spherical Bessels

//#pragma acc routine seq
size_t mchn(int n, int k);
#pragma acc routine seq
size_t mchn_v2(int l, int m);

#pragma acc routine seq
double spherical_bessel_1(int l, double gr);

#pragma acc routine seq
double spherical_bessel_3(int l, double gr);


//Gamma ftns
#pragma acc routine seq
float igamf(float aa, float xx);
#pragma acc routine seq
float igamcf(float aa, float xx);

#pragma acc routine seq
double igam(double aa, double xx);
#pragma acc routine seq
double igamc(double aa, double xx);
#pragma acc routine seq
double igamn(int a, double x);

#pragma acc routine seq
float vinr_gamf(int n, int l, float r, float z);
#pragma acc routine seq
double vinr_gam(int n, int l, double r, double z);

#pragma acc routine seq
float dvinr_gamf(int n, int l, float r, float z);
#pragma acc routine seq
double dvinr_gam(int n, int l, double r, double z);

#endif // _GAMMA_H_
