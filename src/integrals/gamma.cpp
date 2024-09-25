#include "gamma.h"
#include <cmath>
// #include <accelmath.h>

//no range of a that will be used, can precompute?

using namespace std;

#pragma acc routine seq
double spherical_bessel_1(int l, double gr)
{
 //does not include i^(l+2) term
  double val = 0.;
  for (int m=0;m<=l;m++)
  {
   //inlined mchn_v2
    int f1 = 1;
    for (int i=2;i<l+m;i++)
      f1 *= i;
    for (int i=1;i<m;i++)
      f1 /= i;
    for (int i=1;i<l-m;i++)
      f1 /= i;

    double o2gr = pow(2.*gr,-(m+1));
    double f2 = pow(-1.,l-m+1);
    double egr = exp(-gr);
    double egrp = exp(gr);

    val += f1*o2gr*(egr+f2*egrp);
  }
  return val;
}

#pragma acc routine seq
double spherical_bessel_3(int l, double gr)
{
 //does not include i^(l+2) term
  double val = 0.;
  for (int m=0;m<=l;m++)
  {
    int f1 = 1;
    for (int i=2;i<l+m;i++)
      f1 *= i;
    for (int i=1;i<m;i++)
      f1 /= i;
    for (int i=1;i<l-m;i++)
      f1 /= i;

    double o2gr = pow(2.*gr,-(m+1));
    double egr = exp(-gr);

    val += f1*o2gr*egr;
  }
  return 2.*val;
}

#pragma acc routine seq
float igamf(float aa, float xx) {
  float a = aa;
  float x = xx;

  // x,a <= 0 should never happen for our code
  
  // x > 1 and x > a should be dealt with outside
  // this function. 
  if (x > 1.f && x > a) {
    return 1.f - igamcf(a,x);
  }
  float ax = a * logf(x) - x - lgammaf(a);
  if (ax < -MAXLOGF) {
    return 0.;
  }
  ax = expf(ax);

  float r = a;
  float c = 1.f;
  float ans = 1.f;

  while (c/ans > MACHEPF) {
  // for (int i = 0; i < 15; i++) {
    r += 1.0f;
    c *= x/r;
    ans += c;
  }
  return ans * ax / a;
}

#pragma acc routine seq
float igamcf(float aa, float xx) {
  float a = aa;
  float x = xx;
  float t, yc, pk, qk, r;
  t = 1.;
  
  // x,a <= 0 should never happen for our code
  
  // x < 1 and x < a should be dealt with outside
  // this function. 
  if (x < 1. || x < a) {
    return 1. - igamf(a,x);
  }

  float ax = a * logf(x) - x - lgammaf(a);
  if (ax < -MAXLOGF) {
    return 0.f;
  }
  ax = expf(ax);

  float y = 1.0f - a;
  float z = x + y + 1.f;
  float c = 0.f;
  float pkm2 = 1.f;
  float qkm2 = x;
  float pkm1 = x + 1.f;
  float qkm1 = z*x;
  float ans = pkm1 / qkm1;

  while (t > MACHEPF) {
  // for (int i = 0; i < 15; i++) {
    c += 1.f;
    y += 1.f;
    z += 2.f;
    yc = y*c;
    pk = pkm1 * z - pkm2 * yc;
    qk = qkm1 * z - qkm2 * yc;

    if (qk != 0.f) {
      r = pk / qk;
      t = fabsf((ans - r) / r);
      ans = r;
    }
    else {
      t = 1.f;
    }

    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    if (fabs(pk) > BIGF) {
      pkm2 *= MACHEPF;
      pkm1 *= MACHEPF;
      qkm2 *= MACHEPF;
      qkm1 *= MACHEPF;
    } // if fabs(pk)
  } // while t > MACHEPF
  return ans * ax;
}

#pragma acc routine seq
float vinr_gamf(int n, int l, float r, float z) {
  float ans = powf(r,(float)(-l-1));
  ans *= powf(z,(float)(-l-n-2));
  float rz = r*z;
  float rzl = powf(rz,(float)(2*l+1));
  int lnf = 1;
  for (int i = 1; i < n-l+1; i++) {
    lnf *= i;
  }
  int lnf2 = lnf;
  for (int i = n-l+1; i < n+l+2; i++) {
    lnf2 *= i;
  }

  float g1 = igamf(n-l+1,rz) * (float)lnf;
  float g2 = igamf(n+l+2,rz) * (float)lnf2;
  ans *= rzl * (lnf - g1) + g2;

  if (isnan(ans) || isinf(ans)) {
    ans = 0.;
  }

  return ans;
}

#pragma acc routine seq
float dvinr_gamf(int n, int l, float r, float z) {
  float ans = powf(r,(float)(-l-2));
  ans *= powf(z,(float)(-l-n-2));
  float rz = r*z;
  float l1 = (float)l;
  float rzl = l1 * powf(rz,(float)(2*l+1));
  int lnf = 1.;
  for (int i = 1; i < n-l+1; i++) {
    lnf *= i;
  }
  int lnf2 = lnf;
  for (int i = n-l+1; i < n+l+2; i++) {
    lnf2 *= i;
  }

  float g1 = igamf(n-l+1,rz) * (float)lnf;
  float g2 = igamf(n+l+2,rz) * (float)lnf2;
  ans *= rzl * (lnf - g1) - g2 * (l1 + 1.);

  if (isnan(ans) || isinf(ans)) {
    ans = 0.;
  }

  return ans;
}

#pragma acc routine seq
double igamn(int a, double x)
{
  int f1 = a-1;
  for (int i=2;i<a-1;i++)
    f1 *= i;
  double v1 = 1.;
  double d  = 1.;
  for (int i=1;i<a;i++)
  {
    v1 += pow(x,i)/d;
    d *= i;
  }
  return f1*(1.-v1*exp(-x));
}

#pragma acc routine seq
double igam(double aa, double xx) {
  double a = aa;
  double x = xx;

  // x,a <= 0 should never happen for our code

  // x > 1 and x > a should be dealt with outside
  // this function. 
  if (x > 1. && x > a) {
    return 1. - igamc(a,x);
  }

  double ax = a * log(x) - x - lgamma(a);
  if (ax < -MAXLOG) {
    return 0.;
  }
  ax = exp(ax);

  double r = a;
  double c = 1.;
  double ans = 1.;

  while (c/ans > MACHEP) {
  // for (int i = 0; i < 15; i++) {
    r += 1.0;
    c *= x/r;
    ans += c;
  }
  return ans * ax / a;
}

#pragma acc routine seq
double igamc(double aa, double xx) {
  double a = aa;
  double x = xx;
  double t, yc, pk, qk, r;
  t = 1.;
  // x,a <= 0 should never happen for our code
  
  // x < 1 and x < a should be dealt with outside
  // this function. 
  if (x < 1. || x < a) {
    return 1. - igam(a,x);
  }

  double ax = a * log(x) - x - lgamma(a);
  if (ax < -MAXLOG) {
    return 0.;
  }
  ax = exp(ax);

  double y = 1.0 - a;
  double z = x + y + 1.;
  double c = 0.;
  double pkm2 = 1.;
  double qkm2 = x;
  double pkm1 = x + 1;
  double qkm1 = z*x;
  double ans = pkm1 / qkm1;

  while (t > MACHEP) {
  // for (int i = 0; i < 15; i++) {
    c += 1.;
    y += 1.;
    z += 2.;
    yc = y*c;
    pk = pkm1 * z - pkm2 * yc;
    qk = qkm1 * z - qkm2 * yc;

    if (qk != 0) {
      r = pk / qk;
      t = fabs((ans - r) / r);
      ans = r;
    }
    else {
      t = 1.;
    }

    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    if (fabs(pk) > BIG) {
      pkm2 *= MACHEP;
      pkm1 *= MACHEP;
      qkm2 *= MACHEP;
      qkm1 *= MACHEP;
    } // if fabs(pk)
  } // while t > MACHEP
  return ans * ax;
}

double dvinr_gam(int n, int l, double r, double z) {
  double ans = pow(r,(double)(-l-2));
  ans *= pow(z,(double)(-l-n-2));
  double rz = r*z;
  double l1 = (double)l;
  double rzl = l1 * pow(rz,(double)(2*l+1));
  int lnf = 1;
  for (int i = 1; i < n-l+1; i++) {
    lnf *= i;
  }
  int lnf2 = lnf;
  for (int i = n-l+1; i < n+l+2; i++) {
    lnf2 *= i;
  }

  double g1 = igam(n-l+1,rz) * (double)lnf;
  double g2 = igam(n+l+2,rz) * (double)lnf2;
  ans *= rzl * (lnf - g1) - g2 * (l1 + 1.);

  if (std::isnan(ans) || std::isinf(ans)) {
    ans = 0.;
  }

  return ans;
}

#pragma acc routine seq
double vinr_gam(int n, int l, double r, double z) {
  double ans = pow(r,(double)(-l-1));
  ans *= pow(z,(double)(-l-n-2));
  double rz = r*z;
  double rzl = pow(rz,(double)(2*l+1));
  int lnf = 1;
  for (int i = 1; i < n-l+1; i++) {
    lnf *= i;
  }
  int lnf2 = lnf;
  for (int i = n-l+1; i < n+l+2; i++) {
    lnf2 *= i;
  }

  double g1 = igam(n-l+1,rz) * (double)lnf;
  double g2 = igam(n+l+2,rz) * (double)lnf2;
  ans *= rzl * (lnf - g1) + g2;
  
  if (std::isnan(ans) || std::isinf(ans)) {
    ans = 0.;
  }

  return ans;
}

