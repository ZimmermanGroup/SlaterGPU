#include "scf_util.h"


int count_zero_mo(int N, double* jCA)
{
  int nf = 0;
  for (int i=0;i<N;i++)
  {
    bool iszero = 1;
    for (int j=0;j<N;j++)
    if (fabs(jCA[j*N+i])>1.e-4)
    {
      iszero = 0;
      break;
    }
    if (iszero)
      nf++;
  }
  return nf;
}

int core_count(int natoms, int* atno)
{
  int Nc = 0;
  for (int m=0;m<natoms;m++)
  {
    if (atno[m]>36)
    {
      printf(" WARNING: is the grid updated for Z>36? \n");
      Nc += 19;
    }
    else if (atno[m]>18)
      Nc += 9;
    else if (atno[m]>10)
      Nc += 5;
    else if (atno[m]>2)
      Nc++;
  }

  return Nc;
}

int electron_count(int charge, int natoms, int* atno, int& Nc, int& Na, int& Nb)
{
  int No = 0;
  for (int m=0;m<natoms;m++)
    No += atno[m];
  No -= charge;
  if (No%2==0)
    Na = Nb = No/2;
  else
  {
    Nb = (No-1)/2;
    Na = Nb+1;
  }
  No = min(Na,Nb);

  Nc = core_count(natoms,atno);
  if (Nc>No) Nc = 0;

  return No;
}

void mo_coeff_to_pao_occ_complex(double* occ, int N, double* jCA, double* Pao, double* Pai)
{
  int N2 = N*N;
  //int M2 = 2*N2;
  double* jCAi = &jCA[N2];

  for (int m=0;m<N2;m++) Pao[m] = Pai[m] = 0.;

  for (int i=0;i<N;i++)
  {
    double d1 = occ[i];
    if (fabs(d1)>0.)
    for (int m=0;m<N;m++)
    {
      for (int n=0;n<N;n++)
        Pao[m*N+n] += d1*jCA[m*N+i]*jCA[n*N+i];
      for (int n=0;n<N;n++)
        Pao[m*N+n] += d1*jCAi[m*N+i]*jCAi[n*N+i];
    }
  }

  for (int i=0;i<N;i++)
  {
    double d1 = occ[i];
    if (fabs(d1)>0.)
    for (int m=0;m<N;m++)
    {
      for (int n=0;n<N;n++)
        Pai[m*N+n] += d1*jCA[m*N+i]*jCAi[n*N+i];
      for (int n=0;n<N;n++)
        Pai[m*N+n] -= d1*jCAi[m*N+i]*jCA[n*N+i];
    }
  }

  return;
}

void mo_coeff_to_pao_complex(int No, int N, double* jCA, double* Pao, double* Pai)
{
  int N2 = N*N;
  //int M2 = 2*N2;
  double* jCAi = &jCA[N2];

// mat rP = rC*strans(rC)+iC*strans(iC);
// mat iP = iC*strans(rC)-rC*strans(iC);

  for (int m=0;m<N2;m++) Pao[m] = Pai[m] = 0.;

  for (int i=0;i<No;i++)
  {
    for (int m=0;m<N;m++)
    {
      for (int n=0;n<N;n++)
        Pao[m*N+n] += jCA[m*N+i]*jCA[n*N+i];
      for (int n=0;n<N;n++)
        Pao[m*N+n] += jCAi[m*N+i]*jCAi[n*N+i];
    }
  }
  for (int m=0;m<N2;m++)
    Pao[m] *= 2.;

  for (int i=0;i<No;i++)
  {
    for (int m=0;m<N;m++)
    {
      for (int n=0;n<N;n++)
        Pai[m*N+n] += jCA[m*N+i]*jCAi[n*N+i];
      for (int n=0;n<N;n++)
        Pai[m*N+n] -= jCAi[m*N+i]*jCA[n*N+i];
    }
  }
  for (int m=0;m<N2;m++)
    Pai[m] *= 2.;

  return;
}

void ao_to_mo_complex(int No, int N, double* jAO, double* jMO, double* jCA, double* tmp)
{
  int N2 = N*N;
  int M2 = 2*N2;
  #pragma acc update self(jAO[0:N2],jCA[0:M2])

  int offset = N2;
  double* jCAi = &jCA[offset];
  double* jMOi = &jMO[offset];
  double* tmpi = &tmp[offset];

  for (int i=0;i<M2;i++) tmp[i] = 0.;
  for (int i=0;i<M2;i++) jMO[i] = 0.;

 //real only part
  //mat_times_mat_cpu(tmp,jAO,jCA,N);
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    tmp[i*N+j] += jAO[i*N+k]*jCA[k*N+j];

 //imaginary only part
  //mat_times_mat_cpu(tmpi,jAO,jCAi,N);
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    tmpi[i*N+j] += jAO[i*N+k]*jCAi[k*N+j];

 //real-real
  for (int i=0;i<No;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    jMO[i*N+j] += jCA[k*N+i]*tmp[k*N+j];

 //real-imag
  for (int i=0;i<No;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    jMOi[i*N+j] += jCA[k*N+i]*tmpi[k*N+j];

 //imag-imag (+ due to conjugate*i)
  for (int i=0;i<No;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    jMO[i*N+j] += jCAi[k*N+i]*tmpi[k*N+j];

 //imag-real (- due to conjugate)
  for (int i=0;i<No;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    jMOi[i*N+j] -= jCAi[k*N+i]*tmp[k*N+j];

  #pragma acc update device(jMO[0:M2])

  return;
}

void ao_to_mo_occ(int No, int N, double* jAO, double* jMO, double* jCA, double* tmp)
{
  int N2 = N*N;
  #pragma acc update self(jAO[0:N2],jCA[0:N2])

  for (int i=0;i<N2;i++) tmp[i] = 0.;
  for (int i=0;i<N2;i++) jMO[i] = 0.;

  //mat_times_mat(tmp,jAO,jCA,N);
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    tmp[i*N+j] += jAO[i*N+k]*jCA[k*N+j];

  for (int i=0;i<No;i++)
  for (int j=0;j<N;j++)
  for (int k=0;k<N;k++)
    jMO[i*N+j] += jCA[k*N+i]*tmp[k*N+j];

  #pragma acc update device(jMO[0:N2])

  return;
}

void ao_to_mo(int N, double* jAO, double* jMO, double* jCA, double* tmp)
{
 #if !USE_ACC
  return ao_to_mo_cpu(N,jAO,jMO,jCA,tmp);
 #endif
  mat_times_mat(tmp,jAO,jCA,N);
  mat_times_mat_at(jMO,jCA,tmp,N);

  return;
}

void ao_to_mo_cpu(int N, double* jAO, double* jMO, double* jCA, double* tmp)
{
  mat_times_mat_cpu(tmp,jAO,jCA,N);
  mat_times_mat_at_cpu(jMO,jCA,tmp,N);

  return;
}

void compute_CSC(int N, double* jCA1, double* S, double* jCA2, double* O)
{
  double* tmp = new double[N*N];

  mat_times_mat_at_cpu(tmp,jCA1,S,N);
  mat_times_mat_cpu(O,tmp,jCA2,N);

  delete [] tmp;
}

