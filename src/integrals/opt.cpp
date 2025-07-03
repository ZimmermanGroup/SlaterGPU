#include "opt.h"

double cg_step(float scalar, float maxstep, double& grmsp, int natoms, double* step, double* grad, double* xyz)
{
  int N3 = 3*natoms;
  double grms = 0.;
  for (int i=0;i<N3;i++)
    grms += grad[i]*grad[i];
  grms = sqrt(grms/N3);

  double beta = 0.;
  double max_beta = 0.8;
  if (grmsp>0.)
    beta = grms*grms/grmsp/grmsp;
  if (beta>max_beta) beta = max_beta;

  double mag = 0.;
  for (int i=0;i<N3;i++)
  {
    double step1 = -scalar*grad[i] + beta*step[i];
    step[i] = step1;
    mag += step1*step1;
  }
  mag = sqrt(mag);

  if (mag>maxstep)
  {
    scalar = maxstep/mag;
    mag = 0.;
    for (int i=0;i<N3;i++)
    {
      step[i] *= scalar;
      double step1 = step[i];
      mag += step1*step1;
    }
    mag = sqrt(mag);
  }

  for (int i=0;i<N3;i++)
    xyz[i] += step[i];

  grmsp = grms;

  return mag;
}


double sd_step(float scalar, float maxstep, int natoms, double* grad, double* xyz)
{
  int N3 = 3*natoms;
  double mag = 0.;
  double* step = new double[N3];
  for (int i=0;i<N3;i++)
  {
    double step1 = scalar*grad[i];
    step[i] = step1;
    mag += step1*step1;
  }
  mag = sqrt(mag);

  if (mag>maxstep)
  {
    scalar = maxstep/mag;
    mag = 0.;
    for (int i=0;i<N3;i++)
    {
      double step1 = scalar*step[i];
      step[i] = step1;
      mag += step1*step1;
    }
    mag = sqrt(mag);
  }

  for (int i=0;i<N3;i++)
    xyz[i] -= step[i];

  delete [] step;

  return mag;
}

double sd_step(float scalar, float maxstep, int natoms, double* grad, float* xyz)
{
  int N3 = 3*natoms;
  double mag = 0.; 
  double* step = new double[N3];
  for (int i=0;i<N3;i++)
  {
    double step1 = scalar*grad[i];
    step[i] = step1;
    mag += step1*step1;
  }
  mag = sqrt(mag);

  if (mag>maxstep)
  {
    scalar = maxstep/mag;
    mag = 0.; 
    for (int i=0;i<N3;i++)
    {   
      double step1 = scalar*step[i];
      step[i] = step1;
      mag += step1*step1;
    }   
    mag = sqrt(mag);
  }

  for (int i=0;i<N3;i++)
    xyz[i] -= step[i];

  delete [] step;

  return mag;
}
