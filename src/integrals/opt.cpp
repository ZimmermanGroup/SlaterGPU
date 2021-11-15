#include "opt.h"
#include "fp_def.h"

FP2 cg_step(FP1 scalar, FP1 maxstep, FP2& grmsp, int natoms, FP2* step, FP2* grad, FP1* xyz)
{
  int N3 = 3*natoms;
  FP2 grms = 0.;
  for (int i=0;i<N3;i++)
    grms += grad[i]*grad[i];
  grms = sqrt(grms/N3);

  FP2 beta = 0.;
  FP2 max_beta = 0.8;
  if (grmsp>0.)
    beta = grms*grms/grmsp/grmsp;
  if (beta>max_beta) beta = max_beta;

  FP2 mag = 0.;
  for (int i=0;i<N3;i++)
  {
    FP2 step1 = -scalar*grad[i] + beta*step[i];
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
      FP2 step1 = step[i];
      mag += step1*step1;
    }
    mag = sqrt(mag);
  }

  for (int i=0;i<N3;i++)
    xyz[i] += step[i];

  grmsp = grms;

  return mag;
}


FP2 sd_step(FP1 scalar, FP1 maxstep, int natoms, FP2* grad, FP1* xyz)
{
  int N3 = 3*natoms;
  FP2 mag = 0.;
  FP2* step = new FP2[N3];
  for (int i=0;i<N3;i++)
  {
    FP2 step1 = scalar*grad[i];
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
      FP2 step1 = scalar*step[i];
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
