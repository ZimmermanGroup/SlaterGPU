#include "gpu_util.h"
#include "integrals.h"

void copy_to_all_gpu(int ngpu, int s1, double* A, int include_first)
{
  if (ngpu<1) return;
  if (include_first==0)
  {
    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update self(A[0:s1])
  }

  int start = 1;
  if (include_first>0) start = 0;

 //#pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=start;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);
    #pragma acc update device(A[0:s1])
  }

  acc_set_device_num(0,acc_device_nvidia);
  return;
}

void copy_to_all_gpu(int ngpu, int s1, float* A, int include_first)
{
  if (ngpu<1) return;
  if (include_first==0)
  {
    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update self(A[0:s1])
  }

  int start = 1;
  if (include_first>0) start = 0;

 //#pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=start;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);
    #pragma acc update device(A[0:s1])
  }

  acc_set_device_num(0,acc_device_nvidia);
  return;
}

void copy_to_all_gpu(int ngpu, int s1, int s2, double** A, int include_first)
{
  if (ngpu<1) return;
  if (include_first==0)
  {
    acc_set_device_num(0,acc_device_nvidia);
    #pragma acc update self(A[0:s1][0:s2])
  }

  int start = 1;
  if (include_first>0) start = 0;

 //#pragma omp parallel for schedule(static,1) num_threads(nomp)
  for (int n=start;n<ngpu;n++)
  {
    acc_set_device_num(n,acc_device_nvidia);
    #pragma acc update device(A[0:s1][0:s2])
  }

  acc_set_device_num(0,acc_device_nvidia);
  return;
}