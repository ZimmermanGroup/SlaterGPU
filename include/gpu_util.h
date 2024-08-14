#ifndef GPU_UTIL_H
#define GPU_UTIL_H

void copy_to_all_gpu(int ngpu, int s1, double* A, int include_first);
void copy_to_all_gpu(int ngpu, int s1, float* A, int include_first);
void copy_to_all_gpu(int ngpu, int s1, int s2, double** A, int include_first);

#endif

