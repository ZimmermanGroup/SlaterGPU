#ifndef _GRID_UTIL_H_
#define _GRID_UTIL_H_


void get_angular_grid(int size_ang, double* ang_g, double* ang_w);
int electron_count(int charge, int natoms, int* atno, int& Nc, int& Na, int& Nb);

#endif // _GRID_UTIL_H_
