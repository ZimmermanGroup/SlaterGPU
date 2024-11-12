#ifndef BSRADII
#define BSRADII

 //only programmed first 36 elements for now

const float ang2bohr = 1.8897261;

/* const float bsradii_angstroms[36] = { 
   0.25, 1.40,
   1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 1.50,
   1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.80,
   2.20, 1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.05 
 }; */

const float bsradii[36] = { 
  0.472, 2.646,
  2.740, 1.984, 1.606, 1.323, 1.228, 1.134, 0.945, 2.835,
  3.402, 2.835, 2.362, 2.079, 1.890, 1.890, 1.890, 3.402,
  4.157, 3.402, 3.024, 2.646, 2.551, 2.646, 2.646, 2.646, 2.551, 2.551, 2.551, 2.551, 2.457, 2.362, 2.173, 2.173, 2.173, 1.984
 };

float get_radii(int Z1)
{ 
 //units are Bohr
  float val = 2.;

 //this is odd
  if (Z1<=1) val = 0.661;
  if (Z1<=36) val = bsradii[Z1-1];

 //CPMZ edits
 // if (Z1<=1) val = 0.8;
 // if (Z1==3) val = 2.0;
 //end edits

  return val;
}

//#pragma acc routine seq
float bsf(int a0, int a1, int a2)
{
  float val = 0.35*a0 + 0.85*a1 + a2;
  return val;
}

//#pragma acc routine seq
float get_radii_2(int Z)
{
  if (Z<=0) Z = 1;

  if (Z>36) { printf(" ERROR in Bragg-Slater radii \n"); exit(1); }

  int zn[4] = { 2, 10, 18, 36 };
  float zv[4] = { 1., 2., 3., 3.7 };
  float val0;
  for (int i=0;i<4;i++)
  if (Z <= zn[i])
  {
    val0 = zv[i];
    break;
  }

  float val1 = 0.;
  if (Z<=1) val1 = 0.;
  else if (Z==2) val1 = 0.3;
  else if (Z<=10) val1 = bsf(Z-3,2,0);
  else if (Z<=18) val1 = bsf(Z-11,8,2);
  else if (Z<=20) val1 = bsf(Z-19,8,10);
  else if (Z<=30) val1 = bsf(1,Z-12,10);
  else if (Z<=36) val1 = bsf(Z-29,18,10);

  float r1 = val0*val0/(Z-val1);

  if (Z==3)
  {
    //printf("  brsl, Z==3. val0: %8.5f val1: %8.5f r1: %8.5f \n",val0,val1,r1);
    r1 = 2.0;
  }

  return r1;
}


#endif
