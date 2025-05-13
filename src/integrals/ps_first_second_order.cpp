#include "prosph.h"


//idp: integral over phi DOF
double get_idp_one_atom(bool do_coulomb, int l1, int m1, int l2, int m2)
{
  if (l1!=l2 && do_coulomb) return 0.;

  return get_idp(l1,m1,l2,m2);
}

double get_idp(int l1, int m1, int l2, int m2)
{
 //zero due to cylindrical symmetry
  if (m1!=m2)
    return 0.;

 //z-only terms: 2*pi
  if (m1==0)
    return 2.*PI;

 //this one doesn't make much sense
 //issue crept in due to alignment of xy axes
 // see phi_phase variable
  //if (m1==-3)
  //  return 2.*PI;

  return PI;
}

double get_idp_4b_m3(int m1, int m2, int m3, int m4)
{
 //CPMZ probably incomplete. will need to test
  if (m1<-3 || m1>3) { printf("\n ERROR: m1 out of range: %2i \n",m1); exit(-1); }
  if (m2<-3 || m2>3) { printf("\n ERROR: m2 out of range: %2i \n",m2); exit(-1); }
  if (m3<-3 || m3>3) { printf("\n ERROR: m3 out of range: %2i \n",m3); exit(-1); }

  int m10 = m1 + 3;
  int m20 = m2 + 3;
  int m30 = m3 + 3;
  int m40 = m4 + 3;

  double PI2 = PI/2.;
  double PI4 = PI2/2.;

  double phi_integrals[7][7][7][7] =
    {{{{(3.*PI)/4, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 0}, {0,
    0, PI2, 0, 0, 0, 0}, {0, 0, 0, PI, 0, 0, 0}, {0, 0, 0,
    0, PI2, 0, 0}, {0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0, 0, 0,
    0, PI4}}, {{0, PI2, 0, 0, 0, 0, 0}, {PI2, 0, PI4,
    0, 0, 0, 0}, {0, PI4, 0, 0, 0, 0, 0}, {0, 0, 0, 0, PI2, 0,
     0}, {0, 0, 0, PI2, 0, PI4, 0}, {0, 0, 0, 0, PI4, 0,
    0}, {0, 0, 0, 0, 0, 0, 0}}, {{0, 0, PI2, 0, 0, 0,
    0}, {0, PI4, 0, 0, 0, 0, 0}, {PI2, 0, -(PI4), 0, 0, 0,
     0}, {0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0, 0, PI4, 0, 0}, {0,
     0, 0, PI2, 0, -(PI4), 0}, {0, 0, 0, 0, 0, 0, 0}}, {{0, 0,
     0, PI, 0, 0, 0}, {0, 0, 0, 0, PI2, 0, 0}, {0, 0, 0, 0,
    0, PI2, 0}, {PI, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, 0, 0,
    0, 0}, {0, 0, PI2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}}, {{0,
    0, 0, 0, PI2, 0, 0}, {0, 0, 0, PI2, 0, PI4, 0}, {0, 0,
     0, 0, PI4, 0, 0}, {0, PI2, 0, 0, 0, 0, 0}, {PI2,
    0, PI4, 0, 0, 0, 0}, {0, PI4, 0, 0, 0, 0, 0}, {0, 0, 0, 0,
     0, 0, 0}}, {{0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0, 0, PI4, 0,
     0}, {0, 0, 0, PI2, 0, -(PI4), 0}, {0, 0, PI2, 0, 0,
    0, 0}, {0, PI4, 0, 0, 0, 0, 0}, {PI2, 0, -(PI4), 0, 0,
     0, 0}, {0, 0, 0, 0, 0, 0, 0}}, {{0, 0, 0, 0, 0, 0, PI4}, {0, 
    0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
    0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {PI4, 0, 0, 
    0, 0, 0, 0}}}, {{{0, PI2, 0, 0, 0, 0, 0}, {PI2, 0, PI4,
    0, 0, 0, 0}, {0, PI4, 0, 0, 0, 0, 0}, {0, 0, 0, 0, PI2,
     0, 0}, {0, 0, 0, PI2, 0, PI4, 0}, {0, 0, 0, 0, PI4, 
    0, 0}, {0, 0, 0, 0, 0, 0, 0}}, {{PI2, 0, PI4, 0, 0, 0, 
    0}, {0, (3.*PI)/4, 0, 0, 0, 0, 0}, {PI4, 0, PI2, 0, 0, 
    0, 0}, {0, 0, 0, PI, 0, 0, 0}, {0, 0, 0, 0, PI2,
    0, -(PI4)}, {0, 0, 0, 0, 0, PI4, 0}, {0, 0, 0,
    0, -(PI4), 0, PI2}}, {{0, PI4, 0, 0, 0, 0, 0}, {PI4, 
    0, PI2, 0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 0}, {0, 0, 0,
     0, PI2, 0, -(PI2)}, {0, 0, 0, PI2, 0, 0, 0}, {0, 0, 
    0, 0, 0, 0, PI4}, {0, 0, 0, -(PI2), 0, PI4, 0}}, {{0, 
    0, 0, 0, PI2, 0, 0}, {0, 0, 0, PI, 0, 0, 0}, {0, 0, 0, 
    0, PI2, 0, -(PI2)}, {0, PI, 0, 0, 0, 0, 0}, {PI2, 
    0, PI2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, -(PI2),
     0, 0, 0, 0}}, {{0, 0, 0, PI2, 0, PI4, 0}, {0, 0, 0, 
    0, PI2, 0, -(PI4)}, {0, 0, 0, PI2, 0, 0, 0}, {PI2,
     0, PI2, 0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 0}, {PI4, 
    0, 0, 0, 0, 0, 0}, {0, -(PI4), 0, 0, 0, 0, 0}}, {{0, 0, 0, 
    0, PI4, 0, 0}, {0, 0, 0, 0, 0, PI4, 0}, {0, 0, 0, 0, 0, 
    0, PI4}, {0, 0, 0, 0, 0, 0, 0}, {PI4, 0, 0, 0, 0, 0, 
    0}, {0, PI4, 0, 0, 0, 0, 0}, {0, 0, PI4, 0, 0, 0, 
    0}}, {{0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, -(PI4), 0, PI2}
    , {0, 0, 0, -(PI2), 0, PI4, 0}, {0, 0, -(PI2), 0, 0,
     0, 0}, {0, -(PI4), 0, 0, 0, 0, 0}, {0, 0, PI4, 0, 0, 0, 
    0}, {0, PI2, 0, 0, 0, 0, 0}}}, {{{0, 0, PI2, 0, 0, 0, 
    0}, {0, PI4, 0, 0, 0, 0, 0}, {PI2, 0, -(PI4), 0, 0, 0,
     0}, {0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0, 0, PI4, 0, 0}, {0,
     0, 0, PI2, 0, -(PI4), 0}, {0, 0, 0, 0, 0, 0, 
    0}}, {{0, PI4, 0, 0, 0, 0, 0}, {PI4, 0, PI2, 0, 0, 0, 
    0}, {0, PI2, 0, 0, 0, 0, 0}, {0, 0, 0, 0, PI2, 
    0, -(PI2)}, {0, 0, 0, PI2, 0, 0, 0}, {0, 0, 0, 0, 0, 
    0, PI4}, {0, 0, 0, -(PI2), 0, PI4, 0}}, {{PI2, 
    0, -(PI4), 0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 
    0}, {-(PI4), 0, (3.*PI)/4, 0, 0, 0, 0}, {0, 0, 0, PI, 
    0, -(PI2), 0}, {0, 0, 0, 0, PI4, 0, -(PI4)}, {0, 0, 
    0, -(PI2), 0, PI2, 0}, {0, 0, 0, 0, -(PI4), 0, PI2}
    }, {{0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0, 0, PI2, 
    0, -(PI2)}, {0, 0, 0, PI, 0, -(PI2), 0}, {0, 0, PI, 
    0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 0}, {PI2, 0, -(PI2),
     0, 0, 0, 0}, {0, -(PI2), 0, 0, 0, 0, 0}}, {{0, 0, 0, 
    0, PI4, 0, 0}, {0, 0, 0, PI2, 0, 0, 0}, {0, 0, 0, 
    0, PI4, 0, -(PI4)}, {0, PI2, 0, 0, 0, 0, 0}, {PI4,
     0, PI4, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 
    0, -(PI4), 0, 0, 0, 0}}, {{0, 0, 0, PI2, 0, -(PI4), 
    0}, {0, 0, 0, 0, 0, 0, PI4}, {0, 0, 0, -(PI2), 0, PI2,
     0}, {PI2, 0, -(PI2), 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
    0}, {-(PI4), 0, PI2, 0, 0, 0, 0}, {0, PI4, 0, 0, 0, 0,
     0}}, {{0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, -(PI2), 0, PI4, 
    0}, {0, 0, 0, 0, -(PI4), 0, PI2}, {0, -(PI2), 0, 0, 0,
     0, 0}, {0, 0, -(PI4), 0, 0, 0, 0}, {0, PI4, 0, 0, 0, 0, 
    0}, {0, 0, PI2, 0, 0, 0, 0}}}, {{{0, 0, 0, PI, 0, 0, 
    0}, {0, 0, 0, 0, PI2, 0, 0}, {0, 0, 0, 0, 0, PI2, 
    0}, {PI, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 0}, {0, 
    0, PI2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}}, {{0, 0, 0, 
    0, PI2, 0, 0}, {0, 0, 0, PI, 0, 0, 0}, {0, 0, 0, 0, PI2, 0,
    -(PI2)}, {0, PI, 0, 0, 0, 0, 0}, {PI2, 0, PI2, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0}, {0, 0, -(PI2), 0, 0, 0,
     0}}, {{0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0, 0, PI2, 
    0, -(PI2)}, {0, 0, 0, PI, 0, -(PI2), 0}, {0, 0, PI, 
    0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 0}, {PI2, 0, -(PI2),
     0, 0, 0, 0}, {0, -(PI2), 0, 0, 0, 0, 0}}, {{PI, 0, 0, 0, 
    0, 0, 0}, {0, PI, 0, 0, 0, 0, 0}, {0, 0, PI, 0, 0, 0, 
    0}, {0, 0, 0, 2.*PI, 0, 0, 0}, {0, 0, 0, 0, PI, 0, 0}, {0, 0,
     0, 0, 0, PI, 0}, {0, 0, 0, 0, 0, 0, PI}}, {{0, PI2, 0, 
    0, 0, 0, 0}, {PI2, 0, PI2, 0, 0, 0, 0}, {0, PI2, 0, 0,
     0, 0, 0}, {0, 0, 0, 0, PI, 0, 0}, {0, 0, 0, PI, 0, PI2,
     0}, {0, 0, 0, 0, PI2, 0, PI2}, {0, 0, 0, 0, 0, PI2, 
    0}}, {{0, 0, PI2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {PI2,
    0, -(PI2), 0, 0, 0, 0}, {0, 0, 0, 0, 0, PI, 0}, {0, 0, 
    0, 0, PI2, 0, PI2}, {0, 0, 0, PI, 0, 0, 0}, {0, 0, 0, 
    0, PI2, 0, 0}}, {{0, 0, 0, 0, 0, 0, 0}, {0, 0, -(PI2), 0, 
    0, 0, 0}, {0, -(PI2), 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 
    0, PI}, {0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0, 0, PI2, 0, 
    0}, {0, 0, 0, PI, 0, 0, 0}}}, {{{0, 0, 0, 0, PI2, 0, 
    0}, {0, 0, 0, PI2, 0, PI4, 0}, {0, 0, 0, 0, PI4, 0, 
    0}, {0, PI2, 0, 0, 0, 0, 0}, {PI2, 0, PI4, 0, 0, 0, 
    0}, {0, PI4, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}}, {{0, 0, 
    0, PI2, 0, PI4, 0}, {0, 0, 0, 0, PI2, 
    0, -(PI4)}, {0, 0, 0, PI2, 0, 0, 0}, {PI2, 0, PI2,
     0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 0}, {PI4, 0, 0, 0, 0, 
    0, 0}, {0, -(PI4), 0, 0, 0, 0, 0}}, {{0, 0, 0, 0, PI4, 0, 
    0}, {0, 0, 0, PI2, 0, 0, 0}, {0, 0, 0, 0, PI4, 
    0, -(PI4)}, {0, PI2, 0, 0, 0, 0, 0}, {PI4, 0, PI4,
     0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, -(PI4), 0, 0, 0, 
    0}}, {{0, PI2, 0, 0, 0, 0, 0}, {PI2, 0, PI2, 0, 0, 0, 
    0}, {0, PI2, 0, 0, 0, 0, 0}, {0, 0, 0, 0, PI, 0, 0}, {0, 0,
     0, PI, 0, PI2, 0}, {0, 0, 0, 0, PI2, 0, PI2}, {0, 
    0, 0, 0, 0, PI2, 0}}, {{PI2, 0, PI4, 0, 0, 0, 
    0}, {0, PI2, 0, 0, 0, 0, 0}, {PI4, 0, PI4, 0, 0, 0, 
    0}, {0, 0, 0, PI, 0, PI2, 0}, {0, 0, 0, 0, (3.*PI)/4, 
    0, PI4}, {0, 0, 0, PI2, 0, PI2, 0}, {0, 0, 0, 
    0, PI4, 0, PI2}}, {{0, PI4, 0, 0, 0, 0, 0}, {PI4, 
    0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, PI2, 
    0, PI2}, {0, 0, 0, PI2, 0, PI2, 0}, {0, 0, 0, 
    0, PI2, 0, PI4}, {0, 0, 0, PI2, 0, PI4, 0}}, {{0, 
    0, 0, 0, 0, 0, 0}, {0, -(PI4), 0, 0, 0, 0, 0}, {0, 
    0, -(PI4), 0, 0, 0, 0}, {0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0,
     0, PI4, 0, PI2}, {0, 0, 0, PI2, 0, PI4, 0}, {0, 
    0, 0, 0, PI2, 0, 0}}}, {{{0, 0, 0, 0, 0, PI2, 0}, {0, 0, 
    0, 0, PI4, 0, 0}, {0, 0, 0, PI2, 0, -(PI4), 0}, {0, 
    0, PI2, 0, 0, 0, 0}, {0, PI4, 0, 0, 0, 0, 0}, {PI2, 
    0, -(PI4), 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}}, {{0, 0, 0, 
    0, PI4, 0, 0}, {0, 0, 0, 0, 0, PI4, 0}, {0, 0, 0, 0, 0, 
    0, PI4}, {0, 0, 0, 0, 0, 0, 0}, {PI4, 0, 0, 0, 0, 0, 
    0}, {0, PI4, 0, 0, 0, 0, 0}, {0, 0, PI4, 0, 0, 0, 
    0}}, {{0, 0, 0, PI2, 0, -(PI4), 0}, {0, 0, 0, 0, 0, 
    0, PI4}, {0, 0, 0, -(PI2), 0, PI2, 0}, {PI2, 
    0, -(PI2), 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {-(PI4), 
    0, PI2, 0, 0, 0, 0}, {0, PI4, 0, 0, 0, 0, 0}}, {{0, 
    0, PI2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {PI2, 
    0, -(PI2), 0, 0, 0, 0}, {0, 0, 0, 0, 0, PI, 0}, {0, 0, 0, 
    0, PI2, 0, PI2}, {0, 0, 0, PI, 0, 0, 0}, {0, 0, 0, 
    0, PI2, 0, 0}}, {{0, PI4, 0, 0, 0, 0, 0}, {PI4, 0, 0, 
    0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, PI2, 
    0, PI2}, {0, 0, 0, PI2, 0, PI2, 0}, {0, 0, 0, 
    0, PI2, 0, PI4}, {0, 0, 0, PI2, 0, PI4, 
    0}}, {{PI2, 0, -(PI4), 0, 0, 0, 0}, {0, PI4, 0, 0, 0, 
    0, 0}, {-(PI4), 0, PI2, 0, 0, 0, 0}, {0, 0, 0, PI, 0, 
    0, 0}, {0, 0, 0, 0, PI2, 0, PI4}, {0, 0, 0, 0, 0, (3.*PI)/4, 0}, 
   {0, 0, 0, 0, PI4, 0, PI2}}, {{0, 0, 0, 0, 
    0, 0, 0}, {0, 0, PI4, 0, 0, 0, 0}, {0, PI4, 0, 0, 0, 0, 
    0}, {0, 0, 0, 0, PI2, 0, 0}, {0, 0, 0, PI2, 0, PI4, 
    0}, {0, 0, 0, 0, PI4, 0, PI2}, {0, 0, 0, 0, 0, PI2, 
    0}}}, {{{0, 0, 0, 0, 0, 0, PI4}, {0, 0, 0, 0, 0, 0, 0}, {0, 0,
     0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0,
     0, 0, 0, 0, 0, 0}, {PI4, 0, 0, 0, 0, 0, 0}}, {{0, 0, 0, 0, 0,
     0, 0}, {0, 0, 0, 0, -(PI4), 0, PI2}, {0, 0, 
    0, -(PI2), 0, PI4, 0}, {0, 0, -(PI2), 0, 0, 0, 
    0}, {0, -(PI4), 0, 0, 0, 0, 0}, {0, 0, PI4, 0, 0, 0, 
    0}, {0, PI2, 0, 0, 0, 0, 0}}, {{0, 0, 0, 0, 0, 0, 0}, {0, 0, 
    0, -(PI2), 0, PI4, 0}, {0, 0, 0, 0, -(PI4), 0, PI2},
   {0, -(PI2), 0, 0, 0, 0, 0}, {0, 0, -(PI4), 0, 0, 0, 
    0}, {0, PI4, 0, 0, 0, 0, 0}, {0, 0, PI2, 0, 0, 0, 
    0}}, {{0, 0, 0, 0, 0, 0, 0}, {0, 0, -(PI2), 0, 0, 0, 
    0}, {0, -(PI2), 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, PI}, {0,
     0, 0, 0, 0, PI2, 0}, {0, 0, 0, 0, PI2, 0, 0}, {0, 0, 
    0, PI, 0, 0, 0}}, {{0, 0, 0, 0, 0, 0, 0}, {0, -(PI4), 0, 0,
     0, 0, 0}, {0, 0, -(PI4), 0, 0, 0, 0}, {0, 0, 0, 0, 0, PI2, 0}, 
    {0, 0, 0, 0, PI4, 0, PI2}, {0, 0, 0, PI2, 
    0, PI4, 0}, {0, 0, 0, 0, PI2, 0, 0}}, {{0, 0, 0, 0, 0, 0, 
    0}, {0, 0, PI4, 0, 0, 0, 0}, {0, PI4, 0, 0, 0, 0, 0}, {0, 
    0, 0, 0, PI2, 0, 0}, {0, 0, 0, PI2, 0, PI4, 0}, {0, 0,
     0, 0, PI4, 0, PI2}, {0, 0, 0, 0, 0, PI2, 
    0}}, {{PI4, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 
    0}, {0, 0, PI2, 0, 0, 0, 0}, {0, 0, 0, PI, 0, 0, 0}, {0, 0,
     0, 0, PI2, 0, 0}, {0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0, 0,
    0, 0, (3.*PI)/4}}}};

  {
    int x1 = -2; int x2 = -1; int x3 = 0; int x4 = 1; double v1 = 2.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -2; int x3 = 0; int x4 = 1; double v1 = 4./3.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -2; int x2 = -1; int x3 = 0; int x4 = 3; double v1 = 2.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -1; int x3 = 0; int x4 = 2; double v1 = 4.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -3; int x3 = -2; int x4 = -2; double v1 = 4./3.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -3; int x3 = -1; int x4 = -1; double v1 = 4.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -3; int x3 = 1; int x4 = 1; double v1 = 4./3.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -3; int x3 = 2; int x4 = 2; double v1 = 4.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -2; int x3 = 1; int x4 = 2; double v1 = 8./3.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -1; int x3 = 2; int x4 = 2; double v1 = 8.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -2; int x3 = -2; int x4 = -1; double v1 = 8./3.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -2; int x2 = -2; int x3 = 1; int x4 = 3; double v1 = 2.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -2; int x2 = -2; int x3 = 3; int x4 = 3; double v1 = 2.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -1; int x3 = 1; int x4 = 1; double v1 = 8./3.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }

  {
    int x1 = -3; int x2 = -1; int x3 = -1; int x4 = -1; double v1 = 8.;
    int x10 = x1 + 3; int x20 = x2 + 3; int x30 = x3 + 3; int x40 = x4 + 3;

    vector<int> xs;
    xs.push_back(x10); xs.push_back(x20); xs.push_back(x30); xs.push_back(x40);
    do { phi_integrals[xs[0]][xs[1]][xs[2]][xs[3]] *= v1; }
    while (next_permutation(xs.begin(), xs.end()));
  }


  return fabs(phi_integrals[m10][m20][m30][m40]);
}

double get_idp_m5(int m1, int m2, int m3)
{
  if (m1<-5 || m1>5) { printf("\n ERROR: m1 out of range: %2i \n",m1); exit(-1); }
  if (m2<-5 || m2>5) { printf("\n ERROR: m2 out of range: %2i \n",m2); exit(-1); }
  if (m3<-5 || m3>5) { printf("\n ERROR: m3 out of range: %2i \n",m3); exit(-1); }

  int m10 = m1 + 5;
  int m20 = m2 + 5;
  int m30 = m3 + 5;

  double PI2 = PI/2.;

 //table generated in Mathematica.
 // See "Phi Integrals Phase.nb"

  double phi_integrals [11][11][11] = {{{0, 0, 0, 0, 0, PI, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, PI2, 
   0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, PI2, 0, 0, 0}, {0, 0, 0, 0, 
   0, 0, 0, 0, PI2, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, PI2, 
   0}, {PI, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0,
    0, 0, 0, 0, 0}, {0, 0, PI2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 
   0, PI2, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, PI2, 0, 0, 0, 0,
    0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, {{0, 0, 0, 0, 0, 
   0, PI2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, PI, 0, 0, 0, 0, 0}, {0,
    0, 0, 0, 0, 0, PI2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, PI2,
   0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, PI2, 
   0, -(PI2)}, {0, PI, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {PI2, 
   0, PI2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, PI2, 0, 0, 0, 0,
    0, 0, 0}, {0, 0, 0, 0, PI2, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, -(PI2), 0, 0, 0, 0, 0, 
   0}}, {{0, 0, 0, 0, 0, 0, 0, PI2, 0, 0, 0}, {0, 0, 0, 0, 0, 
   0, PI2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, PI, 0, 0, 0, 0, 0}, {0,
    0, 0, 0, 0, 0, PI2, 0, 0, 0, -(PI2)}, {0, 0, 0, 0, 0, 0, 
   0, PI2, 0, -(PI2), 0}, {0, 0, PI, 0, 0, 0, 0, 0, 0, 0, 
   0}, {0, PI2, 0, PI2, 0, 0, 0, 0, 0, 0, 0}, {PI2, 0, 0, 
   0, PI2, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0}, {0, 0, 0, 0, -(PI2), 0, 0, 0, 0, 0, 0}, {0, 0, 
   0, -(PI2), 0, 0, 0, 0, 0, 0, 0}}, {{0, 0, 0, 0, 0, 0, 0, 
   0, PI2, 0, 0}, {0, 0, 0, 0, 0, 0, 0, PI2, 0, 0, 0}, {0, 0, 
   0, 0, 0, 0, PI2, 0, 0, 0, -(PI2)}, {0, 0, 0, 0, 0, PI, 
   0, 0, 0, -(PI2), 0}, {0, 0, 0, 0, 0, 0, PI2, 0, -(PI2),
    0, 0}, {0, 0, 0, PI, 0, 0, 0, 0, 0, 0, 0}, {0, 0, PI2, 
   0, PI2, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 0, 0, 0, 0,
    0}, {PI2, 0, 0, 0, -(PI2), 0, 0, 0, 0, 0, 0}, {0, 0, 
   0, -(PI2), 0, 0, 0, 0, 0, 0, 0}, {0, 0, -(PI2), 0, 0, 0, 0,
    0, 0, 0, 0}}, {{0, 0, 0, 0, 0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0, 
   0, 0, 0, 0, 0, PI2, 0, -(PI2)}, {0, 0, 0, 0, 0, 0, 
   0, PI2, 0, -(PI2), 0}, {0, 0, 0, 0, 0, 0, PI2, 
   0, -(PI2), 0, 0}, {0, 0, 0, 0, 0, PI, 0, -(PI2), 0, 0, 
   0}, {0, 0, 0, 0, PI, 0, 0, 0, 0, 0, 0}, {0, 0, 0, PI2, 0, 0,
    0, 0, 0, 0, 0}, {0, 0, PI2, 0, -(PI2), 0, 0, 0, 0, 0, 
   0}, {0, PI2, 0, -(PI2), 0, 0, 0, 0, 0, 0, 0}, {PI2, 
   0, -(PI2), 0, 0, 0, 0, 0, 0, 0, 0}, {0, -(PI2), 0, 0, 0, 0,
    0, 0, 0, 0, 0}}, {{PI, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0}, {0, PI, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, PI, 0, 0, 0, 
   0, 0, 0, 0, 0}, {0, 0, 0, PI, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 
   0, PI, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 2.*PI, 0, 0, 0, 0, 
   0}, {0, 0, 0, 0, 0, 0, PI, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
   0, PI, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, PI, 0, 0}, {0, 0, 
   0, 0, 0, 0, 0, 0, 0, PI, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, PI}}, {{0, PI2, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {PI2, 
   0, PI2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, PI2, 0, 0,
    0, 0, 0, 0, 0}, {0, 0, PI2, 0, PI2, 0, 0, 0, 0, 0, 0}, {0,
    0, 0, PI2, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, PI, 0, 
   0, 0, 0}, {0, 0, 0, 0, 0, PI, 0, PI2, 0, 0, 0}, {0, 0, 0, 0,
    0, 0, PI2, 0, PI2, 0, 0}, {0, 0, 0, 0, 0, 0, 0, PI2, 
   0, PI2, 0}, {0, 0, 0, 0, 0, 0, 0, 0, PI2, 0, PI2}, {0, 
   0, 0, 0, 0, 0, 0, 0, 0, PI2, 0}}, {{0, 0, PI2, 0, 0, 0, 0, 
   0, 0, 0, 0}, {0, 0, 0, PI2, 0, 0, 0, 0, 0, 0, 0}, {PI2, 0, 
   0, 0, PI2, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 0, 0, 0,
    0, 0}, {0, 0, PI2, 0, -(PI2), 0, 0, 0, 0, 0, 0}, {0, 0, 0,
    0, 0, 0, 0, PI, 0, 0, 0}, {0, 0, 0, 0, 0, 0, PI2, 0, PI2,
   0, 0}, {0, 0, 0, 0, 0, PI, 0, 0, 0, PI2, 0}, {0, 0, 0, 0,
    0, 0, PI2, 0, 0, 0, PI2}, {0, 0, 0, 0, 0, 0, 0, PI2, 
   0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, PI2, 0, 0}}, {{0, 0, 
   0, PI2, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, PI2, 0, 0, 0, 0,
    0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {PI2, 0, 0, 
   0, -(PI2), 0, 0, 0, 0, 0, 0}, {0, PI2, 0, -(PI2), 0, 0,
    0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, PI, 0, 0}, {0, 0, 0, 
   0, 0, 0, 0, PI2, 0, PI2, 0}, {0, 0, 0, 0, 0, 0, PI2, 0,
    0, 0, PI2}, {0, 0, 0, 0, 0, PI, 0, 0, 0, 0, 0}, {0, 0, 0, 
   0, 0, 0, PI2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, PI2, 0, 0,
    0}}, {{0, 0, 0, 0, PI2, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0}, {0, 0, 0, 0, -(PI2), 0, 0, 0, 0, 0, 0}, {0, 0, 
   0, -(PI2), 0, 0, 0, 0, 0, 0, 0}, {PI2, 0, -(PI2), 0, 0,
    0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, PI, 0}, {0, 0, 
   0, 0, 0, 0, 0, 0, PI2, 0, PI2}, {0, 0, 0, 0, 0, 0, 
   0, PI2, 0, 0, 0}, {0, 0, 0, 0, 0, 0, PI2, 0, 0, 0, 0}, {0, 
   0, 0, 0, 0, PI, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, PI2, 0, 
   0, 0, 0}}, {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 
   0, -(PI2), 0, 0, 0, 0, 0, 0}, {0, 0, 0, -(PI2), 0, 0, 0, 0,
    0, 0, 0}, {0, 0, -(PI2), 0, 0, 0, 0, 0, 0, 0, 
   0}, {0, -(PI2), 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, PI}, {0, 0, 0, 0, 0, 0, 0, 0, 0, PI2, 0}, {0, 0,
    0, 0, 0, 0, 0, 0, PI2, 0, 0}, {0, 0, 0, 0, 0, 0, 0, PI2, 
   0, 0, 0}, {0, 0, 0, 0, 0, 0, PI2, 0, 0, 0, 0}, {0, 0, 0, 0, 
   0, PI, 0, 0, 0, 0, 0}}};

  {
    int x1 = 1; int x2 = -1; int x3 = -2; double v1 = 2.;
    int x10 = x1 + 5;
    int x20 = x2 + 5;
    int x30 = x3 + 5;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = 1; int x2 = -2; int x3 = -3; double v1 = 4./3.;
    int x10 = x1 + 5; int x20 = x2 + 5; int x30 = x3 + 5;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = -1; int x2 = -2; int x3 = 3; double v1 = 2.;
    int x10 = x1 + 5; int x20 = x2 + 5; int x30 = x3 + 5;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = -1; int x2 = 2; int x3 = -3; double v1 = 4.;
    int x10 = x1 + 5; int x20 = x2 + 5; int x30 = x3 + 5;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = -1; int x2 = -3; int x3 = 4; double v1 = 4.;
    int x10 = x1 + 5; int x20 = x2 + 5; int x30 = x3 + 5;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = -2; int x2 = 2; int x3 = -4; double v1 = 2.;
    int x10 = x1 + 5; int x20 = x2 + 5; int x30 = x3 + 5;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = -2; int x2 = -3; int x3 = 5; double v1 = 1./0.75;
    int x10 = x1 + 5; int x20 = x2 + 5; int x30 = x3 + 5;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  return fabs(phi_integrals[m10][m20][m30]);
}

double get_idp_m4(int m1, int m2, int m3)
{
  if (m1<-4 || m1>4) { printf("\n ERROR: m1 out of range: %2i \n",m1); exit(-1); }
  if (m2<-4 || m2>4) { printf("\n ERROR: m2 out of range: %2i \n",m2); exit(-1); }
  if (m3<-4 || m3>4) { printf("\n ERROR: m3 out of range: %2i \n",m3); exit(-1); }

  int m10 = m1 + 4;
  int m20 = m2 + 4;
  int m30 = m3 + 4;

  double PI2 = PI/2.;

 //table generated in Mathematica.
 // See "Phi Integrals Phase.nb"

  double phi_integrals [9][9][9] = {{{0, 0, 0, 0, PI, 0, 0, 0, 0}, {0, 0, 0, 0, 0, PI2, 0, 0, 
   0}, {0, 0, 0, 0, 0, 0, PI2, 0, 0}, {0, 0, 0, 0, 0, 0, 0, PI/
   2, 0}, {PI, 0, 0, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, 0, 0, 0, 0,
    0, 0}, {0, 0, PI2, 0, 0, 0, 0, 0, 0}, {0, 0, 0, PI2, 0, 0,
    0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}}, {{0, 0, 0, 0, 0, PI2, 
   0, 0, 0}, {0, 0, 0, 0, PI, 0, 0, 0, 0}, {0, 0, 0, 0, 0, PI2,
    0, 0, 0}, {0, 0, 0, 0, 0, 0, PI2, 0, -(PI2)}, {0, PI, 
   0, 0, 0, 0, 0, 0, 0}, {PI2, 0, PI2, 0, 0, 0, 0, 0, 0}, {0, 
   0, 0, PI2, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 
   0, -(PI2), 0, 0, 0, 0, 0}}, {{0, 0, 0, 0, 0, 0, PI2, 0, 
   0}, {0, 0, 0, 0, 0, PI2, 0, 0, 0}, {0, 0, 0, 0, PI, 0, 0, 
   0, -(PI2)}, {0, 0, 0, 0, 0, PI2, 0, -(PI2), 0}, {0, 
   0, PI, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, PI2, 0, 0, 0, 0, 
   0}, {PI2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, -(PI2), 0, 0, 
   0, 0, 0}, {0, 0, -(PI2), 0, 0, 0, 0, 0, 0}}, {{0, 0, 0, 0, 0, 
   0, 0, PI2, 0}, {0, 0, 0, 0, 0, 0, PI2, 0, -(PI2)}, {0, 
   0, 0, 0, 0, PI2, 0, -(PI2), 0}, {0, 0, 0, 0, PI, 
   0, -(PI2), 0, 0}, {0, 0, 0, PI, 0, 0, 0, 0, 0}, {0, 
   0, PI2, 0, 0, 0, 0, 0, 0}, {0, PI2, 0, -(PI2), 0, 0, 0,
    0, 0}, {PI2, 0, -(PI2), 0, 0, 0, 0, 0, 0}, {0, -(PI2),
    0, 0, 0, 0, 0, 0, 0}}, {{PI, 0, 0, 0, 0, 0, 0, 0, 
   0}, {0, PI, 0, 0, 0, 0, 0, 0, 0}, {0, 0, PI, 0, 0, 0, 0, 0, 
   0}, {0, 0, 0, PI, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 2.*PI, 0, 0, 0,
    0}, {0, 0, 0, 0, 0, PI, 0, 0, 0}, {0, 0, 0, 0, 0, 0, PI, 0, 
   0}, {0, 0, 0, 0, 0, 0, 0, PI, 0}, {0, 0, 0, 0, 0, 0, 0, 
   0, PI}}, {{0, PI2, 0, 0, 0, 0, 0, 0, 0}, {PI2, 0, PI2,
   0, 0, 0, 0, 0, 0}, {0, PI2, 0, PI2, 0, 0, 0, 0, 0}, {0, 
   0, PI2, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, PI, 0, 0, 0}, {0,
    0, 0, 0, PI, 0, PI2, 0, 0}, {0, 0, 0, 0, 0, PI2, 
   0, PI2, 0}, {0, 0, 0, 0, 0, 0, PI2, 0, PI2}, {0, 0, 0, 
   0, 0, 0, 0, PI2, 0}}, {{0, 0, PI2, 0, 0, 0, 0, 0, 0}, {0, 
   0, 0, PI2, 0, 0, 0, 0, 0}, {PI2, 0, 0, 0, 0, 0, 0, 0, 
   0}, {0, PI2, 0, -(PI2), 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 
   0, PI, 0, 0}, {0, 0, 0, 0, 0, PI2, 0, PI2, 0}, {0, 0, 0,
    0, PI, 0, 0, 0, PI2}, {0, 0, 0, 0, 0, PI2, 0, 0, 
   0}, {0, 0, 0, 0, 0, 0, PI2, 0, 0}}, {{0, 0, 0, PI2, 0, 0, 
   0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, -(PI2), 0, 0, 
   0, 0, 0}, {PI2, 0, -(PI2), 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 
   0, 0, 0, PI, 0}, {0, 0, 0, 0, 0, 0, PI2, 0, PI2}, {0, 0,
    0, 0, 0, PI2, 0, 0, 0}, {0, 0, 0, 0, PI, 0, 0, 0, 0}, {0, 
   0, 0, 0, 0, PI2, 0, 0, 0}}, {{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 
   0, 0, -(PI2), 0, 0, 0, 0, 0}, {0, 0, -(PI2), 0, 0, 0, 0, 0,
    0}, {0, -(PI2), 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 
   0, PI}, {0, 0, 0, 0, 0, 0, 0, PI2, 0}, {0, 0, 0, 0, 0, 
   0, PI2, 0, 0}, {0, 0, 0, 0, 0, PI2, 0, 0, 0}, {0, 0, 0, 
   0, PI, 0, 0, 0, 0}}};

  //vector<vector<int> > list1;
  //vector<int> tuple1;
  //tuple1.push_back(1); tuple1.push_back(-1); tuple1.push_back
  //update_phi_int(1,-1,-2,2.,phi_integrals);

  {
    int x1 = 1; int x2 = -1; int x3 = -2; double v1 = 2.;
    int x10 = x1 + 4;
    int x20 = x2 + 4;
    int x30 = x3 + 4;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = 1; int x2 = -2; int x3 = -3; double v1 = 4./3.;
    int x10 = x1 + 4; int x20 = x2 + 4; int x30 = x3 + 4;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = -1; int x2 = -2; int x3 = 3; double v1 = 2.;
    int x10 = x1 + 4; int x20 = x2 + 4; int x30 = x3 + 4;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = -1; int x2 = 2; int x3 = -3; double v1 = 4.;
    int x10 = x1 + 4; int x20 = x2 + 4; int x30 = x3 + 4;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = -1; int x2 = -3; int x3 = 4; double v1 = 4.;
    int x10 = x1 + 4; int x20 = x2 + 4; int x30 = x3 + 4;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  {
    int x1 = -2; int x2 = 2; int x3 = -4; double v1 = 2.;
    int x10 = x1 + 4; int x20 = x2 + 4; int x30 = x3 + 4;
    phi_integrals[x10][x20][x30] *= v1;
    phi_integrals[x10][x30][x20] *= v1;
    phi_integrals[x20][x10][x30] *= v1;
    phi_integrals[x20][x30][x10] *= v1;
    phi_integrals[x30][x10][x20] *= v1;
    phi_integrals[x30][x20][x10] *= v1;
  }

  //return phi_integrals[m10][m20][m30];
  return fabs(phi_integrals[m10][m20][m30]);
}

#pragma acc routine seq
double ps_dV(double mu1, double mu2, double nu1, double nu2)
{
  double c3n = cos(3*nu1) - cos(3*nu2);
  double cm = cosh(mu1) - cosh(mu2);
  double cn = cos(nu1) - cos(nu2);
  double c3m = cosh(3*mu1) - cosh(3*mu2);
  double val = ( c3n * cm - cn*c3m ) / 12.;
  return val;
}

