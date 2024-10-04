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

//CPMZ marked: nu2 bug in convert function
double second_order_fx2ordV(double Z1, double Z2, double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc)
{
  double f = 1./27648.;
  double pa4 = pow(a,4.);
  double fp = f*pa4;

  double pm02 = pow(mu0,2.);
  double pm12 = pow(mu1,2.);
  double pm22 = pow(mu2,2.);
  double pn02 = pow(nu0,2.);
  double pn12 = pow(nu1,2.);
  double pn22 = pow(nu2,2.);

  double pm012 = pow(mu0-mu1,2.);
  double pm022 = pow(mu0-mu2,2.);
  double pn012 = pow(nu0-nu1,2.);
  double pn022 = pow(nu0-nu2,2.);

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double cosh2m1 = cosh(2.*mu1);
  double cosh2m2 = cosh(2.*mu2);
  double sinh2m1 = sinh(2.*mu1);
  double sinh2m2 = sinh(2.*mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosh4m1 = cosh(4.*mu1);
  double cosh4m2 = cosh(4.*mu2);
  double sinh4m1 = sinh(4.*mu1);
  double sinh4m2 = sinh(4.*mu2);

  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cos2n1 = cos(2.*nu1);
  double cos2n2 = cos(2.*nu2);
  double sin2n1 = sin(2.*nu1);
  double sin2n2 = sin(2.*nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);
  double sin3n1 = sin(3.*nu1);
  double sin3n2 = sin(3.*nu2);
  double cos4n1 = cos(4.*nu1);
  double cos4n2 = cos(4.*nu2);
  double sin4n1 = sin(4.*nu1);
  double sin4n2 = sin(4.*nu2);

  double nnbc = 9.*bc;
  double nnb1 = 9.*b1;
  double nnb2 = 9.*b2;
  double nnc1 = 9.*c1;
  double nnc2 = 9.*c2;

  double sxbc = 16.*bc;
  double sxb1 = 16.*b1;
  double sxb2 = 16.*b2;
  double sxc1 = 16.*c1;
  double sxc2 = 16.*c2;

  double sfbc = 64.*bc;
  double sfb1 = 64.*b1;
  double sfb2 = 64.*b2;
  double sfc1 = 64.*c1;
  double sfc2 = 64.*c2;

  double stbc = 72.*bc;
  double stb1 = 72.*b1;
  double stb2 = 72.*b2;
  double stc1 = 72.*c1;
  double stc2 = 72.*c2;

  double ofbc = 144.*bc;
  double ofb1 = 144.*b1;
  double ofb2 = 144.*b2;
  double ofc1 = 144.*c1;
  double ofc2 = 144.*c2;

  double eebc = 288.*bc;
  double eeb1 = 288.*b1;
  double eeb2 = 288.*b2;
  double eec1 = 288.*c1;
  double eec2 = 288.*c2;

  double fnbc = 576.*bc;
  double fnb1 = 576.*b1;
  double fnb2 = 576.*b2;
  double fnc1 = 576.*c1;
  double fnc2 = 576.*c2;

  double xebc = 648.*bc;
  double xeb1 = 648.*b1;
  double xeb2 = 648.*b2;
  double xec1 = 648.*c1;
  double xec2 = 648.*c2;

  double tnbc = 1296.*bc;
  double tnb1 = 1296.*b1;
  double tnb2 = 1296.*b2;
  double tnc1 = 1296.*c1;
  double tnc2 = 1296.*c2;

  double val1 = fp*(288*a0*cos4n1*cosh3m1 + sfb2*cos4n1*cosh3m1 - ofc2*cos4n1*cosh3m1 - eeb1*mu0*cos4n1*cosh3m1 + 
       eeb2*pm02*cos4n1*cosh3m1 + eeb1*mu1*cos4n1*cosh3m1 - fnb2*mu0*mu1*cos4n1*cosh3m1 + eeb2*pm12*cos4n1*cosh3m1 - 
       eec1*nu0*cos4n1*cosh3m1 + eebc*mu0*nu0*cos4n1*cosh3m1 - eebc*mu1*nu0*cos4n1*cosh3m1 + eec2*pn02*cos4n1*cosh3m1 + 
       eec1*nu1*cos4n1*cosh3m1 - eebc*mu0*nu1*cos4n1*cosh3m1 + eebc*mu1*nu1*cos4n1*cosh3m1 - fnc2*nu0*nu1*cos4n1*cosh3m1 + 
       eec2*pn12*cos4n1*cosh3m1 - 72*a0*cos4n1*cosh3m1 - sxb2*cos4n1*cosh3m1 + nnc2*cos4n1*cosh3m1 + 
       stb1*mu0*cos4n1*cosh3m1 - stb2*pm02*cos4n1*cosh3m1 - stb1*mu1*cos4n1*cosh3m1 + ofb2*mu0*mu1*cos4n1*cosh3m1 - 
       stb2*pm12*cos4n1*cosh3m1 + stc1*nu0*cos4n1*cosh3m1 - stbc*mu0*nu0*cos4n1*cosh3m1 + stbc*mu1*nu0*cos4n1*cosh3m1 - 
       stc2*pn02*cos4n1*cosh3m1 - stc1*nu1*cos4n1*cosh3m1 + stbc*mu0*nu1*cos4n1*cosh3m1 - stbc*mu1*nu1*cos4n1*cosh3m1 + 
       ofc2*nu0*nu1*cos4n1*cosh3m1 - stc2*pn12*cos4n1*cosh3m1 - 288*a0*cos4n2*cosh3m1 - sfb2*cos4n2*cosh3m1 + 
       ofc2*cos4n2*cosh3m1 + eeb1*mu0*cos4n2*cosh3m1 - eeb2*pm02*cos4n2*cosh3m1 - eeb1*mu1*cos4n2*cosh3m1 + 
       fnb2*mu0*mu1*cos4n2*cosh3m1 - eeb2*pm12*cos4n2*cosh3m1 + eec1*nu0*cos4n2*cosh3m1 - eebc*mu0*nu0*cos4n2*cosh3m1 + 
       eebc*mu1*nu0*cos4n2*cosh3m1 - eec2*pn02*cos4n2*cosh3m1 - eec1*nu1*cos4n2*cosh3m1 + eebc*mu0*nu1*cos4n2*cosh3m1 - 
       eebc*mu1*nu1*cos4n2*cosh3m1 + fnc2*nu0*nu1*cos4n2*cosh3m1 - eec2*pn22*cos4n2*cosh3m1 + 72*a0*cos4n2*cosh3m1 + 
       sxb2*cos4n2*cosh3m1 - nnc2*cos4n2*cosh3m1 - stb1*mu0*cos4n2*cosh3m1 + stb2*pm02*cos4n2*cosh3m1 + 
       stb1*mu1*cos4n2*cosh3m1 - ofb2*mu0*mu1*cos4n2*cosh3m1 + stb2*pm12*cos4n2*cosh3m1 - stc1*nu0*cos4n2*cosh3m1 + 
       stbc*mu0*nu0*cos4n2*cosh3m1 - stbc*mu1*nu0*cos4n2*cosh3m1 + stc2*pn02*cos4n2*cosh3m1 + stc1*nu1*cos4n2*cosh3m1 - 
       stbc*mu0*nu1*cos4n2*cosh3m1 + stbc*mu1*nu1*cos4n2*cosh3m1 - ofc2*nu0*nu1*cos4n2*cosh3m1 + stc2*pn22*cos4n2*cosh3m1 - 
       648*a0*cosn1*cosh4m1 - 81*b2*cosn1*cosh4m1 + tnc2*cosn1*cosh4m1 + xeb1*mu0*cosn1*cosh4m1 - 
       xeb2*pm02*cosn1*cosh4m1 - xeb1*mu1*cosn1*cosh4m1 + tnb2*mu0*mu1*cosn1*cosh4m1 - xeb2*pm12*cosn1*cosh4m1 + 
       xec1*nu0*cosn1*cosh4m1 - xebc*mu0*nu0*cosn1*cosh4m1 + xebc*mu1*nu0*cosn1*cosh4m1 - xec2*pn02*cosn1*cosh4m1 - 
       xec1*nu1*cosn1*cosh4m1 + xebc*mu0*nu1*cosn1*cosh4m1 - xebc*mu1*nu1*cosn1*cosh4m1 + tnc2*nu0*nu1*cosn1*cosh4m1 - 
       xec2*pn12*cosn1*cosh4m1 + 72*a0*cos3n1*cosh4m1 + nnb2*cos3n1*cosh4m1 - sxc2*cos3n1*cosh4m1 - 
       stb1*mu0*cos3n1*cosh4m1 + stb2*pm02*cos3n1*cosh4m1 + stb1*mu1*cos3n1*cosh4m1 - ofb2*mu0*mu1*cos3n1*cosh4m1 + 
       stb2*pm12*cos3n1*cosh4m1 - stc1*nu0*cos3n1*cosh4m1 + stbc*mu0*nu0*cos3n1*cosh4m1 - stbc*mu1*nu0*cos3n1*cosh4m1 + 
       stc2*pn02*cos3n1*cosh4m1 + stc1*nu1*cos3n1*cosh4m1 - stbc*mu0*nu1*cos3n1*cosh4m1 + stbc*mu1*nu1*cos3n1*cosh4m1 - 
       ofc2*nu0*nu1*cos3n1*cosh4m1 + stc2*pn12*cos3n1*cosh4m1 + 648*a0*cosn2*cosh4m1 + 81*b2*cosn2*cosh4m1 - 
       tnc2*cosn2*cosh4m1 - xeb1*mu0*cosn2*cosh4m1 + xeb2*pm02*cosn2*cosh4m1 + xeb1*mu1*cosn2*cosh4m1 - 
       tnb2*mu0*mu1*cosn2*cosh4m1 + xeb2*pm12*cosn2*cosh4m1 - xec1*nu0*cosn2*cosh4m1 + xebc*mu0*nu0*cosn2*cosh4m1 - 
       xebc*mu1*nu0*cosn2*cosh4m1 + xec2*pn02*cosn2*cosh4m1 + xec1*nu1*cosn2*cosh4m1 - xebc*mu0*nu1*cosn2*cosh4m1 + 
       xebc*mu1*nu1*cosn2*cosh4m1 - tnc2*nu0*nu1*cosn2*cosh4m1 + xec2*pn22*cosn2*cosh4m1 - 72*a0*cos3n2*cosh4m1 - 
       nnb2*cos3n2*cosh4m1 + sxc2*cos3n2*cosh4m1 + stb1*mu0*cos3n2*cosh4m1 - stb2*pm02*cos3n2*cosh4m1 - 
       stb1*mu1*cos3n2*cosh4m1 + ofb2*mu0*mu1*cos3n2*cosh4m1 - stb2*pm12*cos3n2*cosh4m1 + stc1*nu0*cos3n2*cosh4m1 - 
       stbc*mu0*nu0*cos3n2*cosh4m1 + stbc*mu1*nu0*cos3n2*cosh4m1 - stc2*pn02*cos3n2*cosh4m1 - stc1*nu1*cos3n2*cosh4m1 + 
       stbc*mu0*nu1*cos3n2*cosh4m1 - stbc*mu1*nu1*cos3n2*cosh4m1 + ofc2*nu0*nu1*cos3n2*cosh4m1 - stc2*pn22*cos3n2*cosh4m1 + 
       2592*a0*cos4n1*coshm2 + 5184*b2*cos4n1*coshm2 - tnc2*cos4n1*coshm2 - 2592*b1*mu0*cos4n1*coshm2 + 
       2592*b2*pm02*cos4n1*coshm2 + 2592*b1*mu2*cos4n1*coshm2 - 5184*b2*mu0*mu2*cos4n1*coshm2 + 2592*b2*pm22*cos4n1*coshm2 - 
       2592*c1*nu0*cos4n1*coshm2 + 2592*bc*mu0*nu0*cos4n1*coshm2 - 2592*bc*mu2*nu0*cos4n1*coshm2 + 2592*c2*pn02*cos4n1*coshm2 + 
       2592*c1*nu1*cos4n1*coshm2 - 2592*bc*mu0*nu1*cos4n1*coshm2 + 2592*bc*mu2*nu1*cos4n1*coshm2 - 5184*c2*nu0*nu1*cos4n1*coshm2 + 
       2592*c2*pn12*cos4n1*coshm2 - 648*a0*cos4n1*coshm2 - tnb2*cos4n1*coshm2 + 81*c2*cos4n1*coshm2 + 
       xeb1*mu0*cos4n1*coshm2 - xeb2*pm02*cos4n1*coshm2 - xeb1*mu2*cos4n1*coshm2 + tnb2*mu0*mu2*cos4n1*coshm2 - 
       xeb2*pm22*cos4n1*coshm2 + xec1*nu0*cos4n1*coshm2 - xebc*mu0*nu0*cos4n1*coshm2 + xebc*mu2*nu0*cos4n1*coshm2 - 
       xec2*pn02*cos4n1*coshm2 - xec1*nu1*cos4n1*coshm2 + xebc*mu0*nu1*cos4n1*coshm2 - xebc*mu2*nu1*cos4n1*coshm2 + 
       tnc2*nu0*nu1*cos4n1*coshm2 - xec2*pn12*cos4n1*coshm2 - 2592*a0*cos4n2*coshm2 - 5184*b2*cos4n2*coshm2 + 
       tnc2*cos4n2*coshm2 + 2592*b1*mu0*cos4n2*coshm2 - 2592*b2*pm02*cos4n2*coshm2 - 2592*b1*mu2*cos4n2*coshm2 + 
       5184*b2*mu0*mu2*cos4n2*coshm2 - 2592*b2*pm22*cos4n2*coshm2 + 2592*c1*nu0*cos4n2*coshm2 - 2592*bc*mu0*nu0*cos4n2*coshm2 + 
       2592*bc*mu2*nu0*cos4n2*coshm2 - 2592*c2*pn02*cos4n2*coshm2 - 2592*c1*nu1*cos4n2*coshm2 + 2592*bc*mu0*nu1*cos4n2*coshm2 - 
       2592*bc*mu2*nu1*cos4n2*coshm2 + 5184*c2*nu0*nu1*cos4n2*coshm2 - 2592*c2*pn22*cos4n2*coshm2 + 648*a0*cos4n2*coshm2 + 
       tnb2*cos4n2*coshm2 - 81*c2*cos4n2*coshm2 - xeb1*mu0*cos4n2*coshm2 + xeb2*pm02*cos4n2*coshm2 + 
       xeb1*mu2*cos4n2*coshm2 - tnb2*mu0*mu2*cos4n2*coshm2 + xeb2*pm22*cos4n2*coshm2 - xec1*nu0*cos4n2*coshm2 + 
       xebc*mu0*nu0*cos4n2*coshm2 - xebc*mu2*nu0*cos4n2*coshm2 + xec2*pn02*cos4n2*coshm2 + xec1*nu1*cos4n2*coshm2 - 
       xebc*mu0*nu1*cos4n2*coshm2 + xebc*mu2*nu1*cos4n2*coshm2 - tnc2*nu0*nu1*cos4n2*coshm2 + xec2*pn22*cos4n2*coshm2 - 
       2592*a0*cosn1*cosh2m2 - tnb2*cosn1*cosh2m2 + 5184*c2*cosn1*cosh2m2 + 2592*b1*mu0*cosn1*cosh2m2 - 
       2592*b2*pm02*cosn1*cosh2m2 - 2592*b1*mu2*cosn1*cosh2m2 + 5184*b2*mu0*mu2*cosn1*cosh2m2 - 2592*b2*pm22*cosn1*cosh2m2 + 
       2592*c1*nu0*cosn1*cosh2m2 - 2592*bc*mu0*nu0*cosn1*cosh2m2 + 2592*bc*mu2*nu0*cosn1*cosh2m2 - 2592*c2*pn02*cosn1*cosh2m2 - 
       2592*c1*nu1*cosn1*cosh2m2 + 2592*bc*mu0*nu1*cosn1*cosh2m2 - 2592*bc*mu2*nu1*cosn1*cosh2m2 + 5184*c2*nu0*nu1*cosn1*cosh2m2 - 
       2592*c2*pn12*cosn1*cosh2m2 + 288*a0*cos3n1*cosh2m2 + ofb2*cos3n1*cosh2m2 - sfc2*cos3n1*cosh2m2 - 
       eeb1*mu0*cos3n1*cosh2m2 + eeb2*pm02*cos3n1*cosh2m2 + eeb1*mu2*cos3n1*cosh2m2 - fnb2*mu0*mu2*cos3n1*cosh2m2 + 
       eeb2*pm22*cos3n1*cosh2m2 - eec1*nu0*cos3n1*cosh2m2 + eebc*mu0*nu0*cos3n1*cosh2m2 - eebc*mu2*nu0*cos3n1*cosh2m2 + 
       eec2*pn02*cos3n1*cosh2m2 + eec1*nu1*cos3n1*cosh2m2 - eebc*mu0*nu1*cos3n1*cosh2m2 + eebc*mu2*nu1*cos3n1*cosh2m2 - 
       fnc2*nu0*nu1*cos3n1*cosh2m2 + eec2*pn12*cos3n1*cosh2m2 + 2592*a0*cosn2*cosh2m2 + tnb2*cosn2*cosh2m2 - 
       5184*c2*cosn2*cosh2m2 - 2592*b1*mu0*cosn2*cosh2m2 + 2592*b2*pm02*cosn2*cosh2m2 + 2592*b1*mu2*cosn2*cosh2m2 - 
       5184*b2*mu0*mu2*cosn2*cosh2m2 + 2592*b2*pm22*cosn2*cosh2m2 - 2592*c1*nu0*cosn2*cosh2m2 + 2592*bc*mu0*nu0*cosn2*cosh2m2 - 
       2592*bc*mu2*nu0*cosn2*cosh2m2 + 2592*c2*pn02*cosn2*cosh2m2 + 2592*c1*nu1*cosn2*cosh2m2 - 2592*bc*mu0*nu1*cosn2*cosh2m2 + 
       2592*bc*mu2*nu1*cosn2*cosh2m2 - 5184*c2*nu0*nu1*cosn2*cosh2m2 + 2592*c2*pn22*cosn2*cosh2m2 - 288*a0*cos3n2*cosh2m2 - 
       ofb2*cos3n2*cosh2m2 + sfc2*cos3n2*cosh2m2 + eeb1*mu0*cos3n2*cosh2m2 - eeb2*pm02*cos3n2*cosh2m2 - 
       eeb1*mu2*cos3n2*cosh2m2 + fnb2*mu0*mu2*cos3n2*cosh2m2 - eeb2*pm22*cos3n2*cosh2m2 + eec1*nu0*cos3n2*cosh2m2 - 
       eebc*mu0*nu0*cos3n2*cosh2m2 + eebc*mu2*nu0*cos3n2*cosh2m2 - eec2*pn02*cos3n2*cosh2m2 - eec1*nu1*cos3n2*cosh2m2 + 
       eebc*mu0*nu1*cos3n2*cosh2m2 - eebc*mu2*nu1*cos3n2*cosh2m2 + fnc2*nu0*nu1*cos3n2*cosh2m2 - eec2*pn22*cos3n2*cosh2m2 - 
       288*a0*cos4n1*cosh3m2 - sfb2*cos4n1*cosh3m2 + ofc2*cos4n1*cosh3m2 + eeb1*mu0*cos4n1*cosh3m2 - 
       eeb2*pm02*cos4n1*cosh3m2 - eeb1*mu2*cos4n1*cosh3m2 + fnb2*mu0*mu2*cos4n1*cosh3m2 - eeb2*pm22*cos4n1*cosh3m2 + 
       eec1*nu0*cos4n1*cosh3m2 - eebc*mu0*nu0*cos4n1*cosh3m2 + eebc*mu2*nu0*cos4n1*cosh3m2 - eec2*pn02*cos4n1*cosh3m2 - 
       eec1*nu1*cos4n1*cosh3m2 + eebc*mu0*nu1*cos4n1*cosh3m2 - eebc*mu2*nu1*cos4n1*cosh3m2 + fnc2*nu0*nu1*cos4n1*cosh3m2 - 
       eec2*pn12*cos4n1*cosh3m2 + 72*a0*cos4n1*cosh3m2 + sxb2*cos4n1*cosh3m2 - nnc2*cos4n1*cosh3m2 - 
       stb1*mu0*cos4n1*cosh3m2 + stb2*pm02*cos4n1*cosh3m2 + stb1*mu2*cos4n1*cosh3m2 - ofb2*mu0*mu2*cos4n1*cosh3m2 + 
       stb2*pm22*cos4n1*cosh3m2 - stc1*nu0*cos4n1*cosh3m2 + stbc*mu0*nu0*cos4n1*cosh3m2 - stbc*mu2*nu0*cos4n1*cosh3m2 + 
       stc2*pn02*cos4n1*cosh3m2 + stc1*nu1*cos4n1*cosh3m2 - stbc*mu0*nu1*cos4n1*cosh3m2 + stbc*mu2*nu1*cos4n1*cosh3m2 - 
       ofc2*nu0*nu1*cos4n1*cosh3m2 + stc2*pn12*cos4n1*cosh3m2 + 288*a0*cos4n2*cosh3m2 + sfb2*cos4n2*cosh3m2 - 
       ofc2*cos4n2*cosh3m2 - eeb1*mu0*cos4n2*cosh3m2 + eeb2*pm02*cos4n2*cosh3m2 + eeb1*mu2*cos4n2*cosh3m2 - 
       fnb2*mu0*mu2*cos4n2*cosh3m2 + eeb2*pm22*cos4n2*cosh3m2 - eec1*nu0*cos4n2*cosh3m2 + eebc*mu0*nu0*cos4n2*cosh3m2 - 
       eebc*mu2*nu0*cos4n2*cosh3m2 + eec2*pn02*cos4n2*cosh3m2 + eec1*nu1*cos4n2*cosh3m2 - eebc*mu0*nu1*cos4n2*cosh3m2 + 
       eebc*mu2*nu1*cos4n2*cosh3m2 - fnc2*nu0*nu1*cos4n2*cosh3m2 + eec2*pn22*cos4n2*cosh3m2 - 72*a0*cos4n2*cosh3m2 - 
       sxb2*cos4n2*cosh3m2 + nnc2*cos4n2*cosh3m2 + stb1*mu0*cos4n2*cosh3m2 - stb2*pm02*cos4n2*cosh3m2 - 
       stb1*mu2*cos4n2*cosh3m2 + ofb2*mu0*mu2*cos4n2*cosh3m2 - stb2*pm22*cos4n2*cosh3m2 + stc1*nu0*cos4n2*cosh3m2 - 
       stbc*mu0*nu0*cos4n2*cosh3m2 + stbc*mu2*nu0*cos4n2*cosh3m2 - stc2*pn02*cos4n2*cosh3m2 - stc1*nu1*cos4n2*cosh3m2 + 
       stbc*mu0*nu1*cos4n2*cosh3m2 - stbc*mu2*nu1*cos4n2*cosh3m2 + ofc2*nu0*nu1*cos4n2*cosh3m2 - stc2*pn22*cos4n2*cosh3m2 + 
       648*a0*cosn1*cosh4m2 + 81*b2*cosn1*cosh4m2 - tnc2*cosn1*cosh4m2 - xeb1*mu0*cosn1*cosh4m2 + 
       xeb2*pm02*cosn1*cosh4m2 + xeb1*mu2*cosn1*cosh4m2 - tnb2*mu0*mu2*cosn1*cosh4m2 + xeb2*pm22*cosn1*cosh4m2 - 
       xec1*nu0*cosn1*cosh4m2 + xebc*mu0*nu0*cosn1*cosh4m2 - xebc*mu2*nu0*cosn1*cosh4m2 + xec2*pn02*cosn1*cosh4m2 + 
       xec1*nu1*cosn1*cosh4m2 - xebc*mu0*nu1*cosn1*cosh4m2 + xebc*mu2*nu1*cosn1*cosh4m2 - tnc2*nu0*nu1*cosn1*cosh4m2 + 
       xec2*pn12*cosn1*cosh4m2 - 72*a0*cos3n1*cosh4m2 - nnb2*cos3n1*cosh4m2 + sxc2*cos3n1*cosh4m2 + 
       stb1*mu0*cos3n1*cosh4m2 - stb2*pm02*cos3n1*cosh4m2 - stb1*mu2*cos3n1*cosh4m2 + ofb2*mu0*mu2*cos3n1*cosh4m2 - 
       stb2*pm22*cos3n1*cosh4m2 + stc1*nu0*cos3n1*cosh4m2 - stbc*mu0*nu0*cos3n1*cosh4m2 + stbc*mu2*nu0*cos3n1*cosh4m2 - 
       stc2*pn02*cos3n1*cosh4m2 - stc1*nu1*cos3n1*cosh4m2 + stbc*mu0*nu1*cos3n1*cosh4m2 - stbc*mu2*nu1*cos3n1*cosh4m2 + 
       ofc2*nu0*nu1*cos3n1*cosh4m2 - stc2*pn12*cos3n1*cosh4m2 - 648*a0*cosn2*cosh4m2 - 81*b2*cosn2*cosh4m2 + 
       tnc2*cosn2*cosh4m2 + xeb1*mu0*cosn2*cosh4m2 - xeb2*pm02*cosn2*cosh4m2 - xeb1*mu2*cosn2*cosh4m2 + 
       tnb2*mu0*mu2*cosn2*cosh4m2 - xeb2*pm22*cosn2*cosh4m2 + xec1*nu0*cosn2*cosh4m2 - xebc*mu0*nu0*cosn2*cosh4m2 + 
       xebc*mu2*nu0*cosn2*cosh4m2 - xec2*pn02*cosn2*cosh4m2 - xec1*nu1*cosn2*cosh4m2 + xebc*mu0*nu1*cosn2*cosh4m2 - 
       xebc*mu2*nu1*cosn2*cosh4m2 + tnc2*nu0*nu1*cosn2*cosh4m2 - xec2*pn22*cosn2*cosh4m2 + 72*a0*cos3n2*cosh4m2 + 
       nnb2*cos3n2*cosh4m2 - sxc2*cos3n2*cosh4m2 - stb1*mu0*cos3n2*cosh4m2 + stb2*pm02*cos3n2*cosh4m2 + 
       stb1*mu2*cos3n2*cosh4m2 - ofb2*mu0*mu2*cos3n2*cosh4m2 + stb2*pm22*cos3n2*cosh4m2 - stc1*nu0*cos3n2*cosh4m2 + 
       stbc*mu0*nu0*cos3n2*cosh4m2 - stbc*mu2*nu0*cos3n2*cosh4m2 + stc2*pn02*cos3n2*cosh4m2 + stc1*nu1*cos3n2*cosh4m2 - 
       stbc*mu0*nu1*cos3n2*cosh4m2 + stbc*mu2*nu1*cos3n2*cosh4m2 - ofc2*nu0*nu1*cos3n2*cosh4m2 + stc2*pn22*cos3n2*cosh4m2 + 
       xec1*cosh4m1*sinn1 - xebc*mu0*cosh4m1*sinn1 + xebc*mu1*cosh4m1*sinn1 - tnc2*nu0*cosh4m1*sinn1 + 
       tnc2*nu1*cosh4m1*sinn1 + 2592*c1*cosh2m2*sinn1 - 2592*bc*mu0*cosh2m2*sinn1 + 2592*bc*mu2*cosh2m2*sinn1 - 
       5184*c2*nu0*cosh2m2*sinn1 + 5184*c2*nu1*cosh2m2*sinn1 - xec1*cosh4m2*sinn1 + xebc*mu0*cosh4m2*sinn1 - 
       xebc*mu2*cosh4m2*sinn1 + tnc2*nu0*cosh4m2*sinn1 - tnc2*nu1*cosh4m2*sinn1 - ofc1*cosh3m1*sin4n1 + 
       ofbc*mu0*cosh3m1*sin4n1 - ofbc*mu1*cosh3m1*sin4n1 + eec2*nu0*cosh3m1*sin4n1 - eec2*nu1*cosh3m1*sin4n1 - 
       tnc1*coshm2*sin4n1 + tnbc*mu0*coshm2*sin4n1 - tnbc*mu2*coshm2*sin4n1 + 2592*c2*nu0*coshm2*sin4n1 - 
       2592*c2*nu1*coshm2*sin4n1 + ofc1*cosh3m2*sin4n1 - ofbc*mu0*cosh3m2*sin4n1 + ofbc*mu2*cosh3m2*sin4n1 - 
       eec2*nu0*cosh3m2*sin4n1 + eec2*nu1*cosh3m2*sin4n1 - 24*c1*cosh4m1*sin3n1 + 24*bc*mu0*cosh4m1*sin3n1 - 
       24*bc*mu1*cosh4m1*sin3n1 + 48*c2*nu0*cosh4m1*sin3n1 - 48*c2*nu1*cosh4m1*sin3n1 - 96*c1*cosh2m2*sin3n1 + 
       96*bc*mu0*cosh2m2*sin3n1 - 96*bc*mu2*cosh2m2*sin3n1 + 192*c2*nu0*cosh2m2*sin3n1 - 192*c2*nu1*cosh2m2*sin3n1 + 
       24*c1*cosh4m2*sin3n1 - 24*bc*mu0*cosh4m2*sin3n1 + 24*bc*mu2*cosh4m2*sin3n1 - 48*c2*nu0*cosh4m2*sin3n1 + 
       48*c2*nu1*cosh4m2*sin3n1 + 18*c1*cosh3m1*sin4n1 - 18*bc*mu0*cosh3m1*sin4n1 + 18*bc*mu1*cosh3m1*sin4n1 - 
       36*c2*nu0*cosh3m1*sin4n1 + 36*c2*nu1*cosh3m1*sin4n1 + 162*c1*coshm2*sin4n1 - 162*bc*mu0*coshm2*sin4n1 + 
       162*bc*mu2*coshm2*sin4n1 - 324*c2*nu0*coshm2*sin4n1 + 324*c2*nu1*coshm2*sin4n1 - 18*c1*cosh3m2*sin4n1 + 
       18*bc*mu0*cosh3m2*sin4n1 - 18*bc*mu2*cosh3m2*sin4n1 + 36*c2*nu0*cosh3m2*sin4n1 - 36*c2*nu1*cosh3m2*sin4n1 - 
       xec1*cosh4m1*sinn2 + xebc*mu0*cosh4m1*sinn2 - xebc*mu1*cosh4m1*sinn2 + tnc2*nu0*cosh4m1*sinn2 - 
       tnc2*nu1*cosh4m1*sinn2 - 2592*c1*cosh2m2*sinn2 + 2592*bc*mu0*cosh2m2*sinn2 - 2592*bc*mu2*cosh2m2*sinn2 + 
       5184*c2*nu0*cosh2m2*sinn2 - 5184*c2*nu1*cosh2m2*sinn2 + xec1*cosh4m2*sinn2 - xebc*mu0*cosh4m2*sinn2 + 
       xebc*mu2*cosh4m2*sinn2 - tnc2*nu0*cosh4m2*sinn2 + tnc2*nu1*cosh4m2*sinn2 + ofc1*cosh3m1*sin4n2 - 
       ofbc*mu0*cosh3m1*sin4n2 + ofbc*mu1*cosh3m1*sin4n2 - eec2*nu0*cosh3m1*sin4n2 + eec2*nu1*cosh3m1*sin4n2 + 
       tnc1*coshm2*sin4n2 - tnbc*mu0*coshm2*sin4n2 + tnbc*mu2*coshm2*sin4n2 - 2592*c2*nu0*coshm2*sin4n2 + 
       2592*c2*nu1*coshm2*sin4n2 - ofc1*cosh3m2*sin4n2 + ofbc*mu0*cosh3m2*sin4n2 - ofbc*mu2*cosh3m2*sin4n2 + 
       eec2*nu0*cosh3m2*sin4n2 - eec2*nu1*cosh3m2*sin4n2 + 24*c1*cosh4m1*sin3n2 - 24*bc*mu0*cosh4m1*sin3n2 + 
       24*bc*mu1*cosh4m1*sin3n2 - 48*c2*nu0*cosh4m1*sin3n2 + 48*c2*nu1*cosh4m1*sin3n2 + 96*c1*cosh2m2*sin3n2 - 
       96*bc*mu0*cosh2m2*sin3n2 + 96*bc*mu2*cosh2m2*sin3n2 - 192*c2*nu0*cosh2m2*sin3n2 + 192*c2*nu1*cosh2m2*sin3n2 - 
       24*c1*cosh4m2*sin3n2 + 24*bc*mu0*cosh4m2*sin3n2 - 24*bc*mu2*cosh4m2*sin3n2 + 48*c2*nu0*cosh4m2*sin3n2 - 
       48*c2*nu1*cosh4m2*sin3n2 + 16*cosh2m1*(81*(2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn012) + 2*c1*(-nu0 + nu1) - 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 - 18*a0*cos3n1 - nnb2*cos3n1 + 4*c2*cos3n1 + 18*b1*mu0*cos3n1 - 18*b2*pm02*cos3n1 - 
          18*b1*mu1*cos3n1 + 36*b2*mu0*mu1*cos3n1 - 18*b2*pm12*cos3n1 + 18*c1*nu0*cos3n1 - 18*bc*mu0*nu0*cos3n1 + 18*bc*mu1*nu0*cos3n1 - 
          18*c2*pn02*cos3n1 - 18*c1*nu1*cos3n1 + 18*bc*mu0*nu1*cos3n1 - 18*bc*mu1*nu1*cos3n1 + 36*c2*nu0*nu1*cos3n1 - 
          18*c2*pn12*cos3n1 - 81*(2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn022) + 2*c1*(-nu0 + nu2) - 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 + (18*a0 + nnb2*(1 + 2*pm012) + 
             2*(nnb1*(-mu0 + mu1) + c2*(-2 + 9*pn022) + 9*(c1 + bc*(-mu0 + mu1))*(-nu0 + nu2)))*cos3n2 - 162*c1*sinn1 + 162*bc*mu0*sinn1 - 
          162*bc*mu1*sinn1 + 324*c2*nu0*sinn1 - 324*c2*nu1*sinn1 + 6*c1*sin3n1 - 6*bc*mu0*sin3n1 + 6*bc*mu1*sin3n1 - 12*c2*nu0*sin3n1 + 
          12*c2*nu1*sin3n1 + 162*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sinn2 - 6*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin3n2) - 
       18*c1*cosh3m1*sin4n2 + 18*bc*mu0*cosh3m1*sin4n2 - 18*bc*mu1*cosh3m1*sin4n2 + 36*c2*nu0*cosh3m1*sin4n2 - 
       36*c2*nu1*cosh3m1*sin4n2 - 162*c1*coshm2*sin4n2 + 162*bc*mu0*coshm2*sin4n2 - 162*bc*mu2*coshm2*sin4n2 + 324*c2*nu0*coshm2*sin4n2 - 
       324*c2*nu1*coshm2*sin4n2 + 18*c1*cosh3m2*sin4n2 - 18*bc*mu0*cosh3m2*sin4n2 + 18*bc*mu2*cosh3m2*sin4n2 - 
       36*c2*nu0*cosh3m2*sin4n2 + 36*c2*nu1*cosh3m2*sin4n2 + 
       81*coshm1*(-16*(2*a0 + 2*b2*(2 + pm012) + c2*(-1 + 2*pn012) + 2*c1*(-nu0 + nu1) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos4n1 + 
          8*a0*cos4n1 + sxb2*cos4n1 - c2*cos4n1 - 8*b1*mu0*cos4n1 + 8*b2*pm02*cos4n1 + 8*b1*mu1*cos4n1 - sxb2*mu0*mu1*cos4n1 + 
          8*b2*pm12*cos4n1 - 8*c1*nu0*cos4n1 + 8*bc*mu0*nu0*cos4n1 - 8*bc*mu1*nu0*cos4n1 + 8*c2*pn02*cos4n1 + 8*c1*nu1*cos4n1 - 
          8*bc*mu0*nu1*cos4n1 + 8*bc*mu1*nu1*cos4n1 - sxc2*nu0*nu1*cos4n1 + 8*c2*pn12*cos4n1 + 
          16*(2*a0 + 2*b2*(2 + pm012) + c2*(-1 + 2*pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2 + 
          (-8*a0 - 8*b2*(2 + pm012) + c2*(1 - 8*pn022) + 8*c1*(nu0 - nu2) + 8*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2 + sxc1*sin4n1 - 
          sxbc*mu0*sin4n1 + sxbc*mu1*sin4n1 - 32*c2*nu0*sin4n1 + 32*c2*nu1*sin4n1 - 2*c1*sin4n1 + 2*bc*mu0*sin4n1 - 2*bc*mu1*sin4n1 + 
          4*c2*nu0*sin4n1 - 4*c2*nu1*sin4n1 - 16*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n2 + 2*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n2)\
        + 2592*b1*cos4n1*sinhm1 - 5184*b2*mu0*cos4n1*sinhm1 + 5184*b2*mu1*cos4n1*sinhm1 - 2592*bc*nu0*cos4n1*sinhm1 + 
       2592*bc*nu1*cos4n1*sinhm1 - xeb1*cos4n1*sinhm1 + tnb2*mu0*cos4n1*sinhm1 - tnb2*mu1*cos4n1*sinhm1 + 
       xebc*nu0*cos4n1*sinhm1 - xebc*nu1*cos4n1*sinhm1 - 2592*b1*cos4n2*sinhm1 + 5184*b2*mu0*cos4n2*sinhm1 - 
       5184*b2*mu1*cos4n2*sinhm1 + 2592*bc*nu0*cos4n2*sinhm1 - 2592*bc*nu1*cos4n2*sinhm1 + xeb1*cos4n2*sinhm1 - 
       tnb2*mu0*cos4n2*sinhm1 + tnb2*mu1*cos4n2*sinhm1 - xebc*nu0*cos4n2*sinhm1 + xebc*nu1*cos4n2*sinhm1 - 
       tnbc*sin4n1*sinhm1 + 162*bc*sin4n1*sinhm1 + tnbc*sin4n2*sinhm1 - 162*bc*sin4n2*sinhm1 - tnb1*cosn1*sinh2m1 + 
       2592*b2*mu0*cosn1*sinh2m1 - 2592*b2*mu1*cosn1*sinh2m1 + tnbc*nu0*cosn1*sinh2m1 - tnbc*nu1*cosn1*sinh2m1 + 
       ofb1*cos3n1*sinh2m1 - eeb2*mu0*cos3n1*sinh2m1 + eeb2*mu1*cos3n1*sinh2m1 - ofbc*nu0*cos3n1*sinh2m1 + 
       ofbc*nu1*cos3n1*sinh2m1 + tnb1*cosn2*sinh2m1 - 2592*b2*mu0*cosn2*sinh2m1 + 2592*b2*mu1*cosn2*sinh2m1 - 
       tnbc*nu0*cosn2*sinh2m1 + tnbc*nu1*cosn2*sinh2m1 - ofb1*cos3n2*sinh2m1 + eeb2*mu0*cos3n2*sinh2m1 - 
       eeb2*mu1*cos3n2*sinh2m1 + ofbc*nu0*cos3n2*sinh2m1 - ofbc*nu1*cos3n2*sinh2m1 + tnbc*sinn1*sinh2m1 - 
       48*bc*sin3n1*sinh2m1 - tnbc*sinn2*sinh2m1 + 48*bc*sin3n2*sinh2m1 - 96*b1*cos4n1*sinh3m1 + 192*b2*mu0*cos4n1*sinh3m1 - 
       192*b2*mu1*cos4n1*sinh3m1 + 96*bc*nu0*cos4n1*sinh3m1 - 96*bc*nu1*cos4n1*sinh3m1 + 24*b1*cos4n1*sinh3m1 - 
       48*b2*mu0*cos4n1*sinh3m1 + 48*b2*mu1*cos4n1*sinh3m1 - 24*bc*nu0*cos4n1*sinh3m1 + 24*bc*nu1*cos4n1*sinh3m1 + 
       96*b1*cos4n2*sinh3m1 - 192*b2*mu0*cos4n2*sinh3m1 + 192*b2*mu1*cos4n2*sinh3m1 - 96*bc*nu0*cos4n2*sinh3m1 + 
       96*bc*nu1*cos4n2*sinh3m1 - 24*b1*cos4n2*sinh3m1 + 48*b2*mu0*cos4n2*sinh3m1 - 48*b2*mu1*cos4n2*sinh3m1 + 
       24*bc*nu0*cos4n2*sinh3m1 - 24*bc*nu1*cos4n2*sinh3m1 + 48*bc*sin4n1*sinh3m1 - 6*bc*sin4n1*sinh3m1 - 48*bc*sin4n2*sinh3m1 + 
       6*bc*sin4n2*sinh3m1 + 162*b1*cosn1*sinh4m1 - 324*b2*mu0*cosn1*sinh4m1 + 324*b2*mu1*cosn1*sinh4m1 - 162*bc*nu0*cosn1*sinh4m1 + 
       162*bc*nu1*cosn1*sinh4m1 - 18*b1*cos3n1*sinh4m1 + 36*b2*mu0*cos3n1*sinh4m1 - 36*b2*mu1*cos3n1*sinh4m1 + 
       18*bc*nu0*cos3n1*sinh4m1 - 18*bc*nu1*cos3n1*sinh4m1 - 162*b1*cosn2*sinh4m1 + 324*b2*mu0*cosn2*sinh4m1 - 
       324*b2*mu1*cosn2*sinh4m1 + 162*bc*nu0*cosn2*sinh4m1 - 162*bc*nu1*cosn2*sinh4m1 + 18*b1*cos3n2*sinh4m1 - 
       36*b2*mu0*cos3n2*sinh4m1 + 36*b2*mu1*cos3n2*sinh4m1 - 18*bc*nu0*cos3n2*sinh4m1 + 18*bc*nu1*cos3n2*sinh4m1 - 
       162*bc*sinn1*sinh4m1 + 6*bc*sin3n1*sinh4m1 + 162*bc*sinn2*sinh4m1 - 6*bc*sin3n2*sinh4m1 - 2592*b1*cos4n1*sinhm2 + 
       5184*b2*mu0*cos4n1*sinhm2 - 5184*b2*mu2*cos4n1*sinhm2 + 2592*bc*nu0*cos4n1*sinhm2 - 2592*bc*nu1*cos4n1*sinhm2 + 
       xeb1*cos4n1*sinhm2 - tnb2*mu0*cos4n1*sinhm2 + tnb2*mu2*cos4n1*sinhm2 - xebc*nu0*cos4n1*sinhm2 + 
       xebc*nu1*cos4n1*sinhm2 + 2592*b1*cos4n2*sinhm2 - 5184*b2*mu0*cos4n2*sinhm2 + 5184*b2*mu2*cos4n2*sinhm2 - 
       2592*bc*nu0*cos4n2*sinhm2 + 2592*bc*nu1*cos4n2*sinhm2 - xeb1*cos4n2*sinhm2 + tnb2*mu0*cos4n2*sinhm2 - 
       tnb2*mu2*cos4n2*sinhm2 + xebc*nu0*cos4n2*sinhm2 - xebc*nu1*cos4n2*sinhm2 + tnbc*sin4n1*sinhm2 - 162*bc*sin4n1*sinhm2 - 
       tnbc*sin4n2*sinhm2 + 162*bc*sin4n2*sinhm2 + tnb1*cosn1*sinh2m2 - 2592*b2*mu0*cosn1*sinh2m2 + 2592*b2*mu2*cosn1*sinh2m2 - 
       tnbc*nu0*cosn1*sinh2m2 + tnbc*nu1*cosn1*sinh2m2 - ofb1*cos3n1*sinh2m2 + eeb2*mu0*cos3n1*sinh2m2 - 
       eeb2*mu2*cos3n1*sinh2m2 + ofbc*nu0*cos3n1*sinh2m2 - ofbc*nu1*cos3n1*sinh2m2 - tnb1*cosn2*sinh2m2 + 
       2592*b2*mu0*cosn2*sinh2m2 - 2592*b2*mu2*cosn2*sinh2m2 + tnbc*nu0*cosn2*sinh2m2 - tnbc*nu1*cosn2*sinh2m2 + 
       ofb1*cos3n2*sinh2m2 - eeb2*mu0*cos3n2*sinh2m2 + eeb2*mu2*cos3n2*sinh2m2 - ofbc*nu0*cos3n2*sinh2m2 + 
       ofbc*nu1*cos3n2*sinh2m2 - tnbc*sinn1*sinh2m2 + 48*bc*sin3n1*sinh2m2 + tnbc*sinn2*sinh2m2 - 48*bc*sin3n2*sinh2m2 + 
       96*b1*cos4n1*sinh3m2 - 192*b2*mu0*cos4n1*sinh3m2 + 192*b2*mu2*cos4n1*sinh3m2 - 96*bc*nu0*cos4n1*sinh3m2 + 
       96*bc*nu1*cos4n1*sinh3m2 - 24*b1*cos4n1*sinh3m2 + 48*b2*mu0*cos4n1*sinh3m2 - 48*b2*mu2*cos4n1*sinh3m2 + 
       24*bc*nu0*cos4n1*sinh3m2 - 24*bc*nu1*cos4n1*sinh3m2 - 96*b1*cos4n2*sinh3m2 + 192*b2*mu0*cos4n2*sinh3m2 - 
       192*b2*mu2*cos4n2*sinh3m2 + 96*bc*nu0*cos4n2*sinh3m2 - 96*bc*nu1*cos4n2*sinh3m2 + 24*b1*cos4n2*sinh3m2 - 
       48*b2*mu0*cos4n2*sinh3m2 + 48*b2*mu2*cos4n2*sinh3m2 - 24*bc*nu0*cos4n2*sinh3m2 + 24*bc*nu1*cos4n2*sinh3m2 - 
       48*bc*sin4n1*sinh3m2 + 6*bc*sin4n1*sinh3m2 + 48*bc*sin4n2*sinh3m2 - 6*bc*sin4n2*sinh3m2 + 
       6*(-27*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cosn1 + 3*b1*cos3n1 - 6*b2*mu0*cos3n1 + 6*b2*mu2*cos3n1 - 3*bc*nu0*cos3n1 + 
          3*bc*nu1*cos3n1 + 27*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cosn2 - 3*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cos3n2 + 27*bc*sinn1 - 
          bc*sin3n1 - 27*bc*sinn2 + bc*sin3n2)*sinh4m2);

  double val2 = fp*(-288*a0*cos4n1*cosh3m1 - sfb2*cos4n1*cosh3m1 + ofc2*cos4n1*cosh3m1 + eeb1*mu0*cos4n1*cosh3m1 - 
       eeb2*pm02*cos4n1*cosh3m1 - eeb1*mu1*cos4n1*cosh3m1 + fnb2*mu0*mu1*cos4n1*cosh3m1 - eeb2*pm12*cos4n1*cosh3m1 + 
       eec1*nu0*cos4n1*cosh3m1 - eebc*mu0*nu0*cos4n1*cosh3m1 + eebc*mu1*nu0*cos4n1*cosh3m1 - eec2*pn02*cos4n1*cosh3m1 - 
       eec1*nu1*cos4n1*cosh3m1 + eebc*mu0*nu1*cos4n1*cosh3m1 - eebc*mu1*nu1*cos4n1*cosh3m1 + fnc2*nu0*nu1*cos4n1*cosh3m1 - 
       eec2*pn12*cos4n1*cosh3m1 + 72*a0*cos4n1*cosh3m1 + sxb2*cos4n1*cosh3m1 - nnc2*cos4n1*cosh3m1 - 
       stb1*mu0*cos4n1*cosh3m1 + stb2*pm02*cos4n1*cosh3m1 + stb1*mu1*cos4n1*cosh3m1 - ofb2*mu0*mu1*cos4n1*cosh3m1 + 
       stb2*pm12*cos4n1*cosh3m1 - stc1*nu0*cos4n1*cosh3m1 + stbc*mu0*nu0*cos4n1*cosh3m1 - stbc*mu1*nu0*cos4n1*cosh3m1 + 
       stc2*pn02*cos4n1*cosh3m1 + stc1*nu1*cos4n1*cosh3m1 - stbc*mu0*nu1*cos4n1*cosh3m1 + stbc*mu1*nu1*cos4n1*cosh3m1 - 
       ofc2*nu0*nu1*cos4n1*cosh3m1 + stc2*pn12*cos4n1*cosh3m1 + 288*a0*cos4n2*cosh3m1 + sfb2*cos4n2*cosh3m1 - 
       ofc2*cos4n2*cosh3m1 - eeb1*mu0*cos4n2*cosh3m1 + eeb2*pm02*cos4n2*cosh3m1 + eeb1*mu1*cos4n2*cosh3m1 - 
       fnb2*mu0*mu1*cos4n2*cosh3m1 + eeb2*pm12*cos4n2*cosh3m1 - eec1*nu0*cos4n2*cosh3m1 + eebc*mu0*nu0*cos4n2*cosh3m1 - 
       eebc*mu1*nu0*cos4n2*cosh3m1 + eec2*pn02*cos4n2*cosh3m1 + eec1*nu1*cos4n2*cosh3m1 - eebc*mu0*nu1*cos4n2*cosh3m1 + 
       eebc*mu1*nu1*cos4n2*cosh3m1 - fnc2*nu0*nu1*cos4n2*cosh3m1 + eec2*pn22*cos4n2*cosh3m1 - 72*a0*cos4n2*cosh3m1 - 
       sxb2*cos4n2*cosh3m1 + nnc2*cos4n2*cosh3m1 + stb1*mu0*cos4n2*cosh3m1 - stb2*pm02*cos4n2*cosh3m1 - 
       stb1*mu1*cos4n2*cosh3m1 + ofb2*mu0*mu1*cos4n2*cosh3m1 - stb2*pm12*cos4n2*cosh3m1 + stc1*nu0*cos4n2*cosh3m1 - 
       stbc*mu0*nu0*cos4n2*cosh3m1 + stbc*mu1*nu0*cos4n2*cosh3m1 - stc2*pn02*cos4n2*cosh3m1 - stc1*nu1*cos4n2*cosh3m1 + 
       stbc*mu0*nu1*cos4n2*cosh3m1 - stbc*mu1*nu1*cos4n2*cosh3m1 + ofc2*nu0*nu1*cos4n2*cosh3m1 - stc2*pn22*cos4n2*cosh3m1 - 
       648*a0*cosn1*cosh4m1 - 81*b2*cosn1*cosh4m1 + tnc2*cosn1*cosh4m1 + xeb1*mu0*cosn1*cosh4m1 - 
       xeb2*pm02*cosn1*cosh4m1 - xeb1*mu1*cosn1*cosh4m1 + tnb2*mu0*mu1*cosn1*cosh4m1 - xeb2*pm12*cosn1*cosh4m1 + 
       xec1*nu0*cosn1*cosh4m1 - xebc*mu0*nu0*cosn1*cosh4m1 + xebc*mu1*nu0*cosn1*cosh4m1 - xec2*pn02*cosn1*cosh4m1 - 
       xec1*nu1*cosn1*cosh4m1 + xebc*mu0*nu1*cosn1*cosh4m1 - xebc*mu1*nu1*cosn1*cosh4m1 + tnc2*nu0*nu1*cosn1*cosh4m1 - 
       xec2*pn12*cosn1*cosh4m1 + 72*a0*cos3n1*cosh4m1 + nnb2*cos3n1*cosh4m1 - sxc2*cos3n1*cosh4m1 - 
       stb1*mu0*cos3n1*cosh4m1 + stb2*pm02*cos3n1*cosh4m1 + stb1*mu1*cos3n1*cosh4m1 - ofb2*mu0*mu1*cos3n1*cosh4m1 + 
       stb2*pm12*cos3n1*cosh4m1 - stc1*nu0*cos3n1*cosh4m1 + stbc*mu0*nu0*cos3n1*cosh4m1 - stbc*mu1*nu0*cos3n1*cosh4m1 + 
       stc2*pn02*cos3n1*cosh4m1 + stc1*nu1*cos3n1*cosh4m1 - stbc*mu0*nu1*cos3n1*cosh4m1 + stbc*mu1*nu1*cos3n1*cosh4m1 - 
       ofc2*nu0*nu1*cos3n1*cosh4m1 + stc2*pn12*cos3n1*cosh4m1 + 648*a0*cosn2*cosh4m1 + 81*b2*cosn2*cosh4m1 - 
       tnc2*cosn2*cosh4m1 - xeb1*mu0*cosn2*cosh4m1 + xeb2*pm02*cosn2*cosh4m1 + xeb1*mu1*cosn2*cosh4m1 - 
       tnb2*mu0*mu1*cosn2*cosh4m1 + xeb2*pm12*cosn2*cosh4m1 - xec1*nu0*cosn2*cosh4m1 + xebc*mu0*nu0*cosn2*cosh4m1 - 
       xebc*mu1*nu0*cosn2*cosh4m1 + xec2*pn02*cosn2*cosh4m1 + xec1*nu1*cosn2*cosh4m1 - xebc*mu0*nu1*cosn2*cosh4m1 + 
       xebc*mu1*nu1*cosn2*cosh4m1 - tnc2*nu0*nu1*cosn2*cosh4m1 + xec2*pn22*cosn2*cosh4m1 - 72*a0*cos3n2*cosh4m1 - 
       nnb2*cos3n2*cosh4m1 + sxc2*cos3n2*cosh4m1 + stb1*mu0*cos3n2*cosh4m1 - stb2*pm02*cos3n2*cosh4m1 - 
       stb1*mu1*cos3n2*cosh4m1 + ofb2*mu0*mu1*cos3n2*cosh4m1 - stb2*pm12*cos3n2*cosh4m1 + stc1*nu0*cos3n2*cosh4m1 - 
       stbc*mu0*nu0*cos3n2*cosh4m1 + stbc*mu1*nu0*cos3n2*cosh4m1 - stc2*pn02*cos3n2*cosh4m1 - stc1*nu1*cos3n2*cosh4m1 + 
       stbc*mu0*nu1*cos3n2*cosh4m1 - stbc*mu1*nu1*cos3n2*cosh4m1 + ofc2*nu0*nu1*cos3n2*cosh4m1 - stc2*pn22*cos3n2*cosh4m1 - 
       2592*a0*cos4n1*coshm2 - 5184*b2*cos4n1*coshm2 + tnc2*cos4n1*coshm2 + 2592*b1*mu0*cos4n1*coshm2 - 
       2592*b2*pm02*cos4n1*coshm2 - 2592*b1*mu2*cos4n1*coshm2 + 5184*b2*mu0*mu2*cos4n1*coshm2 - 2592*b2*pm22*cos4n1*coshm2 + 
       2592*c1*nu0*cos4n1*coshm2 - 2592*bc*mu0*nu0*cos4n1*coshm2 + 2592*bc*mu2*nu0*cos4n1*coshm2 - 2592*c2*pn02*cos4n1*coshm2 - 
       2592*c1*nu1*cos4n1*coshm2 + 2592*bc*mu0*nu1*cos4n1*coshm2 - 2592*bc*mu2*nu1*cos4n1*coshm2 + 5184*c2*nu0*nu1*cos4n1*coshm2 - 
       2592*c2*pn12*cos4n1*coshm2 + 648*a0*cos4n1*coshm2 + tnb2*cos4n1*coshm2 - 81*c2*cos4n1*coshm2 - 
       xeb1*mu0*cos4n1*coshm2 + xeb2*pm02*cos4n1*coshm2 + xeb1*mu2*cos4n1*coshm2 - tnb2*mu0*mu2*cos4n1*coshm2 + 
       xeb2*pm22*cos4n1*coshm2 - xec1*nu0*cos4n1*coshm2 + xebc*mu0*nu0*cos4n1*coshm2 - xebc*mu2*nu0*cos4n1*coshm2 + 
       xec2*pn02*cos4n1*coshm2 + xec1*nu1*cos4n1*coshm2 - xebc*mu0*nu1*cos4n1*coshm2 + xebc*mu2*nu1*cos4n1*coshm2 - 
       tnc2*nu0*nu1*cos4n1*coshm2 + xec2*pn12*cos4n1*coshm2 + 2592*a0*cos4n2*coshm2 + 5184*b2*cos4n2*coshm2 - 
       tnc2*cos4n2*coshm2 - 2592*b1*mu0*cos4n2*coshm2 + 2592*b2*pm02*cos4n2*coshm2 + 2592*b1*mu2*cos4n2*coshm2 - 
       5184*b2*mu0*mu2*cos4n2*coshm2 + 2592*b2*pm22*cos4n2*coshm2 - 2592*c1*nu0*cos4n2*coshm2 + 2592*bc*mu0*nu0*cos4n2*coshm2 - 
       2592*bc*mu2*nu0*cos4n2*coshm2 + 2592*c2*pn02*cos4n2*coshm2 + 2592*c1*nu1*cos4n2*coshm2 - 2592*bc*mu0*nu1*cos4n2*coshm2 + 
       2592*bc*mu2*nu1*cos4n2*coshm2 - 5184*c2*nu0*nu1*cos4n2*coshm2 + 2592*c2*pn22*cos4n2*coshm2 - 648*a0*cos4n2*coshm2 - 
       tnb2*cos4n2*coshm2 + 81*c2*cos4n2*coshm2 + xeb1*mu0*cos4n2*coshm2 - xeb2*pm02*cos4n2*coshm2 - 
       xeb1*mu2*cos4n2*coshm2 + tnb2*mu0*mu2*cos4n2*coshm2 - xeb2*pm22*cos4n2*coshm2 + xec1*nu0*cos4n2*coshm2 - 
       xebc*mu0*nu0*cos4n2*coshm2 + xebc*mu2*nu0*cos4n2*coshm2 - xec2*pn02*cos4n2*coshm2 - xec1*nu1*cos4n2*coshm2 + 
       xebc*mu0*nu1*cos4n2*coshm2 - xebc*mu2*nu1*cos4n2*coshm2 + tnc2*nu0*nu1*cos4n2*coshm2 - xec2*pn22*cos4n2*coshm2 - 
       2592*a0*cosn1*cosh2m2 - tnb2*cosn1*cosh2m2 + 5184*c2*cosn1*cosh2m2 + 2592*b1*mu0*cosn1*cosh2m2 - 
       2592*b2*pm02*cosn1*cosh2m2 - 2592*b1*mu2*cosn1*cosh2m2 + 5184*b2*mu0*mu2*cosn1*cosh2m2 - 2592*b2*pm22*cosn1*cosh2m2 + 
       2592*c1*nu0*cosn1*cosh2m2 - 2592*bc*mu0*nu0*cosn1*cosh2m2 + 2592*bc*mu2*nu0*cosn1*cosh2m2 - 2592*c2*pn02*cosn1*cosh2m2 - 
       2592*c1*nu1*cosn1*cosh2m2 + 2592*bc*mu0*nu1*cosn1*cosh2m2 - 2592*bc*mu2*nu1*cosn1*cosh2m2 + 5184*c2*nu0*nu1*cosn1*cosh2m2 - 
       2592*c2*pn12*cosn1*cosh2m2 + 288*a0*cos3n1*cosh2m2 + ofb2*cos3n1*cosh2m2 - sfc2*cos3n1*cosh2m2 - 
       eeb1*mu0*cos3n1*cosh2m2 + eeb2*pm02*cos3n1*cosh2m2 + eeb1*mu2*cos3n1*cosh2m2 - fnb2*mu0*mu2*cos3n1*cosh2m2 + 
       eeb2*pm22*cos3n1*cosh2m2 - eec1*nu0*cos3n1*cosh2m2 + eebc*mu0*nu0*cos3n1*cosh2m2 - eebc*mu2*nu0*cos3n1*cosh2m2 + 
       eec2*pn02*cos3n1*cosh2m2 + eec1*nu1*cos3n1*cosh2m2 - eebc*mu0*nu1*cos3n1*cosh2m2 + eebc*mu2*nu1*cos3n1*cosh2m2 - 
       fnc2*nu0*nu1*cos3n1*cosh2m2 + eec2*pn12*cos3n1*cosh2m2 + 2592*a0*cosn2*cosh2m2 + tnb2*cosn2*cosh2m2 - 
       5184*c2*cosn2*cosh2m2 - 2592*b1*mu0*cosn2*cosh2m2 + 2592*b2*pm02*cosn2*cosh2m2 + 2592*b1*mu2*cosn2*cosh2m2 - 
       5184*b2*mu0*mu2*cosn2*cosh2m2 + 2592*b2*pm22*cosn2*cosh2m2 - 2592*c1*nu0*cosn2*cosh2m2 + 2592*bc*mu0*nu0*cosn2*cosh2m2 - 
       2592*bc*mu2*nu0*cosn2*cosh2m2 + 2592*c2*pn02*cosn2*cosh2m2 + 2592*c1*nu1*cosn2*cosh2m2 - 2592*bc*mu0*nu1*cosn2*cosh2m2 + 
       2592*bc*mu2*nu1*cosn2*cosh2m2 - 5184*c2*nu0*nu1*cosn2*cosh2m2 + 2592*c2*pn22*cosn2*cosh2m2 - 288*a0*cos3n2*cosh2m2 - 
       ofb2*cos3n2*cosh2m2 + sfc2*cos3n2*cosh2m2 + eeb1*mu0*cos3n2*cosh2m2 - eeb2*pm02*cos3n2*cosh2m2 - 
       eeb1*mu2*cos3n2*cosh2m2 + fnb2*mu0*mu2*cos3n2*cosh2m2 - eeb2*pm22*cos3n2*cosh2m2 + eec1*nu0*cos3n2*cosh2m2 - 
       eebc*mu0*nu0*cos3n2*cosh2m2 + eebc*mu2*nu0*cos3n2*cosh2m2 - eec2*pn02*cos3n2*cosh2m2 - eec1*nu1*cos3n2*cosh2m2 + 
       eebc*mu0*nu1*cos3n2*cosh2m2 - eebc*mu2*nu1*cos3n2*cosh2m2 + fnc2*nu0*nu1*cos3n2*cosh2m2 - eec2*pn22*cos3n2*cosh2m2 + 
       288*a0*cos4n1*cosh3m2 + sfb2*cos4n1*cosh3m2 - ofc2*cos4n1*cosh3m2 - eeb1*mu0*cos4n1*cosh3m2 + 
       eeb2*pm02*cos4n1*cosh3m2 + eeb1*mu2*cos4n1*cosh3m2 - fnb2*mu0*mu2*cos4n1*cosh3m2 + eeb2*pm22*cos4n1*cosh3m2 - 
       eec1*nu0*cos4n1*cosh3m2 + eebc*mu0*nu0*cos4n1*cosh3m2 - eebc*mu2*nu0*cos4n1*cosh3m2 + eec2*pn02*cos4n1*cosh3m2 + 
       eec1*nu1*cos4n1*cosh3m2 - eebc*mu0*nu1*cos4n1*cosh3m2 + eebc*mu2*nu1*cos4n1*cosh3m2 - fnc2*nu0*nu1*cos4n1*cosh3m2 + 
       eec2*pn12*cos4n1*cosh3m2 - 72*a0*cos4n1*cosh3m2 - sxb2*cos4n1*cosh3m2 + nnc2*cos4n1*cosh3m2 + 
       stb1*mu0*cos4n1*cosh3m2 - stb2*pm02*cos4n1*cosh3m2 - stb1*mu2*cos4n1*cosh3m2 + ofb2*mu0*mu2*cos4n1*cosh3m2 - 
       stb2*pm22*cos4n1*cosh3m2 + stc1*nu0*cos4n1*cosh3m2 - stbc*mu0*nu0*cos4n1*cosh3m2 + stbc*mu2*nu0*cos4n1*cosh3m2 - 
       stc2*pn02*cos4n1*cosh3m2 - stc1*nu1*cos4n1*cosh3m2 + stbc*mu0*nu1*cos4n1*cosh3m2 - stbc*mu2*nu1*cos4n1*cosh3m2 + 
       ofc2*nu0*nu1*cos4n1*cosh3m2 - stc2*pn12*cos4n1*cosh3m2 - 288*a0*cos4n2*cosh3m2 - sfb2*cos4n2*cosh3m2 + 
       ofc2*cos4n2*cosh3m2 + eeb1*mu0*cos4n2*cosh3m2 - eeb2*pm02*cos4n2*cosh3m2 - eeb1*mu2*cos4n2*cosh3m2 + 
       fnb2*mu0*mu2*cos4n2*cosh3m2 - eeb2*pm22*cos4n2*cosh3m2 + eec1*nu0*cos4n2*cosh3m2 - eebc*mu0*nu0*cos4n2*cosh3m2 + 
       eebc*mu2*nu0*cos4n2*cosh3m2 - eec2*pn02*cos4n2*cosh3m2 - eec1*nu1*cos4n2*cosh3m2 + eebc*mu0*nu1*cos4n2*cosh3m2 - 
       eebc*mu2*nu1*cos4n2*cosh3m2 + fnc2*nu0*nu1*cos4n2*cosh3m2 - eec2*pn22*cos4n2*cosh3m2 + 72*a0*cos4n2*cosh3m2 + 
       sxb2*cos4n2*cosh3m2 - nnc2*cos4n2*cosh3m2 - stb1*mu0*cos4n2*cosh3m2 + stb2*pm02*cos4n2*cosh3m2 + 
       stb1*mu2*cos4n2*cosh3m2 - ofb2*mu0*mu2*cos4n2*cosh3m2 + stb2*pm22*cos4n2*cosh3m2 - stc1*nu0*cos4n2*cosh3m2 + 
       stbc*mu0*nu0*cos4n2*cosh3m2 - stbc*mu2*nu0*cos4n2*cosh3m2 + stc2*pn02*cos4n2*cosh3m2 + stc1*nu1*cos4n2*cosh3m2 - 
       stbc*mu0*nu1*cos4n2*cosh3m2 + stbc*mu2*nu1*cos4n2*cosh3m2 - ofc2*nu0*nu1*cos4n2*cosh3m2 + stc2*pn22*cos4n2*cosh3m2 + 
       648*a0*cosn1*cosh4m2 + 81*b2*cosn1*cosh4m2 - tnc2*cosn1*cosh4m2 - xeb1*mu0*cosn1*cosh4m2 + 
       xeb2*pm02*cosn1*cosh4m2 + xeb1*mu2*cosn1*cosh4m2 - tnb2*mu0*mu2*cosn1*cosh4m2 + xeb2*pm22*cosn1*cosh4m2 - 
       xec1*nu0*cosn1*cosh4m2 + xebc*mu0*nu0*cosn1*cosh4m2 - xebc*mu2*nu0*cosn1*cosh4m2 + xec2*pn02*cosn1*cosh4m2 + 
       xec1*nu1*cosn1*cosh4m2 - xebc*mu0*nu1*cosn1*cosh4m2 + xebc*mu2*nu1*cosn1*cosh4m2 - tnc2*nu0*nu1*cosn1*cosh4m2 + 
       xec2*pn12*cosn1*cosh4m2 - 72*a0*cos3n1*cosh4m2 - nnb2*cos3n1*cosh4m2 + sxc2*cos3n1*cosh4m2 + 
       stb1*mu0*cos3n1*cosh4m2 - stb2*pm02*cos3n1*cosh4m2 - stb1*mu2*cos3n1*cosh4m2 + ofb2*mu0*mu2*cos3n1*cosh4m2 - 
       stb2*pm22*cos3n1*cosh4m2 + stc1*nu0*cos3n1*cosh4m2 - stbc*mu0*nu0*cos3n1*cosh4m2 + stbc*mu2*nu0*cos3n1*cosh4m2 - 
       stc2*pn02*cos3n1*cosh4m2 - stc1*nu1*cos3n1*cosh4m2 + stbc*mu0*nu1*cos3n1*cosh4m2 - stbc*mu2*nu1*cos3n1*cosh4m2 + 
       ofc2*nu0*nu1*cos3n1*cosh4m2 - stc2*pn12*cos3n1*cosh4m2 - 648*a0*cosn2*cosh4m2 - 81*b2*cosn2*cosh4m2 + 
       tnc2*cosn2*cosh4m2 + xeb1*mu0*cosn2*cosh4m2 - xeb2*pm02*cosn2*cosh4m2 - xeb1*mu2*cosn2*cosh4m2 + 
       tnb2*mu0*mu2*cosn2*cosh4m2 - xeb2*pm22*cosn2*cosh4m2 + xec1*nu0*cosn2*cosh4m2 - xebc*mu0*nu0*cosn2*cosh4m2 + 
       xebc*mu2*nu0*cosn2*cosh4m2 - xec2*pn02*cosn2*cosh4m2 - xec1*nu1*cosn2*cosh4m2 + xebc*mu0*nu1*cosn2*cosh4m2 - 
       xebc*mu2*nu1*cosn2*cosh4m2 + tnc2*nu0*nu1*cosn2*cosh4m2 - xec2*pn22*cosn2*cosh4m2 + 72*a0*cos3n2*cosh4m2 + 
       nnb2*cos3n2*cosh4m2 - sxc2*cos3n2*cosh4m2 - stb1*mu0*cos3n2*cosh4m2 + stb2*pm02*cos3n2*cosh4m2 + 
       stb1*mu2*cos3n2*cosh4m2 - ofb2*mu0*mu2*cos3n2*cosh4m2 + stb2*pm22*cos3n2*cosh4m2 - stc1*nu0*cos3n2*cosh4m2 + 
       stbc*mu0*nu0*cos3n2*cosh4m2 - stbc*mu2*nu0*cos3n2*cosh4m2 + stc2*pn02*cos3n2*cosh4m2 + stc1*nu1*cos3n2*cosh4m2 - 
       stbc*mu0*nu1*cos3n2*cosh4m2 + stbc*mu2*nu1*cos3n2*cosh4m2 - ofc2*nu0*nu1*cos3n2*cosh4m2 + stc2*pn22*cos3n2*cosh4m2 + 
       xec1*cosh4m1*sinn1 - xebc*mu0*cosh4m1*sinn1 + xebc*mu1*cosh4m1*sinn1 - tnc2*nu0*cosh4m1*sinn1 + 
       tnc2*nu1*cosh4m1*sinn1 + 2592*c1*cosh2m2*sinn1 - 2592*bc*mu0*cosh2m2*sinn1 + 2592*bc*mu2*cosh2m2*sinn1 - 
       5184*c2*nu0*cosh2m2*sinn1 + 5184*c2*nu1*cosh2m2*sinn1 - xec1*cosh4m2*sinn1 + xebc*mu0*cosh4m2*sinn1 - 
       xebc*mu2*cosh4m2*sinn1 + tnc2*nu0*cosh4m2*sinn1 - tnc2*nu1*cosh4m2*sinn1 + ofc1*cosh3m1*sin4n1 - 
       ofbc*mu0*cosh3m1*sin4n1 + ofbc*mu1*cosh3m1*sin4n1 - eec2*nu0*cosh3m1*sin4n1 + eec2*nu1*cosh3m1*sin4n1 + 
       tnc1*coshm2*sin4n1 - tnbc*mu0*coshm2*sin4n1 + tnbc*mu2*coshm2*sin4n1 - 2592*c2*nu0*coshm2*sin4n1 + 
       2592*c2*nu1*coshm2*sin4n1 - ofc1*cosh3m2*sin4n1 + ofbc*mu0*cosh3m2*sin4n1 - ofbc*mu2*cosh3m2*sin4n1 + 
       eec2*nu0*cosh3m2*sin4n1 - eec2*nu1*cosh3m2*sin4n1 - 24*c1*cosh4m1*sin3n1 + 24*bc*mu0*cosh4m1*sin3n1 - 
       24*bc*mu1*cosh4m1*sin3n1 + 48*c2*nu0*cosh4m1*sin3n1 - 48*c2*nu1*cosh4m1*sin3n1 - 96*c1*cosh2m2*sin3n1 + 
       96*bc*mu0*cosh2m2*sin3n1 - 96*bc*mu2*cosh2m2*sin3n1 + 192*c2*nu0*cosh2m2*sin3n1 - 192*c2*nu1*cosh2m2*sin3n1 + 
       24*c1*cosh4m2*sin3n1 - 24*bc*mu0*cosh4m2*sin3n1 + 24*bc*mu2*cosh4m2*sin3n1 - 48*c2*nu0*cosh4m2*sin3n1 + 
       48*c2*nu1*cosh4m2*sin3n1 - 18*c1*cosh3m1*sin4n1 + 18*bc*mu0*cosh3m1*sin4n1 - 18*bc*mu1*cosh3m1*sin4n1 + 
       36*c2*nu0*cosh3m1*sin4n1 - 36*c2*nu1*cosh3m1*sin4n1 - 162*c1*coshm2*sin4n1 + 162*bc*mu0*coshm2*sin4n1 - 
       162*bc*mu2*coshm2*sin4n1 + 324*c2*nu0*coshm2*sin4n1 - 324*c2*nu1*coshm2*sin4n1 + 18*c1*cosh3m2*sin4n1 - 
       18*bc*mu0*cosh3m2*sin4n1 + 18*bc*mu2*cosh3m2*sin4n1 - 36*c2*nu0*cosh3m2*sin4n1 + 36*c2*nu1*cosh3m2*sin4n1 - 
       xec1*cosh4m1*sinn2 + xebc*mu0*cosh4m1*sinn2 - xebc*mu1*cosh4m1*sinn2 + tnc2*nu0*cosh4m1*sinn2 - 
       tnc2*nu1*cosh4m1*sinn2 - 2592*c1*cosh2m2*sinn2 + 2592*bc*mu0*cosh2m2*sinn2 - 2592*bc*mu2*cosh2m2*sinn2 + 
       5184*c2*nu0*cosh2m2*sinn2 - 5184*c2*nu1*cosh2m2*sinn2 + xec1*cosh4m2*sinn2 - xebc*mu0*cosh4m2*sinn2 + 
       xebc*mu2*cosh4m2*sinn2 - tnc2*nu0*cosh4m2*sinn2 + tnc2*nu1*cosh4m2*sinn2 - ofc1*cosh3m1*sin4n2 + 
       ofbc*mu0*cosh3m1*sin4n2 - ofbc*mu1*cosh3m1*sin4n2 + eec2*nu0*cosh3m1*sin4n2 - eec2*nu1*cosh3m1*sin4n2 - 
       tnc1*coshm2*sin4n2 + tnbc*mu0*coshm2*sin4n2 - tnbc*mu2*coshm2*sin4n2 + 2592*c2*nu0*coshm2*sin4n2 - 
       2592*c2*nu1*coshm2*sin4n2 + ofc1*cosh3m2*sin4n2 - ofbc*mu0*cosh3m2*sin4n2 + ofbc*mu2*cosh3m2*sin4n2 - 
       eec2*nu0*cosh3m2*sin4n2 + eec2*nu1*cosh3m2*sin4n2 + 24*c1*cosh4m1*sin3n2 - 24*bc*mu0*cosh4m1*sin3n2 + 
       24*bc*mu1*cosh4m1*sin3n2 - 48*c2*nu0*cosh4m1*sin3n2 + 48*c2*nu1*cosh4m1*sin3n2 + 96*c1*cosh2m2*sin3n2 - 
       96*bc*mu0*cosh2m2*sin3n2 + 96*bc*mu2*cosh2m2*sin3n2 - 192*c2*nu0*cosh2m2*sin3n2 + 192*c2*nu1*cosh2m2*sin3n2 - 
       24*c1*cosh4m2*sin3n2 + 24*bc*mu0*cosh4m2*sin3n2 - 24*bc*mu2*cosh4m2*sin3n2 + 48*c2*nu0*cosh4m2*sin3n2 - 
       48*c2*nu1*cosh4m2*sin3n2 + 16*cosh2m1*(81*(2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn012) + 2*c1*(-nu0 + nu1) - 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 - 18*a0*cos3n1 - nnb2*cos3n1 + 4*c2*cos3n1 + 18*b1*mu0*cos3n1 - 18*b2*pm02*cos3n1 - 
          18*b1*mu1*cos3n1 + 36*b2*mu0*mu1*cos3n1 - 18*b2*pm12*cos3n1 + 18*c1*nu0*cos3n1 - 18*bc*mu0*nu0*cos3n1 + 18*bc*mu1*nu0*cos3n1 - 
          18*c2*pn02*cos3n1 - 18*c1*nu1*cos3n1 + 18*bc*mu0*nu1*cos3n1 - 18*bc*mu1*nu1*cos3n1 + 36*c2*nu0*nu1*cos3n1 - 
          18*c2*pn12*cos3n1 - 81*(2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn022) + 2*c1*(-nu0 + nu2) - 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 + (18*a0 + nnb2*(1 + 2*pm012) + 
             2*(nnb1*(-mu0 + mu1) + c2*(-2 + 9*pn022) + 9*(c1 + bc*(-mu0 + mu1))*(-nu0 + nu2)))*cos3n2 - 162*c1*sinn1 + 162*bc*mu0*sinn1 - 
          162*bc*mu1*sinn1 + 324*c2*nu0*sinn1 - 324*c2*nu1*sinn1 + 6*c1*sin3n1 - 6*bc*mu0*sin3n1 + 6*bc*mu1*sin3n1 - 12*c2*nu0*sin3n1 + 
          12*c2*nu1*sin3n1 + 162*(c1 + bc*(-mu0 + mu1) + 2*c2*(-nu0 + nu2))*sinn2 - 6*(c1 + bc*(-mu0 + mu1) + 2*c2*(-nu0 + nu2))*sin3n2) + 
       18*c1*cosh3m1*sin4n2 - 18*bc*mu0*cosh3m1*sin4n2 + 18*bc*mu1*cosh3m1*sin4n2 - 36*c2*nu0*cosh3m1*sin4n2 + 
       36*c2*nu1*cosh3m1*sin4n2 + 162*c1*coshm2*sin4n2 - 162*bc*mu0*coshm2*sin4n2 + 162*bc*mu2*coshm2*sin4n2 - 324*c2*nu0*coshm2*sin4n2 + 
       324*c2*nu1*coshm2*sin4n2 - 18*c1*cosh3m2*sin4n2 + 18*bc*mu0*cosh3m2*sin4n2 - 18*bc*mu2*cosh3m2*sin4n2 + 
       36*c2*nu0*cosh3m2*sin4n2 - 36*c2*nu1*cosh3m2*sin4n2 + 
       81*coshm1*(16*(2*a0 + 2*b2*(2 + pm012) + c2*(-1 + 2*pn012) + 2*c1*(-nu0 + nu1) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos4n1 - 
          8*a0*cos4n1 - sxb2*cos4n1 + c2*cos4n1 + 8*b1*mu0*cos4n1 - 8*b2*pm02*cos4n1 - 8*b1*mu1*cos4n1 + sxb2*mu0*mu1*cos4n1 - 
          8*b2*pm12*cos4n1 + 8*c1*nu0*cos4n1 - 8*bc*mu0*nu0*cos4n1 + 8*bc*mu1*nu0*cos4n1 - 8*c2*pn02*cos4n1 - 8*c1*nu1*cos4n1 + 
          8*bc*mu0*nu1*cos4n1 - 8*bc*mu1*nu1*cos4n1 + sxc2*nu0*nu1*cos4n1 - 8*c2*pn12*cos4n1 - 
          16*(2*a0 + 2*b2*(2 + pm012) + c2*(-1 + 2*pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2 + 
          (8*a0 + 8*b2*(2 + pm012) + c2*(-1 + 8*pn022) + 8*c1*(-nu0 + nu2) - 8*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2 - 
          sxc1*sin4n1 + sxbc*mu0*sin4n1 - sxbc*mu1*sin4n1 + 32*c2*nu0*sin4n1 - 32*c2*nu1*sin4n1 + 2*c1*sin4n1 - 2*bc*mu0*sin4n1 + 
          2*bc*mu1*sin4n1 - 4*c2*nu0*sin4n1 + 4*c2*nu1*sin4n1 + 16*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n2 - 
          2*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n2) - 2592*b1*cos4n1*sinhm1 + 5184*b2*mu0*cos4n1*sinhm1 - 5184*b2*mu1*cos4n1*sinhm1 + 
       2592*bc*nu0*cos4n1*sinhm1 - 2592*bc*nu1*cos4n1*sinhm1 + xeb1*cos4n1*sinhm1 - tnb2*mu0*cos4n1*sinhm1 + 
       tnb2*mu1*cos4n1*sinhm1 - xebc*nu0*cos4n1*sinhm1 + xebc*nu1*cos4n1*sinhm1 + 2592*b1*cos4n2*sinhm1 - 
       5184*b2*mu0*cos4n2*sinhm1 + 5184*b2*mu1*cos4n2*sinhm1 - 2592*bc*nu0*cos4n2*sinhm1 + 2592*bc*nu1*cos4n2*sinhm1 - 
       xeb1*cos4n2*sinhm1 + tnb2*mu0*cos4n2*sinhm1 - tnb2*mu1*cos4n2*sinhm1 + xebc*nu0*cos4n2*sinhm1 - 
       xebc*nu1*cos4n2*sinhm1 + tnbc*sin4n1*sinhm1 - 162*bc*sin4n1*sinhm1 - tnbc*sin4n2*sinhm1 + 162*bc*sin4n2*sinhm1 - 
       tnb1*cosn1*sinh2m1 + 2592*b2*mu0*cosn1*sinh2m1 - 2592*b2*mu1*cosn1*sinh2m1 + tnbc*nu0*cosn1*sinh2m1 - 
       tnbc*nu1*cosn1*sinh2m1 + ofb1*cos3n1*sinh2m1 - eeb2*mu0*cos3n1*sinh2m1 + eeb2*mu1*cos3n1*sinh2m1 - 
       ofbc*nu0*cos3n1*sinh2m1 + ofbc*nu1*cos3n1*sinh2m1 + tnb1*cosn2*sinh2m1 - 2592*b2*mu0*cosn2*sinh2m1 + 
       2592*b2*mu1*cosn2*sinh2m1 - tnbc*nu0*cosn2*sinh2m1 + tnbc*nu1*cosn2*sinh2m1 - ofb1*cos3n2*sinh2m1 + 
       eeb2*mu0*cos3n2*sinh2m1 - eeb2*mu1*cos3n2*sinh2m1 + ofbc*nu0*cos3n2*sinh2m1 - ofbc*nu1*cos3n2*sinh2m1 + 
       tnbc*sinn1*sinh2m1 - 48*bc*sin3n1*sinh2m1 - tnbc*sinn2*sinh2m1 + 48*bc*sin3n2*sinh2m1 + 96*b1*cos4n1*sinh3m1 - 
       192*b2*mu0*cos4n1*sinh3m1 + 192*b2*mu1*cos4n1*sinh3m1 - 96*bc*nu0*cos4n1*sinh3m1 + 96*bc*nu1*cos4n1*sinh3m1 - 
       24*b1*cos4n1*sinh3m1 + 48*b2*mu0*cos4n1*sinh3m1 - 48*b2*mu1*cos4n1*sinh3m1 + 24*bc*nu0*cos4n1*sinh3m1 - 
       24*bc*nu1*cos4n1*sinh3m1 - 96*b1*cos4n2*sinh3m1 + 192*b2*mu0*cos4n2*sinh3m1 - 192*b2*mu1*cos4n2*sinh3m1 + 
       96*bc*nu0*cos4n2*sinh3m1 - 96*bc*nu1*cos4n2*sinh3m1 + 24*b1*cos4n2*sinh3m1 - 48*b2*mu0*cos4n2*sinh3m1 + 
       48*b2*mu1*cos4n2*sinh3m1 - 24*bc*nu0*cos4n2*sinh3m1 + 24*bc*nu1*cos4n2*sinh3m1 - 48*bc*sin4n1*sinh3m1 + 6*bc*sin4n1*sinh3m1 + 
       48*bc*sin4n2*sinh3m1 - 6*bc*sin4n2*sinh3m1 + 162*b1*cosn1*sinh4m1 - 324*b2*mu0*cosn1*sinh4m1 + 324*b2*mu1*cosn1*sinh4m1 - 
       162*bc*nu0*cosn1*sinh4m1 + 162*bc*nu1*cosn1*sinh4m1 - 18*b1*cos3n1*sinh4m1 + 36*b2*mu0*cos3n1*sinh4m1 - 
       36*b2*mu1*cos3n1*sinh4m1 + 18*bc*nu0*cos3n1*sinh4m1 - 18*bc*nu1*cos3n1*sinh4m1 - 162*b1*cosn2*sinh4m1 + 
       324*b2*mu0*cosn2*sinh4m1 - 324*b2*mu1*cosn2*sinh4m1 + 162*bc*nu0*cosn2*sinh4m1 - 162*bc*nu1*cosn2*sinh4m1 + 18*b1*cos3n2*sinh4m1 - 
       36*b2*mu0*cos3n2*sinh4m1 + 36*b2*mu1*cos3n2*sinh4m1 - 18*bc*nu0*cos3n2*sinh4m1 + 18*bc*nu1*cos3n2*sinh4m1 - 
       162*bc*sinn1*sinh4m1 + 6*bc*sin3n1*sinh4m1 + 162*bc*sinn2*sinh4m1 - 6*bc*sin3n2*sinh4m1 + 2592*b1*cos4n1*sinhm2 - 
       5184*b2*mu0*cos4n1*sinhm2 + 5184*b2*mu2*cos4n1*sinhm2 - 2592*bc*nu0*cos4n1*sinhm2 + 2592*bc*nu1*cos4n1*sinhm2 - 
       xeb1*cos4n1*sinhm2 + tnb2*mu0*cos4n1*sinhm2 - tnb2*mu2*cos4n1*sinhm2 + xebc*nu0*cos4n1*sinhm2 - 
       xebc*nu1*cos4n1*sinhm2 - 2592*b1*cos4n2*sinhm2 + 5184*b2*mu0*cos4n2*sinhm2 - 5184*b2*mu2*cos4n2*sinhm2 + 
       2592*bc*nu0*cos4n2*sinhm2 - 2592*bc*nu1*cos4n2*sinhm2 + xeb1*cos4n2*sinhm2 - tnb2*mu0*cos4n2*sinhm2 + 
       tnb2*mu2*cos4n2*sinhm2 - xebc*nu0*cos4n2*sinhm2 + xebc*nu1*cos4n2*sinhm2 - tnbc*sin4n1*sinhm2 + 162*bc*sin4n1*sinhm2 + 
       tnbc*sin4n2*sinhm2 - 162*bc*sin4n2*sinhm2 + tnb1*cosn1*sinh2m2 - 2592*b2*mu0*cosn1*sinh2m2 + 2592*b2*mu2*cosn1*sinh2m2 - 
       tnbc*nu0*cosn1*sinh2m2 + tnbc*nu1*cosn1*sinh2m2 - ofb1*cos3n1*sinh2m2 + eeb2*mu0*cos3n1*sinh2m2 - 
       eeb2*mu2*cos3n1*sinh2m2 + ofbc*nu0*cos3n1*sinh2m2 - ofbc*nu1*cos3n1*sinh2m2 - tnb1*cosn2*sinh2m2 + 
       2592*b2*mu0*cosn2*sinh2m2 - 2592*b2*mu2*cosn2*sinh2m2 + tnbc*nu0*cosn2*sinh2m2 - tnbc*nu1*cosn2*sinh2m2 + 
       ofb1*cos3n2*sinh2m2 - eeb2*mu0*cos3n2*sinh2m2 + eeb2*mu2*cos3n2*sinh2m2 - ofbc*nu0*cos3n2*sinh2m2 + 
       ofbc*nu1*cos3n2*sinh2m2 - tnbc*sinn1*sinh2m2 + 48*bc*sin3n1*sinh2m2 + tnbc*sinn2*sinh2m2 - 48*bc*sin3n2*sinh2m2 - 
       96*b1*cos4n1*sinh3m2 + 192*b2*mu0*cos4n1*sinh3m2 - 192*b2*mu2*cos4n1*sinh3m2 + 96*bc*nu0*cos4n1*sinh3m2 - 
       96*bc*nu1*cos4n1*sinh3m2 + 24*b1*cos4n1*sinh3m2 - 48*b2*mu0*cos4n1*sinh3m2 + 48*b2*mu2*cos4n1*sinh3m2 - 
       24*bc*nu0*cos4n1*sinh3m2 + 24*bc*nu1*cos4n1*sinh3m2 + 96*b1*cos4n2*sinh3m2 - 192*b2*mu0*cos4n2*sinh3m2 + 
       192*b2*mu2*cos4n2*sinh3m2 - 96*bc*nu0*cos4n2*sinh3m2 + 96*bc*nu1*cos4n2*sinh3m2 - 24*b1*cos4n2*sinh3m2 + 
       48*b2*mu0*cos4n2*sinh3m2 - 48*b2*mu2*cos4n2*sinh3m2 + 24*bc*nu0*cos4n2*sinh3m2 - 24*bc*nu1*cos4n2*sinh3m2 + 
       48*bc*sin4n1*sinh3m2 - 6*bc*sin4n1*sinh3m2 - 48*bc*sin4n2*sinh3m2 + 6*bc*sin4n2*sinh3m2 + 
       6*(-27*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cosn1 + 3*b1*cos3n1 - 6*b2*mu0*cos3n1 + 6*b2*mu2*cos3n1 - 3*bc*nu0*cos3n1 + 
          3*bc*nu1*cos3n1 + 27*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cosn2 - 3*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cos3n2 + 27*bc*sinn1 - 
          bc*sin3n1 - 27*bc*sinn2 + bc*sin3n2)*sinh4m2);

  return -(Z1*val1 + Z2*val2);
}

//CPMZ marked: nu2 bug in convert function
double second_order_fzordV(double Z1, double Z2, double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc)
{
  double f = 1./864.;
  double pa3 = pow(a,3.);
  double fp = f*pa3;

  double pm02 = pow(mu0,2.);
  double pm12 = pow(mu1,2.);
  double pm22 = pow(mu2,2.);
  double pn02 = pow(nu0,2.);
  double pn12 = pow(nu1,2.);
  double pn22 = pow(nu2,2.);

  double pm012 = pow(mu0-mu1,2.);
  double pm022 = pow(mu0-mu2,2.);
  double pn012 = pow(nu0-nu1,2.);
  double pn022 = pow(nu0-nu2,2.);

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double cosh2m1 = cosh(2.*mu1);
  double cosh2m2 = cosh(2.*mu2);
  double sinh2m1 = sinh(2.*mu1);
  double sinh2m2 = sinh(2.*mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosh4m1 = cosh(4.*mu1);
  double cosh4m2 = cosh(4.*mu2);
  double sinh4m1 = sinh(4.*mu1);
  double sinh4m2 = sinh(4.*mu2);

  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cos2n1 = cos(2.*nu1);
  double cos2n2 = cos(2.*nu2);
  double sin2n1 = sin(2.*nu1);
  double sin2n2 = sin(2.*nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);
  double sin3n1 = sin(3.*nu1);
  double sin3n2 = sin(3.*nu2);
  double cos4n1 = cos(4.*nu1);
  double cos4n2 = cos(4.*nu2);
  double sin4n1 = sin(4.*nu1);
  double sin4n2 = sin(4.*nu2);

  double xsbc = 6.*bc;
  double xsb1 = 6.*b1;
  double xsb2 = 6.*b2;
  double xsc1 = 6.*c1;
  double xsc2 = 6.*c2;

  double nnbc = 9.*bc;
  double nnb1 = 9.*b1;
  double nnb2 = 9.*b2;
  double nnc1 = 9.*c1;
  double nnc2 = 9.*c2;

  double twbc = 12.*bc;
  double twb1 = 12.*b2;
  double twb2 = 12.*b1;
  double twc1 = 12.*c1;
  double twc2 = 12.*c2;

  double gtbc = 18.*bc;
  double gtb1 = 18.*b1;
  double gtb2 = 18.*b2;
  double gtc1 = 18.*c1;
  double gtc2 = 18.*c2;

  double txbc = 36.*bc;
  double txb1 = 36.*b1;
  double txb2 = 36.*b2;
  double txc1 = 36.*c1;
  double txc2 = 36.*c2;

  double rfbc = 135.*bc;
  double rfb1 = 135.*b1;
  double rfb2 = 135.*b2;
  double rfc1 = 135.*c1;
  double rfc2 = 135.*c2;

  double tdbc = 270.*bc;
  double tdb1 = 270.*b1;
  double tdb2 = 270.*b2;
  double tdc1 = 270.*c1;
  double tdc2 = 270.*c2;

  double vfbc = 540.*bc;
  double vfb1 = 540.*b1;
  double vfb2 = 540.*b2;
  double vfc1 = 540.*c1;
  double vfc2 = 540.*c2;

  double val1 = fp*(-18*a0*cos4n1*cosh3m1 - 4*b2*cos4n1*cosh3m1 + nnc2*cos4n1*cosh3m1 + gtb1*mu0*cos4n1*cosh3m1 - 
       gtb2*pm02*cos4n1*cosh3m1 - gtb1*mu1*cos4n1*cosh3m1 + txb2*mu0*mu1*cos4n1*cosh3m1 - gtb2*pm12*cos4n1*cosh3m1 + 
       gtc1*nu0*cos4n1*cosh3m1 - gtbc*mu0*nu0*cos4n1*cosh3m1 + gtbc*mu1*nu0*cos4n1*cosh3m1 - gtc2*pn02*cos4n1*cosh3m1 - 
       gtc1*nu1*cos4n1*cosh3m1 + gtbc*mu0*nu1*cos4n1*cosh3m1 - gtbc*mu1*nu1*cos4n1*cosh3m1 + txc2*nu0*nu1*cos4n1*cosh3m1 - 
       gtc2*pn12*cos4n1*cosh3m1 + 18*a0*cos4n2*cosh3m1 + 4*b2*cos4n2*cosh3m1 - nnc2*cos4n2*cosh3m1 - 
       gtb1*mu0*cos4n2*cosh3m1 + gtb2*pm02*cos4n2*cosh3m1 + gtb1*mu1*cos4n2*cosh3m1 - txb2*mu0*mu1*cos4n2*cosh3m1 + 
       gtb2*pm12*cos4n2*cosh3m1 - gtc1*nu0*cos4n2*cosh3m1 + gtbc*mu0*nu0*cos4n2*cosh3m1 - gtbc*mu1*nu0*cos4n2*cosh3m1 + 
       gtc2*pn02*cos4n2*cosh3m1 + gtc1*nu1*cos4n2*cosh3m1 - gtbc*mu0*nu1*cos4n2*cosh3m1 + gtbc*mu1*nu1*cos4n2*cosh3m1 - 
       txc2*nu0*nu1*cos4n2*cosh3m1 + gtc2*pn22*cos4n2*cosh3m1 + 270*a0*cos4n1*coshm2 + vfb2*cos4n1*coshm2 - 
       rfc2*cos4n1*coshm2 - tdb1*mu0*cos4n1*coshm2 + tdb2*pm02*cos4n1*coshm2 + tdb1*mu2*cos4n1*coshm2 - 
       vfb2*mu0*mu2*cos4n1*coshm2 + tdb2*pm22*cos4n1*coshm2 - tdc1*nu0*cos4n1*coshm2 + tdbc*mu0*nu0*cos4n1*coshm2 - 
       tdbc*mu2*nu0*cos4n1*coshm2 + tdc2*pn02*cos4n1*coshm2 + tdc1*nu1*cos4n1*coshm2 - tdbc*mu0*nu1*cos4n1*coshm2 + 
       tdbc*mu2*nu1*cos4n1*coshm2 - vfc2*nu0*nu1*cos4n1*coshm2 + tdc2*pn12*cos4n1*coshm2 - 270*a0*cos4n2*coshm2 - 
       vfb2*cos4n2*coshm2 + rfc2*cos4n2*coshm2 + tdb1*mu0*cos4n2*coshm2 - tdb2*pm02*cos4n2*coshm2 - 
       tdb1*mu2*cos4n2*coshm2 + vfb2*mu0*mu2*cos4n2*coshm2 - tdb2*pm22*cos4n2*coshm2 + tdc1*nu0*cos4n2*coshm2 - 
       tdbc*mu0*nu0*cos4n2*coshm2 + tdbc*mu2*nu0*cos4n2*coshm2 - tdc2*pn02*cos4n2*coshm2 - tdc1*nu1*cos4n2*coshm2 + 
       tdbc*mu0*nu1*cos4n2*coshm2 - tdbc*mu2*nu1*cos4n2*coshm2 + vfc2*nu0*nu1*cos4n2*coshm2 - tdc2*pn22*cos4n2*coshm2 + 
       270*a0*cosn1*cosh2m2 + rfb2*cosn1*cosh2m2 - vfc2*cosn1*cosh2m2 - tdb1*mu0*cosn1*cosh2m2 + 
       tdb2*pm02*cosn1*cosh2m2 + tdb1*mu2*cosn1*cosh2m2 - vfb2*mu0*mu2*cosn1*cosh2m2 + tdb2*pm22*cosn1*cosh2m2 - 
       tdc1*nu0*cosn1*cosh2m2 + tdbc*mu0*nu0*cosn1*cosh2m2 - tdbc*mu2*nu0*cosn1*cosh2m2 + tdc2*pn02*cosn1*cosh2m2 + 
       tdc1*nu1*cosn1*cosh2m2 - tdbc*mu0*nu1*cosn1*cosh2m2 + tdbc*mu2*nu1*cosn1*cosh2m2 - vfc2*nu0*nu1*cosn1*cosh2m2 + 
       tdc2*pn12*cosn1*cosh2m2 + 18*a0*cos3n1*cosh2m2 + nnb2*cos3n1*cosh2m2 - 4*c2*cos3n1*cosh2m2 - 
       gtb1*mu0*cos3n1*cosh2m2 + gtb2*pm02*cos3n1*cosh2m2 + gtb1*mu2*cos3n1*cosh2m2 - txb2*mu0*mu2*cos3n1*cosh2m2 + 
       gtb2*pm22*cos3n1*cosh2m2 - gtc1*nu0*cos3n1*cosh2m2 + gtbc*mu0*nu0*cos3n1*cosh2m2 - gtbc*mu2*nu0*cos3n1*cosh2m2 + 
       gtc2*pn02*cos3n1*cosh2m2 + gtc1*nu1*cos3n1*cosh2m2 - gtbc*mu0*nu1*cos3n1*cosh2m2 + gtbc*mu2*nu1*cos3n1*cosh2m2 - 
       txc2*nu0*nu1*cos3n1*cosh2m2 + gtc2*pn12*cos3n1*cosh2m2 - 270*a0*cosn2*cosh2m2 - rfb2*cosn2*cosh2m2 + 
       vfc2*cosn2*cosh2m2 + tdb1*mu0*cosn2*cosh2m2 - tdb2*pm02*cosn2*cosh2m2 - tdb1*mu2*cosn2*cosh2m2 + 
       vfb2*mu0*mu2*cosn2*cosh2m2 - tdb2*pm22*cosn2*cosh2m2 + tdc1*nu0*cosn2*cosh2m2 - tdbc*mu0*nu0*cosn2*cosh2m2 + 
       tdbc*mu2*nu0*cosn2*cosh2m2 - tdc2*pn02*cosn2*cosh2m2 - tdc1*nu1*cosn2*cosh2m2 + tdbc*mu0*nu1*cosn2*cosh2m2 - 
       tdbc*mu2*nu1*cosn2*cosh2m2 + vfc2*nu0*nu1*cosn2*cosh2m2 - tdc2*pn22*cosn2*cosh2m2 - 18*a0*cos3n2*cosh2m2 - 
       nnb2*cos3n2*cosh2m2 + 4*c2*cos3n2*cosh2m2 + gtb1*mu0*cos3n2*cosh2m2 - gtb2*pm02*cos3n2*cosh2m2 - 
       gtb1*mu2*cos3n2*cosh2m2 + txb2*mu0*mu2*cos3n2*cosh2m2 - gtb2*pm22*cos3n2*cosh2m2 + gtc1*nu0*cos3n2*cosh2m2 - 
       gtbc*mu0*nu0*cos3n2*cosh2m2 + gtbc*mu2*nu0*cos3n2*cosh2m2 - gtc2*pn02*cos3n2*cosh2m2 - gtc1*nu1*cos3n2*cosh2m2 + 
       gtbc*mu0*nu1*cos3n2*cosh2m2 - gtbc*mu2*nu1*cos3n2*cosh2m2 + txc2*nu0*nu1*cos3n2*cosh2m2 - gtc2*pn22*cos3n2*cosh2m2 + 
       18*a0*cos4n1*cosh3m2 + 4*b2*cos4n1*cosh3m2 - nnc2*cos4n1*cosh3m2 - gtb1*mu0*cos4n1*cosh3m2 + 
       gtb2*pm02*cos4n1*cosh3m2 + gtb1*mu2*cos4n1*cosh3m2 - txb2*mu0*mu2*cos4n1*cosh3m2 + gtb2*pm22*cos4n1*cosh3m2 - 
       gtc1*nu0*cos4n1*cosh3m2 + gtbc*mu0*nu0*cos4n1*cosh3m2 - gtbc*mu2*nu0*cos4n1*cosh3m2 + gtc2*pn02*cos4n1*cosh3m2 + 
       gtc1*nu1*cos4n1*cosh3m2 - gtbc*mu0*nu1*cos4n1*cosh3m2 + gtbc*mu2*nu1*cos4n1*cosh3m2 - txc2*nu0*nu1*cos4n1*cosh3m2 + 
       gtc2*pn12*cos4n1*cosh3m2 - 18*a0*cos4n2*cosh3m2 - 4*b2*cos4n2*cosh3m2 + nnc2*cos4n2*cosh3m2 + 
       gtb1*mu0*cos4n2*cosh3m2 - gtb2*pm02*cos4n2*cosh3m2 - gtb1*mu2*cos4n2*cosh3m2 + txb2*mu0*mu2*cos4n2*cosh3m2 - 
       gtb2*pm22*cos4n2*cosh3m2 + gtc1*nu0*cos4n2*cosh3m2 - gtbc*mu0*nu0*cos4n2*cosh3m2 + gtbc*mu2*nu0*cos4n2*cosh3m2 - 
       gtc2*pn02*cos4n2*cosh3m2 - gtc1*nu1*cos4n2*cosh3m2 + gtbc*mu0*nu1*cos4n2*cosh3m2 - gtbc*mu2*nu1*cos4n2*cosh3m2 + 
       txc2*nu0*nu1*cos4n2*cosh3m2 - gtc2*pn22*cos4n2*cosh3m2 - tdc1*cosh2m2*sinn1 + tdbc*mu0*cosh2m2*sinn1 - 
       tdbc*mu2*cosh2m2*sinn1 + vfc2*nu0*cosh2m2*sinn1 - vfc2*nu1*cosh2m2*sinn1 + nnc1*cosh3m1*sin4n1 - nnbc*mu0*cosh3m1*sin4n1 + 
       nnbc*mu1*cosh3m1*sin4n1 - gtc2*nu0*cosh3m1*sin4n1 + gtc2*nu1*cosh3m1*sin4n1 - rfc1*coshm2*sin4n1 + 
       rfbc*mu0*coshm2*sin4n1 - rfbc*mu2*coshm2*sin4n1 + tdc2*nu0*coshm2*sin4n1 - tdc2*nu1*coshm2*sin4n1 - nnc1*cosh3m2*sin4n1 + 
       nnbc*mu0*cosh3m2*sin4n1 - nnbc*mu2*cosh3m2*sin4n1 + gtc2*nu0*cosh3m2*sin4n1 - gtc2*nu1*cosh3m2*sin4n1 - 
       xsc1*cosh2m2*sin3n1 + xsbc*mu0*cosh2m2*sin3n1 - xsbc*mu2*cosh2m2*sin3n1 + twc2*nu0*cosh2m2*sin3n1 - 
       twc2*nu1*cosh2m2*sin3n1 + tdc1*cosh2m2*sinn2 - tdbc*mu0*cosh2m2*sinn2 + tdbc*mu2*cosh2m2*sinn2 - vfc2*nu0*cosh2m2*sinn2 + 
       vfc2*nu1*cosh2m2*sinn2 - nnc1*cosh3m1*sin4n2 + nnbc*mu0*cosh3m1*sin4n2 - nnbc*mu1*cosh3m1*sin4n2 + gtc2*nu0*cosh3m1*sin4n2 - 
       gtc2*nu1*cosh3m1*sin4n2 + rfc1*coshm2*sin4n2 - rfbc*mu0*coshm2*sin4n2 + rfbc*mu2*coshm2*sin4n2 - tdc2*nu0*coshm2*sin4n2 + 
       tdc2*nu1*coshm2*sin4n2 + nnc1*cosh3m2*sin4n2 - nnbc*mu0*cosh3m2*sin4n2 + nnbc*mu2*cosh3m2*sin4n2 - gtc2*nu0*cosh3m2*sin4n2 + 
       gtc2*nu1*cosh3m2*sin4n2 + 135*coshm1*((-2*a0 - 2*b2*(2 + pm012) + c2*(1 - 2*pn012) + 2*c1*(nu0 - nu1) + 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos4n1 + (2*a0 + 2*b2*(2 + pm012) + c2*(-1 + 2*pn022) + 2*c1*(-nu0 + nu2) - 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2 + (c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n1 - 
          (c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n2) + xsc1*cosh2m2*sin3n2 - xsbc*mu0*cosh2m2*sin3n2 + xsbc*mu2*cosh2m2*sin3n2 - 
       twc2*nu0*cosh2m2*sin3n2 + twc2*nu1*cosh2m2*sin3n2 + 
       cosh2m1*(-135*(2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn012) + 2*c1*(-nu0 + nu1) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 - 
          18*a0*cos3n1 - nnb2*cos3n1 + 4*c2*cos3n1 + gtb1*mu0*cos3n1 - gtb2*pm02*cos3n1 - gtb1*mu1*cos3n1 + txb2*mu0*mu1*cos3n1 - 
          gtb2*pm12*cos3n1 + gtc1*nu0*cos3n1 - gtbc*mu0*nu0*cos3n1 + gtbc*mu1*nu0*cos3n1 - gtc2*pn02*cos3n1 - gtc1*nu1*cos3n1 + 
          gtbc*mu0*nu1*cos3n1 - gtbc*mu1*nu1*cos3n1 + txc2*nu0*nu1*cos3n1 - gtc2*pn12*cos3n1 + 
          135*(2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 + 
          (18*a0 + nnb2*(1 + 2*pm012) + 2*(nnb1*(-mu0 + mu1) + c2*(-2 + 9*pn022) + 9*(c1 + bc*(-mu0 + mu1))*(-nu0 + nu2)))*cos3n2 + 
          tdc1*sinn1 - tdbc*mu0*sinn1 + tdbc*mu1*sinn1 - vfc2*nu0*sinn1 + vfc2*nu1*sinn1 + xsc1*sin3n1 - xsbc*mu0*sin3n1 + 
          xsbc*mu1*sin3n1 - twc2*nu0*sin3n1 + twc2*nu1*sin3n1 - 270*(c1 + bc*(-mu0 + mu1) + 2*c2*(-nu0 + nu2))*sinn2 - 
          6*(c1 + bc*(-mu0 + mu1) + 2*c2*(-nu0 + nu2))*sin3n2) + tdb1*cos4n1*sinhm1 - vfb2*mu0*cos4n1*sinhm1 + vfb2*mu1*cos4n1*sinhm1 - 
       tdbc*nu0*cos4n1*sinhm1 + tdbc*nu1*cos4n1*sinhm1 - tdb1*cos4n2*sinhm1 + vfb2*mu0*cos4n2*sinhm1 - vfb2*mu1*cos4n2*sinhm1 + 
       tdbc*nu0*cos4n2*sinhm1 - tdbc*nu1*cos4n2*sinhm1 - rfbc*sin4n1*sinhm1 + rfbc*sin4n2*sinhm1 + rfb1*cosn1*sinh2m1 - 
       tdb2*mu0*cosn1*sinh2m1 + tdb2*mu1*cosn1*sinh2m1 - rfbc*nu0*cosn1*sinh2m1 + rfbc*nu1*cosn1*sinh2m1 + nnb1*cos3n1*sinh2m1 - 
       gtb2*mu0*cos3n1*sinh2m1 + gtb2*mu1*cos3n1*sinh2m1 - nnbc*nu0*cos3n1*sinh2m1 + nnbc*nu1*cos3n1*sinh2m1 - 
       rfb1*cosn2*sinh2m1 + tdb2*mu0*cosn2*sinh2m1 - tdb2*mu1*cosn2*sinh2m1 + rfbc*nu0*cosn2*sinh2m1 - rfbc*nu1*cosn2*sinh2m1 - 
       nnb1*cos3n2*sinh2m1 + gtb2*mu0*cos3n2*sinh2m1 - gtb2*mu1*cos3n2*sinh2m1 + nnbc*nu0*cos3n2*sinh2m1 - 
       nnbc*nu1*cos3n2*sinh2m1 - rfbc*sinn1*sinh2m1 - 3*bc*sin3n1*sinh2m1 + rfbc*sinn2*sinh2m1 + 3*bc*sin3n2*sinh2m1 + 
       xsb1*cos4n1*sinh3m1 - twb2*mu0*cos4n1*sinh3m1 + twb2*mu1*cos4n1*sinh3m1 - xsbc*nu0*cos4n1*sinh3m1 + 
       xsbc*nu1*cos4n1*sinh3m1 - xsb1*cos4n2*sinh3m1 + twb2*mu0*cos4n2*sinh3m1 - twb2*mu1*cos4n2*sinh3m1 + 
       xsbc*nu0*cos4n2*sinh3m1 - xsbc*nu1*cos4n2*sinh3m1 - 3*bc*sin4n1*sinh3m1 + 3*bc*sin4n2*sinh3m1 - tdb1*cos4n1*sinhm2 + 
       vfb2*mu0*cos4n1*sinhm2 - vfb2*mu2*cos4n1*sinhm2 + tdbc*nu0*cos4n1*sinhm2 - tdbc*nu1*cos4n1*sinhm2 + tdb1*cos4n2*sinhm2 - 
       vfb2*mu0*cos4n2*sinhm2 + vfb2*mu2*cos4n2*sinhm2 - tdbc*nu0*cos4n2*sinhm2 + tdbc*nu1*cos4n2*sinhm2 + rfbc*sin4n1*sinhm2 - 
       rfbc*sin4n2*sinhm2 - rfb1*cosn1*sinh2m2 + tdb2*mu0*cosn1*sinh2m2 - tdb2*mu2*cosn1*sinh2m2 + rfbc*nu0*cosn1*sinh2m2 - 
       rfbc*nu1*cosn1*sinh2m2 - nnb1*cos3n1*sinh2m2 + gtb2*mu0*cos3n1*sinh2m2 - gtb2*mu2*cos3n1*sinh2m2 + 
       nnbc*nu0*cos3n1*sinh2m2 - nnbc*nu1*cos3n1*sinh2m2 + rfb1*cosn2*sinh2m2 - tdb2*mu0*cosn2*sinh2m2 + tdb2*mu2*cosn2*sinh2m2 - 
       rfbc*nu0*cosn2*sinh2m2 + rfbc*nu1*cosn2*sinh2m2 + nnb1*cos3n2*sinh2m2 - gtb2*mu0*cos3n2*sinh2m2 + 
       gtb2*mu2*cos3n2*sinh2m2 - nnbc*nu0*cos3n2*sinh2m2 + nnbc*nu1*cos3n2*sinh2m2 + rfbc*sinn1*sinh2m2 + 3*bc*sin3n1*sinh2m2 - 
       rfbc*sinn2*sinh2m2 - 3*bc*sin3n2*sinh2m2 - 3*
        (2*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cos4n1 - 2*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cos4n2 + bc*(-sin4n1 + sin4n2))*sinh3m2);

  double val2 = fp*(-18*a0*cos4n1*cosh3m1 - 4*b2*cos4n1*cosh3m1 + nnc2*cos4n1*cosh3m1 + gtb1*mu0*cos4n1*cosh3m1 - 
       gtb2*pm02*cos4n1*cosh3m1 - gtb1*mu1*cos4n1*cosh3m1 + txb2*mu0*mu1*cos4n1*cosh3m1 - gtb2*pm12*cos4n1*cosh3m1 + 
       gtc1*nu0*cos4n1*cosh3m1 - gtbc*mu0*nu0*cos4n1*cosh3m1 + gtbc*mu1*nu0*cos4n1*cosh3m1 - gtc2*pn02*cos4n1*cosh3m1 - 
       gtc1*nu1*cos4n1*cosh3m1 + gtbc*mu0*nu1*cos4n1*cosh3m1 - gtbc*mu1*nu1*cos4n1*cosh3m1 + txc2*nu0*nu1*cos4n1*cosh3m1 - 
       gtc2*pn12*cos4n1*cosh3m1 + 18*a0*cos4n2*cosh3m1 + 4*b2*cos4n2*cosh3m1 - nnc2*cos4n2*cosh3m1 - 
       gtb1*mu0*cos4n2*cosh3m1 + gtb2*pm02*cos4n2*cosh3m1 + gtb1*mu1*cos4n2*cosh3m1 - txb2*mu0*mu1*cos4n2*cosh3m1 + 
       gtb2*pm12*cos4n2*cosh3m1 - gtc1*nu0*cos4n2*cosh3m1 + gtbc*mu0*nu0*cos4n2*cosh3m1 - gtbc*mu1*nu0*cos4n2*cosh3m1 + 
       gtc2*pn02*cos4n2*cosh3m1 + gtc1*nu1*cos4n2*cosh3m1 - gtbc*mu0*nu1*cos4n2*cosh3m1 + gtbc*mu1*nu1*cos4n2*cosh3m1 - 
       txc2*nu0*nu1*cos4n2*cosh3m1 + gtc2*pn22*cos4n2*cosh3m1 - 162*a0*cos4n1*coshm2 - 324*b2*cos4n1*coshm2 + 
       81*c2*cos4n1*coshm2 + 162*b1*mu0*cos4n1*coshm2 - 162*b2*pm02*cos4n1*coshm2 - 162*b1*mu2*cos4n1*coshm2 + 
       324*b2*mu0*mu2*cos4n1*coshm2 - 162*b2*pm22*cos4n1*coshm2 + 162*c1*nu0*cos4n1*coshm2 - 162*bc*mu0*nu0*cos4n1*coshm2 + 
       162*bc*mu2*nu0*cos4n1*coshm2 - 162*c2*pn02*cos4n1*coshm2 - 162*c1*nu1*cos4n1*coshm2 + 162*bc*mu0*nu1*cos4n1*coshm2 - 
       162*bc*mu2*nu1*cos4n1*coshm2 + 324*c2*nu0*nu1*cos4n1*coshm2 - 162*c2*pn12*cos4n1*coshm2 + 162*a0*cos4n2*coshm2 + 
       324*b2*cos4n2*coshm2 - 81*c2*cos4n2*coshm2 - 162*b1*mu0*cos4n2*coshm2 + 162*b2*pm02*cos4n2*coshm2 + 
       162*b1*mu2*cos4n2*coshm2 - 324*b2*mu0*mu2*cos4n2*coshm2 + 162*b2*pm22*cos4n2*coshm2 - 162*c1*nu0*cos4n2*coshm2 + 
       162*bc*mu0*nu0*cos4n2*coshm2 - 162*bc*mu2*nu0*cos4n2*coshm2 + 162*c2*pn02*cos4n2*coshm2 + 162*c1*nu1*cos4n2*coshm2 - 
       162*bc*mu0*nu1*cos4n2*coshm2 + 162*bc*mu2*nu1*cos4n2*coshm2 - 324*c2*nu0*nu1*cos4n2*coshm2 + 162*c2*pn22*cos4n2*coshm2 + 
       162*a0*cosn1*cosh2m2 + 81*b2*cosn1*cosh2m2 - 324*c2*cosn1*cosh2m2 - 162*b1*mu0*cosn1*cosh2m2 + 
       162*b2*pm02*cosn1*cosh2m2 + 162*b1*mu2*cosn1*cosh2m2 - 324*b2*mu0*mu2*cosn1*cosh2m2 + 162*b2*pm22*cosn1*cosh2m2 - 
       162*c1*nu0*cosn1*cosh2m2 + 162*bc*mu0*nu0*cosn1*cosh2m2 - 162*bc*mu2*nu0*cosn1*cosh2m2 + 162*c2*pn02*cosn1*cosh2m2 + 
       162*c1*nu1*cosn1*cosh2m2 - 162*bc*mu0*nu1*cosn1*cosh2m2 + 162*bc*mu2*nu1*cosn1*cosh2m2 - 324*c2*nu0*nu1*cosn1*cosh2m2 + 
       162*c2*pn12*cosn1*cosh2m2 - 18*a0*cos3n1*cosh2m2 - nnb2*cos3n1*cosh2m2 + 4*c2*cos3n1*cosh2m2 + 
       gtb1*mu0*cos3n1*cosh2m2 - gtb2*pm02*cos3n1*cosh2m2 - gtb1*mu2*cos3n1*cosh2m2 + txb2*mu0*mu2*cos3n1*cosh2m2 - 
       gtb2*pm22*cos3n1*cosh2m2 + gtc1*nu0*cos3n1*cosh2m2 - gtbc*mu0*nu0*cos3n1*cosh2m2 + gtbc*mu2*nu0*cos3n1*cosh2m2 - 
       gtc2*pn02*cos3n1*cosh2m2 - gtc1*nu1*cos3n1*cosh2m2 + gtbc*mu0*nu1*cos3n1*cosh2m2 - gtbc*mu2*nu1*cos3n1*cosh2m2 + 
       txc2*nu0*nu1*cos3n1*cosh2m2 - gtc2*pn12*cos3n1*cosh2m2 - 162*a0*cosn2*cosh2m2 - 81*b2*cosn2*cosh2m2 + 
       324*c2*cosn2*cosh2m2 + 162*b1*mu0*cosn2*cosh2m2 - 162*b2*pm02*cosn2*cosh2m2 - 162*b1*mu2*cosn2*cosh2m2 + 
       324*b2*mu0*mu2*cosn2*cosh2m2 - 162*b2*pm22*cosn2*cosh2m2 + 162*c1*nu0*cosn2*cosh2m2 - 162*bc*mu0*nu0*cosn2*cosh2m2 + 
       162*bc*mu2*nu0*cosn2*cosh2m2 - 162*c2*pn02*cosn2*cosh2m2 - 162*c1*nu1*cosn2*cosh2m2 + 162*bc*mu0*nu1*cosn2*cosh2m2 - 
       162*bc*mu2*nu1*cosn2*cosh2m2 + 324*c2*nu0*nu1*cosn2*cosh2m2 - 162*c2*pn22*cosn2*cosh2m2 + 18*a0*cos3n2*cosh2m2 + 
       nnb2*cos3n2*cosh2m2 - 4*c2*cos3n2*cosh2m2 - gtb1*mu0*cos3n2*cosh2m2 + gtb2*pm02*cos3n2*cosh2m2 + 
       gtb1*mu2*cos3n2*cosh2m2 - txb2*mu0*mu2*cos3n2*cosh2m2 + gtb2*pm22*cos3n2*cosh2m2 - gtc1*nu0*cos3n2*cosh2m2 + 
       gtbc*mu0*nu0*cos3n2*cosh2m2 - gtbc*mu2*nu0*cos3n2*cosh2m2 + gtc2*pn02*cos3n2*cosh2m2 + gtc1*nu1*cos3n2*cosh2m2 - 
       gtbc*mu0*nu1*cos3n2*cosh2m2 + gtbc*mu2*nu1*cos3n2*cosh2m2 - txc2*nu0*nu1*cos3n2*cosh2m2 + gtc2*pn22*cos3n2*cosh2m2 + 
       18*a0*cos4n1*cosh3m2 + 4*b2*cos4n1*cosh3m2 - nnc2*cos4n1*cosh3m2 - gtb1*mu0*cos4n1*cosh3m2 + 
       gtb2*pm02*cos4n1*cosh3m2 + gtb1*mu2*cos4n1*cosh3m2 - txb2*mu0*mu2*cos4n1*cosh3m2 + gtb2*pm22*cos4n1*cosh3m2 - 
       gtc1*nu0*cos4n1*cosh3m2 + gtbc*mu0*nu0*cos4n1*cosh3m2 - gtbc*mu2*nu0*cos4n1*cosh3m2 + gtc2*pn02*cos4n1*cosh3m2 + 
       gtc1*nu1*cos4n1*cosh3m2 - gtbc*mu0*nu1*cos4n1*cosh3m2 + gtbc*mu2*nu1*cos4n1*cosh3m2 - txc2*nu0*nu1*cos4n1*cosh3m2 + 
       gtc2*pn12*cos4n1*cosh3m2 - 18*a0*cos4n2*cosh3m2 - 4*b2*cos4n2*cosh3m2 + nnc2*cos4n2*cosh3m2 + 
       gtb1*mu0*cos4n2*cosh3m2 - gtb2*pm02*cos4n2*cosh3m2 - gtb1*mu2*cos4n2*cosh3m2 + txb2*mu0*mu2*cos4n2*cosh3m2 - 
       gtb2*pm22*cos4n2*cosh3m2 + gtc1*nu0*cos4n2*cosh3m2 - gtbc*mu0*nu0*cos4n2*cosh3m2 + gtbc*mu2*nu0*cos4n2*cosh3m2 - 
       gtc2*pn02*cos4n2*cosh3m2 - gtc1*nu1*cos4n2*cosh3m2 + gtbc*mu0*nu1*cos4n2*cosh3m2 - gtbc*mu2*nu1*cos4n2*cosh3m2 + 
       txc2*nu0*nu1*cos4n2*cosh3m2 - gtc2*pn22*cos4n2*cosh3m2 - 162*c1*cosh2m2*sinn1 + 162*bc*mu0*cosh2m2*sinn1 - 
       162*bc*mu2*cosh2m2*sinn1 + 324*c2*nu0*cosh2m2*sinn1 - 324*c2*nu1*cosh2m2*sinn1 + nnc1*cosh3m1*sin4n1 - nnbc*mu0*cosh3m1*sin4n1 + 
       nnbc*mu1*cosh3m1*sin4n1 - gtc2*nu0*cosh3m1*sin4n1 + gtc2*nu1*cosh3m1*sin4n1 + 81*c1*coshm2*sin4n1 - 81*bc*mu0*coshm2*sin4n1 + 
       81*bc*mu2*coshm2*sin4n1 - 162*c2*nu0*coshm2*sin4n1 + 162*c2*nu1*coshm2*sin4n1 - nnc1*cosh3m2*sin4n1 + nnbc*mu0*cosh3m2*sin4n1 - 
       nnbc*mu2*cosh3m2*sin4n1 + gtc2*nu0*cosh3m2*sin4n1 - gtc2*nu1*cosh3m2*sin4n1 + xsc1*cosh2m2*sin3n1 - 
       xsbc*mu0*cosh2m2*sin3n1 + xsbc*mu2*cosh2m2*sin3n1 - twc2*nu0*cosh2m2*sin3n1 + twc2*nu1*cosh2m2*sin3n1 + 
       162*c1*cosh2m2*sinn2 - 162*bc*mu0*cosh2m2*sinn2 + 162*bc*mu2*cosh2m2*sinn2 - 324*c2*nu0*cosh2m2*sinn2 + 324*c2*nu1*cosh2m2*sinn2 - 
       nnc1*cosh3m1*sin4n2 + nnbc*mu0*cosh3m1*sin4n2 - nnbc*mu1*cosh3m1*sin4n2 + gtc2*nu0*cosh3m1*sin4n2 - 
       gtc2*nu1*cosh3m1*sin4n2 - 81*c1*coshm2*sin4n2 + 81*bc*mu0*coshm2*sin4n2 - 81*bc*mu2*coshm2*sin4n2 + 162*c2*nu0*coshm2*sin4n2 - 
       162*c2*nu1*coshm2*sin4n2 + nnc1*cosh3m2*sin4n2 - nnbc*mu0*cosh3m2*sin4n2 + nnbc*mu2*cosh3m2*sin4n2 - gtc2*nu0*cosh3m2*sin4n2 + 
       gtc2*nu1*cosh3m2*sin4n2 + 81*coshm1*((2*a0 + 2*b2*(2 + pm012) + c2*(-1 + 2*pn012) + 2*c1*(-nu0 + nu1) - 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos4n1 + (-2*a0 - 2*b2*(2 + pm012) + c2*(1 - 2*pn022) + 2*c1*(nu0 - nu2) + 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2 - (c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n1 + 
          (c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n2) - xsc1*cosh2m2*sin3n2 + xsbc*mu0*cosh2m2*sin3n2 - xsbc*mu2*cosh2m2*sin3n2 + 
       twc2*nu0*cosh2m2*sin3n2 - twc2*nu1*cosh2m2*sin3n2 + 
       cosh2m1*(-81*(2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn012) + 2*c1*(-nu0 + nu1) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 + 
          18*a0*cos3n1 + nnb2*cos3n1 - 4*c2*cos3n1 - gtb1*mu0*cos3n1 + gtb2*pm02*cos3n1 + gtb1*mu1*cos3n1 - txb2*mu0*mu1*cos3n1 + 
          gtb2*pm12*cos3n1 - gtc1*nu0*cos3n1 + gtbc*mu0*nu0*cos3n1 - gtbc*mu1*nu0*cos3n1 + gtc2*pn02*cos3n1 + gtc1*nu1*cos3n1 - 
          gtbc*mu0*nu1*cos3n1 + gtbc*mu1*nu1*cos3n1 - txc2*nu0*nu1*cos3n1 + gtc2*pn12*cos3n1 + 
          81*(2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 - 
          (18*a0 + nnb2*(1 + 2*pm012) + 2*(nnb1*(-mu0 + mu1) + c2*(-2 + 9*pn022) + 9*(c1 + bc*(-mu0 + mu1))*(-nu0 + nu2)))*cos3n2 + 
          162*c1*sinn1 - 162*bc*mu0*sinn1 + 162*bc*mu1*sinn1 - 324*c2*nu0*sinn1 + 324*c2*nu1*sinn1 - xsc1*sin3n1 + xsbc*mu0*sin3n1 - 
          xsbc*mu1*sin3n1 + twc2*nu0*sin3n1 - twc2*nu1*sin3n1 - 162*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sinn2 + 
          6*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin3n2) - 162*b1*cos4n1*sinhm1 + 324*b2*mu0*cos4n1*sinhm1 - 324*b2*mu1*cos4n1*sinhm1 + 
       162*bc*nu0*cos4n1*sinhm1 - 162*bc*nu1*cos4n1*sinhm1 + 162*b1*cos4n2*sinhm1 - 324*b2*mu0*cos4n2*sinhm1 + 324*b2*mu1*cos4n2*sinhm1 - 
       162*bc*nu0*cos4n2*sinhm1 + 162*bc*nu1*cos4n2*sinhm1 + 81*bc*sin4n1*sinhm1 - 81*bc*sin4n2*sinhm1 + 81*b1*cosn1*sinh2m1 - 
       162*b2*mu0*cosn1*sinh2m1 + 162*b2*mu1*cosn1*sinh2m1 - 81*bc*nu0*cosn1*sinh2m1 + 81*bc*nu1*cosn1*sinh2m1 - nnb1*cos3n1*sinh2m1 + 
       gtb2*mu0*cos3n1*sinh2m1 - gtb2*mu1*cos3n1*sinh2m1 + nnbc*nu0*cos3n1*sinh2m1 - nnbc*nu1*cos3n1*sinh2m1 - 81*b1*cosn2*sinh2m1 + 
       162*b2*mu0*cosn2*sinh2m1 - 162*b2*mu1*cosn2*sinh2m1 + 81*bc*nu0*cosn2*sinh2m1 - 81*bc*nu1*cosn2*sinh2m1 + nnb1*cos3n2*sinh2m1 - 
       gtb2*mu0*cos3n2*sinh2m1 + gtb2*mu1*cos3n2*sinh2m1 - nnbc*nu0*cos3n2*sinh2m1 + nnbc*nu1*cos3n2*sinh2m1 - 81*bc*sinn1*sinh2m1 + 
       3*bc*sin3n1*sinh2m1 + 81*bc*sinn2*sinh2m1 - 3*bc*sin3n2*sinh2m1 + xsb1*cos4n1*sinh3m1 - twb2*mu0*cos4n1*sinh3m1 + 
       twb2*mu1*cos4n1*sinh3m1 - xsbc*nu0*cos4n1*sinh3m1 + xsbc*nu1*cos4n1*sinh3m1 - xsb1*cos4n2*sinh3m1 + 
       twb2*mu0*cos4n2*sinh3m1 - twb2*mu1*cos4n2*sinh3m1 + xsbc*nu0*cos4n2*sinh3m1 - xsbc*nu1*cos4n2*sinh3m1 - 
       3*bc*sin4n1*sinh3m1 + 3*bc*sin4n2*sinh3m1 + 162*b1*cos4n1*sinhm2 - 324*b2*mu0*cos4n1*sinhm2 + 324*b2*mu2*cos4n1*sinhm2 - 
       162*bc*nu0*cos4n1*sinhm2 + 162*bc*nu1*cos4n1*sinhm2 - 162*b1*cos4n2*sinhm2 + 324*b2*mu0*cos4n2*sinhm2 - 324*b2*mu2*cos4n2*sinhm2 + 
       162*bc*nu0*cos4n2*sinhm2 - 162*bc*nu1*cos4n2*sinhm2 - 81*bc*sin4n1*sinhm2 + 81*bc*sin4n2*sinhm2 - 81*b1*cosn1*sinh2m2 + 
       162*b2*mu0*cosn1*sinh2m2 - 162*b2*mu2*cosn1*sinh2m2 + 81*bc*nu0*cosn1*sinh2m2 - 81*bc*nu1*cosn1*sinh2m2 + nnb1*cos3n1*sinh2m2 - 
       gtb2*mu0*cos3n1*sinh2m2 + gtb2*mu2*cos3n1*sinh2m2 - nnbc*nu0*cos3n1*sinh2m2 + nnbc*nu1*cos3n1*sinh2m2 + 81*b1*cosn2*sinh2m2 - 
       162*b2*mu0*cosn2*sinh2m2 + 162*b2*mu2*cosn2*sinh2m2 - 81*bc*nu0*cosn2*sinh2m2 + 81*bc*nu1*cosn2*sinh2m2 - nnb1*cos3n2*sinh2m2 + 
       gtb2*mu0*cos3n2*sinh2m2 - gtb2*mu2*cos3n2*sinh2m2 + nnbc*nu0*cos3n2*sinh2m2 - nnbc*nu1*cos3n2*sinh2m2 + 81*bc*sinn1*sinh2m2 - 
       3*bc*sin3n1*sinh2m2 - 81*bc*sinn2*sinh2m2 + 3*bc*sin3n2*sinh2m2 - 
       3*(2*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cos4n1 - 2*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cos4n2 + bc*(-sin4n1 + sin4n2))*sinh3m2);

  return -(Z1*val1 + Z2*val2);
}

//CPMZ marked: nu2 bug in convert function
double second_order_fzzordV(double Z1, double Z2, double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc)
{
  double f = 1./55296.;
  double pa4 = f*pow(a,4.);
  double fp = pa4;

  double pm02 = pow(mu0,2.);
  double pm12 = pow(mu1,2.);
  double pm22 = pow(mu2,2.);
  double pn02 = pow(nu0,2.);
  double pn12 = pow(nu1,2.);
  double pn22 = pow(nu2,2.);

  double pm012 = pow(mu0-mu1,2.);
  double pm022 = pow(mu0-mu2,2.);
  double pn012 = pow(nu0-nu1,2.);
  double pn022 = pow(nu0-nu2,2.);

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double cosh2m1 = cosh(2.*mu1);
  double cosh2m2 = cosh(2.*mu2);
  double sinh2m1 = sinh(2.*mu1);
  double sinh2m2 = sinh(2.*mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosh4m1 = cosh(4.*mu1);
  double cosh4m2 = cosh(4.*mu2);
  double sinh4m1 = sinh(4.*mu1);
  double sinh4m2 = sinh(4.*mu2);

  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cos2n1 = cos(2.*nu1);
  double cos2n2 = cos(2.*nu2);
  double sin2n1 = sin(2.*nu1);
  double sin2n2 = sin(2.*nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);
  double sin3n1 = sin(3.*nu1);
  double sin3n2 = sin(3.*nu2);
  double cos4n1 = cos(4.*nu1);
  double cos4n2 = cos(4.*nu2);
  double sin4n1 = sin(4.*nu1);
  double sin4n2 = sin(4.*nu2);
  double cosn13 = cosn1*cosn1*cosn1;

  double exp2m1 = exp(2.*mu1);
  double exp3m1 = exp(3.*mu1);
  double exp4m1 = exp(4.*mu1);
  double exp5m1 = exp(5.*mu1);
  double exp6m1 = exp(6.*mu1);
  double exp7m1 = exp(7.*mu1);
  double exp8m1 = exp(8.*mu1);

  double xsbc = 6.*bc;
  double xsb1 = 6.*b1;
  double xsb2 = 6.*b2;
  double xsc1 = 6.*c1;
  double xsc2 = 6.*c2;

  double nnbc = 9.*bc;
  double nnb1 = 9.*b1;
  double nnb2 = 9.*b2;
  double nnc1 = 9.*c1;
  double nnc2 = 9.*c2;

  double sxbc = 16.*bc;
  double sxb1 = 16.*b1;
  double sxb2 = 16.*b2;
  double sxc1 = 16.*c1;
  double sxc2 = 16.*c2;

  double gtbc = 18.*bc;
  double gtb1 = 18.*b1;
  double gtb2 = 18.*b2;
  double gtc1 = 18.*c1;
  double gtc2 = 18.*c2;

  double txbc = 36.*bc;
  double txb1 = 36.*b1;
  double txb2 = 36.*b2;
  double txc1 = 36.*c1;
  double txc2 = 36.*c2;

  double sfbc = 64.*bc;
  double sfb1 = 64.*b1;
  double sfb2 = 64.*b2;
  double sfc1 = 64.*c1;
  double sfc2 = 64.*c2;

  double stbc = 72.*bc;
  double stb1 = 72.*b1;
  double stb2 = 72.*b2;
  double stc1 = 72.*c1;
  double stc2 = 72.*c2;

  double ofbc = 144.*bc;
  double ofb1 = 144.*b1;
  double ofb2 = 144.*b2;
  double ofc1 = 144.*c1;
  double ofc2 = 144.*c2;

  double tsbc = 216.*bc;
  double tsb1 = 216.*b1;
  double tsb2 = 216.*b2;
  double tsc1 = 216.*c1;
  double tsc2 = 216.*c2;

  double eebc = 288.*bc;
  double eeb1 = 288.*b1;
  double eeb2 = 288.*b2;
  double eec1 = 288.*c1;
  double eec2 = 288.*c2;

  double ftbc = 432.*bc;
  double ftb1 = 432.*b1;
  double ftb2 = 432.*b2;
  double ftc1 = 432.*c1;
  double ftc2 = 432.*c2;

  double fnbc = 576.*bc;
  double fnb1 = 576.*b1;
  double fnb2 = 576.*b2;
  double fnc1 = 576.*c1;
  double fnc2 = 576.*c2;

  double xebc = 648.*bc;
  double xeb1 = 648.*b1;
  double xeb2 = 648.*b2;
  double xec1 = 648.*c1;
  double xec2 = 648.*c2;

  double tnbc = 1296.*bc;
  double tnb1 = 1296.*b1;
  double tnb2 = 1296.*b2;
  double tnc1 = 1296.*c1;
  double tnc2 = 1296.*c2;

  double val1 = fp*(16*exp(2*mu1 + 4*mu2)*(1 + exp4m1)*(189*(2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn012) + 2*c1*(-nu0 + nu1) - 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 - (18*a0 + nnb2*(1 + 2*pm012) + 
             2*(c2*(-2 + 9*pn012) + nnc1*(-nu0 + nu1) - 9*(mu0 - mu1)*(b1 - bc*nu0 + bc*nu1)))*cos3n1 - 
          189*(2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 + 
          (18*a0 + nnb2*(1 + 2*pm012) + 2*(c2*(-2 + 9*pn022) + nnc1*(-nu0 + nu2) - 9*(mu0 - mu1)*(b1 - bc*nu0 + bc*nu1)))*cos3n2) + 
       27*exp(3*mu1 + 4*mu2)*(1 + exp2m1)*(112*(2*a0 + 2*b2*(2 + pm012) + c2*(-1 + 2*pn012) + 2*c1*(-nu0 + nu1) - 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos4n1 + (-8*a0 - 8*b2*(2 + pm012) + c2*(1 - 8*pn012) + 8*c1*(nu0 - nu1) + 
             8*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos4n1 - 112*(2*a0 + 2*b2*(2 + pm012) + c2*(-1 + 2*pn022) + 2*c1*(-nu0 + nu2) - 
             2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2 + (8*a0 + 8*b2*(2 + pm012) + c2*(-1 + 8*pn022) + 8*c1*(-nu0 + nu2) - 
             8*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2) + exp4m1*
        (216*a0*cosn1 + 54*b1*cosn1 + 27*b2*cosn1 - ftc2*cosn1 - tsb1*mu0*cosn1 - 108*b2*mu0*cosn1 + tsb2*pm02*cosn1 + 
          tsb1*mu2*cosn1 + 108*b2*mu2*cosn1 - ftb2*mu0*mu2*cosn1 + tsb2*pm22*cosn1 - 54*bc*nu0*cosn1 - tsc1*nu0*cosn1 + 
          tsbc*mu0*nu0*cosn1 - tsbc*mu2*nu0*cosn1 + tsc2*pn02*cosn1 + 54*bc*nu1*cosn1 + tsc1*nu1*cosn1 - tsbc*mu0*nu1*cosn1 + 
          tsbc*mu2*nu1*cosn1 - ftc2*nu0*nu1*cosn1 + tsc2*pn12*cosn1 + 72*a0*cos3n1 + gtb1*cos3n1 + nnb2*cos3n1 - sxc2*cos3n1 - 
          stb1*mu0*cos3n1 - txb2*mu0*cos3n1 + stb2*pm02*cos3n1 + stb1*mu2*cos3n1 + txb2*mu2*cos3n1 - ofb2*mu0*mu2*cos3n1 + 
          stb2*pm22*cos3n1 - gtbc*nu0*cos3n1 - stc1*nu0*cos3n1 + stbc*mu0*nu0*cos3n1 - stbc*mu2*nu0*cos3n1 + stc2*pn02*cos3n1 + 
          gtbc*nu1*cos3n1 + stc1*nu1*cos3n1 - stbc*mu0*nu1*cos3n1 + stbc*mu2*nu1*cos3n1 - ofc2*nu0*nu1*cos3n1 + stc2*pn12*cos3n1 - 
          216*a0*cosn2 - 54*b1*cosn2 - 27*b2*cosn2 + ftc2*cosn2 + tsb1*mu0*cosn2 + 108*b2*mu0*cosn2 - tsb2*pm02*cosn2 - 
          tsb1*mu2*cosn2 - 108*b2*mu2*cosn2 + ftb2*mu0*mu2*cosn2 - tsb2*pm22*cosn2 + 54*bc*nu0*cosn2 + tsc1*nu0*cosn2 - 
          tsbc*mu0*nu0*cosn2 + tsbc*mu2*nu0*cosn2 - tsc2*pn02*cosn2 - 54*bc*nu1*cosn2 - tsc1*nu1*cosn2 + tsbc*mu0*nu1*cosn2 - 
          tsbc*mu2*nu1*cosn2 + ftc2*nu0*nu1*cosn2 - tsc2*pn22*cosn2 - 72*a0*cos3n2 - gtb1*cos3n2 - nnb2*cos3n2 + sxc2*cos3n2 + 
          stb1*mu0*cos3n2 + txb2*mu0*cos3n2 - stb2*pm02*cos3n2 - stb1*mu2*cos3n2 - txb2*mu2*cos3n2 + ofb2*mu0*mu2*cos3n2 - 
          stb2*pm22*cos3n2 + gtbc*nu0*cos3n2 + stc1*nu0*cos3n2 - stbc*mu0*nu0*cos3n2 + stbc*mu2*nu0*cos3n2 - stc2*pn02*cos3n2 - 
          gtbc*nu1*cos3n2 - stc1*nu1*cos3n2 + stbc*mu0*nu1*cos3n2 - stbc*mu2*nu1*cos3n2 + ofc2*nu0*nu1*cos3n2 - stc2*pn22*cos3n2 + 
          exp(8*mu2)*(27*(8*a0 + b1*(-2 - 8*mu0 + 8*mu2) + b2*(1 + 4*mu0 + 8*pm02 - 4*(1 + 4*mu0)*mu2 + 8*pm22) + 8*c2*(-2 + pn012) + 
                2*(bc - 4*c1 + 4*bc*mu0 - 4*bc*mu2)*(nu0 - nu1))*cosn1 + 
             (72*a0 - gtb1*(1 + 4*mu0 - 4*mu2) + nnb2*(1 + 4*mu0 + 8*pm02 - 4*(1 + 4*mu0)*mu2 + 8*pm22) + 8*c2*(-2 + 9*pn012) + 
                18*(bc - 4*c1 + 4*bc*mu0 - 4*bc*mu2)*(nu0 - nu1))*cos3n1 - 
             27*(8*a0 + b1*(-2 - 8*mu0 + 8*mu2) + b2*(1 + 4*mu0 + 8*pm02 - 4*(1 + 4*mu0)*mu2 + 8*pm22) + 8*c2*(-2 + pn022) + 
                2*(bc - 4*c1 + 4*bc*mu0 - 4*bc*mu2)*(nu0 - nu2))*cosn2 - 
             (72*a0 - gtb1*(1 + 4*mu0 - 4*mu2) + nnb2*(1 + 4*mu0 + 8*pm02 - 4*(1 + 4*mu0)*mu2 + 8*pm22) + 8*c2*(-2 + 9*pn022) + 
                18*(bc - 4*c1 + 4*bc*mu0 - 4*bc*mu2)*(nu0 - nu2))*cos3n2) + 
          16*exp(6*mu2)*(-189*(2*a0 + b1*(-1 - 2*mu0 + 2*mu2) + b2*(1 - 2*mu2 + 2*(mu0 + pm02 - 2*mu0*mu2 + pm22)) + 2*c2*(-2 + pn012) + 
                (bc - 2*c1 + 2*bc*mu0 - 2*bc*mu2)*(nu0 - nu1))*cosn1 + 
             (18*a0 - nnb1*(1 + 2*mu0 - 2*mu2) + nnb2*(1 - 2*mu2 + 2*(mu0 + pm02 - 2*mu0*mu2 + pm22)) + 2*c2*(-2 + 9*pn012) + 
                9*(bc - 2*c1 + 2*bc*mu0 - 2*bc*mu2)*(nu0 - nu1))*cos3n1 + 
             189*(2*a0 + b1*(-1 - 2*mu0 + 2*mu2) + b2*(1 - 2*mu2 + 2*(mu0 + pm02 - 2*mu0*mu2 + pm22)) + 2*c2*(-2 + pn022) + 
                (bc - 2*c1 + 2*bc*mu0 - 2*bc*mu2)*(nu0 - nu2))*cosn2 - 
             (18*a0 - nnb1*(1 + 2*mu0 - 2*mu2) + nnb2*(1 - 2*mu2 + 2*(mu0 + pm02 - 2*mu0*mu2 + pm22)) + 2*c2*(-2 + 9*pn022) + 
                9*(bc - 2*c1 + 2*bc*mu0 - 2*bc*mu2)*(nu0 - nu2))*cos3n2) + 
          16*exp(2*mu2)*(-189*(2*a0 + b1*(1 - 2*mu0 + 2*mu2) + b2*(1 + 2*pm02 + 2*mu2*(1 + mu2) - 2*mu0*(1 + 2*mu2)) + 2*c2*(-2 + pn012) - 
                (bc + 2*c1 - 2*bc*mu0 + 2*bc*mu2)*(nu0 - nu1))*cosn1 + 
             (18*a0 + nnb1*(1 - 2*mu0 + 2*mu2) + nnb2*(1 + 2*pm02 + 2*mu2*(1 + mu2) - 2*mu0*(1 + 2*mu2)) + 2*c2*(-2 + 9*pn012) - 
                9*(bc + 2*c1 - 2*bc*mu0 + 2*bc*mu2)*(nu0 - nu1))*cos3n1 + 
             189*(2*a0 + b1*(1 - 2*mu0 + 2*mu2) + b2*(1 + 2*pm02 + 2*mu2*(1 + mu2) - 2*mu0*(1 + 2*mu2)) + 2*c2*(-2 + pn022) - 
                (bc + 2*c1 - 2*bc*mu0 + 2*bc*mu2)*(nu0 - nu2))*cosn2 - 
             (18*a0 + nnb1*(1 - 2*mu0 + 2*mu2) + nnb2*(1 + 2*pm02 + 2*mu2*(1 + mu2) - 2*mu0*(1 + 2*mu2)) + 2*c2*(-2 + 9*pn022) - 
                9*(bc + 2*c1 - 2*bc*mu0 + 2*bc*mu2)*(nu0 - nu2))*cos3n2) + 
          exp(7*mu2)*(16*(18*a0 - xsb1*(1 + 3*mu0 - 3*mu2) + 2*b2*(2 + 6*mu0 + 9*pm02 - 6*(1 + 3*mu0)*mu2 + 9*pm22) + 
                nnc2*(-1 + 2*pn012) + 6*(bc - 3*c1 + 3*bc*mu0 - 3*bc*mu2)*(nu0 - nu1))*cos4n1 + 
             (72*a0 - 24*b1*(1 + 3*mu0 - 3*mu2) + 8*b2*(2 + 6*mu0 + 9*pm02 - 6*(1 + 3*mu0)*mu2 + 9*pm22) + nnc2*(-1 + 8*pn012) + 
                24*(bc - 3*c1 + 3*bc*mu0 - 3*bc*mu2)*(nu0 - nu1))*cos4n1 - 
             16*(18*a0 - xsb1*(1 + 3*mu0 - 3*mu2) + 2*b2*(2 + 6*mu0 + 9*pm02 - 6*(1 + 3*mu0)*mu2 + 9*pm22) + nnc2*(-1 + 2*pn022) + 
                6*(bc - 3*c1 + 3*bc*mu0 - 3*bc*mu2)*(nu0 - nu2))*cos4n2 - 
             (72*a0 - 24*b1*(1 + 3*mu0 - 3*mu2) + 8*b2*(2 + 6*mu0 + 9*pm02 - 6*(1 + 3*mu0)*mu2 + 9*pm22) + nnc2*(-1 + 8*pn022) + 
                24*(bc - 3*c1 + 3*bc*mu0 - 3*bc*mu2)*(nu0 - nu2))*cos4n2) + 
          exp(mu2)*(16*(18*a0 + xsb1*(1 - 3*mu0 + 3*mu2) + 2*b2*(2 + 9*pm02 + 6*mu2 + 9*pm22 - 6*mu0*(1 + 3*mu2)) + 
                nnc2*(-1 + 2*pn012) - 6*(bc + 3*c1 - 3*bc*mu0 + 3*bc*mu2)*(nu0 - nu1))*cos4n1 + 
             (72*a0 + 24*b1*(1 - 3*mu0 + 3*mu2) + 8*b2*(2 + 9*pm02 + 6*mu2 + 9*pm22 - 6*mu0*(1 + 3*mu2)) + nnc2*(-1 + 8*pn012) - 
                24*(bc + 3*c1 - 3*bc*mu0 + 3*bc*mu2)*(nu0 - nu1))*cos4n1 - 
             16*(18*a0 + xsb1*(1 - 3*mu0 + 3*mu2) + 2*b2*(2 + 9*pm02 + 6*mu2 + 9*pm22 - 6*mu0*(1 + 3*mu2)) + nnc2*(-1 + 2*pn022) - 
                6*(bc + 3*c1 - 3*bc*mu0 + 3*bc*mu2)*(nu0 - nu2))*cos4n2 - 
             (72*a0 + 24*b1*(1 - 3*mu0 + 3*mu2) + 8*b2*(2 + 9*pm02 + 6*mu2 + 9*pm22 - 6*mu0*(1 + 3*mu2)) + nnc2*(-1 + 8*pn022) - 
                24*(bc + 3*c1 - 3*bc*mu0 + 3*bc*mu2)*(nu0 - nu2))*cos4n2) + 
          27*exp(5*mu2)*(-112*(2*a0 - 2*b1*(1 + mu0 - mu2) + 2*b2*(2 + 2*mu0 + pm02 - 2*(1 + mu0)*mu2 + pm22) + c2*(-1 + 2*pn012) + 
                2*(-c1 + bc*(1 + mu0 - mu2))*(nu0 - nu1))*cos4n1 + 
             (8*a0 - 8*b1*(1 + mu0 - mu2) + 8*b2*(2 + 2*mu0 + pm02 - 2*(1 + mu0)*mu2 + pm22) + c2*(-1 + 8*pn012) + 
                8*(-c1 + bc*(1 + mu0 - mu2))*(nu0 - nu1))*cos4n1 + 
             112*(2*a0 - 2*b1*(1 + mu0 - mu2) + 2*b2*(2 + 2*mu0 + pm02 - 2*(1 + mu0)*mu2 + pm22) + c2*(-1 + 2*pn022) + 
                2*(-c1 + bc*(1 + mu0 - mu2))*(nu0 - nu2))*cos4n2 + 
             (-8*a0 + 8*b1*(1 + mu0 - mu2) - 8*b2*(2 + 2*mu0 + pm02 - 2*(1 + mu0)*mu2 + pm22) + c2*(1 - 8*pn022) + 
                8*(c1 + bc*(-1 - mu0 + mu2))*(nu0 - nu2))*cos4n2) + 
          27*exp(3*mu2)*(-112*(2*a0 + 2*b1*(1 - mu0 + mu2) + 2*b2*(2 + pm02 - 2*mu0*(1 + mu2) + mu2*(2 + mu2)) + c2*(-1 + 2*pn012) - 
                2*(c1 + bc*(1 - mu0 + mu2))*(nu0 - nu1))*cos4n1 + 
             (8*a0 + 8*b1*(1 - mu0 + mu2) + 8*b2*(2 + pm02 - 2*mu0*(1 + mu2) + mu2*(2 + mu2)) + c2*(-1 + 8*pn012) - 
                8*(c1 + bc*(1 - mu0 + mu2))*(nu0 - nu1))*cos4n1 + 
             112*(2*a0 + 2*b1*(1 - mu0 + mu2) + 2*b2*(2 + pm02 - 2*mu0*(1 + mu2) + mu2*(2 + mu2)) + c2*(-1 + 2*pn022) - 
                2*(c1 + bc*(1 - mu0 + mu2))*(nu0 - nu2))*cos4n2 - 
             (8*a0 + 8*b1*(1 - mu0 + mu2) + 8*b2*(2 + pm02 - 2*mu0*(1 + mu2) + mu2*(2 + mu2)) + c2*(-1 + 8*pn022) - 
                8*(c1 + bc*(1 - mu0 + mu2))*(nu0 - nu2))*cos4n2)) - 54*bc*exp4m1*sinn1 - tsc1*exp4m1*sinn1 + 
       tsbc*exp4m1*mu0*sinn1 - tsbc*exp4m1*mu2*sinn1 + ftc2*exp4m1*nu0*sinn1 - ftc2*exp4m1*nu1*sinn1 - 
       xsbc*exp4m1*sin3n1 - 24*c1*exp4m1*sin3n1 + 24*bc*exp4m1*mu0*sin3n1 - 24*bc*exp4m1*mu2*sin3n1 + 
       48*c2*exp4m1*nu0*sin3n1 - 48*c2*exp4m1*nu1*sin3n1 + 54*bc*exp4m1*sinn2 + tsc1*exp4m1*sinn2 - 
       tsbc*exp4m1*mu0*sinn2 + tsbc*exp4m1*mu2*sinn2 - ftc2*exp4m1*nu0*sinn2 + ftc2*exp4m1*nu1*sinn2 + 
       xsbc*exp4m1*sin3n2 + 24*c1*exp4m1*sin3n2 - 24*bc*exp4m1*mu0*sin3n2 + 24*bc*exp4m1*mu2*sin3n2 - 
       48*c2*exp4m1*nu0*sin3n2 + 48*c2*exp4m1*nu1*sin3n2 + 
       48*exp(4*mu1 + 6*mu2)*((bc - 2*c1 + 2*bc*mu0 - 2*bc*mu2 + 4*c2*nu0 - 4*c2*nu1)*(-63*sinn1 + sin3n1) + 
          63*(bc - 2*c1 + 2*bc*mu0 - 2*bc*mu2 + 4*c2*nu0 - 4*c2*nu1)*sinn2 - (bc - 2*c1 + 2*bc*mu0 - 2*bc*mu2 + 4*c2*nu0 - 4*c2*nu1)*sin3n2) + 
       48*exp(4*mu1 + 2*mu2)*(-((bc + 2*c1 - 2*bc*mu0 + 2*bc*mu2 - 4*c2*nu0 + 4*c2*nu1)*(-63*sinn1 + sin3n1)) - 
          63*(bc + 2*c1 - 2*bc*mu0 + 2*bc*mu2 - 4*c2*nu0 + 4*c2*nu1)*sinn2 + (bc + 2*c1 - 2*bc*mu0 + 2*bc*mu2 - 4*c2*nu0 + 4*c2*nu1)*sin3n2) + 
       6*exp(4*mu1 + 8*mu2)*(9*(bc - 4*c1 + 4*bc*mu0 - 4*bc*mu2 + 8*c2*nu0 - 8*c2*nu1)*sinn1 + (bc - 4*c1 + 4*bc*mu0 - 4*bc*mu2 + 8*c2*nu0 - 8*c2*nu1)*sin3n1 - 
          (bc - 4*c1 + 4*bc*mu0 - 4*bc*mu2 + 8*c2*nu0 - 8*c2*nu1)*(9*sinn2 + sin3n2)) + 
       54*exp(4*mu1 + 3*mu2)*(-2*(c1 + bc*(1 - mu0 + mu2) - 2*c2*nu0 + 2*c2*nu1)*(-28 + cos4n1)*sin4n1 - 
          (c1 + bc*(1 - mu0 + mu2) - 2*c2*nu0 + 2*c2*nu1)*(56*sin4n2 - sin4n2)) + 
       6*exp(4*mu1 + mu2)*(-((bc + 3*c1 - 3*bc*mu0 + 3*bc*mu2 - xsc2*nu0 + xsc2*nu1)*(8*sin4n1 + sin4n1)) + 
          8*(bc + 3*c1 - 3*bc*mu0 + 3*bc*mu2 - xsc2*nu0 + xsc2*nu1)*sin4n2 + (bc + 3*c1 - 3*bc*mu0 + 3*bc*mu2 - xsc2*nu0 + xsc2*nu1)*sin4n2) + 
       54*exp(4*mu1 + 5*mu2)*(2*(-c1 + bc*(1 + mu0 - mu2) + 2*c2*nu0 - 2*c2*nu1)*(-28 + cos4n1)*sin4n1 + 
          (c1 + bc*(-1 - mu0 + mu2) - 2*c2*nu0 + 2*c2*nu1)*(-56*sin4n2 + sin4n2)) + 
       6*exp(4*mu1 + 7*mu2)*(2*(bc - 3*c1 + 3*bc*mu0 - 3*bc*mu2 + xsc2*nu0 - xsc2*nu1)*(4 + cos4n1)*sin4n1 - 
          (bc - 3*c1 + 3*bc*mu0 - 3*bc*mu2 + xsc2*nu0 - xsc2*nu1)*(8*sin4n2 + sin4n2)) + 
       2*exp(4*mu2)*(-3*(bc*(-1 + 4*mu0 - 4*mu1) - 4*(c1 + 2*c2*(-nu0 + nu1)))*(9*sinn1 + sin3n1) + 
          27*(bc*(-1 + 4*mu0 - 4*mu1) - 4*(c1 + 2*c2*(-nu0 + nu2)))*sinn2 + 3*(bc*(-1 + 4*mu0 - 4*mu1) - 4*(c1 + 2*c2*(-nu0 + nu2)))*sin3n2 + 
          3*exp8m1*(-((bc - 4*c1 + 4*bc*mu0 - 4*bc*mu1 + 8*c2*nu0 - 8*c2*nu1)*(9*sinn1 + sin3n1)) + 
             9*(bc - 4*c1 + 4*bc*mu0 - 4*bc*mu1 + 8*c2*nu0 - 8*c2*nu1)*sinn2 + (bc - 4*c1 + 4*bc*mu0 - 4*bc*mu1 + 8*c2*nu0 - 8*c2*nu1)*sin3n2) + 
          24*exp2m1*((bc + 2*c1 - 2*bc*mu0 + 2*bc*mu1 - 4*c2*nu0 + 4*c2*nu1)*(-63*sinn1 + sin3n1) + 
             63*(bc + 2*c1 - 2*bc*mu0 + 2*bc*mu1 - 4*c2*nu0 + 4*c2*nu1)*sinn2 - (bc + 2*c1 - 2*bc*mu0 + 2*bc*mu1 - 4*c2*nu0 + 4*c2*nu1)*sin3n2) + 
          24*exp6m1*(63*(bc - 2*c1 + 2*bc*mu0 - 2*bc*mu1 + 4*c2*nu0 - 4*c2*nu1)*sinn1 - (bc - 2*c1 + 2*bc*mu0 - 2*bc*mu1 + 4*c2*nu0 - 4*c2*nu1)*sin3n1 + 
             (bc - 2*c1 + 2*bc*mu0 - 2*bc*mu1 + 4*c2*nu0 - 4*c2*nu1)*(-63*sinn2 + sin3n2)) + 
          27*exp5m1*((c1 + bc*(-1 - mu0 + mu1) + 2*c2*(-nu0 + nu1))*(-56*sin4n1 + sin4n1) + 56*(c1 + bc*(-1 - mu0 + mu1) + 2*c2*(-nu0 + nu2))*sin4n2 + 
             (-c1 + bc*(1 + mu0 - mu1) + 2*c2*(nu0 - nu2))*sin4n2) + 
          3*exp(mu1)*((bc + 3*c1 - 3*bc*mu0 + 3*bc*mu1 - xsc2*nu0 + xsc2*nu1)*(8*sin4n1 + sin4n1) - 
             8*(bc + 3*c1 - 3*bc*mu0 + 3*bc*mu1 - xsc2*nu0 + xsc2*nu1)*sin4n2 - (bc + 3*c1 - 3*bc*mu0 + 3*bc*mu1 - xsc2*nu0 + xsc2*nu1)*sin4n2) + 
          27*exp3m1*((c1 + bc*(1 - mu0 + mu1) + 2*c2*(-nu0 + nu1))*(-56*sin4n1 + sin4n1) + 56*(c1 + bc*(1 - mu0 + mu1) + 2*c2*(-nu0 + nu2))*sin4n2 - 
             (c1 + bc*(1 - mu0 + mu1) + 2*c2*(-nu0 + nu2))*sin4n2) + 
          3*exp7m1*(-2*(bc - 3*c1 + 3*bc*mu0 - 3*bc*mu1 + xsc2*nu0 - xsc2*nu1)*(4 + cos4n1)*sin4n1 + 
             (bc - 3*c1 + 3*bc*mu0 - 3*bc*mu1 + xsc2*nu0 - xsc2*nu1)*(8*sin4n2 + sin4n2)) + 
          exp4m1*((-16*(18*a0 + 2*b2*(2 + 9*pm012) - 9*(c2*(1 - 2*pn012) + 2*c1*(nu0 - nu1) + 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1))))*
                 cos4n1 - (72*a0 + 8*b2*(2 + 9*pm012) - 9*(c2*(1 - 8*pn012) + 8*c1*(nu0 - nu1) + 8*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1))))*
                 cos4n1 + 16*(18*a0 + 2*b2*(2 + 9*pm012) - 9*(c2*(1 - 2*pn022) + 2*c1*(nu0 - nu2) + 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2))))*
                 cos4n2 + (72*a0 + 8*b2*(2 + 9*pm012) - 9*(c2*(1 - 8*pn022) + 8*c1*(nu0 - nu2) + 8*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2))))*
                 cos4n2)*cosh3m1 + (-27*(8*a0 + b2*(1 + 8*pm012) + 8*c2*(-2 + pn012) + 8*c1*(-nu0 + nu1) - 
                   8*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 - (72*a0 + nnb2*(1 + 8*pm012) + 
                   8*(c2*(-2 + 9*pn012) + nnc1*(-nu0 + nu1) - 9*(mu0 - mu1)*(b1 - bc*nu0 + bc*nu1)))*cos3n1 + 
                27*(8*a0 + b2*(1 + 8*pm012) + 8*c2*(-2 + pn022) + 8*c1*(-nu0 + nu2) - 8*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 + 
                (72*a0 + nnb2*(1 + 8*pm012) + 8*(c2*(-2 + 9*pn022) + nnc1*(-nu0 + nu2) - 9*(mu0 - mu1)*(b1 - bc*nu0 + bc*nu1)))*cos3n2)*
              cosh4m1 + 6*(-32*(b1 - 2*b2*mu0 + 2*b2*mu1 - bc*nu0 + bc*nu1)*cos4n2*(-31 + cosh2m1)*sinhm1 - 
                8*(b1 - 2*b2*mu0 + 2*b2*mu1 - bc*nu0 + bc*nu1)*cos4n2*(5 + cosh2m1)*sinhm1 - 
                9*(b1 - 2*b2*mu0 + 2*b2*mu1 - bc*nu0 + bc*nu1)*cosn2*(-56*sinh2m1 + sinh4m1) - 
                3*(b1 - 2*b2*mu0 + 2*b2*mu1 - bc*nu0 + bc*nu1)*cos3n2*(8*sinh2m1 + sinh4m1) - 
                (b1 - 2*b2*mu0 + 2*b2*mu1 - bc*nu0 + bc*nu1)*(8*(-4*cos4n1*(-31 + cosh2m1) - cos4n1*(5 + cosh2m1))*sinhm1 - 
                   24*(-21*cosn1 + cos3n1)*sinh2m1 - 12*cosn13*sinh4m1)))))/exp(4*(mu1 + mu2));

  double val2 = 0.;

  return -(Z1*val1 + Z2*val2);
}


//CPMZ marked: nu2 bug in convert function
double second_order_fordV(double Z1, double Z2, double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc)
{
 //zc determines whether center 1 or center 2 has the z ftn
  //b1 = b2 = 0.;
  //c1 = c2 = 0.;
  //bc = 0.;

  double f = 1./8.;
  double pa2 = pow(a,2.);
  double fp = f*pa2;

  double pm02 = pow(mu0,2.);
  double pm12 = pow(mu1,2.);
  double pm22 = pow(mu2,2.);
  double pn02 = pow(nu0,2.);
  double pn12 = pow(nu1,2.);
  double pn22 = pow(nu2,2.);

  double pm012 = pow(mu0-mu1,2.);
  double pm022 = pow(mu0-mu2,2.);
  double pn012 = pow(nu0-nu1,2.);
  double pn022 = pow(nu0-nu2,2.);

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double cosh2m1 = cosh(2.*mu1);
  double cosh2m2 = cosh(2.*mu2);
  double sinh2m1 = sinh(2.*mu1);
  double sinh2m2 = sinh(2.*mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosh4m1 = cosh(4.*mu1);
  double cosh4m2 = cosh(4.*mu2);
  double sinh4m1 = sinh(4.*mu1);
  double sinh4m2 = sinh(4.*mu2);

  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cos2n1 = cos(2.*nu1);
  double cos2n2 = cos(2.*nu2);
  double sin2n1 = sin(2.*nu1);
  double sin2n2 = sin(2.*nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);
  double sin3n1 = sin(3.*nu1);
  double sin3n2 = sin(3.*nu2);
  double cos4n1 = cos(4.*nu1);
  double cos4n2 = cos(4.*nu2);
  double sin4n1 = sin(4.*nu1);
  double sin4n2 = sin(4.*nu2);

  double val1 = fp*(-(cosh2m1*((2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn012) + 2*c1*(-nu0 + nu1) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 - 
            (2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 - 
            2*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sinn1 + 2*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sinn2)) + 
       cosh2m2*((2*a0 + b2*(1 + 2*pm022) + 2*c2*(-2 + pn012) + 2*c1*(-nu0 + nu1) - 2*(mu0 - mu2)*(b1 + bc*(-nu0 + nu1)))*cosn1 - 
          (2*a0 + b2*(1 + 2*pm022) + 2*c2*(-2 + pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu2)*(b1 + bc*(-nu0 + nu2)))*cosn2 - 
          2*(c1 - bc*mu0 + bc*mu2 - 2*c2*nu0 + 2*c2*nu1)*sinn1 + 2*(c1 - bc*mu0 + bc*mu2 - 2*c2*nu0 + 2*c2*nu1)*sinn2) - 
       coshm1*((2*a0 + 2*b2*(2 + pm012) + c2*(-1 + 2*pn012) + 2*c1*(-nu0 + nu1) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos4n1 + 
          (-2*a0 - 2*b2*(2 + pm012) + c2*(1 - 2*pn022) + 2*c1*(nu0 - nu2) + 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2 - 
          (c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n1 + (c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n2) + 
       coshm2*((2*a0 + 2*b2*(2 + pm022) + 2*b1*(-mu0 + mu2) + c2*(-1 + 2*pn012) + 2*(c1 + bc*(-mu0 + mu2))*(-nu0 + nu1))*cos4n1 + 
          (-2*a0 - 2*b2*(2 + pm022) + c2*(1 - 2*pn022) + 2*c1*(nu0 - nu2) + 2*(mu0 - mu2)*(b1 - bc*nu0 + bc*nu1))*cos4n2 - 
          (c1 + bc*(-mu0 + mu2) + 2*c2*(-nu0 + nu1))*sin4n1 + (c1 + bc*(-mu0 + mu2) + 2*c2*(-nu0 + nu2))*sin4n2) + 2*b1*cos4n1*sinhm1 - 
       4*b2*mu0*cos4n1*sinhm1 + 4*b2*mu1*cos4n1*sinhm1 - 2*bc*nu0*cos4n1*sinhm1 + 2*bc*nu1*cos4n1*sinhm1 - 2*b1*cos4n2*sinhm1 + 
       4*b2*mu0*cos4n2*sinhm1 - 4*b2*mu1*cos4n2*sinhm1 + 2*bc*nu0*cos4n2*sinhm1 - 2*bc*nu1*cos4n2*sinhm1 - bc*sin4n1*sinhm1 + 
       bc*sin4n2*sinhm1 + b1*cosn1*sinh2m1 - 2*b2*mu0*cosn1*sinh2m1 + 2*b2*mu1*cosn1*sinh2m1 - bc*nu0*cosn1*sinh2m1 + 
       bc*nu1*cosn1*sinh2m1 - b1*cosn2*sinh2m1 + 2*b2*mu0*cosn2*sinh2m1 - 2*b2*mu1*cosn2*sinh2m1 + bc*nu0*cosn2*sinh2m1 - 
       bc*nu1*cosn2*sinh2m1 - bc*sinn1*sinh2m1 + bc*sinn2*sinh2m1 - 2*b1*cos4n1*sinhm2 + 4*b2*mu0*cos4n1*sinhm2 - 
       4*b2*mu2*cos4n1*sinhm2 + 2*bc*nu0*cos4n1*sinhm2 - 2*bc*nu1*cos4n1*sinhm2 + 2*b1*cos4n2*sinhm2 - 4*b2*mu0*cos4n2*sinhm2 + 
       4*b2*mu2*cos4n2*sinhm2 - 2*bc*nu0*cos4n2*sinhm2 + 2*bc*nu1*cos4n2*sinhm2 + bc*sin4n1*sinhm2 - bc*sin4n2*sinhm2 - 
       b1*cosn1*sinh2m2 + 2*b2*mu0*cosn1*sinh2m2 - 2*b2*mu2*cosn1*sinh2m2 + bc*nu0*cosn1*sinh2m2 - bc*nu1*cosn1*sinh2m2 + 
       b1*cosn2*sinh2m2 - 2*b2*mu0*cosn2*sinh2m2 + 2*b2*mu2*cosn2*sinh2m2 - bc*nu0*cosn2*sinh2m2 + bc*nu1*cosn2*sinh2m2 + 
       bc*sinn1*sinh2m2 - bc*sinn2*sinh2m2);

  double val2 = fp*(-(cosh2m1*((2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn012) + 2*c1*(-nu0 + nu1) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 - 
            (2*a0 + b2*(1 + 2*pm012) + 2*c2*(-2 + pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 - 
            2*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sinn1 + 2*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sinn2)) + 
       cosh2m2*((2*a0 + b2*(1 + 2*pm022) + 2*c2*(-2 + pn012) + 2*c1*(-nu0 + nu1) - 2*(mu0 - mu2)*(b1 + bc*(-nu0 + nu1)))*cosn1 - 
          (2*a0 + b2*(1 + 2*pm022) + 2*c2*(-2 + pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu2)*(b1 + bc*(-nu0 + nu2)))*cosn2 - 
          2*(c1 - bc*mu0 + bc*mu2 - 2*c2*nu0 + 2*c2*nu1)*sinn1 + 2*(c1 - bc*mu0 + bc*mu2 - 2*c2*nu0 + 2*c2*nu1)*sinn2) - 
       coshm1*((-2*a0 - 2*b2*(2 + pm012) + c2*(1 - 2*pn012) + 2*c1*(nu0 - nu1) + 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos4n1 + 
          (2*a0 + 2*b2*(2 + pm012) + c2*(-1 + 2*pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2 + 
          (c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n1 - (c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n2) + 
       coshm2*((-2*a0 - 2*b2*(2 + pm022) + c2*(1 - 2*pn012) + 2*c1*(nu0 - nu1) + 2*(mu0 - mu2)*(b1 + bc*(-nu0 + nu1)))*cos4n1 + 
          (2*a0 + 2*b2*(2 + pm022) + c2*(-1 + 2*pn022) + 2*c1*(-nu0 + nu2) - 2*(mu0 - mu2)*(b1 + bc*(-nu0 + nu2)))*cos4n2 + 
          (c1 - bc*mu0 + bc*mu2 - 2*c2*nu0 + 2*c2*nu1)*sin4n1 - (c1 - bc*mu0 + bc*mu2 - 2*c2*nu0 + 2*c2*nu1)*sin4n2) - 2*b1*cos4n1*sinhm1 + 
       4*b2*mu0*cos4n1*sinhm1 - 4*b2*mu1*cos4n1*sinhm1 + 2*bc*nu0*cos4n1*sinhm1 - 2*bc*nu1*cos4n1*sinhm1 + 2*b1*cos4n2*sinhm1 - 
       4*b2*mu0*cos4n2*sinhm1 + 4*b2*mu1*cos4n2*sinhm1 - 2*bc*nu0*cos4n2*sinhm1 + 2*bc*nu1*cos4n2*sinhm1 + bc*sin4n1*sinhm1 - 
       bc*sin4n2*sinhm1 + b1*cosn1*sinh2m1 - 2*b2*mu0*cosn1*sinh2m1 + 2*b2*mu1*cosn1*sinh2m1 - bc*nu0*cosn1*sinh2m1 + 
       bc*nu1*cosn1*sinh2m1 - b1*cosn2*sinh2m1 + 2*b2*mu0*cosn2*sinh2m1 - 2*b2*mu1*cosn2*sinh2m1 + bc*nu0*cosn2*sinh2m1 - 
       bc*nu1*cosn2*sinh2m1 - bc*sinn1*sinh2m1 + bc*sinn2*sinh2m1 + 2*b1*cos4n1*sinhm2 - 4*b2*mu0*cos4n1*sinhm2 + 
       4*b2*mu2*cos4n1*sinhm2 - 2*bc*nu0*cos4n1*sinhm2 + 2*bc*nu1*cos4n1*sinhm2 - 2*b1*cos4n2*sinhm2 + 4*b2*mu0*cos4n2*sinhm2 - 
       4*b2*mu2*cos4n2*sinhm2 + 2*bc*nu0*cos4n2*sinhm2 - 2*bc*nu1*cos4n2*sinhm2 - bc*sin4n1*sinhm2 + bc*sin4n2*sinhm2 - 
       b1*cosn1*sinh2m2 + 2*b2*mu0*cosn1*sinh2m2 - 2*b2*mu2*cosn1*sinh2m2 + bc*nu0*cosn1*sinh2m2 - bc*nu1*cosn1*sinh2m2 + 
       b1*cosn2*sinh2m2 - 2*b2*mu0*cosn2*sinh2m2 + 2*b2*mu2*cosn2*sinh2m2 - bc*nu0*cosn2*sinh2m2 + bc*nu1*cosn2*sinh2m2 + 
       bc*sinn1*sinh2m2 - bc*sinn2*sinh2m2);


  return -(Z1*val1 + Z2*val2);
}

//CPMZ marked: nu2 bug in convert function
double second_order_fzzdV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc)
{
  //b1 = b2 = 0.;
  //c1 = c2 = 0.;
  //bc = 0.;

  double f = 1./216000.;
  double pa5 = pow(a,5.);
  double fp = f*pa5;

  double pm02 = pow(mu0,2.);
  double pm12 = pow(mu1,2.);
  double pm22 = pow(mu2,2.);
  double pn02 = pow(nu0,2.);
  double pn12 = pow(nu1,2.);
  double pn22 = pow(nu2,2.);

  double pm012 = pow(mu0-mu1,2.);
  double pm022 = pow(mu0-mu2,2.);
  double pn012 = pow(nu0-nu1,2.);
  double pn022 = pow(nu0-nu2,2.);

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double cosh2m1 = cosh(2.*mu1);
  double cosh2m2 = cosh(2.*mu2);
  double sinh2m1 = sinh(2.*mu1);
  double sinh2m2 = sinh(2.*mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosh4m1 = cosh(4.*mu1);
  double cosh4m2 = cosh(4.*mu2);
  double sinh4m1 = sinh(4.*mu1);
  double sinh4m2 = sinh(4.*mu2);
  double cosh5m1 = cosh(5.*mu1);
  double cosh5m2 = cosh(5.*mu2);
  double sinh5m1 = sinh(5.*mu1);
  double sinh5m2 = sinh(5.*mu2);
  double coshm13 = coshm1*coshm1*coshm1;
  double coshm23 = coshm2*coshm2*coshm2;

  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cos2n1 = cos(2.*nu1);
  double cos2n2 = cos(2.*nu2);
  double sin2n1 = sin(2.*nu1);
  double sin2n2 = sin(2.*nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);
  double sin3n1 = sin(3.*nu1);
  double sin3n2 = sin(3.*nu2);
  double cos4n1 = cos(4.*nu1);
  double cos4n2 = cos(4.*nu2);
  double sin4n1 = sin(4.*nu1);
  double sin4n2 = sin(4.*nu2);
  double cos5n1 = cos(5.*nu1);
  double cos5n2 = cos(5.*nu2);
  double sin5n1 = sin(5.*nu1);
  double sin5n2 = sin(5.*nu2);

  double ttbc = 225.*bc;
  double ttb1 = 225.*b1;
  double ttb2 = 225.*b2;
  double ttc1 = 225.*c1;
  double ttc2 = 225.*c2;

  double fsbc = 5625.*bc;
  double fsb1 = 5625.*b1;
  double fsb2 = 5625.*b2;
  double fsc1 = 5625.*c1;
  double fsc2 = 5625.*c2;

  double val = fp*(-675*a0*cosn1*cosh5m1 - 54*b2*cosn1*cosh5m1 + 1350*c2*cosn1*cosh5m1 + 675*b1*mu0*cosn1*cosh5m1 - 
       675*b2*pm02*cosn1*cosh5m1 - 675*b1*mu1*cosn1*cosh5m1 + 1350*b2*mu0*mu1*cosn1*cosh5m1 - 675*b2*pm12*cosn1*cosh5m1 + 
       675*c1*nu0*cosn1*cosh5m1 - 675*bc*mu0*nu0*cosn1*cosh5m1 + 675*bc*mu1*nu0*cosn1*cosh5m1 - 675*c2*pn02*cosn1*cosh5m1 - 
       675*c1*nu1*cosn1*cosh5m1 + 675*bc*mu0*nu1*cosn1*cosh5m1 - 675*bc*mu1*nu1*cosn1*cosh5m1 + 1350*c2*nu0*nu1*cosn1*cosh5m1 - 
       675*c2*pn12*cosn1*cosh5m1 - 225*a0*cos3n1*cosh5m1 - 18*b2*cos3n1*cosh5m1 + 50*c2*cos3n1*cosh5m1 + 
       ttb1*mu0*cos3n1*cosh5m1 - ttb2*pm02*cos3n1*cosh5m1 - ttb1*mu1*cos3n1*cosh5m1 + 450*b2*mu0*mu1*cos3n1*cosh5m1 - 
       ttb2*pm12*cos3n1*cosh5m1 + ttc1*nu0*cos3n1*cosh5m1 - ttbc*mu0*nu0*cos3n1*cosh5m1 + ttbc*mu1*nu0*cos3n1*cosh5m1 - 
       ttc2*pn02*cos3n1*cosh5m1 - ttc1*nu1*cos3n1*cosh5m1 + ttbc*mu0*nu1*cos3n1*cosh5m1 - ttbc*mu1*nu1*cos3n1*cosh5m1 + 
       450*c2*nu0*nu1*cos3n1*cosh5m1 - ttc2*pn12*cos3n1*cosh5m1 + 675*a0*cosn2*cosh5m1 + 54*b2*cosn2*cosh5m1 - 
       1350*c2*cosn2*cosh5m1 - 675*b1*mu0*cosn2*cosh5m1 + 675*b2*pm02*cosn2*cosh5m1 + 675*b1*mu1*cosn2*cosh5m1 - 
       1350*b2*mu0*mu1*cosn2*cosh5m1 + 675*b2*pm12*cosn2*cosh5m1 - 675*c1*nu0*cosn2*cosh5m1 + 675*bc*mu0*nu0*cosn2*cosh5m1 - 
       675*bc*mu1*nu0*cosn2*cosh5m1 + 675*c2*pn02*cosn2*cosh5m1 + 675*c1*nu1*cosn2*cosh5m1 - 675*bc*mu0*nu1*cosn2*cosh5m1 + 
       675*bc*mu1*nu1*cosn2*cosh5m1 - 1350*c2*nu0*nu1*cosn2*cosh5m1 + 675*c2*pn22*cosn2*cosh5m1 + 225*a0*cos3n2*cosh5m1 + 
       18*b2*cos3n2*cosh5m1 - 50*c2*cos3n2*cosh5m1 - ttb1*mu0*cos3n2*cosh5m1 + ttb2*pm02*cos3n2*cosh5m1 + 
       ttb1*mu1*cos3n2*cosh5m1 - 450*b2*mu0*mu1*cos3n2*cosh5m1 + ttb2*pm12*cos3n2*cosh5m1 - ttc1*nu0*cos3n2*cosh5m1 + 
       ttbc*mu0*nu0*cos3n2*cosh5m1 - ttbc*mu1*nu0*cos3n2*cosh5m1 + ttc2*pn02*cos3n2*cosh5m1 + ttc1*nu1*cos3n2*cosh5m1 - 
       ttbc*mu0*nu1*cos3n2*cosh5m1 + ttbc*mu1*nu1*cos3n2*cosh5m1 - 450*c2*nu0*nu1*cos3n2*cosh5m1 + ttc2*pn22*cos3n2*cosh5m1 + 
       16875*a0*cos3n1*coshm2 + 33750*b2*cos3n1*coshm2 - 3750*c2*cos3n1*coshm2 - 16875*b1*mu0*cos3n1*coshm2 + 
       16875*b2*pm02*cos3n1*coshm2 + 16875*b1*mu2*cos3n1*coshm2 - 33750*b2*mu0*mu2*cos3n1*coshm2 + 16875*b2*pm22*cos3n1*coshm2 - 
       16875*c1*nu0*cos3n1*coshm2 + 16875*bc*mu0*nu0*cos3n1*coshm2 - 16875*bc*mu2*nu0*cos3n1*coshm2 + 16875*c2*pn02*cos3n1*coshm2 + 
       16875*c1*nu1*cos3n1*coshm2 - 16875*bc*mu0*nu1*cos3n1*coshm2 + 16875*bc*mu2*nu1*cos3n1*coshm2 - 33750*c2*nu0*nu1*cos3n1*coshm2 + 
       16875*c2*pn12*cos3n1*coshm2 - 1350*b2*cos5n1*coshm2 - 16875*a0*cos3n2*coshm2 - 33750*b2*cos3n2*coshm2 + 
       3750*c2*cos3n2*coshm2 + 16875*b1*mu0*cos3n2*coshm2 - 16875*b2*pm02*cos3n2*coshm2 - 16875*b1*mu2*cos3n2*coshm2 + 
       33750*b2*mu0*mu2*cos3n2*coshm2 - 16875*b2*pm22*cos3n2*coshm2 + 16875*c1*nu0*cos3n2*coshm2 - 16875*bc*mu0*nu0*cos3n2*coshm2 + 
       16875*bc*mu2*nu0*cos3n2*coshm2 - 16875*c2*pn02*cos3n2*coshm2 - 16875*c1*nu1*cos3n2*coshm2 + 16875*bc*mu0*nu1*cos3n2*coshm2 - 
       16875*bc*mu2*nu1*cos3n2*coshm2 + 33750*c2*nu0*nu1*cos3n2*coshm2 - 16875*c2*pn22*cos3n2*coshm2 + 1350*b2*cos5n2*coshm2 - 
       900*a0*cos5n1*coshm23 + 72*c2*cos5n1*coshm23 + 900*b1*mu0*cos5n1*coshm23 - 
       900*b2*pm02*cos5n1*coshm23 - 900*b1*mu2*cos5n1*coshm23 + 1800*b2*mu0*mu2*cos5n1*coshm23 - 
       900*b2*pm22*cos5n1*coshm23 + 900*c1*nu0*cos5n1*coshm23 - 900*bc*mu0*nu0*cos5n1*coshm23 + 
       900*bc*mu2*nu0*cos5n1*coshm23 - 900*c2*pn02*cos5n1*coshm23 - 900*c1*nu1*cos5n1*coshm23 + 
       900*bc*mu0*nu1*cos5n1*coshm23 - 900*bc*mu2*nu1*cos5n1*coshm23 + 1800*c2*nu0*nu1*cos5n1*coshm23 - 
       900*c2*pn12*cos5n1*coshm23 + 900*a0*cos5n2*coshm23 - 72*c2*cos5n2*coshm23 - 
       900*b1*mu0*cos5n2*coshm23 + 900*b2*pm02*cos5n2*coshm23 + 900*b1*mu2*cos5n2*coshm23 - 
       1800*b2*mu0*mu2*cos5n2*coshm23 + 900*b2*pm22*cos5n2*coshm23 - 900*c1*nu0*cos5n2*coshm23 + 
       900*bc*mu0*nu0*cos5n2*coshm23 - 900*bc*mu2*nu0*cos5n2*coshm23 + 900*c2*pn02*cos5n2*coshm23 + 
       900*c1*nu1*cos5n2*coshm23 - 900*bc*mu0*nu1*cos5n2*coshm23 + 900*bc*mu2*nu1*cos5n2*coshm23 - 
       1800*c2*nu0*nu1*cos5n2*coshm23 + 900*c2*pn22*cos5n2*coshm23 - 16875*a0*cosn1*cosh3m2 - 3750*b2*cosn1*cosh3m2 + 
       33750*c2*cosn1*cosh3m2 + 16875*b1*mu0*cosn1*cosh3m2 - 16875*b2*pm02*cosn1*cosh3m2 - 16875*b1*mu2*cosn1*cosh3m2 + 
       33750*b2*mu0*mu2*cosn1*cosh3m2 - 16875*b2*pm22*cosn1*cosh3m2 + 16875*c1*nu0*cosn1*cosh3m2 - 16875*bc*mu0*nu0*cosn1*cosh3m2 + 
       16875*bc*mu2*nu0*cosn1*cosh3m2 - 16875*c2*pn02*cosn1*cosh3m2 - 16875*c1*nu1*cosn1*cosh3m2 + 16875*bc*mu0*nu1*cosn1*cosh3m2 - 
       16875*bc*mu2*nu1*cosn1*cosh3m2 + 33750*c2*nu0*nu1*cosn1*cosh3m2 - 16875*c2*pn12*cosn1*cosh3m2 - 50*b2*cos5n1*cosh3m2 + 
       16875*a0*cosn2*cosh3m2 + 3750*b2*cosn2*cosh3m2 - 33750*c2*cosn2*cosh3m2 - 16875*b1*mu0*cosn2*cosh3m2 + 
       16875*b2*pm02*cosn2*cosh3m2 + 16875*b1*mu2*cosn2*cosh3m2 - 33750*b2*mu0*mu2*cosn2*cosh3m2 + 16875*b2*pm22*cosn2*cosh3m2 - 
       16875*c1*nu0*cosn2*cosh3m2 + 16875*bc*mu0*nu0*cosn2*cosh3m2 - 16875*bc*mu2*nu0*cosn2*cosh3m2 + 16875*c2*pn02*cosn2*cosh3m2 + 
       16875*c1*nu1*cosn2*cosh3m2 - 16875*bc*mu0*nu1*cosn2*cosh3m2 + 16875*bc*mu2*nu1*cosn2*cosh3m2 - 33750*c2*nu0*nu1*cosn2*cosh3m2 + 
       16875*c2*pn22*cosn2*cosh3m2 + 50*b2*cos5n2*cosh3m2 + 675*a0*cosn1*cosh5m2 + 54*b2*cosn1*cosh5m2 - 
       1350*c2*cosn1*cosh5m2 - 675*b1*mu0*cosn1*cosh5m2 + 675*b2*pm02*cosn1*cosh5m2 + 675*b1*mu2*cosn1*cosh5m2 - 
       1350*b2*mu0*mu2*cosn1*cosh5m2 + 675*b2*pm22*cosn1*cosh5m2 - 675*c1*nu0*cosn1*cosh5m2 + 675*bc*mu0*nu0*cosn1*cosh5m2 - 
       675*bc*mu2*nu0*cosn1*cosh5m2 + 675*c2*pn02*cosn1*cosh5m2 + 675*c1*nu1*cosn1*cosh5m2 - 675*bc*mu0*nu1*cosn1*cosh5m2 + 
       675*bc*mu2*nu1*cosn1*cosh5m2 - 1350*c2*nu0*nu1*cosn1*cosh5m2 + 675*c2*pn12*cosn1*cosh5m2 + 225*a0*cos3n1*cosh5m2 + 
       18*b2*cos3n1*cosh5m2 - 50*c2*cos3n1*cosh5m2 - ttb1*mu0*cos3n1*cosh5m2 + ttb2*pm02*cos3n1*cosh5m2 + 
       ttb1*mu2*cos3n1*cosh5m2 - 450*b2*mu0*mu2*cos3n1*cosh5m2 + ttb2*pm22*cos3n1*cosh5m2 - ttc1*nu0*cos3n1*cosh5m2 + 
       ttbc*mu0*nu0*cos3n1*cosh5m2 - ttbc*mu2*nu0*cos3n1*cosh5m2 + ttc2*pn02*cos3n1*cosh5m2 + ttc1*nu1*cos3n1*cosh5m2 - 
       ttbc*mu0*nu1*cos3n1*cosh5m2 + ttbc*mu2*nu1*cos3n1*cosh5m2 - 450*c2*nu0*nu1*cos3n1*cosh5m2 + ttc2*pn12*cos3n1*cosh5m2 - 
       675*a0*cosn2*cosh5m2 - 54*b2*cosn2*cosh5m2 + 1350*c2*cosn2*cosh5m2 + 675*b1*mu0*cosn2*cosh5m2 - 
       675*b2*pm02*cosn2*cosh5m2 - 675*b1*mu2*cosn2*cosh5m2 + 1350*b2*mu0*mu2*cosn2*cosh5m2 - 675*b2*pm22*cosn2*cosh5m2 + 
       675*c1*nu0*cosn2*cosh5m2 - 675*bc*mu0*nu0*cosn2*cosh5m2 + 675*bc*mu2*nu0*cosn2*cosh5m2 - 675*c2*pn02*cosn2*cosh5m2 - 
       675*c1*nu1*cosn2*cosh5m2 + 675*bc*mu0*nu1*cosn2*cosh5m2 - 675*bc*mu2*nu1*cosn2*cosh5m2 + 1350*c2*nu0*nu1*cosn2*cosh5m2 - 
       675*c2*pn22*cosn2*cosh5m2 - 225*a0*cos3n2*cosh5m2 - 18*b2*cos3n2*cosh5m2 + 50*c2*cos3n2*cosh5m2 + 
       ttb1*mu0*cos3n2*cosh5m2 - ttb2*pm02*cos3n2*cosh5m2 - ttb1*mu2*cos3n2*cosh5m2 + 450*b2*mu0*mu2*cos3n2*cosh5m2 - 
       ttb2*pm22*cos3n2*cosh5m2 + ttc1*nu0*cos3n2*cosh5m2 - ttbc*mu0*nu0*cos3n2*cosh5m2 + ttbc*mu2*nu0*cos3n2*cosh5m2 - 
       ttc2*pn02*cos3n2*cosh5m2 - ttc1*nu1*cos3n2*cosh5m2 + ttbc*mu0*nu1*cos3n2*cosh5m2 - ttbc*mu2*nu1*cos3n2*cosh5m2 + 
       450*c2*nu0*nu1*cos3n2*cosh5m2 - ttc2*pn22*cos3n2*cosh5m2 + 675*c1*cosh5m1*sinn1 - 675*bc*mu0*cosh5m1*sinn1 + 
       675*bc*mu1*cosh5m1*sinn1 - 1350*c2*nu0*cosh5m1*sinn1 + 1350*c2*nu1*cosh5m1*sinn1 + 16875*c1*cosh3m2*sinn1 - 
       16875*bc*mu0*cosh3m2*sinn1 + 16875*bc*mu2*cosh3m2*sinn1 - 33750*c2*nu0*cosh3m2*sinn1 + 33750*c2*nu1*cosh3m2*sinn1 - 
       675*c1*cosh5m2*sinn1 + 675*bc*mu0*cosh5m2*sinn1 - 675*bc*mu2*cosh5m2*sinn1 + 1350*c2*nu0*cosh5m2*sinn1 - 
       1350*c2*nu1*cosh5m2*sinn1 + 75*c1*cosh5m1*sin3n1 - 75*bc*mu0*cosh5m1*sin3n1 + 75*bc*mu1*cosh5m1*sin3n1 - 
       150*c2*nu0*cosh5m1*sin3n1 + 150*c2*nu1*cosh5m1*sin3n1 - fsc1*coshm2*sin3n1 + fsbc*mu0*coshm2*sin3n1 - 
       fsbc*mu2*coshm2*sin3n1 + 11250*c2*nu0*coshm2*sin3n1 - 11250*c2*nu1*coshm2*sin3n1 - 75*c1*cosh5m2*sin3n1 + 
       75*bc*mu0*cosh5m2*sin3n1 - 75*bc*mu2*cosh5m2*sin3n1 + 150*c2*nu0*cosh5m2*sin3n1 - 150*c2*nu1*cosh5m2*sin3n1 + 
       180*c1*coshm23*sin5n1 - 180*bc*mu0*coshm23*sin5n1 + 180*bc*mu2*coshm23*sin5n1 - 
       360*c2*nu0*coshm23*sin5n1 + 360*c2*nu1*coshm23*sin5n1 - 675*c1*cosh5m1*sinn2 + 675*bc*mu0*cosh5m1*sinn2 - 
       675*bc*mu1*cosh5m1*sinn2 + 1350*c2*nu0*cosh5m1*sinn2 - 1350*c2*nu1*cosh5m1*sinn2 - 16875*c1*cosh3m2*sinn2 + 
       16875*bc*mu0*cosh3m2*sinn2 - 16875*bc*mu2*cosh3m2*sinn2 + 33750*c2*nu0*cosh3m2*sinn2 - 33750*c2*nu1*cosh3m2*sinn2 + 
       675*c1*cosh5m2*sinn2 - 675*bc*mu0*cosh5m2*sinn2 + 675*bc*mu2*cosh5m2*sinn2 - 1350*c2*nu0*cosh5m2*sinn2 + 
       1350*c2*nu1*cosh5m2*sinn2 + 25*cosh3m1*(75*(9*a0 + b2*(2 + 9*pm012) + 9*c2*(-2 + pn012) + 9*c1*(-nu0 + nu1) - 
             9*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 + 2*b2*cos5n1 - 
          75*(9*a0 + b2*(2 + 9*pm012) + 9*c2*(-2 + pn022) + 9*c1*(-nu0 + nu2) - 9*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 - 
          2*b2*cos5n2 - 675*(c1 + bc*(-mu0 + mu1) + 2*c2*(-nu0 + nu1))*sinn1 + 675*(c1 + bc*(-mu0 + mu1) + 2*c2*(-nu0 + nu2))*sinn2) - 
       75*c1*cosh5m1*sin3n2 + 75*bc*mu0*cosh5m1*sin3n2 - 75*bc*mu1*cosh5m1*sin3n2 + 150*c2*nu0*cosh5m1*sin3n2 - 
       150*c2*nu1*cosh5m1*sin3n2 + fsc1*coshm2*sin3n2 - fsbc*mu0*coshm2*sin3n2 + fsbc*mu2*coshm2*sin3n2 - 
       11250*c2*nu0*coshm2*sin3n2 + 11250*c2*nu1*coshm2*sin3n2 + 75*c1*cosh5m2*sin3n2 - 75*bc*mu0*cosh5m2*sin3n2 + 
       75*bc*mu2*cosh5m2*sin3n2 - 150*c2*nu0*cosh5m2*sin3n2 + 150*c2*nu1*cosh5m2*sin3n2 + 
       75*coshm1*(-25*(9*a0 + 9*b2*(2 + pm012) + c2*(-2 + 9*pn012) + 9*c1*(-nu0 + nu1) - 9*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos3n1 + 
          18*b2*cos5n1 + 25*(9*a0 + 9*b2*(2 + pm012) + c2*(-2 + 9*pn022) + 9*c1*(-nu0 + nu2) - 9*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*
           cos3n2 - 18*b2*cos5n2 + 75*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin3n1 - 75*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin3n2) - 
       180*c1*coshm23*sin5n2 + 180*bc*mu0*coshm23*sin5n2 - 180*bc*mu2*coshm23*sin5n2 + 
       360*c2*nu0*coshm23*sin5n2 - 360*c2*nu1*coshm23*sin5n2 + 
       36*coshm13*((25*a0 + c2*(-2 + 25*pn012) - 25*(c1*(nu0 - nu1) + (mu0 - mu1)*(b1 - b2*mu0 + b2*mu1 - bc*nu0 + bc*nu1)))*cos5n1 - 
          (25*a0 + c2*(-2 + 25*pn022) - 25*(c1*(nu0 - nu2) + (mu0 - mu1)*(b1 - b2*mu0 + b2*mu1 - bc*nu0 + bc*nu1)))*cos5n2 - 
          5*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin5n1 + 5*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin5n2) + 16875*b1*cos3n1*sinhm1 - 
       33750*b2*mu0*cos3n1*sinhm1 + 33750*b2*mu1*cos3n1*sinhm1 - 16875*bc*nu0*cos3n1*sinhm1 + 16875*bc*nu1*cos3n1*sinhm1 - 
       675*b1*cos5n1*sinhm1 + 1350*b2*mu0*cos5n1*sinhm1 - 1350*b2*mu1*cos5n1*sinhm1 + 675*bc*nu0*cos5n1*sinhm1 - 
       675*bc*nu1*cos5n1*sinhm1 - 16875*b1*cos3n2*sinhm1 + 33750*b2*mu0*cos3n2*sinhm1 - 33750*b2*mu1*cos3n2*sinhm1 + 
       16875*bc*nu0*cos3n2*sinhm1 - 16875*bc*nu1*cos3n2*sinhm1 + 675*b1*cos5n2*sinhm1 - 1350*b2*mu0*cos5n2*sinhm1 + 
       1350*b2*mu1*cos5n2*sinhm1 - 675*bc*nu0*cos5n2*sinhm1 + 675*bc*nu1*cos5n2*sinhm1 - fsbc*sin3n1*sinhm1 + 135*bc*sin5n1*sinhm1 + 
       fsbc*sin3n2*sinhm1 - 135*bc*sin5n2*sinhm1 - fsb1*cosn1*sinh3m1 + 11250*b2*mu0*cosn1*sinh3m1 - 11250*b2*mu1*cosn1*sinh3m1 + 
       fsbc*nu0*cosn1*sinh3m1 - fsbc*nu1*cosn1*sinh3m1 - 75*b1*cos5n1*sinh3m1 + 150*b2*mu0*cos5n1*sinh3m1 - 
       150*b2*mu1*cos5n1*sinh3m1 + 75*bc*nu0*cos5n1*sinh3m1 - 75*bc*nu1*cos5n1*sinh3m1 + fsb1*cosn2*sinh3m1 - 
       11250*b2*mu0*cosn2*sinh3m1 + 11250*b2*mu1*cosn2*sinh3m1 - fsbc*nu0*cosn2*sinh3m1 + fsbc*nu1*cosn2*sinh3m1 + 
       75*b1*cos5n2*sinh3m1 - 150*b2*mu0*cos5n2*sinh3m1 + 150*b2*mu1*cos5n2*sinh3m1 - 75*bc*nu0*cos5n2*sinh3m1 + 
       75*bc*nu1*cos5n2*sinh3m1 + fsbc*sinn1*sinh3m1 + 15*bc*sin5n1*sinh3m1 - fsbc*sinn2*sinh3m1 - 15*bc*sin5n2*sinh3m1 + 
       135*b1*cosn1*sinh5m1 - 270*b2*mu0*cosn1*sinh5m1 + 270*b2*mu1*cosn1*sinh5m1 - 135*bc*nu0*cosn1*sinh5m1 + 135*bc*nu1*cosn1*sinh5m1 + 
       45*b1*cos3n1*sinh5m1 - 90*b2*mu0*cos3n1*sinh5m1 + 90*b2*mu1*cos3n1*sinh5m1 - 45*bc*nu0*cos3n1*sinh5m1 + 
       45*bc*nu1*cos3n1*sinh5m1 - 135*b1*cosn2*sinh5m1 + 270*b2*mu0*cosn2*sinh5m1 - 270*b2*mu1*cosn2*sinh5m1 + 135*bc*nu0*cosn2*sinh5m1 - 
       135*bc*nu1*cosn2*sinh5m1 - 45*b1*cos3n2*sinh5m1 + 90*b2*mu0*cos3n2*sinh5m1 - 90*b2*mu1*cos3n2*sinh5m1 + 
       45*bc*nu0*cos3n2*sinh5m1 - 45*bc*nu1*cos3n2*sinh5m1 - 135*bc*sinn1*sinh5m1 - 15*bc*sin3n1*sinh5m1 + 135*bc*sinn2*sinh5m1 + 
       15*bc*sin3n2*sinh5m1 - 16875*b1*cos3n1*sinhm2 + 33750*b2*mu0*cos3n1*sinhm2 - 33750*b2*mu2*cos3n1*sinhm2 + 
       16875*bc*nu0*cos3n1*sinhm2 - 16875*bc*nu1*cos3n1*sinhm2 + 675*b1*cos5n1*sinhm2 - 1350*b2*mu0*cos5n1*sinhm2 + 
       1350*b2*mu2*cos5n1*sinhm2 - 675*bc*nu0*cos5n1*sinhm2 + 675*bc*nu1*cos5n1*sinhm2 + 16875*b1*cos3n2*sinhm2 - 
       33750*b2*mu0*cos3n2*sinhm2 + 33750*b2*mu2*cos3n2*sinhm2 - 16875*bc*nu0*cos3n2*sinhm2 + 16875*bc*nu1*cos3n2*sinhm2 - 
       675*b1*cos5n2*sinhm2 + 1350*b2*mu0*cos5n2*sinhm2 - 1350*b2*mu2*cos5n2*sinhm2 + 675*bc*nu0*cos5n2*sinhm2 - 
       675*bc*nu1*cos5n2*sinhm2 + fsbc*sin3n1*sinhm2 - 135*bc*sin5n1*sinhm2 - fsbc*sin3n2*sinhm2 + 135*bc*sin5n2*sinhm2 + 
       fsb1*cosn1*sinh3m2 - 11250*b2*mu0*cosn1*sinh3m2 + 11250*b2*mu2*cosn1*sinh3m2 - fsbc*nu0*cosn1*sinh3m2 + 
       fsbc*nu1*cosn1*sinh3m2 + 75*b1*cos5n1*sinh3m2 - 150*b2*mu0*cos5n1*sinh3m2 + 150*b2*mu2*cos5n1*sinh3m2 - 
       75*bc*nu0*cos5n1*sinh3m2 + 75*bc*nu1*cos5n1*sinh3m2 - fsb1*cosn2*sinh3m2 + 11250*b2*mu0*cosn2*sinh3m2 - 
       11250*b2*mu2*cosn2*sinh3m2 + fsbc*nu0*cosn2*sinh3m2 - fsbc*nu1*cosn2*sinh3m2 - 75*b1*cos5n2*sinh3m2 + 
       150*b2*mu0*cos5n2*sinh3m2 - 150*b2*mu2*cos5n2*sinh3m2 + 75*bc*nu0*cos5n2*sinh3m2 - 75*bc*nu1*cos5n2*sinh3m2 - 
       fsbc*sinn1*sinh3m2 - 15*bc*sin5n1*sinh3m2 + fsbc*sinn2*sinh3m2 + 15*bc*sin5n2*sinh3m2 - 
       15*(9*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cosn1 + 3*b1*cos3n1 - 6*b2*mu0*cos3n1 + 6*b2*mu2*cos3n1 - 3*bc*nu0*cos3n1 + 
          3*bc*nu1*cos3n1 - 9*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cosn2 - 3*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cos3n2 - 9*bc*sinn1 - 
          bc*sin3n1 + 9*bc*sinn2 + bc*sin3n2)*sinh5m2);
;

  return val;
}

//CPMZ marked: nu2 bug in convert function
double second_order_fzdV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double zc, double a0, double b1, double b2, double c1, double c2, double bc)
{
 //zc determines whether center 1 or center 2 has the z ftn
  double f = 1./27648.;
  double pa4 = pow(a,4.);
  double fp = f*pa4;

  double pm02 = pow(mu0,2.);
  double pm12 = pow(mu1,2.);
  double pm22 = pow(mu2,2.);
  double pn02 = pow(nu0,2.);
  double pn12 = pow(nu1,2.);
  double pn22 = pow(nu2,2.);

  double pm012 = pow(mu0-mu1,2.);
  double pm022 = pow(mu0-mu2,2.);
  double pn012 = pow(nu0-nu1,2.);
  double pn022 = pow(nu0-nu2,2.);

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double cosh2m1 = cosh(2.*mu1);
  double cosh2m2 = cosh(2.*mu2);
  double sinh2m1 = sinh(2.*mu1);
  double sinh2m2 = sinh(2.*mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosh4m1 = cosh(4.*mu1);
  double cosh4m2 = cosh(4.*mu2);
  double sinh4m1 = sinh(4.*mu1);
  double sinh4m2 = sinh(4.*mu2);

  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cos2n1 = cos(2.*nu1);
  double cos2n2 = cos(2.*nu2);
  double sin2n1 = sin(2.*nu1);
  double sin2n2 = sin(2.*nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);
  double sin3n1 = sin(3.*nu1);
  double sin3n2 = sin(3.*nu2);
  double cos4n1 = cos(4.*nu1);
  double cos4n2 = cos(4.*nu2);
  double sin4n1 = sin(4.*nu1);
  double sin4n2 = sin(4.*nu2);

  double tsbc = 216.*bc;
  double tsb1 = 216.*b1;
  double tsb2 = 216.*b2;
  double tsc1 = 216.*c1;
  double tsc2 = 216.*c2;

  double ftbc = 432.*bc;
  double ftb1 = 432.*b1;
  double ftb2 = 432.*b2;
  double ftc1 = 432.*c1;
  double ftc2 = 432.*c2;

  double tobc = 2304.*bc;
  double tob1 = 2304.*b1;
  double tob2 = 2304.*b2;
  double toc1 = 2304.*c1;
  double toc2 = 2304.*c2;

  double febc = 4608.*bc;
  double feb1 = 4608.*b1;
  double feb2 = 4608.*b2;
  double fec1 = 4608.*c1;
  double fec2 = 4608.*c2;

  double val = fp*(2304*a0*zc*cosn1*cosh3m1 + 512*b2*zc*cosn1*cosh3m1 - 
       fec2*zc*cosn1*cosh3m1 - tob1*mu0*zc*cosn1*cosh3m1 + 
       tob2*pm02*zc*cosn1*cosh3m1 + tob1*mu1*zc*cosn1*cosh3m1 - 
       feb2*mu0*mu1*zc*cosn1*cosh3m1 + 
       tob2*pm12*zc*cosn1*cosh3m1 - toc1*nu0*zc*cosn1*cosh3m1 + 
       tobc*mu0*nu0*zc*cosn1*cosh3m1 - tobc*mu1*nu0*zc*cosn1*cosh3m1 + 
       toc2*pn02*zc*cosn1*cosh3m1 + toc1*nu1*zc*cosn1*cosh3m1 - 
       tobc*mu0*nu1*zc*cosn1*cosh3m1 + tobc*mu1*nu1*zc*cosn1*cosh3m1 - 
       fec2*nu0*nu1*zc*cosn1*cosh3m1 + 
       toc2*pn12*zc*cosn1*cosh3m1 - 2304*a0*zc*cosn2*cosh3m1 - 
       512*b2*zc*cosn2*cosh3m1 + fec2*zc*cosn2*cosh3m1 + 
       tob1*mu0*zc*cosn2*cosh3m1 - tob2*pm02*zc*cosn2*cosh3m1 - 
       tob1*mu1*zc*cosn2*cosh3m1 + feb2*mu0*mu1*zc*cosn2*cosh3m1 - 
       tob2*pm12*zc*cosn2*cosh3m1 + toc1*nu0*zc*cosn2*cosh3m1 - 
       tobc*mu0*nu0*zc*cosn2*cosh3m1 + tobc*mu1*nu0*zc*cosn2*cosh3m1 - 
       toc2*pn02*zc*cosn2*cosh3m1 - toc1*nu1*zc*cosn2*cosh3m1 + 
       tobc*mu0*nu1*zc*cosn2*cosh3m1 - tobc*mu1*nu1*zc*cosn2*cosh3m1 + 
       fec2*nu0*nu1*zc*cosn2*cosh3m1 - 
       toc2*pn22*zc*cosn2*cosh3m1 - 216*a0*cos4n1*cosh4m1 - 
       27*b2*cos4n1*cosh4m1 + 108*c2*cos4n1*cosh4m1 + 
       tsb1*mu0*cos4n1*cosh4m1 - tsb2*pm02*cos4n1*cosh4m1 - 
       tsb1*mu1*cos4n1*cosh4m1 + ftb2*mu0*mu1*cos4n1*cosh4m1 - 
       tsb2*pm12*cos4n1*cosh4m1 + tsc1*nu0*cos4n1*cosh4m1 - 
       tsbc*mu0*nu0*cos4n1*cosh4m1 + tsbc*mu1*nu0*cos4n1*cosh4m1 - 
       tsc2*pn02*cos4n1*cosh4m1 - tsc1*nu1*cos4n1*cosh4m1 + 
       tsbc*mu0*nu1*cos4n1*cosh4m1 - tsbc*mu1*nu1*cos4n1*cosh4m1 + 
       ftc2*nu0*nu1*cos4n1*cosh4m1 - tsc2*pn12*cos4n1*cosh4m1 + 
       216*a0*cos4n2*cosh4m1 + 27*b2*cos4n2*cosh4m1 - 
       108*c2*cos4n2*cosh4m1 - tsb1*mu0*cos4n2*cosh4m1 + 
       tsb2*pm02*cos4n2*cosh4m1 + tsb1*mu1*cos4n2*cosh4m1 - 
       ftb2*mu0*mu1*cos4n2*cosh4m1 + tsb2*pm12*cos4n2*cosh4m1 - 
       tsc1*nu0*cos4n2*cosh4m1 + tsbc*mu0*nu0*cos4n2*cosh4m1 - 
       tsbc*mu1*nu0*cos4n2*cosh4m1 + tsc2*pn02*cos4n2*cosh4m1 + 
       tsc1*nu1*cos4n2*cosh4m1 - tsbc*mu0*nu1*cos4n2*cosh4m1 + 
       tsbc*mu1*nu1*cos4n2*cosh4m1 - ftc2*nu0*nu1*cos4n2*cosh4m1 + 
       tsc2*pn22*cos4n2*cosh4m1 + 2304*a0*zc*cos3n1*coshm2 + 
       feb2*zc*cos3n1*coshm2 - 512*c2*zc*cos3n1*coshm2 - 
       tob1*mu0*zc*cos3n1*coshm2 + tob2*pm02*zc*cos3n1*coshm2 + 
       tob1*mu2*zc*cos3n1*coshm2 - feb2*mu0*mu2*zc*cos3n1*coshm2 + 
       tob2*pm22*zc*cos3n1*coshm2 - toc1*nu0*zc*cos3n1*coshm2 + 
       tobc*mu0*nu0*zc*cos3n1*coshm2 - tobc*mu2*nu0*zc*cos3n1*coshm2 + 
       toc2*pn02*zc*cos3n1*coshm2 + toc1*nu1*zc*cos3n1*coshm2 - 
       tobc*mu0*nu1*zc*cos3n1*coshm2 + tobc*mu2*nu1*zc*cos3n1*coshm2 - 
       fec2*nu0*nu1*zc*cos3n1*coshm2 + 
       toc2*pn12*zc*cos3n1*coshm2 - 2304*a0*zc*cos3n2*coshm2 - 
       feb2*zc*cos3n2*coshm2 + 512*c2*zc*cos3n2*coshm2 + 
       tob1*mu0*zc*cos3n2*coshm2 - tob2*pm02*zc*cos3n2*coshm2 - 
       tob1*mu2*zc*cos3n2*coshm2 + feb2*mu0*mu2*zc*cos3n2*coshm2 - 
       tob2*pm22*zc*cos3n2*coshm2 + toc1*nu0*zc*cos3n2*coshm2 - 
       tobc*mu0*nu0*zc*cos3n2*coshm2 + tobc*mu2*nu0*zc*cos3n2*coshm2 - 
       toc2*pn02*zc*cos3n2*coshm2 - toc1*nu1*zc*cos3n2*coshm2 + 
       tobc*mu0*nu1*zc*cos3n2*coshm2 - tobc*mu2*nu1*zc*cos3n2*coshm2 + 
       fec2*nu0*nu1*zc*cos3n2*coshm2 - 
       toc2*pn22*zc*cos3n2*coshm2 - 216*a0*cos4n1*cosh2m2 - 
       108*b2*cos4n1*cosh2m2 + 27*c2*cos4n1*cosh2m2 + 
       tsb1*mu0*cos4n1*cosh2m2 - tsb2*pm02*cos4n1*cosh2m2 - 
       tsb1*mu2*cos4n1*cosh2m2 + ftb2*mu0*mu2*cos4n1*cosh2m2 - 
       tsb2*pm22*cos4n1*cosh2m2 + tsc1*nu0*cos4n1*cosh2m2 - 
       tsbc*mu0*nu0*cos4n1*cosh2m2 + tsbc*mu2*nu0*cos4n1*cosh2m2 - 
       tsc2*pn02*cos4n1*cosh2m2 - tsc1*nu1*cos4n1*cosh2m2 + 
       tsbc*mu0*nu1*cos4n1*cosh2m2 - tsbc*mu2*nu1*cos4n1*cosh2m2 + 
       ftc2*nu0*nu1*cos4n1*cosh2m2 - tsc2*pn12*cos4n1*cosh2m2 + 
       216*a0*cos4n2*cosh2m2 + 108*b2*cos4n2*cosh2m2 - 
       27*c2*cos4n2*cosh2m2 - tsb1*mu0*cos4n2*cosh2m2 + 
       tsb2*pm02*cos4n2*cosh2m2 + tsb1*mu2*cos4n2*cosh2m2 - 
       ftb2*mu0*mu2*cos4n2*cosh2m2 + tsb2*pm22*cos4n2*cosh2m2 - 
       tsc1*nu0*cos4n2*cosh2m2 + tsbc*mu0*nu0*cos4n2*cosh2m2 - 
       tsbc*mu2*nu0*cos4n2*cosh2m2 + tsc2*pn02*cos4n2*cosh2m2 + 
       tsc1*nu1*cos4n2*cosh2m2 - tsbc*mu0*nu1*cos4n2*cosh2m2 + 
       tsbc*mu2*nu1*cos4n2*cosh2m2 - ftc2*nu0*nu1*cos4n2*cosh2m2 + 
       tsc2*pn22*cos4n2*cosh2m2 - 2304*a0*zc*cosn1*cosh3m2 - 
       512*b2*zc*cosn1*cosh3m2 + fec2*zc*cosn1*cosh3m2 + 
       tob1*mu0*zc*cosn1*cosh3m2 - tob2*pm02*zc*cosn1*cosh3m2 - 
       tob1*mu2*zc*cosn1*cosh3m2 + feb2*mu0*mu2*zc*cosn1*cosh3m2 - 
       tob2*pm22*zc*cosn1*cosh3m2 + toc1*nu0*zc*cosn1*cosh3m2 - 
       tobc*mu0*nu0*zc*cosn1*cosh3m2 + tobc*mu2*nu0*zc*cosn1*cosh3m2 - 
       toc2*pn02*zc*cosn1*cosh3m2 - toc1*nu1*zc*cosn1*cosh3m2 + 
       tobc*mu0*nu1*zc*cosn1*cosh3m2 - tobc*mu2*nu1*zc*cosn1*cosh3m2 + 
       fec2*nu0*nu1*zc*cosn1*cosh3m2 - 
       toc2*pn12*zc*cosn1*cosh3m2 + 2304*a0*zc*cosn2*cosh3m2 + 
       512*b2*zc*cosn2*cosh3m2 - fec2*zc*cosn2*cosh3m2 - 
       tob1*mu0*zc*cosn2*cosh3m2 + tob2*pm02*zc*cosn2*cosh3m2 + 
       tob1*mu2*zc*cosn2*cosh3m2 - feb2*mu0*mu2*zc*cosn2*cosh3m2 + 
       tob2*pm22*zc*cosn2*cosh3m2 - toc1*nu0*zc*cosn2*cosh3m2 + 
       tobc*mu0*nu0*zc*cosn2*cosh3m2 - tobc*mu2*nu0*zc*cosn2*cosh3m2 + 
       toc2*pn02*zc*cosn2*cosh3m2 + toc1*nu1*zc*cosn2*cosh3m2 - 
       tobc*mu0*nu1*zc*cosn2*cosh3m2 + tobc*mu2*nu1*zc*cosn2*cosh3m2 - 
       fec2*nu0*nu1*zc*cosn2*cosh3m2 + 
       toc2*pn22*zc*cosn2*cosh3m2 + 216*a0*cos4n1*cosh4m2 + 
       27*b2*cos4n1*cosh4m2 - 108*c2*cos4n1*cosh4m2 - 
       tsb1*mu0*cos4n1*cosh4m2 + tsb2*pm02*cos4n1*cosh4m2 + 
       tsb1*mu2*cos4n1*cosh4m2 - ftb2*mu0*mu2*cos4n1*cosh4m2 + 
       tsb2*pm22*cos4n1*cosh4m2 - tsc1*nu0*cos4n1*cosh4m2 + 
       tsbc*mu0*nu0*cos4n1*cosh4m2 - tsbc*mu2*nu0*cos4n1*cosh4m2 + 
       tsc2*pn02*cos4n1*cosh4m2 + tsc1*nu1*cos4n1*cosh4m2 - 
       tsbc*mu0*nu1*cos4n1*cosh4m2 + tsbc*mu2*nu1*cos4n1*cosh4m2 - 
       ftc2*nu0*nu1*cos4n1*cosh4m2 + tsc2*pn12*cos4n1*cosh4m2 - 
       216*a0*cos4n2*cosh4m2 - 27*b2*cos4n2*cosh4m2 + 
       108*c2*cos4n2*cosh4m2 + tsb1*mu0*cos4n2*cosh4m2 - 
       tsb2*pm02*cos4n2*cosh4m2 - tsb1*mu2*cos4n2*cosh4m2 + 
       ftb2*mu0*mu2*cos4n2*cosh4m2 - tsb2*pm22*cos4n2*cosh4m2 + 
       tsc1*nu0*cos4n2*cosh4m2 - tsbc*mu0*nu0*cos4n2*cosh4m2 + 
       tsbc*mu2*nu0*cos4n2*cosh4m2 - tsc2*pn02*cos4n2*cosh4m2 - 
       tsc1*nu1*cos4n2*cosh4m2 + tsbc*mu0*nu1*cos4n2*cosh4m2 - 
       tsbc*mu2*nu1*cos4n2*cosh4m2 + ftc2*nu0*nu1*cos4n2*cosh4m2 - 
       tsc2*pn22*cos4n2*cosh4m2 - toc1*zc*cosh3m1*sinn1 + 
       tobc*mu0*zc*cosh3m1*sinn1 - tobc*mu1*zc*cosh3m1*sinn1 + 
       fec2*nu0*zc*cosh3m1*sinn1 - fec2*nu1*zc*cosh3m1*sinn1 + 
       toc1*zc*cosh3m2*sinn1 - tobc*mu0*zc*cosh3m2*sinn1 + 
       tobc*mu2*zc*cosh3m2*sinn1 - fec2*nu0*zc*cosh3m2*sinn1 + 
       fec2*nu1*zc*cosh3m2*sinn1 + 108*c1*cosh4m1*sin4n1 - 
       108*bc*mu0*cosh4m1*sin4n1 + 108*bc*mu1*cosh4m1*sin4n1 - 
       tsc2*nu0*cosh4m1*sin4n1 + tsc2*nu1*cosh4m1*sin4n1 - 
       108*c1*cosh4m2*sin4n1 + 108*bc*mu0*cosh4m2*sin4n1 - 
       108*bc*mu2*cosh4m2*sin4n1 + tsc2*nu0*cosh4m2*sin4n1 - 
       tsc2*nu1*cosh4m2*sin4n1 - 768*c1*zc*coshm2*sin3n1 + 
       768*bc*mu0*zc*coshm2*sin3n1 - 768*bc*mu2*zc*coshm2*sin3n1 + 
       1536*c2*nu0*zc*coshm2*sin3n1 - 1536*c2*nu1*zc*coshm2*sin3n1 + 
       54*c1*cosh2m2*sin4n1 - 54*bc*mu0*cosh2m2*sin4n1 + 
       54*bc*mu2*cosh2m2*sin4n1 - 108*c2*nu0*cosh2m2*sin4n1 + 
       108*c2*nu1*cosh2m2*sin4n1 + toc1*zc*cosh3m1*sinn2 - 
       tobc*mu0*zc*cosh3m1*sinn2 + tobc*mu1*zc*cosh3m1*sinn2 - 
       fec2*nu0*zc*cosh3m1*sinn2 + fec2*nu1*zc*cosh3m1*sinn2 - 
       toc1*zc*cosh3m2*sinn2 + tobc*mu0*zc*cosh3m2*sinn2 - 
       tobc*mu2*zc*cosh3m2*sinn2 + fec2*nu0*zc*cosh3m2*sinn2 - 
       fec2*nu1*zc*cosh3m2*sinn2 - 108*c1*cosh4m1*sin4n2 + 
       108*bc*mu0*cosh4m1*sin4n2 - 108*bc*mu1*cosh4m1*sin4n2 + 
       tsc2*nu0*cosh4m1*sin4n2 - tsc2*nu1*cosh4m1*sin4n2 + 
       108*c1*cosh4m2*sin4n2 - 108*bc*mu0*cosh4m2*sin4n2 + 
       108*bc*mu2*cosh4m2*sin4n2 - tsc2*nu0*cosh4m2*sin4n2 + 
       tsc2*nu1*cosh4m2*sin4n2 + 768*c1*zc*coshm2*sin3n2 - 
       768*bc*mu0*zc*coshm2*sin3n2 + 768*bc*mu2*zc*coshm2*sin3n2 - 
       1536*c2*nu0*zc*coshm2*sin3n2 + 1536*c2*nu1*zc*coshm2*sin3n2 + 
       256*zc*coshm1*(-((9*a0 + 9*b2*(2 + pm012) + 
               c2*(-2 + 9*pn012) + 9*c1*(-nu0 + nu1) - 
               9*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos3n1) + 
          (9*a0 + 9*b2*(2 + pm012) + c2*(-2 + 9*pn022) + 
             9*c1*(-nu0 + nu2) - 9*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos3n2 + 
          3*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin3n1 - 
          3*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin3n2) - 
       54*c1*cosh2m2*sin4n2 + 54*bc*mu0*cosh2m2*sin4n2 - 
       54*bc*mu2*cosh2m2*sin4n2 + 108*c2*nu0*cosh2m2*sin4n2 - 
       108*c2*nu1*cosh2m2*sin4n2 + 
       27*cosh2m1*((8*a0 + 4*b2*(1 + 2*pm012) + 
             c2*(-1 + 8*pn012) + 8*c1*(-nu0 + nu1) - 
             8*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos4n1 + 
          (-8*a0 - 4*b2*(1 + 2*pm012) + c2*(1 - 8*pn022) + 
             8*c1*(nu0 - nu2) + 8*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos4n2 - 
          2*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n1 + 
          2*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin4n2) + 
       tob1*zc*cos3n1*sinhm1 - feb2*mu0*zc*cos3n1*sinhm1 + 
       feb2*mu1*zc*cos3n1*sinhm1 - tobc*nu0*zc*cos3n1*sinhm1 + 
       tobc*nu1*zc*cos3n1*sinhm1 - tob1*zc*cos3n2*sinhm1 + 
       feb2*mu0*zc*cos3n2*sinhm1 - feb2*mu1*zc*cos3n2*sinhm1 + 
       tobc*nu0*zc*cos3n2*sinhm1 - tobc*nu1*zc*cos3n2*sinhm1 - 
       768*bc*zc*sin3n1*sinhm1 + 768*bc*zc*sin3n2*sinhm1 - 
       108*b1*cos4n1*sinh2m1 + tsb2*mu0*cos4n1*sinh2m1 - 
       tsb2*mu1*cos4n1*sinh2m1 + 108*bc*nu0*cos4n1*sinh2m1 - 
       108*bc*nu1*cos4n1*sinh2m1 + 108*b1*cos4n2*sinh2m1 - 
       tsb2*mu0*cos4n2*sinh2m1 + tsb2*mu1*cos4n2*sinh2m1 - 
       108*bc*nu0*cos4n2*sinh2m1 + 108*bc*nu1*cos4n2*sinh2m1 + 
       27*bc*sin4n1*sinh2m1 - 27*bc*sin4n2*sinh2m1 - 
       768*b1*zc*cosn1*sinh3m1 + 1536*b2*mu0*zc*cosn1*sinh3m1 - 
       1536*b2*mu1*zc*cosn1*sinh3m1 + 768*bc*nu0*zc*cosn1*sinh3m1 - 
       768*bc*nu1*zc*cosn1*sinh3m1 + 768*b1*zc*cosn2*sinh3m1 - 
       1536*b2*mu0*zc*cosn2*sinh3m1 + 1536*b2*mu1*zc*cosn2*sinh3m1 - 
       768*bc*nu0*zc*cosn2*sinh3m1 + 768*bc*nu1*zc*cosn2*sinh3m1 + 
       768*bc*zc*sinn1*sinh3m1 - 768*bc*zc*sinn2*sinh3m1 + 
       54*b1*cos4n1*sinh4m1 - 108*b2*mu0*cos4n1*sinh4m1 + 
       108*b2*mu1*cos4n1*sinh4m1 - 54*bc*nu0*cos4n1*sinh4m1 + 
       54*bc*nu1*cos4n1*sinh4m1 - 54*b1*cos4n2*sinh4m1 + 
       108*b2*mu0*cos4n2*sinh4m1 - 108*b2*mu1*cos4n2*sinh4m1 + 
       54*bc*nu0*cos4n2*sinh4m1 - 54*bc*nu1*cos4n2*sinh4m1 - 
       27*bc*sin4n1*sinh4m1 + 27*bc*sin4n2*sinh4m1 - 
       tob1*zc*cos3n1*sinhm2 + feb2*mu0*zc*cos3n1*sinhm2 - 
       feb2*mu2*zc*cos3n1*sinhm2 + tobc*nu0*zc*cos3n1*sinhm2 - 
       tobc*nu1*zc*cos3n1*sinhm2 + tob1*zc*cos3n2*sinhm2 - 
       feb2*mu0*zc*cos3n2*sinhm2 + feb2*mu2*zc*cos3n2*sinhm2 - 
       tobc*nu0*zc*cos3n2*sinhm2 + tobc*nu1*zc*cos3n2*sinhm2 + 
       768*bc*zc*sin3n1*sinhm2 - 768*bc*zc*sin3n2*sinhm2 + 
       108*b1*cos4n1*sinh2m2 - tsb2*mu0*cos4n1*sinh2m2 + 
       tsb2*mu2*cos4n1*sinh2m2 - 108*bc*nu0*cos4n1*sinh2m2 + 
       108*bc*nu1*cos4n1*sinh2m2 - 108*b1*cos4n2*sinh2m2 + 
       tsb2*mu0*cos4n2*sinh2m2 - tsb2*mu2*cos4n2*sinh2m2 + 
       108*bc*nu0*cos4n2*sinh2m2 - 108*bc*nu1*cos4n2*sinh2m2 - 
       27*bc*sin4n1*sinh2m2 + 27*bc*sin4n2*sinh2m2 + 
       768*b1*zc*cosn1*sinh3m2 - 1536*b2*mu0*zc*cosn1*sinh3m2 + 
       1536*b2*mu2*zc*cosn1*sinh3m2 - 768*bc*nu0*zc*cosn1*sinh3m2 + 
       768*bc*nu1*zc*cosn1*sinh3m2 - 768*b1*zc*cosn2*sinh3m2 + 
       1536*b2*mu0*zc*cosn2*sinh3m2 - 1536*b2*mu2*zc*cosn2*sinh3m2 + 
       768*bc*nu0*zc*cosn2*sinh3m2 - 768*bc*nu1*zc*cosn2*sinh3m2 - 
       768*bc*zc*sinn1*sinh3m2 + 768*bc*zc*sinn2*sinh3m2 + 
       27*(-2*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cos4n1 + 
          2*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cos4n2 + 
          bc*(sin4n1 - sin4n2))*sinh4m2);

  return val;
}


//CPMZ marked: nu2 bug in convert function
//x2 or y2 terms
double second_order_fx2dV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc)
{
  double f = 1./216000.;
  double pa5 = pow(a,5.);
  double fp = f*pa5;

  double pm02 = pow(mu0,2.);
  double pm12 = pow(mu1,2.);
  double pm22 = pow(mu2,2.);
  double pn02 = pow(nu0,2.);
  double pn12 = pow(nu1,2.);
  double pn22 = pow(nu2,2.);

  double pm012 = pow(mu0-mu1,2.);
  double pm022 = pow(mu0-mu2,2.);
  double pn012 = pow(nu0-nu1,2.);
  double pn022 = pow(nu0-nu2,2.);

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosh5m1 = cosh(5.*mu1);
  double cosh5m2 = cosh(5.*mu2);
  double sinh5m1 = sinh(5.*mu1);
  double sinh5m2 = sinh(5.*mu2);

  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);
  double sin3n1 = sin(3.*nu1);
  double sin3n2 = sin(3.*nu2);
  double cos5n1 = cos(5.*nu1);
  double cos5n2 = cos(5.*nu2);
  double sin5n1 = sin(5.*nu1);
  double sin5n2 = sin(5.*nu2);


  double tfbc = 2025.*bc;
  double tfb1 = 2025.*b1;
  double tfb2 = 2025.*b2;
  double tfc1 = 2025.*c1;
  double tfc2 = 2025.*c2;

  double ffbc = 4050.*bc;
  double ffb1 = 4050.*b1;
  double ffb2 = 4050.*b2;
  double ffc1 = 4050.*c1;
  double ffc2 = 4050.*c2;

  double ttbc = 225.*bc;
  double ttb1 = 225.*b1;
  double ttb2 = 225.*b2;
  double ttc1 = 225.*c1;
  double ttc2 = 225.*c2;

  double fsbc = 5625.*bc;
  double fsb1 = 5625.*b1;
  double fsb2 = 5625.*b2;
  double fsc1 = 5625.*c1;
  double fsc2 = 5625.*c2;

  double val = fp*(-2025*a0*cosn1*cosh5m1 - 162*b2*cosn1*cosh5m1 + 
       ffc2*cosn1*cosh5m1 + tfb1*mu0*cosn1*cosh5m1 - 
       tfb2*pm02*cosn1*cosh5m1 - tfb1*mu1*cosn1*cosh5m1 + 
       ffb2*mu0*mu1*cosn1*cosh5m1 - tfb2*pm12*cosn1*cosh5m1 + 
       tfc1*nu0*cosn1*cosh5m1 - tfbc*mu0*nu0*cosn1*cosh5m1 + 
       tfbc*mu1*nu0*cosn1*cosh5m1 - tfc2*pn02*cosn1*cosh5m1 - 
       tfc1*nu1*cosn1*cosh5m1 + tfbc*mu0*nu1*cosn1*cosh5m1 - 
       tfbc*mu1*nu1*cosn1*cosh5m1 + ffc2*nu0*nu1*cosn1*cosh5m1 - 
       tfc2*pn12*cosn1*cosh5m1 + 225*a0*cos3n1*cosh5m1 + 
       18*b2*cos3n1*cosh5m1 - 50*c2*cos3n1*cosh5m1 - 
       ttb1*mu0*cos3n1*cosh5m1 + ttb2*pm02*cos3n1*cosh5m1 + 
       ttb1*mu1*cos3n1*cosh5m1 - 450*b2*mu0*mu1*cos3n1*cosh5m1 + 
       ttb2*pm12*cos3n1*cosh5m1 - ttc1*nu0*cos3n1*cosh5m1 + 
       ttbc*mu0*nu0*cos3n1*cosh5m1 - ttbc*mu1*nu0*cos3n1*cosh5m1 + 
       ttc2*pn02*cos3n1*cosh5m1 + ttc1*nu1*cos3n1*cosh5m1 - 
       ttbc*mu0*nu1*cos3n1*cosh5m1 + ttbc*mu1*nu1*cos3n1*cosh5m1 - 
       450*c2*nu0*nu1*cos3n1*cosh5m1 + ttc2*pn12*cos3n1*cosh5m1 + 
       2025*a0*cosn2*cosh5m1 + 162*b2*cosn2*cosh5m1 - 
       ffc2*cosn2*cosh5m1 - tfb1*mu0*cosn2*cosh5m1 + 
       tfb2*pm02*cosn2*cosh5m1 + tfb1*mu1*cosn2*cosh5m1 - 
       ffb2*mu0*mu1*cosn2*cosh5m1 + tfb2*pm12*cosn2*cosh5m1 - 
       tfc1*nu0*cosn2*cosh5m1 + tfbc*mu0*nu0*cosn2*cosh5m1 - 
       tfbc*mu1*nu0*cosn2*cosh5m1 + tfc2*pn02*cosn2*cosh5m1 + 
       tfc1*nu1*cosn2*cosh5m1 - tfbc*mu0*nu1*cosn2*cosh5m1 + 
       tfbc*mu1*nu1*cosn2*cosh5m1 - ffc2*nu0*nu1*cosn2*cosh5m1 + 
       tfc2*pn22*cosn2*cosh5m1 - 225*a0*cos3n2*cosh5m1 - 
       18*b2*cos3n2*cosh5m1 + 50*c2*cos3n2*cosh5m1 + 
       ttb1*mu0*cos3n2*cosh5m1 - ttb2*pm02*cos3n2*cosh5m1 - 
       ttb1*mu1*cos3n2*cosh5m1 + 450*b2*mu0*mu1*cos3n2*cosh5m1 - 
       ttb2*pm12*cos3n2*cosh5m1 + ttc1*nu0*cos3n2*cosh5m1 - 
       ttbc*mu0*nu0*cos3n2*cosh5m1 + ttbc*mu1*nu0*cos3n2*cosh5m1 - 
       ttc2*pn02*cos3n2*cosh5m1 - ttc1*nu1*cos3n2*cosh5m1 + 
       ttbc*mu0*nu1*cos3n2*cosh5m1 - ttbc*mu1*nu1*cos3n2*cosh5m1 + 
       450*c2*nu0*nu1*cos3n2*cosh5m1 - ttc2*pn22*cos3n2*cosh5m1 + 
       5625*a0*cos3n1*coshm2 + 11250*b2*cos3n1*coshm2 - 
       1250*c2*cos3n1*coshm2 - fsb1*mu0*cos3n1*coshm2 + 
       fsb2*pm02*cos3n1*coshm2 + fsb1*mu2*cos3n1*coshm2 - 
       11250*b2*mu0*mu2*cos3n1*coshm2 + fsb2*pm22*cos3n1*coshm2 - 
       fsc1*nu0*cos3n1*coshm2 + fsbc*mu0*nu0*cos3n1*coshm2 - 
       fsbc*mu2*nu0*cos3n1*coshm2 + fsc2*pn02*cos3n1*coshm2 + 
       fsc1*nu1*cos3n1*coshm2 - fsbc*mu0*nu1*cos3n1*coshm2 + 
       fsbc*mu2*nu1*cos3n1*coshm2 - 11250*c2*nu0*nu1*cos3n1*coshm2 + 
       fsc2*pn12*cos3n1*coshm2 - 2025*a0*cos5n1*coshm2 - 
       ffb2*cos5n1*coshm2 + 162*c2*cos5n1*coshm2 + 
       tfb1*mu0*cos5n1*coshm2 - tfb2*pm02*cos5n1*coshm2 - 
       tfb1*mu2*cos5n1*coshm2 + ffb2*mu0*mu2*cos5n1*coshm2 - 
       tfb2*pm22*cos5n1*coshm2 + tfc1*nu0*cos5n1*coshm2 - 
       tfbc*mu0*nu0*cos5n1*coshm2 + tfbc*mu2*nu0*cos5n1*coshm2 - 
       tfc2*pn02*cos5n1*coshm2 - tfc1*nu1*cos5n1*coshm2 + 
       tfbc*mu0*nu1*cos5n1*coshm2 - tfbc*mu2*nu1*cos5n1*coshm2 + 
       ffc2*nu0*nu1*cos5n1*coshm2 - tfc2*pn12*cos5n1*coshm2 - 
       5625*a0*cos3n2*coshm2 - 11250*b2*cos3n2*coshm2 + 
       1250*c2*cos3n2*coshm2 + fsb1*mu0*cos3n2*coshm2 - 
       fsb2*pm02*cos3n2*coshm2 - fsb1*mu2*cos3n2*coshm2 + 
       11250*b2*mu0*mu2*cos3n2*coshm2 - fsb2*pm22*cos3n2*coshm2 + 
       fsc1*nu0*cos3n2*coshm2 - fsbc*mu0*nu0*cos3n2*coshm2 + 
       fsbc*mu2*nu0*cos3n2*coshm2 - fsc2*pn02*cos3n2*coshm2 - 
       fsc1*nu1*cos3n2*coshm2 + fsbc*mu0*nu1*cos3n2*coshm2 - 
       fsbc*mu2*nu1*cos3n2*coshm2 + 11250*c2*nu0*nu1*cos3n2*coshm2 - 
       fsc2*pn22*cos3n2*coshm2 + 2025*a0*cos5n2*coshm2 + 
       ffb2*cos5n2*coshm2 - 162*c2*cos5n2*coshm2 - 
       tfb1*mu0*cos5n2*coshm2 + tfb2*pm02*cos5n2*coshm2 + 
       tfb1*mu2*cos5n2*coshm2 - ffb2*mu0*mu2*cos5n2*coshm2 + 
       tfb2*pm22*cos5n2*coshm2 - tfc1*nu0*cos5n2*coshm2 + 
       tfbc*mu0*nu0*cos5n2*coshm2 - tfbc*mu2*nu0*cos5n2*coshm2 + 
       tfc2*pn02*cos5n2*coshm2 + tfc1*nu1*cos5n2*coshm2 - 
       tfbc*mu0*nu1*cos5n2*coshm2 + tfbc*mu2*nu1*cos5n2*coshm2 - 
       ffc2*nu0*nu1*cos5n2*coshm2 + tfc2*pn22*cos5n2*coshm2 - 
       5625*a0*cosn1*cosh3m2 - 1250*b2*cosn1*cosh3m2 + 
       11250*c2*cosn1*cosh3m2 + fsb1*mu0*cosn1*cosh3m2 - 
       fsb2*pm02*cosn1*cosh3m2 - fsb1*mu2*cosn1*cosh3m2 + 
       11250*b2*mu0*mu2*cosn1*cosh3m2 - fsb2*pm22*cosn1*cosh3m2 + 
       fsc1*nu0*cosn1*cosh3m2 - fsbc*mu0*nu0*cosn1*cosh3m2 + 
       fsbc*mu2*nu0*cosn1*cosh3m2 - fsc2*pn02*cosn1*cosh3m2 - 
       fsc1*nu1*cosn1*cosh3m2 + fsbc*mu0*nu1*cosn1*cosh3m2 - 
       fsbc*mu2*nu1*cosn1*cosh3m2 + 11250*c2*nu0*nu1*cosn1*cosh3m2 - 
       fsc2*pn12*cosn1*cosh3m2 + 225*a0*cos5n1*cosh3m2 + 
       50*b2*cos5n1*cosh3m2 - 18*c2*cos5n1*cosh3m2 - 
       ttb1*mu0*cos5n1*cosh3m2 + ttb2*pm02*cos5n1*cosh3m2 + 
       ttb1*mu2*cos5n1*cosh3m2 - 450*b2*mu0*mu2*cos5n1*cosh3m2 + 
       ttb2*pm22*cos5n1*cosh3m2 - ttc1*nu0*cos5n1*cosh3m2 + 
       ttbc*mu0*nu0*cos5n1*cosh3m2 - ttbc*mu2*nu0*cos5n1*cosh3m2 + 
       ttc2*pn02*cos5n1*cosh3m2 + ttc1*nu1*cos5n1*cosh3m2 - 
       ttbc*mu0*nu1*cos5n1*cosh3m2 + ttbc*mu2*nu1*cos5n1*cosh3m2 - 
       450*c2*nu0*nu1*cos5n1*cosh3m2 + ttc2*pn12*cos5n1*cosh3m2 + 
       5625*a0*cosn2*cosh3m2 + 1250*b2*cosn2*cosh3m2 - 
       11250*c2*cosn2*cosh3m2 - fsb1*mu0*cosn2*cosh3m2 + 
       fsb2*pm02*cosn2*cosh3m2 + fsb1*mu2*cosn2*cosh3m2 - 
       11250*b2*mu0*mu2*cosn2*cosh3m2 + fsb2*pm22*cosn2*cosh3m2 - 
       fsc1*nu0*cosn2*cosh3m2 + fsbc*mu0*nu0*cosn2*cosh3m2 - 
       fsbc*mu2*nu0*cosn2*cosh3m2 + fsc2*pn02*cosn2*cosh3m2 + 
       fsc1*nu1*cosn2*cosh3m2 - fsbc*mu0*nu1*cosn2*cosh3m2 + 
       fsbc*mu2*nu1*cosn2*cosh3m2 - 11250*c2*nu0*nu1*cosn2*cosh3m2 + 
       fsc2*pn22*cosn2*cosh3m2 - 225*a0*cos5n2*cosh3m2 - 
       50*b2*cos5n2*cosh3m2 + 18*c2*cos5n2*cosh3m2 + 
       ttb1*mu0*cos5n2*cosh3m2 - ttb2*pm02*cos5n2*cosh3m2 - 
       ttb1*mu2*cos5n2*cosh3m2 + 450*b2*mu0*mu2*cos5n2*cosh3m2 - 
       ttb2*pm22*cos5n2*cosh3m2 + ttc1*nu0*cos5n2*cosh3m2 - 
       ttbc*mu0*nu0*cos5n2*cosh3m2 + ttbc*mu2*nu0*cos5n2*cosh3m2 - 
       ttc2*pn02*cos5n2*cosh3m2 - ttc1*nu1*cos5n2*cosh3m2 + 
       ttbc*mu0*nu1*cos5n2*cosh3m2 - ttbc*mu2*nu1*cos5n2*cosh3m2 + 
       450*c2*nu0*nu1*cos5n2*cosh3m2 - ttc2*pn22*cos5n2*cosh3m2 + 
       2025*a0*cosn1*cosh5m2 + 162*b2*cosn1*cosh5m2 - 
       ffc2*cosn1*cosh5m2 - tfb1*mu0*cosn1*cosh5m2 + 
       tfb2*pm02*cosn1*cosh5m2 + tfb1*mu2*cosn1*cosh5m2 - 
       ffb2*mu0*mu2*cosn1*cosh5m2 + tfb2*pm22*cosn1*cosh5m2 - 
       tfc1*nu0*cosn1*cosh5m2 + tfbc*mu0*nu0*cosn1*cosh5m2 - 
       tfbc*mu2*nu0*cosn1*cosh5m2 + tfc2*pn02*cosn1*cosh5m2 + 
       tfc1*nu1*cosn1*cosh5m2 - tfbc*mu0*nu1*cosn1*cosh5m2 + 
       tfbc*mu2*nu1*cosn1*cosh5m2 - ffc2*nu0*nu1*cosn1*cosh5m2 + 
       tfc2*pn12*cosn1*cosh5m2 - 225*a0*cos3n1*cosh5m2 - 
       18*b2*cos3n1*cosh5m2 + 50*c2*cos3n1*cosh5m2 + 
       ttb1*mu0*cos3n1*cosh5m2 - ttb2*pm02*cos3n1*cosh5m2 - 
       ttb1*mu2*cos3n1*cosh5m2 + 450*b2*mu0*mu2*cos3n1*cosh5m2 - 
       ttb2*pm22*cos3n1*cosh5m2 + ttc1*nu0*cos3n1*cosh5m2 - 
       ttbc*mu0*nu0*cos3n1*cosh5m2 + ttbc*mu2*nu0*cos3n1*cosh5m2 - 
       ttc2*pn02*cos3n1*cosh5m2 - ttc1*nu1*cos3n1*cosh5m2 + 
       ttbc*mu0*nu1*cos3n1*cosh5m2 - ttbc*mu2*nu1*cos3n1*cosh5m2 + 
       450*c2*nu0*nu1*cos3n1*cosh5m2 - ttc2*pn12*cos3n1*cosh5m2 - 
       2025*a0*cosn2*cosh5m2 - 162*b2*cosn2*cosh5m2 + 
       ffc2*cosn2*cosh5m2 + tfb1*mu0*cosn2*cosh5m2 - 
       tfb2*pm02*cosn2*cosh5m2 - tfb1*mu2*cosn2*cosh5m2 + 
       ffb2*mu0*mu2*cosn2*cosh5m2 - tfb2*pm22*cosn2*cosh5m2 + 
       tfc1*nu0*cosn2*cosh5m2 - tfbc*mu0*nu0*cosn2*cosh5m2 + 
       tfbc*mu2*nu0*cosn2*cosh5m2 - tfc2*pn02*cosn2*cosh5m2 - 
       tfc1*nu1*cosn2*cosh5m2 + tfbc*mu0*nu1*cosn2*cosh5m2 - 
       tfbc*mu2*nu1*cosn2*cosh5m2 + ffc2*nu0*nu1*cosn2*cosh5m2 - 
       tfc2*pn22*cosn2*cosh5m2 + 225*a0*cos3n2*cosh5m2 + 
       18*b2*cos3n2*cosh5m2 - 50*c2*cos3n2*cosh5m2 - 
       ttb1*mu0*cos3n2*cosh5m2 + ttb2*pm02*cos3n2*cosh5m2 + 
       ttb1*mu2*cos3n2*cosh5m2 - 450*b2*mu0*mu2*cos3n2*cosh5m2 + 
       ttb2*pm22*cos3n2*cosh5m2 - ttc1*nu0*cos3n2*cosh5m2 + 
       ttbc*mu0*nu0*cos3n2*cosh5m2 - ttbc*mu2*nu0*cos3n2*cosh5m2 + 
       ttc2*pn02*cos3n2*cosh5m2 + ttc1*nu1*cos3n2*cosh5m2 - 
       ttbc*mu0*nu1*cos3n2*cosh5m2 + ttbc*mu2*nu1*cos3n2*cosh5m2 - 
       450*c2*nu0*nu1*cos3n2*cosh5m2 + ttc2*pn22*cos3n2*cosh5m2 + 
       tfc1*cosh5m1*sinn1 - tfbc*mu0*cosh5m1*sinn1 + 
       tfbc*mu1*cosh5m1*sinn1 - ffc2*nu0*cosh5m1*sinn1 + 
       ffc2*nu1*cosh5m1*sinn1 + fsc1*cosh3m2*sinn1 - 
       fsbc*mu0*cosh3m2*sinn1 + fsbc*mu2*cosh3m2*sinn1 - 
       11250*c2*nu0*cosh3m2*sinn1 + 11250*c2*nu1*cosh3m2*sinn1 - 
       tfc1*cosh5m2*sinn1 + tfbc*mu0*cosh5m2*sinn1 - 
       tfbc*mu2*cosh5m2*sinn1 + ffc2*nu0*cosh5m2*sinn1 - 
       ffc2*nu1*cosh5m2*sinn1 - 75*c1*cosh5m1*sin3n1 + 
       75*bc*mu0*cosh5m1*sin3n1 - 75*bc*mu1*cosh5m1*sin3n1 + 
       150*c2*nu0*cosh5m1*sin3n1 - 150*c2*nu1*cosh5m1*sin3n1 - 
       1875*c1*coshm2*sin3n1 + 1875*bc*mu0*coshm2*sin3n1 - 
       1875*bc*mu2*coshm2*sin3n1 + 3750*c2*nu0*coshm2*sin3n1 - 
       3750*c2*nu1*coshm2*sin3n1 + 75*c1*cosh5m2*sin3n1 - 
       75*bc*mu0*cosh5m2*sin3n1 + 75*bc*mu2*cosh5m2*sin3n1 - 
       150*c2*nu0*cosh5m2*sin3n1 + 150*c2*nu1*cosh5m2*sin3n1 + 
       405*c1*coshm2*sin5n1 - 405*bc*mu0*coshm2*sin5n1 + 
       405*bc*mu2*coshm2*sin5n1 - 810*c2*nu0*coshm2*sin5n1 + 
       810*c2*nu1*coshm2*sin5n1 - 45*c1*cosh3m2*sin5n1 + 
       45*bc*mu0*cosh3m2*sin5n1 - 45*bc*mu2*cosh3m2*sin5n1 + 
       90*c2*nu0*cosh3m2*sin5n1 - 90*c2*nu1*cosh3m2*sin5n1 - 
       tfc1*cosh5m1*sinn2 + tfbc*mu0*cosh5m1*sinn2 - 
       tfbc*mu1*cosh5m1*sinn2 + ffc2*nu0*cosh5m1*sinn2 - 
       ffc2*nu1*cosh5m1*sinn2 - fsc1*cosh3m2*sinn2 + 
       fsbc*mu0*cosh3m2*sinn2 - fsbc*mu2*cosh3m2*sinn2 + 
       11250*c2*nu0*cosh3m2*sinn2 - 11250*c2*nu1*cosh3m2*sinn2 + 
       tfc1*cosh5m2*sinn2 - tfbc*mu0*cosh5m2*sinn2 + 
       tfbc*mu2*cosh5m2*sinn2 - ffc2*nu0*cosh5m2*sinn2 + 
       ffc2*nu1*cosh5m2*sinn2 + 75*c1*cosh5m1*sin3n2 - 
       75*bc*mu0*cosh5m1*sin3n2 + 75*bc*mu1*cosh5m1*sin3n2 - 
       150*c2*nu0*cosh5m1*sin3n2 + 150*c2*nu1*cosh5m1*sin3n2 + 
       1875*c1*coshm2*sin3n2 - 1875*bc*mu0*coshm2*sin3n2 + 
       1875*bc*mu2*coshm2*sin3n2 - 3750*c2*nu0*coshm2*sin3n2 + 
       3750*c2*nu1*coshm2*sin3n2 - 75*c1*cosh5m2*sin3n2 + 
       75*bc*mu0*cosh5m2*sin3n2 - 75*bc*mu2*cosh5m2*sin3n2 + 
       150*c2*nu0*cosh5m2*sin3n2 - 150*c2*nu1*cosh5m2*sin3n2 - 
       405*c1*coshm2*sin5n2 + 405*bc*mu0*coshm2*sin5n2 - 
       405*bc*mu2*coshm2*sin5n2 + 810*c2*nu0*coshm2*sin5n2 - 
       810*c2*nu1*coshm2*sin5n2 + 45*c1*cosh3m2*sin5n2 - 
       45*bc*mu0*cosh3m2*sin5n2 + 45*bc*mu2*cosh3m2*sin5n2 - 
       90*c2*nu0*cosh3m2*sin5n2 + 90*c2*nu1*cosh3m2*sin5n2 + 
       cosh3m1*(625*(9*a0 + b2*(2 + 9*pm012) + 
             9*c2*(-2 + pn012) + 9*c1*(-nu0 + nu1) - 
             9*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 - 225*a0*cos5n1 - 
          50*b2*cos5n1 + 18*c2*cos5n1 + ttb1*mu0*cos5n1 - 
          ttb2*pm02*cos5n1 - ttb1*mu1*cos5n1 + 
          450*b2*mu0*mu1*cos5n1 - ttb2*pm12*cos5n1 + 
          ttc1*nu0*cos5n1 - ttbc*mu0*nu0*cos5n1 + ttbc*mu1*nu0*cos5n1 - 
          ttc2*pn02*cos5n1 - ttc1*nu1*cos5n1 + 
          ttbc*mu0*nu1*cos5n1 - ttbc*mu1*nu1*cos5n1 + 450*c2*nu0*nu1*cos5n1 - 
          ttc2*pn12*cos5n1 - 
          625*(9*a0 + b2*(2 + 9*pm012) + 9*c2*(-2 + pn022) + 
             9*c1*(-nu0 + nu2) - 9*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 + 
          (225*a0 + 25*b2*(2 + 9*pm012) - 
             225*(b1*(mu0 - mu1) + (c1 - bc*mu0 + bc*mu1)*(nu0 - nu2)) + 
             9*c2*(-2 + 25*pn022))*cos5n2 - fsc1*sinn1 + 
          fsbc*mu0*sinn1 - fsbc*mu1*sinn1 + 11250*c2*nu0*sinn1 - 
          11250*c2*nu1*sinn1 + 45*c1*sin5n1 - 45*bc*mu0*sin5n1 + 
          45*bc*mu1*sin5n1 - 90*c2*nu0*sin5n1 + 90*c2*nu1*sin5n1 + 
          5625*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sinn2 - 
          45*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin5n2) + 
       coshm1*(-625*(9*a0 + 9*b2*(2 + pm012) + 
             c2*(-2 + 9*pn012) + 9*c1*(-nu0 + nu1) - 
             9*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos3n1 + 
          81*(25*a0 + 25*b2*(2 + pm012) - 
             25*(b1*(mu0 - mu1) + (c1 - bc*mu0 + bc*mu1)*(nu0 - nu1)) + 
             c2*(-2 + 25*pn012))*cos5n1 + 
          625*(9*a0 + 9*b2*(2 + pm012) + c2*(-2 + 9*pn022) + 
             9*c1*(-nu0 + nu2) - 9*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos3n2 - 
          81*(25*a0 + 25*b2*(2 + pm012) - 
             25*(b1*(mu0 - mu1) + (c1 - bc*mu0 + bc*mu1)*(nu0 - nu2)) + 
             c2*(-2 + 25*pn022))*cos5n2 + 
          1875*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin3n1 - 
          405*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin5n1 - 
          1875*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin3n2 + 
          405*(c1 - bc*mu0 + bc*mu1 - 2*c2*nu0 + 2*c2*nu1)*sin5n2) + 
       fsb1*cos3n1*sinhm1 - 11250*b2*mu0*cos3n1*sinhm1 + 
       11250*b2*mu1*cos3n1*sinhm1 - fsbc*nu0*cos3n1*sinhm1 + 
       fsbc*nu1*cos3n1*sinhm1 - tfb1*cos5n1*sinhm1 + 
       ffb2*mu0*cos5n1*sinhm1 - ffb2*mu1*cos5n1*sinhm1 + 
       tfbc*nu0*cos5n1*sinhm1 - tfbc*nu1*cos5n1*sinhm1 - 
       fsb1*cos3n2*sinhm1 + 11250*b2*mu0*cos3n2*sinhm1 - 
       11250*b2*mu1*cos3n2*sinhm1 + fsbc*nu0*cos3n2*sinhm1 - 
       fsbc*nu1*cos3n2*sinhm1 + tfb1*cos5n2*sinhm1 - 
       ffb2*mu0*cos5n2*sinhm1 + ffb2*mu1*cos5n2*sinhm1 - 
       tfbc*nu0*cos5n2*sinhm1 + tfbc*nu1*cos5n2*sinhm1 - 
       1875*bc*sin3n1*sinhm1 + 405*bc*sin5n1*sinhm1 + 
       1875*bc*sin3n2*sinhm1 - 405*bc*sin5n2*sinhm1 - 
       1875*b1*cosn1*sinh3m1 + 3750*b2*mu0*cosn1*sinh3m1 - 
       3750*b2*mu1*cosn1*sinh3m1 + 1875*bc*nu0*cosn1*sinh3m1 - 
       1875*bc*nu1*cosn1*sinh3m1 + 75*b1*cos5n1*sinh3m1 - 
       150*b2*mu0*cos5n1*sinh3m1 + 150*b2*mu1*cos5n1*sinh3m1 - 
       75*bc*nu0*cos5n1*sinh3m1 + 75*bc*nu1*cos5n1*sinh3m1 + 
       1875*b1*cosn2*sinh3m1 - 3750*b2*mu0*cosn2*sinh3m1 + 
       3750*b2*mu1*cosn2*sinh3m1 - 1875*bc*nu0*cosn2*sinh3m1 + 
       1875*bc*nu1*cosn2*sinh3m1 - 75*b1*cos5n2*sinh3m1 + 
       150*b2*mu0*cos5n2*sinh3m1 - 150*b2*mu1*cos5n2*sinh3m1 + 
       75*bc*nu0*cos5n2*sinh3m1 - 75*bc*nu1*cos5n2*sinh3m1 + 
       1875*bc*sinn1*sinh3m1 - 15*bc*sin5n1*sinh3m1 - 
       1875*bc*sinn2*sinh3m1 + 15*bc*sin5n2*sinh3m1 + 
       405*b1*cosn1*sinh5m1 - 810*b2*mu0*cosn1*sinh5m1 + 
       810*b2*mu1*cosn1*sinh5m1 - 405*bc*nu0*cosn1*sinh5m1 + 
       405*bc*nu1*cosn1*sinh5m1 - 45*b1*cos3n1*sinh5m1 + 
       90*b2*mu0*cos3n1*sinh5m1 - 90*b2*mu1*cos3n1*sinh5m1 + 
       45*bc*nu0*cos3n1*sinh5m1 - 45*bc*nu1*cos3n1*sinh5m1 - 
       405*b1*cosn2*sinh5m1 + 810*b2*mu0*cosn2*sinh5m1 - 
       810*b2*mu1*cosn2*sinh5m1 + 405*bc*nu0*cosn2*sinh5m1 - 
       405*bc*nu1*cosn2*sinh5m1 + 45*b1*cos3n2*sinh5m1 - 
       90*b2*mu0*cos3n2*sinh5m1 + 90*b2*mu1*cos3n2*sinh5m1 - 
       45*bc*nu0*cos3n2*sinh5m1 + 45*bc*nu1*cos3n2*sinh5m1 - 
       405*bc*sinn1*sinh5m1 + 15*bc*sin3n1*sinh5m1 + 
       405*bc*sinn2*sinh5m1 - 15*bc*sin3n2*sinh5m1 - 
       fsb1*cos3n1*sinhm2 + 11250*b2*mu0*cos3n1*sinhm2 - 
       11250*b2*mu2*cos3n1*sinhm2 + fsbc*nu0*cos3n1*sinhm2 - 
       fsbc*nu1*cos3n1*sinhm2 + tfb1*cos5n1*sinhm2 - 
       ffb2*mu0*cos5n1*sinhm2 + ffb2*mu2*cos5n1*sinhm2 - 
       tfbc*nu0*cos5n1*sinhm2 + tfbc*nu1*cos5n1*sinhm2 + 
       fsb1*cos3n2*sinhm2 - 11250*b2*mu0*cos3n2*sinhm2 + 
       11250*b2*mu2*cos3n2*sinhm2 - fsbc*nu0*cos3n2*sinhm2 + 
       fsbc*nu1*cos3n2*sinhm2 - tfb1*cos5n2*sinhm2 + 
       ffb2*mu0*cos5n2*sinhm2 - ffb2*mu2*cos5n2*sinhm2 + 
       tfbc*nu0*cos5n2*sinhm2 - tfbc*nu1*cos5n2*sinhm2 + 
       1875*bc*sin3n1*sinhm2 - 405*bc*sin5n1*sinhm2 - 
       1875*bc*sin3n2*sinhm2 + 405*bc*sin5n2*sinhm2 + 
       1875*b1*cosn1*sinh3m2 - 3750*b2*mu0*cosn1*sinh3m2 + 
       3750*b2*mu2*cosn1*sinh3m2 - 1875*bc*nu0*cosn1*sinh3m2 + 
       1875*bc*nu1*cosn1*sinh3m2 - 75*b1*cos5n1*sinh3m2 + 
       150*b2*mu0*cos5n1*sinh3m2 - 150*b2*mu2*cos5n1*sinh3m2 + 
       75*bc*nu0*cos5n1*sinh3m2 - 75*bc*nu1*cos5n1*sinh3m2 - 
       1875*b1*cosn2*sinh3m2 + 3750*b2*mu0*cosn2*sinh3m2 - 
       3750*b2*mu2*cosn2*sinh3m2 + 1875*bc*nu0*cosn2*sinh3m2 - 
       1875*bc*nu1*cosn2*sinh3m2 + 75*b1*cos5n2*sinh3m2 - 
       150*b2*mu0*cos5n2*sinh3m2 + 150*b2*mu2*cos5n2*sinh3m2 - 
       75*bc*nu0*cos5n2*sinh3m2 + 75*bc*nu1*cos5n2*sinh3m2 - 
       1875*bc*sinn1*sinh3m2 + 15*bc*sin5n1*sinh3m2 + 
       1875*bc*sinn2*sinh3m2 - 15*bc*sin5n2*sinh3m2 + 
       15*(-27*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cosn1 + 3*b1*cos3n1 - 
          6*b2*mu0*cos3n1 + 6*b2*mu2*cos3n1 - 3*bc*nu0*cos3n1 + 
          3*bc*nu1*cos3n1 + 27*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cosn2 - 
          3*(b1 - 2*b2*mu0 + 2*b2*mu2 - bc*nu0 + bc*nu1)*cos3n2 + 27*bc*sinn1 - 
          bc*sin3n1 - 27*bc*sinn2 + bc*sin3n2)*sinh5m2);

  return val;
}

//CPMZ marked: nu2 bug in convert function
double second_order_r1fdV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc)
{
  double f = 1./108.;
  double pa3 = pow(a,3.);
  double pm012 = pow(mu0-mu1,2.);
  double pm022 = pow(mu0-mu2,2.);
  double pn012 = pow(nu0-nu1,2.);
  double pn022 = pow(nu0-nu2,2.);
  double cm01n01 = c1 - bc*mu0 + bc*mu1 - 2.*c2*nu0 + 2.*c2*nu1;
  double cm01n02 = c1 - bc*mu0 + bc*mu1 - 2.*c2*nu0 + 2.*c2*nu2;
  double cm02n01 = c1 - bc*mu0 + bc*mu2 - 2.*c2*nu0 + 2.*c2*nu1;
  double bm02n01 = b1 - 2.*b2*mu0 + 2.*b2*mu2 - bc*nu0 + bc*nu1;
  double bm01n01 = b1 - 2.*b2*mu0 + 2.*b2*mu1 - bc*nu0 + bc*nu1;
  double bm01n02 = b1 - 2.*b2*mu0 + 2.*b2*mu1 - bc*nu0 + bc*nu2;
  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);
  double sin3n1 = sin(3.*nu1);
  double sin3n2 = sin(3.*nu2);

  double val = 0.;

  return val;
}

//CPMZ marked: nu2 bug in convert function
double second_order_fdV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a, double a0, double b1, double b2, double c1, double c2, double bc)
{
  //b1 = b2 = 0.;
  //c1 = c2 = 0.;
  //bc = 0.;

  double f = 1./108.;
  double pa3 = pow(a,3.);
  double pm012 = pow(mu0-mu1,2.);
  double pm022 = pow(mu0-mu2,2.);
  double pn012 = pow(nu0-nu1,2.);
  double pn022 = pow(nu0-nu2,2.);
  double cm01n01 = c1 - bc*mu0 + bc*mu1 - 2.*c2*nu0 + 2.*c2*nu1;
  double cm01n02 = c1 - bc*mu0 + bc*mu1 - 2.*c2*nu0 + 2.*c2*nu2;
  double cm02n01 = c1 - bc*mu0 + bc*mu2 - 2.*c2*nu0 + 2.*c2*nu1;
  double bm02n01 = b1 - 2.*b2*mu0 + 2.*b2*mu2 - bc*nu0 + bc*nu1;
  double bm01n01 = b1 - 2.*b2*mu0 + 2.*b2*mu1 - bc*nu0 + bc*nu1;
  double bm01n02 = b1 - 2.*b2*mu0 + 2.*b2*mu1 - bc*nu0 + bc*nu2;
  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);
  double sin3n1 = sin(3.*nu1);
  double sin3n2 = sin(3.*nu2);

  double val = f*(pa3*(-(cosh3m1*((9.*a0 + b2*(2. + 9.*pm012) + 
               9.*c2*(-2. + pn012) + 9.*c1*(-nu0 + nu1) - 
               9.*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cosn1 - 
            (9.*a0 + b2*(2. + 9.*pm012) + 9.*c2*(-2. + pn022) + 
               9.*c1*(-nu0 + nu2) - 9.*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cosn2 - 
            9.*(cm01n01)*sinn1 + 
            9.*(cm01n02)*sinn2)) + 
       cosh3m2*((9.*a0 + b2*(2. + 9.*pm022) + 
             9.*c2*(-2. + pn012) + 9.*c1*(-nu0 + nu1) - 
             9.*(mu0 - mu2)*(b1 + bc*(-nu0 + nu1)))*cosn1 - 
          (9.*a0 + b2*(2. + 9.*pm022) + 9.*c2*(-2. + pn022) + 
             9.*c1*(-nu0 + nu2) - 9.*(mu0 - mu2)*(b1 + bc*(-nu0 + nu2)))*cosn2 - 
          9.*(c1 - bc*mu0 + bc*mu2 - 2.*c2*nu0 + 2.*c2*nu1)*sinn1 + 
          9.*(c1 - bc*mu0 + bc*mu2 - 2.*c2*nu0 + 2.*c2*nu2)*sinn2) - 
       coshm1*(-((9.*a0 + 9.*b2*(2. + pm012) + 
               c2*(-2. + 9.*pn012) + 9.*c1*(-nu0 + nu1) - 
               9.*(mu0 - mu1)*(b1 + bc*(-nu0 + nu1)))*cos3n1) + 
          (9.*a0 + 9.*b2*(2. + pm012) + c2*(-2. + 9.*pn022) + 
             9.*c1*(-nu0 + nu2) - 9.*(mu0 - mu1)*(b1 + bc*(-nu0 + nu2)))*cos3n2 + 
          3.*(cm01n01)*sin3n1 - 
          3.*(cm01n02)*sin3n2) + 
       coshm2*(-((9.*a0 + 9.*b2*(2. + pm022) + 
               c2*(-2. + 9.*pn012) + 9.*c1*(-nu0 + nu1) - 
               9.*(mu0 - mu2)*(b1 + bc*(-nu0 + nu1)))*cos3n1) + 
          (9.*a0 + 9.*b2*(2. + pm022) + c2*(-2. + 9.*pn022) + 
             9.*c1*(-nu0 + nu2) - 9.*(mu0 - mu2)*(b1 + bc*(-nu0 + nu2)))*cos3n2 + 
          3.*(c1 - bc*mu0 + bc*mu2 - 2.*c2*nu0 + 2.*c2*nu1)*sin3n1 - 
          3.*(c1 - bc*mu0 + bc*mu2 - 2.*c2*nu0 + 2.*c2*nu2)*sin3n2) - 
       3.*((3.*(bm01n01)*cos3n1 - 
             3.*(bm01n02)*cos3n2 + 
             bc*(-sin3n1 + sin3n2))*sinhm1 + 
          (-((bm01n01)*cosn1) + 
             (bm01n02)*cosn2 + bc*(sinn1 - sinn2))
            *sinh3m1) + 3.*((3.*(b1 - 2.*b2*mu0 + 2.*b2*mu2 - bc*nu0 + bc*nu1)*cos3n1 - 
             3.*(b1 - 2.*b2*mu0 + 2.*b2*mu2 - bc*nu0 + bc*nu2)*cos3n2 + 
             bc*(-sin3n1 + sin3n2))*sinhm2 + 
          (-((b1 - 2.*b2*mu0 + 2.*b2*mu2 - bc*nu0 + bc*nu1)*cosn1) + 
             (b1 - 2.*b2*mu0 + 2.*b2*mu2 - bc*nu0 + bc*nu2)*cosn2 + bc*(sinn1 - sinn2))
            *sinh3m2)));

  return val;
}

#pragma acc routine seq
double first_order_frdV(double Z1, double Z2, double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double phi0, double phi1, double phi2, double a, double b, double c, double d)
{
  double f0 = 0.125;

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double coshm12 = coshm1*coshm1;
  double coshm22 = coshm2*coshm2;
  double sinhm12 = sinhm1*sinhm1;
  double sinhm22 = sinhm2*sinhm2;

  double sinh2m1 = sinh(2.*mu1);
  double sinh2m2 = sinh(2.*mu2);
  double cosh2m1 = cosh(2.*mu1);
  double cosh2m2 = cosh(2.*mu2);

  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double cos2n1 = cos(2.*nu1);
  double cos2n2 = cos(2.*nu2);

  double sinn1n2 = sin(nu1)-sin(nu2);
  double sinn2n1 = -sinn1n2;
  double cosn2n1 = cos(nu2)-cos(nu1);
  double cosn1n2 = -cosn2n1; 
  double sin2n1n2 = sin(2.*nu1)-sin(2.*nu2);
  double sin2n2n1 = -sin2n1n2;
  double cos2n1n2 = cos(2.*nu1)-cos(2.*nu2);
  double cos2n2n1 = -cos2n1n2;

  double abm1 = 2.*a + 2.*b*(mu1-mu0);
  double abm2 = 2.*a + 2.*b*(mu2-mu0);
  double cn1 = 2.*c*(nu1-nu0);
  double cn2 = 2.*c*(nu2-nu0);
  double dp12 = d*(phi1+phi2-2.*phi0);

  double abm1n1p12 = abm1 + cn1 + dp12;
  double abm2n1p12 = abm2 + cn1 + dp12;
  double abm1n2p12 = abm1 + cn2 + dp12;
  double abm2n2p12 = abm2 + cn2 + dp12;
  double abmn1p12 = 2.*a - 2.*b*mu0 + cn1 + dp12;
  double abmn2p12 = 2.*a - 2.*b*mu0 + cn2 + dp12;

  double val1 = 2.*coshm12*(abmn1p12*cosn1 - abmn2p12*cosn2 + 2.*c*sinn2n1);
  double val2 = -2.*coshm22*(abmn1p12*cosn1 - abmn2p12*cosn2 + 2.*c*sinn2n1);
  double val3 = -coshm1*(abm1n1p12*cos2n1-abm1n2p12*cos2n2 + c*sin2n2n1);
  double val4 = coshm2*(abm2n1p12*cos2n1-abm2n2p12*cos2n2 + c*sin2n2n1);

  double val5 =   b*(2.*cos2n1n2*sinhm1 + cosn1n2*(2.*mu1*cosh2m1-sinh2m1));
  double val5b = -b*(2.*cos2n2n1*sinhm1 + cosn1n2*(2.*mu1*cosh2m1-sinh2m1));

  double val6  = b*(2.*cos2n2n1*sinhm2 - cosn1n2*(2.*mu2*cosh2m2-sinh2m2));
  double val6b = b*(2.*cos2n2n1*sinhm2 + cosn1n2*(2.*mu2*cosh2m2-sinh2m2));

  double vt1 = f0*Z1*(phi1-phi2)*(-val1 - val2 + val3 + val4 + val5b + val6b);
  double vt2 = -f0*Z2*(phi1-phi2)*(val1 + val2 + val3 + val4 + val5 + val6);

  return vt1 + vt2;
}

#pragma acc routine seq
double first_order_frdV(double Z1, double Z2, double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a0, double b1, double c1)
{
  double fp = 0.125;

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double coshm12 = coshm1*coshm1;
  double coshm22 = coshm2*coshm2;
  double sinhm12 = sinhm1*sinhm1;
  double sinhm22 = sinhm2*sinhm2;

  double sinh2m1 = sinh(2.*mu1);
  double sinh2m2 = sinh(2.*mu2);
  double cosh2m1 = cosh(2.*mu1);
  double cosh2m2 = cosh(2.*mu2);

  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double cosn12 = cosn1*cosn1;
  double cosn22 = cosn2*cosn2;

  double sin2n1 = sin(2.*nu1);
  double sin2n2 = sin(2.*nu2);
  double cos2n1 = cos(2.*nu1);
  double cos2n2 = cos(2.*nu2);

  double val1 = fp*(2*(a0 + b1*(-mu0 + mu1) + c1*(-nu0 + nu2))*cos2n2*coshm1 + 2*a0*cos2n1*coshm2 - 2*b1*mu0*cos2n1*coshm2 + 
       2*b1*mu2*cos2n1*coshm2 - 2*c1*nu0*cos2n1*coshm2 + 2*c1*nu1*cos2n1*coshm2 - 2*a0*cos2n2*coshm2 + 2*b1*mu0*cos2n2*coshm2 - 
       2*b1*mu2*cos2n2*coshm2 + 2*c1*nu0*cos2n2*coshm2 - 2*c1*nu2*cos2n2*coshm2 + 2*a0*cosn1*cosh2m2 - 2*b1*mu0*cosn1*cosh2m2 + 
       2*b1*mu2*cosn1*cosh2m2 - 2*c1*nu0*cosn1*cosh2m2 + 2*c1*nu1*cosn1*cosh2m2 - 2*a0*cosn2*cosh2m2 + 2*b1*mu0*cosn2*cosh2m2 - 
       2*b1*mu2*cosn2*cosh2m2 + 2*c1*nu0*cosn2*cosh2m2 - 2*c1*nu2*cosn2*cosh2m2 - 2*c1*cosh2m2*sinn1 - c1*coshm2*sin2n1 + 
       coshm1*(-2*(a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu1)*cos2n1 + c1*sin2n1) + 2*c1*cosh2m2*sinn2 - 
       2*cosh2m1*((a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu1)*cosn1 - (a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu2)*cosn2 + c1*(-sinn1 + sinn2)) - 
       c1*coshm1*sin2n2 + c1*coshm2*sin2n2 + 4*b1*cosn12*sinhm1 - 4*b1*cosn22*sinhm1 - 2*b1*cosn2*coshm1*sinhm1 + 
       b1*cosn1*sinh2m1 - 4*b1*cosn12*sinhm2 + 4*b1*cosn22*sinhm2 + b1*(-cosn1 + cosn2)*sinh2m2);

  double val2 = fp*(-2*a0*cos2n1*coshm2 + 2*b1*mu0*cos2n1*coshm2 - 2*b1*mu2*cos2n1*coshm2 + 2*c1*nu0*cos2n1*coshm2 - 
       2*c1*nu1*cos2n1*coshm2 + 2*a0*cos2n2*coshm2 - 2*b1*mu0*cos2n2*coshm2 + 2*b1*mu2*cos2n2*coshm2 - 2*c1*nu0*cos2n2*coshm2 + 
       2*c1*nu2*cos2n2*coshm2 + 2*a0*cosn1*cosh2m2 - 2*b1*mu0*cosn1*cosh2m2 + 2*b1*mu2*cosn1*cosh2m2 - 2*c1*nu0*cosn1*cosh2m2 + 
       2*c1*nu1*cosn1*cosh2m2 - 2*a0*cosn2*cosh2m2 + 2*b1*mu0*cosn2*cosh2m2 - 2*b1*mu2*cosn2*cosh2m2 + 2*c1*nu0*cosn2*cosh2m2 - 
       2*c1*nu2*cosn2*cosh2m2 - 2*c1*cosh2m2*sinn1 + c1*coshm2*sin2n1 + 2*c1*cosh2m2*sinn2 - 
       2*cosh2m1*((a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu1)*cosn1 - (a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu2)*cosn2 + c1*(-sinn1 + sinn2)) - 
       c1*coshm2*sin2n2 + coshm1*(2*(a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu1)*cos2n1 - 2*(a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu2)*cos2n2 + 
          c1*(-sin2n1 + sin2n2)) - 4*b1*cosn12*sinhm1 + 4*b1*cosn22*sinhm1 + b1*cosn1*sinh2m1 - 
       b1*cosn2*sinh2m1 + 4*b1*cosn12*sinhm2 - 4*b1*cosn22*sinhm2 + b1*(-cosn1 + cosn2)*sinh2m2);

  return -(Z1*val1 + Z2*val2);
}

#pragma acc routine seq
double first_order_frdV(double Z1, double Z2, double mu0, double mu1, double mu2, double nu1, double nu2, double a, double b)
{
  double cosnpn = cos(nu1) + cos(nu2);
  double cosnmn = cos(nu1) - cos(nu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double cosh2m1 = cosh(2.*mu1);
  double cosh2m2 = cosh(2.*mu2);
  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);

  double abm1m0 = a+b*(mu1-mu0);
  double abm2m0 = a+b*(mu2-mu0);

 //nucleus 1 + 2 (first at z0, other at -z0)
  double val1 = -0.25*Z1*cosnmn*(-2.*abm1m0*cosnpn*coshm1 - abm1m0*cosh2m1 + 2.*abm2m0*cosnpn*coshm2 + abm2m0*cosh2m2
                    + b*(2.*cosnpn+coshm1)*sinhm1 - b*(2.*cosnpn + coshm2)*sinhm2);
  double val2 =  0.25*Z2*cosnmn*(-2.*abm1m0*cosnpn*coshm1 + abm1m0*cosh2m1 + 2.*abm2m0*cosnpn*coshm2 - abm2m0*cosh2m2 
                    - b*(-2.*cosnpn+coshm1)*sinhm1 + b*(-2.*cosnpn + coshm2)*sinhm2);

  return val1+val2;
}

//last 4 values are the Taylor expansion coefficients for f
//analytic integral over mu,nu,phi to 1st order
#pragma acc routine seq
double first_order_fdV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double phi0, double phi1, double phi2, double a, double b, double c, double d)
{
  double f0 = 1./72.;

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double coshm13 = coshm1*coshm1*coshm1;
  double coshm23 = coshm2*coshm2*coshm2;
  double sinhm13 = sinhm1*sinhm1*sinhm1;
  double sinhm23 = sinhm2*sinhm2*sinhm2;

  double sinh2m1 = sinh(2.*mu1);
  double sinh2m2 = sinh(2.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);

  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);

  double sinn1n2 = sin(nu1)-sin(nu2);
  double sinn2n1 = -sinn1n2;
  double cosn2n1 = cos(nu2)-cos(nu1);
  double sin3n1n2 = sin(3.*nu1)-sin(3.*nu2);
  double sin3n2n1 = -sin3n1n2;

  double abm1 = 2.*a + 2.*b*(mu1-mu0);
  double abm2 = 2.*a + 2.*b*(mu2-mu0);
  double cn1 = 2.*c*(nu1-nu0);
  double cn2 = 2.*c*(nu2-nu0);
  double dp12 = d*(phi1+phi2-2.*phi0);
  double abm1n1p12 = abm1 + cn1 + dp12;
  double abm2n1p12 = abm2 + cn1 + dp12;
  double abm1n2p12 = abm1 + cn2 + dp12;
  double abm2n2p12 = abm2 + cn2 + dp12;

  double val1 = 3.*coshm23*(-abm2n1p12*cosn1 + abm2n2p12*cosn2 + 2.*c*sinn1n2);
  double val2 = 3.*coshm13*(abm1n1p12*cosn1 - abm1n2p12*cosn2 + 2.*c*sinn2n1);
  double val3 = 6.*b*cos3n1*sinhm1 - 6.*b*cos3n2*sinhm1 + 6.*b*cosn2n1*coshm1*coshm1*sinhm1 - 2.*b*cosn1*sinhm13 + 2.*b*cosn2*sinhm13;
  double val4 = coshm1*(-3.*abm1n1p12*cos3n1 + 3.*abm1n2p12*cos3n2 + 2.*c*sin3n1n2 + 9.*(abm1n1p12*cosn1 - abm1n2p12*cosn2)*sinhm1*sinhm1);
  double val5 = -9.*c*sinn1*sinhm1*sinh2m1 + 9.*c*sinn2*sinhm1*sinh2m1 - 6.*b*cos3n1*sinhm2 + 6.*b*cos3n2*sinhm2 - 6.*b*cosn2n1*coshm2*coshm2*sinhm2 + 2.*b*cosn1*sinhm23 - 2.*b*cosn2*sinhm23;
  double val6 = coshm2*(3.*abm2n1p12*cos3n1 - 3.*abm2n2p12*cos3n2 + 2.*c*sin3n2n1 + 9.*(-abm2n1p12*cosn1 + abm2n2p12*cosn2)*sinhm2*sinhm2);
  double val7 = 9.*c*sinn1n2*sinhm2*sinh2m2;

  double val = f0*(phi1-phi2)*(val1 + val2 + val3 + val4 + val5 + val6 + val7);

  return val;
}

//analytic integral over mu,nu to 1st order
#pragma acc routine seq
double first_order_fdV(double mu0, double mu1, double mu2, double nu0, double nu1, double nu2, double a0, double b1, double c1)
{
 //a3 multiplied outside of this ftn
  double fp = 1./36.;

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);

  double sinh2m1 = sinh(2.*mu1);
  double sinh2m2 = sinh(2.*mu2);
  double sinh3m1 = sinh(3.*mu1);
  double sinh3m2 = sinh(3.*mu2);
  double cosh3m1 = cosh(3.*mu1);
  double cosh3m2 = cosh(3.*mu2);

  double sinn1 = sin(nu1);
  double sinn2 = sin(nu2);
  double cosn1 = cos(nu1);
  double cosn2 = cos(nu2);
  double sin3n1 = sin(3.*nu1);
  double sin3n2 = sin(3.*nu2);
  double cos3n1 = cos(3.*nu1);
  double cos3n2 = cos(3.*nu2);

  double ta0 = 3.*a0;
  double tb1 = 3.*b1;
  double tc1 = 3.*c1;
  double abm01n01 = a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu1;

  double val = fp*(-ta0*cos3n1*coshm2 + tb1*mu0*cos3n1*coshm2 - tb1*mu2*cos3n1*coshm2 + tc1*nu0*cos3n1*coshm2 - 
       tc1*nu1*cos3n1*coshm2 + ta0*cos3n2*coshm2 - tb1*mu0*cos3n2*coshm2 + tb1*mu2*cos3n2*coshm2 - 
       tc1*nu0*cos3n2*coshm2 + tc1*nu2*cos3n2*coshm2 + ta0*cosn1*cosh3m2 - tb1*mu0*cosn1*cosh3m2 + 
       tb1*mu2*cosn1*cosh3m2 - tc1*nu0*cosn1*cosh3m2 + tc1*nu1*cosn1*cosh3m2 - ta0*cosn2*cosh3m2 + 
       tb1*mu0*cosn2*cosh3m2 - tb1*mu2*cosn2*cosh3m2 + tc1*nu0*cosn2*cosh3m2 - tc1*nu2*cosn2*cosh3m2 - 
       tc1*cosh3m2*sinn1 + c1*coshm2*sin3n1 + tc1*cosh3m2*sinn2 - 
       3*cosh3m1*((a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu1)*cosn1 - (a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu2)*cosn2 + 
          c1*(-sinn1 + sinn2)) - c1*coshm2*sin3n2 + 
       coshm1*(3*(a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu1)*cos3n1 - 3*(a0 - b1*mu0 + b1*mu1 - c1*nu0 + c1*nu2)*cos3n2 + 
          c1*(-sin3n1 + sin3n2)) - tb1*cos3n1*sinhm1 + tb1*cos3n2*sinhm1 + b1*cosn1*sinh3m1 - 
       b1*cosn2*sinh3m1 + tb1*cos3n1*sinhm2 - tb1*cos3n2*sinhm2 + b1*(-cosn1 + cosn2)*sinh3m2);

  return val;
}

//analytic integral over mu to 1st order
#pragma acc routine seq
double first_order_fdV(double mu0, double mu1, double mu2, double nu1, double nu2, double a, double b)
{
  double p = (cos(nu2)-cos(nu1))/36.;

  double mu13 = 3.*mu1;
  double mu23 = 3.*mu2;
  double tabm1m0 = 3.*(a+b*(mu1-mu0));
  double tabm2m0 = 3.*(a+b*(mu2-mu0));

  double sinhm1 = sinh(mu1);
  double sinhm2 = sinh(mu2);
  double coshm1 = cosh(mu1);
  double coshm2 = cosh(mu2);
  double sinh3m1 = sinh(mu13);
  double sinh3m2 = sinh(mu23);
  double cosh3m1 = cosh(mu13);
  double cosh3m2 = cosh(mu23);
  double cosn1mn2 = cos(nu1-nu2);
  double cosn1pn2 = cos(nu1+nu2);

  double opcn = 1.+2.*cosn1pn2;
  double omcn = 1.+2.*cosn1mn2;
  double opmcn = opcn*omcn;

  double val = p * ( -tabm1m0*opmcn*coshm1 + tabm1m0*cosh3m1 + tabm2m0*opmcn*coshm2 - tabm2m0*cosh3m2 
                   - b*(-3.*opmcn*sinhm1 + sinh3m1) + b*(-3.*opmcn*sinhm2+sinh3m2) );
  return val;
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


void get_ab_r1r2(int gs, const double z0, const double Z1, const double Z2, double* val1, double* val2, double* val1p, double* val2p, double* grid, double* gridm, double* wt, double* val)
{
  const double a2 = z0*z0; //a^3/a, where 1/r1 and 1/r2 have 1/a factors in them
  const int gs3 = 3*gs; const int gs6 = 6*gs;

 #pragma acc parallel loop present(val[0:gs],val1[0:gs],val2[0:gs],val1p[0:gs3],val2p[0:gs3],grid[0:gs6],gridm[0:gs6])
  for (int j=0;j<gs;j++)
  {
    double mu = gridm[6*j+0];  double nu = gridm[6*j+1];  double phi = gridm[6*j+2];
    double dmu = gridm[6*j+3]; double dnu = gridm[6*j+4]; double dphi = gridm[6*j+5];

    double mu1  = mu-0.5*dmu;   double mu2  = mu+0.5*dmu; 
    double nu1  = nu-0.5*dnu;   double nu2  = nu+0.5*dnu;
    double phi1 = phi-0.5*dphi; double phi2 = phi+0.5*dphi;

    double sinhm = sinh(mu); double coshm = cosh(mu);
    double sinn  = sin(nu);  double cosn  = cos(nu);
    double sinp  = sin(phi); double cosp  = cos(phi);

    double dxdmu = z0*sinn*cosp*coshm;
    double dydmu = z0*sinn*sinp*coshm;
    double dzdmu = z0*cosn*sinhm;

    double dxdnu =  z0*sinhm*cosp*cosn;
    double dydnu =  z0*sinhm*sinp*cosn;
    double dzdnu = -z0*coshm*sinn;

    double dxdphi = -z0*sinhm*sinn*sinp; double dydphi = z0*sinhm*sinn*cosp; double dzdphi = 0.;

    double v1 = val1[j]; double v2 = val2[j];
    double dfdx1 = val1p[3*j+0]; double dfdy1 = val1p[3*j+1]; double dfdz1 = val1p[3*j+2];
    double dfdx2 = val2p[3*j+0]; double dfdy2 = val2p[3*j+1]; double dfdz2 = val2p[3*j+2];
    dfdx1 *= v2; dfdy1 *= v2; dfdz1 *= v2;
    dfdx2 *= v1; dfdy2 *= v1; dfdz2 *= v1;
    double dfx = dfdx1+dfdx2; double dfy = dfdy1+dfdy2; double dfz = dfdz1+dfdz2;
    double dfdmu  = dfx*dxdmu + dfy*dydmu + dfz*dzdmu;
    double dfdnu  = dfx*dxdnu + dfy*dydnu + dfz*dzdnu;
    //double dfdphi = dfx*dxdphi + dfy*dydphi;

    double a = v1*v2;
    double b = dfdmu;
    double c = dfdnu;
    //double d = dfdphi;

    //val[j] = a2*first_order_frdV(Z1,Z2,mu,mu1,mu2,nu1,nu2,a,b)*idp;
    val[j] = a2*first_order_frdV(Z1,Z2,mu,mu1,mu2,nu,nu1,nu2,a,b,c);
    //val[j] = -a2*first_order_frdV(Z1,Z2,mu,mu1,mu2,nu,nu1,nu2,phi,1.,0.,a,b,c,d)*idp;
    //printf("  val/n: %11.8f %11.8f diff: %5.1e \n",val[j],vn,fabs(val[j]-vn));
  }
  return;
}

void get_ab_2d(int gs, const double z0, double* val1, double* val2, double* val1p, double* val2p, double* grid, double* gridm, double* wt, double* val)
{
 //note: wt not needed

  const double a3 = z0*z0*z0;
  const int gs3 = 3*gs; const int gs6 = 6*gs;
 #pragma acc parallel loop present(val[0:gs],val1[0:gs],val2[0:gs],val1p[0:gs3],val2p[0:gs3],grid[0:gs6],gridm[0:gs6])
  for (int j=0;j<gs;j++)
  {
    double mu = gridm[6*j+0];  double nu = gridm[6*j+1];  double phi = gridm[6*j+2];
    double dmu = gridm[6*j+3]; double dnu = gridm[6*j+4]; double dphi = gridm[6*j+5];

    double mu1  = mu-0.5*dmu;   double mu2  = mu+0.5*dmu; 
    double nu1  = nu-0.5*dnu;   double nu2  = nu+0.5*dnu;
    double phi1 = phi-0.5*dphi; double phi2 = phi+0.5*dphi;

    double sinhm = sinh(mu); double coshm = cosh(mu);
    double sinn  = sin(nu);  double cosn  = cos(nu);
    double sinp  = sin(phi); double cosp  = cos(phi);

    double dxdmu = z0*sinn*cosp*coshm;
    double dydmu = z0*sinn*sinp*coshm;
    double dzdmu = z0*cosn*sinhm;

    double dxdnu =  z0*sinhm*cosp*cosn;
    double dydnu =  z0*sinhm*sinp*cosn;
    double dzdnu = -z0*coshm*sinn;

    double dxdphi = -z0*sinhm*sinn*sinp; double dydphi = z0*sinhm*sinn*cosp; double dzdphi = 0.;

   //could calc val1p starting from val2, vice versa
    double v1 = val1[j]; double v2 = val2[j];
    double dfdx1 = val1p[3*j+0]; double dfdy1 = val1p[3*j+1]; double dfdz1 = val1p[3*j+2];
    double dfdx2 = val2p[3*j+0]; double dfdy2 = val2p[3*j+1]; double dfdz2 = val2p[3*j+2];
    dfdx1 *= v2; dfdy1 *= v2; dfdz1 *= v2;
    dfdx2 *= v1; dfdy2 *= v1; dfdz2 *= v1;
    double dfx = dfdx1+dfdx2; double dfy = dfdy1+dfdy2; double dfz = dfdz1+dfdz2;
    double dfdmu  = dfx*dxdmu + dfy*dydmu + dfz*dzdmu;
    double dfdnu  = dfx*dxdnu + dfy*dydnu + dfz*dzdnu;
    double dfdphi = dfx*dxdphi + dfy*dydphi;

    double a = v1*v2;
    double b = dfdmu;
    double c = dfdnu;
    //double d = dfdphi;

    //CPMZ debug
    //a = 1.;
    //b = c = 0.;

    double vp = a3*first_order_fdV(mu,mu1,mu2,nu,nu1,nu2,a,b,c);
    //val[j] = a3*first_order_fdV(mu,mu1,mu2,nu,nu1,nu2,phi,phi1,phi2,a,b,c,d);
    //  printf("    vp/vn: %5.1e %5.1e diff: %4.1e \n",vp,val[j],fabs(val[j]-vp));

    val[j] = vp;
  }
  return;
}

void get_ab_3d(int gs, const double z0, double* val1, double* val2, double* val1p, double* val2p, double* grid, double* gridm, double* wt, double* val)
{
  double a3 = z0*z0*z0;
  int gs3 = 3*gs; int gs6 = 6*gs;
 //#pragma acc parallel loop present(val[0:gs],val1[0:gs],val2[0:gs],val1p[0:gs3],val2p[0:gs3],grid[0:gs6],gridm[0:gs6],wt[0:gs])
  for (int j=0;j<gs;j++)
  {
    double mu = gridm[6*j+0];  double nu = gridm[6*j+1];  double phi = gridm[6*j+2];
    double dmu = gridm[6*j+3]; double dnu = gridm[6*j+4]; double dphi = gridm[6*j+5];

    double mu1  = mu-0.5*dmu;   double mu2  = mu+0.5*dmu; 
    double nu1  = nu-0.5*dnu;   double nu2  = nu+0.5*dnu;
    double phi1 = phi-0.5*dphi; double phi2 = phi+0.5*dphi;

    double sinhm = sinh(mu); double coshm = cosh(mu);
    double sinn  = sin(nu);  double cosn  = cos(nu);
    double sinp  = sin(phi); double cosp  = cos(phi);

    double dxdmu = z0*sinn*cosp*coshm;
    double dydmu = z0*sinn*sinp*coshm;
    double dzdmu = z0*cosn*sinhm;

    double dxdnu =  z0*sinhm*cosp*cosn;
    double dydnu =  z0*sinhm*sinp*cosn;
    double dzdnu = -z0*coshm*sinn;

    double dxdphi = -z0*sinhm*sinn*sinp;
    double dydphi =  z0*sinhm*sinn*cosp;
    //double dzdphi = 0.;

   //could calc val1p starting from val2, vice versa
    double v1 = val1[j]; double v2 = val2[j];
    double dfdx1 = val1p[3*j+0]; double dfdy1 = val1p[3*j+1]; double dfdz1 = val1p[3*j+2];
    double dfdx2 = val2p[3*j+0]; double dfdy2 = val2p[3*j+1]; double dfdz2 = val2p[3*j+2];
    dfdx1 *= v2; dfdy1 *= v2; dfdz1 *= v2;
    dfdx2 *= v1; dfdy2 *= v1; dfdz2 *= v1;
    double dfx = dfdx1+dfdx2; double dfy = dfdy1+dfdy2; double dfz = dfdz1+dfdz2;
    double dfdmu  = dfx*dxdmu + dfy*dydmu + dfz*dzdmu;
    double dfdnu  = dfx*dxdnu + dfy*dydnu + dfz*dzdnu;
    double dfdphi = dfx*dxdphi + dfy*dydphi;

    double a = v1*v2;
    double b = dfdmu;
    double c = dfdnu;
    double d = dfdphi;

    //double vp = a3*first_order_fdV(mu,mu1,mu2,nu1,nu2,a,b)*dphi;
    val[j] = a3*first_order_fdV(mu,mu1,mu2,nu,nu1,nu2,phi,phi1,phi2,a,b,c,d);
    //if (fabs(a)<1.e-10)
    //if (fabs(val[j]-vp)>1.e-10)
    //  printf("    vp/vn: %5.1e %5.1e diff: %4.1e \n",vp,val[j],fabs(val[j]-vp));
  }
  return;
}




//1s/2p(partial) ftns only for the moment
void get_ab_mnp(int gs, const double a, bool do_r12, double Z1, double Z2, int n1, int l1, int m1, int n2, int l2, int m2, double norm1, double norm2, double zt1, double zt2, double* grid, double* gridm, double* wt, double* val)
{
  printf("\n ERROR: need to rework the mnp functions (nu2 bug in convert) \n"); exit(-1);

  int gs6 = 6*gs;

  double n12 = norm1*norm2;
  double zt12 = zt1+zt2;
  double zm12 = zt1-zt2;

  int nr1 = n1-l1-1;
  int nr2 = n2-l2-1;

  int nx = 0;
  int ny = 0;
  int nz1 = 0;
  int nz2 = 0;

  if (l1==1)
  {
    if (m1==-1) ny = 1;
    if (m1== 0) nz1 = 1;
    if (m1== 1) nx = 1;
  }
  if (l2==1)
  {
    if (m2==-1) ny += 1;
    if (m2== 0) nz2  = 1;
    if (m2== 1) nx += 1;
  }
  int nxy = nx + ny;
  int nz = nz1 + nz2;

  if (n1-l1-1>0) { printf("\n WARNING: get_ab_mnp not ready \n"); exit(-1); }
  if (n2-l2-1>0) { printf("\n WARNING: get_ab_mnp not ready \n"); exit(-1); }

  int wz = 0;
  if (nz1>0) wz = 1;
  else if (nz2>0) wz = -1;

  printf("    nx/y: %i %i  z1/2: %i %i \n",nx,ny,nz1,nz2);

  if (nx%2==1 || ny%2==1)
  {
    printf("  odd x/y: set integral to zero \n");
   //odd in x or y: goes to zero
    for (int j=0;j<gs;j++)
      val[j] = 0.;
    return;
  }

  //double idp = get_idp(nx,ny);
  double idp = get_idp(l1,m1,l2,m2);

 //eventually move if statements out of the loop
 //#pragma acc parallel loop present(val[0:gs],grid[0:gs6],gridm[0:gs6],wt[0:gs])
  for (int j=0;j<gs;j++)
  {
    double mu = gridm[6*j+0];  double nu = gridm[6*j+1];  double phi = gridm[6*j+2];
    double dmu = gridm[6*j+3]; double dnu = gridm[6*j+4]; double dphi = gridm[6*j+5];

    double mu1  = mu-0.5*dmu;   double mu2  = mu+0.5*dmu; 
    double nu1  = nu-0.5*dnu;   double nu2  = nu+0.5*dnu;
    double phi1 = phi-0.5*dphi; double phi2 = phi+0.5*dphi;

    double sinhm = sinh(mu); double coshm = cosh(mu);
    double sinn  = sin(nu);  double cosn  = cos(nu);
    //double sinp  = sin(phi); double cosp  = cos(phi);
 
    double r1 = a*(coshm-cosn); double r2 = a*(coshm+cosn);
    double f1 = exp(-zt1*r1);   double f2 = exp(-zt2*r2);
    double f12 = f1*f2;

    double dfdmu = a*zt12*sinhm;
    double dfdnu = a*zm12*sinn;
    double df2dmu = a*zt12*(coshm+dfdmu*sinhm);
    double df2dnu = a*zm12*(cosn+dfdnu*sinn);
    double dfdmn = dfdmu*dfdnu;

    dfdmu  *= f12;
    dfdnu  *= f12;
    df2dmu *= f12;
    df2dnu *= f12;
    dfdmn  *= f12;

    //double v0 = val[j];

    if (do_r12)
    {
      if (nxy==0 && nz==0)
        val[j] = n12*second_order_fordV(Z1,Z2,mu,mu1,mu2,nu,nu1,nu2,a,f12,dfdmu,df2dmu,dfdnu,df2dnu,dfdmn)*idp;
      else if (nx==2 || ny==2)
        val[j] = n12*second_order_fx2ordV(Z1,Z2,mu,mu1,mu2,nu,nu1,nu2,a,f12,dfdmu,df2dmu,dfdnu,df2dnu,dfdmn)*idp;
      else if ((nz1==1 || nz2==1) && nz<2)
        val[j] = n12*second_order_fzordV(Z1,Z2,mu,mu1,mu2,nu,nu1,nu2,a,f12,dfdmu,df2dmu,dfdnu,df2dnu,dfdmn)*idp;
      else if (nz1==1 && nz2==1)
        val[j] = n12*second_order_fzzordV(Z1,Z2,mu,mu1,mu2,nu,nu1,nu2,a,f12,dfdmu,df2dmu,dfdnu,df2dnu,dfdmn)*idp;
    }
    else
    {
      if (nxy==0 && nz==0)
      {
        if (nr1==1)
          val[j] = n12*second_order_r1fdV(mu,mu1,mu2,nu,nu1,nu2,a,f12,dfdmu,df2dmu,dfdnu,df2dnu,dfdmn)*idp;
        //else if (nr2==1)
        //  val[j] = n12*second_order_r2fdV(mu,mu1,mu2,nu,nu1,nu2,a,f12,dfdmu,df2dmu,dfdnu,df2dnu,dfdmn)*idp;
        else
          val[j] = n12*second_order_fdV(mu,mu1,mu2,nu,nu1,nu2,a,f12,dfdmu,df2dmu,dfdnu,df2dnu,dfdmn)*idp;
      }
      else if (nx==2 || ny==2)
        val[j] = n12*second_order_fx2dV(mu,mu1,mu2,nu,nu1,nu2,a,f12,dfdmu,df2dmu,dfdnu,df2dnu,dfdmn)*idp;
      else if ((nz1==1 || nz2==1) && nz<2)
        val[j] = n12*second_order_fzdV(mu,mu1,mu2,nu,nu1,nu2,a,wz,f12,dfdmu,df2dmu,dfdnu,df2dnu,dfdmn)*idp;
      else if (nz1==1 && nz2==1)
        val[j] = n12*second_order_fzzdV(mu,mu1,mu2,nu,nu1,nu2,a,f12,dfdmu,df2dmu,dfdnu,df2dnu,dfdmn)*idp;
    }

    //printf("   val/0: %11.8f %11.8f  diff: %11.8f \n",val[j],v0,fabs(val[j]-v0));
  }

  return;
}
