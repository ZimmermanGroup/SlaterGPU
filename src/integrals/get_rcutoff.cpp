#include "get_rcutoff.h"
#include "fp_def.h"

FP1 get_rcutoff(int n, int l, FP1 z) {
  switch(n) {
    case 1 :
      switch(l) {
        case 0 :
          return  7.000000000000000146e-03/z;
        default:
          return 0.;
      } // switch l
    case 2 :
      switch(l) {
        case 0 :
          return  1.200000000000000025e-02/z;
        case 1 :
          return  2.999999999999999889e-01/z;
        default:
          return 0.;
      } // switch l
    case 3 :
      switch(l) {
        case 0 :
          return  1.400000000000000029e-02/z;
        case 1 :
          return  2.999999999999999889e-01/z;
        case 2 :
          return  2.999999999999999889e-01/z;
        default:
          return 0.;
      } // switch l
    case 4 :
      switch(l) {
        case 0 :
          return  2.000000000000000042e-02/z;
        case 1 :
          return  2.999999999999999889e-01/z;
        case 2 :
          return  5.000000000000000000e-01/z;
        case 3 :
          return  8.000000000000000444e-01/z;
        default:
          return 0.;
      } // switch l
    case 5 :
      switch(l) {
        case 0 :
          return  2.300000000000000100e-01/z;
        case 1 :
          return  5.500000000000000444e-01/z;
        case 2 :
          return  6.999999999999999556e-01/z;
        case 3 :
          return  7.500000000000000000e-01/z;
        case 4 :
          return  1.000000000000000000e+00/z;
        default:
          return 0.;
      } // switch l
    case 6 :
      switch(l) {
        case 0 :
          return  2.500000000000000000e-01/z;
        case 1 :
          return  8.000000000000000444e-01/z;
        case 2 :
          return  9.000000000000000222e-01/z;
        case 3 :
          return  1.000000000000000000e+00/z;
        case 4 :
          return  1.500000000000000000e+00/z;
        case 5 :
          return  2.000000000000000000e+00/z;
        default:
          return 0.;
      } // switch l
    case 7 :
      switch(l) {
        case 0 :
          return  2.999999999999999889e-02/z;
        case 1 :
          return  9.000000000000000222e-01/z;
        case 2 :
          return  1.199999999999999956e+00/z;
        case 3 :
          return  1.399999999999999911e+00/z;
        case 4 :
          return  1.800000000000000044e+00/z;
        case 5 :
          return  2.399999999999999911e+00/z;
        case 6 :
          return  2.700000000000000178e+00/z;
        default:
          return 0.;
      } // switch l
    case 8 :
      switch(l) {
        case 0 :
          return  2.999999999999999889e-02/z;
        case 1 :
          return  1.199999999999999956e+00/z;
        case 2 :
          return  2.000000000000000000e+00/z;
        case 3 :
          return  2.000000000000000000e+00/z;
        case 4 :
          return  2.700000000000000178e+00/z;
        case 5 :
          return  3.000000000000000000e+00/z;
        case 6 :
          return  3.100000000000000089e+00/z;
        case 7 :
          return  3.200000000000000178e+00/z;
        default:
          return 0.;
      } // switch l
    case 9 :
      switch(l) {
        case 0 :
          return  2.999999999999999889e-02/z;
        case 1 :
          return  1.399999999999999911e+00/z;
        case 2 :
          return  1.699999999999999956e+00/z;
        case 3 :
          return  2.000000000000000000e+00/z;
        case 4 :
          return  2.500000000000000000e+00/z;
        case 5 :
          return  3.000000000000000000e+00/z;
        case 6 :
          return  3.100000000000000089e+00/z;
        case 7 :
          return  4.000000000000000000e+00/z;
        case 8 :
          return  5.000000000000000000e+00/z;
        default:
          return 0.;
      } // switch l
    default :
      return 0.;
  } // switch n
}

FP1 get_drcutoff(int n, int l, FP1 z) {
  switch(n) {
    case 1 :
      switch(l) {
        case 0 :
          return  1.499999999999999944e-01/z;
        default:
          return 0.;
      } // switch l
    case 2 :
      switch(l) {
        case 0 :
          return  5.000000000000000000e-01/z;
        case 1 :
          return  5.000000000000000000e-01/z;
        default:
          return 0.;
      } // switch l
    case 3 :
      switch(l) {
        case 0 :
          return  5.000000000000000000e-01/z;
        case 1 :
          return  5.999999999999999778e-01/z;
        case 2 :
          return  6.999999999999999556e-01/z;
        default:
          return 0.;
      } // switch l
    case 4 :
      switch(l) {
        case 0 :
          return  9.000000000000000222e-01/z;
        case 1 :
          return  1.000000000000000000e+00/z;
        case 2 :
          return  1.000000000000000000e+00/z;
        case 3 :
          return  1.199999999999999956e+00/z;
        default:
          return 0.;
      } // switch l
    case 5 :
      switch(l) {
        case 0 :
          return  8.000000000000000444e-01/z;
        case 1 :
          return  1.000000000000000000e+00/z;
        case 2 :
          return  1.500000000000000000e+00/z;
        case 3 :
          return  1.500000000000000000e+00/z;
        case 4 :
          return  1.800000000000000044e+00/z;
        default:
          return 0.;
      } // switch l
    case 6 :
      switch(l) {
        case 0 :
          return  1.800000000000000044e+00/z;
        case 1 :
          return  1.399999999999999911e+00/z;
        case 2 :
          return  1.500000000000000000e+00/z;
        case 3 :
          return  1.899999999999999911e+00/z;
        case 4 :
          return  2.000000000000000000e+00/z;
        case 5 :
          return  2.299999999999999822e+00/z;
        default:
          return 0.;
      } // switch l
    case 7 :
      switch(l) {
        case 0 :
          return  1.800000000000000044e+00/z;
        case 1 :
          return  1.399999999999999911e+00/z;
        case 2 :
          return  1.800000000000000044e+00/z;
        case 3 :
          return  2.299999999999999822e+00/z;
        case 4 :
          return  3.000000000000000000e+00/z;
        case 5 :
          return  3.200000000000000178e+00/z;
        case 6 :
          return  3.500000000000000000e+00/z;
        default:
          return 0.;
      } // switch l
    case 8 :
      switch(l) {
        case 0 :
          return  2.000000000000000000e+00/z;
        case 1 :
          return  1.500000000000000000e+00/z;
        case 2 :
          return  2.500000000000000000e+00/z;
        case 3 :
          return  2.799999999999999822e+00/z;
        case 4 :
          return  3.399999999999999911e+00/z;
        case 5 :
          return  3.799999999999999822e+00/z;
        case 6 :
          return  4.000000000000000000e+00/z;
        case 7 :
          return  4.099999999999999645e+00/z;
        default:
          return 0.;
      } // switch l
    case 9 :
      switch(l) {
        case 0 :
          return  2.500000000000000000e+00/z;
        case 1 :
          return  2.000000000000000000e+00/z;
        case 2 :
          return  2.299999999999999822e+00/z;
        case 3 :
          return  3.399999999999999911e+00/z;
        case 4 :
          return  3.500000000000000000e+00/z;
        case 5 :
          return  4.000000000000000000e+00/z;
        case 6 :
          return  4.000000000000000000e+00/z;
        case 7 :
          return  4.299999999999999822e+00/z;
        case 8 :
          return  4.500000000000000000e+00/z;
        default:
          return 0.;
      } // switch l
    default :
      return 0.;
  } // switch n
}

