#ifndef _FP_DEF_H_
#define _FP_DEF_H_

#ifndef EVL64
#define EVL64 1 // eval grid points in fp64 if true, else fp32
#endif
#define RED64 1 // reduce grid in fp64 if true, else fp32

// Changed by Soumi
// ----------------
#if EVL64
#define FP1 float
#else
#define FP1 float
#endif

#if RED64
#define FP2 double
#else
#define FP2 double
#endif
// -----------------
/*
#if EVL64
#define FP1 double
#else
#define FP1 float
#endif

#if RED64
#define FP2 double
#else
#define FP2 float
#endif
*/


#endif // _FP_DEF_H_
