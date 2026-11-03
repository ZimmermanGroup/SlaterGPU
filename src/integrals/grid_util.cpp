#include "grid_util.h"
#include "lebedev2.h"
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdlib>

using namespace std;

void compare_pao_12(int ngpu, bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis, cusolverDnHandle_t cu_hdl)
{
  int gs = nrad*nang;
  int gsa = natoms*gs;
  int gsa6 = gsa*6;
  int dim = gsa;

  int N = basis.size();
  int N2 = N*N;

  double* Pao1 = new double[N2];
  double* Pao2 = new double[N2];
  bool found1 = read_square(N,Pao1,"Paoci");
  bool found2 = read_square(N,Pao2,"Paoks");
  if (!found1 || !found2) { printf("\n\n  couldn't find Paoci or Paoks file \n"); return; }

  double* rho1 = new double[gsa];
  double* rho2 = new double[gsa];
  double* grid = new double[gsa6];
  float* gridf = new float[gsa6];
  double* wt = new double[gsa];

  printf("\n gbasis: %i natoms: %i nrad: %3i nang: %4i gsa: %4i \n",(int)gbasis,natoms,nrad,nang,gsa);

  #pragma acc enter data create(rho1[0:gsa],rho2[0:gsa])
  #pragma acc enter data create(grid[0:gsa6],wt[0:gsa],gridf[0:gsa6])
  #pragma acc enter data copyin(coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  printf("\n comparing Pao_ci to Pao_ks \n");

  if (read_int("GRID_SQUARE"))
  {
    int nsq = 1;
    for (int nq=0;nq<1000;nq++)
    {
      int gsize = nq*nq*nq;
      if (gsize>gsa) { nsq = nq-1; break; }
    }
    int gsize = nsq*nsq*nsq;
    printf("\n using square grid (size: %4i) \n",gsize);
    double a = 3;
    double volume = a*a*a;
    double sp = a/(nsq+1);
    double wt1 = volume/gsize;
    for (int j=0;j<gsize;j++)
      wt[j] = wt1;

    int ng = 0;
    int nsq2 = nsq/2;
    for (int i=0;i<nsq;i++)
    {
      double x1 = sp*(i-nsq2);
      for (int j=0;j<nsq;j++)
      {
        double y1 = sp*(j-nsq2);
        for (int k=0;k<nsq;k++)
        {
          double z1 = sp*(k-nsq2);

          grid[6*ng+0] = x1;
          grid[6*ng+1] = y1;
          grid[6*ng+2] = z1;
          ng++;
        }
      }
    }
    for (int i=6*ng;i<gsa6;i++)
      grid[i] = 0.;
    for (int i=ng;i<gsa;i++)
      wt[i] = 0.;
    dim = gsize;

    #pragma acc update device(grid[0:gsa6],wt[0:gs])
  }
  else
    get_becke_grid_full(natoms,atno,coords,nrad,nang,ang_g,ang_w,6,grid,wt);

  if (gbasis)
  {
    for (int j=0;j<gsa6;j++) gridf[j] = grid[j];
    #pragma acc update device(gridf[0:gsa6])

    compute_rho(1,natoms,atno,coords,basis,Pao1,nrad,dim,gridf,rho1,NULL,1);
    compute_rho(1,natoms,atno,coords,basis,Pao2,nrad,dim,gridf,rho2,NULL,1);
  }
  else
  {
    compute_rhod(natoms,atno,coords,basis,Pao1,nrad,dim,grid,rho1,NULL,NULL,1);
    compute_rhod(natoms,atno,coords,basis,Pao2,nrad,dim,grid,rho2,NULL,NULL,1);
  }

  #pragma acc update self(rho1[0:gsa],rho2[0:gsa],grid[0:gsa6],wt[0:gsa])

  if (natoms==1)
  {
    printf("       r        rho1      rho2 \n");
    for (int j=0;j<nrad;j++)
      printf(" %8.5f %12.6e %12.6e\n",grid[6*j*nang+3],rho1[j*nang],rho2[j*nang]);
    printf("\n");
  }

  double diff1 = 0.;
  double diff2 = 0.;
 #pragma acc parallel loop present(rho1[0:gsa],rho2[0:gsa],wt[0:gsa]) reduction(+:diff1,diff2)
  for (int j=0;j<gsa;j++)
  {
    double d12 = rho1[j]-rho2[j];
    diff1 += abs(d12)*wt[j];
    diff2 += d12*d12*wt[j];
  }

  diff2 = sqrt(diff2);
  printf("\n L1/L2: %10.6f %10.6f \n",diff1,diff2);

  // Extra XC100 metrics: real-space shape functions on the SAME quadrature.
  // rho1 = CI, rho2 = CI-derived KS. No grid or density files are written.
  // The #pragma acc update self above already copied all arrays to the host.
  double ne_ci = 0.0, ne_ks = 0.0, norm2_ci = 0.0;
  const double negative_tolerance = 1.e-12; // density, e/bohr^3
  for (int j=0; j<gsa; ++j)
  {
    const double w = wt[j];
    if (w == 0.0) continue;
    if (!std::isfinite(w) || w < 0.0 ||
        !std::isfinite(rho1[j]) || !std::isfinite(rho2[j]) ||
        rho1[j] < -negative_tolerance || rho2[j] < -negative_tolerance)
    {
      fprintf(stderr, "Invalid density/weight at grid point %d: rho_ci=%g rho_ks=%g w=%g\n",
              j,rho1[j],rho2[j],w);
      std::exit(EXIT_FAILURE);
    }
    // Clip ONLY tiny negative roundoff artifacts to zero for entropy metrics.
    const double a = std::max(0.0,rho1[j]);
    const double b = std::max(0.0,rho2[j]);
    ne_ci += w*a;
    ne_ks += w*b;
    norm2_ci += w*a*a;
  }
  if (!(ne_ci > 0.) || !(ne_ks > 0.) || !(norm2_ci > 0.))
  {
    fprintf(stderr,"Density integrals must be positive. CI=%g KS=%g\n",ne_ci,ne_ks);
    std::exit(EXIT_FAILURE);
  }

  // For KL and Jeffreys ONLY, discard grid points where BOTH electronic
  // densities are extremely small. Floating-point underflow in Gaussian
  // tails otherwise gives spurious infinite KL due to rho_CI>0,rho_KS=0
  // (or vice versa). Renormalize the retained shapes to unit integral.
  // Default 1e-10 e/bohr^3. Optional KL_RHO_CUTOFF input for sensitivity tests.
  double rho_cutoff = read_float("KL_RHO_CUTOFF");
  if (rho_cutoff == 0.) rho_cutoff = 1.e-10;
  if (!std::isfinite(rho_cutoff) || rho_cutoff <= 0.)
  {
    fprintf(stderr, "KL_RHO_CUTOFF must be positive and finite.\n");
    std::exit(EXIT_FAILURE);
  }
  double dropped_ci = 0., dropped_ks = 0.;
  int dropped_points = 0;
  for (int j=0; j<gsa; ++j)
  {
    const double w = wt[j];
    if (w == 0.) continue;
    const double a = std::max(0.,rho1[j]);
    const double b = std::max(0.,rho2[j]);
    if (std::max(a,b) < rho_cutoff)
    {
      dropped_ci += w*a;
      dropped_ks += w*b;
      ++dropped_points;
    }
  }
  const double ne_ci_kl = ne_ci-dropped_ci;
  const double ne_ks_kl = ne_ks-dropped_ks;
  if (!(ne_ci_kl > 0.) || !(ne_ks_kl > 0.))
  {
    fprintf(stderr, "KL tail cutoff eliminated all density.\n");
    std::exit(EXIT_FAILURE);
  }
  // Warn when tail removal is not negligible relative to the target density.
  if (dropped_ci/ne_ci > 1.e-8 || dropped_ks/ne_ks > 1.e-8)
    fprintf(stderr, "WARNING: KL cutoff removed non-negligible density; test cutoff convergence.\n");

  // KL is directional; Jeffreys is the sum of both KL directions.
  // JS distance and Hellinger distance are symmetric and bounded [0,1].
  double kl_ci_ks = 0., kl_ks_ci = 0.;
  double j_div = 0., js_div = 0., hellinger_sq = 0., tv = 0.;
  bool singular_ci_ks = false, singular_ks_ci = false;
  for (int j=0; j<gsa; ++j)
  {
    const double w = wt[j];
    if (w == 0.0) continue;
    const double p = std::max(0.0,rho1[j])/ne_ci;
    const double q = std::max(0.0,rho2[j])/ne_ks;
    const double m = 0.5*(p+q);
    if (m == 0.0) continue;

    tv += 0.5*w*std::abs(p-q);
    const double root_difference = std::sqrt(p)-std::sqrt(q);
    hellinger_sq += 0.5*w*root_difference*root_difference;
    // Jensen-Shannon is bounded and is evaluated on the FULL grid.
    if (p > 0.0 && q > 0.0)
    {
      const double t = (p-q)/(p+q);
      const double t2 = t*t;
      const double f = (std::abs(t) < 1.e-3)
        ? t2*(1.0 + t2/6.0 + t2*t2/15.0)
        : (1.0+t)*std::log1p(t) + (1.0-t)*std::log1p(-t);
      js_div += 0.5*w*m*f;
    }
    else
      js_div += w*m*std::log(2.0);

    // KL/Jeffreys: work on the retained, re-normalized shape functions.
    // Crucially, skipping the tail without re-normalizing could yield
    // negative directional KL values, despite KL being nonnegative.
    const double a = std::max(0.,rho1[j]);
    const double b = std::max(0.,rho2[j]);
    if (std::max(a,b) >= rho_cutoff)
    {
      const double pk = a/ne_ci_kl;
      const double qk = b/ne_ks_kl;
      if (pk > 0.0 && qk > 0.0)
      {
        const double lratio = (pk > 0.5*qk && pk < 2.0*qk)
          ? std::log1p((pk-qk)/qk) : (std::log(pk)-std::log(qk));
        kl_ci_ks += w*pk*lratio;
        kl_ks_ci -= w*qk*lratio;
        j_div += w*(pk-qk)*lratio;  // nonnegative pointwise
      }
      else if (pk > 0.0) singular_ci_ks = true;
      else if (qk > 0.0) singular_ks_ci = true;
    }
  }
  if (singular_ci_ks) kl_ci_ks = std::numeric_limits<double>::infinity();
  if (singular_ks_ci) kl_ks_ci = std::numeric_limits<double>::infinity();
  if (singular_ci_ks || singular_ks_ci)
    j_div = std::numeric_limits<double>::infinity();

  // The normalized KL integrals must be >=0, apart from tiny roundoff.
  if (std::isfinite(kl_ci_ks) && std::isfinite(kl_ks_ci) &&
      (kl_ci_ks < -1.e-12 || kl_ks_ci < -1.e-12 || j_div < -1.e-12))
  {
    fprintf(stderr, "Invalid negative KL/Jeffreys: KL1=%g KL2=%g J=%g\n",
            kl_ci_ks,kl_ks_ci,j_div);
    std::exit(EXIT_FAILURE);
  }
  if (std::isfinite(kl_ci_ks)) kl_ci_ks = std::max(0.,kl_ci_ks);
  if (std::isfinite(kl_ks_ci)) kl_ks_ci = std::max(0.,kl_ks_ci);

  const double js_dist = std::sqrt(std::max(0.0,js_div)/std::log(2.0));
  const double hellinger = std::sqrt(std::max(0.0,hellinger_sq));
  const double rel_l2 = diff2/std::sqrt(norm2_ci);
  // For matching electron counts, this is the Jeffreys value obtained
  // with electron densities rather than unit-normalized shape functions.
  const double j_electron_scaled = j_div*ne_ci;
  const double l2_per_e = diff2/ne_ci;
  const double l1_per_e = diff1/ne_ci;

  printf("\n XC100 KL tail cutoff=% .3e dropped CI=% .12e KS=% .12e (%d points)\n",
         rho_cutoff,dropped_ci,dropped_ks,dropped_points);
  printf(" XC100 density integrals: CI=% .12e KS=% .12e\n",ne_ci,ne_ks);
  printf(" XC100 L1/Ne(CI)=% .12e L2/Ne(CI)=% .12e relative_L2=% .12e\n",
         l1_per_e,l2_per_e,rel_l2);
  printf(" XC100 shape L1=% .12e TV=% .12e Hellinger=% .12e\n",
         2.0*tv,tv,hellinger);
  printf(" XC100 KL(CI||KS)=% .12e KL(KS||CI)=% .12e Jeffreys(normalized)=% .12e\n",
         kl_ci_ks,kl_ks_ci,j_div);
  printf(" XC100 Jeffreys(electron_scaled)=% .12e\n",j_electron_scaled);
  printf(" XC100 JS_divergence=% .12e JS_distance(normalized_0_to_1)=% .12e\n",
         js_div,js_dist);
  // Small output suitable for collecting across all XC100 species:
  printf("XC100_METRICS_CSV_HEADER,NeCI,NeKS,L1,L2,L1perNe,L2perNe,relativeL2,shapeL1,TV,KL_CI_KS,KL_KS_CI,JeffreysNormalized,JeffreysElectronScaled,JSdivergence,JSdistance,Hellinger\n");
  printf("XC100_METRICS_CSV,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
         ne_ci,ne_ks,diff1,diff2,l1_per_e,l2_per_e,rel_l2,2.0*tv,tv,
         kl_ci_ks,kl_ks_ci,j_div,j_electron_scaled,js_div,js_dist,hellinger);

  #pragma acc exit data delete(rho1[0:gsa],rho2[0:gsa],grid[0:gsa6],wt[0:gsa],gridf[0:gsa],coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  delete [] rho1;
  delete [] rho2;
  delete [] grid;
  delete [] gridf;
  delete [] wt;
  delete [] Pao1;
  delete [] Pao2;

  return;
}

void save_grid_rho(bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis)
{
  int gs = nrad*nang;
  int gsa = natoms*gs;
  int gsa6 = gsa*6;

  int N = basis.size();
  int N2 = N*N;

  double* Pao = new double[N2];
  bool found = read_square(N,Pao,"Pao_ci");
  if (!found) read_square(N,Pao,"Pao");
  if (!found) { printf("  couldn't find Pao or Pao_ci file \n"); return; }

  double* rho = new double[gsa];
  double* grid = new double[gsa6];
  float* gridf = new float[gsa6];
  double* wt = new double[gsa];

  printf("\n gbasis: %i natoms: %i nrad: %3i nang: %4i gsa: %4i \n",(int)gbasis,natoms,nrad,nang,gsa);

  #pragma acc enter data create(rho[0:gsa])
  #pragma acc enter data create(grid[0:gsa6],wt[0:gsa],gridf[0:gsa6])
  #pragma acc enter data copyin(coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  get_becke_grid_full(natoms,atno,coords,nrad,nang,ang_g,ang_w,6,grid,wt);
  if (gbasis)
  {
    //for (int j=0;j<gsa6;j++) gridf[j] = grid[j];
    //#pragma acc update device(gridf[0:gsa6])

    compute_rhodg(1,natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,NULL,1);
  }
  else
    compute_rhod(natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,NULL,NULL,1);

  #pragma acc update self(rho[0:gsa],grid[0:gsa6],wt[0:gsa])

  if (natoms==1)
  {
    printf("       r        rho \n");
    for (int j=0;j<nrad;j++)
      printf("  %8.5f  %8.5f \n",grid[6*j*nang+3],rho[j*nang]);
    printf("\n");
  }

  double dent = 0.;
 #pragma acc parallel loop present(rho[0:gsa],wt[0:gsa]) reduction(+:dent)
  for (int j=0;j<gsa;j++)
    dent += rho[j]*wt[j];
  printf("\n total e-s: %8.5f \n",dent);

  printf("\n  writing RHO_WF to disk \n");
  write_vector(gsa,rho,"RHO_WF");
  printf("\n  writing GRID_WTS to disk \n");
  write_grid(natoms,nrad,nang,grid,wt);

  #pragma acc exit data delete(rho[0:gsa],grid[0:gsa6],wt[0:gsa],gridf[0:gsa],coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  delete [] rho;
  delete [] grid;
  delete [] gridf;
  delete [] wt;
  delete [] Pao;

  return;
}

void save_grid_ao_basis(bool gbasis, int natoms, int* atno, double* coords, int nrad, int nang, double* ang_g, double* ang_w, vector<vector<double> >& basis)
{
  int gs = nrad*nang;
  int gsa = natoms*gs;
  int gsa6 = gsa*6;

  int N = basis.size();
  int N2 = N*N;

  double* Pao = new double[N2];
  double* rho = new double[gsa];
  double* grid = new double[gsa6];
  double* wt = new double[gsa];

  printf("\n gbasis: %i natoms: %i nrad: %3i nang: %4i gsa: %4i \n",(int)gbasis,natoms,nrad,nang,gsa);

  #pragma acc enter data create(rho[0:gsa])
  #pragma acc enter data create(grid[0:gsa6],wt[0:gsa])
  #pragma acc enter data copyin(coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  get_becke_grid_full(natoms,atno,coords,nrad,nang,ang_g,ang_w,6,grid,wt);

  for (int i=0;i<N;i++)
  for (int j=0;j<=i;j++)
  {
    for (int k=0;k<N2;k++) Pao[k] = 0.;
    Pao[i*N+j] = Pao[j*N+i] = 0.5;
    if (i==j) Pao[i*N+j] = 1.;

    if (gbasis)
      compute_rhodg(1,natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,NULL,1);
    else
      compute_rhod(natoms,atno,coords,basis,Pao,nrad,gsa,grid,rho,NULL,NULL,1);
    #pragma acc update self(rho[0:gsa])

    double dent = 0.;
   #pragma acc parallel loop present(rho[0:gsa],wt[0:gsa]) reduction(+:dent)
    for (int k=0;k<gsa;k++)
      dent += rho[k]*wt[k];
    printf("  %2i-%2i total e-s: %8.5f \n",i,j,dent);

    string filename = "RHO_AO_"+to_string(i+1)+"_"+to_string(j+1);
    write_vector(gsa,rho,filename);
  }

  #pragma acc exit data delete(rho[0:gsa],grid[0:gsa6],wt[0:gsa],coords[0:3*natoms],ang_g[0:3*nang],ang_w[0:nang])

  delete [] rho;
  delete [] grid;
  delete [] wt;
  delete [] Pao;

  return;
}

