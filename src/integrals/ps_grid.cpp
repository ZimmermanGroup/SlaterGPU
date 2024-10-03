#include "prosph.h"
#include "lebedev2.h"

//just edited  get_munuphi and munuphi_to_xyz


//enable zeta-dependent maximum volume (over mu)
#define VOLUME_RESIZE 0

#define RMAX 1.
//was 2.3 until recently
#define CF0 2.1




void shift_grid(const double shift_z, int gs, double* grid)
{
  printf(" shifting grid by %8.5f in z \n",shift_z);

  int gs6 = 6*gs;
 #pragma acc parallel loop present(grid[0:gs6])
  for (int j=0;j<gs;j++)
    grid[6*j+2] += shift_z;

  return;
}

//everything already on gpu
void reorient_grid(const double z0, int gs, double* grid, double* grid2, double* rot)
{
  int gs6 = 6*gs;
  //printf("  z0 in reorient grid: %8.5f \n",z0);

#if 0
 //already done in initialize_ps_coords/quad_grid_munu

 //put first focus on atom 1 (should be at 0,0,0)
  const double shift_z = -z0;
 #pragma acc parallel loop present(grid[0:gs6])
  for (int j=0;j<gs;j++)
    grid[6*j+2] += shift_z;
#endif

 #pragma acc parallel loop collapse(2) present(grid[0:gs6], rot[0:9], grid2[0:gs6])
  for (int i=0;i<gs;i++)
  {
    for(int j=0;j<3;j++)
    {
      double tsum = 0.;
     #pragma acc loop reduction(+:tsum)
      for(int k=0;k<3;k++)
        tsum += grid[6*i+k] * rot[3*k+j];

      grid2[6*i+j] = tsum;
    }
  }

  return;
}

void rotate_3x3(double* rot_mat, double* A, double* B)
{
  #pragma acc parallel loop collapse(2) present(rot_mat[0:9],A[0:9],B[0:9])
  for (int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      double tsum = 0.;
      #pragma acc loop reduction(+:tsum)
      for(int k=0;k<3;k++)
        tsum += A[3*i+k] * rot_mat[3*k+j];

      B[3*i+j] = tsum;
    }
  }
}

void rotate_3x3t_cpu(double* rot_mat, double* A, double* B)
{
  for (int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      double tsum = 0.;
      for(int k=0;k<3;k++)
        tsum += A[3*i+k] * rot_mat[3*j+k];

      B[3*i+j] = tsum;
    }
  }
}

void rotate_4x3t_cpu(double* rot_mat, double* A, double* B)
{
  for (int i=0;i<4;i++)
  {
    for(int j=0;j<3;j++)
    {
      double tsum = 0.;
      for(int k=0;k<3;k++)
        tsum += A[3*i+k] * rot_mat[3*j+k];

      B[3*i+j] = tsum;
    }
  }
}

void rotate_3xyz(double& A1, double& B1, double& C1, double* rot)
{
  double A0 = A1; double B0 = B1; double C0 = C1;
 /*
  A1 = A0*rot[0] + B0*rot[3] + C0*rot[6];
  B1 = A0*rot[1] + B0*rot[4] + C0*rot[7];
  C1 = A0*rot[2] + B0*rot[5] + C0*rot[8];
 */
  A1 = A0*rot[0] + B0*rot[1] + C0*rot[2];
  B1 = A0*rot[3] + B0*rot[4] + C0*rot[5];
  C1 = A0*rot[6] + B0*rot[7] + C0*rot[8];
}

void set_up_axes(int natoms, double* coords, double* coord1, double* first_rot)
{
  if (coords[4]==0. && coords[5]==0.)
  {
   #pragma acc serial present(first_rot[0:9])
    {
      bool fs = 1; if (coords[3]<0.) fs = 0;
      first_rot[4] = 1.;
      if (!fs) { first_rot[6] = -1.; first_rot[2] = 1.; }
      else { first_rot[6] = 1.; first_rot[2] = -1.; }
    }
  }
  else if (coords[5]==0.)
  {
   #pragma acc serial present(first_rot[0:9])
    {
      bool fs = 1; if (coords[4]<0.) fs = 0;
      first_rot[0] = 1.;
      if (!fs) { first_rot[7] = -1.; first_rot[5] = 1.; }
      else { first_rot[7] = 1.; first_rot[5] = -1.; }   
    }
  }
  else
  {
   #pragma acc serial present(first_rot[0:9])
    {
      first_rot[0] = first_rot[4] = first_rot[8] = 1.;
      if (coords[5]>0.) first_rot[8] = -1.;
    }
  }
  rotate_3x3(first_rot,coords,coord1);

 #if 0
  #pragma acc update self(coord1[0:9])
  printf(" after setup. xyz: \n");
  for (int j=0;j<natoms;j++)
    printf(" %8.5f %8.5f %8.5f \n",coord1[3*j+0],coord1[3*j+1],coord1[3*j+2]);
 #endif

  return;
}

//coords/trot already on gpu
void gen_total_rot_n2(int natoms, double* coords, double* trot)
{
  if (natoms==1)
  {
    #pragma acc serial
      trot[0] = trot[4] = trot[8] = 1.;
    return;
  }

 //first atom must be at (0,0,0)
  if (coords[0]!=0. || coords[1]!=0. || coords[2]!=0.)
  {
    printf("\n ERROR: gen_total_rot_n2. first atom must be at (0,0,0) \n");
    exit(-1);
  }

  double coord1[9];
  double first_rot[9];
  double rot_mat_1[9];
  double rot_mat_2[9];
  double tmp[9];

  #pragma acc enter data create(first_rot[0:9],rot_mat_1[0:9],rot_mat_2[0:9],tmp[0:9])
  #pragma acc enter data create(coord1[0:9])

 #pragma acc parallel loop present(tmp[0:9])
  for(int i=0;i<9;i++)
    tmp[i] = 0.;

  #pragma acc parallel loop present(rot_mat_1[0:9],rot_mat_2[0:9],trot[0:9])
  for(int i=0;i<9;i++)
  {
    first_rot[i] = 0.;
    rot_mat_1[i] = 0.;
    rot_mat_2[i] = 0.;
    trot[i] = 0.;
  }

  set_up_axes(natoms,coords,coord1,first_rot);

  #pragma acc serial present(rot_mat_1[0:9],coord1[0:9])
  {
    double angle_1 = 0.;
    if (coord1[5]!=0.)
      angle_1 = atan(-coord1[3*1+0]/coord1[3*1+2]);
    //printf("  angle_1: %8.5f \n",angle_1);
   //xz rotation
    rot_mat_1[0] = cos(angle_1);
    rot_mat_1[2] = -sin(angle_1);
    rot_mat_1[4] = 1.;
    rot_mat_1[6] = sin(angle_1);
    rot_mat_1[8] = cos(angle_1);
  }

  rotate_3x3(rot_mat_1,coord1,tmp);

 #if 0
  #pragma acc update self(tmp[0:9])
  printf("   after rot1 \n");
  for (int j=0;j<natoms;j++)
    printf(" %8.5f %8.5f %8.5f \n",tmp[3*j+0],tmp[3*j+1],tmp[3*j+2]);
 #endif

  #pragma acc serial present(rot_mat_2[0:9],tmp[0:9])
  {
    double angle_2 = 0.;
    if (tmp[5]!=0.)
      angle_2 = atan(-tmp[3*1+1]/tmp[3*1+2]);
    //printf("  angle_2: %8.5f \n",angle_2);
    rot_mat_2[0] = 1.;
    rot_mat_2[4] = cos(angle_2);
    rot_mat_2[5] = -sin(angle_2);
    rot_mat_2[7] = sin(angle_2);
    rot_mat_2[8] = cos(angle_2);
  }


 //construct full rotation matrix
  rotate_3x3(rot_mat_1,first_rot,tmp);

  #pragma acc parallel loop present(rot_mat_1[0:9],tmp[0:9])
  for (int j=0;j<9;j++)
    rot_mat_1[j] = tmp[j];

  rotate_3x3(rot_mat_2,rot_mat_1,trot);

 #if 0
  #pragma acc update self(trot[0:9])
  printf("\n trot: \n");
  for (int j=0;j<3;j++)
    printf("   %8.5f %8.5f %8.5f \n",trot[3*j+0],trot[3*j+1],trot[3*j+2]);
 #endif

 //check the result
 #if 0
  rotate_3x3(trot,coords,tmp);
  #pragma acc update self(tmp[0:9])
  printf("\n final orientation: \n");
  for (int j=0;j<natoms;j++)
    printf(" %8.5f %8.5f %8.5f \n",tmp[3*j+0],tmp[3*j+1],tmp[3*j+2]);
  if (fabs(tmp[4])>1.e-10) printf(" WARNING: y axis off \n"); 
 #endif

  #pragma acc update self(trot[0:9])
  #pragma acc exit data delete(coord1[0:9])
  #pragma acc exit data delete(first_rot[0:9],rot_mat_1[0:9],rot_mat_2[0:9],tmp[0:9])

  return;    
}

void gen_rot_n3z1(double* coords, double* rot_mat_1, double* rot_mat_2, double* rot_mat_3, double* tmp1, double* tmp2)
{
  #pragma acc serial present(rot_mat_1[0:9],coords[0:9])
  {
    double angle_1 = atan(-coords[3*1+0]/coords[3*1+2]);
   //rotate around y axis
    rot_mat_1[0] = cos(angle_1);
    rot_mat_1[2] = -sin(angle_1);
    rot_mat_1[4] = 1.;
    rot_mat_1[6] = sin(angle_1);
    rot_mat_1[8] = cos(angle_1);
  }

  rotate_3x3(rot_mat_1,coords,tmp1);

  #pragma acc serial present(rot_mat_2[0:9],tmp1[0:9])
  {
    double angle_2 = 0.;
    if (tmp1[5]!=0.)
      angle_2 = atan(-(tmp1[3*1+1]/tmp1[3*1+2]));
   //rotate around x axis
    rot_mat_2[0] = 1.;
    rot_mat_2[4] = cos(angle_2);
    rot_mat_2[5] = -sin(angle_2);
    rot_mat_2[7] = sin(angle_2);
    rot_mat_2[8] = cos(angle_2);
  }

  rotate_3x3(rot_mat_2,tmp1,tmp2);

  #pragma acc serial present(tmp2[0:9],rot_mat_3[0:9])
  {
    double angle_3 = 0.;
    if (tmp2[6]!=0.)
      angle_3 = atan(tmp2[3*2+1]/tmp2[3*2+0]);
    else
      angle_3 = PI*0.5;
    //printf(" angle_3: %8.5f \n",angle_3);
   //rotate around z axis
    rot_mat_3[0] = cos(angle_3);
    rot_mat_3[1] = -sin(angle_3);
    rot_mat_3[3] = sin(angle_3);
    rot_mat_3[4] = cos(angle_3);
    rot_mat_3[8] = 1.;
  }

  rotate_3x3(rot_mat_3,tmp2,tmp1);

  #pragma acc serial present(tmp1[0:9],rot_mat_3[0:9])
  {
   //x positive
    if (tmp1[6]<0.)
    {
      rot_mat_3[0] *= -1.;
      rot_mat_3[3] *= -1.;
    }
  }

  return;
}


//coords/trot already on gpu
void gen_total_rot_n3(int natoms, double* coords, double* trot)
{
  //int N3 = 3*natoms;

  if (natoms==2)
    return gen_total_rot_n2(natoms,coords,trot);

  if (natoms!=3) 
  {
    printf("\n ERROR: gen_total_rot_n3. natoms must be equal to 3 \n");
    exit(-1);
  }
 //first atom must be at (0,0,0)
  if (coords[0]!=0. || coords[1]!=0. || coords[2]!=0.)
  {
    printf("\n ERROR: gen_total_rot_n3. first atom must be at (0,0,0) \n");
    exit(-1);
  }

  double coord1[9];
  double first_rot[9];
  double part_rot[9];
  double rot_mat_1[9];
  double rot_mat_2[9];
  double rot_mat_3[9];
  double tmp1[9];
  double tmp2[9];
  #pragma acc enter data create(first_rot[0:9],part_rot[0:9],rot_mat_1[0:9],rot_mat_2[0:9],rot_mat_3[0:9],tmp1[0:9],tmp2[0:9])
  #pragma acc enter data create(coord1[0:9])

  #pragma acc parallel loop present(tmp1[0:9],tmp2[0:9],first_rot[0:9],part_rot[0:9],rot_mat_1[0:9],rot_mat_2[0:9],rot_mat_3[0:9],trot[0:9])
  for(int i=0;i<9;i++)
  {
    tmp1[i] = 0.;
    tmp2[i] = 0.;

    first_rot[i] = 0.;
    part_rot[i] = 0.;
    rot_mat_1[i] = 0.;
    rot_mat_2[i] = 0.;
    rot_mat_3[i] = 0.;
    trot[i] = 0.;
  }

  set_up_axes(natoms,coords,coord1,first_rot);
  //#pragma acc enter data copyin(coord1[0:9])

  gen_rot_n3z1(coord1,rot_mat_1,rot_mat_2,rot_mat_3,tmp1,tmp2);

 //construct full rotation matrix
  rotate_3x3(rot_mat_1,first_rot,tmp1);
  rotate_3x3(rot_mat_2,tmp1,part_rot);
  rotate_3x3(rot_mat_3,part_rot,trot);

 #if 0
 //check the result
  rotate_3x3(trot,coords,tmp1);

  #pragma acc update self(tmp1[0:9])
  printf("\n final orientation: \n");
  for (int j=0;j<3;j++)
    printf(" %8.5f %8.5f %8.5f \n",tmp1[3*j+0],tmp1[3*j+1],tmp1[3*j+2]);
  if (fabs(tmp1[4])>1.e-10) printf(" WARNING: y2 axis off \n"); 
  if (fabs(tmp1[7])>1.e-10) printf(" WARNING: y3 axis off \n"); 
  if (tmp1[6]<0.) printf(" WARNING: x3 axis negative \n");
 #endif

  #pragma acc exit data delete(first_rot[0:9],part_rot[0:9],rot_mat_1[0:9],rot_mat_2[0:9],rot_mat_3[0:9],tmp1[0:9],tmp2[0:9])
  #pragma acc exit data delete(coord1[0:9])

  #pragma acc update self(trot[0:9])

  return;    
}

int sign(double x)
  {
    if (x>0) return 1;
    else if (x<=0) return -1;
    return 0;
  }

void munuphi_to_xyz(double a, double mu, double nu, double phi, double& x, double& y, double& z)
{
  double sinhm = sinh(mu);
  double coshm = cosh(mu);
  double sinn = sin(nu);
  double cosn = cos(nu);
  x = a*sinhm*sinn*cos(phi);
  y = a*sinhm*sinn*sin(phi);
  z = a*coshm*cosn - a;
}

#pragma acc routine seq
void munuphi_to_xyz_gpu(double a, double mu, double nu, double phi, double& x, double& y, double& z)
{
  double sinhm = sinh(mu);
  double coshm = cosh(mu);
  double sinn = sin(nu);
  double cosn = cos(nu);
  x = a*sinhm*sinn*cos(phi);
  y = a*sinhm*sinn*sin(phi);
  z = a*coshm*cosn - a;
}

void munuphi_to_xyz(double a, int gs, double* gridm, double* grid, double* wt)
{
  double a3 = a*a*a;
  if (grid!=NULL && wt!=NULL)
  for (int j=0;j<gs;j++)
  {
    double mu  = gridm[6*j+0];
    double nu  = gridm[6*j+1];
    double phi = gridm[6*j+2];
    double dmu = gridm[6*j+3];
    double dnu = gridm[6*j+4];
    double dphi = gridm[6*j+5];

    double sinhm = sinh(mu);
    double coshm = cosh(mu);
    double sinn = sin(nu);
    double cosn = cos(nu);

    double x = a*sinhm*sinn*cos(phi);
    double y = a*sinhm*sinn*sin(phi);
    double z = a*coshm*cosn;

    double wt1 = a3*ps_dV(mu-0.5*dmu,mu+0.5*dmu,nu-0.5*dnu,nu+0.5*dnu) * dphi;

    //if (x>xmax) xmax = x;
    //if (y>ymax) ymax = y;
    //if (z>zmax) zmax = z;

    grid[6*j+0] = x;
    grid[6*j+1] = y;
    grid[6*j+2] = z-a;
    wt[j] = wt1;
  }

  //if (prl>1) printf("  xyz max: %9.5f %9.5f %9.5f \n",xmax,ymax,zmax);
}

//determines range of mu values
double get_c0(double a, double nmu, double cf, double rm)
{
  double ca = CF0*cf*pow(a,-0.25);
  // rm *= cf;
  //double cb = 2./log(3.) * acosh(1.+rm/a);

  return ca;
}

void initialize_ps_coords_1c(double a, double cf, const double c2, int nmu, int nnu, int nphi, double phi_zero, double* grid, double* gridm, double* wt, int prl)
{
  if (a<=0.) { printf(" ERROR: initialize_ps_coords with a<0 \n",a); exit(-1); }
  //if (nnu<2 || nphi<2) { printf(" ERROR: cannot initialize PS coordinates with nnu<2 or nphi<2 \n"); exit(-1); }

  double c0 = get_c0(a,nmu,cf,RMAX);
  if (cf<0.) c0 = -cf; //c1 passed to ftn
  const double c1 = c0;

  int gs = nmu*nnu*abs(nphi);
  int gs6 = 6*gs;

  const double dnx = PI/nnu;
  double dphi0 = 2.*PI/nphi;
  if (nphi==-1) dphi0 = 1.;
  nphi = abs(nphi);
  const double dphi = dphi0;
  const double dx1 = 1./(nmu+1);
  const double dx2 = 1/nnu;
  double a3 = a*a*a;

  double rmax = a*cosh(c1*atanh(nmu*dx1));
  printf("    cf: %5.1f  c1: %8.5f  rmax: %8.3f  c2: %8.5f \n",cf,c1,rmax,c2);

  //double xmax = 0.; double ymax = 0.; double zmax = 0.;
  for (int i=0;i<nmu;i++)
  {
    double x0  =  i*dx1;
    double x1  = x0+dx1;
    double mu0 = c1*atanh(x0);
    double mu1 = c1*atanh(x1);
    double mu = 0.5*(mu1+mu0); //"average" mu value

    double sinhm = sinh(mu);
    double coshm = cosh(mu);
    double dmu = mu1-mu0;

    for (int j=0;j<nnu;j++)
    {
      double x20 = j*dnx;
      double x21 = (j+1)*dnx;

     //removed 2.*x2
      double nu0 = x20 + c2*sin(x20);
      double nu1 = x21 + c2*sin(x21);

      double nu = 0.5*(nu0+nu1);
      double dnu = nu1 - nu0;

      if (nnu<12) printf("  x2: %8.5f %8.5f  nu: %8.5f %8.5f  nu(m): %8.5f  dnu: %8.5f \n",x20,x21,nu0,nu1,nu,dnu);

      double sinn = sin(nu);
      double cosn = cos(nu);

      int i1 = i*nnu*nphi + j*nphi;
     #pragma acc parallel loop present(grid[0:gs6],gridm[0:gs6],wt[0:gs])
      for (int k=0;k<nphi;k++)
      {
        int i2 = i1+k;
        double phi = k*dphi + phi_zero;

        gridm[6*i2+0] = mu; gridm[6*i2+1] = nu; gridm[6*i2+2] = phi;
        gridm[6*i2+3] = dmu; gridm[6*i2+4] = dnu; gridm[6*i2+5] = dphi;

        //if (prl>2) printf("    mu/nu/phi: %8.5f %8.5f %8.5f  wt: %8.5f \n",mu,nu,phi,wt1);

      } //loop phi
    } //loop nu
  } //loop mu

  munuphi_to_xyz(a,gs,gridm,grid,wt);

 #if 0
  #pragma acc update self(gridm[0:gs6])
  for (int i=0;i<nmu;i++)
  {
    double mu1 = gridm[6*i*nnu*nphi+0];
    double r1 = a*cosh(mu1);
    printf("  mu: %8.5f  r: %5.2f \n",mu1,r1);
  }
 #endif

  return;
}

void initialize_ps_coords_2c(double a, double cf, int nmu, int nnu, int nphi, double phi_zero, double* grid, double* gridm, double* wt, int prl)
{
  if (a<=0.) { printf(" ERROR: initialize_ps_coords with a<0 \n",a); exit(-1); }
  //if (nnu<2 || nphi<2) { printf(" ERROR: cannot initialize PS coordinates with nnu<2 or nphi<2 \n"); exit(-1); }

  double c0 = get_c0(a,nmu,cf,RMAX);
  if (cf<0.) c0 = -cf; //c1 passed to ftn
  const double c1 = c0;

  int gs = nmu*nnu*abs(nphi);
  int gs6 = 6*gs;

 //c1==2   -> about 20 Bohr
 //c1==2.3 -> about 30 Bohr
 //c1==3   -> about 100 Bohr
  const double dnu = PI/nnu;
  double dphi0 = 2.*PI/nphi;
  if (nphi==-1) dphi0 = 1.;
  nphi = abs(nphi);
  const double dphi = dphi0;
  const double dx1 = 1./(nmu+1);
  double a3 = a*a*a;

  if (prl>0)
    printf("\n   init_ps_coords. a: %8.5f nmu/nu/phi: %2i %2i %2i  dnu/dphi: %8.5f %8.5f  dx1: %8.5f \n",a,nmu,nnu,nphi,dnu,dphi,dx1);
  double rmax = a*cosh(c1*atanh(nmu*dx1));
  printf("    cf: %5.1f  c1: %8.5f  rmax: %8.3f \n",cf,c1,rmax);

  //double xmax = 0.; double ymax = 0.; double zmax = 0.;
  for (int i=0;i<nmu;i++)
  {
    double x0  =  i*dx1;
    double x1  = x0+dx1;
    double mu0 = c1*atanh(x0);
    double mu1 = c1*atanh(x1);
    double mu = 0.5*(mu1+mu0); //"average" mu value

    double sinhm = sinh(mu);
    double coshm = cosh(mu);
    double dmu = mu1-mu0;

    //if (i+1==nmu) printf("  x1/0: %8.5f %8.5f  mu/0: %8.5f %8.5f  dmu: %8.5f \n",x1,x0,mu,mu0,dmu);
    //if (i+1==nmu) printf("    c1: %8.5f  mu/0: %8.5f %8.5f %8.5f \n",c1,mu,mu0,mu1);

    for (int j=0;j<nnu;j++)
    {
      double nu = (j+0.5)*dnu;

      double sinn = sin(nu);
      double cosn = cos(nu);

      int i1 = i*nnu*nphi + j*nphi;
     #pragma acc parallel loop present(grid[0:gs6],gridm[0:gs6],wt[0:gs])
      for (int k=0;k<nphi;k++)
      {
        int i2 = i1+k;
        double phi = k*dphi + phi_zero;

        gridm[6*i2+0] = mu; gridm[6*i2+1] = nu; gridm[6*i2+2] = phi;
        gridm[6*i2+3] = dmu; gridm[6*i2+4] = dnu; gridm[6*i2+5] = dphi;

        //if (prl>2) printf("    mu/nu/phi: %8.5f %8.5f %8.5f  wt: %8.5f \n",mu,nu,phi,wt1);

      } //loop phi
    } //loop nu
  } //loop mu

  munuphi_to_xyz(a,gs,gridm,grid,wt);

 #if 0
  #pragma acc update self(gridm[0:gs6])
  for (int i=0;i<nmu;i++)
  {
    double mu1 = gridm[6*i*nnu*nphi+0];
    double r1 = a*cosh(mu1);
    printf("  mu: %8.5f  r: %5.2f \n",mu1,r1);
  }
 #endif

 #if 0
  if (grid!=NULL && wt!=NULL)
  for (int j=0;j<gs;j++)
  {
    double mu  = gridm[6*j+0];
    double nu  = gridm[6*j+1];
    double phi = gridm[6*j+2];
    double dmu = gridm[6*j+3];
    double dnu = gridm[6*j+4];

    double sinhm = sinh(mu);
    double coshm = cosh(mu);
    double sinn = sin(nu);
    double cosn = cos(nu);

    double x = a*sinhm*sinn*cos(phi);
    double y = a*sinhm*sinn*sin(phi);
    double z = a*coshm*cosn;

    double wt1 = a3*ps_dV(mu-0.5*dmu,mu+0.5*dmu,nu-0.5*dnu,nu+0.5*dnu) * dphi;

    //if (x>xmax) xmax = x;
    //if (y>ymax) ymax = y;
    //if (z>zmax) zmax = z;

    grid[6*j+0] = x;
    grid[6*j+1] = y;
    grid[6*j+2] = z-a;
    wt[j] = wt1;
  }

  //if (prl>1) printf("  xyz max: %9.5f %9.5f %9.5f \n",xmax,ymax,zmax);
 #endif

  return;
}

void initialize_ps_coords_2c_phi(double a, double cf, int nmu, int nnu, int nphi, double phi_zero, double phi_add, double* grid, double* gridm, double* wt, int prl)
{
  if (a<=0.) { printf(" ERROR: initialize_ps_coords with a<0 \n",a); exit(-1); }
  if (nphi<3) { printf(" ERROR: initialize_ps_coords_phi with nphi<3 \n"); exit(-1); }

  printf("  testing initialize_ps_coords_2c_phi \n");

  double c0 = get_c0(a,nmu,cf,RMAX);
  const double c1 = c0;

  int gs = nmu*nnu*abs(nphi);
  int gs6 = 6*gs;

 //c1==2   -> about 20 Bohr
 //c1==2.3 -> about 30 Bohr
 //c1==3   -> about 100 Bohr
  const double dnu = PI/nnu;
  double dphi0 = 2.*PI/nphi;
  if (nphi==-1) dphi0 = 1.;
  nphi = abs(nphi);
  const double dphi = dphi0;
  const double dx1 = 1./(nmu+1);
  double a3 = a*a*a;

  double phia[nphi];
  double phiad[nphi];
  for (int k=0;k<nphi;k++)
  {
    double phi = k*dphi + phi_zero;
    phia[k] = phi;
  }
  int wk0 = -1;
  for (int k=0;k<nphi;k++)
  {
    double phi1 = phia[k];
    double phi2 = phi1+dphi;
    if (phi_add>phi1 && phi_add<phi2)
    { wk0 = k; break; }
  }
  int wk = wk0; if (wk0<2) wk = 2; if (wk0==nphi-1) wk = nphi-2;
  //printf("   wk0: %2i wk: %2i \n",wk0,wk);

  double dphik1 = phi_add/wk;
  for (int k=0;k<wk;k++)
    phiad[k] = dphik1;
  double dphik2 = (2.*PI-phi_add)/(nphi-wk);
  for (int k=wk;k<nphi;k++)
    phiad[k] = dphik2;

  for (int k=1;k<nphi;k++)
  {
    double phi = phia[k-1] + phiad[k-1];
    phia[k] = phi;
  }

  if (prl>0)
  {
    printf("  phia: ");
    for (int k=0;k<nphi;k++)
      printf(" %8.5f",phia[k]);
    printf("\n");

    printf("  phiad: ");
    for (int k=0;k<nphi;k++)
      printf(" %8.5f",phiad[k]);
    printf("\n");
  }

  #pragma acc enter data copyin(phia[0:nphi],phiad[0:nphi])

  if (prl>0)
  printf("\n init_ps_coords. a: %8.5f nmu/nu/phi: %2i %2i %2i  dnu/dphi: %8.5f %8.5f  dx1: %8.5f \n",a,nmu,nnu,nphi,dnu,dphi,dx1);

  //double xmax = 0.; double ymax = 0.; double zmax = 0.;
  for (int i=0;i<nmu;i++)
  {
    double x0  =  i*dx1;
    double x1  = x0+dx1;
    double mu0 = c1*atanh(x0);
    double mu1 = c1*atanh(x1);
    double mu = 0.5*(mu1+mu0); //"average" mu value

    double sinhm = sinh(mu);
    double coshm = cosh(mu);
    double dmu = mu1-mu0;

    //printf("  x1/0: %8.5f %8.5f  mu/0: %8.5f %8.5f  dmu: %8.5f \n",x1,x0,mu,mu0,dmu);
    //if (i+1==nmu) printf("  mu/0: %8.5f %8.5f \n",mu,mu0);

    for (int j=0;j<nnu;j++)
    {
      double nu = (j+0.5)*dnu;

      double sinn = sin(nu);
      double cosn = cos(nu);

      int i1 = i*nnu*nphi + j*nphi;
     #pragma acc parallel loop present(grid[0:gs6],gridm[0:gs6],wt[0:gs],phia[0:nphi],phiad[0:nphi])
      for (int k=0;k<nphi;k++)
      {
        int i2 = i1+k;
        double phi = phia[k];
        double dphi = phiad[k];

        gridm[6*i2+0] = mu; gridm[6*i2+1] = nu; gridm[6*i2+2] = phi;
        gridm[6*i2+3] = dmu; gridm[6*i2+4] = dnu; gridm[6*i2+5] = dphi;

        //if (prl>2) printf("    mu/nu/phi: %8.5f %8.5f %8.5f  wt: %8.5f \n",mu,nu,phi,wt1);

      } //loop phi
    } //loop nu
  } //loop mu

  #pragma acc exit data delete(phia[0:nphi],phiad[0:nphi])

  munuphi_to_xyz(a,gs,gridm,grid,wt);

  return;
}

void get_munuphi(double a, double A3, double B3, double C3, double& mu, double& nu, double& phi)
{
 //just changed this, due to first atom being at zero
  C3 += a;
  double ab2 = A3*A3+B3*B3;
  double rr = a*a + ab2 + C3*C3;
  double rrm = a*a + ab2 - C3*C3;
  double r1 = sqrt(ab2+(C3-a)*(C3-a));
  double r2 = sqrt(ab2+(C3+a)*(C3+a));
  double ra = rr+r1*r2;
 //avoid near-miss of zero
  double rb = fabs(rrm+r1*r2);
  double sq2 = sqrt(2.);

  double mar = sqrt(ra)/(sq2*a);
  mu = acosh(mar); if (fabs(mar-1.)<1.e-10) mu = 0.;
  nu = 0.5*PI - atan(sq2*C3*sqrt(1./rb));
  phi = 0.; if (A3!=0.0) { if (B3==0. && A3<0.) phi = -PI; else phi = atan(B3/A3); } else phi = PI*0.5*sign(B3);
  if (phi<-1.e-10) phi += 2.*PI;

  if (std::isnan(mu) || std::isnan(nu))
  {
    printf("  ra: %8.5e sqrt(ra): %8.5e  sqra/2sqa: %8.5e \n",ra,sqrt(ra),sqrt(ra)/(sq2*a));
    printf("  ab2: %8.5e rr: %8.5e r12: %8.5e %8.5e ra: %8.5e rb: %8.5e C3: %8.5e \n",ab2,rr,r1,r2,ra,rb,C3);
  }
  return;
}

void shift_boundary_munu(double mu, double nu, int ix, int jx, int ib, int ia, int jb, int ja, int kb, int ka, int gs, double* gridm)
{
  if (ix<0 || jx<0 || ib<0 || ia<0 || jb<0 || ja<0 || kb<0 || ka<0)
  { printf(" ERROR: shift_boundary_munu \n"); exit(-1); }

  int gs6 = 6*gs;

  double mu_ib  = gridm[ix*ib+0];
  double mu_ia  = gridm[ix*ia+0];
  double dmu_ib = gridm[ix*ib+3];
  double dmu_ia = gridm[ix*ia+3];

  double nu_jb  = gridm[jx*jb+1];
  double nu_ja  = gridm[jx*ja+1];
  double dnu_jb = gridm[jx*jb+4];
  double dnu_ja = gridm[jx*ja+4];

  //double phi_kb  = gridm[6*kb+2];
  //double phi_ka  = gridm[6*ka+2];
  //double dphi    = gridm[6*kb+5];

 //new mu division
  double mu0 = mu_ib-dmu_ib*0.5;
  double mu1 = mu;
  double mu2 = mu_ia+dmu_ia*0.5;

 //new nu division
  double nu0 = nu_jb-dnu_jb*0.5;
  double nu1 = nu;
  double nu2 = nu_ja+dnu_ja*0.5;

 //recenters mu/nu so 3rd atom is at boundary of cells, leave phi alone
  #pragma acc serial present(gridm[0:gs6])
  for (int k=0;k<2;k++)
  {
    int kn = kb; if (k==1) kn = ka;
    //double phin = phi_kb; if (k==1) phin = phi_ka;
    int in = ib*ix+jb*jx+6*kn;

   //first cell
    double dmu = mu1-mu0;
    double dnu = nu1-nu0;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
    gridm[in+0] = mu0+0.5*dmu;
    gridm[in+1] = nu0+0.5*dnu;
    gridm[in+3] = dmu;
    gridm[in+4] = dnu;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);

   //second cell
    in = ib*ix+ja*jx+6*kn;
    dnu = nu2-nu1;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
    gridm[in+0] = mu0+0.5*dmu;
    gridm[in+1] = nu1+0.5*dnu;
    gridm[in+3] = dmu;
    gridm[in+4] = dnu;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);

   //third cell
    in = ia*ix+jb*jx+6*kn;
    dmu = mu2-mu1;
    dnu = nu1-nu0;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
    gridm[in+0] = mu1+0.5*dmu;
    gridm[in+1] = nu0+0.5*dnu;
    gridm[in+3] = dmu;
    gridm[in+4] = dnu;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);

   //fourth cell
    in = ia*ix+ja*jx+6*kn;
    dmu = mu2-mu1;
    dnu = nu2-nu1;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
    gridm[in+0] = mu1+0.5*dmu;
    gridm[in+1] = nu1+0.5*dnu;
    gridm[in+3] = dmu;
    gridm[in+4] = dnu;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
  }

  return;
}

void shift_boundary_nu(double nu, int ix, int jx, int jb, int ja, int kb, int ka, int gs, double* gridm)
{
  if (ix<0 || jx<0 || jb<0 || ja<0 || kb<0 || ka<0)
  { printf(" ERROR: shift_boundary_nu \n"); exit(-1); }

 //mu is close to zero
  int ib = 0;

  int gs6 = 6*gs;

  double nu_jb  = gridm[jx*jb+1];
  double nu_ja  = gridm[jx*ja+1];
  double dnu_jb = gridm[jx*jb+4];
  double dnu_ja = gridm[jx*ja+4];

 //new nu division
  double nu0 = nu_jb-dnu_jb*0.5;
  double nu1 = nu;
  double nu2 = nu_ja+dnu_ja*0.5;

 //recenters nu so 3rd atom is at boundary of cells, leave mu/phi alone
  #pragma acc serial present(gridm[0:gs6])
  for (int k=0;k<2;k++)
  {
    int kn = kb; if (k==1) kn = ka;
    int in = ib*ix+jb*jx+6*kn;

   //first cell
    double dnu = nu1-nu0;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
    gridm[in+1] = nu0+0.5*dnu;
    gridm[in+4] = dnu;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f (dnu) \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);

   //second cell
    in = ib*ix+ja*jx+6*kn;
    dnu = nu2-nu1;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
    gridm[in+1] = nu1+0.5*dnu;
    gridm[in+4] = dnu;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f (dnu) \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
  }

  return;
}

void shift_boundary_mu(double mu, int ix, int jx, int ib, int ia, int jb, int kb, int ka, int gs, double* gridm)
{
  if (ix<0 || jx<0 || ib<0 || ia<0 || jb<0 || kb<0 || ka<0)
  { printf(" ERROR: shift_boundary_mu \n"); exit(-1); }

  int gs6 = 6*gs;

  double mu_ib  = gridm[ix*ib+0];
  double mu_ia  = gridm[ix*ia+0];
  double dmu_ib = gridm[ix*ib+3];
  double dmu_ia = gridm[ix*ia+3];

 //new mu division
  double mu0 = mu_ib-dmu_ib*0.5;
  double mu1 = mu;
  double mu2 = mu_ia+dmu_ia*0.5;

 //recenters mu so 3rd atom is at boundary of cells, leave nu/phi alone
  #pragma acc serial present(gridm[0:gs6])
  for (int k=0;k<2;k++)
  {
    int kn = kb; if (k==1) kn = ka;
    int in = ib*ix+jb*jx+6*kn;

   //first cell
    double dmu = mu1-mu0;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
    gridm[in+0] = mu0+0.5*dmu;
    gridm[in+3] = dmu;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f (dmu) \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);

   //second cell
    in = ia*ix+jb*jx+6*kn;
    dmu = mu2-mu1;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
    gridm[in+0] = mu1+0.5*dmu;
    gridm[in+3] = dmu;
    //printf("   in: %3i  munuphi: %8.5f %8.5f %8.5f (dmu) \n",in,gridm[in+0],gridm[in+1],gridm[in+2]);
  }

  return;
}


bool handle_linear_nu(int& jb, int& ja, int nnu, double nu, double munu_thresh)
{
  bool skip_nu = 0;
  if (jb==ja)
  {
    if (jb==0 && nu>munu_thresh)
    {
      printf("   near-linear angle found in 3-atom setup (nu: %8.5f) \n",nu);
      ja = 1;
    }
    else
      skip_nu = 1;
  }
  if (jb==nnu-1)
  {
    if (nu<PI-munu_thresh)
    {
      printf("   near-linear angle found in 3-atom setup (nu3: %8.5f) \n",nu);
      jb = nnu-2;
      ja = nnu-1;
    }
    else
    {
      jb = ja = nnu-1;
      skip_nu = 1;
    }
  }

  return skip_nu;
}

void initialize_ps_coords_3c(double cf, int nmu, int nnu, int nphi, double phi_zero, double* coordn, double* grid, double* gridm, double* wt, double* rot, int prl)
{
 //this code assumes first two atoms aligned along z axis, after the rotation
  double coordm[9];
  rotate_3x3t_cpu(rot,coordn,coordm);

  double A2 = coordm[3]; double B2 = coordm[4]; double C2 = coordm[5];
  double A3 = coordm[6]; double B3 = coordm[7]; double C3 = coordm[8];
  double a = sqrt(A2*A2+B2*B2+C2*C2)*0.5;

  initialize_ps_coords_2c(a,cf,nmu,nnu,nphi,phi_zero,grid,gridm,wt,prl);

  if (prl>1) printf("\n initialize_ps_coords with third center \n");

  nphi = abs(nphi);
  int gs = nmu*nnu*nphi;
  int gs6 = 6*gs;

  double mu = 0.; double nu = 0.; double phi = 0.;
  get_munuphi(a,A3,B3,C3,mu,nu,phi);

  if (prl>1) printf("   mu/nu/phi: %8.5f %8.5f %8.5f \n",mu,nu,phi);

  if (prl>1)
  {
    double x; double y; double z;
    munuphi_to_xyz(a,mu,nu,phi,x,y,z);
    printf("  ABC: %8.5f %8.5f %8.5f  xyz: %8.5f %8.5f %8.5f \n",A3,B3,C3,x,y,z);

    if (fabs(A3-x)>1.e-6) { printf("\n ERROR: x doesn't match \n"); exit(-1); }
    if (fabs(B3-y)>1.e-6) { printf("\n ERROR: y doesn't match \n"); exit(-1); }
    if (fabs(C3-z)>1.e-6) { printf("\n ERROR: z doesn't match \n"); exit(-1); }
  }

  #pragma acc update self(gridm[0:gs6])

  int ix = 6*nnu*nphi;
  int jx = 6*nphi;

 //get the points that surround the third center
  int ib = 0; int ia = 0;
  for (int i=0;i<nmu;i++) if (gridm[ix*i+0]<mu) ib = i; else break;
  for (int i=nmu-1;i>=0;i--) if (gridm[ix*i+0]>mu) ia = i;  else break;

  int jb = 0; int ja = 0;
  for (int j=0;j<nnu;j++) if (gridm[jx*j+1]<=nu) jb = j; else break;
  for (int j=nnu-1;j>=0;j--) if (gridm[jx*j+1]>=nu) ja = j;  else break;

 //xyz is (x>0,0,z), so second line not needed
  int kb = 0; int ka = nphi-1;
  //if (fabs(phi-phi_zero-PI)<1.e-6) kb = ka = nphi / 2;

 //check if munu are near the boundaries already
  double munu_thresh = 1.e-3;

  bool skip_mu = 0;
  if (mu<munu_thresh)
    skip_mu = 1;
  else if (ib==ia)
    ia = 1;

  bool skip_nu = handle_linear_nu(jb,ja,nnu,nu,munu_thresh);

  if (skip_mu && prl>0)
  {
    printf("   nu:");
    for (int j=0;j<nnu;j++)
      printf(" %8.5f",gridm[j*jx+1]);
    printf("\n");
  }

  if (skip_nu && prl>0)
  {
    printf("   mu:");
    for (int i=0;i<nmu;i++)
      printf(" %8.5f",gridm[i*ix]);
    printf("\n");
  }


  if (prl>1) printf("   iba: %i %i  jba: %i %i  kba: %i %i \n",ib,ia,jb,ja,kb,ka);

  if ((ib==ia && !skip_mu) || (jb==ja && !skip_nu))
  {
    printf("\n ERROR: couldn't determine interval in 3-center PS \n");
    printf("   mu/nu/phi: %8.5f %8.5f %8.5f \n",mu,nu,phi);
    double x; double y; double z;
    munuphi_to_xyz(a,mu,nu,phi,x,y,z);
    printf("   ABC: %8.5f %8.5f %8.5f  xyz: %8.5f %8.5f %8.5f \n",A3,B3,C3,x,y,z);
    printf("   iba: %i %i jba: %i %i \n",ib,ia,jb,ja);

    printf("   AC1: %8.5f %8.5f %8.5f \n",coordm[0],coordm[1],coordm[2]);
    printf("   AC2: %8.5f %8.5f %8.5f \n",coordm[3],coordm[4],coordm[5]);
    printf("   AC3: %8.5f %8.5f %8.5f \n",coordm[6],coordm[7],coordm[8]);

    printf("   mu:");
    for (int i=0;i<nmu;i++)
      printf(" %8.5f",gridm[i*ix]);
    printf("\n");

    printf("   nu:");
    for (int j=0;j<nnu;j++)
      printf(" %8.5f",gridm[j*jx+1]);
    printf("\n");

    printf("   phi:");
    for (int k=0;k<nphi;k++)
      printf(" %8.5f",gridm[6*k+2]);
    printf("\n");

    exit(-1);
  }

 #if 0
  if (skip_mu || skip_nu)
  {
    printf("   skipmu/nu: %i %i \n",skip_mu,skip_nu);
    printf("   mu: %8.5f is between %8.5f and %8.5f \n",mu,gridm[ix*ib+0],gridm[ix*ia+0]);
    printf("   nu: %8.5f is between %8.5f and %8.5f \n",nu,gridm[jx*jb+1],gridm[jx*ja+1]);
    printf("  phi: %8.5f is between %8.5f and %8.5f \n",phi,gridm[6*kb+2],gridm[6*ka+2]);
  }
 #endif
  

  if (skip_mu && skip_nu)
  {
    if (prl>0)
      printf("   mu/nu are small \n");
  }
  else if (skip_mu && !skip_nu)
    shift_boundary_nu(nu,ix,jx,jb,ja,kb,ka,gs,gridm);
  else if (!skip_mu && skip_nu)
    shift_boundary_mu(mu,ix,jx,ib,ia,jb,kb,ka,gs,gridm);
  else
    shift_boundary_munu(mu,nu,ix,jx,ib,ia,jb,ja,kb,ka,gs,gridm);


  //CPMZ be careful: should uniquely identify all 8 cells
  //need to watch the near-linear angles

  //move shifted grid points to end
  {
    double tmp[6];
    #pragma acc enter data create(tmp[0:6])

    #pragma acc serial present(gridm[0:gs6],tmp[0:6])
    for (int k=0;k<2;k++)
    {
      int ws = gs-8+k*4;
      int kn = kb; if (k==1) kn = ka;

      int in = ib*ix+jb*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;

      in = ib*ix+ja*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;

      in = ia*ix+jb*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;

      in = ia*ix+ja*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;
    }

    #pragma acc exit data delete(tmp[0:6])
  }

 //order munuphi so third atom is "first"
  {
    if (prl>1) printf("   reordering end points of 3c grid (mnp: %8.5f %8.5f %8.5f) \n",mu,nu,phi);

    int i0 = gs6-48;
   #pragma acc parallel loop independent present(gridm[0:gs6])
    for (int j=0;j<8;j++)
    {
      int i1 = i0+6*j;
      double mup = gridm[i1+0];
      double nup = gridm[i1+1];
      double phip = gridm[i1+2];
      double dmup = gridm[i1+3];
      double dnup = gridm[i1+4];
      double dphip = gridm[i1+5];

      double mu1 = mup-0.5*dmup; double mu2 = mu1 + dmup;
      double nu1 = nup-0.5*dnup; double nu2 = nu1 + dnup;
      double phi1 = phip; double phi2 = phi1 + dphip;

      int sm = 1;
      if (fabs(mu1-mu)>fabs(mu2-mu))
      {
        double t1 = mu1;
        mu1 = mu2;
        mu2 = t1;
        sm = -1;
      }
      int sn = 1;
      if (fabs(nu1-nu)>fabs(nu2-nu))
      {
        double t1 = nu1;
        nu1 = nu2;
        nu2 = t1;
        sn = -1;
      }
      int sp = 1;
      double p1d = fabs(phi1-phi); if (p1d>3.14) p1d = fabs(phi1-phi-2.*PI);
      double p2d = fabs(phi2-phi); if (p2d>3.14) p2d = fabs(phi2-phi-2.*PI);
      if (p1d>p2d)
      {
        double t1 = phi1;
        phi1 = phi2;
        phi2 = t1;
        sp = -1;
      }

      double mun = mup;
      double dmun = sm*dmup;
      double nun = nup;
      double dnun = sn*dnup;
      double phin = phi1;
      double dphin = sp*dphip;

      //printf("  munuphi: %8.5f %8.5f %8.5f / %8.5f %8.5f %8.5f  ordered: %8.5f %8.5f %8.5f / %8.5f %8.5f %8.5f \n",mup-0.5*dmup,nup-0.5*dnup,phip,mup+0.5*dmup,nup+0.5*dnup,phip+dphip,mun-0.5*dmun,nun-0.5*dnun,phin,mun+0.5*dmun,nun+0.5*dnun,phin+dphin);

      gridm[i1+0] = mun;
      gridm[i1+1] = nun;
      gridm[i1+2] = phin;
      gridm[i1+3] = dmun;
      gridm[i1+4] = dnun;
      gridm[i1+5] = dphin;
    }

  }

 #if 0
  {
    #pragma acc update self(gridm[0:gs6])
    printf("\n new grid xyz: \n");
    for (int k=0;k<2;k++)
    {
      int kn = kb; if (k==1) kn = ka;

      int in = ib*ix+jb*jx+6*kn;
      mu  = gridm[in+0];
      nu  = gridm[in+1];
      phi = gridm[in+2];

      munuphi_to_xyz(a,mu,nu,phi,x,y,z);
      printf("   %8.5f %8.5f %8.5f \n",x,y,z);

      in = ib*ix+ja*jx+6*kn;
      mu  = gridm[in+0];
      nu  = gridm[in+1];
      phi = gridm[in+2];

      munuphi_to_xyz(a,mu,nu,phi,x,y,z);
      printf("   %8.5f %8.5f %8.5f \n",x,y,z);

      in = ia*ix+jb*jx+6*kn;
      mu  = gridm[in+0];
      nu  = gridm[in+1];
      phi = gridm[in+2];

      munuphi_to_xyz(a,mu,nu,phi,x,y,z);
      printf("   %8.5f %8.5f %8.5f \n",x,y,z);

      in = ia*ix+ja*jx+6*kn;
      mu  = gridm[in+0];
      nu  = gridm[in+1];
      phi = gridm[in+2];

      munuphi_to_xyz(a,mu,nu,phi,x,y,z);
      printf("   %8.5f %8.5f %8.5f \n",x,y,z);
    }
  }
 #endif

 //wt will be handled by quadrature into xyz

  return;
}

void initialize_ps_coords_4c(double cf, int nmu, int nnu, int nphi, double phi_zero, double* coordn, double* grid, double* gridm, double* wt, double* rot, int prl)
{
 //this code assumes first two atoms aligned along z axis, after the rotation
  double coordm[12];
  rotate_4x3t_cpu(rot,coordn,coordm);

  double A2 = coordm[3]; double B2 = coordm[4]; double C2 = coordm[5];
  double A3 = coordm[6]; double B3 = coordm[7]; double C3 = coordm[8];
  double A4 = coordm[9]; double B4 = coordm[10]; double C4 = coordm[11];
  double a = sqrt(A2*A2+B2*B2+C2*C2)*0.5;

  double mu3 = 0.; double nu3 = 0.; double phi3 = 0.;
  double mu4 = 0.; double nu4 = 0.; double phi4 = 0.;
  get_munuphi(a,A3,B3,C3,mu3,nu3,phi3);
  get_munuphi(a,A4,B4,C4,mu4,nu4,phi4);

  if (phi4!=0.)
    initialize_ps_coords_2c_phi(a,cf,nmu,nnu,nphi,phi_zero,phi4,grid,gridm,wt,prl);
  else
    initialize_ps_coords_2c(a,cf,nmu,nnu,nphi,phi_zero,grid,gridm,wt,prl);

  if (prl>1) printf("\n initialize_ps_coords with third center \n");

  nphi = abs(nphi);
  int gs = nmu*nnu*nphi;
  int gs6 = 6*gs;

  if (prl>1)
  {
    printf("   mu/nu/phi(3): %8.5f %8.5f %8.5f \n",mu3,nu3,phi3);
    printf("   mu/nu/phi(4): %8.5f %8.5f %8.5f \n",mu4,nu4,phi4);
  }

  if (prl>0)
  {
    double x; double y; double z;
    munuphi_to_xyz(a,mu3,nu3,phi3,x,y,z);
    printf("  ABC3: %8.5f %8.5f %8.5f  xyz: %8.5f %8.5f %8.5f \n",A3,B3,C3,x,y,z);
    munuphi_to_xyz(a,mu4,nu4,phi4,x,y,z);
    printf("  ABC4: %8.5f %8.5f %8.5f  xyz: %8.5f %8.5f %8.5f \n",A4,B4,C4,x,y,z);
  }

  #pragma acc update self(gridm[0:gs6])

  int ix = 6*nnu*nphi;
  int jx = 6*nphi;
  int kx = 6;

 //get the points that surround the third center
  int ib3 = 0; int ia3 = 0;
  for (int i=0;i<nmu;i++) if (gridm[ix*i+0]<mu3) ib3 = i; else break;
  for (int i=nmu-1;i>=0;i--) if (gridm[ix*i+0]>mu3) ia3 = i;  else break;

  int jb3 = 0; int ja3 = 0;
  for (int j=0;j<nnu;j++) if (gridm[jx*j+1]<=nu3) jb3 = j; else break;
  for (int j=nnu-1;j>=0;j--) if (gridm[jx*j+1]>=nu3) ja3 = j;  else break;

  int kb3 = 0; int ka3 = nphi-1;


 //get the points that surround the fourth center
  int ib4 = 0; int ia4 = 0;
  for (int i=0;i<nmu;i++) if (gridm[ix*i+0]<mu4) ib4 = i; else break;
  for (int i=nmu-1;i>=0;i--) if (gridm[ix*i+0]>mu4) ia4 = i;  else break;

  int jb4 = 0; int ja4 = 0;
  for (int j=0;j<nnu;j++) if (gridm[jx*j+1]<=nu4) jb4 = j; else break;

  for (int j=nnu-1;j>=0;j--) if (gridm[jx*j+1]>=nu4) ja4 = j;  else break;

  int kba = 0;
  for (int k=0;k<nphi;k++) if (fabs(gridm[kx*k+2]-phi4)<1.e-6) { kba = k; break; }
  int kb4 = kba-1; int ka4 = kba;
  if (kb4<0) kb4 = nphi-1;

 //check if munu are near the boundaries already
  double munu_thresh = 1.e-3;

  bool skip_mu3 = 0;
  if (mu3<munu_thresh)
    skip_mu3 = 1;
  else if (ib3==ia3)
    ia3 = 1;

  bool skip_nu3 = handle_linear_nu(jb3,ja3,nnu,nu3,munu_thresh);

  if (skip_mu3 && prl>0)
  {
    printf("   nu:");
    for (int j=0;j<nnu;j++)
      printf(" %8.5f",gridm[j*jx+1]);
    printf("\n");
  }

  if (skip_nu3 && prl>0)
  {
    printf("   mu:");
    for (int i=0;i<nmu;i++)
      printf(" %8.5f",gridm[i*ix]);
    printf("\n");
  }


  if (prl>0)
  {
    printf("   iba(3): %2i %2i  jba: %2i %2i  kba: %2i %2i \n",ib3,ia3,jb3,ja3,kb3,ka3);
    printf("   iba(4): %2i %2i  jba: %2i %2i  kba: %2i %2i \n",ib4,ia4,jb4,ja4,kb4,ka4);
  }

  bool i34m = 0; if (ib3==ib4 && ia3==ia4) i34m = 1;
  bool j34m = 0; if (jb3==jb4 && ja3==ja4) j34m = 1;
  bool k34m = 0; if (kb3==kb4 && ka3==ka4) k34m = 1;

  if ((ib3==ia3 && !skip_mu3) || (jb3==ja3 && !skip_nu3) || (i34m && j34m && k34m) )
  {
    printf("\n ERROR: couldn't determine interval in 3-center PS \n");
    printf("   mu/nu/phi(3): %8.5f %8.5f %8.5f \n",mu3,nu3,phi3);
    double x; double y; double z;
    munuphi_to_xyz(a,mu3,nu3,phi3,x,y,z);
    printf("   ABC: %8.5f %8.5f %8.5f  xyz: %8.5f %8.5f %8.5f \n",A3,B3,C3,x,y,z);
    printf("   iba3: %i %i jba3: %i %i \n",ib3,ia3,jb3,ja3);
    printf("   iba4: %i %i jba4: %i %i \n",ib4,ia4,jb4,ja4);

    printf("   AC1: %8.5f %8.5f %8.5f \n",coordm[0],coordm[1],coordm[2]);
    printf("   AC2: %8.5f %8.5f %8.5f \n",coordm[3],coordm[4],coordm[5]);
    printf("   AC3: %8.5f %8.5f %8.5f \n",coordm[6],coordm[7],coordm[8]);

    printf("   mu:");
    for (int i=0;i<nmu;i++)
      printf(" %8.5f",gridm[i*ix]);
    printf("\n");

    printf("   nu:");
    for (int j=0;j<nnu;j++)
      printf(" %8.5f",gridm[j*jx+1]);
    printf("\n");

    printf("   phi:");
    for (int k=0;k<nphi;k++)
      printf(" %8.5f",gridm[6*k+2]);
    printf("\n");

    exit(-1);
  }

 #if 0
  printf("   mu3: %8.5f is between %8.5f and %8.5f \n",mu3,gridm[ix*ib+0],gridm[ix*ia+0]);
  printf("   nu3: %8.5f is between %8.5f and %8.5f \n",nu3,gridm[jx*jb+1],gridm[jx*ja+1]);
  printf("  phi3: %8.5f is between %8.5f and %8.5f \n",phi3,gridm[6*kb+2],gridm[6*ka+2]);
  printf("   mu4: %8.5f is between %8.5f and %8.5f \n",mu4,gridm[ix*ib+0],gridm[ix*ia+0]);
  printf("   nu4: %8.5f is between %8.5f and %8.5f \n",nu4,gridm[jx*jb+1],gridm[jx*ja+1]);
  printf("  phi4: %8.5f is between %8.5f and %8.5f \n",phi4,gridm[6*kb+2],gridm[6*ka+2]);
 #endif

  if (skip_mu3 && skip_nu3)
  {
    if (prl>0)
      printf("   mu/nu are small \n");
  }
  else if (skip_mu3 && !skip_nu3)
    shift_boundary_nu(nu3,ix,jx,jb3,ja3,kb3,ka3,gs,gridm);
  else if (!skip_mu3 && skip_nu3)
    shift_boundary_mu(mu3,ix,jx,ib3,ia3,jb3,kb3,ka3,gs,gridm);
  else
    shift_boundary_munu(mu3,nu3,ix,jx,ib3,ia3,jb3,ja3,kb3,ka3,gs,gridm);

  bool skip_mu4 = 0;
  if (mu4<munu_thresh)
    skip_mu4 = 1;
  else if (ib4==ia4)
    ia4 = 1;

  bool skip_nu4 = handle_linear_nu(jb4,ja4,nnu,nu4,munu_thresh);

  if (skip_mu4 && skip_nu4)
  {
    if (prl>0)
      printf("   mu/nu are small \n");
  }
  else if (skip_mu4 && !skip_nu4)
    shift_boundary_nu(nu4,ix,jx,jb4,ja4,kb4,ka4,gs,gridm);
  else if (!skip_mu4 && skip_nu4)
    shift_boundary_mu(mu4,ix,jx,ib4,ia4,jb4,kb4,ka4,gs,gridm);
  else
    shift_boundary_munu(mu4,nu4,ix,jx,ib4,ia4,jb4,ja4,kb4,ka4,gs,gridm);

  //move shifted grid points to end
  {
    double tmp[6];
    #pragma acc enter data create(tmp[0:6])

    #pragma acc serial present(gridm[0:gs6],tmp[0:6])
    for (int k=0;k<2;k++)
    {
      int ws = gs-16+k*4;
      int kn = kb3; if (k==1) kn = ka3;

      int in = ib3*ix+jb3*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;

      in = ib3*ix+ja3*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;

      in = ia3*ix+jb3*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;

      in = ia3*ix+ja3*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;
    }

    #pragma acc serial present(gridm[0:gs6],tmp[0:6])
    for (int k=0;k<2;k++)
    {
      int ws = gs-8+k*4;
      int kn = kb4; if (k==1) kn = ka4;

      int in = ib4*ix+jb4*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;

      in = ib4*ix+ja4*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;

      in = ia4*ix+jb4*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;

      in = ia4*ix+ja4*jx+6*kn;
      for (int j=0;j<6;j++)
      {
        tmp[j] = gridm[6*ws+j];
        gridm[6*ws+j] = gridm[in+j];
        gridm[in+j] = tmp[j];
      }
      ws++;
    }

    #pragma acc exit data delete(tmp[0:6])
  }


  return;
}

void initialize_ps_coords_batch(int wb, int nbatch, double a, double cf, int nmu, int nnu, int nphi, double phi_zero, double* grid, double* gridm, double* wt, int prl)
{
  if (a<=0.) { printf(" ERROR: initialize_ps_coords with a<0 \n",a); exit(-1); }
  //if (nnu<2 || nphi<2) { printf(" ERROR: cannot initialize PS coordinates with nnu<2 or nphi<2 \n"); exit(-1); }

  double c0 = get_c0(a,nmu,cf,RMAX);
  if (cf<0.) c0 = -cf; //c1 passed to ftn
  const double c1 = c0;

  int gs = nmu*nnu*abs(nphi)/nbatch;
  int gs6 = 6*gs;

  const double dnu = PI/nnu;
  double dphi0 = 2.*PI/nphi;
  if (nphi==-1) dphi0 = 1.;
  nphi = abs(nphi);
  const double dphi = dphi0;
  const double dx1 = 1./(nmu+1);
  double a3 = a*a*a;

  if (prl>0)
    printf("\n init_ps_coords. wb/nbatch: %i %i -- a: %8.5f nmu/nu/phi: %2i %2i %2i  dnu/dphi: %8.5f %8.5f  dx1: %8.5f \n",wb,nbatch,a,nmu,nnu,nphi,dnu,dphi,dx1);
  double rmax = a*cosh(c1*atanh(nmu*dx1));
  printf("    cf: %5.1f  c1: %8.5f  rmax: %8.3f \n",cf,c1,rmax);

  int ic = 0;
  for (int i=0;i<nmu;i++)
  if (i%nbatch==wb)
  {
    double x0  =  i*dx1;
    double x1  = x0+dx1;
    double mu0 = c1*atanh(x0);
    double mu1 = c1*atanh(x1);
    double mu = 0.5*(mu1+mu0); //"average" mu value

    double sinhm = sinh(mu);
    double coshm = cosh(mu);
    double dmu = mu1-mu0;

    //printf("  wb: %i  mu: %8.5f \n",wb,mu);

    for (int j=0;j<nnu;j++)
    {
      double nu = (j+0.5)*dnu;

      double sinn = sin(nu);
      double cosn = cos(nu);

      int i1 = ic*nnu*nphi + j*nphi;
     #pragma acc parallel loop present(grid[0:gs6],gridm[0:gs6],wt[0:gs])
      for (int k=0;k<nphi;k++)
      {
        int i2 = i1+k;
        double phi = k*dphi + phi_zero;

        gridm[6*i2+0] = mu; gridm[6*i2+1] = nu; gridm[6*i2+2] = phi;
        gridm[6*i2+3] = dmu; gridm[6*i2+4] = dnu; gridm[6*i2+5] = dphi;
      } //loop phi
    } //loop nu
    ic++;
  } //loop mu

  //munuphi_to_xyz(a,gs,gridm,grid,wt);

 #if 0
  #pragma acc update self(gridm[0:gs6])
  for (int i=0;i<nmu;i++)
  {
    double mu1 = gridm[6*i*nnu*nphi+0];
    double r1 = a*cosh(mu1);
    printf("  mu: %8.5f  r: %5.2f \n",mu1,r1);
  }
 #endif

 #if 0
  ic = 0;
  if (grid!=NULL && wt!=NULL)
  for (int i=0;i<nmu;i++)
  if (i%nbatch==wb)
  {
    int i1 = ic*nnu*nphi;
    for (int j=0;j<nnu*nphi;j++)
    {
      double mu  = gridm[i1+6*j+0];
      double nu  = gridm[i1+6*j+1];
      double phi = gridm[i1+6*j+2];
      double dmu = gridm[i1+6*j+3];
      double dnu = gridm[i1+6*j+4];

      double sinhm = sinh(mu);
      double coshm = cosh(mu);
      double sinn = sin(nu);
      double cosn = cos(nu);

      double x = a*sinhm*sinn*cos(phi);
      double y = a*sinhm*sinn*sin(phi);
      double z = a*coshm*cosn;

      double wt1 = a3*ps_dV(mu-0.5*dmu,mu+0.5*dmu,nu-0.5*dnu,nu+0.5*dnu) * dphi;

      //if (x>xmax) xmax = x;
      //if (y>ymax) ymax = y;
      //if (z>zmax) zmax = z;

      grid[i1+6*j+0] = x;
      grid[i1+6*j+1] = y;
      grid[i1+6*j+2] = z-a;
      wt[j] = wt1;
    }
    ic++;
  }
 #endif

  return;
}

double get_rsp(int nsplit, double ztm, double* rsp, int prl)
{
  int c0 = 0; if (nsplit==2) c0 = 1; 
  for (int j=0;j<nsplit-1;j++)
  {
    double cut = 0.5*pow(10.,-(c0+j));
    double rc = -log(cut)/ztm;

    rsp[j] = rc;

    if (prl>1)
      printf("  cut: %8.5f  rc: %8.5f \n",cut,rc);
  }

  return rsp[nsplit-2];
}

void refine_3c_grid(int nsplit, double ztm, double a, double A3, double B3, double C3, int gs, int nsg, double* gridm, double* rot, int prl)
{
  //if (nsplit<2) return;
  if (nsg<1) return;

  if (prl>1) printf("   refine_3c_grid (gs: %4i  nsplit: %2i  nsg: %2i) \n",gs,nsplit,nsg);

 //put third center in PS plane
  rotate_3xyz(A3,B3,C3,rot);
  //printf("   ABC: %8.5f %8.5f %8.5f \n",A3,B3,C3);

  int gs0 = gs-nsg;
  int gs6 = 6*gs;
  int i0 = 6*gs0;
  int ns3 = nsplit*nsplit*nsplit;

  double* gride = new double[48];
  #pragma acc enter data create(gride[0:48])

  #pragma acc parallel loop collapse(2) present(gridm[0:gs6],gride[0:48])
  for (int j=0;j<8;j++)
  for (int k=0;k<6;k++)
    gride[6*j+k] = gridm[i0+6*j+k];

 //not using rsp/rsmax
  //double* rsp = new double[nsplit-1];
  //double rsmax = get_rsp(nsplit,ztm,rsp,prl);
  //#pragma acc enter data copyin(rsp[0:nsplit-1])

 #pragma acc parallel loop independent present(gridm[0:gs6],gride[0:48])
  for (int j=0;j<8;j++)
  {
    double   mu = gride[6*j+0];
    double   nu = gride[6*j+1];
    double  phi = gride[6*j+2];
    double  dmu = gride[6*j+3];
    double  dnu = gride[6*j+4];
    double dphi = gride[6*j+5];

    double mu1 = mu-0.5*dmu; double mu2 = mu1+dmu;
    double nu1 = nu-0.5*dnu; double nu2 = nu1+dnu;
    double phi1 = phi; double phi2 = phi+dphi;
    if (phi1>=2.*PI) phi1 -= 2.*PI;

    double x,y,z;
    munuphi_to_xyz_gpu(a,mu1,nu1,phi1,x,y,z);
    double x1 = x-A3; double y1 = y-B3; double z1 = z-C3;

   //distance to third center (first grid point)
    double R1 = sqrt(x1*x1+y1*y1+z1*z1);

    munuphi_to_xyz_gpu(a,mu2,nu1,phi1,x,y,z);
    double x2 = x-A3; double y2 = y-B3; double z2 = z-C3;
    munuphi_to_xyz_gpu(a,mu1,nu2,phi1,x,y,z);
    double x3 = x-A3; double y3 = y-B3; double z3 = z-C3;
    munuphi_to_xyz_gpu(a,mu1,nu1,phi2,x,y,z);
    double x4 = x-A3; double y4 = y-B3; double z4 = z-C3;
    double Rmu  = sqrt(x2*x2+y2*y2+z2*z2);
    double Rnu  = sqrt(x3*x3+y3*y3+z3*z3);
    double Rphi = sqrt(x4*x4+y4*y4+z4*z4);

    munuphi_to_xyz_gpu(a,mu1,nu1,phi1,x,y,z);
    if (prl>1)
      printf("  munuphi1: %8.5f %8.5f %8.5f  xyz: %8.5f %8.5f %8.5f  Rs: %8.5f %8.5f %8.5f %8.5f \n",mu1,nu1,phi1,x,y,z,R1,Rmu,Rnu,Rphi);

    //if (Rmu<rsmax || Rnu<rsmax || Rphi<rsmax)
    {
      double fs = 1./nsplit;
      double dmun = fs*dmu;
      double dnun = fs*dnu;
      double dphin = fs*dphi;

      int i1 = i0+6*j*ns3;
     #pragma acc loop independent
      for (int m=0;m<nsplit;m++)
      {
        double fm = m*fs;
        double mun = mu1+fm*dmu;
        int i2 = i1 + 6*m*nsplit*nsplit;

       #pragma acc loop independent
        for (int n=0;n<nsplit;n++)
        {
          double fn = n*fs;
          double nun = nu1+fn*dnu;
          int i3 = i2 + 6*n*nsplit;

         #pragma acc loop independent
          for (int p=0;p<nsplit;p++)
          {
            double fp = p*fs;
            double phin = phi1+fp*dphi;
            int i4 = i3 + 6*p;

            //printf("  fmnp: %8.5f %8.5f %8.5f  mnp: %8.5f %8.5f %8.5f \n",fm,fn,fp,mun,nun,phin);

            gridm[i4+0] = mun+0.5*dmun;
            gridm[i4+1] = nun+0.5*dnun;
            gridm[i4+2] = phin;
            gridm[i4+3] = dmun;
            gridm[i4+4] = dnun;
            gridm[i4+5] = dphin;
          }
        }
      }
    } //cell is small, divide up evenly
   #if 0
    else
    {
      //cell division based on size of 3rd center
      printf("\n ERROR: cell division type 2 not available \n");
    }
   #endif

  } //loop j over octants

  if (prl>1)
  {
    #pragma acc update self(gridm[0:gs6])
    printf("  split grid pts: \n");
    for (int j=0;j<8*ns3;j++)
    {
      double   mu = gridm[i0+6*j+0];
      double   nu = gridm[i0+6*j+1];
      double  phi = gridm[i0+6*j+2];
      double  dmu = gridm[i0+6*j+3];
      double  dnu = gridm[i0+6*j+4];
      double dphi = gridm[i0+6*j+5];

      double mu1 = mu-0.5*dmu; double mu2 = mu1+dmu;
      double nu1 = nu-0.5*dnu; double nu2 = nu1+dnu;
      double phi1 = phi; double phi2 = phi+dphi;
      if (phi1>=2.*PI) phi1 -= 2.*PI;

      double x,y,z;
      munuphi_to_xyz(a,mu1,nu1,phi1,x,y,z);
      printf(" xyz(1): %8.5f %8.5f %8.5f ",x,y,z);

      double x1 = x-A3; double y1 = y-B3; double z1 = z-C3;
      double R1 = sqrt(x1*x1+y1*y1+z1*z1);

      munuphi_to_xyz(a,mu2,nu1,phi1,x,y,z);
      printf(" xyz(mu2): %8.5f %8.5f %8.5f ",x,y,z);

      double x2 = x-A3; double y2 = y-B3; double z2 = z-C3;
      munuphi_to_xyz(a,mu1,nu2,phi1,x,y,z);
      printf(" xyz(nu2): %8.5f %8.5f %8.5f ",x,y,z);

      double x3 = x-A3; double y3 = y-B3; double z3 = z-C3;
      munuphi_to_xyz(a,mu1,nu1,phi2,x,y,z);
      printf(" xyz(phi2): %8.5f %8.5f %8.5f ",x,y,z);

      double x4 = x-A3; double y4 = y-B3; double z4 = z-C3;
      double Rmu  = sqrt(x2*x2+y2*y2+z2*z2);
      double Rnu  = sqrt(x3*x3+y3*y3+z3*z3);
      double Rphi = sqrt(x4*x4+y4*y4+z4*z4);

      //munuphi_to_xyz(a,mu1,nu1,phi1,x,y,z);
      if (R1<0.001)
        printf("  mnp1: %8.5f %8.5f %8.5f  mnp2: %8.5f %8.5f %8.5f  Rs: %8.5f %8.5f %8.5f %8.5f \n",mu1,nu1,phi1,mu2,nu2,phi2,R1,Rmu,Rnu,Rphi);
    }
  }

 #if 0
 #pragma acc parallel loop present(gridm[0:gs6])
  for (int j=gs0+8;j<gs;j++)
  {
    for (int k=0;k<6;k++)
      gridm[6*j+k] = 0.;
  }
 #endif

  #pragma acc exit data delete(gride[0:48])
  //#pragma acc exit data delete(rsp[0:nsplit-1])

  //delete [] rsp;
  delete [] gride;

  return;
}

double get_cfn_3c(double ztm1, double ztm2, int nmu, double* coordn)
{
  double A2 = coordn[3]; double B2 = coordn[4]; double C2 = coordn[5];
  double A3 = coordn[6]; double B3 = coordn[7]; double C3 = coordn[8];
  double A23 = A3-A2; double B23 = B3-B2; double C23 = C3-C2;
  double z0 = sqrt(A2*A2+B2*B2+C2*C2)*0.5;

  double r12 = 2.*z0;
  double r13 = sqrt(A3*A3+B3*B3+C3*C3);
  double r23 = sqrt(A23*A23+B23*B23+C23*C23);

 //parameters
  double ln14 = log(1.e-14);
  double RMIN = 4.; //lower limit

  double rm = -ln14/(ztm1+ztm2); //effective size of 2 center Slaters
  if (r13<r23) //center 3 closer to center 1
  {
    if (rm<r13) rm = r13+1.;
  }
  else //center 3 closer to center 2
  {
    if (rm<r23) rm = r23+1.;
  }
  if (rm<RMIN)
    rm = RMIN;

  double mu = acosh(rm/z0);
  double fm = nmu/(nmu+1.);
  double cfn = mu/atanh(fm);

  //printf("    get_cfn_3c.  r12/13/23: %8.5f %8.5f %8.5f  ztm: %8.5f %8.5f  rm: %8.5f  mu: %8.5f  fm: %8.5f  cfn: %8.5f \n",r12,r13,r23,ztm1,ztm2,rm,mu,fm,cfn);

  //debug
  //cfn = 1.1;

  return -cfn;
}


void generate_ps_quad_grid_3c_refine(double cfn, int wb, int nb, double ztm1, double ztm2, int nsplit, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt)
{
  if (natoms<3) { printf("\n ERROR: cannot use refined grid with natoms<3 \n"); exit(-1); }
  int prl = 0;
  if (prl>1) printf("\n in generate_ps_quad_grid_3c_refine (natoms: %i) \n",natoms);
  if (wb>=nb) { printf("\n ERROR: wb cannot be larger than nb \n"); exit(-1); }

  //double cfn = 1.1; //3-atom specific
 #if VOLUME_RESIZE && 0
  cfn = get_cfn_3c(ztm1,ztm2,nmu,coordn);
 #endif

  int qo = quad_order;
  int qoh = quad_r_order;

  int qos = qo*qo*qo;
  int qosh = qoh*qoh*qoh;
  int qo2 = 2*quad_order;
  int qoh2 = 2*qoh;

  double* Qx = new double[qo2];
  double* Qy = new double[qo2];
  double* Qz = new double[qo2];

  get_quad(qo,Qx);
  get_quad(qo,Qy);
  get_quad(qo,Qz);

  #pragma acc enter data copyin(Qx[0:qo2],Qy[0:qo2],Qz[0:qo2])

  int gs0 = nmu*nnu*nphi-8;
  int nsg = 8*nsplit*nsplit*nsplit;
  int gs = gs0+nsg;
  int gsq = gs0*qos+nsg*qosh; //3-atom grid size
  if (gsq%nb>0) { printf("\n ERROR: couldn't divide up grid in 3c_refine (remainder: %i) \n",gsq%nb); exit(-1); }
  gsq /= nb;
  int gs6 = 6*gs; // w/o quad
  int gsq6 = 6*gsq; // w/ quad
  //if (nsg==8) nsg = 0; //call quad only once (debug)

  //printf("  gs: %4i  gsq: %5i  wb/nb: %i %i \n",gs,gsq,wb,nb);

  double* gridq = new double[gsq6];
  double* gridm = new double[gs6];

  #pragma acc enter data create(gridq[0:gsq6],gridm[0:gs6])

  double rot[9];
  #pragma acc enter data create(rot[0:9])
  #pragma acc parallel loop present(rot[0:9])
  for (int j=0;j<9;j++) 
    rot[j] = 0.;
  #pragma acc parallel loop present(rot[0:9])
  for (int j=0;j<3;j++) 
    rot[j*3+j] = 1.;

  double A2 = coordn[3]; double B2 = coordn[4]; double C2 = coordn[5];
  double z0 = sqrt(A2*A2+B2*B2+C2*C2)*0.5;

  {
    get_3c_position(coordn,rot);
    double A3 = coordn[6]; double B3 = coordn[7]; double C3 = coordn[8];

    initialize_ps_coords_3c(cfn,nmu,nnu,nphi,0.,coordn,NULL,gridm,NULL,rot,prl);

    refine_3c_grid(nsplit,0.,z0,A3,B3,C3,gs,nsg,gridm,rot,prl);

   //now w/batching
    quad_grid_munuphi(wb,nb,qo,qo,qo,0,z0,Qx,Qy,Qz,gs-nsg,0,gridm,gridq,wt);

   //if splitting, create rest of grid
    if (nsg>0)
    {
      double Qxh[qoh2]; double Qyh[qoh2]; double Qzh[qoh2];
      get_quad(qoh,Qxh); get_quad(qoh,Qyh); get_quad(qoh,Qzh);
      #pragma acc enter data copyin(Qxh[0:qoh2],Qyh[0:qoh2],Qzh[0:qoh2])

     //nsg pts surrounding third nucleus w/refined grid
      quad_grid_munuphi(wb,nb,qoh,qoh,qoh,qos,z0,Qxh,Qyh,Qzh,gs,gs-nsg,gridm,gridq,wt);

      #pragma acc exit data delete(Qxh[0:qoh2],Qyh[0:qoh2],Qzh[0:qoh2])
    }

    reorient_grid(z0,gsq,gridq,grid,rot);
  }

  #pragma acc exit data delete(Qx[0:qo2],Qy[0:qo2],Qz[0:qo2])
  #pragma acc exit data delete(rot[0:9])
  #pragma acc exit data delete(gridq[0:gsq6],gridm[0:gs6])

  if (0)
  {
    printf("  coordn(1/2): %8.5f %8.5f %8.5f   %8.5f %8.5f %8.5f \n",coordn[0],coordn[1],coordn[2],coordn[3],coordn[4],coordn[5]);
    printf("  grid around atom 3 (%8.5f %8.5f %8.5f) \n",coordn[6],coordn[7],coordn[8]);
    #pragma acc update self(grid[0:gsq6],wt[0:gsq])
    for (int j=gsq-nsg*qosh;j<gsq;j++)
      printf("  xyz: %8.5f %8.5f %8.5f  wt: %9.6f \n",grid[6*j+0],grid[6*j+1],grid[6*j+2],wt[j]);
  }
 
  delete [] Qx;
  delete [] Qy;
  delete [] Qz;
  delete [] gridq;
  delete [] gridm;

  return;
}

void generate_ps_quad_grid_3c_refine(int wb, int nb, double ztm1, double ztm2, int nsplit, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt)
{
  double cfn = 1.;
  return generate_ps_quad_grid_3c_refine(cfn,wb,nb,ztm1,ztm2,nsplit,natoms,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
}

void generate_ps_quad_grid_3c_refine(double ztm1, double ztm2, int nsplit, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt)
{
  double cfn = 1.;
  int wb = 0; int nb = 1;
  return generate_ps_quad_grid_3c_refine(cfn,wb,nb,ztm1,ztm2,nsplit,natoms,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
}

void get_nrad_nang_order(int& nrad, int& nang, int ang_order, int nb, int gsq)
{
  while (nrad<500)
  {
    ang_order--;
    nang = get_npts_from_order(ang_order);
    if (nang<1) { printf("\n ERROR: couldn't construct an atomic grid \n"); exit(-1); }
    nrad = gsq/nang;
  }
  if (nrad<1) { printf("\n ERROR: couldn't construct an atomic grid \n"); exit(-1); }

  if (nrad>2000)
    nrad = 2000;
  nrad -= nrad%nb;

  return;
}

//wb,nb --> batches over grid
void generate_ps_quad_grid(double cfn, int wb, int nb, double Z1, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt)
{
  int prl = 0;
  if (prl>1) printf("\n in generate_ps_quad_grid (natoms: %i) \n",natoms);
  if (nb>1 && natoms>2) printf(" TESTING: batching for 3c systems \n");
  if (wb>=nb) { printf("\n ERROR: wb cannot be larger than nb \n"); exit(-1); }

  int qo = quad_order;
  int qoh = quad_r_order;
  //qoh = quad_order;

  int qos = qo*qo*qo;
  int qosh = qoh*qoh*qoh;
  int qo2 = 2*quad_order;
  int qoh2 = 2*qoh;

  double* Qx = new double[qo2];
  double* Qy = new double[qo2];
  double* Qz = new double[qo2];

  get_quad(qo,Qx);
  get_quad(qo,Qy);
  get_quad(qo,Qz);

  #pragma acc enter data copyin(Qx[0:qo2],Qy[0:qo2],Qz[0:qo2])

  int gs = nmu*nnu*nphi;
  int gsq = qos*gs;
  int nqa = min(2,natoms-2);
  if (natoms>=3)
    gsq = (gs-nqa*8)*qos+nqa*8*qosh;

  if (natoms<=2 && nmu%nb>0) { printf("\n ERROR: batch must divide nmu \n"); exit(-1); }
  if (natoms>2 && nb>1) { printf("\n ERROR: cannot use batching here when natoms>2 \n"); exit(-1); }
  if (gs%nb>0) { printf("\n ERROR: batch must divide gs \n"); exit(-1); }
  if (gsq%nb>0) { printf("\n ERROR: batch must divide gsq \n"); exit(-1); }

  gs /= nb;
  gsq /= nb;
  int gs6 = 6*gs;
  int gsq6 = 6*gsq;

  if (natoms>=3)
    printf("    natoms: %i  gs/q: %6i %6i \n",natoms,gs,gsq);

  double* gridq = new double[gsq6];
  double* gridm = new double[gs6];

  #pragma acc enter data create(gridq[0:gsq6],gridm[0:gs6])

  double rot[9];
  #pragma acc enter data create(rot[0:9])
  #pragma acc parallel loop present(rot[0:9])
  for (int j=0;j<9;j++) 
    rot[j] = 0.;
  #pragma acc parallel loop present(rot[0:9])
  for (int j=0;j<3;j++) 
    rot[j*3+j] = 1.;

  double A2 = coordn[3]; double B2 = coordn[4]; double C2 = coordn[5];
  double z0 = sqrt(A2*A2+B2*B2+C2*C2)*0.5;
  if (natoms==1) z0 = 2.;

  #define ATOMIC_GRID 1
 
  if (natoms==1)
  {
   #if ATOMIC_GRID
    int nrad = 0;
    int nang = 0;

   //maximum degree of Lebedev grid
    int ang_order = 24;
    get_nrad_nang_order(nrad,nang,ang_order,nb,gsq*nb);
    int gsb = nrad*nang/nb;

   //zero weight on unused data pts. also need r1 to be nonzero
    #pragma acc parallel loop present(wt[0:gsq])
    for (int j=0;j<gsq;j++)
      wt[j] = 0.;
    #pragma acc parallel loop present(grid[0:gsq6])
    for (int j=0;j<gsq6;j++)
      grid[j] = 1.e4;
    if (nrad<1) { printf("\n ERROR: could not form atomic grid \n"); exit(-1); }

    printf("    using atomic grid (Z: %2i).  nrad: %4i  nang: %4i  gs: %7i  gsq: %7i \n",(int)Z1,nrad,nang,gsb,gsq);
    generate_central_grid_2d(wb,nb,1,grid,wt,Z1,nrad,nang);

   #else
   //realigned PS grid
    //double c2 = read_float_2("C2");
    //if (c2<=-1.) c2 = 0.;
    //if (c2<-0.5) { c2 = -0.5; printf("\n WARNING: C2 cannot be less than -0.5 \n"); } 
    //if (c2>0.5)  { c2 =  0.5; printf("\n WARNING: C2 cannot be larger than 0.5 \n"); } 

    //initialize_ps_coords_1c(z0,cfn,c2,nmu,nnu,nphi,0.,NULL,gridm,NULL,prl);
    initialize_ps_coords_batch(wb,nb,z0,cfn,nmu,nnu,nphi,0.,NULL,gridm,NULL,prl);

    quad_grid_munuphi(qo,qo,qo,0,z0,Qx,Qy,Qz,gs,0,gridm,gridq,wt);

    #pragma acc parallel loop present(gridq[0:gsq6],grid[0:gsq6])
    for (int j=0;j<gsq6;j++)
      grid[j] = gridq[j];
    //shift_grid(z0,gsq,grid);
   #endif
  }
  else if (natoms==2)
  {
    get_2c_position(coordn,rot);

    initialize_ps_coords_batch(wb,nb,z0,cfn,nmu,nnu,nphi,0.,NULL,gridm,NULL,prl);
    //initialize_ps_coords_2c(z0,cfn,nmu,nnu,nphi,0.,NULL,gridm,NULL,prl);

    quad_grid_munuphi(qo,qo,qo,0,z0,Qx,Qy,Qz,gs,0,gridm,gridq,wt);

    reorient_grid(z0,gsq,gridq,grid,rot);
  }
  else if (natoms==3)
  {
    get_3c_position(coordn,rot);

    initialize_ps_coords_3c(cfn,nmu,nnu,nphi,0.,coordn,NULL,gridm,NULL,rot,prl);

    quad_grid_munuphi(qo,qo,qo,0,z0,Qx,Qy,Qz,gs-8,0,gridm,gridq,wt);

    double Qxh[qoh2]; double Qyh[qoh2]; double Qzh[qoh2];
    get_quad(qoh,Qxh); get_quad(qoh,Qyh); get_quad(qoh,Qzh);
    #pragma acc enter data copyin(Qxh[0:qoh2],Qyh[0:qoh2],Qzh[0:qoh2])

   //8 pts surrounding third nucleus w/refined grid
    quad_grid_munuphi(wb,nb,qoh,qoh,qoh,qos,z0,Qxh,Qyh,Qzh,gs,gs-8,gridm,gridq,wt);

    #pragma acc exit data delete(Qxh[0:qoh2],Qyh[0:qoh2],Qzh[0:qoh2])

    reorient_grid(z0,gsq,gridq,grid,rot);
  }
  else if (natoms==4)
  {
    if (nb>1) { printf("\n ERROR: shouldn't be here (4-atom grid) with nbatch>1 \n"); exit(-1); }

   //first rotation matrix that takes three atoms into xz plane
    get_4c_position(coordn,rot);

    //initialize_ps_coords_3c(cfn,nmu,nnu,nphi,0.,coordn,NULL,gridm,NULL,rot,prl);
    initialize_ps_coords_4c(cfn,nmu,nnu,nphi,0.,coordn,NULL,gridm,NULL,rot,prl);

    quad_grid_munuphi(qo,qo,qo,0,z0,Qx,Qy,Qz,gs-16,0,gridm,gridq,wt);

    double Qxh[qoh2]; double Qyh[qoh2]; double Qzh[qoh2];
    get_quad(qoh,Qxh); get_quad(qoh,Qyh); get_quad(qoh,Qzh);
    #pragma acc enter data copyin(Qxh[0:qoh2],Qyh[0:qoh2],Qzh[0:qoh2])

   //16 pts surrounding latter two nuclei w/refined grid
    quad_grid_munuphi(qoh,qoh,qoh,qos,z0,Qxh,Qyh,Qzh,gs,gs-16,gridm,gridq,wt);

    #pragma acc exit data delete(Qxh[0:qoh2],Qyh[0:qoh2],Qzh[0:qoh2])

    reorient_grid(z0,gsq,gridq,grid,rot);
  }
  else
  {
    printf("\n ERROR: generate_ps_grid_quad supports natoms<5 \n"); exit(-1);
  }

  #pragma acc exit data delete(Qx[0:qo2],Qy[0:qo2],Qz[0:qo2])
  #pragma acc exit data delete(rot[0:9])
  #pragma acc exit data delete(gridq[0:gsq6],gridm[0:gs6])

 #if 1
  if (natoms==1 && gsq<10000 && prl > 1)
  {
    #pragma acc update self(grid[0:gsq6],wt[0:gsq])
    for (int j=0;j<gsq;j++)
      printf("  xyz: %8.5f %8.5f %8.5f  wt: %9.6f \n",grid[6*j+0],grid[6*j+1],grid[6*j+2],wt[j]);
  }
 #endif

  delete [] Qx;
  delete [] Qy;
  delete [] Qz;
  delete [] gridq;
  delete [] gridm;

  return;
}

void generate_ps_quad_grid(double cfn, double Z1, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt)
{
  int wb = 0; int nb = 1;
  return generate_ps_quad_grid(cfn,wb,nb,Z1,natoms,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
}

void generate_ps_quad_grid(double Z1, int natoms, double* coordn, int quad_order, int quad_r_order, int nmu, int nnu, int nphi, double* grid, double* wt)
{
  double cfn = 1.;
  int wb = 0; int nb = 1;
  return generate_ps_quad_grid(cfn,wb,nb,Z1,natoms,coordn,quad_order,quad_r_order,nmu,nnu,nphi,grid,wt);
}
