#include "symm.h"

#define THRESH_SYMM 1.e-4

//Notes:
// 1. need C2v
// 2. check for closeness between +/- values
// 3. 4H tzp --> one orbital has no symmetry

void generate_mo(int i1, int gsa, int iN, int N, float* norm, float* tm1, float** val1, double* jCA1)
{
  int N2 = N*N;

  #pragma acc parallel loop present(tm1[0:gsa])
  for (int j=0;j<gsa;j++)
    tm1[j] = 0.f;

   //tm1 is MO on grid
 #pragma acc parallel loop present(tm1[0:gsa],val1[0:iN][0:gsa],norm[0:N],jCA1[0:N2])
  for (int j=0;j<gsa;j++)
  {
    float v1 = 0.f;
    #pragma acc loop reduction(+:v1)
    for (int k=0;k<N;k++)
      v1 += jCA1[k*N+i1]*norm[k]*val1[k][j];
    tm1[j] = v1;
  }

  #pragma acc update self(tm1[0:gsa])

  return;
}

void eval_cs_mos(int natoms, int gsa, int N, float* norm, float* tm1, float** val1, double* jCA1, int prl, int* symm)
{
  for (int i1=0;i1<N;i1++)
  {
    generate_mo(i1,gsa,N,N,norm,tm1,val1,jCA1);

    const double thresh = THRESH_SYMM*THRESH_SYMM;

    short z_symm = 0;
    short n_symm = 0;
    for (int n=0;n<natoms;n++)
    {
      float vz1 = tm1[2*n+0]; float vz2 = tm1[2*n+1]; // +/- z on atom n
      //printf("  orb: %2i  vz: %9.6f %9.6f \n",i1+1,vz1,vz2);
      double vz12 = vz1*vz2;
      if (fabs(vz12)>thresh)
      {
        n_symm++;
        if (vz12>0.)
          z_symm++;
        else
          z_symm--;
      }
    }
    if (prl>1) printf("  orbital %2i. parity? %i \n",i1+1,z_symm);

    if (n_symm>0)
    {
      if (z_symm==-n_symm)
        symm[i1] = 2;
      else
        symm[i1] = 1;
    }
    else
      symm[i1] = 99;
  }

  return;
}

void eval_c2v_ish_mos(int gsa, int iN, int N, float* norm, float* tm1, float** val1, double* jCA1, int prl, int* symm)
{
 //currently only useful for sigma bonds
 //can update to include π bonds
  for (int i1=0;i1<N;i1++)
  {
    if (prl>1) printf("\n orbital %2i \n",i1+1);

    generate_mo(i1,gsa,iN,N,norm,tm1,val1,jCA1);

    const double thresh = THRESH_SYMM;

    double v1 = tm1[0]; double v2 = tm1[3];
    double v1m = tm1[1]; double v2m = tm1[4];
    double v1p = tm1[2]; double v2p = tm1[5];

    int symm_z = 1; //value nonzero along z axis
    if (fabs(v1m)<thresh && fabs(v2m)<thresh && fabs(v1p)<thresh && fabs(v2p)<thresh)
      symm_z = 0;

    if (symm_z)
      symm[i1] = 1;
    else
      symm[i1] = 2;
  }

  return;
}

void eval_c2v_mos(int natoms, int gsa, int iN, int N, float* norm, float* tm1, float** val1, double* jCA1, int prl, int* symm)
{
  for (int i1=0;i1<N;i1++)
  {
    if (prl>1) printf("\n orbital %2i \n",i1+1);

    generate_mo(i1,gsa,iN,N,norm,tm1,val1,jCA1);

    const double thresh = THRESH_SYMM;

    short y_symm = 0;
    float vy1 = tm1[0]; float vy2 = tm1[1]; // +/- y on atom 1
    if (fabs(vy1)<thresh || fabs(vy2)<thresh)
    { 
      //printf("    using alt vy \n");
      vy1 = tm1[4]; vy2 = tm1[5];
    }
    if (fabs(vy1)>thresh)
    {
      if (vy1*vy2<0.) y_symm = -1; else y_symm = 1;
    }

    if (prl>1) printf("   vy: %8.5f %8.5f  y_symm: %2i \n",vy1,vy2,y_symm);

    short x_symm = 0;
    float vyz1 = tm1[2]; float vyz2 = tm1[3];
    if (fabs(vyz1)<thresh || fabs(vyz2)<thresh)
    {
      vyz1 = tm1[6]; vyz2 = tm1[7];
      //printf("    using alt vyz \n");
    }
    if (fabs(vyz1)>thresh)
    {
      if (vyz1*vyz2<0.) x_symm = -1; else x_symm = 1;
    }

    if (prl>1) printf("  vyz: %8.5f %8.5f  x_symm: %2i \n",vyz1,vyz2,x_symm);

    if (y_symm==-1)
    {
      if (x_symm==-1)
        symm[i1] = 2;
      else
        symm[i1] = 4;
    }
    else
    {
      if (x_symm==-1)
        symm[i1] = 3;
      else
        symm[i1] = 1;
    }
  }


  return;
}

void eval_d2h_ish_two_atom_mos(const int pgs, int gsa, int iN, int N, float* norm, float* tm1, float** val1, double* jCA1, int prl, int* symm)
{
  const double thresh = THRESH_SYMM/10.;
  int a1 = 0; int a2 = 1;

  for (int i1=0;i1<N;i1++)
  {
    if (prl>1) printf("\n orbital %2i \n",i1+1);

    generate_mo(i1,gsa,iN,N,norm,tm1,val1,jCA1);

    int zs = 0;
    if (fabs(tm1[a1*pgs+5])>thresh)
    {
      if (tm1[a1*pgs+5]*tm1[a2*pgs+4]<0.) zs = -1;
      else zs = 1;
    }
    else if (fabs(tm1[a1*pgs+4])>thresh)
    {
      if (tm1[a1*pgs+4]*tm1[a2*pgs+5]<0.) zs = -1;
      else zs = 1;
    }
    else
      zs = 0;

    if (prl>1) printf("  zs: %i \n",zs);

   //just do z symmetry, then whether val is zero along z direction
    if (zs==0)  symm[i1] = 3; //not cylindrical
    if (zs==1)  symm[i1] = 1; //reflects + over z
    if (zs==-1) symm[i1] = 2; //reflects - over z
  }

  return;
}

void eval_d2h_mos(const int pgs, int natoms, int gsa, int iN, int N, float* norm, float* tm1, float** val1, double* jCA1, int prl, int* symm)
{
  //prl+=2;

  for (int i1=0;i1<N;i1++)
  {
    if (prl>1) printf("\n orbital %2i \n",i1+1);

    generate_mo(i1,gsa,iN,N,norm,tm1,val1,jCA1);

    const double thresh = THRESH_SYMM;

    short xsall[natoms];
    short ysall[natoms];
    short zsall[natoms];
    for (int n=0;n<natoms;n++)
    {
      int ind1 = n*pgs;

      int xs = 0; int ys = 0; int zs = 0;
      if (fabs(tm1[ind1+0])>thresh)
      {
        if (tm1[ind1+0]*tm1[ind1+1]<0.) xs = -1;
        else xs = 1;
      }
      if (fabs(tm1[ind1+2])>thresh)
      {
        if (tm1[ind1+2]*tm1[ind1+3]<0.) ys = -1;
        else ys = 1;
      }
      if (fabs(tm1[ind1+4])>thresh)
      {
        if (tm1[ind1+4]*tm1[ind1+5]<0.) zs = -1;
        else zs = 1;
      }

      xsall[n] = xs; ysall[n] = ys; zsall[n] = zs;

      if (prl>1) printf("  atom %i  orbital %2i  xyz symm: %2i %2i %2i \n",n,i1,xs,ys,zs);
    }

   //reflection across xy axis
    short z_symm = 0;
    for (int n=0;n<natoms;n++)
      z_symm += zsall[n];
    if (prl>1) printf("   z reflect: %i \n",z_symm);

    //check for xz/yz reflections
    int x_symm = 0; int y_symm = 0;
    {
      int xz = 0; //reflection across xz plane
      float v1 = tm1[0*pgs+2]; //y direction of atom 1
      float v2;
      if (natoms==2)
        v2 = tm1[0*pgs+3]; //-y direction of atom 1
      else
        v2 = tm1[3*pgs+3]; //-y direction of atom 4
      if (fabs(v1)<thresh || fabs(v2)<thresh)
      {
        v1 = tm1[0*pgs+8];
        v2 = tm1[3*pgs+9];
      }
      if (fabs(v1)>thresh)
      {
        if (v1*v2<0.) xz = -1;
        else xz = 1;
      }
      y_symm = xz;

      int xz2 = 0;
      if (natoms==4)
      {
        float v3 = tm1[1*pgs+2];
        float v4 = tm1[2*pgs+3];
        if (fabs(v3)<thresh || fabs(v4)<thresh)
        {
          v3 = tm1[1*pgs+8];
          v4 = tm1[2*pgs+9];
        }
        if (fabs(v3)>thresh)
        {
          if (v3*v4<0.)
          xz2 = -1; else xz2 = 1;
        }

        if (xz==xz2) y_symm = xz;
        else y_symm = 0;
      }

      int yz = 0; //reflection across yz plane
      float v5 = tm1[0*pgs+1]; //-x direction of atom 1
      float v6 = tm1[1*pgs+0]; //x direction of atom 2
      if (fabs(v5)<thresh || fabs(v6)<thresh)
      {
        v5 = tm1[0*pgs+7];
        v6 = tm1[1*pgs+6];
      }
      if (fabs(v5)>thresh)
      {
        if (v5*v6<0.) yz = -1;
        else yz = 1;
      }
      x_symm = yz;

      int yz2 = 0;
      if (natoms==4)
      {
        float v7 = tm1[3*pgs+1];
        float v8 = tm1[2*pgs+0];
        if (fabs(v7)<thresh || fabs(v8)<thresh)
        {
          v7 = tm1[3*pgs+7];
          v8 = tm1[2*pgs+6];
        }
        if (fabs(v7)>thresh)
        { 
          if (v7*v8<0.) yz2 = -1;
          else yz2 = 1;
        }
        if (yz==yz2) x_symm = yz;
        else x_symm = 0;
      }

      if (prl>1) printf("  x reflect: %i %i  y reflect: %i %i \n",yz,yz2,xz,xz2);
    }

    if (z_symm>=0)
    {
     //ordered according to http://symmetry.jacobs-university.de/cgi-bin/group.cgi?group=602&option=4
      if (x_symm>=0 && y_symm>=0)
        symm[i1] = 1; //ag
      else if (x_symm==-1 && y_symm>=0)
        symm[i1] = 8; //b3u
      else if (x_symm>=0  && y_symm==-1)
        symm[i1] = 7; //b2u
      else if (x_symm==-1 && y_symm==-1)
        symm[i1] = 4; //b3g
    }
    else if (z_symm<0)
    {
      if (x_symm>=0 && y_symm>=0)
        symm[i1] = 6; //b1u
      else if (x_symm==-1 && y_symm>=0)
        symm[i1] = 3; //b2g
      else if (x_symm>=0  && y_symm==-1)
        symm[i1] = 2; //b1g
      else if (x_symm==-1 && y_symm==-1)
        symm[i1] = 5; //au
    }

  }
  return;
}

bool determine_mo_symmetry(const int point_group, int natoms, int* atno, double* coords, vector<vector<double> > basis, double* jCA, double* jS, int* symm, int prl)
{
  int N = basis.size();
  if (point_group==9)
  {
   //read from file
    return read_mo_symmetry(N,symm);
  }
  if (point_group<=1) return 0;

  printf("\n setting orbital symmetries for group: %i \n",point_group);
  printf("   available point groups: Cs (2), C2v (4), D2h (8) \n");

 //handle empty dummy atoms
  natoms = get_natoms_with_basis(natoms,atno,basis);

  bool two_atom_c2v = 0;
  bool two_atom_d2h = 0;
  int pgs = 6; //max grid pts per atom
  if (point_group==2)
  {
    bool in_z_plane = 1;
    for (int n=0;n<natoms;n++)
    {
      if (coords[3*n+2]!=0.)
        in_z_plane = 0;
    }
    if (!in_z_plane)
    {
      printf("  cannot use Cs symmetry because atoms are not in xy plane (z=0) \n");
      exit(-1);
    }
    pgs = 2;
  }
  else if (point_group==4)
  {
    if (natoms==1)
    {
      printf("\n cannot use C2v symmetry with only 1 atom \n");
      exit(-1);
    }
    else if (natoms==2)
    {
      two_atom_c2v = 1;
    }
    else if (natoms!=3)
    {
      printf("  WARNING: C2v little tested with >3 atoms \n");
      //printf("  C2v only for 2 or 3 atoms \n");
      //for (int j=0;j<N;j++) symm[j] = 0;
      //return 0;
    }
    pgs = 8;

    printf("   using C2v symmetry \n");

   //assuming atoms in this order, in xz plane (y=0)
   //      *1
   //
   //  *2      *3
   //
   // or for two atom case:
   //  *1        *2
  }
  else if (point_group==8)
  {
    if (natoms==2)
      two_atom_d2h = 1;
    else if (natoms==4)
    {
    }
    else if (natoms>4)
    {
      //printf("  checking D2h point group. Coords: %9.6f %9.6f %9.6f  %9.6f %9.6f %9.6f \n",coords[0],coords[1],coords[2],coords[3],coords[4],coords[5]);
      double dthresh = 1.e-6;

      bool two_symm = 0;
      if (coords[0]+coords[3]<dthresh && fabs(coords[0])>dthresh)
      {
        if (coords[1]==0. && coords[4]==0. && coords[2]==0. && coords[5]==0.)
          two_symm = 1;
      }
     //not clear this is usable
      //else if (coords[1]+coords[4]<dthresh && coords[1]!=0.)
      //{
      //  if (coords[0]==0. && coords[3]==0. && coords[2]==0. && coords[5]==0.)
      //    two_symm = 1;
      //}

      if (!two_symm)
      {
        printf("  D2h only for 2 or 4 atoms \n");
        for (int j=0;j<N;j++) symm[j] = 0;
        return 0;
      }

      two_atom_d2h = 1;
    }
    printf("   using D2h symmetry \n");
    pgs = 10;

   //assuming atoms in this order, all in xy plane
   //  *1        *2
   //
   //  *4        *3
   //
   // or for two atom case:
   //  *1        *2
  }
  else
  {
    printf("\n point group not supported \n");
    return 0;
  }


  int gsa = natoms*pgs;
  int gsa6 = gsa*6;

  //printf("  pgs: %2i  gsa: %3i \n",pgs,gsa);

  float* grid = new float[gsa6]();
  float* grid1 = new float[gsa6];

  int iN = N;
  float** val1 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gsa];

  int N2 = N*N;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  double* jCA1 = new double[N2];
  for (int j=0;j<N2;j++)
    jCA1[j] = jCA[j];

  float norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];

  float* tm1 = new float[gsa];

  #pragma acc enter data create(grid[0:gsa6],grid1[0:gsa6])
  #pragma acc enter data copyin(norm[0:N])
  #pragma acc enter data copyin(jCA1[0:N2])
  #pragma acc enter data create(val1[0:iN][0:gsa])
  #pragma acc enter data create(tm1[0:gsa])

  //create grid
  int wg = 0;
  if (point_group==2)
  {
    for (int n=0;n<natoms;n++)
    {
     //above and below the xy plane
      double A1 = coords[3*n+0]; double B1 = coords[3*n+1]; double C1 = 0.;
      grid[6*wg+0] = A1; grid[6*wg+1] = B1; grid[6*wg+2] = C1+1.; wg++;
      grid[6*wg+0] = A1; grid[6*wg+1] = B1; grid[6*wg+2] = C1-1.; wg++;
    }
  }
  else if (point_group==4 && !two_atom_c2v)
  {
    int a1 = 0; int a2 = 1; int a3 = 2;
    double A1 = coords[3*a1+0]; double B1 = coords[3*a1+1]; double C1 = coords[3*a1+2];
    double A2 = coords[3*a2+0]; double B2 = coords[3*a2+1]; double C2 = coords[3*a2+2];
    double A3 = coords[3*a3+0]; double B3 = coords[3*a3+1]; double C3 = coords[3*a3+2];
    //printf("  XYZ: \n %8.5f %8.5f %8.5f \n %8.5f %8.5f %8.5f \n %8.5f %8.5f %8.5f \n",A1,B1,C1,A2,B2,C2,A3,B3,C3);

   //to get y symm
    grid[6*wg+0] = A1; grid[6*wg+1] = B1+1.; grid[6*wg+2] = C1; wg++;
    grid[6*wg+0] = A1; grid[6*wg+1] = B1-1.; grid[6*wg+2] = C1; wg++;
   //atom 2-3 along x
    grid[6*wg+0] = A2-1.; grid[6*wg+1] = B2; grid[6*wg+2] = C2; wg++;
    grid[6*wg+0] = A3+1.; grid[6*wg+1] = B3; grid[6*wg+2] = C3; wg++;

   //alterative to get y symm
    grid[6*wg+0] = A2; grid[6*wg+1] = B2+1.; grid[6*wg+2] = C2; wg++;
    grid[6*wg+0] = A2; grid[6*wg+1] = B2-1.; grid[6*wg+2] = C2; wg++;
   //alternative to get x symm
    grid[6*wg+0] = A2; grid[6*wg+1] = B2+1.; grid[6*wg+2] = C2; wg++;
    grid[6*wg+0] = A3; grid[6*wg+1] = B3+1.; grid[6*wg+2] = C3; wg++;

   #if 0
   //alterative to get y symm
    grid[6*wg+0] = (A1+A2+A3)/3.; grid[6*wg+1] = (B1+B2+B3)/3.+1.; grid[6*wg+2] = (C1+C2+C3)/3.; wg++;
    grid[6*wg+0] = (A1+A2+A3)/3.; grid[6*wg+1] = (B1+B2+B3)/3.-1.; grid[6*wg+2] = (C1+C2+C3)/3.; wg++;
   //alternative to get x symm
    grid[6*wg+0] = (A1+A2+A3)/3.+1.; grid[6*wg+1] = (B1+B2+B3)/3.; grid[6*wg+2] = (C1+C2+C3)/3.; wg++;
    grid[6*wg+0] = (A1+A2+A3)/3.-1.; grid[6*wg+1] = (B1+B2+B3)/3.; grid[6*wg+2] = (C1+C2+C3)/3.; wg++;
   #endif
  }
  else if (point_group==4 && two_atom_c2v)
  {
    int a1 = 0; int a2 = 1;
    double A1 = coords[3*a1+0]; double B1 = coords[3*a1+1]; double C1 = coords[3*a1+2];
    double A2 = coords[3*a2+0]; double B2 = coords[3*a2+1]; double C2 = coords[3*a2+2];

   //points only along z axis
    grid[6*wg+0] = A1; grid[6*wg+1] = B1; grid[6*wg+2] = C1;    wg++;
    grid[6*wg+0] = A1; grid[6*wg+1] = B1; grid[6*wg+2] = C1-1.; wg++;
    grid[6*wg+0] = A1; grid[6*wg+1] = B1; grid[6*wg+2] = C1+1.; wg++;
    grid[6*wg+0] = A2; grid[6*wg+1] = B2; grid[6*wg+2] = C2;    wg++;
    grid[6*wg+0] = A2; grid[6*wg+1] = B2; grid[6*wg+2] = C2-1.; wg++;
    grid[6*wg+0] = A2; grid[6*wg+1] = B2; grid[6*wg+2] = C2+1.; wg++;
  }
  else if (point_group==8)
  {
    //1 Bohr above/below each atom
    for (int n=0;n<natoms;n++)
    {
      double A1 = coords[3*n+0]; double B1 = coords[3*n+1]; double C1 = coords[3*n+2];

      grid[6*wg+0] = A1+1.; grid[6*wg+1] = B1; grid[6*wg+2] = C1; wg++;
      grid[6*wg+0] = A1-1.; grid[6*wg+1] = B1; grid[6*wg+2] = C1; wg++;
      grid[6*wg+0] = A1; grid[6*wg+1] = B1+1.; grid[6*wg+2] = C1; wg++;
      grid[6*wg+0] = A1; grid[6*wg+1] = B1-1.; grid[6*wg+2] = C1; wg++;
      grid[6*wg+0] = A1; grid[6*wg+1] = B1; grid[6*wg+2] = C1+1.; wg++;
      grid[6*wg+0] = A1; grid[6*wg+1] = B1; grid[6*wg+2] = C1-1.; wg++;
     //+1 in z
      grid[6*wg+0] = A1+1.; grid[6*wg+1] = B1; grid[6*wg+2] = C1+1.; wg++;
      grid[6*wg+0] = A1-1.; grid[6*wg+1] = B1; grid[6*wg+2] = C1+1.; wg++;
      grid[6*wg+0] = A1; grid[6*wg+1] = B1+1.; grid[6*wg+2] = C1+1.; wg++;
      grid[6*wg+0] = A1; grid[6*wg+1] = B1-1.; grid[6*wg+2] = C1+1.; wg++;
    }
  }
  printf("   found %2i grid points \n",wg);

  #pragma acc update device(grid[0:gsa6])

  //evaluate AOs on grid
  for (int m=0;m<natoms;m++)
  {
   //working on this block of the matrix
    int s1 = 0; if (m>0) s1 = n2i[m-1]; int s2 = n2i[m];

    float Z1 = (float)atno[m];
    double A1 = coords[3*m+0]; double B1 = coords[3*m+1]; double C1 = coords[3*m+2];

    copy_grid(gsa,grid1,grid);
    recenter_grid_zero(gsa,grid1,-A1,-B1,-C1);

    #pragma acc parallel loop collapse(2) present(val1[0:iN][0:gsa])
    for (int i1=s1;i1<s2;i1++)
    for (int j=0;j<gsa;j++)
      val1[i1][j] = 1.f;

    for (int i1=s1;i1<s2;i1++)
    {
      int ii1 = i1;

      vector<double> basis1 = basis[i1];
      int n1 = basis1[0]; int l1 = basis1[1]; int m1 = basis1[2]; float zeta1 = basis1[3];

      eval_sh(ii1,gsa,grid1,val1[ii1],n1,l1,m1,zeta1);
    }

  } //loop m over natoms

  for (int i1=0;i1<N;i1++)
    symm[i1] = 0;

 //put MOs on grid and evaluate their symmetry
  if (point_group==2)
    eval_cs_mos(natoms,gsa,N,norm,tm1,val1,jCA1,prl,symm);
  else if (point_group==4 && !two_atom_c2v)
    eval_c2v_mos(natoms,gsa,iN,N,norm,tm1,val1,jCA1,prl,symm);
  else if (point_group==4 && two_atom_c2v)
   //rough c2v function
    eval_c2v_ish_mos(gsa,iN,N,norm,tm1,val1,jCA1,prl,symm);
  else if (point_group==8 && two_atom_d2h)
   //rough d2h functions
    eval_d2h_ish_two_atom_mos(pgs,gsa,iN,N,norm,tm1,val1,jCA1,prl,symm);
  else if (point_group==8)
   //should either fix or remove diatomics here:
    eval_d2h_mos(pgs,natoms,gsa,iN,N,norm,tm1,val1,jCA1,prl,symm);


  if (prl>0)
  {
    printf("\n symm:");
    for (int i1=0;i1<N;i1++)
      printf(" %i",symm[i1]);
    printf("\n");
  }

  for (int i1=0;i1<N;i1++)
  if (symm[i1]==0)
    printf(" WARNING: Missing symmetry assignment for orbital %i \n",i1+1);

 //clean up
  #pragma acc exit data delete(grid[0:gsa6],grid1[0:gsa6])
  #pragma acc exit data delete(val1[0:iN][0:gsa])
  #pragma acc exit data delete(norm[0:N],jCA1[0:N2])
  #pragma acc exit data delete(tm1[0:gsa])


  for (int i=0;i<iN;i++)
    delete [] val1[i];
  delete [] val1;

  delete [] grid;
  delete [] grid1;
  delete [] jCA1;
  delete [] tm1;
  delete [] n2i;

  return 1;
}


bool determine_mo_symmetry(const int point_group, int natoms, int* atno, float* coordsf, vector<vector<double> > basis, double* jCA, double* jS, int* symm, int prl)
{
  double coords[3*natoms];
  for (int j=0;j<3*natoms;j++)
    coords[j] = coordsf[j];

  bool result = determine_mo_symmetry(point_group,natoms,atno,coords,basis,jCA,jS,symm,prl);
  return result;
}

void determine_symmetry_atomic_sz(int N, int* symmblocks, vector<vector<double> >& basis, double* jCA, int prl)
{
 //symmetries only for partly symmetric atoms
 //z direction is priviledged

 //label "all others"
  for (int n=0;n<N;n++)
    symmblocks[n] = 42;

  double thresh = 1.e-4;
  for (int m=0;m<N;m++)
  {
    vector<double> basis1 = basis[m];

   //m==0 quantum numbers
    if (basis1[2]==0)
    for (int j=0;j<N;j++)
    if(fabs(jCA[m*N+j])>thresh)
      symmblocks[j] = 1;

    if (basis1[2]!=0)
    for (int j=0;j<N;j++)
    if(fabs(jCA[m*N+j])>thresh)
      symmblocks[j] = 2;

  }

  if (prl>0)
  {
    printf("  symm blocks: ");
    for (int n=0;n<N;n++)
      printf(" %i",symmblocks[n]);
    printf("\n");
  }

  return;
}

void determine_symmetry_atomic(int N, int* symmblocks, vector<vector<double> >& basis, double* jCA, int prl)
{
 //symmetries only for fully symmetric atoms

 //label "all others"
  for (int n=0;n<N;n++)
    symmblocks[n] = 42;

  double thresh = 1.e-4;
  for (int m=0;m<N;m++)
  {
    vector<double> basis1 = basis[m];

   //s symm ao
    if (basis1[1]==0)
    for (int j=0;j<N;j++)
    if(fabs(jCA[m*N+j])>thresh)
      symmblocks[j] = 1;

   //p symm ao (2-4)
    if (basis1[1]==1)
    for (int j=0;j<N;j++)
    if(fabs(jCA[m*N+j])>thresh)
    {
      symmblocks[j] = 3+basis1[2];
    }

   //d symm ao (5-9)
    if (basis1[1]==2)
    for (int j=0;j<N;j++)
    if(fabs(jCA[m*N+j])>thresh)
      symmblocks[j] = 7+basis1[2];

   //f symm ao (10-16)
    if (basis1[1]==3)
    for (int j=0;j<N;j++)
    if(fabs(jCA[m*N+j])>thresh)
      symmblocks[j] = 13+basis1[2];

   //g symm ao (17-25)
    if (basis1[1]==4)
    for (int j=0;j<N;j++)
    if(fabs(jCA[m*N+j])>thresh)
      symmblocks[j] = 21+basis1[2];

   //h symm ao (26-36)
    if (basis1[1]==5)
    for (int j=0;j<N;j++)
    if(fabs(jCA[m*N+j])>thresh)
      symmblocks[j] = 31+basis1[2];
  }

  if (prl>0)
  {
    printf("  symm blocks: ");
    for (int n=0;n<N;n++)
      printf(" %i",symmblocks[n]);
    printf("\n");
  }

  return;
}

void determine_symmetry_aos(int N, int* symb, double* jCA, int prl)
{
  double thresh = 1.e-4;
  int N2 = N*N;
  bool aoa[N2];
  for (int j=0;j<N2;j++) aoa[j] = 0;

  //printf("\n jCA: \n");
  //print_square(N,jCA);

  for (int i=0;i<N;i++)
  {
    for (int j=0;j<N;j++)
    if (fabs(jCA[j*N+i])>thresh)
      aoa[i*N+j] = 1;
  }

  if (prl>2)
  {
    printf("\n  aoa: \n");
    for (int i=0;i<N;i++)
    {
      printf("   ");
      for (int j=0;j<N;j++)
        printf(" %i",(int)aoa[i*N+j]);
      printf("\n");
    }
  }

  for (int i=0;i<N;i++) symb[i] = 0;

  int ncut[N];
  for (int i1=0;i1<N;i1++)
  {
    int ct1 = 0;
    for (int i2=0;i2<N;i2++)
      ct1 += aoa[i1*N+i2];

    if (ct1>1) ct1--;
    ncut[i1] = ct1;
  }

  int sn = 1; //symmetry block #
  symb[0] = sn;
  for (int i1=0;i1<N;i1++)
  {
    if (symb[i1]>0)
    for (int i2=i1+1;i2<N;i2++)
    {
      int nsame = 0;
      for (int j=0;j<N;j++)
      if (aoa[i1*N+j] && aoa[i2*N+j])
        nsame++;

      if (nsame>=ncut[i1])
        symb[i2] = symb[i1];
    }
    else
    {
      sn++;
      symb[i1] = sn;
      i1--;
    }
  }

  if (prl>0)
  {
    int Npr = min(N,24);
    if (prl>1) Npr = N;
    printf("  symm blocks: ");
    for (int i=0;i<Npr;i++)
      printf(" %i",symb[i]);
    printf("\n");
  }

  return;
}
