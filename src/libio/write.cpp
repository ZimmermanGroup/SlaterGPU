#include "write.h"
#include "qr.h"


#include <sstream>

#include "becke.h"

#define A2B 1.8897261


string get_iarray_name(short type1, short type2, short i1)
{
  string ft;
  if (type1==1)
    ft = "gao_";
  else if (type1==2)
    ft = "gpq_";
  else if (type1==3)
    ft = "gft_";
  else if (type1==4)
    ft = "rho_";
  else if (type1==5)
    ft = "tmp_";
  else
  {
    printf(" ERROR: type1 not understood \n");
    exit(1);
  }
  string fn;
  if (type2==1)
    fn = "S_";
  else if (type2==2)
    fn = "D_";
  else if (type2==3)
    fn = "T_";
  else if (type2==4)
    fn = "Q_";
  else if (type2==5)
    fn = "P_";
  else
  {
    printf(" ERROR: type2 not understood \n");
    exit(1);
  }

  string fc = to_string(i1);
  string filename = "scratch/"+ft+fn+fc;

  return filename;
}

void write_iarray(short type1, short type2, short i1, int s1, int s2, float* A)
{
  printf("\n shouldn't be here (for now) \n"); 

  string filename = get_iarray_name(type1,type2,i1);
  //printf("  writing %8s \n",filename.c_str());

  ofstream outfile;
  outfile.open(filename.c_str());
  for (int i=0;i<s1;i++)
  {
    for (int j=0;j<s2;j++)
      outfile << SSTRF(A[i*s2+j]) + " ";
    outfile << endl;
  }
  outfile.close();

  return;
}

void write_iarray(short type1, short type2, short i1, int s1, int s2, double* A)
{
  printf("\n shouldn't be here (for now) \n"); 

  string filename = get_iarray_name(type1,type2,i1);
  //printf("  writing %8s \n",filename.c_str());

  ofstream outfile;
  outfile.open(filename.c_str());
  for (int i=0;i<s1;i++)
  {
    for (int j=0;j<s2;j++)
      outfile << SSTRF2(A[i*s2+j]) + " ";
    outfile << endl;
  }
  outfile.close();

  return;
}

void save_geoms(int natoms, int* atno, vector<float> E, vector<float*> geoms, string fname)
{
  int ng = geoms.size();

  double B2A = 1./A2B;

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(8);

  for (int i=0;i<ng;i++)
  {
    float* coords = geoms[i];
    outfile << natoms << endl << E[i] << endl;
    for (int j=0;j<natoms;j++)
    {
      string aname = get_aname(atno[j]);
      double x1 = coords[3*j]*B2A; double y1 = coords[3*j+1]*B2A; double z1 = coords[3*j+2]*B2A;
      outfile << " " << aname << " " << x1 << " " << y1 << " " << z1 << endl;
    }
  }

  outfile.close();
  return;
}

void get_nxyzr(int n1, int l1, int m1, int& nx, int& ny, int& nz, int& nr)
{
  nx=ny=nz=nr=0;
  if (n1==1)
  {
  }
  else if (n1==2)
  {
    if (l1==0)
      nr = 1;
    else if (l1==1)
    {
      if (m1==0)      nz = 1;
      else if (m1==1) nx = 1;
      else            ny = 1;
    }
  }
  else if (n1==3)
  {
    if (l1==0)
      nr = 2;
    else if (l1==1)
    {
      nr = 1;
      if (m1==0)      nz = 1;
      else if (m1==1) nx = 1;
      else            ny = 1;
    }
    else if (l1==2)
    {
      nr = 0;
     #if CART_D
      if (m1==0) { nx = 2; }
      if (m1==1) { ny = 2; }
      if (m1==2) { nz = 2; }
      if (m1==3) { nx = 1; ny = 1; }
      if (m1==4) { nx = 1; nz = 1; }
      if (m1==5) { ny = 1; nz = 1; }
     #else
      if (m1==-2) { nx = 1; ny = 1; }
      if (m1==-1) { ny = 1; nz = 1; }
      if (m1== 0) { nz = 2; }
      if (m1== 1) { nx = 1; nz = 1; }
      if (m1== 2) { nx = 2; ny = 2; }
     #endif
    }
  }
  else if (n1==4)
  {
    //printf(" WARNING: n=4 MO print incomplete \n");
    if (l1==0)
      nr = 3;
    else if (l1==1)
    {
      nr = 2;
      if (m1==0)      nz = 1;
      else if (m1==1) nx = 1;
      else            ny = 1;
    }
    else if (l1==2)
    {
      nr = 1;
     #if CART_D
      if (m1==0) { nx = 2; }
      if (m1==1) { ny = 2; }
      if (m1==2) { nz = 2; }
      if (m1==3) { nx = 1; ny = 1; }
      if (m1==4) { nx = 1; nz = 1; }
      if (m1==5) { ny = 1; nz = 1; }
     #else
      if (m1==-2) { nx = 1; ny = 1; }
      if (m1==-1) { ny = 1; nz = 1; }
      if (m1== 0) { nz = 2; }
      if (m1== 1) { nx = 1; nz = 1; }
      if (m1== 2) { nx = 2; ny = 2; }
     #endif
    }
    else if (l1==3)
    {
     //INCOMPLETE
      if (m1==-2) { nx = 1; ny = 1; nz = 1; }
      else if (m1==0) { nz = 3; }
      else
        nx = nx = nz = 444;
    }

  }
  return;
}

void write_molden_g(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, double* eig, string fname)
{
  int N = basis.size();

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());

  double B2A = 1./A2B;

  outfile << fixed << setprecision(8);

  outfile << "[Molden Format]" << endl;
  outfile << "[Atoms] (Angs)" << endl;
  for (int i=0;i<natoms;i++)
  {
    string aname = get_aname(atno[i]);
    double x1 = coords[3*i]*B2A; double y1 = coords[3*i+1]*B2A; double z1 = coords[3*i+2]*B2A;
    outfile << " " << aname << "     " << i+1 << "   " << atno[i] << "   ";
    outfile << x1 << "   " << y1 << "   " << z1 << endl;
  }


  outfile << "[GTO]" << endl;

  for (int i=0;i<natoms;i++)
  {
    //double x1 = coords[3*i]; double y1 = coords[3*i+1]; double z1 = coords[3*i+2];

    string aname = get_aname(atno[i]);
    string b1 = read_basis_text(aname);
    outfile << i+1 << "    0 " << endl;
    outfile << b1 << endl;
      //outfile << nr << " " << zeta << " " << norm << " " << endl;
  }

  outfile << "[MO]" << endl;
  for (int i=0;i<N;i++)
  {
    int occ = 0; if (i<No) occ = 1;

    outfile << "Sym=X" << endl;
    outfile << "Ene=" << eig[i] << endl;
    outfile << "Spin=Alpha" << endl;
    outfile << "Occup=" << occ << endl;

    for (int j=0;j<N;j++)
      outfile << " " << j+1 << "  " << jCA[j*N+i] << endl;
  }

  outfile << " [5D7F] " << endl;


  outfile.close();

  return;
}

bool close_val(double v1, double v2)
{
  if (fabs(v1-v2)<1.e-4) return true;
  return false;
}

void write_molden_ss(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, string fname)
{
  int N = basis.size();

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());

  double B2A = 1./A2B;

  outfile << fixed << setprecision(8);

  outfile << "[Molden Format]" << endl;
  outfile << "[Atoms] (Angs)" << endl;
  for (int i=0;i<natoms;i++)
  {
    string aname = get_aname(atno[i]);
    double x1 = coords[3*i]*B2A; double y1 = coords[3*i+1]*B2A; double z1 = coords[3*i+2]*B2A;
    outfile << " " << aname << "     " << i+1 << "   " << atno[i] << "   ";
    outfile << x1 << "   " << y1 << "   " << z1 << endl;
  }


  outfile << "[STO]" << endl;

  int nss = 4;
  for (int i=0;i<natoms;i++)
  {
    double x1 = coords[3*i]; double y1 = coords[3*i+1]; double z1 = coords[3*i+2];

    for (int j=0;j<N;j++)
    {
      if (close_val(x1,basis[j][5]) && close_val(y1,basis[j][6]) && close_val(z1,basis[j][7]))
      {
        int n1 = basis[j][0]; int l1 = basis[j][1]; int m1 = basis[j][2]; double zeta = basis[j][3]*B2A; double norm = basis[j][4];
        int nx,ny,nz,nr;
        get_nxyzr(n1,l1,m1,nx,ny,nz,nr);

        if (l1<2)
        {
          for (int n=0;n<nss;n++)
          {
           //nr+n-l1?
            double norm1 = norm*pow(zeta,n); if (n==2) norm1 *= 0.4; if (n==3) norm1 /= 15.;
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr+n << " " << zeta << " " << norm1 << " " << endl;
          }
        }
        else if (l1==2)
        {
         //for 2 of the 5 3d functions, divide into parts
          if (nz==2)
          {
           //3dz2 --> 2.z2 - x2 - y2
            for (int n=0;n<nss;n++)
            {
              double norm1 = norm*pow(zeta,n); if (n==2) norm1 *= 0.4; if (n==3) norm1 /= 15.;
              nz = 2; nx = ny = 0;
              outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
              outfile << nr+n << " " << zeta << " " << 2.*norm1 << " " << endl;
              nz = 0; nx = 2; ny = 0;
              outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
              outfile << nr+n << " " << zeta << " " << norm1 << " " << endl;
              nz = 0; nx = 0; ny = 2;
              outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
              outfile << nr+n << " " << zeta << " " << norm1 << " " << endl;
            }
          }
          else if (nx==2 && ny==2)
          {
            //x2-y2 --> x2 - y2
            for (int n=0;n<nss;n++)
            {
              double norm1 = norm*pow(zeta,n); if (n==2) norm1 *= 0.4; if (n==3) norm1 /= 15.;
              nx = 2; ny = 0;
              outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
              outfile << nr+n << " " << zeta << " " << norm1 << " " << endl;
              nx = 0; ny = 2;
              outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
              outfile << nr+n << " " << zeta << " " << norm1 << " " << endl;
            }
          }
          else
          {
            for (int n=0;n<nss;n++)
            {
              double norm1 = norm*pow(zeta,n); if (n==2) norm1 *= 0.4; if (n==3) norm1 /= 15.;
              outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
              outfile << nr+n << " " << zeta << " " << norm1 << " " << endl;
            }
          }
        } //if l1==2

       #if 0
       //just drop most f or higher functions
        else if (l1==3)
        {
         //m==-2
          if (nx==1 && ny==1 && nz==1) //fxyz is the simplest 4f orbital
          {
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << norm << " " << endl;
          }
         //m==0
          if (nz==3)
          {
           //5z3
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << 5.*norm << " " << endl;
           //-3zr2
            nz = 1; nr = 2;
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << -3.*norm << " " << endl;
          }
          else
          {
            //incomplete
          }
        } //if l1==3
       #endif

      }
    }
  }

  outfile << "[MO]" << endl;
  for (int i=0;i<N;i++)
  {
    int occ = 0; if (i<No) occ = 1;

    outfile << "Sym=X" << endl;
    outfile << "Ene=0.0" << endl;
    outfile << "Spin=Alpha" << endl;
    outfile << "Occup=" << occ << endl;

    int wj = 0;
    for (int j=0;j<N;j++)
    {
      int n1 = basis[j][0]; int l1 = basis[j][1]; int m1 = basis[j][2];
      int nx,ny,nz,nr;
      get_nxyzr(n1,l1,m1,nx,ny,nz,nr);

      if (l1<2)
      {
        for (int n=0;n<nss;n++)
          outfile << " " << 1+wj++ << "  " << jCA[j*N+i] << endl;
      }
      else if (l1==2)
      {
        if (nz==2)
        {
          for (int n=0;n<nss;n++)
          {
            outfile << " " << 1+wj++ << "  " <<  jCA[j*N+i] << endl;
            outfile << " " << 1+wj++ << "  " << -jCA[j*N+i] << endl;
            outfile << " " << 1+wj++ << "  " << -jCA[j*N+i] << endl;
          }
        }
        else if (nx==2 && ny==2)
        {
          for (int n=0;n<nss;n++)
          {
            outfile << " " << 1+wj++ << "  " <<  jCA[j*N+i] << endl;
            outfile << " " << 1+wj++ << "  " << -jCA[j*N+i] << endl;
          }
        }
        else
        {
          for (int n=0;n<nss;n++)
            outfile << " " << 1+wj++ << "  " << jCA[j*N+i] << endl;
        }
      }
     #if 0
      else if (l1==3)
      {
       //incomplete
        if (nx==1 && ny==1 && nz==1) //m==-2
          outfile << " " << 1+wj++ << "  " << jCA[j*N+i] << endl;
        else if (nz==3) //m==0
          outfile << " " << 1+wj++ << "  " << jCA[j*N+i] << endl;
      }
     #endif
    }
  }

  outfile.close();

  return;
}

void write_molden(bool gbasis, int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, double* eig, string fname)
{
  if (gbasis) return write_molden_g(natoms,atno,coords,basis,jCA,No,eig,fname);
  if (basis[0].size()>10) return write_molden_ss(natoms,atno,coords,basis,jCA,No,fname);

  int N = basis.size();

  bool missing_ftns = 0;
  for (int j=0;j<N;j++)
  if (basis[j][0]>3)
    missing_ftns = 1;
  if (missing_ftns)
    printf(" WARNING: n=4 MO print incomplete \n");
  //printf("\n WARNING: no f function printing to molden \n");

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());

  double B2A = 1./A2B;

  outfile << fixed << setprecision(8);

  outfile << "[Molden Format]" << endl;
  outfile << "[Atoms] (Angs)" << endl;
  for (int i=0;i<natoms;i++)
  {
    string aname = get_aname(atno[i]);
    double x1 = coords[3*i]*B2A; double y1 = coords[3*i+1]*B2A; double z1 = coords[3*i+2]*B2A;
    outfile << " " << aname << "     " << i+1 << "   " << atno[i] << "   ";
    outfile << x1 << "   " << y1 << "   " << z1 << endl;
  }


  outfile << "[STO]" << endl;

  for (int i=0;i<natoms;i++)
  {
    double x1 = coords[3*i]; double y1 = coords[3*i+1]; double z1 = coords[3*i+2];

    for (int j=0;j<N;j++)
    {
      if (close_val(x1,basis[j][5]) && close_val(y1,basis[j][6]) && close_val(z1,basis[j][7]))
      {
        int n1 = basis[j][0]; int l1 = basis[j][1]; int m1 = basis[j][2]; double zeta = basis[j][3]*B2A; double norm = basis[j][4];
        int nx,ny,nz,nr;
        get_nxyzr(n1,l1,m1,nx,ny,nz,nr);

        if (l1<2)
        {
          outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
          outfile << nr << " " << zeta << " " << norm << " " << endl;
        }
        else if (l1==2)
        {
         //for 2 of the 5 3d functions, divide into parts
          if (nz==2)
          {
           //3dz2 --> 2.z2 - x2 - y2
            nz = 2; nx = ny = 0;
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << 2.*norm << " " << endl;
            nz = 0; nx = 2; ny = 0;
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << norm << " " << endl;
            nz = 0; nx = 0; ny = 2;
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << norm << " " << endl;
          }
          else if (nx==2 && ny==2)
          {
            //x2-y2 --> x2 - y2
            nx = 2; ny = 0;
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << norm << " " << endl;
            nx = 0; ny = 2;
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << norm << " " << endl;
          }
          else
          {
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << norm << " " << endl;
          }
        } //if l1==2

       #if 0
       //just drop most f or higher functions
        else if (l1==3)
        {
         //m==-2
          if (nx==1 && ny==1 && nz==1) //fxyz is the simplest 4f orbital
          {
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << norm << " " << endl;
          }
         //m==0
          if (nz==3)
          {
           //5z3
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << 5.*norm << " " << endl;
           //-3zr2
            nz = 1; nr = 2;
            outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
            outfile << nr << " " << zeta << " " << -3.*norm << " " << endl;
          }
          else
          {
            //incomplete
          }
        } //if l1==3
       #endif

      }
    }
  }

  outfile << "[MO]" << endl;
  for (int i=0;i<N;i++)
  {
    int occ = 0; if (i<No) occ = 1;

    outfile << "Sym=X" << endl;
    outfile << "Ene=" << eig[i] << endl;
    outfile << "Spin=Alpha" << endl;
    outfile << "Occup=" << occ << endl;

    int wj = 0;
    for (int j=0;j<N;j++)
    {
      int n1 = basis[j][0]; int l1 = basis[j][1]; int m1 = basis[j][2];
      int nx,ny,nz,nr;
      get_nxyzr(n1,l1,m1,nx,ny,nz,nr);

      if (l1<2)
        outfile << " " << 1+wj++ << "  " << jCA[j*N+i] << endl;
      else if (l1==2)
      {
        if (nz==2)
        {
          outfile << " " << 1+wj++ << "  " <<  jCA[j*N+i] << endl;
          outfile << " " << 1+wj++ << "  " << -jCA[j*N+i] << endl;
          outfile << " " << 1+wj++ << "  " << -jCA[j*N+i] << endl;
        }
        else if (nx==2 && ny==2)
        {
          outfile << " " << 1+wj++ << "  " <<  jCA[j*N+i] << endl;
          outfile << " " << 1+wj++ << "  " << -jCA[j*N+i] << endl;
        }
        else
          outfile << " " << 1+wj++ << "  " << jCA[j*N+i] << endl;
      }
     #if 0
      else if (l1==3)
      {
       //incomplete
        if (nx==1 && ny==1 && nz==1) //m==-2
          outfile << " " << 1+wj++ << "  " << jCA[j*N+i] << endl;
        else if (nz==3) //m==0
          outfile << " " << 1+wj++ << "  " << jCA[j*N+i] << endl;
      }
     #endif
    }
  }


  outfile.close();

  return;
}

void write_molden(bool gbasis, int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, string fname)
{
  int N = basis.size();
  double eig[N];
  for (int i=0;i<N;i++)
    eig[i] = 0.;
  return write_molden(gbasis,natoms,atno,coords,basis,jCA,No,eig,fname);
}

#if 0
void write_molden_vcf(int natoms, int* atno, double* coords, vector<vector<vcf> > &vcfs, string fname)
{
  for (int n=0;n<natoms;n++)
  for (int j=0;j<vcfs[n].size();j++)
  if (vcfs[n][j].numerical)
    return;

  int nvcft = 0;
  for (int n=0;n<natoms;n++)
    nvcft += vcfs[n].size();
  bool found_mu = 0;
  for (int n=0;n<natoms;n++)
  for (int j=0;j<vcfs[n].size();j++)
  if (vcfs[n][j].mu)
  { found_mu = 1; nvcft++; }
  if (found_mu) printf(" WARNING: mu lobes from p orbitals in molden output \n");

  int N = nvcft;
  int N2 = N*N;

  double* jCA = new double[N2]();
  double R12 = 1.;
  if (natoms==2) { double a12 = coords[0]-coords[3]; double b12 = coords[1]-coords[4]; double c12 = coords[2]-coords[5]; R12 = sqrt(a12*a12+b12*b12+c12*c12); }

  vector<vector<double> > vbasis;
  int s1 = 0;
  for (int n=0;n<natoms;n++)
  {
    float A1 = coords[3*n+0]; float B1 = coords[3*n+1]; float C1 = coords[3*n+2];

    vector<vcf> vcf1 = vcfs[n];
    int nvcf = vcf1.size();

    int jc = 0;
    for (int j=0;j<nvcf;j++)
    {
      float v1 = vcf1[j].v1;

      int n1 = vcf1[j].n; int l1 = vcf1[j].l; int m1 = vcf1[j].m;
      vector<double> b1;
      for (int m=0;m<10;m++) b1.push_back(0);
      b1[5] = A1; b1[6] = B1; b1[7] = C1;
      b1[8] = atno[n]; b1[9] = atno[n];

      b1[0] = n1; b1[1] = l1; b1[2] = m1;
      b1[3] = vcf1[j].zeta;
      b1[4] = vcf1[j].norm;

      vbasis.push_back(b1);

      if (vcf1[j].mu>0)
      {
       //ftns on 2 atoms
        vbasis.back()[3] = b1[3] = vcf1[j].zeta/R12;
        if (n==0) { b1[5] = coords[3]; b1[6] = coords[4]; b1[7] = coords[5]; }
        if (n==1) { b1[5] = coords[0]; b1[6] = coords[1]; b1[7] = coords[2]; }
        vbasis.push_back(b1);

        jCA[(s1+jc)*N+0] = v1;
        jCA[(s1+jc+1)*N+0] = -v1;

        jc++;
      }
      else if (vcf1[j].Vatom==0) //can't write Vatom
      {
        jCA[(s1+jc)*N+0] = v1;
       //to visualize individual components (all but first component)
        jCA[(s1+jc)*N+s1+jc] = v1;
      }
      jc++;

    } //loop j for atom n
    s1 += jc;
  }

  write_molden(0,natoms,atno,coords,vbasis,jCA,1,fname);

  delete [] jCA;

  return;
}
#endif

#if 0
void write_molden_vcf(int natoms, int* atno, float* coordsf, vector<vector<vcf> > &vcfs, string fname)
{
  double coords[3*natoms];
  for (int j=0;j<3*natoms;j++)
    coords[j] = coordsf[j];
  return write_molden_vcf(natoms,atno,coords,vcfs,fname);
}
#endif

void write_graph(int size, double* h, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  outfile << "<gexf xmlns=\"http://gexf.net/1.2\" version=\"1.2\">" << endl;

  outfile << "  <graph mode=\"static\" defaultedgetype=\"directed\">" << endl;

  outfile << "    <attributes class=\"node\">" << endl;
  outfile << "      <attribute id=\"0\" title=\"twos\" type=\"float\"/>" << endl;
  outfile << "    </attributes>" << endl;
  outfile << "    <nodes>" << endl;
  for (int i=0;i<size;i++)
  {
    double v1 = h[i*size+i];
    double v2 = -(v1-h[0])+2.;
    if (v2<0.) v2 = 0.;
    //outfile << "      <node id=\"" << i << "\" label=\"" << v1 << "\" twos=\"" << v2 << "\" />" << endl;
    outfile << "      <node id=\"" << i << "\" label=\"" << v1 << "\" >" << endl; 
    outfile << "        <attvalues>" << endl;
    outfile << "          <attvalue for=\"0\" value=\"" << v2 << "\" />" << endl;
    outfile << "        </attvalues>" << endl;
    outfile << "      </node>" << endl;
  }
  outfile << "    </nodes>" << endl;

  double threshw = 1.e-3;

  outfile << "    <edges>" << endl;
  int wi = 0;
  for (int i=0;i<size;i++)
  for (int j=0;j<size;j++)
  if (i!=j)
  {
    double wt = fabs(h[i*size+j]);
    if (wt>=threshw)
    {
      outfile << "      <edge id=\"" << wi << "\" source=\"" << i << "\" target=\"" << j << "\" weight=\"" << wt << "\" />" << endl;
      wi++;
    }
  }
  outfile << "    </edges>" << endl;

  outfile << "  </graph>" << endl;
  outfile << "</gexf>" << endl;

 #if 0
 //CSV file
  for (int i=0;i<size;i++)
  {
    outfile << ";" << i+1;
  }
  outfile << endl;

  for (int i=0;i<size;i++)
  {
    outfile << i+1 << ";";
    for (int j=0;j<size;j++)
    {
      if (i==j)
        outfile << "0";
      else
        outfile << fabs(h[i*size+j]);
      if (j!=size-1)
        outfile << ";";
    }
    outfile << endl;
  }
 #endif

  outfile.close();
  return;
}

#if 0
//use save_xyzv instead
void write_vecvec_csv(vector<vector<double> >& vec, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());
  outfile << fixed << setprecision(14);

  int s1 = vec.size();
  for (int j=0;j<s1;j++)
  {
    vector<double>& v1 = vec[j];

    int v1s = v1.size();
    for (int k=0;k<v1s-1;v1++)
      line += SSTRF2(v1[k]) + ",";
    line += SSTRF2(v1[v1s-1]);
    outfile << line << endl;
  }

  outfile.close();
  return;
}
#endif

void write_vector(int size, double* vals, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());
  outfile << fixed << setprecision(14);

  for (int j=0;j<size;j++)
  {
    string line = SSTRF2(vals[j]);

    outfile << line << " ";
  }
  outfile << endl;

  outfile.close();

  return;
}

void write_vector(int size, float* vals, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());
  outfile << fixed << setprecision(14);

  for (int j=0;j<size;j++)
  {
    string line = SSTRF2(vals[j]);

    outfile << line << " ";
  }
  outfile << endl;

  outfile.close();

  return;
}

void write_int(int val, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());

  string line = to_string(val);

  outfile << line << endl;

  outfile.close();

  return;
}

void write_float(double val, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());
  outfile << fixed << setprecision(10);

  string line = SSTRF(val);

  outfile << line << endl;

  outfile.close();

  return;
}

void write_double(double val, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());
  outfile << fixed << setprecision(10);

  string line = SSTRF2(val);

  outfile << line << endl;

  outfile.close();

  return;
}

void save_xyzv(int size1, int size2, double* vals, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(12);

  for (int j=0;j<size1;j++)
  {
    double* val1 = &vals[size2*j];
    for (int k=0;k<size2-1;k++)
      outfile << val1[k] << ",";
    outfile << val1[size2-1] << endl;
  }

  outfile.close();

  return;
}

void save_xyzv(vector<vector<double> >& vals, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(12);

  for (int j=0;j<vals.size();j++)
  {
    vector<double> val1 = vals[j];
    int nv = val1.size();

    for (int k=0;k<nv-1;k++)
      outfile << val1[k] << ",";
    outfile << val1[nv-1] << endl;
  }

  outfile.close();

  return;
}

void save_dft_exc(int save_type, double thresh, int natoms, int nrad, int nang, double* grid, double* rho, double* exc, string filename)
{
 //cannot yet use thresh/rho to select pts to print
 // instead, just print everything
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(12);

  if (save_type==1)
  {
    outfile << "DFT.  nrad: " << nrad << "  exc" << endl;

   //radial component only
    for (int j=0;j<nrad;j++)
    {
      int j1 = j*nang;
      //if (rho[j1]>thresh)
      {
        double v1 = exc[j1];
        string line = SSTRF2(v1);
        outfile << line << endl;
      }
    }
  }
  else if (save_type==2)
  {
    outfile << "DFT.  yz: " << 2*nrad << "   exc" << endl;

   //yz plane only
    int gs = nrad*nang;
    int gsa = natoms*gs;

    for (int j=0;j<gsa;j++)
    {
      if (fabs(grid[6*j+0])<1.e-12 && grid[6*j+1]>=0.)
      //if (fabs(grid[6*j+0])<1.e-12 && grid[6*j+1]>=0. && rho[j]>thresh)
      {
        //double z1 = grid[6*j+2];
        double v1 = exc[j];
        //string line = SSTRF(z1) + "," + SSTRF2(v1);
        string line = SSTRF2(v1);
        outfile << line << endl;
      }
    }
  }
  else if (save_type==3)
  {
    outfile << "DFT.  z: " << 2*nrad << "   exc" << endl;

   //z axis only
    int gs = nrad*nang;
    int gsa = natoms*gs;

    for (int j=0;j<gsa;j++)
    {
      if (fabs(grid[6*j+0])<1.e-12 && fabs(grid[6*j+1])<1.e-12)
      //if (fabs(grid[6*j+0])<1.e-12 && fabs(grid[6*j+1])<1.e-12 && rho[j]>thresh)
      {
        //double z1 = grid[6*j+2];
        double v1 = exc[j];
        //string line = SSTRF(z1) + "," + SSTRF2(v1);
        string line = SSTRF2(v1);
        outfile << line << endl;
      }
    }
  }

  outfile.close();
  return;
}

void save_dft_exc(bool save_type, int natoms, int nrad, int nang, double* grid, double* exc, string filename)
{
  return save_dft_exc(save_type,natoms,nrad,nang,grid,exc,filename);
}

void get_atom_neighbors(int natoms, int* atno, double* coords, int* atneighbors, double* atnqr)
{
  for (int n=0;n<natoms;n++)
  {
    double xa1 = coords[3*n+0]; double ya1 = coords[3*n+1]; double za1 = coords[3*n+2];

    double d12min = 1000.;
    int i1 = 0;
    for (int m=0;m<natoms;m++)
    if (n!=m)
    {
      double xa2 = coords[3*m+0]; double ya2 = coords[3*m+1]; double za2 = coords[3*m+2];
      double x12 = xa2-xa1; double y12 = ya2-ya1; double z12 = za2-za1;
      double d12 = x12*x12 + y12*y12 + z12*z12;

      if (d12<d12min) { d12min = d12; i1 = m; }
    }
    atneighbors[n] = i1;
    atnqr[n] = QR[atno[i1]-1];
  }

  return;
}

void dist_to_atom_pair(int natoms, int* atno, double* coords, const double x1, const double y1, const double z1, double& ra1, double& ra2)
{
  for (int n=0;n<natoms;n++) if (atno[n]>18) { printf("\n ERROR: quantum radii not defined for Z>18 \n"); exit(-1); }

  int nmin = 0;
  double r2min = 10000.;
  for (int n=0;n<natoms;n++)
  {
    //double qrn = QR[atno[n]];
    double qrn = 1./atno[n]; //r -> Zr

    double x10 = x1-coords[3*n+0];
    double y10 = y1-coords[3*n+1];
    double z10 = z1-coords[3*n+2];

    double r2 = x10*x10+y10*y10+z10*z10;
    r2 /= qrn*qrn;
    if (r2<r2min)
    {
      r2min = r2;
      nmin = n;
    }
  }
  ra1 = sqrt(r2min);

  r2min = 10000.;
  for (int n=0;n<natoms;n++)
  if (n!=nmin)
  {
    //double qrn = QR[atno[n]];
    double qrn = 1./atno[n]; //r -> Zr

    double x10 = x1-coords[3*n+0];
    double y10 = y1-coords[3*n+1];
    double z10 = z1-coords[3*n+2];

    double r2 = x10*x10+y10*y10+z10*z10;
    r2 /= qrn*qrn;
    if (r2<r2min)
      r2min = r2;
  }
  ra2 = sqrt(r2min);
  if (natoms==1) ra2 = ra1;

  return;
}

double dist_to_hz(int natoms, int* atno, double* coords, const double x1, const double y1, const double z1)
{
 //returns distance to closest H atom
 // or -distance to closest non-H atom
 // distances are scaled by approximate quantum radii

  for (int n=0;n<natoms;n++) if (atno[n]>18) { printf("\n ERROR: quantum radii not defined for Z>18 \n"); exit(-1); }

  int* atneighbors = new int[natoms];
  double* atnqr = new double[natoms];
  get_atom_neighbors(natoms,atno,coords,atneighbors,atnqr);

 //dist to H rescaled by radius of its neighboring atom
  double rh2min = 10000.;
  for (int n=0;n<natoms;n++)
  if (atno[n]==1)
  {
    double qrn = atnqr[n];

    double x10 = x1-coords[3*n+0];
    double y10 = y1-coords[3*n+1];
    double z10 = z1-coords[3*n+2];

    double r2 = x10*x10+y10*y10+z10*z10;
    r2 /= qrn*qrn;
    if (r2<rh2min) rh2min = r2;
  }

  double rz2min = 10000.;
  for (int n=0;n<natoms;n++)
  if (atno[n]>1)
  {
    int Z = atno[n];
    double qr1 = QR[Z-1];

    double x10 = x1-coords[3*n+0];
    double y10 = y1-coords[3*n+1];
    double z10 = z1-coords[3*n+2];

    double r2 = x10*x10+y10*y10+z10*z10;
    r2 /= qr1*qr1;
    if (r2<rz2min) rz2min = r2;
  }

  double rh = sqrt(rh2min);
  double rz = sqrt(rz2min);

  double r = -rz;
  if (rh<fabs(rz)) r = rh;

  delete [] atneighbors;
  delete [] atnqr;

  return r;
}

vector<vector<double> > process_rdtve(vector<vector<double> >& rdtve1, int type, const double rthresh, const double dthresh)
{
  int offset = 0;
  if (type==2) offset = 6; //if this changes also change the sort command below

  vector<vector<double> > rdtve2;
  if (type==1)
    sort(rdtve1.begin(),rdtve1.end());
  else if (type==2)
    sort(rdtve1.begin(),rdtve1.end(),
        [](const vector<double>& a, const vector<double>& b)
         { return a[6] < b[6]; });
  else
  {
    printf("\n ERROR: cannot sort due to incorrect type \n");
  }

 //sorted list with rho>rthresh
  vector<vector<double> > tmpr;
  for (int j=0;j<rdtve1.size();j++)
  if (rdtve1[j][offset+0]>rthresh)
    tmpr.push_back(rdtve1[j]);

 //delete similar cases
  double den = 1.e-8;
  for (int j=0;j<tmpr.size()-1;j++)
  {
    double rho1 = tmpr[j][offset+0];
    double drho1 = tmpr[j][offset+1];
    double Td1 = tmpr[j][offset+2];

    bool keep = 0;
    for (int k=j+1;k<tmpr.size();k++)
    {
      double rho2 = tmpr[k][offset+0];
      double rf = (rho2-rho1)/rho1;

      if (rf>dthresh) //unique
      {
        keep = 1;
        break;
      }

      double drho2 = tmpr[k][offset+1];
      double Td2 = tmpr[k][offset+2];
      double df = fabs(drho1-drho2)/(drho1+den);
      double tf = fabs(Td1-Td2)/(Td1+den);

      if (df<dthresh && tf<dthresh) //duplicate
      {
        //printf("  rdt:  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f   rdtf:  %6.1e  %6.1e  %6.1e  \n",rho1,drho1,Td1,rho2,drho2,Td2,rf,df,tf);
        break;
      }
    } //loop k over partial upper triangle

    if (keep)
    {
      //printf("  rdt(%i):  %9.6f  %9.6f  %9.6f \n",1*keep,rho1,drho1,Td1);
      rdtve2.push_back(tmpr[j]);
    }
  } //loop j over tmpr

  return rdtve2;
}

void save_dft_vals_rh(double thresh, int natoms, int* atno, double* coords, int nrad, int nang, double* grid, double* rho, double* drho, double* Td, string filename)
{
  int gs = nrad*nang;
  int gsa = gs*natoms;

  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(8);

  double rthresh = 1.e-8;

  vector<vector<double> > rdtve;
  for (int j=0;j<gsa;j++)
  if (rho[j]>rthresh)
  {
    vector<double> rdt1;
    // = { grid[6*j+0], grid[6*j+1], grid[6*j+2] };
    rdt1.push_back(rho[j]);
    rdt1.push_back(drho[j]);
    rdt1.push_back(Td[j]);
    rdt1.push_back(0.);
    rdt1.push_back(0.);
    rdt1.push_back(0.);
    rdt1.push_back(grid[6*j+0]);
    rdt1.push_back(grid[6*j+1]);
    rdt1.push_back(grid[6*j+2]);

    rdtve.push_back(rdt1);
  }

  rdtve = process_rdtve(rdtve,1,rthresh,thresh);

  int npts = rdtve.size();
  outfile << "DFT.  rho/|drho|/Td/rh" << endl;

  for (int j=0;j<npts;j++)
  {
    double x1 = rdtve[j][6]; double y1 = rdtve[j][7]; double z1 = rdtve[j][8];
    double rh = dist_to_hz(natoms,atno,coords,x1,y1,z1);

    string line = SSTRF2(rdtve[j][0]) + "," + SSTRF2(rdtve[j][1]) + "," + SSTRF2(rdtve[j][2]) + "," + SSTRF2(rh);
    outfile << line << endl;
  }

  outfile.close();

  return;
}


void save_dft_vals(int save_type, double thresh, int natoms, int* atno, double* coords, int nrad, int nang, double* grid, double* wt,
                   double* rho, double* drho, double* Td, double* lapl, double* hessw, double* hessp,
                   double* epsi, double* vc, double* ec, int zpos, string filename)
{
  int gs = nrad*nang;
  int gsa = gs*natoms;

  if (natoms==1) printf("\n WARNING: save_dft_vals: atomic r is not rescaled \n");

 //just print vc twice if ec is not available
  if (ec==NULL) ec = vc;

  if (epsi==NULL)
  { printf(" ERROR: epsi now required in save_dft_vals \n"); exit(-1); }

 //remove pts with small integration weight
  const double wthresh = 1.e-6; //low weights away from nuclei -> Becke partitioned
  int gscut = gs/4; //outside of core region

  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(12);

  if (save_type==1)
  {
    outfile << "DFT.  nrad: " << nrad << "  grid/rho/|drho|/Td/l/w/hp,eps/vc/ec" << endl;

   //approximating density at nucleus by the nearest radial point
    double atden = rho[0]*pow(atno[0],-3.);

   //radial component only (single atom)
    for (int j=0;j<nrad;j++)
    {
      int j1 = j*nang;

      double x1 = grid[6*j1+0]; double y1 = grid[6*j1+1]; double z1 = grid[6*j1+2];
      double rh = dist_to_hz(natoms,atno,coords,x1,y1,z1);
      double ra1,ra2; dist_to_atom_pair(natoms,atno,coords,x1,y1,z1,ra1,ra2);

      {
        string line = SSTRF2(rh) + "," + SSTRF2(ra1) + "," + SSTRF2(ra2) + ",0.0," + SSTRF2(grid[6*j1+1]) + "," + SSTRF2(grid[6*j1+2])
              + "," + SSTRF2(rho[j1]) + "," + SSTRF2(drho[j1])
              + "," + SSTRF2(Td[j1]) + "," + SSTRF2(lapl[j1]) + "," + SSTRF2(hessw[j1])
              + "," + SSTRF2(hessp[3*j1+0]) + "," + SSTRF2(hessp[3*j1+1]) + "," + SSTRF2(hessp[3*j1+2])
              + "," + SSTRF2(atden)
              + "," + SSTRF2(epsi[j1]) + "," + SSTRF2(vc[j1]) + "," + SSTRF2(ec[j1]);
      /*  string line = SSTRF2(grid[6*j1+3]) + ",0.0,0.0," + SSTRF2(grid[6*j1+3]) + "," + SSTRF2(rho[j1]) + "," + SSTRF2(drho[j1]) + "," + SSTRF2(Td[j1])
         + "," + SSTRF2(lapl[j1]) + "," + SSTRF2(hessw[j1])
         + "," + SSTRF2(hessp[3*j+0]) + "," + SSTRF2(hessp[3*j+1]) + "," + SSTRF2(hessp[3*j+2])
         + "," + SSTRF2(epsi[j1]) + "," + SSTRF2(vc[j1]) + "," + SSTRF2(ec[j1]); */
        outfile << line << endl;
      }
    }
  }
  else if (save_type==2)
  {
    outfile << "DFT. rh/r1/r2/xyz: " << 2*nrad << " grid/rho/|drho|/Td/l/w/hp/eps/vc/ec" << endl;
   //within yz plane (of diatomic)
    for (int n=0;n<natoms;n++)
    {
      double atden = rho[n*gs]*pow(atno[n],-3.); //approximating density at nucleus by the nearest radial point

      for (int k=0;k<gs;k++)
      {
        int j = n*gs+k;

        if (k<gscut || wt[j]>wthresh)
        if (fabs(grid[6*j+0])<1.e-12 && grid[6*j+1]>=0.)
        {
          double x1 = grid[6*j+0]; double y1 = grid[6*j+1]; double z1 = grid[6*j+2];
          double rh = dist_to_hz(natoms,atno,coords,x1,y1,z1);
          double ra1,ra2; dist_to_atom_pair(natoms,atno,coords,x1,y1,z1,ra1,ra2);

          string line = SSTRF2(rh) + "," + SSTRF2(ra1) + "," + SSTRF2(ra2) + ",0.0," + SSTRF2(grid[6*j+1]) + "," + SSTRF2(grid[6*j+2])
                + "," + SSTRF2(rho[j]) + "," + SSTRF2(drho[j])
                + "," + SSTRF2(Td[j]) + "," + SSTRF2(lapl[j]) + "," + SSTRF2(hessw[j])
                + "," + SSTRF2(hessp[3*j+0]) + "," + SSTRF2(hessp[3*j+1]) + "," + SSTRF2(hessp[3*j+2])
                + "," + SSTRF2(atden)
                + "," + SSTRF2(epsi[j]) + "," + SSTRF2(vc[j]) + "," + SSTRF2(ec[j]);
          outfile << line << endl;
        }
      } //loop k over atomic grid
    } //loop n over natoms
  }
  else if (save_type==3)
  {
    vector<vector<double> > rdtve;
    for (int n=0;n<natoms;n++)
    {
      double atden = rho[n*gs]*pow(atno[n],-3.); //approximating density at nucleus by the nearest radial point

      for (int k=0;k<gs;k++)
      {
        int j = n*gs+k;
        if (k<gscut || wt[j]>wthresh)
        {
          vector<double> rdt1;
          rdt1.push_back(rho[j]);
          rdt1.push_back(drho[j]);
          rdt1.push_back(Td[j]);
          rdt1.push_back(lapl[j]);
          rdt1.push_back(hessw[j]); //5th item
          for (int k=0;k<3;k++)
            rdt1.push_back(hessp[3*j+k]);
          rdt1.push_back(atden);
          rdt1.push_back(epsi[j]); //10th item
          rdt1.push_back(vc[j]);
          rdt1.push_back(ec[j]);
          rdt1.push_back(grid[6*j+0]); //13th item (index 12)
          rdt1.push_back(grid[6*j+1]);
          rdt1.push_back(grid[6*j+2]);

          rdtve.push_back(rdt1);
        }

      } //loop over atomic grid
    } //loop over natoms

    double rthresh = 1.e-8;
    rdtve = process_rdtve(rdtve,1,rthresh,thresh);

    int npts = rdtve.size();
    int ndata = rdtve[0].size();
    int offset = ndata-3;
    outfile << "DFT.  rh/r1/r2/xyz: " << npts << " grid/rho/|drho|/Td/l/w/hp,eps/vc/ec" << endl;

    for (int j=0;j<npts;j++)
    {
      string line = "";

      double x1 = rdtve[j][offset]; double y1 = rdtve[j][offset+1]; double z1 = rdtve[j][offset+2];
      double rh = dist_to_hz(natoms,atno,coords,x1,y1,z1);
      double ra1,ra2; dist_to_atom_pair(natoms,atno,coords,x1,y1,z1,ra1,ra2);
      line += SSTRF2(rh) + "," + SSTRF2(ra1) + "," + SSTRF2(ra2) + ",";

      for (int k=offset;k<ndata;k++) //XYZ
        line += SSTRF2(rdtve[j][k]) + ",";
      for (int k=0;k<offset-1;k++) //rho,drho,tau,lapl,hessw,atden,eps,vc
        line += SSTRF2(rdtve[j][k]) + ",";
      line += SSTRF2(rdtve[j][offset-1]); //ec

      outfile << line << endl;
    }
  }
  else if (save_type==4)
  {
    printf("\n WARNING: type=4 save_dft_vals not up to date \n\n");
    outfile << "DFT.  z: " << 2*nrad*natoms << " grid/rho/|drho|/Td/vc" << endl;
   //along z axis
    for (int j=0;j<gsa;j++)
    {
      if (fabs(grid[6*j+0])<1.e-12 && fabs(grid[6*j+1])<1.e-12)
      {
        string line = SSTRF2(grid[6*j+2]) + "," + SSTRF2(rho[j]) + "," + SSTRF2(drho[j]) + "," + SSTRF2(Td[j])
           + "," + SSTRF2(lapl[j]) + "," + SSTRF2(hessw[j]) + "," + SSTRF2(epsi[j]) + "," + SSTRF2(vc[j]) + "," + SSTRF2(ec[j]);
        outfile << line << endl;
      }
    }
  }
  else
  {
    printf(" ERROR: wtype %i not supported in save_dft_vals \n",save_type);
  }

  outfile.close();
  return;
}

void save_dft_vals(bool save_type, double thresh, int natoms, int nrad, int nang, double* grid, double* rho, double* drho, double* Td, double* epsi, double* vc, int zpos, string filename)
{
  return save_dft_vals(save_type,thresh,natoms,nrad,nang,grid,rho,drho,Td,epsi,vc,zpos,filename);
}

void write_grid(int natoms, int nrad, int nang, double* grid, double* wt)
{
  int gs = nrad*nang;
  string filename = "GRID_WTS";
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(10);

  outfile << filename << ":" << endl;
  for (int n=0;n<natoms;n++)
  {
    int i0 = n*gs;
    for (int i=0;i<nrad;i++)
    {
      for (int j=0;j<nang;j++)
      {
        int i1 = i0 + i*nang+j;
        string line = "r: " + SSTRF(grid[6*i1+3]) + " xyz: ";
        line += " " + SSTRF(grid[6*i1+0]) + " " + SSTRF(grid[6*i1+1]) + " " + SSTRF(grid[6*i1+2]);
        line += " wt: " + SSTRF(wt[i1]);
        outfile << line << endl;
      }
    }
  }

  outfile.close();
  return;
}

void write_gridpts(int s1, int s2, float* A, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(10);

  outfile << filename << ":" << endl;
  for (int i=0;i<s1;i++)
  {
    string line = "";
    for (int j=0;j<s2;j++)
      line += " " + SSTRF(A[i*s2+j]);
    outfile << line << endl;
  }

  outfile.close();
  return;
}

void write_gridpts(int s1, int s2, double* A, string filename)
{
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(15);

  outfile << filename << ":" << endl;
  for (int i=0;i<s1;i++)
  {
    string line = "";
    for (int j=0;j<s2;j++)
      line += " " + SSTRF2(A[i*s2+j]);
    outfile << line << endl;
  }

  outfile.close();
  return;
}

void write_square_clean(int N, double* A, string fname, double thresh, int prl)
{
  if (prl>1) printf("  writing %s to file \n",fname.c_str());

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(14);

  outfile << fname << ":" << endl;
  for (int i=0;i<N;i++)
  {
    string line = "";
    for (int j=0;j<N;j++)
    if (fabs(A[i*N+j])>thresh)
      line += " " + SSTRF2(A[i*N+j]);
    else
      line += " 0.";
    outfile << line << endl;
  }

  outfile.close();
  return;
}

void write_square(int N, int M, double** A, string fname, int prl)
{
  if (prl>1) printf("  writing %s to file \n",fname.c_str());

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(14);

  outfile << fname << ":" << endl;
  for (int i=0;i<N;i++)
  {
    string line = "";
    for (int j=0;j<M;j++)
      line += " " + SSTRF2(A[i][j]);
    outfile << line << endl;
  }

  outfile.close();
  return;
}

void write_square(int N, double* A, string fname, int prl)
{
  if (prl>1) printf("  writing %s to file \n",fname.c_str());

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(14);

  outfile << fname << ":" << endl;
  for (int i=0;i<N;i++)
  {
    string line = "";
    for (int j=0;j<N;j++)
      line += " " + SSTRF2(A[i*N+j]);
    outfile << line << endl;
  }

  outfile.close();
  return;
}

void write_square(int N, float* A, string fname, int prl)
{
  if (prl>1) printf(" writing %s to file \n",fname.c_str());

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(10);

  outfile << fname << ":" << endl;
  for (int i=0;i<N;i++)
  {
    string line = "";
    for (int j=0;j<N;j++)
      line += " " + SSTRF(A[i*N+j]);
    outfile << line << endl;
  }

  outfile.close();
  return;
}

void write_rect(int N2, int Naux, double* C, string filename)
{
  printf("  writing %s to file \n",filename.c_str());

  ofstream outfile;
  outfile.open(filename.c_str());

  //const double thresh = 1.e-15;
  outfile << fixed << setprecision(16);

  outfile << filename << endl;
  for (int i=0;i<N2;i++)
  {
    string line = "";
    for (int j=0;j<Naux;j++)
    //if (fabs(C[i*Naux+j])>thresh)
      line += " " + SSTRF2(C[i*Naux+j]);
    //else
    //  line += " 0.0";
    outfile << line << endl;
  }

  outfile.close();
  return;
}

void write_C(int Naux, int N2, float* C)
{
  printf("  writing C to file \n");

  string filename = "Ciap";
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(10);

  outfile << "Ciap:" << endl;
#if 1
  for (int i=0;i<N2;i++)
  {
    string line = "";
    for (int j=0;j<Naux;j++)
      line += " " + SSTRF(C[i*Naux+j]);
    outfile << line << endl;
  }
#else
  int nna = Naux*sqrt(N2);
  int na = Naux;
  for (int p=0;p<Naux;p++)
  {
    string line = "";
    for (int m=0;m<N;m++)
    for (int n=0;n<N;n++)
      line += " " + SSTRF(C[m*nna+n*na+p]);
    outfile << line << endl;
  }
#endif

  outfile.close();
  return;
}

void write_C(int Naux, int N2, double* C)
{
  printf("  writing C to file \n");

  string filename = "Ciap";
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(16);

  outfile << "Ciap:" << endl;
#if 1
  for (int i=0;i<N2;i++)
  {
    string line = "";
    for (int j=0;j<Naux;j++)
      line += " " + SSTRF2(C[i*Naux+j]);
    outfile << line << endl;
  }
#else
  int nna = Naux*sqrt(N2);
  int na = Naux;
  for (int p=0;p<Naux;p++)
  {
    string line = "";
    for (int m=0;m<N;m++)
    for (int n=0;n<N;n++)
      line += " " + SSTRF2(C[m*nna+n*na+p]);
    outfile << line << endl;
  }
#endif

  outfile.close();
  return;
}

void write_Cy(int Naux, int N2, double* C)
{
  printf("  writing Cy to file \n");

  string filename = "Cyiap";
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(16);

  outfile << "Cyiap:" << endl;
  for (int i=0;i<N2;i++)
  {
    string line = "";
    for (int j=0;j<Naux;j++)
      line += " " + SSTRF2(C[i*Naux+j]);
    outfile << line << endl;
  }

  outfile.close();
  return;
}

void write_Col(int Naux, int N2, double* C)
{
  printf("  writing Col to file \n");

  string filename = "Coiap";
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(16);

  outfile << "Coiap:" << endl;
  for (int i=0;i<N2;i++)
  {
    string line = "";
    for (int j=0;j<Naux;j++)
      line += " " + SSTRF2(C[i*Naux+j]);
    outfile << line << endl;
  }

  outfile.close();
  return;
}

void write_S_En_T(int N, float* S, float* En, float* T)
{
  printf("  writing S/En/T to file \n");

  string filename = "SENT";
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(10);

  outfile << "S:" << endl;
  for (int i=0;i<N;i++)
  {
    string line = "";
    for (int j=0;j<N;j++)
      line += " " + SSTRF(S[i*N+j]);
    outfile << line << endl;
  }

  outfile << "En:" << endl;
  for (int i=0;i<N;i++)
  {
    string line = "";
    for (int j=0;j<N;j++)
      line += " " + SSTRF(En[i*N+j]);
    outfile << line << endl;
  }

  outfile << "T:" << endl;
  for (int i=0;i<N;i++)
  {
    string line = "";
    for (int j=0;j<N;j++)
      line += " " + SSTRF(T[i*N+j]);
    outfile << line << endl;
  }

  outfile.close();
  return;
}


void write_S_En_T(int N, double* S, double* En, double* T)
{
  printf("  writing S/En/T to file \n");

  string filename = "SENT";
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(14);

  outfile << "S:" << endl;
  for (int i=0;i<N;i++)
  {
    string line = "";
    for (int j=0;j<N;j++)
      line += " " + SSTRF2(S[i*N+j]);
    outfile << line << endl;
  }

  outfile << "En:" << endl;
  for (int i=0;i<N;i++)
  {
    string line = "";
    for (int j=0;j<N;j++)
      line += " " + SSTRF2(En[i*N+j]);
    outfile << line << endl;
  }

  outfile << "T:" << endl;
  for (int i=0;i<N;i++)
  {
    string line = "";
    for (int j=0;j<N;j++)
      line += " " + SSTRF2(T[i*N+j]);
    outfile << line << endl;
  }

  outfile.close();
  return;
}

void write_mo_grid(int natoms, int* atno, double* coords, int nrad, int gsa, vector<vector<double> > basis,
                   int No, double* Pao, double* jCA, double* jS, float* grid, float* wt, int prl)
{
  if (nrad<1) nrad = 1;
  if (natoms>3) return;

 //assumes jCA on gpu

  //int gs = gsa/natoms;
  //int gsa3 = 3*gsa;
  int gsa6 = 6*gsa;

  printf("\n getting MOs on the grid (size: %6i) \n",gsa);
  if (gsa<1) { exit(-1); }

  int nang = gsa/nrad/natoms;
  
  bool sgs_basis = 0;
  bool ss_basis = 0;
  if (basis[0].size()>10) ss_basis = 1;
  //if (basis[0].size()>10) sgs_basis = 1;
  if (sgs_basis && prl>1) printf("  this is an SGS basis \n");
  if (ss_basis && prl>1)  printf("  this is an SS basis \n");

  double Zg = read_float("ZG");
  double ztg = read_float("ZTG");
  double Rc = read_float("RC");
  //int jorder = 2;

  int gs1 = 0;
  if (Rc>0.)
    gs1 = get_gs1(nrad,nang,grid,Rc);

  //#include "jsetup.cpp"
  int gs2 = gsa;

  int N = basis.size();
  int N2 = N*N;
  int* n2i = new int[natoms];
  int imaxN = get_imax_n2i(natoms,N,basis,n2i);

  float norm[N];
  for (int i=0;i<N;i++)
    norm[i] = basis[i][4];

 #if 0
  float* Paon = new float[N2];
  for (int i=0;i<N;i++)
  for (int j=0;j<N;j++)
    Paon[i*N+j] = Pao[i*N+j]*norm[i]*norm[j];

  if (prl>1)
  {
    printf("\n Pao: \n");
    print_square(N,Pao);
    printf("\n");
  }
 #endif

  int iN = N;
  float** val1 = new float*[iN];
  for (int i=0;i<iN;i++)
    val1[i] = new float[gsa];

  float* grid1 = new float[gsa6];

  #pragma acc enter data copyin(norm[0:N])
  #pragma acc enter data create(val1[0:iN][0:gsa])
  #pragma acc enter data create(grid1[0:gsa6])

 #if 0
  #pragma acc parallel loop present(rho[0:gsa])
  for (int j=0;j<gsa;j++)
    rho[j] = 0.f;
 #endif

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

      if (ss_basis)
        eval_ss(ii1,gs2,grid,val1[ii1],n1,l1,m1,zeta1,Rc);
      else if (sgs_basis)
        eval_sgs(ii1,gs1,gs2,grid,val1[ii1],n1,l1,m1,zeta1,Rc);
      else
        eval_sh(ii1,gsa,grid1,val1[ii1],n1,l1,m1,zeta1);
    }

  } //loop m over natoms


 #if 0
 //rho
  for (int i1=0;i1<N;i1++)
  {
    int ii1 = i1;
    float* valn = val1[ii1];

    for (int i2=0;i2<N;i2++)
    {
      float d1 = Paon[i1*N+i2];

      int ii2 = i2;
      float* valm = val1[ii2];

      #pragma acc parallel loop present(rho[0:gsa],gf[0:gsa],TL[0:gsa],valm[0:gsa],valn[0:gsa],valL[0:gsa])
      for (int j=0;j<gsa;j++)
        rho[j] += d1*valn[j]*valm[j];

    } //loop i2
  } //loop i1

  double den = 0.;
  #pragma acc parallel loop present(rho[0:gsa],wt[0:gsa]) reduction(+:den)
  for (int j=0;j<gsa;j++)
    den += rho[j]*wt[j];
  printf("  total den: %12.8f \n",den);

  if (prl>1 && gsa<20000)
  {
    #pragma acc update self(rho[0:gsa])
    printf("\n rho: \n");
    for (int n=0;n<natoms;n++)
    {
      for (int j=0;j<gs;j+=pr_inc)
        printf(" %8.5f",rho[n*gs+j]);
      printf("\n");
    }
  }
 #endif

  //float thresh = 1.e-20f;

  double* tm1 = new double[gsa];
  #pragma acc enter data create(tm1[0:gsa])

  printf("  z:");
  if (natoms>1)
  for (int k=0;k<gsa;k++)
  {
    float x1 = grid[6*k+0]; float y1 = grid[6*k+1]; float z1 = grid[6*k+2];
    if (x1==0. && y1==0.)
      printf(" %9.5f",z1);
  }
  else
  for (int k=0;k<gsa;k++)
  {
    float x1 = grid[6*k+0]; float y1 = grid[6*k+1]; float z1 = grid[6*k+2];
    if (x1==0. && y1==0. && z1>0.)
      printf(" %9.5f",z1);
  }
  printf("\n");

 //MOs
  string filename = "MO_ON_GRID";
  ofstream outfile;
  outfile.open(filename.c_str());
  //outfile << right << setw(6) << "MO" << setw(12) << "x" 
  //        << setw(12) << "y" << setw(12) << "z" << setw(12) << "value" << '\n';
  
  outfile << fixed << setprecision(5);
  
  for (int i1=0;i1<N;i1++)
  {
    #pragma acc parallel loop present(tm1[0:gsa])
    for (int j=0;j<gsa;j++)
      tm1[j] = 0.f;

    #pragma acc parallel loop present(tm1[0:gsa],val1[0:iN][0:gsa],norm[0:N],jCA[0:N2])
    for (int j=0;j<gsa;j++)
    {
      float v1 = 0.f;
      #pragma acc loop reduction(+:v1)
      for (int k=0;k<N;k++)
        v1 += jCA[k*N+i1]*norm[k]*val1[k][j]; 
      tm1[j] = v1;                            
    }

    float sign1 = 1.;
   #pragma acc serial present(tm1[0:gsa])
    if (tm1[0]<0.)
      sign1 *= -1.;

    #pragma acc update self(tm1[0:gsa])
    // printf("  MO%i:",i1);
    outfile << "MO " << i1 << ":\n";
    //if (natoms>1)
    for (int k=0;k<gsa;k++)
    {
      float x1 = grid[6*k+0]; float y1 = grid[6*k+1]; float z1 = grid[6*k+2];
      float mo_val = sign1*tm1[k];
      //if (x1==0. && y1==0.)

      outfile << setw(12) << x1 << setw(12) << y1 << setw(12) << z1
              << setw(12) << mo_val << '\n';
    }
  }
  outfile.close();
  printf("\n");

  #pragma acc exit data delete(tm1[0:gsa])

  delete [] tm1;


 //cleanup
  #pragma acc exit data delete(norm[0:N])
  #pragma acc exit data delete(val1[0:iN][0:gsa])
  #pragma acc exit data delete(grid1[0:gsa6])

  //delete [] Paon;
  delete [] n2i;

  delete [] grid1;

  for (int i=0;i<iN;i++)
    delete [] val1[i];
  delete [] val1;

  return;
}

void write_mo_grid(int natoms, int* atno, double* coords, int nrad, int gsa, vector<vector<double> > basis,
                   int No, double* Pao, double* jCA, double* jS, double* grid, double* wt, int prl)
{
  int gs6 = 6*gsa;
  float* gridf = new float[gs6];
  float* wtf = new float[gsa];
  #pragma acc enter data create(gridf[0:gs6],wtf[0:gsa])

  #pragma acc parallel loop present(gridf[0:gs6],grid[0:gs6])
  for (int j=0;j<gs6;j++)
    gridf[j] = grid[j];
  #pragma acc parallel loop present(wtf[0:gsa],wt[0:gsa])
  for (int j=0;j<gsa;j++)
    wtf[j] = wt[j];
  #pragma acc update self(gridf[0:gs6])

  write_mo_grid(natoms,atno,coords,nrad,gsa,basis,No,Pao,jCA,jS,gridf,wtf,prl);

  #pragma acc exit data delete(gridf[0:gs6],wtf[0:gsa])
  delete [] gridf;
  delete [] wtf;
}