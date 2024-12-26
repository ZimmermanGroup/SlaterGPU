#include "write.h"

#include <sstream>
#define SSTRF( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << fixed << setprecision(8) << scientific << x ) ).str()
#define SSTRF2( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << fixed << setprecision(14) << scientific << x ) ).str()
//#define SSTRF( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << fixed << setprecision(8) << x ) ).str()
//#define SSTRF2( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << fixed << setprecision(14) << x ) ).str()

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

void write_molden_g(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, string fname)
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
    outfile << "Ene=0.0" << endl;
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

void write_molden(bool gbasis, int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, string fname)
{
  if (basis[0].size()>10) return write_molden_ss(natoms,atno,coords,basis,jCA,No,fname);
  if (gbasis) return write_molden_g(natoms,atno,coords,basis,jCA,No,fname);

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

void save_dft_vals(int save_type, double thresh, int natoms, int nrad, int nang, double* grid, double* rho, double* drho, double* Td, double* epsi, double* vc, int zpos, string filename)
{
  //if (!save_radial) { printf(" ERROR: save_dft_vals not ready for angular component \n"); return; }

  int gs = nrad*nang;
  int gsa = gs*natoms;

  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(12);

  if (save_type==1)
  {
    outfile << "DFT.  nrad: " << nrad << "  grid/rho/|drho|/Td/vc" << endl;

   //radial component only (single atom)
    if (epsi!=NULL)
    for (int j=0;j<nrad;j++)
    {
      int j1 = j*nang;
      //if (rho[j1]>thresh)
      {
        string line = SSTRF2(grid[6*j1+3]) + "," + SSTRF2(rho[j1]) + "," + SSTRF2(drho[j1]) + "," + SSTRF2(Td[j1]) + "," + SSTRF2(epsi[j1]) + "," + SSTRF2(vc[j1]);
        outfile << line << endl;
      }
    }
    else
    for (int j=0;j<nrad;j++)
    {
      int j1 = j*nang;
      //if (rho[j1]>thresh)
      {
        string line = SSTRF2(grid[6*j1+3]) + "," + SSTRF2(rho[j1]) + "," + SSTRF2(drho[j1]) + "," + SSTRF2(Td[j1]) + "," + SSTRF2(vc[j1]);
        outfile << line << endl;
      }
    }
  }
  else if (save_type==2)
  {
    outfile << "DFT.  yz: " << 2*nrad << " grid/rho/|drho|/Td/vc" << endl;
   //within yz plane
    if (epsi!=NULL)
    for (int j=0;j<gsa;j++)
    {
      if (fabs(grid[6*j+0])<1.e-12 && grid[6*j+1]>=0.)
      //if (fabs(grid[6*j+0])<1.e-12 && grid[6*j+1]>=0. && rho[j]>thresh)
      {
        string line = SSTRF2(grid[6*j+2]) + "," + SSTRF2(rho[j]) + "," + SSTRF2(drho[j]) + "," + SSTRF2(Td[j]) + "," + SSTRF2(epsi[j]) + "," + SSTRF2(vc[j]);
        outfile << line << endl;
      }
    }
    else
    for (int j=0;j<gsa;j++)
    {
      if (fabs(grid[6*j+0])<1.e-12 && grid[6*j+1]>=0.)
      //if (fabs(grid[6*j+0])<1.e-12 && grid[6*j+1]>=0. && rho[j]>thresh)
      {
        string line = SSTRF2(grid[6*j+2]) + "," + SSTRF2(rho[j]) + "," + SSTRF2(drho[j]) + "," + SSTRF2(Td[j]) + "," + SSTRF2(vc[j]);
        outfile << line << endl;
      }
    }
  }
  else if (save_type==3)
  {
    outfile << "DFT.  z: " << 2*nrad*natoms << " grid/rho/|drho|/Td/vc" << endl;
   //along z axis
    if (epsi!=NULL)
    for (int j=0;j<gsa;j++)
    {
      if (fabs(grid[6*j+0])<1.e-12 && fabs(grid[6*j+1])<1.e-12)
      //if (fabs(grid[6*j+0])<1.e-12 && fabs(grid[6*j+1])<1.e-12 && rho[j]>thresh)
      {
        string line = SSTRF2(grid[6*j+2]) + "," + SSTRF2(rho[j]) + "," + SSTRF2(drho[j]) + "," + SSTRF2(Td[j]) + "," + SSTRF2(epsi[j]) + "," + SSTRF2(vc[j]);
        outfile << line << endl;
      }
    }
    else
    for (int j=0;j<gsa;j++)
    {
      if (fabs(grid[6*j+0])<1.e-12 && fabs(grid[6*j+1])<1.e-12)
      //if (fabs(grid[6*j+0])<1.e-12 && fabs(grid[6*j+1])<1.e-12 && rho[j]>thresh)
      {
        string line = SSTRF2(grid[6*j+2]) + "," + SSTRF2(rho[j]) + "," + SSTRF2(drho[j]) + "," + SSTRF2(Td[j]) + "," + SSTRF2(vc[j]);
        outfile << line << endl;
      }
    }
  }

  outfile.close();
  return;
}

void save_dft_vals(bool save_type, double thresh, int natoms, int nrad, int nang, double* grid, double* rho, double* drho, double* Td, double* epsi, double* vc, int zpos, string filename)
{
  return save_dft_vals(save_type,thresh,natoms,nrad,nang,grid,rho,drho,Td,epsi,vc,zpos,filename);
}

void write_grid(int wc, int nrad, int nang, double* grid, double* wt)
{
  string filename = "GRID_WTS";
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(10);

  outfile << filename << ":" << endl;
  for (int i=0;i<nrad;i++)
  {
    for (int j=0;j<nang;j++)
    {
      int i1 = i*nang+j;
      string line = "r: " + SSTRF(grid[wc*i1+3]) + " xyz: ";
      line += " " + SSTRF(grid[wc*i1+0]) + " " + SSTRF(grid[wc*i1+1]) + " " + SSTRF(grid[wc*i1+2]);
      line += " wt: " + SSTRF(wt[i1]);
      outfile << line << endl;
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

