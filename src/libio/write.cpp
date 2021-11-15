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

void write_iarray(short type1, short type2, short i1, int s1, int s2, FP1* A)
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

#if !EVL64
void write_iarray(short type1, short type2, short i1, int s1, int s2, FP2* A)
{
  printf("\n shouldn't be here (for now) \n"); 

  string filename = get_iarray_name(type1,type2,i1);
  printf("  writing %8s \n",filename.c_str());

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
#endif

void save_geoms(int natoms, int* atno, vector<FP1> E, vector<FP1*> geoms, string fname)
{
  int ng = geoms.size();

  FP2 B2A = 1./A2B;

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());

  outfile << fixed << setprecision(8);

  for (int i=0;i<ng;i++)
  {
    FP1* coords = geoms[i];
    outfile << natoms << endl << E[i] << endl;
    for (int j=0;j<natoms;j++)
    {
      string aname = get_aname(atno[j]);
      FP2 x1 = coords[3*j]*B2A; FP2 y1 = coords[3*j+1]*B2A; FP2 z1 = coords[3*j+2]*B2A;
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
      if (m1==0)      nx = 1;
      else if (m1==1) ny = 1;
      else            nz = 1;
    }
  }
  else if (n1==3)
  {
    if (l1==0)
      nr = 2;
    else if (l1==1)
    {
      nr = 1;
      if (m1==0)      nx = 1;
      else if (m1==1) ny = 1;
      else            nz = 1;
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
    printf(" WARNING: n=4 MO print incomplete \n");
    if (l1==0)
      nr = 3;
    else if (l1==1)
    {
      nr = 2;
      if (m1==0)      nx = 1;
      else if (m1==1) ny = 1;
      else            nz = 1;
    }
    else if (l1==3)
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
  }
  return;
}


void write_molden(int natoms, int* atno, FP2* coords, vector<vector<FP2> > &basis, FP2* jCA, int No, string fname)
{
  int N = basis.size();

  for (int j=0;j<N;j++)
  if (basis[j][0]>3)
    return;

  printf(" WARNING: cannot print spherical d functions to molden \n"); 

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());

  FP2 B2A = 1./A2B;

  outfile << fixed << setprecision(8);

  outfile << "[Molden Format]" << endl;
  outfile << "[Atoms] (Angs)" << endl;
  for (int i=0;i<natoms;i++)
  {
    string aname = get_aname(atno[i]);
    FP2 x1 = coords[3*i]*B2A; FP2 y1 = coords[3*i+1]*B2A; FP2 z1 = coords[3*i+2]*B2A;
    outfile << " " << aname << "     " << i+1 << "   " << atno[i] << "   ";
    outfile << x1 << "   " << y1 << "   " << z1 << endl;
  }
  

  outfile << "[STO]" << endl;

  for (int i=0;i<natoms;i++)
  {
    FP2 x1 = coords[3*i]; FP2 y1 = coords[3*i+1]; FP2 z1 = coords[3*i+2];

    for (int j=0;j<N;j++)
    {
      if (x1==basis[j][5] && y1==basis[j][6] && z1==basis[j][7])
      {
        int n1 = basis[j][0]; int l1 = basis[j][1]; int m1 = basis[j][2]; FP2 zeta = basis[j][3]*B2A; FP2 norm = basis[j][4];
        int nx,ny,nz,nr;
        get_nxyzr(n1,l1,m1,nx,ny,nz,nr);

        outfile << i+1 << " " << nx << " " << ny << " " << nz << " ";
        outfile << nr << " " << zeta << " " << norm << " " << endl;
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

    for (int j=0;j<N;j++)
      outfile << " " << j+1 << "  " << jCA[j*N+i] << endl;
  }


  outfile.close();
  
  return;
}

#if !EVL64
void write_square(int N, FP2* A, string fname, int prl)
{
  if (prl>1) printf(" writing %s to file \n",fname.c_str());

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
#endif

void write_square(int N, FP1* A, string fname, int prl)
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

void write_C(int Naux, int N2, FP1* C)
{
  printf(" writing C to file \n");

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

#if !EVL64
void write_C(int Naux, int N2, FP2* C)
{
  printf(" writing C to file \n");

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
#endif

void write_S_En_T(int N, FP1* S, FP1* En, FP1* T)
{
  printf(" writing S/En/T to file \n");
   
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

#if !EVL64
void write_S_En_T(int N, FP2* S, FP2* En, FP2* T)
{
  printf(" writing S/En/T to file \n");
   
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
#endif
