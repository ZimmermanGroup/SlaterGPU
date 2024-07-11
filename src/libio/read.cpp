//use only m=0 components of spherical harmonics
//do this for l>#
//#define CYLINDER_Z 1
#define CYLINDER_Z 6
//#define CYLINDER_Z 10

#define GEOM_DEBUG 0
#define PDEBUG 0
#define DDEBUG 0
#define FDEBUG 0
#define GDEBUG 0
#define HDEBUG 0

#include "read.h"

#include <sstream>
#define SSTRF( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << fixed << setprecision(14) << x ) ).str()
#define SSTRF2( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << fixed << setprecision(4) << x ) ).str()

double norm(int n, int l, int m, double zeta);

#define A2B 1.8897261

   
int check_file(string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;
     
  infile.close();
   
  return 1;   
}  

void print_coords(int natoms, float* coordsf)
{
  float B2A = 1.f/A2B;
  for (int n=0;n<natoms;n++)
    printf(" %10.5f %10.5f %10.5f \n",coordsf[3*n+0]*B2A,coordsf[3*n+1]*B2A,coordsf[3*n+2]*B2A);
}

void print_gradient(int natoms, double* grad)
{
  for (int n=0;n<natoms;n++)
    printf(" %10.5f %10.5f %10.5f \n",grad[3*n+0],grad[3*n+1],grad[3*n+2]);
}

void print_square_diff(int N, double* S1, double* S2) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S1[n*N+m],S1[n*N+m]);
    printf("\n");
  }
}

void print_square_fine(int N, float* S) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %12.8f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_fine(int N, double* S) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %12.8f",S[n*N+m]);
    printf("\n");
  }
}

void print_square(int N, double* S) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square(int N, float* S) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square(int M, int N, float* S) 
{
  for (int n=0;n<M;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square(int M, int N, double* S) 
{
  for (int n=0;n<M;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_mos_col(int M, int N, vector<vector<double> > basis, double* jCA)
{
  if (M<1) return;
  for (int j=0;j<N;j++)
  {
    int n1 = basis[j][0]; int l1 = basis[j][1]; int m1 = basis[j][2];

    printf("  %i%i%2i ",n1,l1,m1);
    for (int k=0;k<M;k++)
      printf(" %10.5f",jCA[j*N+k]);
    printf("\n");
  }
}

void print_square_col(int M, int N, double* S) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<M;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_col_sm(int M, int N, double* S) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<M;m++)
      printf(" %6.3f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_ss(int N, double* S)
{
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %8.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_sm(int N, double* S)
{
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %5.2f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_sm(int N, float* S)
{
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %5.2f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_ss_sm(int N, double* S)
{ 
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %6.3f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_nxn(int No, int N, float* S)
{
  for (int n=0;n<No;n++)
  {
    for (int m=0;m<No;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_nxn(int No, int N, double* S)
{
  for (int n=0;n<No;n++)
  {
    for (int m=0;m<No;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_rectangle(int N1, int N2, double* S)
{
  for (int n=0;n<N1;n++)
  {
    printf("  ");
    for (int m=0;m<N2;m++)
      printf(" %10.5f",S[n*N2+m]);
    printf("\n");
  }
}

void print_rectangle_e(int N1, int N2, double* S)
{
  for (int n=0;n<N1;n++)
  {
    printf("  ");
    for (int m=0;m<N2;m++)
      printf(" %5.1e",S[n*N2+m]);
    printf("\n");
  }
}

void print_rectangle_sm(int N1, int N2, double* S)
{
  for (int n=0;n<N1;n++)
  {
    printf("  ");
    for (int m=0;m<N2;m++)
      printf(" %7.4f",S[n*N2+m]);
    printf("\n");
  }
}


vector<string> split1(const string &s, char delim)
{
  stringstream ss(s);
  string item;  
  vector<string> tokens;
  while ((bool)getline(ss, item, delim))
    tokens.push_back(item);
  return tokens;
}

string read_basis_text(string aname)
{
  string filename = "basis";
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    printf("  couldn't open file: %s \n",filename.c_str());
    return 0;
  }

  bool done = 0;
  string text = "";
  while (!infile.eof() && !done)
  {  
    string line;
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>1 && tok_line[0].compare(aname)==0)
    {
      //printf("  found elem %2s \n",aname.c_str());
      while (!infile.eof())
      {
        (bool)getline(infile, line);
        vector<string> tok2 = split1(line,' ');
        if (tok2.size() & tok2[0].compare("****")==0)
        {
          done = 1;
          break;
        }
        text += line;
        text += "\n";
      }
    }
  }

  infile.close();

  return text;
}

double read_float(string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0.;

  string line;
  bool success = (bool)getline(infile, line);

  double val = 0.;
  if (success)
    val = atof(line.c_str());

  infile.close();

  return val;
}

double read_float_2(string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return -1.;

  string line;
  bool success = (bool)getline(infile, line);

  double val = -1.;
  if (success)
    val = atof(line.c_str());

  infile.close();

  return val;
}

int read_int(string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);

  int val = 0;
  if (success)
    val = atoi(line.c_str());

  infile.close();

  return val;
}

int read_int_2(string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return -1;

  string line;
  bool success = (bool)getline(infile, line);

  int val = 0.;
  if (success)
    val = atoi(line.c_str());

  infile.close();

  return val;
}

void read_thresh(float& no_thresh, float& occ_thresh)
{
  string filename = "THRESH";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return;

  string line;
  bool success = (bool)getline(infile, line);

  if (success)
  {
    no_thresh = atof(line.c_str());

    success = (bool)getline(infile, line);
    if (success)
      occ_thresh = atof(line.c_str());
  }

  infile.close();

  return;
}

void read_eps(double& eps1, double& eps2, double& eps1s)
{
  string filename = "EPS";

  eps1 = 0.001;
  eps2 = 0.000001;
  eps1s = 0.0001;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return;

  string line;
  bool success = (bool)getline(infile, line);

  if (success)
  {
    eps1 = atof(line.c_str());

    success = (bool)getline(infile, line);
    if (success)
      eps2 = atof(line.c_str());

    success = (bool)getline(infile, line);
    if (success)
      eps1s = atof(line.c_str());
  }

  infile.close();

  if (eps1<1.e-6)
    eps1 = 1.e-6;
  if (eps2<1.e-12)
    eps2 = 1.e-12;
  if (eps1s<0.)
    eps1s = 0.;

  return;
}

void read_eps(double& eps1, double& eps2)
{
  string filename = "EPS";

  eps1 = 0.001;
  eps2 = 0.000001;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return;

  string line;
  bool success = (bool)getline(infile, line);

  if (success)
  {
    eps1 = atof(line.c_str());

    success = (bool)getline(infile, line);
    if (success)
      eps2 = atof(line.c_str());
  }

  infile.close();

  if (eps1<1.e-6)
    eps1 = 1.e-6;
  if (eps2<1.e-12)
    eps2 = 1.e-12;

  return;
}

bool read_dft(int& dt1, int& dt2, int type)
{
  string filename = "DFT";
  if (type==2) filename = "DFTS";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);

  dt1 = dt2 = -1;
  if (success)
  {
   //exchange
    dt1 = atoi(line.c_str());

    success = (bool)getline(infile, line);
    if (success)
     //correlation
      dt2 = atoi(line.c_str());
  }
  else
  {
    infile.close();
    return 0;
  }

  infile.close();

  return 1;
}

int read_esci()
{
  string filename = "ESCI";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);

  int default_val = 1;
  int do_esci;
  if (success)
    do_esci = atoi(line.c_str());
  else
    do_esci = default_val;

  infile.close();

  return do_esci;
}

int read_cusp()
{
  string filename = "CUSP";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 1;

  string line;
  bool success = (bool)getline(infile, line);

  int default_val = 1;
  int cusp;
  if (success)
    cusp = atoi(line.c_str());
  else
    cusp = default_val;

  infile.close();

  return cusp;
}

int read_hfx()
{
  string filename = "HFX";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);

  int default_val = 1;
  int hfx;
  if (success)
    hfx = atoi(line.c_str());
  else
    hfx = default_val;

  infile.close();

  return hfx;
}

int read_pp()
{
  string filename = "PP";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);

  int default_val = 1;
  int pp;
  if (success)
    pp = atoi(line.c_str());
  else
    pp = default_val;

  infile.close();

  return pp;
}

int read_symm()
{
  string filename = "SYMM";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);

  int default_val = 1;
  int pp;
  if (success)
    pp = atoi(line.c_str());
  else
    pp = default_val;

  infile.close();

  return pp;
}

int read_vnuc()
{
  string filename = "VNUC";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);

  int default_val = 1;
  int vnuc;
  if (success)
    vnuc = atoi(line.c_str());
  else
    vnuc = default_val;

  infile.close();

  return vnuc;
}

int read_ri()
{
  string filename = "RI";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);

  int default_val = 1;
  int do_ri;
  if (success)
    do_ri = atoi(line.c_str());
  else
    do_ri = default_val;

  infile.close();

  return do_ri;
}

int read_opt()
{
  string filename = "OPT";

  int default_val = 1;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return default_val;

  string line;
  bool success = (bool)getline(infile, line);

  int nopt;
  if (success)
    nopt = atoi(line.c_str());
  else
    nopt = default_val;

  infile.close();

  return nopt;
}

int read_mbe()
{
  string filename = "MBE";

  int default_val = 3;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return default_val;

  string line;
  bool success = (bool)getline(infile, line);

  int mbe;
  if (success)
    mbe = atoi(line.c_str());
  else
    mbe = default_val;

  infile.close();

  return mbe;
}


int read_basis()
{
  string filename = "basis";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  infile.close();

  return 1;
}

int read_hf()
{
  string filename = "HF";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  infile.close();

  return 1;
}

int read_lag()
{
  string filename = "LAG";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  infile.close();

  return 1;
}

void read_T(int& use_td, int& use_tl, int& use_erf, int& use_th, int& use_mixed_t)
{
  string filename = "T";

  int default_val = 0;
  use_td = use_tl = use_th = 0;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return;

  string line;
  bool success = (bool)getline(infile, line);

  int type;
  if (success)
    type = atoi(line.c_str());
  else
    type = default_val;

  if (type==1)
    use_td = 1;
  else if (type==2)
    use_tl = 1;
  else if (type==3)
  {
    use_tl = 1;
    use_erf = 1;
  }
  else if (type==4)
    use_td = 2;
  else if (type==5)
  {
    use_td = 1;
    use_mixed_t = 1;
  }
  else if (type==6)
    use_th = 1;

  infile.close();

  return;
}

int read_cas()
{
  string filename = "CAS";

  int default_val = 0;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return default_val;

  string line;
  bool success = (bool)getline(infile, line);

  int cas_type;
  if (success)
    cas_type = atoi(line.c_str());
  else
    cas_type = default_val;

  infile.close();

  return cas_type;
}

void read_cas_size(int& Nc, int& Na, int& Nb)
{
  string filename = "CAS_2";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return;

  string line;
  bool success; 

  success = (bool)getline(infile, line);
  if (success)
    Nc = atoi(line.c_str());

  success = (bool)getline(infile, line);
  if (success)
    Na = Nc+atoi(line.c_str());

  success = (bool)getline(infile, line);
  if (success)
    Nb = Nc+atoi(line.c_str());

  infile.close();

  return;
}

bool read_array(int size, double* A, string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);

  if (success)
  {
    vector<string> tok_line = split1(line,' ');
    //printf("  size: %2i tok_line.size: %2i \n",size,tok_line.size());
    if (tok_line.size()<size) { infile.close(); return 0; }

    for (int i=0;i<tok_line.size();i++)
      A[i] = atof(tok_line[i].c_str());
  }

  infile.close();

  return success;
}

vector<double> read_vector(string filename)
{
  vector<double> vec;
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return vec;

  string line;
  bool success = (bool)getline(infile, line);

  if (success)
  {
    vector<string> tok_line = split1(line,' ');
    for (int i=0;i<tok_line.size();i++)
      vec.push_back(atof(tok_line[i].c_str()));
  }

  infile.close();

  return vec;
}

bool read_cas_act(int& N, int& M)
{
  string filename = "CAS_ACT";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);

  if (success)
  {
    vector<string> tok_line = split1(line,' ');
    N = atoi(tok_line[0].c_str());
    M = atoi(tok_line[1].c_str());
  }
  else
  {
    printf("\n ERROR reading CAS_ACT \n");
    exit(1);
  }

  infile.close();

  return 1;
}

int read_rotate(int N, double* jCA)
{
  string filename = "ROT";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  printf("\n jCA(in) \n");
  print_square_sm(N,jCA);

  int nrot = 0;
  while (!infile.eof())
  {
    string line;
    bool success = (bool)getline(infile, line);

    if (success)
    {
      double vec1[N];
      double vec2[N];
      vector<string> tok_line = split1(line,' ');

      if (tok_line.size()==2)
      {
        int i1 = atoi(tok_line[0].c_str())-1;
        int i2 = atoi(tok_line[1].c_str())-1;

        printf("   rotate %2i->%2i \n",i1+1,i2+1);
        for (int i=0;i<N;i++) vec1[i] = jCA[i*N+i1];
        for (int i=0;i<N;i++) vec2[i] = jCA[i*N+i2];
        for (int i=0;i<N;i++) jCA[i*N+i1] = vec2[i];
        for (int i=0;i<N;i++) jCA[i*N+i2] = vec1[i];
        nrot++;
      }
    }
  }

  infile.close();

  printf("\n  rotated %2i orbitals \n",nrot);

  printf("\n jCA(out) \n");
  print_square_sm(N,jCA);

  return 1;
}

int read_nsteps()
{
  string filename = "NSTEPS";

  int default_nsteps = 1;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return default_nsteps;

  string line;
  bool success = (bool)getline(infile, line);

  int nsteps;
  if (success)
    nsteps = atoi(line.c_str());
  else
    nsteps = default_nsteps;

  infile.close();

  return nsteps;
}

int read_spinref()
{
  string filename = "SPINREF";

  int spinref;
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return -1;
  else
  {
    string line;
    bool success = (bool)getline(infile, line);

    if (success)
      spinref = atoi(line.c_str());
    else
      spinref = 0;
  }
  infile.close();
 
  return spinref;
}

int read_group()
{
  string filename = "GROUP";

  int group_size;
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;
  else
  {
    string line;
    bool success = (bool)getline(infile, line);

    if (success)
      group_size = atoi(line.c_str());
    else
      group_size = 1;
  }
  infile.close();
 
  return group_size;
}

int read_restart()
{
  string filename = "RESTART";

  int is_restart = 1;
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;
  else
  {
    string line;
    bool success = (bool)getline(infile, line);

    if (success)
      is_restart = atoi(line.c_str());
    else
      is_restart = 1;

    //printf(" found restart: %i \n",is_restart);
  }
  infile.close();
 
  return is_restart;
}

void read_gridps(int& nmu, int& nnu, int& nphi, int type)
{
  string filename = "GRIDPS";
  if (type==2) filename = "GRIDPS2";

  nmu  = 8; //default
  nnu  = 8; //default
  nphi = 8; //default

  ifstream infile;
  infile.open(filename.c_str());    
  if (!infile)
  {
    printf("  couldn't open GRIDPS file. \n");
    nmu = nnu = nphi = 0;
    return;
  }
  
  string line;
  bool success = (bool)getline(infile, line);
  if (success)
    nmu = atoi(line.c_str());

  success = (bool)getline(infile, line);
  if (success)
    nnu = atoi(line.c_str());

  success = (bool)getline(infile, line);
  if (success)
    nphi = atoi(line.c_str());

  infile.close();

  return;
}


void read_quad(int& qo1, int& qo2)
{
  string filename = "QUAD";
  //if (type==2) filename = "QUAD2";

  qo1 = 8;
  qo2 = 8;

  ifstream infile;
  infile.open(filename.c_str());    
  if (!infile)
    return;
  
  string line;
  bool success = (bool)getline(infile, line);
  if (success)
    qo1 = atoi(line.c_str());

  success = (bool)getline(infile, line);
  if (success)
    qo2 = atoi(line.c_str());
  else
    qo2 = qo1;

  if (qo2<=0) qo2 = qo1;

  infile.close();

  return;
}

void read_nrad_nang(int& nrad, int& nang, int type)
{
  string filename = "GRID";
  if (type==2) filename = "GRID2";

  ifstream infile;
  infile.open(filename.c_str());    
  if (!infile)
  {
    printf("  couldn't open GRID file. please provide GRID \n");
    exit(1);
  }
  
  string line;
  bool success = (bool)getline(infile, line);
  if (success)
    nrad = atoi(line.c_str());

  success = (bool)getline(infile, line);
  if (success)
    nang = atoi(line.c_str());

  infile.close();

  return;
}

int read_tuples(vector<vector<int> >& tuples)
{
  string filename = "TUPLES";
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  printf(" TUPLES only reading one line (for now) \n");

  string line;
  bool success = (bool)getline(infile, line);
  vector<string> tok_line = split1(line,' ');
  if (tok_line.size()>0)
  {
    vector<int> tuple1;
    for (int m=0;m<tok_line.size();m++)
      tuple1.push_back(atoi(tok_line[m].c_str()));
    tuples.push_back(tuple1);
  }

  infile.close();

  return tuples.size();
}


string get_iarray_name(short type1, short type2, short i1);

int read_iarray(short type1, short type2, short i1, int s1, int s2, float* A)
{
  string filename = get_iarray_name(type1,type2,i1);
  //printf("  reading %8s \n",filename.c_str());

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    printf("   couldn't open file (%s) \n",filename.c_str());
    return 0;
  }

  string line;
  int found = 1;

  for (int i=0;i<s1;i++)
  {
    bool success = (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');

    int tks = tok_line.size();
    if (tks>=s2)
    {
      for (int m=0;m<s2;m++)
        A[i*s2+m] = atof(tok_line[m].c_str());
    }
    else
      found = 0;
  }

  infile.close();

  return found;
}

int read_iarray(short type1, short type2, short i1, int s1, int s2, double* A)
{
  string filename = get_iarray_name(type1,type2,i1);
  //printf("  reading %8s \n",filename.c_str());

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    printf("   couldn't open file (%s) \n",filename.c_str());
    return 0;
  }

  string line;
  int found = 1;

  for (int i=0;i<s1;i++)
  {
    bool success = (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');

    int tks = tok_line.size();
    if (tks>=s2)
    {
      for (int m=0;m<s2;m++)
        A[i*s2+m] = atof(tok_line[m].c_str());
    }
    else
      found = 0;
  }

  infile.close();

  return found;
}

int read_gridpts(int s1, int s2, float* A, string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool)getline(infile, line);
  int wi = 0;
  while (wi<s1)
  {
    success = (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<s2) { printf(" ERROR: file size incorrect (read_gridpts. file: %s) \n",filename.c_str()); exit(1); }
      for (int m=0;m<s2;m++)
        A[wi*s2+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }
  
  infile.close();

  return 1;
}

int read_square_check(int N, double* A, string filename)
{

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    printf("  couldn't open file %s \n",filename.c_str());
    return 0;
  }

  string line;
  bool success = (bool)getline(infile, line);
  int wi = 0;
  while (wi<N)
  {
    success = (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { infile.close(); return 0; }
      for (int m=0;m<N;m++)
        A[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
    else
      return 0;
  }
  
  infile.close();

  return 1;
}

int read_square(int N, double* A, string filename)
{

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    //printf("  couldn't open file \n");
    return 0;
  }

  string line;
  bool success = (bool)getline(infile, line);
  int wi = 0;
  while (wi<N)
  {
    success = (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file (%s) size not square \n",filename.c_str()); exit(1); }
      for (int m=0;m<N;m++)
        A[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }
  
  infile.close();

  return 1;
}

int read_square(vector<vector<double> > basis, double* A, string filename)
{
  int N = basis.size();
  return read_square(N,A,filename);
}

void read_SENT(string dirname, int N, double* S, double* T, double* En, int prl)
{
  string filename = "SENT"; if (dirname!="0") filename = dirname+"/"+"SENT";
  printf("  attempting file read: %s \n",filename.c_str());

  ifstream infile;
  infile.open(filename.c_str());    
  if (!infile)
  {
    printf("  couldn't open file \n");
    return;
  }

  string line;  
  (bool)getline(infile, line);
  int wi = 0;
  while (wi<N)
  {
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size wrong in S \n"); exit(1); }
      for (int m=0;m<N;m++)
        S[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }

  if (En==NULL) { infile.close(); return; }

  (bool)getline(infile, line);
  wi = 0;
  while (wi<N)
  {
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size wrong in En \n"); exit(1); }
      for (int m=0;m<N;m++)
        En[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }

  if (T==NULL) { infile.close(); return; }

  (bool)getline(infile, line);
  wi = 0;
  while (wi<N)
  {
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size wrong in T \n"); exit(1); }
      for (int m=0;m<N;m++)
        T[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }

  infile.close();

  if (prl>1)
  {
    printf("\n S: \n"); print_square(N,S);
    printf(" En: \n"); print_square(N,En);
    printf(" T: \n"); print_square(N,T);
    printf("\n");
  }

}

bool read_yukawa_potentials(int N, int Naux, double*& Ayd, double*& Cyd)
{
  int N2 = N*N;
  int N2a = N2*Naux;
  int Na2 = Naux*Naux;

  bool yukawa_available = check_file("Ay");
  if (yukawa_available)
  {
    printf("  reading Yukawa potentials \n");
    Ayd = new double[Na2];
    Cyd = new double[N2a];

    read_square(Naux,Ayd,"Ay");
    read_Ciap(N,Naux,Cyd,"Cyiap");
  }

  return yukawa_available;
}

void read_Ciap(int N, int Naux, double* Ciap, string filename)
{
  printf("  attempting file read: %s \n",filename.c_str());
  int prl = 0;

  ifstream infile;
  infile.open(filename.c_str());    
  if (!infile)
  {
    printf("  couldn't open file \n");
    return;
  }
  
  string line;

  (bool)getline(infile, line);
  int wi = 0;
  while (!infile.eof())
  {
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<Naux+1) { printf(" ERROR: file size wrong in Ciap (%2i vs %2i) \n",tok_line.size(),Naux); exit(1); }
      for (int m=0;m<Naux;m++)
        Ciap[wi*Naux+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }
  int N2 = N*N;
  if (wi==N2) { if (prl>1) printf("   found all lines of Ciap \n"); }
  else printf(" Ciap missing lines \n");
  infile.close();

  return;
}

bool read_integrals(string dirname, int N, int Naux, double* S, double* T, double* En, double* A, double* Ciap, int prl)
{
  int N2 = N*N;
  //int N3 = N2*N;
  //int na = Naux;
  //int Naux2 = Naux*Naux;

 //first file
  string filename = "Ciap"; if (dirname!="0") filename = dirname+"/"+"Ciap";
  read_Ciap(N,Naux,Ciap,filename);

 //second file
  filename = "A"; if (dirname!="0") filename = dirname+"/"+"A";
  printf("  attempting file read: %s \n",filename.c_str());

  ifstream infile;
  infile.open(filename.c_str());    
  if (!infile)
  {
    printf("  couldn't open file \n");
    return 0;
  }
  
  string line;
  (bool)getline(infile, line);
  int wi = 0;
  while (!infile.eof())
  {
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<Naux+1) { printf(" ERROR: file size wrong in A \n"); exit(1); }
      for (int m=0;m<Naux;m++)
        A[wi*Naux+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }
  if (wi==Naux) { if (prl>1) printf("   found all lines of A \n"); }
  else printf(" A missing lines \n");
  infile.close();


 //third file
  filename = "SENT"; if (dirname!="0") filename = dirname+"/"+"SENT";
  printf("  attempting file read: %s \n",filename.c_str());

  infile.open(filename.c_str());    
  if (!infile)
  {
    printf("  couldn't open file \n");
    return 0;
  }
  
  (bool)getline(infile, line);
  wi = 0;
  while (wi<N)
  {
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size wrong in S \n"); exit(1); }
      for (int m=0;m<N;m++)
        S[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }

  (bool)getline(infile, line);
  wi = 0;
  while (wi<N)
  {
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size wrong in En \n"); exit(1); }
      for (int m=0;m<N;m++)
        En[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }

  (bool)getline(infile, line);
  wi = 0;
  while (wi<N)
  {
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size wrong in T \n"); exit(1); }
      for (int m=0;m<N;m++)
        T[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }

  infile.close();

  if (prl>1)
  {
    printf("\n S: \n"); print_square(N,S);
    printf(" En: \n"); print_square(N,En);
    printf(" T: \n"); print_square(N,T);
    printf("\n");
  }
  if (prl>1)
  {
    printf(" A: \n");
    print_square_sm(Naux,A);
    printf("\n");
  }

  if (prl>1)
  {
    printf(" C: \n");
    for (int i=0;i<N2;i++) 
    {
      for (int j=0;j<Naux;j++)
        printf(" %6.3f",Ciap[i*Naux+j]);
      printf("\n");
    }
    printf("\n");
  }

  return 1;
}

double nuclear_repulsion(int natoms, int* atno, double* coords)
{
  double Enn = 0;

  for (int n=0;n<natoms;n++)
  for (int m=0;m<n;m++)
  {
    double x12 = coords[3*n+0] - coords[3*m+0];
    double y12 = coords[3*n+1] - coords[3*m+1];
    double z12 = coords[3*n+2] - coords[3*m+2];
    double zz = atno[n]*atno[m];
    Enn += zz/sqrt(x12*x12+y12*y12+z12*z12);
  }

  return Enn;
}

double nuclear_repulsion(int natoms, int* atno, float* coordsf)
{
  double Enn = 0;

  for (int n=0;n<natoms;n++)
  for (int m=0;m<n;m++)
  {
    float x12 = coordsf[3*n+0] - coordsf[3*m+0];
    float y12 = coordsf[3*n+1] - coordsf[3*m+1];
    float z12 = coordsf[3*n+2] - coordsf[3*m+2];
    float zz = atno[n]*atno[m];
    Enn += zz/sqrt(x12*x12+y12*y12+z12*z12);
  }

  return Enn;
}

void add_s(vector<vector<double> > &basis1, vector<double> ao1, int n, double zeta)
{
 //ns function
  ao1[0] = n; ao1[1] = 0; ao1[2] = 0;
  ao1[3] = zeta;
  ao1[4] = norm(n,0,0,zeta);
  basis1.push_back(ao1);

  return;
}

void add_p(vector<vector<double> > &basis1, vector<double> ao1, int n, double zeta, int np)
{
 //2px
  ao1[0] = n; ao1[1] = 1; ao1[2] = 1;
  ao1[3] = zeta;
  ao1[4] = norm(n,1,0,zeta);

 #if PDEBUG || CYLINDER_Z<=1
  //z function only
  ao1[2] = 0;
  basis1.push_back(ao1);
  return;
 #endif

  basis1.push_back(ao1);

  if (np<2) return;

 //2py
  ao1[2] = -1; 
  basis1.push_back(ao1);
 //2pz
  ao1[2] = 0;
  basis1.push_back(ao1);

  return;
}

void add_dz2(vector<vector<double> > &basis1, vector<double> ao1, int n, double zeta)
{
 //3d
  ao1[0] = n; ao1[1] = 2;
  ao1[3] = zeta;

  int m = 0;
  {
    ao1[2] = m;
    ao1[4] = norm(n,2,m,zeta);

    basis1.push_back(ao1);
  }

  return;
}

void add_d(vector<vector<double> > &basis1, vector<double> ao1, int n, double zeta, int nd)
{
  if (nd<2) { printf(" ERROR: cannot use nd anymore \n"); exit(1); } 

 //3d
  ao1[0] = n; ao1[1] = 2;
  ao1[3] = zeta;

 #if CART_D
  #if 0
  ao1[4] = norm(n,2,1,zeta);
  ao1[2] = 0;
  basis1.push_back(ao1);
  ao1[2] = 1;
  basis1.push_back(ao1);
  ao1[2] = 2;
  basis1.push_back(ao1);
  #endif

  ao1[4] = norm(n,2,3,zeta);
  ao1[2] = 3;
  basis1.push_back(ao1);
  ao1[2] = 4;
  basis1.push_back(ao1); 
  ao1[2] = 5;
  basis1.push_back(ao1); 
 #else
  for (int m=-2;m<=2;m++)
  {
    ao1[2] = m;
    ao1[4] = norm(n,2,m,zeta);

   #if DDEBUG || CYLINDER_Z<=2
    if (m==0)
   #endif
    basis1.push_back(ao1);
  }
 #endif

  return;
}

void add_f(vector<vector<double> > &basis1, vector<double> ao1, int n, double zeta, int nf)
{
 //4f
  ao1[0] = n; ao1[1] = 3;
  ao1[3] = zeta;

#if FDEBUG
  printf("   WARNING: restricted f functions \n");
#endif

  int bf = 0;
  for (int m=-3;m<=3;m++)
  {
    ao1[2] = m;
    ao1[4] = norm(n,3,m,zeta);

   #if FDEBUG || CYLINDER_Z<=3
    if (m==0)
   #endif
    {
      basis1.push_back(ao1);
      bf++;
    }

    if (bf>nf) break;
  }

  return;
}

void add_g(vector<vector<double> > &basis1, vector<double> ao1, int n, double zeta, int ng)
{
 //5g
  ao1[0] = n; ao1[1] = 4;
  ao1[3] = zeta;

 //high angular momentum will integrate poorly
 // when becke grid has near-zero weights (around < 10^-8)

 #if GDEBUG
  printf("   WARNING: restricted g functions \n");
 #endif

  int bf = 0;
  for (int m=-4;m<=4;m++)
  {
    ao1[2] = m;
    ao1[4] = norm(n,4,m,zeta);

   #if GDEBUG || CYLINDER_Z<=4
    if (m==0)
   #endif
    {
      basis1.push_back(ao1);
    }
            
    bf++; if (bf>ng) break;
  }
       
  return;
}

void add_h(vector<vector<double> > &basis1, vector<double> ao1, int n, double zeta, int nh)
{
 //6h
  ao1[0] = n; ao1[1] = 5;
  ao1[3] = zeta;

 #if CYLINDER_Z<=5
  ao1[2] = 0;
  ao1[4] = norm(n,5,0,zeta);
  basis1.push_back(ao1);
  return;
 #endif
 #if HDEBUG
  ao1[2] = 0;
  ao1[4] = norm(n,5,0,zeta);
  basis1.push_back(ao1);
  return;
 #endif
 
  int bf = 0;
  for (int m=-5;m<=5;m++)
  {
    ao1[2] = m;
    ao1[4] = norm(n,5,m,zeta);
    basis1.push_back(ao1);
 
    bf++; if (bf>nh) break;
  }
 
  return;
}

void get_n_l(string aotype, int& n, int& l)
{
  int size = aotype.length(); 
  char char_array[size+1]; 

  strcpy(char_array, aotype.c_str()); 

  n = (int)char_array[0]-48;
  if (char_array[1]=='S')
    l = 0;
  else if (char_array[1]=='P')
    l = 1;
  else if (char_array[1]=='D')
    l = 2;
  else if (char_array[1]=='F')
    l = 3;
  else if (char_array[1]=='G')
    l = 4;
  else if (char_array[1]=='H')
    l = 5;
  else if (char_array[1]=='I')
    l = 6;
  else if (char_array[1]=='J')
    l = 7;

  if (char_array[1]=='D' && char_array[2]=='Z' && char_array[3]=='2')
    l = 22;
  else if (n<l+1)
  {
    printf(" ERROR: basis invalid. n: %i l: %i \n",n,l);
    exit(1);
  }

  return;
}

bool get_basis(string aname, int atnum, double Zeff, double* coords1, vector<vector<double> >& basis1, vector<vector<double> >& basis_aux1)
{
  int prl = 1;

  string bfilename = aname + ".basis";
  ifstream bfile;
  bfile.open(bfilename.c_str());
  if (!bfile)
  {
    return 0;
  }

  string line;
  int on_basis = 0;
  while (!bfile.eof())
  {  
    (bool)getline(bfile, line);
   //this split function is sensitive to number of spaces
    vector<string> tok_line = split1(line,' ');

    int size = tok_line.size();
    if (size>0)
    {
     //generic basis function centered at some XYZ
      vector<double> b1; 
      for (int m=0;m<10;m++) b1.push_back(0);
      b1[5] = coords1[0]; b1[6] = coords1[1]; b1[7] = coords1[2];
      b1[8] = Zeff; b1[9] = atnum;

      if (on_basis)
      {
        if (tok_line[0]=="END")
          on_basis = 0;
        else
        {
          double zeta = atof(tok_line[size-1].c_str());
          if (zeta<=0.) { printf("\n ERROR: zeta cannot be zero \n"); cout << " LINE: " << line << endl; exit(-1); }
          int n,l;
          get_n_l(tok_line[0],n,l);
          if (prl>1) printf("  found n/l/zeta: %i %i %5.2f (%s) \n",n,l,zeta,tok_line[size-1].c_str());

          if (on_basis==1)
          {
            if (l==0) add_s(basis1,b1,n,zeta);
            if (l==1) add_p(basis1,b1,n,zeta,3);
            if (l==2) add_d(basis1,b1,n,zeta,5);
            if (l==3) add_f(basis1,b1,n,zeta,7);
            if (l==4) add_g(basis1,b1,n,zeta,9);
            if (l==5) add_h(basis1,b1,n,zeta,11);
            //if (l==6) add_i(basis1,b1,n,zeta,13);

           //add just dz2 function
            if (l==22) add_dz2(basis1,b1,n,zeta);
          }
          else if (on_basis==2)
          {
            if (l==0) add_s(basis_aux1,b1,n,zeta);
            if (l==1) add_p(basis_aux1,b1,n,zeta,3);
            if (l==2) add_d(basis_aux1,b1,n,zeta,5);
            if (l==3) add_f(basis_aux1,b1,n,zeta,7);
            if (l==4) add_g(basis_aux1,b1,n,zeta,9);
            if (l==5) add_h(basis_aux1,b1,n,zeta,11);
            //if (l==6) add_i(basis_aux1,b1,n,zeta,13);
          }
        }
      }
      if (tok_line[0]=="BASIS")
      {
        if (prl>1) printf(" found BASIS \n");
        on_basis = 1;
      }
      else if (tok_line[0]=="FIT")
      {
        if (prl>1) printf(" found FIT \n");
        on_basis = 2;
      }
    }
  }


  bfile.close();

  return 1;
}

string get_aname(int Z)
{
  if (Z==0) return "X";
  if (Z==1) return "H";
  if (Z==2) return "He";
  if (Z==3) return "Li";
  if (Z==4) return "Be";
  if (Z==5) return "B";
  if (Z==6) return "C";
  if (Z==7) return "N";
  if (Z==8) return "O";
  if (Z==9) return "F";
  if (Z==10) return "Ne";
  if (Z==11) return "Na";
  if (Z==12) return "Mg";
  if (Z==13) return "Al";
  if (Z==14) return "Si";
  if (Z==15) return "P";
  if (Z==16) return "S";
  if (Z==17) return "Cl";
  if (Z==18) return "Ar";
  if (Z==19) return "K";
  if (Z==20) return "Ca";
  if (Z==21) return "Sc";
  if (Z==22) return "Ti";
  if (Z==23) return "V";
  if (Z==24) return "Cr";
  if (Z==25) return "Mn";
  if (Z==26) return "Fe";
  if (Z==27) return "Co";
  if (Z==28) return "Ni";
  if (Z==29) return "Cu";
  if (Z==30) return "Zn";
  if (Z==31) return "Ga";
  if (Z==32) return "Ge";
  if (Z==33) return "As";
  if (Z==34) return "Se";
  if (Z==35) return "Br";
  if (Z==36) return "Kr";
  return "n/a";
}

int get_Z(string atname)
{
  int Zeff = 0;
  if (atname=="H") Zeff = 1;
  if (atname=="He") Zeff = 2;
  if (atname=="Li") Zeff = 3;
  if (atname=="Be") Zeff = 4;
  if (atname=="B") Zeff = 5;
  if (atname=="C") Zeff = 6;
  if (atname=="N") Zeff = 7;
  if (atname=="O") Zeff = 8;
  if (atname=="F") Zeff = 9;
  if (atname=="Ne") Zeff = 10;
  if (atname=="Na") Zeff = 11;
  if (atname=="Mg") Zeff = 12;
  if (atname=="Al") Zeff = 13;
  if (atname=="Si") Zeff = 14;
  if (atname=="P") Zeff = 15;
  if (atname=="S") Zeff = 16;
  if (atname=="Cl") Zeff = 17;
  if (atname=="Ar") Zeff = 18;
  if (atname=="K") Zeff = 19;
  if (atname=="Ca") Zeff = 20;
  if (atname=="Sc") Zeff = 21;
  if (atname=="Ti") Zeff = 22;
  if (atname=="V") Zeff = 23;
  if (atname=="Cr") Zeff = 24;
  if (atname=="Mn") Zeff = 25;
  if (atname=="Fe") Zeff = 26;
  if (atname=="Co") Zeff = 27;
  if (atname=="Ni") Zeff = 28;
  if (atname=="Cu") Zeff = 29;
  if (atname=="Zn") Zeff = 30;
  if (atname=="Kr") Zeff = 36;
  return Zeff;
}

int is_valid_atom(string atname)
{
  if (get_Z(atname)>0)
    return 1;
  return 0;
}

bool get_secondary_basis(string name, int natoms, int* atno, double* coords, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, int prl)
{
  if (prl>0)
    printf(" reading secondary basis set \n");
  //vector<vector<double> > basis_aux;

  double coords1[3];
  for (int i=0;i<natoms;i++)
  {
    double atno1 = atno[i];
    string aname = get_aname(atno1)+name;
    coords1[0] = coords[3*i+0]; 
    coords1[1] = coords[3*i+1]; 
    coords1[2] = coords[3*i+2]; 

    bool found = get_basis(aname,i,atno1,coords1,basis,basis_aux);
    if (!found)
    {
      printf(" no secondary basis found for element %s \n",aname.c_str());
      return 0;
    }
  }

  //printf("  found %2i basis functions and %2i auxiliary functions \n",basis.size(),basis_aux.size());
  print_basis(natoms,basis,basis_aux,prl);
 
  return 1;
}


int read_input(string filename, bool gbasis, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, int& charge, int& nup, int* atno, double* &coords)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    printf("  couldn't open file: %s \n",filename.c_str());
    return 0;
  }

  string line;

  int natoms = 0;
  int geom_in_ang = 0;
  charge = 0;
  while (!infile.eof())
  {  
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');

    //cout << " 1READ: " << line << endl;
    if (tok_line.size()>0 && tok_line[0]=="Angstrom")
      geom_in_ang = 1;
    else if (tok_line.size()==1)
      charge = atoi(tok_line[0].c_str());
    else if (tok_line.size()==2)
    {
      charge = atoi(tok_line[0].c_str());
      nup = atoi(tok_line[1].c_str());
     //interpreting these as multiplicities
      if (nup==1) nup = 0;
      else if (nup==3) nup = 1;
      else if (nup==5) nup = 2;
      else if (nup==2 || nup==4) 
      {
        printf(" ERROR: doublet/quartet not available \n");
        exit(1);
      }
      else if (nup>5)
      {
        printf(" ERROR: spin>5 not available \n");
        exit(1);
      }
    }
    if (tok_line.size()>3) natoms++;

    //printf("   size: %2i \n",tok_line.size());
  }
  infile.clear(); infile.seekg(0);

  coords = new double[3*natoms]();
  int wa = 0;
  while (!infile.eof())
  {  
    (bool)getline(infile, line);
    vector<string> tok_line = split1(line,' ');

    //cout << " 2READ: " << line << endl;
    if (tok_line.size()>3)
    {
      string tkl0 = tok_line[0];
      if (tkl0=="X")
      {
       //dummy atom with (an integer) charge
        double Zeff = 0;
        if (tok_line.size()>4)
          Zeff = atof(tok_line[4].c_str());
        atno[wa] = 0;
        //atno[wa] = Zeff; //was this
        coords[3*wa+0] = atof(tok_line[1].c_str());
        coords[3*wa+1] = atof(tok_line[2].c_str());
        coords[3*wa+2] = atof(tok_line[3].c_str());
        if (geom_in_ang) { coords[3*wa+0] *= A2B;  coords[3*wa+1] *= A2B; coords[3*wa+2] *= A2B; }
        if (!gbasis && Zeff==0) get_basis(tkl0,wa,Zeff,&coords[3*wa],basis,basis_aux);
        wa++;
      }
      else if (is_valid_atom(tkl0))
      {
        double Zeff = get_Z(tkl0);
        atno[wa] = Zeff;
        if (tok_line.size()>4)
          Zeff = atof(tok_line[4].c_str());
        if (Zeff==0.) Zeff = atno[wa];

        coords[3*wa+0] = atof(tok_line[1].c_str());
        coords[3*wa+1] = atof(tok_line[2].c_str());
        coords[3*wa+2] = atof(tok_line[3].c_str());
        if (geom_in_ang) { coords[3*wa+0] *= A2B;  coords[3*wa+1] *= A2B; coords[3*wa+2] *= A2B; }
        //printf(" coords[%i]: %8.5f %8.5f %8.5f \n",wa,coords[3*wa+0],coords[3*wa+1],coords[3*wa+2]);
        if (!gbasis)
        {
          bool found = get_basis(tkl0,wa,Zeff,&coords[3*wa],basis,basis_aux);
          if (!found)
          {
            printf("\n ERROR: could not find basis for %s \n",tkl0.c_str());
            exit(-1);
          }
        }
        wa++;
      }
    }
  }

 #if GEOM_DEBUG
  {
    double coordn[3*natoms];
    double rot[9];
    double PI = 3.141592653;

    for (int j=0;j<9;j++) rot[j] = 0.;
    double angle_1 = randomf(0,PI);
   //xz
    rot[0] = cos(angle_1);
    rot[2] = -sin(angle_1);
    rot[4] = 1.;
    rot[6] = sin(angle_1);
    rot[8] = cos(angle_1);

    for (int i=0;i<natoms;i++)
    for (int j=0;j<3;j++)
    {
      double tr = 0.;
      for (int k=0;k<3;k++)
        tr += coords[i*natoms+k] * rot[k*3+j];
      coordn[i*natoms+j] = tr;
    }

    
    for (int j=0;j<9;j++) rot[j] = 0.;
    double angle_2 = randomf(0,PI);
   //xy
    rot[0] = cos(angle_2);
    rot[1] = -sin(angle_2);
    rot[3] = sin(angle_2);
    rot[4] = cos(angle_2);
    rot[8] = 1.;

    for (int i=0;i<natoms;i++)
    for (int j=0;j<3;j++)
    {
      double tr = 0.;
      for (int k=0;k<3;k++)
        tr += coordn[i*natoms+k] * rot[k*3+j];
      coords[i*natoms+j] = tr;
    }

    for (int n=0;n<basis.size();n++)
    {
      int n1 = basis[n][9];
      basis[n][5] = coords[3*n1+0];
      basis[n][6] = coords[3*n1+1];
      basis[n][7] = coords[3*n1+2];
    }
  }

 #endif

  return natoms;
}

void print_basis(int natoms, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, int prl)
{
  if (prl<0) return;

  int N = basis.size();
  printf("  found %2i atoms and %2i basis functions \n",natoms,N);
  //if (prl>0)
  if (prl>-1)
  {
    for (int n=0;n<N;n++)
      printf("   basis(%3i):  %i %i %2i  %6.3f  norm: %7.2f  (XYZ: %8.5f %8.5f %8.5f Z: %3.2f) \n",n,(int)basis[n][0],(int)basis[n][1],(int)basis[n][2],basis[n][3],basis[n][4],basis[n][5],basis[n][6],basis[n][7],basis[n][8]);
  }

  int Naux = basis_aux.size();
  printf("  found %2i auxiliary basis functions \n",Naux);
#if DDEBUG || FDEBUG || GDEBUG
  if (prl>-1 && Naux<250)
#else
  if (prl>1)
#endif
  {
    for (int n=0;n<Naux;n++)
      printf("   basis_aux(%3i):  %i %i %2i  %6.3f  norm: %7.2f  (XYZ: %8.5f %8.5f %8.5f) \n",n,(int)basis_aux[n][0],(int)basis_aux[n][1],(int)basis_aux[n][2],basis_aux[n][3],basis_aux[n][4],basis_aux[n][5],basis_aux[n][6],basis_aux[n][7]);
  }

}

bool check_valid_basis(int n1, int l1)
{
 //based on currently programmed AO types
  if (n1>12) return 0;
  if (l1==1 && n1>9) return 0;
  if (l1==2 && n1>7) return 0;
  if (l1==3 && n1>6) return 0;
  if (l1==4 && n1>5) return 0;
  if (l1==5 && n1>6) return 0;
  if (l1>=6) return 0;

  return 1;
}

int get_n12(int n1, int n2, int l1, int l2)
{
  int nr1 = n1-1;
  int nr2 = n2-1;
  //int nr1 = n1-l1-1;
  //int nr2 = n2-l2-1;
  return 1+nr1+nr2;
}

void screen_basis_aux(vector<vector<double> >& basis_aux, int prl)
{
  vector<vector<double> > basis1;
  for (int j=0;j<basis_aux.size();j++)
  {
    vector<double> ao1;
   //atom # first for sorting purposes
    ao1.push_back(basis_aux[j][9]);
    for (int k=0;k<9;k++)
      ao1.push_back(basis_aux[j][k]);

    basis1.push_back(ao1);
  }

  sort(basis1.begin(),basis1.end());
  basis1.erase(unique(basis1.begin(),basis1.end()), basis1.end());

  for (int j=0;j<basis1.size();j++)
  {
    vector<double> ao1;
   //restore regular order
    for (int k=0;k<9;k++)
      ao1.push_back(basis1[j][k+1]);
    ao1.push_back(basis1[j][0]);
    basis1[j] = ao1;
  }

  int N = basis1.size();

  int Nkeep = N;
  bool keep[N];
  for (int i1=0;i1<N;i1++) keep[i1] = 1;

  const float thresh = 0.05; //% minimum difference
  for (int i1=0;i1<N-1;i1++)
  {
    vector<double> ao1 = basis1[i1];
    vector<double> ao2 = basis1[i1+1];

    int at1 = ao1[9]; int at2 = ao2[9];
    if (at1==at2)
    {
      int n1 = ao1[0]; int n2 = ao2[0];
      int l1 = ao1[1]; int l2 = ao2[1];
      int m1 = ao1[2]; int m2 = ao2[2];

      if (n1==n2 && l1==l2 && m1==m2)
      {
        float zt1 = ao1[3]; float zt2 = ao2[3];
        if (fabs(zt2/zt1-1.f) < thresh)
        {
          if (prl>2)
            printf("   screening aux: %i %i %2i  zt: %8.5f %8.5f \n",n1,l1,m1,zt1,zt2);
          keep[i1] = 0;
          Nkeep--;
        }
      }
    }
  }

  if (prl>0)
    printf("   screen basis aux. N: %2i --> %2i \n",basis_aux.size(),Nkeep);

  basis_aux.clear();
  for (int j=0;j<N;j++)
  if (keep[j])
    basis_aux.push_back(basis1[j]);

  return;
}

#if 0
int get_bsize(int n1, int l1, int n2, int l2, int bs1, int bs2)
{
  int bsize = 1;
  if (l1==0 && l2==0) //s ftns
  {
    bsize = 8;
  }
  if (l1+l2==1) //p ftns
    bsize = 6;
  if (l2==2)
  {
    if (bs12
    bsize = 1;
    bsize = 4;
  }
}
#endif

int get_bsize(double ztmin, double ztmax)
{
  //double dz = ztmax-ztmin;
  double dd = ztmax/ztmin;

 //wider zeta range, more ftns
  if (dd>4.8) return 6;
  if (dd>3.8) return 5;
  if (dd>2.8) return 4;
  if (dd>1.8) return 3;
  if (dd>1.) return 2;

 //dd==1.
  return 1;
}

void create_basis_aux_v3(int natoms, vector<vector<double> >& basis_std, vector<vector<double> >& basis_aux)
{
  if (basis_std.size()<1) { printf("\n ERROR: couldn't create auxiliary basis set \n"); return; }

  vector<vector<double> > basis_max;
  vector<vector<double> > basis_min;
  vector<int> basis_size;

 //gather basis set, skipping m sequence
 //taking only the min and max values of zeta
  for (int a=0;a<natoms;a++)
  for (int n=1;n<10;n++)
  for (int l=0;l<6;l++)
  {
    double ztmin = 1000.; int minz = -1;
    double ztmax = 0.;    int maxz = -1;
    int bsize = 0;

   //find basis ftns with specific nl(m=0)
    for (int i1=0;i1<basis_std.size();i1++)
    if (basis_std[i1][9]==a && basis_std[i1][0]==n && basis_std[i1][1]==l && basis_std[i1][2]==0)
    {
      double zt1 = basis_std[i1][3];
      if (zt1<ztmin) { ztmin = zt1; minz = i1; }
      if (zt1>ztmax) { ztmax = zt1; maxz = i1; }
      bsize++;
    }
    if (minz!=-1 && maxz!=-1 && l!=5)
    {
      basis_min.push_back(basis_std[minz]);
      basis_max.push_back(basis_std[maxz]);
      basis_size.push_back(bsize);
    }
    else if (minz!=-1 && l==5)
    {
      printf("\n WARNING: H functions in standard basis not supported by RI==3 \n");
      exit(-1);
    }
  }

 //using the truncated, no m basis
  int N = basis_max.size();

  int lmax = 0;
  for (int n=0;n<natoms;n++)
  for (int i1=0;i1<N;i1++)
  if (basis_max[i1][9]==n)
  {
    vector<double> basis1 = basis_min[i1];
    vector<double> basis2 = basis_max[i1];
    int bs1 = basis_size[i1];

    int n1 = basis1[0]; int l1 = basis1[1]; double zt1 = basis1[3];
    int n2 = basis2[0]; int l2 = basis2[1]; double zt2 = basis2[3];

    for (int i2=0;i2<=i1;i2++)
    if (basis_max[i2][9]==n)
    {
      vector<double> basis3 = basis_min[i2];
      vector<double> basis4 = basis_max[i2];
      int bs2 = basis_size[i2];

      int n3 = basis3[0]; int l3 = basis3[1]; double zt3 = basis3[3];
      int n4 = basis4[0]; int l4 = basis4[1]; double zt4 = basis4[3];

      int n13 = get_n12(n1,n3,l1,l3);
      int l13 = l1+l3;
      double zt13 = zt1+zt3;
      double zt24 = zt2+zt4;

      int n24 = get_n12(n2,n4,l2,l4); 
      int l24 = l2+l4;

      int nmax = get_bsize(zt13,zt24);
      if (l24>=4 && nmax>2) nmax--; //trimming

      double zf = zt24/zt13;
      double B = 1;
      if (nmax>1) B = exp(log(zf)/(nmax-1));
      for (int ns=0;ns<nmax;ns++)
      {
        double ztn = zt13*pow(B,ns);
        vector<double> b1 = basis1;
        b1[0] = n13; b1[1] = l13; b1[3] = ztn;

        bool valid_basis = check_valid_basis(n13,l13);
        if (valid_basis)
        {
          //if (n==0)
          //  printf("  nl: %i %i / %i %i  zt: %8.5f \n",n13,l13,n24,l24,ztn);
          for (int m=-l13;m<=l13;m++)
          {
            b1[2] = m;
            b1[4] = norm(n13,l13,m,ztn);
            basis_aux.push_back(b1);
          }
        }

        if (valid_basis && l13>lmax) lmax = l13;

      } //adding ftns

    } //pairs of ftns on same atom
  } //loop n

  //screen_basis_aux(basis_aux,1);
  printf("  lmax in auxiliary basis: %i \n",lmax);

  if (basis_aux.size()<1) printf(" WARNING: didn't find any auxiliary basis ftns \n");

  return;
}

void create_basis_aux(int natoms, vector<vector<double> >& basis_std, vector<vector<double> >& basis_aux)
{
  if (basis_std.size()<1) { printf("\n ERROR: couldn't create auxiliary basis set \n"); return; }

  vector<vector<double> > basis;

 //gather basis set, skipping m sequence
  int natp = basis_std[0][9]; int np = basis_std[0][0]; int lp = basis_std[0][1]; double zetap = basis_std[0][3];
  basis.push_back(basis_std[0]);
  for (int i1=1;i1<basis_std.size();i1++)
  {
    vector<double> basis1 = basis_std[i1]; 
    int nat1 = basis1[9]; int n1 = basis1[0]; int l1 = basis1[1]; double zeta1 = basis1[3];
 
    if (fabs(zeta1-zetap)>0.01 || n1!=np || l1!=lp || nat1!=natp)
    {
      zetap = zeta1; natp = nat1; np = n1; lp = l1;
      basis.push_back(basis1);
    }
  }

 //using the truncated, no m basis
  int N = basis.size();

  int lmax = 0;
  for (int n=0;n<natoms;n++)
  {
    for (int i1=0;i1<N;i1++)
    {
      vector<double> basis1 = basis[i1];

      if (basis1[9]==n)
      {
        vector<double> b1 = basis1; int n1 = basis1[0]; int l1 = basis1[1]; double zeta1 = basis1[3];
        for (int i2=0;i2<=i1;i2++)
        {
          vector<double> basis2 = basis[i2]; int n2 = basis2[0]; int l2 = basis2[1]; double zeta2 = basis2[3];

          if (basis2[9]==n) //both basis ftns on same atom
          {
            int n12 = get_n12(n1,n2,l1,l2);
            int l12 = l1+l2;
            double z12 = zeta1+zeta2;
            b1[0] = n12;
            b1[1] = l12;
            b1[3] = z12;


            bool valid_basis = check_valid_basis(n12,l12);
            if (valid_basis)
            {
              for (int m=-l12;m<=l12;m++)
              {
                b1[2] = m;
                b1[4] = norm(n12,l12,m,z12);
                basis_aux.push_back(b1);
              }
            }

            if (valid_basis && l12>lmax) lmax = l12;

          } //basis on same atom
        } //pairs of ftns
      } //basis on atom n
    } //loop i1
  } //loop n

  screen_basis_aux(basis_aux,1);
  printf("  lmax in auxiliary basis: %i \n",lmax);

  if (basis_aux.size()<1) printf(" WARNING: didn't find any auxiliary basis ftns \n");

  return;
}

//puts first atom at 0.,0.,0.
void set_zero_coords_basis(int natoms, double* coords, vector<vector<double> >& basis, vector<vector<double> >& basis_aux)
{
  double A1 = coords[0]; double B1 = coords[1]; double C1 = coords[2];
  if (A1==0. && B1==0. && C1==0.)
    return;

  for (int n=0;n<natoms;n++)
  {
    coords[3*n+0] -= A1;
    coords[3*n+1] -= B1;
    coords[3*n+2] -= C1;
  }

  int N = basis.size();
  for (int j=0;j<N;j++)
  {
    basis[j][5] -= A1;
    basis[j][6] -= B1;
    basis[j][7] -= C1;
  }

  int Naux = basis_aux.size();
  for (int j=0;j<Naux;j++)
  {
    basis_aux[j][5] -= A1;
    basis_aux[j][6] -= B1;
    basis_aux[j][7] -= C1;
  }

  return;
}

int initialize(bool gbasis, vector<vector<double> >& basis, vector<vector<double> >& basis_aux, int* atno, double* &coords, int& charge, int& unpaired, double& Enn, int prl)
{
  //int prl = 1;

  //int charge = 0;
  //double* coords;
  //int* atno = new int[100]();

  //vector<vector<double> > basis;
  //vector<vector<double> > basis_aux;

  string geomfile = "GEOM";
  int natoms = read_input(geomfile,gbasis,basis,basis_aux,charge,unpaired,atno,coords);
  //void set_zero_coords_basis(int natoms, double* coords, vector<vector<double> >& basis, vector<vector<double> >& basis_aux)

  bool do_ps_integrals = read_int("PS");
  if (!do_ps_integrals) do_ps_integrals = read_int("QUAD");
  if (do_ps_integrals)
    set_zero_coords_basis(natoms,coords,basis,basis_aux);

  if (natoms<1)
  {
    printf("\n  no geometry found \n"); exit(1);
  }

  if (prl>-1)
  {
    printf("  Geometry: \n");
    for (int m=0;m<natoms;m++)
      printf("  %2i %8.5f %8.5f %8.5f \n",atno[m],coords[3*m+0],coords[3*m+1],coords[3*m+2]);
  }
  Enn = nuclear_repulsion(natoms,atno,coords);
  if (prl>0) printf("  nuclear repulsion: %12.8f \n",Enn);


  if (gbasis) return natoms;

  int auto_ri = read_ri();

  if (auto_ri>=2)
  {
    basis_aux.clear();
    if (auto_ri==3)
      create_basis_aux_v3(natoms,basis,basis_aux);
    else
      create_basis_aux(natoms,basis,basis_aux);
    prl++;
  }

  print_basis(natoms,basis,basis_aux,prl);

  return natoms;
}

