#define DDEBUG 0
#define FDEBUG 0
#define GDEBUG 0

#include "read.h"
#include "fp_def.h"

#include <sstream>
#define SSTRF( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << fixed << setprecision(14) << x ) ).str()
#define SSTRF2( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << fixed << setprecision(4) << x ) ).str()

FP2 norm(int n, int l, int m, FP2 zeta);

#define A2B 1.8897261


void print_coords(int natoms, FP1* coordsf)
{
  FP1 B2A = 1.f/A2B;
  for (int n=0;n<natoms;n++)
    printf(" %10.5f %10.5f %10.5f \n",coordsf[3*n+0]*B2A,coordsf[3*n+1]*B2A,coordsf[3*n+2]*B2A);
}

void print_gradient(int natoms, FP2* grad)
{
  for (int n=0;n<natoms;n++)
    printf(" %10.5f %10.5f %10.5f \n",grad[3*n+0],grad[3*n+1],grad[3*n+2]);
}

void print_square_diff(int N, FP2* S1, FP2* S2) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S1[n*N+m],S1[n*N+m]);
    printf("\n");
  }
}

void print_square_fine(int N, FP2* S) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %12.8f",S[n*N+m]);
    printf("\n");
  }
}

#if !EVL64
void print_square(int N, FP2* S) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}
#endif

void print_square(int N, FP1* S) 
{
  for (int n=0;n<N;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

#if !EVL64
void print_square(int M, int N, FP1* S) 
{
  for (int n=0;n<M;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}
#endif

void print_square(int M, int N, FP2* S) 
{
  for (int n=0;n<M;n++)
  {
    for (int m=0;m<N;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_ss(int N, FP2* S)
{
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %8.5f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_sm(int N, FP2* S)
{
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %5.2f",S[n*N+m]);
    printf("\n");
  }
}

#if !EVL64
void print_square_sm(int N, FP1* S)
{
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %5.2f",S[n*N+m]);
    printf("\n");
  }
}
#endif

void print_square_ss_sm(int N, FP2* S)
{ 
  for (int n=0;n<N;n++)
  {
    printf("  ");
    for (int m=0;m<N;m++)
      printf(" %6.3f",S[n*N+m]);
    printf("\n");
  }
}

void print_square_nxn(int No, int N, FP1* S)
{
  for (int n=0;n<No;n++)
  {
    for (int m=0;m<No;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}

#if !EVL64
void print_square_nxn(int No, int N, FP2* S)
{
  for (int n=0;n<No;n++)
  {
    for (int m=0;m<No;m++)
      printf(" %10.5f",S[n*N+m]);
    printf("\n");
  }
}
#endif

void print_rectangle_sm(int N1, int N2, FP2* S)
{
  for (int n=0;n<N1;n++)
  {
    printf("  ");
    for (int m=0;m<N2;m++)
      printf(" %5.2f",S[n*N2+m]);
    printf("\n");
  }
}


vector<string> split1(const string &s, char delim)
{
  stringstream ss(s);
  string item;  
  vector<string> tokens;
  while (getline(ss, item, delim))
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
    printf("   couldn't open file: %s \n",filename.c_str());
    return 0;
  }

  bool done = 0;
  string text = "";
  while (!infile.eof() && !done)
  {  
    string line;
    getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>1 && tok_line[0].compare(aname)==0)
    {
      //printf("  found elem %2s \n",aname.c_str());
      while (!infile.eof())
      {
        getline(infile, line);
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

FP2 read_FP1(string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0.;

  string line;
  bool success = (bool) getline(infile, line);

  FP2 val = 0.;
  if (success)
    val = atof(line.c_str());

  infile.close();

  return val;
}

void read_thresh(FP1& no_thresh, FP1& occ_thresh)
{
  string filename = "THRESH";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return;

  string line;
  bool success = (bool) getline(infile, line);

  if (success)
  {
    no_thresh = atof(line.c_str());

    success = (bool) getline(infile, line);
    if (success)
      occ_thresh = atof(line.c_str());
  }

  infile.close();

  return;
}

void read_eps(FP2& eps1, FP2& eps2, FP2& eps1s)
{
  string filename = "EPS";

  eps1 = 0.001;
  eps2 = 0.0000001;
  eps1s = 0.00005;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return;

  string line;
  bool success = (bool) getline(infile, line);

  if (success)
  {
    eps1 = atof(line.c_str());

    success = (bool) getline(infile, line);
    if (success)
      eps2 = atof(line.c_str());

    success = (bool) getline(infile, line);
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

void read_eps(FP2& eps1, FP2& eps2)
{
  string filename = "EPS";

  eps1 = 0.001;
  eps2 = 0.000001;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return;

  string line;
  bool success = (bool) getline(infile, line);

  if (success)
  {
    eps1 = atof(line.c_str());

    success = (bool) getline(infile, line);
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

bool read_dft(int& dt1, int& dt2)
{
  string filename = "DFT";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool) getline(infile, line);

  dt1 = dt2 = -1;
  if (success)
  {
   //exchange
    dt1 = atoi(line.c_str());

    success = (bool) getline(infile, line);
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
  bool success = (bool) getline(infile, line);

  int default_val = 1;
  int do_esci;
  if (success)
    do_esci = atoi(line.c_str());
  else
    do_esci = default_val;

  infile.close();

  return do_esci;
}

int read_ri()
{
  string filename = "RI";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool) getline(infile, line);

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
  bool success = (bool) getline(infile, line);

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
  bool success = (bool) getline(infile, line);

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

int read_cas()
{
  string filename = "CAS";

  int default_val = 0;

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return default_val;

  string line;
  bool success = (bool) getline(infile, line);

  int cas_type;
  if (success)
    cas_type = atoi(line.c_str());
  else
    cas_type = default_val;

  infile.close();

  return cas_type;
}

bool read_cas_act(int& N, int& M)
{
  string filename = "CAS_ACT";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool) getline(infile, line);

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

void read_rotate(int N, FP2* jCA)
{
  string filename = "ROT";

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return;


  int nrot = 0;
  while (!infile.eof())
  {
    string line;
    bool success = (bool) getline(infile, line);

    if (success)
    {
      FP2 vec1[N];
      FP2 vec2[N];
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

  return;
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
  bool success = (bool) getline(infile, line);

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
    bool success = (bool) getline(infile, line);

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
    bool success = (bool) getline(infile, line);

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
    bool success = (bool) getline(infile, line);

    if (success)
      is_restart = atoi(line.c_str());
    else
      is_restart = 1;

    //printf(" found restart: %i \n",is_restart);
  }
  infile.close();
 
  return is_restart;
}

void read_nrad_nang(int& nrad, int& nang, int type)
{
  string filename = "GRID";
  if (type==2) filename = "GRID2";

  ifstream infile;
  infile.open(filename.c_str());    
  if (!infile)
  {
    printf("   couldn't open GRID file. please provide GRID \n");
    exit(1);
  }
  
  string line;
  bool success = (bool) getline(infile, line);
  if (success)
    nrad = atoi(line.c_str());

  success = (bool) getline(infile, line);
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
  bool success = (bool) getline(infile, line);
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

int read_iarray(short type1, short type2, short i1, int s1, int s2, FP1* A)
{
  string filename = get_iarray_name(type1,type2,i1);
  //printf("  reading %8s \n",filename.c_str());

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    printf("    couldn't open file (%s) \n",filename.c_str());
    return 0;
  }

  string line;
  int found = 1;

  for (int i=0;i<s1;i++)
  {
    bool success = (bool) getline(infile, line);
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

#if !EVL64
int read_iarray(short type1, short type2, short i1, int s1, int s2, FP2* A)
{
  string filename = get_iarray_name(type1,type2,i1);
  //printf("  reading %8s \n",filename.c_str());

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    printf("    couldn't open file (%s) \n",filename.c_str());
    return 0;
  }

  string line;
  int found = 1;

  for (int i=0;i<s1;i++)
  {
    bool success = (bool) getline(infile, line);
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
#endif

int read_gridpts(int s1, int s2, FP1* A, string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
    return 0;

  string line;
  bool success = (bool) getline(infile, line);
  int wi = 0;
  while (wi<s1)
  {
    success = (bool) getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<s2) { printf(" ERROR: file size incorrect (read_gridpts) \n"); exit(1); }
      for (int m=0;m<s2;m++)
        A[wi*s2+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }
  
  infile.close();

  return 1;
}

int read_square(int N, FP2* Pao, string filename)
{

  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    //printf("   couldn't open file \n");
    return 0;
  }

  string line;
  bool success = (bool) getline(infile, line);
  int wi = 0;
  while (wi<N)
  {
    success = (bool) getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size not square \n"); exit(1); }
      for (int m=0;m<N;m++)
        Pao[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }
  
  infile.close();

  return 1;
}

int read_square(vector<vector<FP2> > basis, FP2* Pao, string filename)
{
  int N = basis.size();
  return read_square(N,Pao,filename);
}

void read_SENT(string dirname, int N, FP2* S, FP2* T, FP2* En, int prl)
{
  string filename = "SENT"; if (dirname!="0") filename = dirname+"/"+"SENT";
  printf("  attempting file read: %s \n",filename.c_str());

  ifstream infile;
  infile.open(filename.c_str());    
  if (!infile)
  {
    printf("   couldn't open file \n");
    return;
  }

  string line;  
  getline(infile, line);
  int wi = 0;
  while (wi<N)
  {
    getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size wrong in S \n"); exit(1); }
      for (int m=0;m<N;m++)
        S[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }

  getline(infile, line);
  wi = 0;
  while (wi<N)
  {
    getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size wrong in En \n"); exit(1); }
      for (int m=0;m<N;m++)
        En[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }

  getline(infile, line);
  wi = 0;
  while (wi<N)
  {
    getline(infile, line);
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

void read_integrals(string dirname, int N, int Naux, FP2* S, FP2* T, FP2* En, FP2* A, FP2* Ciap, int prl)
{
  int N2 = N*N;
  //int N3 = N2*N;
  //int na = Naux;
  //int Naux2 = Naux*Naux;

 //first file
  string filename = "Ciap"; if (dirname!="0") filename = dirname+"/"+"Ciap";
  printf("  attempting file read: %s \n",filename.c_str());

  ifstream infile;
  infile.open(filename.c_str());    
  if (!infile)
  {
    printf("   couldn't open file \n");
    return;
  }
  
  string line;

  getline(infile, line);
  int wi = 0;
  while (!infile.eof())
  {
    getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<Naux+1) { printf(" ERROR: file size wrong in Ciap (%2i vs %2i) \n",tok_line.size(),Naux); exit(1); }
      for (int m=0;m<Naux;m++)
        Ciap[wi*Naux+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }
  if (wi==N2) { if (prl>1) printf("   found all lines of Ciap \n"); }
  else printf(" Ciap missing lines \n");
  infile.close();


 //second file
  filename = "A"; if (dirname!="0") filename = dirname+"/"+"A";
  printf("  attempting file read: %s \n",filename.c_str());

  infile.open(filename.c_str());    
  if (!infile)
  {
    printf("   couldn't open file \n");
    return;
  }
  
  getline(infile, line);
  wi = 0;
  while (!infile.eof())
  {
    getline(infile, line);
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
    printf("   couldn't open file \n");
    return;
  }
  
  getline(infile, line);
  wi = 0;
  while (wi<N)
  {
    getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size wrong in S \n"); exit(1); }
      for (int m=0;m<N;m++)
        S[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }

  getline(infile, line);
  wi = 0;
  while (wi<N)
  {
    getline(infile, line);
    vector<string> tok_line = split1(line,' ');
    if (tok_line.size()>0)
    {
      if (tok_line.size()<N) { printf(" ERROR: file size wrong in En \n"); exit(1); }
      for (int m=0;m<N;m++)
        En[wi*N+m] = atof(tok_line[m+1].c_str());
      wi++;
    }
  }

  getline(infile, line);
  wi = 0;
  while (wi<N)
  {
    getline(infile, line);
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

#if 0
 //now compute g matrix
  FP2* Ai = new FP2[Naux2];
  for (int m=0;m<Naux2;m++) Ai[m] = A[m];
  printf(" cannot compute g, need invert function \n");
  //Invert_stable(Ai,Naux,inv_cutoff);

  FP2* AC = new FP2[Naux*N2]();
  int nanan2 = Naux2+N2;

#if 0
  FP2* tmp = new FP2[nomp*nanan2]();

  compute_AC(N,na,Ai,Ciap,AC,tmp);
  compute_g(N,Naux,Ciap,AC,tmp,g);
  delete [] tmp;
#endif

  delete [] Ai;
  delete [] AC;

  if (prl>2)
  {
    printf("\n printing gmnls \n");
    for (int m=0;m<N;m++)
    for (int n=0;n<N;n++)
    {
      printf("  g%i%i:\n",m,n);
      print_square_ss(N,&g[m*N3+n*N2]);
    }
  }
  printf("\n");
#endif

  return;
}

FP2 nuclear_repulsion(int natoms, int* atno, FP2* coords)
{
  FP2 Enn = 0;

  for (int n=0;n<natoms;n++)
  for (int m=0;m<n;m++)
  {
    FP2 x12 = coords[3*n+0] - coords[3*m+0];
    FP2 y12 = coords[3*n+1] - coords[3*m+1];
    FP2 z12 = coords[3*n+2] - coords[3*m+2];
    FP2 zz = atno[n]*atno[m];
    Enn += zz/sqrt(x12*x12+y12*y12+z12*z12);
  }

  return Enn;
}

#if !EVL64
FP2 nuclear_repulsion(int natoms, int* atno, FP1* coordsf)
{
  FP2 Enn = 0;

  for (int n=0;n<natoms;n++)
  for (int m=0;m<n;m++)
  {
    FP1 x12 = coordsf[3*n+0] - coordsf[3*m+0];
    FP1 y12 = coordsf[3*n+1] - coordsf[3*m+1];
    FP1 z12 = coordsf[3*n+2] - coordsf[3*m+2];
    FP1 zz = atno[n]*atno[m];
    Enn += zz/sqrt(x12*x12+y12*y12+z12*z12);
  }

  return Enn;
}
#endif

void add_s(vector<vector<FP2> > &basis1, vector<FP2> ao1, int n, FP2 zeta)
{
 //ns function
  ao1[0] = n; ao1[1] = 0; ao1[2] = 0;
  ao1[3] = zeta;
  ao1[4] = norm(n,0,0,zeta);
  basis1.push_back(ao1);

  return;
}

void add_p(vector<vector<FP2> > &basis1, vector<FP2> ao1, int n, FP2 zeta, int np)
{
 //2px
  ao1[0] = n; ao1[1] = 1; ao1[2] = 0;
  ao1[3] = zeta;
  ao1[4] = norm(n,1,0,zeta);
  basis1.push_back(ao1);

  if (np<2) return;

 //2py
  ao1[2] = 1; 
  basis1.push_back(ao1);
 //2pz
  ao1[2] = -1;
  basis1.push_back(ao1);

  return;
}

void add_d(vector<vector<FP2> > &basis1, vector<FP2> ao1, int n, FP2 zeta, int nd)
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

   #if DDEBUG
    //if (n==3 || (n>3 && m==0))
   #endif
    basis1.push_back(ao1);
  }
 #endif

  return;
}

void add_f(vector<vector<FP2> > &basis1, vector<FP2> ao1, int n, FP2 zeta, int nf)
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

#if FDEBUG
    if (m==2)
#endif
    {
      basis1.push_back(ao1);
      bf++;
    }

    if (bf>nf) break;
  }

  return;
}

void add_g(vector<vector<FP2> > &basis1, vector<FP2> ao1, int n, FP2 zeta, int ng)
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

#if GDEBUG
 //m==2 problem?
 //m==-1,0 problem?
  //  if (m!=2 && m!=-1 && m!=0)
   // if (m==2 || m==-1 || m==0)
    if (m<=0)
#endif
    {
      basis1.push_back(ao1);
    }
            
    bf++; if (bf>ng) break;
  }
       
  return;
}

void add_h(vector<vector<FP2> > &basis1, vector<FP2> ao1, int n, FP2 zeta, int nh)
{
 //6h
  ao1[0] = n; ao1[1] = 5;
  ao1[3] = zeta;
 
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

  if (n<l+1)
  {
    printf(" ERROR: basis invalid. n: %i l: %i \n",n,l);
    exit(1);
  }

  return;
}

void get_basis(string aname, int atnum, int Zeff, FP2* coords1, vector<vector<FP2> >& basis1, vector<vector<FP2> >& basis_aux1)
{
  int prl = 1;

  string bfilename = aname + ".basis";
  ifstream bfile;
  bfile.open(bfilename.c_str());
  if (!bfile)
  {
    printf(" ERROR: couldn't find %s \n",bfilename.c_str());
    exit(1);
  }

  string line;
  int on_basis = 0;
  while (!bfile.eof())
  {  
    getline(bfile, line);
   //this split function is sensitive to number of spaces
    vector<string> tok_line = split1(line,' ');

    int size = tok_line.size();
    if (size>0)
    {
     //generic basis function centered at some XYZ
      vector<FP2> b1; 
      for (int m=0;m<10;m++) b1.push_back(0);
      b1[5] = coords1[0]; b1[6] = coords1[1]; b1[7] = coords1[2];
      b1[8] = Zeff; b1[9] = atnum;

      if (on_basis)
      {
        if (tok_line[0]=="END")
          on_basis = 0;
        else
        {
          FP2 zeta = atof(tok_line[size-1].c_str());
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

}

string get_aname(int Z)
{
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

int read_input(string filename, bool gbasis, vector<vector<FP2> >& basis, vector<vector<FP2> >& basis_aux, int& charge, int& nup, int* atno, FP2* &coords)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
    printf("   couldn't open file: %s \n",filename.c_str());
    return 0;
  }

  string line;

  int natoms = 0;
  int geom_in_ang = 0;
  charge = 0;
  while (!infile.eof())
  {  
    getline(infile, line);
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
        printf(" ERROR: FP2t/quartet not available \n");
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

  coords = new FP2[3*natoms]();
  int wa = 0;
  while (!infile.eof())
  {  
    getline(infile, line);
    vector<string> tok_line = split1(line,' ');

    //cout << " 2READ: " << line << endl;
    if (tok_line.size()>3)
    {
      string tkl0 = tok_line[0];
      if (tkl0=="X")
      {
       //dummy atom with an integer charge
        int Zeff = atoi(tok_line[4].c_str());
        atno[wa] = Zeff;
        coords[3*wa+0] = atof(tok_line[1].c_str());
        coords[3*wa+1] = atof(tok_line[2].c_str());
        coords[3*wa+2] = atof(tok_line[3].c_str());
        wa++;
      }
      else if (is_valid_atom(tkl0))
      {
        int Zeff = get_Z(tkl0);
        atno[wa] = Zeff;
        coords[3*wa+0] = atof(tok_line[1].c_str());
        coords[3*wa+1] = atof(tok_line[2].c_str());
        coords[3*wa+2] = atof(tok_line[3].c_str());
        if (geom_in_ang) { coords[3*wa+0] *= A2B;  coords[3*wa+1] *= A2B; coords[3*wa+2] *= A2B; }
        //printf(" coords[%i]: %8.5f %8.5f %8.5f \n",wa,coords[3*wa+0],coords[3*wa+1],coords[3*wa+2]);
        if (!gbasis) get_basis(tkl0,wa,Zeff,&coords[3*wa],basis,basis_aux);
        wa++;
      }
    }
  }

  return natoms;
}

int initialize(bool gbasis, vector<vector<FP2> >& basis, vector<vector<FP2> >& basis_aux, int* atno, FP2* &coords, int& charge, int& unpaired, FP2& Enn, int prl)
{
  //int prl = 1;

  //int charge = 0;
  //FP2* coords;
  //int* atno = new int[100]();

  //vector<vector<FP2> > basis;
  //vector<vector<FP2> > basis_aux;

  string geomfile = "GEOM";
  int natoms = read_input(geomfile,gbasis,basis,basis_aux,charge,unpaired,atno,coords);

  if (prl>0)
  {
    printf("  Geometry: \n");
    for (int m=0;m<natoms;m++)
      printf("  %2i %8.5f %8.5f %8.5f \n",atno[m],coords[3*m+0],coords[3*m+1],coords[3*m+2]);
  }
  Enn = nuclear_repulsion(natoms,atno,coords);
  if (prl>0) printf("  nuclear repulsion: %12.8f \n",Enn);


  if (gbasis) return natoms;

  int N = basis.size();
  printf("  found %2i atoms and %2i basis functions \n",natoms,N);
  prl = 1;
  if (prl>0)
  {
    for (int n=0;n<N;n++)
      printf("   basis(%3i):  %i %i %2i  %6.3f  norm: %7.2f  (XYZ: %8.5f %8.5f %8.5f) \n",n,(int)basis[n][0],(int)basis[n][1],(int)basis[n][2],basis[n][3],basis[n][4],basis[n][5],basis[n][6],basis[n][7]);
  }

  int Naux = basis_aux.size();
  printf("  found %2i auxiliary basis functions \n",Naux);
#if DDEBUG || FDEBUG || GDEBUG
  if (prl>0 && Naux<250)
#else
  if ((prl>0 && Naux<20) || prl>1)
#endif
  {
    for (int n=0;n<Naux;n++)
      printf("   basis_aux(%3i):  %i %i %2i  %6.3f  norm: %7.2f  (XYZ: %8.5f %8.5f %8.5f) \n",n,(int)basis_aux[n][0],(int)basis_aux[n][1],(int)basis_aux[n][2],basis_aux[n][3],basis_aux[n][4],basis_aux[n][5],basis_aux[n][6],basis_aux[n][7]);
  }

  return natoms;
}
