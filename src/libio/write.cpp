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

void write_molden(bool gbasis, int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, string fname)
{
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

