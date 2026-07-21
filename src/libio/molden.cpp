#include "write.h"
#include <sstream>

bool close_val(double v1, double v2);

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

void write_molden_mo_g(ofstream& outfile, vector<vector<double> > &basis, double* jCA, int No, double* eig)
{
  int N = basis.size();

  outfile << "[MO]" << endl;
  for (int i=0;i<N;i++)
  {
    int occ = 0; if (i<No) occ = 1;

    outfile << "Sym=X" << endl;
    if (eig!=NULL)
      outfile << "Ene=" << eig[i] << endl;
    else
      outfile << "Ene=0." << endl;
    outfile << "Spin=Alpha" << endl;
    outfile << "Occup=" << occ << endl;

    for (int j=0;j<N;j++)
    {
      int l1 = basis[j][1];
      if (l1<2)
        outfile << " " << j+1 << "  " << jCA[j*N+i] << endl;
      else if (l1==2) //d orbs
      {
        int j1 = j+2; //0
        outfile << " " << j+1 << "  " << jCA[j1*N+i] << endl;
        j1 = j+3; //+1
        outfile << " " << j+2 << "  " << jCA[j1*N+i] << endl;
        j1 = j+1; //-1
        outfile << " " << j+3 << "  " << jCA[j1*N+i] << endl;
        j1 = j+4; //+2
        outfile << " " << j+4 << "  " << jCA[j1*N+i] << endl;
        j1 = j; //-2
        outfile << " " << j+5 << "  " << jCA[j1*N+i] << endl;

        j += 4; //jump ahead
      }
      else if (l1==3) //f orbs
      {
        int j1 = j+3; //0
        outfile << " " << j+1 << "  " << jCA[j1*N+i] << endl;
        j1 = j+4; //1
        outfile << " " << j+2 << "  " << jCA[j1*N+i] << endl;
        j1 = j+2; //-1
        outfile << " " << j+3 << "  " << jCA[j1*N+i] << endl;
        j1 = j+5; //2
        outfile << " " << j+4 << "  " << jCA[j1*N+i] << endl;
        j1 = j+1; //-2
        outfile << " " << j+5 << "  " << jCA[j1*N+i] << endl;
        j1 = j+6; //3
        outfile << " " << j+6 << "  " << jCA[j1*N+i] << endl;
        j1 = j; //-3
        outfile << " " << j+7 << "  " << jCA[j1*N+i] << endl;

        j += 6;
      }
      else if (l1==4) //g orbs
      {
        int j1 = j+4;
        outfile << " " << j+1 << "  " << jCA[j1*N+i] << endl;
        j1 = j+5;
        outfile << " " << j+2 << "  " << jCA[j1*N+i] << endl;
        j1 = j+3;
        outfile << " " << j+3 << "  " << jCA[j1*N+i] << endl;
        j1 = j+6;
        outfile << " " << j+4 << "  " << jCA[j1*N+i] << endl;
        j1 = j+2;
        outfile << " " << j+5 << "  " << jCA[j1*N+i] << endl;
        j1 = j+7;
        outfile << " " << j+6 << "  " << jCA[j1*N+i] << endl;
        j1 = j+1;
        outfile << " " << j+7 << "  " << jCA[j1*N+i] << endl;
        j1 = j+8;
        outfile << " " << j+8 << "  " << jCA[j1*N+i] << endl;
        j1 = j;
        outfile << " " << j+9 << "  " << jCA[j1*N+i] << endl;

        j += 8;
      }
      else //not printing other orbs at the moment
        outfile << " " << j+1 << "  0.0" << endl;
    }
  }

  int lmax = 0;
  for (int i=0;i<N;i++)
    lmax = max(lmax,(int)basis[i][1]);

  if (lmax==2)
    outfile << " [5D] " << endl;
  else if (lmax==3)
    outfile << " [5D7F] " << endl;
  else
    outfile << " [5D7F9G] " << endl;

  return;
}

void write_molden_g(int btype, int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, double* eig, string fname)
{
 //this ftn expects orbitals to be spherical, not Cartesian. May crash for Cartesians
  int N = basis.size();

  No = abs(No);

  const double n0 = norm_sh(0,0);
  double B2A = 1./A2B;

  string filename = fname;
  ofstream outfile;
  outfile.open(filename.c_str());
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

 //get molden formatted basis from basis vector
  for (int i=0;i<natoms;i++)
  {
    string aname = get_aname(atno[i]);

    int wb1 = -1;
    for (int j=0;j<basis.size();j++)
    if (basis[j][9]==i)
    { wb1 = j; break; }
    int wb2 = wb1;
    for (int j=wb1+1;j<basis.size();j++)
    if (basis[j][9]!=i)
    { wb2 = j; break; }
    if (i+1==natoms) wb2 = N;
    //printf(" wb12: %2i %2i \n",wb1,wb2);

    string b1 = "";
    for (int j=wb1;j<wb2;)
    {
      int l1 = (int)basis[j][1];
      int ng = (int)basis[j][3];
      int row0 = j + l1; // the m==0 row of this shell
      double nsh0 = norm_sh(l1,0);

      string lv = "S"; if (l1==1) lv = "P"; else if (l1==2) lv = "D"; else if (l1==3) lv = "F"; else if (l1==4) lv = "G"; else if (l1==5) lv = "H";
      b1 += lv+" "+to_string(ng)+" 1.0 \n";
      for (int k=0;k<ng;k++)
      {
        double g1 = basis[row0][10+k];
        double n1 = basis[row0][10+ng+k] * n0 / ( get_gto_norm(l1,g1) * nsh0 );
        b1 += " " + to_string(g1) + "  " + to_string(n1) + "\n";
      }

      j += (2*l1+1);
    }
    outfile << i+1 << "    0 " << endl;
    outfile << b1 << endl;
  }

  write_molden_mo_g(outfile,basis,jCA,No,eig);

  outfile.close();

  return;
}

void write_molden_g(int natoms, int* atno, double* coords, vector<vector<double> > &basis, double* jCA, int No, double* eig, string fname)
{
  return write_molden_g(0,natoms,atno,coords,basis,jCA,No,eig,fname);
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
