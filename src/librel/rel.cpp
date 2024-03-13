#include <cmath>
#include <fstream>
#include "fp_def.h"
#include "utils.h"
#include "rel.h"

using namespace std;

//this function is based on the pyscf x2c code to find hcore
void compute_hcore(int N, FP2* t, FP2* v, FP2* pVp, FP2* s, FP2 c, FP2* hcore, int prl)
{
  /*
    ofstream fileS;
    fileS.open("S");
    for(int i = 0; i < N; i++) 
    {    
      for(int j = 0; j < N; j++) 
      {    
        fileS << s[i*N+j] << " ";
      }    
      fileS << "\n";
    }    
    fileS.close();

    ofstream fileT;
    fileT.open("T");
    for(int i = 0; i < N; i++) 
    {    
      for(int j = 0; j < N; j++) 
      {    
        fileT << t[i*N+j] << " "; 
      }    
      fileT << "\n";
    }    
    fileT.close();

    ofstream fileV;
    fileV.open("V");
    for(int i = 0; i < N; i++) 
    {    
      for(int j = 0; j < N; j++) 
      {    
        fileV << v[i*N+j] << " "; 
      }    
      fileV << "\n";
    }    
    fileV.close();
  */

    ofstream filepVp;
    filepVp.open("pVp");
    for(int i = 0; i < N; i++) 
    {    
      for(int j = 0; j < N; j++) 
      {    
        filepVp << pVp[i*N+j] << " "; 
      }    
      filepVp << "\n";
    }    
    filepVp.close();
  

    int idx1 = 0;
    int idx2 = 0;
    int idx3 = 0;
    int Nsq = N * N;
    int Ntim2 = N * 2;
    int N2 = Ntim2 * Ntim2;
    int N2over2 = N2/2;
    FP2* h = new FP2[N2];
    FP2* m = new FP2[N2];
    
    for(int i = 0; i < N2; i++)
    {
      h[i] = 0.;
    }
    for(int i = 0; i < N2; i++)
    {
      m[i] = 0.;
    }  
  
    //acc_assign(N2,h,0.);
    //acc_assign(N2,m,0.);
    
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < N; j++)
      {
        h[Ntim2*i+j] = v[N*i+j];
        h[Ntim2*i+j+N] = t[N*i+j];
        h[Ntim2*i+j+N2over2] = t[N*i+j];
        h[Ntim2*i+j+N2over2+N] = pVp[N*i+j]/(c * c * 4.) - t[N*i+j];
        m[Ntim2*i+j] = s[N*i+j];
        m[Ntim2*i+j+N2over2+N] = t[N*i+j]/(c * c * 2.);
      }
    }

    if(prl>0)
    {
      printf("S: \n");
      for(int i = 0; i < N; i++) 
      {    
        for(int j = 0; j < N; j++) 
        {    
          printf("%6.4f ", s[i*N+j]);
        }    
        printf("\n");
      }    
      printf("\n");

      printf("hmod: \n");
      for(int i = 0; i < Ntim2; i++)
      {
        for(int j = 0; j < Ntim2; j++)
        {
          printf("%6.4f ", h[i*Ntim2+j]);
        }
        printf("\n");
      }
      printf("\n");

      printf("m before diag: \n");
      for(int i = 0; i < Ntim2; i++) 
      {    
        for(int j = 0; j < Ntim2; j++) 
        {    
          printf("%10.10f ", m[i*Ntim2+j]);
        }    
        printf("\n");
      }
      printf("\n"); 
    }

    FP2* eigvalm = new FP2[Ntim2];
    
    for(int i = 0; i < Ntim2; i++)
    {
      eigvalm[i] = 0.;
    }
    
    //acc_assign(Ntim2,eigvalm,0.);

    DiagonalizeP(m,eigvalm,Ntim2);
    for(int i = 0; i < Ntim2; i++)
    {
      if(eigvalm[i] > 1e-9)                               //Should be the case that the eigenvalues are sorted from largest to smallest but make sure
      idx1++;
    }

    if(prl>0)
    {
      printf("m after diag: \n");
      for(int i = 0; i < Ntim2; i++) 
      {    
        for(int j = 0; j < Ntim2; j++) 
        {    
          printf("%10.10f ", m[i*Ntim2+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }
    
    FP2* mconj = new FP2[N2];

    for(int i = 0; i < N2; i++)
    {
      mconj[i] = 0.;
    }

    for(int i = 0; i < Ntim2; i++)
    {
      for(int j = 0; j < Ntim2; j++)
      {
        mconj[i*Ntim2+j] = m[j*Ntim2+i];
      }
    }

    if(prl>0)
    {
      printf("mconj: \n");
      for(int i = 0; i < Ntim2; i++) 
      {    
        for(int j = 0; j < Ntim2; j++) 
        {    
          printf("%10.10f ", mconj[i*Ntim2+j]);
        }    
        printf("\n");
      }    
      printf("\n");

      printf("eigvalm after: \n");
      for(int i = 0; i < Ntim2; i++) 
      {
        printf("%10.10f ", eigvalm[i]); 
      }
      printf("\n");
      printf("\n");
    }

    int idx1dim = Ntim2 * idx1;
    FP2* mnxi = new FP2[idx1dim];
    
    for(int i = 0; i < idx1dim; i++)
    {
      mnxi[i] = 0.;
    }

    //acc_assign(idx1dim,mnxi,0.);

    for(int i = Ntim2 - idx1; i < Ntim2; i++)
    {
      for(int j = 0; j < Ntim2; j++)
      {
        mnxi[idx1*j+i-Ntim2+idx1] = mconj[Ntim2*j+i] / sqrt(eigvalm[i]);
      }
    }
    
    if(prl>0)
    {
      printf("mnxi: \n");
      for(int i = 0; i < Ntim2; i++) 
      {    
        for(int j = 0; j < idx1; j++) 
        {    
          printf("%6.4f ", mnxi[i*idx1+j]);
        }    
        printf("\n");
      }
      printf("\n");
    }

    FP2* eigvecconjm = new FP2[idx1dim];
    
    for(int i = 0; i < idx1dim; i++)
    {
      eigvecconjm[i] = 0.;
    }

    //acc_assign(idx1dim,eigvecconjm,0.);

    for(int i = 0; i < idx1; i++)
    {
      for(int j = 0; j < Ntim2; j++)
      {
        eigvecconjm[Ntim2*i+j] = mnxi[idx1*j+i];
      }
    }

    if(prl>0)
    {
      printf("eigvecconjm: \n");
      for(int i = 0; i < idx1; i++) 
      {    
        for(int j = 0; j < Ntim2; j++) 
        {    
          printf("%6.4f ", eigvecconjm[i*Ntim2+j]);
        }    
        printf("\n");
      }
      printf("\n");
    }

    FP2* mconjh = new FP2[idx1dim];
    
    for(int i = 0; i < idx1dim; i++)
    {
      mconjh[i] = 0.;
    }

    //acc_assign(idx1dim,mconjh,0.);

    for(int i = 0; i < idx1; i++)
    {
      for(int j = 0; j < Ntim2; j++)
      {
        for(int k = 0; k < Ntim2; k++)
        {
          mconjh[Ntim2*i+j] += eigvecconjm[Ntim2*i+k] * h[Ntim2*k+j];
        }
      }
    }

    if(prl>0)
    {
      printf("mconjh: \n");
      for(int i = 0; i < idx1; i++) 
      {    
        for(int j = 0; j < Ntim2; j++) 
        {    
          printf("%6.4f ", mconjh[i*Ntim2+j]);
        }    
        printf("\n");
      }
      printf("\n");
    }

    int idx1sq = idx1 * idx1;
    FP2* mhm = new FP2[idx1sq];
    
    for(int i = 0; i < idx1sq; i++)
    {
      mhm[i] = 0.;
    }

    //acc_assign(idx1sq,mhm,0.);

    for(int i = 0; i < idx1; i++)
    {
      for(int j = 0; j < idx1; j++)
      {
        for(int k = 0; k < Ntim2; k++)
        {
          mhm[idx1*i+j] += mconjh[Ntim2*i+k] * mnxi[idx1*k+j];
        }
      }
    }

    if(prl>0)
    {
      printf("mhm before diag: \n");
      for(int i = 0; i < idx1; i++) 
      {
        for(int j = 0; j < idx1; j++)
        {
          printf("%6.4f ", mhm[i*idx1+j]);
        }
        printf("\n");
      }
      printf("\n");
    }

    FP2* eigvalmhm = new FP2[idx1];
    for(int i = 0; i < idx1; i++)
    {
      eigvalmhm[i] = 0.;
    }

    //acc_assign(idx1,eigvalmhm,0.);

    DiagonalizeP(mhm,eigvalmhm,idx1);

    if(prl>0)
    {
      printf("mhm after diag: \n");
      for(int i = 0; i < idx1; i++) 
      {    
        for(int j = 0; j < idx1; j++) 
        {    
          printf("%6.4f ", mhm[i*idx1+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }

    FP2* mhmconj = new FP2[idx1sq];
    
    for(int i = 0; i < idx1sq; i++)
    {
      mhmconj[i] = 0.;
    }

    for(int i = 0; i < idx1; i++)
    {
      for(int j = 0; j < idx1; j++)
      {
        mhmconj[i*idx1+j] = mhm[j*idx1+i];
      }
    }
   
    FP2* a = new FP2[idx1dim];
    for(int i = 0; i < idx1dim; i++)
    {   
      a[i] = 0.; 
    }   
    
    for(int i = 0; i < Ntim2; i++)
    {
      for(int j = 0; j < idx1; j++)
      {
        for(int k = 0; k < idx1; k++)
        {
          a[idx1*i+j] += mnxi[idx1*i+k] * mhmconj[idx1*k+j];
        }
      }
    }

    if(prl>0)
    {
      printf("mhm after diag and mult: \n");
      for(int i = 0; i < Ntim2; i++) 
      {
        for(int j = 0; j < idx1; j++)
        {
          printf("%6.4f ", mhm[i*idx1+j]);
        }
        printf("\n");
      }
      printf("\n");
    }

    for(int i = 0; i < idx1; i++)
    {
      if(eigvalmhm[i] > -(c*c))
      idx2++;
    }

    if(prl>0)
    {
      printf("eigvalmhm: \n");
      for(int i = 0; i < idx1; i++) 
      {    
        printf("%6.4f ",eigvalmhm[i]); 
      }    
      printf("\n");
    }

    FP2* eigvalmhm2 = new FP2[idx2];
    for(int i = 0; i < idx2; i++)
    {
      eigvalmhm2[i] = 0.;
    }

    for(int i = idx1 - idx2; i < idx1; i++)
    {
      eigvalmhm2[i-idx1+idx2] = eigvalmhm[i];   
    }

    int Nidx2dim = N * idx2;
    FP2* cl = new FP2[Nidx2dim];
    
    for(int i = 0; i < Nidx2dim; i++)
    {
      cl[i] = 0.;
    }
    
    //acc_assign(Nidx2dim,cl,0.);

    for(int i = idx1 - idx2; i < idx1; i++)
    {
      for(int j = 0; j < N; j++)
      {
        cl[idx2*j+i-idx1+idx2] = mhm[idx1*j+i];
      }
    }

    if(prl>0)
    {
      printf("cl: \n");
      for(int i = 0; i < N; i++) 
      {    
        for(int j = 0; j < idx2; j++) 
        {    
          printf("%6.4f ", cl[i*idx2+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }

    FP2* clconj = new FP2[Nidx2dim];

    for(int i = 0; i < Nidx2dim; i++)
    {
      clconj[i] = 0.;
    }

    //acc_assign(Nidx2dim,clconj,0.);

    for(int i = 0; i < idx2; i++)
    {
      for(int j = 0; j < N; j++)
      {
        clconj[N*i+j] = cl[idx2*j+i];
      }
    }

    if(prl>0)
    {
      printf("clconj: \n");
      for(int i = 0; i < idx2; i++) 
      {    
        for(int j = 0; j < N; j++) 
        {    
          printf("%6.4f ", clconj[i*N+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }
    
    FP2* clconjs = new FP2[Nidx2dim];
    
    for(int i = 0; i < Nidx2dim; i++)
    {
      clconjs[i] = 0.;
    }

    //acc_assign(Nidx2dim,clconjs,0.);

    for(int i = 0; i < idx2; i++)
    {
      for(int j = 0; j < N; j++)
      {
        for(int k = 0; k < N; k++)
        {
          clconjs[N*i+j] += clconj[N*i+k] * s[N*k+j];
        }
      }
    }

    if(prl>0)
    {
      printf("clconjs: \n");
      for(int i = 0; i < idx2; i++) 
      {    
        for(int j = 0; j < N; j++) 
        {    
          printf("%6.4f ", clconjs[i*N+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }

    int idx2sq = idx2 * idx2;
    FP2* clscl = new FP2[idx2sq];
    
    for(int i = 0; i < idx2sq; i++)
    {
      clscl[i] = 0.;
    }

    //acc_assign(idx2sq,clscl,0.);

    for(int i = 0; i < idx2; i++)
    {
      for(int j = 0; j < idx2; j++)
      {
        for(int k = 0; k < N; k++)
        {
          clscl[idx2*i+j] += clconjs[N*i+k] * cl[idx2*k+j];
        }
      }
    }

    if(prl>0)
    {
      printf("clscl before diag: \n");
      for(int i = 0; i < idx2; i++) 
      {    
        for(int j = 0; j < idx2; j++) 
        {    
          printf("%6.4f ", clscl[i*idx2+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }

    FP2* eigvalclscl = new FP2[idx2];
    
    for(int i = 0; i < idx2; i++)
    {
      eigvalclscl[i] = 0.;
    }

    //acc_assign(idx2,eigvalclscl,0.);

    DiagonalizeP(clscl,eigvalclscl,idx2);

    if(prl>0)
    {
      printf("clscl after diag: \n");
      for(int i = 0; i < idx2; i++) 
      {    
        for(int j = 0; j < idx2; j++) 
        {    
          printf("%6.4f ", clscl[i*idx2+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }

    FP2* clscltra = new FP2[idx2sq];

    for(int i = 0; i < idx2sq; i++)
    {
      clscltra[i] = 0.;
    }

    for(int i = 0; i < idx2; i++)
    {
      for(int j = 0; j < idx2; j++)
      {
        clscltra[i*idx2+j] = clscl[j*idx2+i];
      }
    }

    for(int i = 0; i < idx2; i++)
    {
      if(eigvalclscl[i] > 1e-14)
      idx3++;
    }

    if(prl>0)
    {
      printf("eigvalclscl: \n");
      for(int i = 0; i < idx2; i++) 
      {    
        printf("%20.20f ",eigvalclscl[i]);
      }    
      printf("\n");
    }

    int idx2idx3dim = idx2 * idx3;
    FP2* clsclnorm = new FP2[idx2idx3dim];
    
    for(int i = 0; i < idx2idx3dim; i++)
    {
      clsclnorm[i] = 0.;
    }

    //acc_assign(idx2idx3dim,clsclnorm,0.);

    for(int i = idx2 - idx3; i < idx3; i++)
    {
      for(int j = 0; j < idx2; j++)
      {
        clsclnorm[idx3*j+i-idx2+idx3] = clscltra[idx2*j+i] / sqrt(eigvalclscl[i]);
      }
    }

    if(prl>0)
    {
      printf("clsclnorm: \n");
      for(int i = 0; i < idx2; i++) 
      {    
        for(int j = 0; j < idx3; j++) 
        {    
          printf("%6.4f ", clsclnorm[i*idx3+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }

    FP2* clsclre = new FP2[idx2idx3dim];
    
    for(int i = 0; i < idx2idx3dim; i++)
    {
      clsclre[i] = 0.;
    }

    //acc_assign(idx2idx3dim,clsclre,0.);

    for(int i = 0; i < idx3; i++)
    {
      for(int j = 0; j < idx2; j++)
      {
        clsclre[idx3*j+i] = clscltra[idx2*j+i];
      }
    }

    if(prl>0)
    {
      printf("clsclre: \n");
      for(int i = 0; i < idx2; i++) 
      {    
        for(int j = 0; j < idx3; j++) 
        {    
          printf("%6.4f ", clsclre[i*idx3+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }

    FP2* clsclconj = new FP2[idx2idx3dim];
    
    for(int i = 0; i < idx2idx3dim; i++)
    {
      clsclconj[i] = 0.;
    }

    //acc_assign(idx2idx3dim,clsclconj,0.);

    for(int i = 0; i < idx3; i++)
    {
      for(int j = 0; j < idx2; j++)
      {
        clsclconj[idx2*i+j] = clsclre[idx3*j+i];
      }
    }

    if(prl>0)
    {
      printf("clsclconj: \n");
      for(int i = 0; i < idx3; i++) 
      {    
        for(int j = 0; j < idx2; j++) 
        {    
          printf("%6.4f ", clsclconj[i*idx2+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }

    FP2* n1 = new FP2[idx2sq]; 
    for(int i = 0; i < idx2sq; i++)
    {
      n1[i] = 0.;
    }

    //acc_assign(idx2sq,n1,0.);

    for(int i = 0; i < idx2; i++)
    {
      for(int j = 0; j < idx2; j++)
      {
        for(int k = 0; k < idx3; k++)
        {
          n1[idx2*i+j] += clsclnorm[idx3*i+k] * clsclconj[idx2*k+j];
        }
      }
    }

    if(prl>0)
    {
      printf("n1: \n");
      for(int i = 0; i < idx2; i++) 
      {    
        for(int j = 0; j < idx2; j++) 
        {    
          printf("%6.4f ", n1[i*idx2+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }

    FP2* n2 = new FP2[Nidx2dim];    
    for(int i = 0; i < Nidx2dim; i++)
    {
      n2[i] = 0.;
    }

    //acc_assign(Nidx2dim,n2,0.);

    for(int i = 0; i < idx2; i++)
    {
      for(int j = 0; j < N; j++)
      {
        for(int k = 0; k < idx2; k++)
        {
          n2[N*i+j] += n1[idx2*i+k] * clconj[N*k+j];
        }
      }
    }

    if(prl>0)
    {
      printf("n2: \n");
      for(int i = 0; i < idx2; i++) 
      {    
        for(int j = 0; j < N; j++) 
        {    
          printf("%6.4f ", n2[i*N+j]);
        }    
        printf("\n");
      }    
      printf("\n");
    }

    FP2* n3 = new FP2[Nidx2dim];    
    for(int i = 0; i < Nidx2dim; i++)
    {
      n3[i] = 0.;
    }

    //acc_assign(Nidx2dim,n3,0.);

    for(int i = 0; i < idx2; i++)
    {
      for(int j = 0; j < N; j++)
      {
        for(int k = 0; k < N; k++)
        {
          n3[N*i+j] += n2[N*i+k] * s[N*k+j];
        }
      }
    }

    FP2* n3conj = new FP2[Nidx2dim];
    
    for(int i = 0; i < Nidx2dim; i++)
    {
      n3conj[i] = 0.;
    }

    //acc_assign(Nidx2dim,n3conj,0.);

    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < idx2; j++)
      {
        n3conj[idx2*i+j] = n3[N*j+i];
      }
    }

    FP2* n3conje = new FP2[Nidx2dim];
    
    for(int i = 0; i < Nidx2dim; i++)
    {
      n3conje[i] = 0.;
    }

    //acc_assign(Nidx2dim,n3conje,0.);

    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < idx2; j++)
      {
        n3conje[idx2*i+j] = n3conj[idx2*i+j] * eigvalmhm2[j];
      }
    }
    
    if(prl>0)
    {
      printf("n3conje: \n");
      for(int i = 0; i < N; i++)
      {
        for(int j = 0; j < idx2; j++)
        {
          printf("%6.4f ", n3conje[i*N+j]);
        }
        printf("\n");
      }

      printf("n3: \n");
      for(int i = 0; i < idx2; i++)
      {
        for(int j = 0; j < N; j++)
        {
          printf("%6.4f ", n3[i*N+j]);
        }
        printf("\n");
      }
    }

    for(int i = 0; i < Nsq; i++)
    {
      hcore[i] = 0.;
    }

    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < N; j++)
      {
        for(int k = 0; k < idx2; k++)
        {
          hcore[N*i+j] += n3conje[idx2*i+k] * n3[N*k+j];
        }
      }
    }

    ofstream filehcore;
    filehcore.open("Hcore");
    for(int i = 0; i < N; i++) 
    {    
      for(int j = 0; j < N; j++) 
      {    
        filehcore << hcore[i*N+j] << " "; 
      }    
      filehcore << "\n";
    }    
    filehcore.close();

    printf("hcore: \n");
    for(int i = 0; i < N; i++) 
    {    
      for(int j = 0; j < N; j++) 
      {    
        printf("%10.8f ", hcore[i*N+j]);
      }    
      printf("\n");
    }

    if(prl>0)
    {
      printf("idx1: %i idx2: %i idx3: %i \n",idx1,idx2,idx3);
    }

    delete [] h;
    delete [] m;
    delete [] eigvalm;
    delete [] mconj;
    delete [] mnxi;
    delete [] eigvecconjm;
    delete [] mconjh;
    delete [] mhm;
    delete [] eigvalmhm;
    delete [] eigvalmhm2;
    delete [] mhmconj;
    delete [] a;
    delete [] cl;
    delete [] clconj;
    delete [] clconjs;
    delete [] clscl;
    delete [] clscltra;
    delete [] eigvalclscl;
    delete [] clsclnorm;
    delete [] clsclconj;
    delete [] clsclre;
    delete [] n1;
    delete [] n2;
    delete [] n3;
    delete [] n3conj;
    delete [] n3conje;

    return;
}

