// Included by becke.cpp: bounded observation tiles, reused across density frames.
#include <stdexcept>
#include <climits>
#include <memory>
#include <exception>
#include <thread>
#include <cerrno>
#if USE_ACC
#include <cublas_v2.h>
#include <cuda_runtime.h>

namespace {
// Device selection is local to a host thread. Keep it selected until all of
// that worker's OpenACC buffers and cuBLAS handles have been destroyed.
class BFieldDevice {
  int previous_;
 public:
  explicit BFieldDevice(int device) : previous_(acc_get_device_num(acc_device_nvidia)) {
    acc_set_device_num(device,acc_device_nvidia);
  }
  ~BFieldDevice() { acc_set_device_num(previous_,acc_device_nvidia); }
};

// Join even if allocating/starting a later worker throws.
class BFieldThreads {
 public:
  vector<thread> workers;
  ~BFieldThreads() { join(); }
  void join() {
    for (auto& worker:workers) if (worker.joinable()) worker.join();
  }
};

int bfield_device_count(int points)
{
  int count=acc_get_num_devices(acc_device_nvidia);
  if (const char* value=getenv("SGPU_B_NGPU"))
  {
    char* end=NULL;
    errno=0;
    long requested=strtol(value,&end,10);
    if (errno || end==value || *end || requested<0 || requested>INT_MAX)
      throw invalid_argument("SGPU_B_NGPU must be a nonnegative integer (0 means all visible GPUs)");
    if (requested) count=min(count,(int)requested);
  }
  return min(count,points);
}

// Angular tables may have a newer device value than their host mirror.
// Snapshot them on the caller's device before starting independent workers.
void bfield_host_snapshot(double* source,size_t count,double* destination)
{
  if (acc_is_present(source,count*sizeof(double)))
    acc_memcpy_from_device(destination,acc_deviceptr(source),count*sizeof(double));
  else copy(source,source+count,destination);
}

void bfield_cuda_check(cudaError_t status)
{
  if (status != cudaSuccess) throw runtime_error(cudaGetErrorString(status));
}
void bfield_blas_check(cublasStatus_t status)
{
  if (status != CUBLAS_STATUS_SUCCESS) throw runtime_error("B-field cuBLAS failure");
}

// Each allocation owns its host and device storage, including exception paths.
class BFieldBuffer {
  unique_ptr<double[]> host_;
  size_t count_;
 public:
  explicit BFieldBuffer(size_t count) : host_(new double[count]), count_(count) {
    if (!acc_create(host_.get(),count_*sizeof(double))) throw bad_alloc();
  }
  ~BFieldBuffer() { acc_delete(host_.get(),count_*sizeof(double)); }
  double* data() { return host_.get(); }
  BFieldBuffer(const BFieldBuffer&)=delete;
  BFieldBuffer& operator=(const BFieldBuffer&)=delete;
};
class BFieldBlas {
 public:
  cublasHandle_t value;
  BFieldBlas() { bfield_blas_check(cublasCreate(&value)); }
  ~BFieldBlas() { cublasDestroy(value); }
  BFieldBlas(const BFieldBlas&)=delete;
  BFieldBlas& operator=(const BFieldBlas&)=delete;
};

// Phi and Q are column-major M x N; the host density is row-major N x N.
// The transpose below converts its memory interpretation, not its physical indices.
void bfield_contract(cublasHandle_t handle, int M, int N, double* phi, double* p, double* q)
{
  acc_wait_all();
  const double one = 1., zero = 0.;
  #pragma acc host_data use_device(phi,p,q)
  {
    bfield_blas_check(cublasDgemm(handle,CUBLAS_OP_N,CUBLAS_OP_T,
        M,N,N,&one,phi,M,p,N,&zero,q,M));
  }
  // The OpenACC evaluators currently use synchronous launches. Explicitly join
  // the cuBLAS stream before an OpenACC consumer accesses Q.
  bfield_cuda_check(cudaDeviceSynchronize());
}
}
#endif

static void compute_B_field_frames_single(int natoms, int* atno, double* coords, bool gbasis,
    vector<vector<double> >& basis, int nframes, const double* const* pao,
    const double* const* paodt, int nrad, int nang, double* ang_g, double* ang_w,
    int gsa, double* grid, double* B, bool device_output, int prl)
{
  if (nframes<=0 || gsa<=0) return;
  const size_t nb = (size_t)nframes*gsa*3;
#if USE_ACC
  if (gbasis)
  {
    if (device_output && !acc_is_present(B,nb*sizeof(double)))
      throw invalid_argument("Device B output must already be present");
    const bool skip_dt = read_int("SKIP_DJ")!=0;
    if (skip_dt) printf(" -- WARNING: skipping dJ/dt term -- \n");
    bool have_dt = false;
    if (!skip_dt && paodt)
      for (int f=0;f<nframes;f++) have_dt = have_dt || paodt[f]!=NULL;
    const int N = basis.size();
    const size_t N2 = (size_t)N*N;
    const int gs = nrad*nang;
    int tile = 16;
    if (const char* value = getenv("SGPU_B_TILE")) tile = max(1,min(64,atoi(value)));
    // Bound AO/Q tables; a single observation is the irreducible minimum.
    const size_t per_point = sizeof(double)*(size_t)gs*N*(have_dt ? 6 : 5);
    if (per_point>0) tile=min(tile,(int)max((size_t)1,(size_t)(512*1024*1024)/per_point));
    tile = min(tile,gsa);
    if (gs<=0 || (size_t)tile*gs>INT_MAX)
      throw invalid_argument("Invalid B-field quadrature/tile size");
    const int M = tile*gs;
    const size_t mn = (size_t)M*N, grad_size = 3*mn;
    vector<int> ends(natoms);
    get_imax_n2i(natoms,N,basis,ends.data());
    const size_t q2_size=have_dt ? mn : 1, p2_size=have_dt ? N2 : 1;
    // Apply precisely the original element screening once per input frame.
    vector<double> screened((size_t)nframes*N2), screened_dt(have_dt ? (size_t)nframes*N2 : 0);
    vector<char> live(nframes,0);
    for (int frame=0;frame<nframes;frame++)
    {
      if (!pao || !pao[frame]) throw invalid_argument("Missing B-field density frame");
      const bool dt=have_dt && paodt[frame]!=NULL;
      for (size_t k=0;k<N2;k++)
      {
        double value=pao[frame][k];
        screened[(size_t)frame*N2+k]=fabs(value)>1.e-10 ? value : 0.;
        live[frame]=live[frame] || fabs(value)>1.e-10;
        if (dt)
        {
          value=paodt[frame][k];
          screened_dt[(size_t)frame*N2+k]=fabs(value)>1.e-10 ? value : 0.;
          live[frame]=live[frame] || fabs(value)>1.e-10;
        }
      }
    }
    if (find(live.begin(),live.end(),1)==live.end())
    {
      if (device_output)
      {
        #pragma acc parallel loop present(B[0:nb])
        for (size_t k=0;k<nb;k++) B[k]=0.;
      }
      else fill(B,B+nb,0.);
      return;
    }
    BFieldBuffer obs_store(6*(size_t)gsa),base_store(6*(size_t)gs),weight_store(gs);
    BFieldBuffer shift_store(6*(size_t)M),tmp_store(3*(size_t)M);
    BFieldBuffer phi_store(mn),grad_store(grad_size),q1_store(mn),q2_store(q2_size);
    BFieldBuffer p1_store(N2),p2_store(p2_size),contribution_store(3*(size_t)M),out_store(nb);
    double *obs=obs_store.data(),*base=base_store.data(),*weights=weight_store.data();
    double *shifted=shift_store.data(),*tmp=tmp_store.data();
    double *phi=phi_store.data(),*grad=grad_store.data(),*q1=q1_store.data(),*q2=q2_store.data();
    double *p1=p1_store.data(),*p2=p2_store.data(),*contribution=contribution_store.data(),*bout=out_store.data();
    copy(grid,grid+6*(size_t)gsa,obs);
    #pragma acc update device(obs[0:6*(size_t)gsa])
    BFieldBlas blas;
    generate_central_grid_2d(-1,1,base,weights,1.,nrad,nang,ang_g,ang_w);
    // Preserve the original observation-centered quadrature, regularization,
    // weights, and 1/c contribution. AO radii are independently recentered.
    #pragma acc parallel loop present(base[0:6*gs])
    for (int q=0;q<gs;q++)
    {
      double r = base[6*q+3];
      base[6*q+3] = r+1.e-8;
      base[6*q+4] = r*r+1.e-8;
      base[6*q+5] = r*r*r+1.e-8;
    }
    for (int first=0;first<gsa;first+=tile)
    {
      const int nt = min(tile,gsa-first), used = nt*gs;
      #pragma acc parallel loop present(phi[0:mn],grad[0:grad_size])
      for (size_t k=0;k<grad_size;k++)
      {
        grad[k] = 0.;
        if (k<mn) phi[k] = 0.;
      }
      // Every Gaussian AO value and gradient is evaluated once for this tile,
      // independent of which atom supplied the other member of an AO pair.
      for (int atom=0;atom<natoms;atom++)
      {
        const double ax=coords[3*atom], ay=coords[3*atom+1], az=coords[3*atom+2];
        #pragma acc parallel loop present(shifted[0:6*M],base[0:6*gs],obs[0:6*(size_t)gsa])
        for (int q=0;q<used;q++)
        {
          int source=q%gs, observation=first+q/gs;
          double x=base[6*source]+obs[6*observation]-ax;
          double y=base[6*source+1]+obs[6*observation+1]-ay;
          double z=base[6*source+2]+obs[6*observation+2]-az;
          shifted[6*q]=x; shifted[6*q+1]=y; shifted[6*q+2]=z;
          shifted[6*q+3]=sqrt(x*x+y*y+z*z);
          shifted[6*q+4]=base[6*source+4];
          shifted[6*q+5]=base[6*source+5];
        }
        int start=atom ? ends[atom-1] : 0;
        for (int ao=start;ao<ends[atom];ao++)
        {
          const vector<double>& b=basis[ao];
          int l=b[1], m=b[2], ng=b[3], norm_offset=10+ng;
          for (int primitive=0;primitive<ng;primitive++)
          {
            const double exponent=b[10+primitive], coefficient=b[norm_offset+primitive];
            eval_ghd(used,shifted,tmp,l,m,coefficient,exponent);
            #pragma acc parallel loop present(phi[0:mn],tmp[0:3*M])
            for (int q=0;q<used;q++) phi[(size_t)ao*M+q] += tmp[q];
            eval_pd_gh(used,shifted,tmp,l+1,l,m,coefficient,exponent);
            #pragma acc parallel loop present(grad[0:grad_size],tmp[0:3*M])
            for (int q=0;q<3*used;q++) grad[(size_t)ao*3*M+q] += tmp[q];
          }
        }
      }
      // The same bounded AO tile is reused across all frames. Workspaces have
      // call-local ownership, so changing geometry/basis/grid cannot reuse stale data.
      for (int frame=0;frame<nframes;frame++)
      {
        const bool dt = have_dt && paodt[frame]!=NULL;
        if (!live[frame])
        {
          #pragma acc parallel loop present(bout[0:nb])
          for (int point=0;point<3*nt;point++)
            bout[((size_t)frame*gsa+first)*3+point]=0.;
          continue;
        }
        copy(screened.data()+(size_t)frame*N2,screened.data()+(size_t)(frame+1)*N2,p1);
        if (dt) copy(screened_dt.data()+(size_t)frame*N2,screened_dt.data()+(size_t)(frame+1)*N2,p2);
        #pragma acc update device(p1[0:N2])
        bfield_contract(blas.value,M,N,phi,p1,q1);
        if (dt)
        {
          #pragma acc update device(p2[0:N2])
          bfield_contract(blas.value,M,N,phi,p2,q2);
        }
        #pragma acc parallel loop present(q1[0:mn],q2[0:q2_size],grad[0:grad_size],base[0:6*gs],weights[0:gs],contribution[0:3*M])
        for (int q=0;q<used;q++)
        {
          double jx=0.,jy=0.,jz=0.,dx=0.,dy=0.,dz=0.;
          #pragma acc loop seq
          for (int ao=0;ao<N;ao++)
          {
            const size_t index=(size_t)ao*M+q, gi=(size_t)ao*3*M+3*q;
            double v=q1[index];
            jx+=v*grad[gi]; jy+=v*grad[gi+1]; jz+=v*grad[gi+2];
            if (dt)
            {
              v=q2[index];
              dx+=v*grad[gi]; dy+=v*grad[gi+1]; dz+=v*grad[gi+2];
            }
          }
          int s=q%gs;
          double x=base[6*s],y=base[6*s+1],z=base[6*s+2],w=weights[s];
          double inv3=w/base[6*s+5], inv2=dt ? w/(ALPHAFSC*base[6*s+4]) : 0.;
          contribution[3*q]=(y*jz-z*jy)*inv3+(y*dz-z*dy)*inv2;
          contribution[3*q+1]=-(x*jz-z*jx)*inv3-(x*dz-z*dx)*inv2;
          contribution[3*q+2]=(x*jy-y*jx)*inv3+(x*dy-y*dx)*inv2;
        }
        #pragma acc parallel loop gang present(contribution[0:3*M],bout[0:nb])
        for (int point=0;point<nt;point++)
        {
          double bx=0.,by=0.,bz=0.;
          #pragma acc loop vector reduction(+:bx,by,bz)
          for (int source=0;source<gs;source++)
          {
            int q=point*gs+source;
            bx+=contribution[3*q]; by+=contribution[3*q+1]; bz+=contribution[3*q+2];
          }
          size_t out=((size_t)frame*gsa+first+point)*3;
          bout[out]=-bx; bout[out+1]=-by; bout[out+2]=-bz;
        }
      }
    }
    if (device_output)
    {
      #pragma acc parallel loop present(B[0:nb],bout[0:nb])
      for (size_t k=0;k<nb;k++) B[k]=bout[k];
    }
    else
    {
      #pragma acc update self(bout[0:nb])
      copy(bout,bout+nb,B);
    }
    if (prl>1) printf("  Gaussian B: frames=%i tile=%i workspace=%zu bytes\n",nframes,tile,
        sizeof(double)*(6*(size_t)gsa+7*(size_t)gs+12*(size_t)M+2*mn+grad_size+q2_size+N2+p2_size+nb));
    return;
  }
#endif
  // Retain the existing STO and CPU implementation without changing its physics.
  for (int f=0;f<nframes;f++)
    compute_B_field(natoms,atno,coords,gbasis,basis,const_cast<double*>(pao[f]),
        paodt ? const_cast<double*>(paodt[f]) : NULL,
        nrad,nang,ang_g,ang_w,gsa,grid,B+(size_t)f*gsa*3,prl);
#if USE_ACC
  if (device_output)
  {
    #pragma acc update device(B[0:nb])
  }
#endif
}

void compute_B_field_frames(int natoms, int* atno, double* coords, bool gbasis,
    vector<vector<double> >& basis, int nframes, const double* const* pao,
    const double* const* paodt, int nrad, int nang, double* ang_g, double* ang_w,
    int gsa, double* grid, double* B, bool device_output, int prl)
{
  if (nframes<=0 || gsa<=0) return;
#if USE_ACC
  if (gbasis)
  {
    const int ngpu=bfield_device_count(gsa);
    if (ngpu>1)
    {
      if (nrad<=0 || nang<=0 || nrad>INT_MAX/nang)
        throw invalid_argument("Invalid B-field quadrature size");
      const size_t nb=(size_t)nframes*gsa*3;
      if (device_output && !acc_is_present(B,nb*sizeof(double)))
        throw invalid_argument("Device B output must already be present");
      const int caller=acc_get_device_num(acc_device_nvidia);
      const int visible=acc_get_num_devices(acc_device_nvidia);
      acc_wait_all();
      vector<double> angular(3*(size_t)nang),weights(nang),result(nb);
      bfield_host_snapshot(ang_g,angular.size(),angular.data());
      bfield_host_snapshot(ang_w,weights.size(),weights.data());
      vector<exception_ptr> errors(ngpu);
      if (prl>0) printf("  Gaussian B: GPUs=%i frames=%i observations=%i\n",ngpu,nframes,gsa);
      BFieldThreads threads;
      threads.workers.reserve(ngpu);
      for (int worker=0;worker<ngpu;worker++)
      {
        const int device=(caller+worker)%visible;
        const int first=(size_t)gsa*worker/ngpu;
        const int count=(size_t)gsa*(worker+1)/ngpu-first;
        if (prl>1) printf("  Gaussian B: device=%i observations=[%i,%i)\n",device,first,first+count);
        threads.workers.emplace_back([&,worker,device,first,count]() {
          try
          {
            BFieldDevice selected(device);
            // Each GPU owns its angular tables, cuBLAS handle, bounded AO
            // tiles, and output. No OpenACC host mapping is shared by workers.
            BFieldBuffer angular_store(angular.size()),weight_store(weights.size());
            double* ag=angular_store.data();
            double* aw=weight_store.data();
            copy(angular.begin(),angular.end(),ag);
            copy(weights.begin(),weights.end(),aw);
            acc_update_device(ag,angular.size()*sizeof(double));
            acc_update_device(aw,weights.size()*sizeof(double));
            vector<double> local((size_t)nframes*count*3);
            compute_B_field_frames_single(natoms,atno,coords,gbasis,basis,
                nframes,pao,paodt,nrad,nang,ag,aw,count,grid+6*(size_t)first,
                local.data(),false,prl);
            for (int frame=0;frame<nframes;frame++)
              copy(local.data()+(size_t)frame*count*3,
                  local.data()+(size_t)(frame+1)*count*3,
                  result.data()+((size_t)frame*gsa+first)*3);
          }
          catch (...) { errors[worker]=current_exception(); }
        });
      }
      threads.join();
      for (const auto& error:errors) if (error) rethrow_exception(error);
      // The caller thread never switches device. Device-output mode preserves
      // the host mirror, including when the caller started on a nonzero GPU.
      if (device_output)
        acc_memcpy_to_device(acc_deviceptr(B),result.data(),nb*sizeof(double));
      else copy(result.begin(),result.end(),B);
      return;
    }
  }
#endif
  compute_B_field_frames_single(natoms,atno,coords,gbasis,basis,nframes,pao,paodt,
      nrad,nang,ang_g,ang_w,gsa,grid,B,device_output,prl);
}
