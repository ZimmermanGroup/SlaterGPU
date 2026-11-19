#include "cint_soc.h"
extern "C" {
#include "cint_funcs.h"
}
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using namespace std;

namespace cintw {
namespace {

struct Optimizer
{
  CINTOpt* p = nullptr;
  ~Optimizer() { if (p) CINTdel_optimizer(&p); }
};
struct Shell { int n, offset; };
struct Pair { int i, j; };

// Only a bounded collection of output shell tiles is private to each thread.
// A flush holds one row-shell lock; there are no nested locks or AO atomics.
class Tiles {
  struct Tile { int row, col; list<size_t>::iterator position; vector<double> v; };
  const vector<Shell>& shells;
  int n;
  size_t limit, used = 0;
  list<size_t> recent;
  vector<double>& output;
  vector<omp_lock_t>& locks;
  unordered_map<size_t, Tile> tiles;
  void flush(Tile& t)
  {
    const Shell a = shells[t.row], b = shells[t.col];
    omp_set_lock(&locks[t.row]);
    for (int k=0;k<3;k++)
    for (int i=0;i<a.n;i++)
    for (int j=0;j<b.n;j++)
      output[(size_t(k)*n+a.offset+i)*n+b.offset+j] += t.v[(size_t(k)*a.n+i)*b.n+j];
    omp_unset_lock(&locks[t.row]);
  }

public:
  Tiles(const vector<Shell>& s, int nao, size_t bytes, vector<double>& out, vector<omp_lock_t>& l)
    : shells(s), n(nao), limit(bytes), output(out), locks(l) {}

  void add(int row, int col, const vector<double>& block, double factor)
  {
    const size_t key=size_t(row)*shells.size()+col;
    auto it=tiles.find(key);
    if (it==tiles.end())
    {
      const size_t bytes=block.size()*sizeof(double);
      // One tile is allowed even if it exceeds the requested cache capacity.
      while (!tiles.empty() && used+bytes>limit)
      {
        auto victim=tiles.find(recent.front());
        flush(victim->second);
        used-=victim->second.v.size()*sizeof(double);
        recent.pop_front(); tiles.erase(victim);
      }
      recent.push_back(key);
      Tile t{row,col,prev(recent.end()),vector<double>(block.size(),0.)};
      it=tiles.emplace(key,move(t)).first;
      used+=bytes;
    }
    recent.splice(recent.end(),recent,it->second.position);
    for (size_t i=0;i<block.size();++i) it->second.v[i]+=factor*block[i];
  }
  void finish() { for (auto& p:tiles) flush(p.second); tiles.clear(); }
};

double antisymmetry_error(const vector<double>& a, int n)
{
  double error=0.;
  for (int k=0;k<3;k++)
  for (int p=0;p<n;p++)
  for (int q=0;q<n;q++)
    error=max(error,abs(a[(size_t(k)*n+p)*n+q]+a[(size_t(k)*n+q)*n+p]));
  return error;
}
}

SOCIntegrals breit_pauli_soc(int natm, int nbas, int* atm, int* bas, double* env,
                            bool cartesian, const double* density, const SOCOptions& options)
{
  if (natm<=0 || nbas<=0 || !atm || !bas || !env || options.threads<0)
    throw invalid_argument("breit_pauli_soc: invalid basis or thread count");
  if (env[PTR_RANGE_OMEGA]!=0.)
    throw invalid_argument("breit_pauli_soc: range-separated Coulomb is unsupported");
  for (int a=0;a<natm;a++)
  if (atm[a*ATM_SLOTS+CHARGE_OF]<0 || atm[a*ATM_SLOTS+NUC_MOD_OF]==3)
    throw invalid_argument("breit_pauli_soc: ECP/fractional nuclear charges are unsupported");

  vector<Shell> shells;
  vector<Pair> pairs;
  int n=0;

  for (int i=0;i<nbas;i++)
  {
    const int d=cartesian ? CINTcgto_cart(i,bas) : CINTcgto_spheric(i,bas);
    if (d<=0 || d>numeric_limits<int>::max()-n)
      throw invalid_argument("breit_pauli_soc: invalid shell size");

    shells.push_back(Shell{d,n}); n+=d;
    for (int j=0;j<=i;j++)
      pairs.push_back(Pair{i,j});
  }
  const size_t n2=size_t(n)*n;
  if (n2>numeric_limits<size_t>::max()/3)
    throw overflow_error("breit_pauli_soc: AO matrix size overflow");
  if (density)
  {
    double scale  = 1.;
    double defect = 0.;
    for (int p=0;p<n;p++)
    for (int q=0;q<n;q++)
    {
      const double v=density[size_t(p)*n+q];
      if (!isfinite(v)) throw invalid_argument("breit_pauli_soc: nonfinite reference density");

      scale = max(scale,abs(v));
      defect = max(defect,abs(v-density[size_t(q)*n+p]));
    }
    if (defect>1.e-10*scale)
      throw invalid_argument("breit_pauli_soc: reference density must be symmetric and spin-summed");
  }
  SOCIntegrals result;
  result.nao=n;
  result.nuclear.assign(3*n2,0.);

  const int threads=options.threads ? options.threads : omp_get_max_threads();
  auto* nuclear=cartesian ? int1e_pnucxp_cart : int1e_pnucxp_sph;
  Optimizer opt1;
  int1e_pnucxp_optimizer(&opt1.p,atm,natm,bas,nbas,env);

  #pragma omp parallel num_threads(threads)
  {
    vector<double> buf,cache;
    #pragma omp for schedule(dynamic,8)
    for (size_t t=0;t<pairs.size();t++)
    {
      const Pair pair = pairs[t];
      const Shell a = shells[pair.i];
      const Shell b = shells[pair.j];
      int shls[2]={pair.i,pair.j};
      buf.resize(3*size_t(a.n)*b.n);
      const auto need=nuclear(nullptr,nullptr,shls,atm,natm,bas,nbas,env,opt1.p,nullptr);
      if (cache.size()<size_t(need)) cache.resize(need);
      if (!nuclear(buf.data(),nullptr,shls,atm,natm,bas,nbas,env,opt1.p,cache.data())) continue;
      for (int k=0;k<3;k++)
      for (int i=0;i<a.n;i++)
      for (int j=0;j<b.n;j++)
      {
        const double v=buf[(size_t(k)*b.n+j)*a.n+i];
        result.nuclear[(size_t(k)*n+a.offset+i)*n+b.offset+j]=v;
        if (pair.i!=pair.j)
          result.nuclear[(size_t(k)*n+b.offset+j)*n+a.offset+i]=-v;
      }
    }
  }
  result.nuclear_antisymmetry_error=antisymmetry_error(result.nuclear,n);
  if (!density) return result;

  result.somf.assign(3*n2,0.);
  if (pairs.size()>size_t(numeric_limits<long long>::max())/pairs.size())
    throw overflow_error("breit_pauli_soc: shell quartet count overflow");

  result.shell_quartets=pairs.size()*pairs.size();
  auto* two=cartesian ? int2e_p1vxp1_cart : int2e_p1vxp1_sph;
  Optimizer opt2;

  int2e_p1vxp1_optimizer(&opt2.p,atm,natm,bas,nbas,env);
  vector<omp_lock_t> locks(nbas);
  for (auto& l:locks) omp_init_lock(&l);
  #pragma omp parallel num_threads(threads)
  {
    Tiles tiles(shells,n,options.tile_cache_bytes,result.somf,locks);
    vector<double> buf,cache,jblock,k1block,k2block;
    #pragma omp for schedule(dynamic,8)
    for (long long t=0;t<static_cast<long long>(result.shell_quartets);t++)
    {
      const Pair ab=pairs[t/pairs.size()],cd=pairs[t%pairs.size()];
      int shls[4]={ab.i,ab.j,cd.i,cd.j};
      const Shell oi=shells[ab.i],oj=shells[ab.j],ok=shells[cd.i],ol=shells[cd.j];
      const size_t block=size_t(oi.n)*oj.n*ok.n*ol.n;
      buf.resize(3*block);
      const auto need=two(nullptr,nullptr,shls,atm,natm,bas,nbas,env,opt2.p,nullptr);
      if (cache.size()<size_t(need)) cache.resize(need);
      if (!two(buf.data(),nullptr,shls,atm,natm,bas,nbas,env,opt2.p,cache.data())) continue;
      // p1vxp1 is antisymmetric in the first pair and symmetric in the second.
      // It has no pair-exchange symmetry. Same-shell blocks already contain
      // both AO orientations, so only distinct SHELL pairs are expanded here.
      for (int swap1=0;swap1<(ab.i==ab.j ? 1 : 2);swap1++)
      for (int swap2=0;swap2<(cd.i==cd.j ? 1 : 2);swap2++)
      {
        const int sa=swap1 ? ab.j:ab.i, sb=swap1 ? ab.i:ab.j;
        const int sc=swap2 ? cd.j:cd.i, sd=swap2 ? cd.i:cd.j;
        const Shell a=shells[sa],b=shells[sb],c=shells[sc],d=shells[sd];
        const double sign=swap1 ? -1.:1.;

        jblock.assign(3*size_t(a.n)*b.n,0.);
        k1block.assign(3*size_t(a.n)*d.n,0.);
        k2block.assign(3*size_t(c.n)*b.n,0.);

        for (int xyz=0;xyz<3;xyz++)
        for (int l=0;l<ol.n;l++)
        for (int k=0;k<ok.n;k++)
        for (int j=0;j<oj.n;j++)
        for (int i=0;i<oi.n;i++)
        {
          const int ia=swap1 ? j:i, ib=swap1 ? i:j;
          const int ic=swap2 ? l:k, id=swap2 ? k:l;
          const double v=sign*buf[xyz*block+((size_t(l)*ok.n+k)*oj.n+j)*oi.n+i];
          jblock[(size_t(xyz)*a.n+ia)*b.n+ib]+=v*density[size_t(d.offset+id)*n+c.offset+ic];
          k1block[(size_t(xyz)*a.n+ia)*d.n+id]+=v*density[size_t(b.offset+ib)*n+c.offset+ic];
          k2block[(size_t(xyz)*c.n+ic)*b.n+ib]+=v*density[size_t(d.offset+id)*n+a.offset+ia];
        }
        tiles.add(sa,sb,jblock,1.);
        tiles.add(sa,sd,k1block,-1.5);
        tiles.add(sc,sb,k2block,-1.5);
      }
    }
    tiles.finish();
  }
  for (auto& l:locks) omp_destroy_lock(&l);
  result.somf_antisymmetry_error=antisymmetry_error(result.somf,n);
  return result;
}
} // namespace cintw
