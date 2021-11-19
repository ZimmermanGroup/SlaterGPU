#ifndef _CINT_WRAPPER_H_
#define _CINT_WRAPPER_H_

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include "cintwrapper.h"

extern "C" {
  #include "cint_funcs.h"
}

using namespace std;

class CINTPrep {
  public: //should probably make this private...
    string xyzfile;
    string basfile;
    string basfile_ri;
    vector< int > atoms;
    vector< double > coord;
    map< int, basis_t > basmap;
    int var_dim;
    int nbas;
    int nenv;
    int *atm;
    int *bas;
    double *env;

    //RI variables
    bool do_ri;
    map< int, basis_t > basmap_ri;
    int var_dim_ri;
    int nbas_ri;
    int nenv_ri;
  
  public:
    CINTPrep(bool doing_ri = false);
    ~CINTPrep();
    int get_var_dim();
    int get_nbas();
    int get_nenv();
    int *get_atm();
    int *get_bas();
    double *get_env();
    void set_atm(int *&atm_in);
    void set_bas(int *&bas_in);
    void set_env(double *&env_in);
    void assign_coords(int natoms, int *atomlist, double *coords, bool in_bohr = true); 

    //RI funcs
    int get_var_dim_ri();
    int get_nbas_ri();
    int get_nenv_ri();

    void read_xyz(string inpxyz);
    void read_bas(string inpbas);
    void read_bas_ri(string inpbas);
    void prep_env();

    void copy_atoms(vector< int > &atoms_copy);
    void copy_coord(vector< double > &coord_copy);

    unordered_map<short, short> anum_to_N;
};

// TODO 
// 1. do ri stuff
// 2. probably should have read_xyz/read_bas take in filenames
//    and manage memory associated with atm/env
// 3. should allow setting of atm, bas, env via int*& using concept
//    of ownership from rust

#endif
