#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "cintprep.h"
#include "elements.h"
#include <algorithm>
#include <cctype>

CINTPrep::CINTPrep(bool doing_ri) {
  atm = NULL;
  bas = NULL;
  env = NULL;

  var_dim = 0;
  nbas = 0;
  nenv = 0;

  do_ri = doing_ri;
  //printf("doing ri: %d\n", do_ri);
  //RI prep
  var_dim_ri = 0;
  nbas_ri = 0;
  nenv_ri = 0;
}

CINTPrep::~CINTPrep() {
  if (atm != NULL) {
    delete [] atm;
  }
  if (bas != NULL) {
    delete [] bas;
  }
  if (env != NULL) {
    delete [] env;
  }
}

int CINTPrep::get_var_dim() {
  return var_dim;
}

int CINTPrep::get_nbas() {
  return nbas;
}

int CINTPrep::get_nenv() {
  return nenv;
}

int *CINTPrep::get_atm() {
  int *tmp = atm;
  atm = NULL;
  return tmp;
}

int *CINTPrep::get_bas() {
  int *tmp = bas;
  bas = NULL;
  return tmp;
}

double *CINTPrep::get_env() {
  double *tmp = env;
  env = NULL;
  return tmp;
}

void CINTPrep::set_atm(int *&atm_in) {
  if (atm != NULL) {
    delete [] atm;
  }
  atm = atm_in;
  atm_in = NULL;
}

void CINTPrep::set_bas(int *&bas_in) {
  if (bas != NULL) {
    delete [] bas;
  }
  bas = bas_in;
  bas_in = NULL;
}

void CINTPrep::set_env(double *&env_in) {
  if (env != NULL) {
    delete [] env;
  }
  env = env_in;
  env_in = NULL;
}

//RI methods
int CINTPrep::get_var_dim_ri() {
  return var_dim_ri;
}

int CINTPrep::get_nbas_ri() {
  return nbas_ri;
}

int CINTPrep::get_nenv_ri() {
  return nenv_ri;
}
//end ri methods

void CINTPrep::copy_atoms(vector< int > &atoms_copy) {
  vector< int >().swap(atoms_copy);
  int size = atoms.size();
  atoms_copy.reserve(size);
  for (int i = 0; i < size; i++) {
    atoms_copy.push_back(atoms[i]);
  }
}

void CINTPrep::copy_coord(vector< double > &coord_copy) {
  vector< double >().swap(coord_copy);
  int size = coord.size();
  coord_copy.reserve(size);
  for (int i = 0; i < size; i++) {
    coord_copy.push_back(coord[i]);
  }
}

void CINTPrep::assign_coords(int natoms, int *atomlist, double *coords, bool in_bohr) {
  vector< double >().swap(coord);
  vector< int >().swap(atoms);
  atoms.reserve(natoms);
  coord.reserve(natoms*3);

  double scale = (in_bohr) ? 1. : ANG2BOHR;

  for (int i = 0; i < natoms; i++) {
    atoms.push_back(atomlist[i]);
  }
  for (int i = 0; i < natoms * 3; i++) {
    coord.push_back(coords[i]*scale);
  }
}

void CINTPrep::read_xyz(string inpxyz) {
  xyzfile = inpxyz;
  ifstream inpfile(xyzfile.c_str());
  vector< double >().swap(coord);
  vector< int >().swap(atoms);
  int natoms;

  if (inpfile.fail()) {
    printf("ERROR: could not open %s\n", xyzfile.c_str());
  }
  printf("Reading coordinates from: %s\n", xyzfile.c_str());

  string line;
  istringstream iss;
  getline(inpfile, line);
  iss.str(line);

  iss >> natoms;
  if (iss.fail()) {
    printf("ERROR: invalid format for %s\n", xyzfile.c_str());
    exit(1);
  }
  if (natoms <= 0) {
    printf("ERROR: atom count should be > 0\n");
    exit(1);
  }
  iss.clear();
  getline(inpfile, line);

  atoms.reserve(natoms);
  coord.reserve(natoms*3);

  for (int i = 0; i < natoms; i++) {
    string symbol;
    int atomic_num;
    double x, y, z;

    getline(inpfile, line);
    iss.str(line);
    iss >> symbol >> x >> y >> z;
    if (iss.fail()) {
      printf("ERROR: at line %d in file %s\n", i+2, xyzfile.c_str());
      exit(1);
    }
    iss.clear();
    if (elem_2_int.count(symbol) == 0) {
      printf("ERROR: Unknown symbol %s at line %d in file %s\n",
              symbol.c_str(), i, xyzfile.c_str());
      exit(1);
    }

    atomic_num = elem_2_int.at(symbol);
    atoms.push_back(atomic_num);
    coord.push_back(x * ANG2BOHR);
    coord.push_back(y * ANG2BOHR);
    coord.push_back(z * ANG2BOHR);
  } // for i

  if (inpfile.peek() != EOF) {
    printf("ERROR: unexpected lines after last coordinate\n");
    exit(1);
  }
  inpfile.close();
  return;
}

void CINTPrep::read_bas(string inpbas) {
  if (basmap.size() > 0) {
    basmap.clear();
  }
  basfile = inpbas;
  ifstream inpfile(basfile.c_str());
  unordered_set< int > present_atoms;

  if (inpfile.fail()) {
    printf("ERROR: could not open basis file %s\n", basfile.c_str());
    exit(1);
  }

  for (int i = 0; i < atoms.size(); i++) {
    present_atoms.insert(atoms[i]);
  }

  vector< string > lines;
  string line;
  istringstream iss;

  while (getline(inpfile, line)) {
    lines.push_back(line);
  }

  int i = 0;
  while (i < lines.size()) {
    if (!all_of(lines[i].begin(), lines[i].end(), [](char c) {return isspace(c);}) 
        && lines[i].c_str()[0] != '*') {
      iss.str(lines[i]);
      string atom_symbol;
      iss >> atom_symbol;
      if (iss.fail() || elem_2_int.count(atom_symbol) == 0) {
        printf("ERROR reading basis set %s at line %d\n", 
              basfile.c_str(), i);
        exit(1);
      }
      iss.clear();
      i++;

      int atom_num = elem_2_int.at(atom_symbol);
      if (present_atoms.count(atom_num)) {
        present_atoms.erase(atom_num);
        basis_t basis;
        basis.nuc = atom_num;
        while (i < lines.size()) {
          if (!all_of(lines[i].begin(), lines[i].end(), [](char c) {return isspace(c);}) 
              && lines[i].c_str()[0] != '*') {
            iss.str(lines[i]);
            string ang;
            int n_prim;
            iss >> ang >> n_prim;
            int n_ang = ang.length();

            if (iss.fail()) {
              printf("%s\n",lines[i].c_str());
              printf("ERROR reading basis set at line %d\n",i);
              exit(1);
            }

            iss.clear();
            i++;

            for (int j = 0; j < n_ang; j++) {
              int ang_int = angular_mom.at(ang[j]);
              basis.shells.push_back(ang_int);
            } // for j

            vector< double > exps[n_ang];
            vector< double > coefs[n_ang];

            for (int j = 0; j < n_prim; j++) {
              istringstream tmp;
              iss.str(lines[i]);
              string num;
              char skip;
              int power;
              double exp;
              double coeff;

              iss >> num;
              replace(num.begin(),num.end(),'D','e');
              tmp.str(num);
              tmp >> exp;
              tmp.clear();
              // iss >> exp >> skip >> power;
              // exp *= pow(10., power);
              for (int k = 0; k < n_ang; k++) {
                iss >> num;
                replace(num.begin(),num.end(),'D','e');
                tmp.str(num);
                tmp >> coeff;
                tmp.clear();
                // iss >> coeff >> skip >> power;
                // coeff *= pow(10., power);
                exps[k].push_back(exp);
                coefs[k].push_back(coeff);
              } // for k
              iss.clear();
              i++;
            } // for j

            for (int j = 0; j < n_ang; j++) {
              basis.exps.push_back(exps[j]);
              basis.coef.push_back(coefs[j]);
            } // for j

            if (lines[i].c_str()[0] == '*') {
              i++;
              break;
            }
          } // if line.find_first_not_of... != string npos
          else {
            i++;
          } // else
        } // while i < lines.size (inner)

        if (basmap.count(atom_num)) {
          printf("ERROR basis defined twice for atom %s\n", elem_arr[atom_num-1].c_str());
          exit(1);
        }

        basmap[atom_num] = basis;
      } // if present_atoms.count
      else {
        while (i < lines.size()) {
          i++;
          if (lines[i].c_str()[0] == '*') {
            i++;
            break;
          } // if lines[i]
        } // while i < lines.size
      } // else

      iss.clear();
    }
    else {
      i++;
    }
  } // while i < lines.size (outer)
  if (present_atoms.size() != 0) {
    printf("ERROR: Missing basis set for some atoms\n");
    exit(1);
  }
}

void CINTPrep::read_bas_ri(string auxbas) {
  if (basmap_ri.size() > 0) {
    basmap_ri.clear();
  }
  basfile_ri = auxbas;
  ifstream inpfile(basfile_ri.c_str());
  unordered_set< int > present_atoms;

  if (inpfile.fail()) {
    printf("ERROR: could not open aux basis file %s\n", basfile_ri.c_str());
    exit(1);
  }

  for (int i = 0; i < atoms.size(); i++) {
    present_atoms.insert(atoms[i]);
  }

  vector< string > lines;
  string line;
  istringstream iss;

  while (getline(inpfile, line)) {
    lines.push_back(line);
  }

  int i = 0;
  while (i < lines.size()) {
    if (!all_of(lines[i].begin(), lines[i].end(), [](char c) {return isspace(c);}) 
        && lines[i].c_str()[0] != '*') {
      iss.str(lines[i]);
      string atom_symbol;
      iss >> atom_symbol;
      if (iss.fail() || elem_2_int.count(atom_symbol) == 0) {
        printf("ERROR reading aux basis set %s\n", basfile_ri.c_str());
        exit(1);
      }
      iss.clear();
      i++;

      int atom_num = elem_2_int.at(atom_symbol);
      if (present_atoms.count(atom_num)) {
        present_atoms.erase(atom_num);
        basis_t basis;
        basis.nuc = atom_num;
        while (i < lines.size()) {
          if (!all_of(lines[i].begin(), lines[i].end(), [](char c) {return isspace(c);}) 
              && lines[i].c_str()[0] != '*') {
            iss.str(lines[i]);
            string ang;
            int n_prim;
            iss >> ang >> n_prim;
            int n_ang = ang.length();

            if (iss.fail()) {
              printf("ERROR reading basis set\n");
              exit(1);
            }

            iss.clear();
            i++;

            for (int j = 0; j < n_ang; j++) {
              int ang_int = angular_mom.at(ang[j]);
              basis.shells.push_back(ang_int);
              // printf("RI: %s ang: %d\n",ang.c_str(),ang_int);
            } // for j

            vector< double > exps[n_ang];
            vector< double > coefs[n_ang];

            for (int j = 0; j < n_prim; j++) {
              istringstream tmp;
              iss.str(lines[i]);
              string num;
              char skip;
              int power;
              double exp;
              double coeff;

              iss >> num;
              replace(num.begin(),num.end(),'D','e');
              tmp.str(num);
              tmp >> exp;
              tmp.clear();
              // iss >> exp >> skip >> power;
              // exp *= pow(10., power);
              for (int k = 0; k < n_ang; k++) {
                iss >> num;
                replace(num.begin(),num.end(),'D','e');
                tmp.str(num);
                tmp >> coeff;
                tmp.clear();
                // iss >> coeff >> skip >> power;
                // coeff *= pow(10., power);
                exps[k].push_back(exp);
                coefs[k].push_back(coeff);
              } // for k
              iss.clear();
              i++;
            } // for j

            for (int j = 0; j < n_ang; j++) {
              basis.exps.push_back(exps[j]);
              basis.coef.push_back(coefs[j]);
            } // for j

            if (lines[i].c_str()[0] == '*') {
              i++;
              break;
            }
          } // if line.find_first_not_of... != string npos
          else {
            i++;
          } // else
        } // while i < lines.size (inner)

        if (basmap_ri.count(atom_num)) {
          printf("ERROR basis defined twice for atom %s\n", elem_arr[atom_num-1].c_str());
          exit(1);
        }

        basmap_ri[atom_num] = basis;
      } // if present_atoms.count
      else {
        while (i < lines.size()) {
          i++;
          if (lines[i].c_str()[0] == '*') {
            i++;
            break;
          } // if lines[i]
        } // while i < lines.size
      } // else

      iss.clear();
    }
    else {
      i++;
    }
  } // while i < lines.size (outer)
  if (present_atoms.size() != 0) {
    printf("ERROR: Missing basis set for some atoms\n");
    exit(1);
  }
}

void CINTPrep::prep_env() {
  if (atm != NULL) {
    delete [] atm;
  }
  if (bas != NULL) {
    delete [] bas;
  }
  if (env != NULL) {
    delete [] env;
  }

  int natm = atoms.size();
  nbas = 0;
  nenv = PTR_ENV_START;
  int offset = PTR_ENV_START;
  var_dim = 0;

  for (int i = 0; i < natm; i++) {
    int nshls = basmap.at(atoms[i]).shells.size();
    for (int j = 0; j < nshls; j++) {
      nbas++;
    } // for j
  } // for i

  nbas_ri = 0;
  if (do_ri) {
    for (int i = 0; i < natm; i++) {
      int nshls = basmap_ri.at(atoms[i]).shells.size();
      for (int j = 0; j < nshls; j++) {
        nbas_ri++;
      } // for j
    } // for i
  }

  int nbas_all = nbas + nbas_ri;
  if (do_ri) {
    nbas_all++;
  }

  bas = new int[nbas_all * BAS_SLOTS];
  atm = new int[natm * ATM_SLOTS];

  for (int i = 0; i < natm; i++) {
    int atom_num = atoms[i];
    atm[i * ATM_SLOTS + 0] = atom_num;
    atm[i * ATM_SLOTS + 1] = offset;
    atm[i * ATM_SLOTS + 2] = 0;
    atm[i * ATM_SLOTS + 3] = 0;
    offset += 3;
  } // for i

  int bas_num = 0;
  for (int i = 0; i < natm; i++) {
    int atom_num = atoms[i];

    int N_at = 0;
    for (int j = 0; j < basmap.at(atom_num).shells.size(); j++) {
      int nshls = basmap.at(atom_num).shells[j];
      int nprim = basmap.at(atom_num).exps[j].size();
      bas[ATOM_OF  +BAS_SLOTS*bas_num] = i;
      bas[ANG_OF   +BAS_SLOTS*bas_num] = nshls;
      bas[NPRIM_OF +BAS_SLOTS*bas_num] = nprim;
      bas[NCTR_OF  +BAS_SLOTS*bas_num] = 1;
      bas[PTR_EXP  +BAS_SLOTS*bas_num] = offset;
      bas[PTR_COEFF+BAS_SLOTS*bas_num] = offset + nprim;
      offset += nprim*2;
      //var_dim += (BT::DO_CART) ? CINTcgto_cart(bas_num,bas) : CINTcgto_spheric(bas_num,bas);
      N_at += (BT::DO_CART) ? CINTcgto_cart(bas_num,bas) : CINTcgto_spheric(bas_num,bas);
      bas_num++;
    } // for j

    var_dim += N_at;
    anum_to_N[atom_num] = N_at;
  } // for i

  if (do_ri) {
    for (int i = 0; i < natm; i++) {
      int atom_num = atoms[i];
      for (int j = 0; j < basmap_ri.at(atom_num).shells.size(); j++) {
        int nshls = basmap_ri.at(atom_num).shells[j];
        int nprim = basmap_ri.at(atom_num).exps[j].size();
        bas[ATOM_OF  +BAS_SLOTS*bas_num] = i;
        bas[ANG_OF   +BAS_SLOTS*bas_num] = nshls;
        bas[NPRIM_OF +BAS_SLOTS*bas_num] = nprim;
        bas[NCTR_OF  +BAS_SLOTS*bas_num] = 1;
        bas[PTR_EXP  +BAS_SLOTS*bas_num] = offset;
        bas[PTR_COEFF+BAS_SLOTS*bas_num] = offset + nprim;
        offset += nprim*2;
        var_dim_ri += (BT::DO_CART) ? CINTcgto_cart(bas_num,bas) : CINTcgto_spheric(bas_num,bas);
        bas_num++;
      } // for j
    } // for i
    bas[ATOM_OF  +BAS_SLOTS*bas_num] = 0;
    bas[ANG_OF   +BAS_SLOTS*bas_num] = 0;
    bas[NPRIM_OF +BAS_SLOTS*bas_num] = 1;
    bas[NCTR_OF  +BAS_SLOTS*bas_num] = 1;
    bas[PTR_EXP  +BAS_SLOTS*bas_num] = offset;
    bas[PTR_COEFF+BAS_SLOTS*bas_num] = offset + 1;
    offset += 2;
  } // if do_ri

  nenv = offset;
  env = new double[nenv]();
  offset = PTR_ENV_START;
  for (int i = 0; i < natm*3; i++) {
    env[offset++] = coord[i];
  }

  for (int i = 0; i < natm; i++) {
    int atom_num = atoms[i];
    for (int j = 0; j < basmap.at(atom_num).shells.size(); j++) {
      int nexp = basmap.at(atom_num).exps[j].size();
      int ncof = basmap.at(atom_num).coef[j].size();
      int shll = basmap.at(atom_num).shells[j];
      for (int k = 0; k < nexp; k++) {
        env[offset+k] = basmap.at(atom_num).exps[j][k];
      } // for k
      for (int k = 0; k < ncof; k++) {
        env[offset + nexp + k] = basmap.at(atom_num).coef[j][k] * CINTgto_norm(shll, env[offset+k]);
      } // for k
      offset += nexp + ncof;
    } // for j
  } // for i
  if (do_ri) {
    for (int i = 0; i < natm; i++) {
      int atom_num = atoms[i];
      for (int j = 0; j < basmap_ri.at(atom_num).shells.size(); j++) {
        int nexp = basmap_ri.at(atom_num).exps[j].size();
        int ncof = basmap_ri.at(atom_num).coef[j].size();
        int shll = basmap_ri.at(atom_num).shells[j];
        for (int k = 0; k < nexp; k++) {
          env[offset+k] = basmap_ri.at(atom_num).exps[j][k];
        } // for k
        for (int k = 0; k < ncof; k++) {
          env[offset + nexp + k] = basmap_ri.at(atom_num).coef[j][k] * CINTgto_norm(shll, env[offset+k]);
        } // for k
        offset += nexp + ncof;
      } // for j
    } // for i
    env[offset] = 0.;
    env[offset + 1] = 1.;
  } // if do_ri
}
