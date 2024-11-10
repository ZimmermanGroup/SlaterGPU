  bool sgs_basis = 0;
  bool ss_basis = 0;
  if (basis[0].size()>10) ss_basis = 1;
  //if (basis[0].size()>10) sgs_basis = 1;
  if (sgs_basis && prl>1) printf("  this is an SGS basis \n");
  if (ss_basis && prl>1)  printf("  this is an SS basis \n");

  double Zg = read_float("ZG");
  double ztg = read_float("ZTG");
  double Rc = read_float("RC");
  //int jorder = 2;

  int gs1 = 0;
  if (Rc>0.)
    gs1 = get_gs1(nrad,nang,grid,Rc);

