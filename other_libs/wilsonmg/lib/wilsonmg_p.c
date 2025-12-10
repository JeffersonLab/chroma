// =================================================================
// QDP Clover Multigrid Chroma Interface
//   Precision-dependent code

// -----------------------------------------------------------------
// Includes
#define USE_MG

#if QDP_Precision == 'F' || QDP_Precision == 1
#include "generic_f.h"
#else
#include "generic_d.h"
#endif

#include "wilsonmg_p.h"
#define printf0 if(QDP_this_node==0) printf

// The Dirac operator, multigrid structure and QOP structures are global
static QOP(FermionLinksWilson) *wil = NULL;
static QOP_F3_FermionLinksWilson *wilf = NULL;
			       //static QOP_WilsonMg *wilmg = NULL;
static QOP_info_t info = {0,0,QOP_SUCCESS,0,0};
static QOP_invert_arg_t inv = {40,40,1,QOP_EVENODD};
static QOP_resid_arg_t res = {1e-12,0,0,0,0,0};

struct MGP(Clover_Params) PC(g_param);
// -----------------------------------------------------------------
// Takes a global average (only used for profiling)
static double
global_average(double x)
{
  QMP_sum_double(&x);
  return x/QMP_get_number_of_nodes();
}

void* MGP(create_subspace)(int *latsize)
{
  QOP_WilsonMg* ret_val = NULL;
  double timer=0;
  int ndim=4;
  printf0("QDP: Creating multigrid structure\n");
  timer = -QDP_time();
  ret_val = QOP_wilsonMgNew();
  QOP_wilsonMgSetLinks(ret_val, wilf);
  // Set meta-global parameter
  // Set global multigrid parameters
  QOP_wilsonMgSet(ret_val, -1, "nlevels", PC(g_param).levels);
  QOP_wilsonMgSet(ret_val, -2, "verbose", PC(g_param).verb);
  QOP_wilsonMgSet(ret_val, -1, "profile", 1);
  QOP_wilsonMgSet(ret_val, -1, "kappanv", PC(g_param).kappac);
  QOP_wilsonMgSet(ret_val, -1, "kappa",   PC(g_param).kappac);
  QOP_wilsonMgSet(ret_val, -1, "itmax",   PC(g_param).maxiter);
  QOP_wilsonMgSet(ret_val, -1, "ngcr",    PC(g_param).ngcr);

  for (int l=0; l<PC(g_param).levels; l++) {
    printf0("QDP:   Creating %i%s multigrid level with size",l+1,
              (l>2?"th":l>1?"rd":l>0?"nd":"st"));
    // Set lth-level parameters
    double *dlat = malloc(ndim*sizeof(double));
    for (int d=0; d<ndim; d++) {
      dlat[d] = latsize[d];
      for (int ll=0; ll<=l; ll++) dlat[d] /= PC(g_param).block[ll][d];
      printf0(" %g",dlat[d]);
    }
    printf0("\n");
    QOP_wilsonMgSetArray(ret_val, l, "lattice", dlat, ndim);
    QOP_wilsonMgSet(ret_val, l, "nvecs", PC(g_param).nNullVecs[l]);
    // Set lth-level parameters (null-vector generation)
    QOP_wilsonMgSet(ret_val, l, "setup_nvecs", PC(g_param).nNullVecs[l] + PC(g_param).nExtraVecs[l]);
    QOP_wilsonMgSet(ret_val, l, "setup_res", PC(g_param).nullRes[l]);
    QOP_wilsonMgSet(ret_val, l, "setup_maxit", PC(g_param).nullMaxIter[l]);
    QOP_wilsonMgSet(ret_val, l, "setup_change_fac", PC(g_param).nullConv[l]);
    // Set lth-level parameters (V-cycle)
    QOP_wilsonMgSet(ret_val, l, "npre", PC(g_param).npre[l]);
    QOP_wilsonMgSet(ret_val, l, "npost", PC(g_param).npost[l]);
    QOP_wilsonMgSet(ret_val, l, "scale", PC(g_param).urelax[l]);
    // Set lth-level parameters (solver)
    QOP_wilsonMgSet(ret_val, l, "cres", PC(g_param).cres[l]);
    QOP_wilsonMgSet(ret_val, l, "itmax", PC(g_param).cmaxiter[l]);
    QOP_wilsonMgSet(ret_val, l, "ngcr", PC(g_param).cngcr[l]);
  }

  // Testing. Dont set up the MG 
  printf0("QDP:   Solving null vectors for multigrid\n");
  QOP_wilsonMgSetup(ret_val);
  timer += QDP_time(); QMP_max_double(&timer);
  printf0("QDP: multigrid structure setup time = %g secs\n", timer);
  return (void *)ret_val;
}

void MGP(reset_subspace)(int *latsize, void *subspace_in)
{
  QOP_WilsonMg* ret_val = (QOP_WilsonMg*)(subspace_in);
  double timer=0;
  int ndim=4;
  printf0("QDP: Resetting multigrid structure\n");
  timer = -QDP_time();
  QOP_wilsonMgSetLinks(ret_val, wilf);
  // Set meta-global parameter
  // Set global multigrid parameters
  QOP_wilsonMgSet(ret_val, -1, "nlevels", PC(g_param).levels);
  QOP_wilsonMgSet(ret_val, -2, "verbose", PC(g_param).verb);
  QOP_wilsonMgSet(ret_val, -1, "profile", 1);
  QOP_wilsonMgSet(ret_val, -1, "kappanv", PC(g_param).kappac);
  QOP_wilsonMgSet(ret_val, -1, "kappa",   PC(g_param).kappac);
  QOP_wilsonMgSet(ret_val, -1, "itmax",   PC(g_param).maxiter);
  QOP_wilsonMgSet(ret_val, -1, "ngcr",    PC(g_param).ngcr);

  for (int l=0; l<PC(g_param).levels; l++) {
    printf0("QDP:   Creating %i%s multigrid level with size",l+1,
              (l>2?"th":l>1?"rd":l>0?"nd":"st"));
    // Set lth-level parameters
    double *dlat = malloc(ndim*sizeof(double));
    for (int d=0; d<ndim; d++) {
      dlat[d] = latsize[d];
      for (int ll=0; ll<=l; ll++) dlat[d] /= PC(g_param).block[ll][d];
      printf0(" %g",dlat[d]);
    }
    printf0("\n");
    QOP_wilsonMgSetArray(ret_val, l, "lattice", dlat, ndim);
    QOP_wilsonMgSet(ret_val, l, "nvecs", PC(g_param).nNullVecs[l]);
    // Set lth-level parameters (null-vector generation)
    QOP_wilsonMgSet(ret_val, l, "setup_nvecs", PC(g_param).nNullVecs[l] + PC(g_param).nExtraVecs[l]);
    QOP_wilsonMgSet(ret_val, l, "setup_res", PC(g_param).nullRes[l]);
    QOP_wilsonMgSet(ret_val, l, "setup_maxit", PC(g_param).nullMaxIter[l]);
    QOP_wilsonMgSet(ret_val, l, "setup_change_fac", PC(g_param).nullConv[l]);
    // Set lth-level parameters (V-cycle)
    QOP_wilsonMgSet(ret_val, l, "npre", PC(g_param).npre[l]);
    QOP_wilsonMgSet(ret_val, l, "npost", PC(g_param).npost[l]);
    QOP_wilsonMgSet(ret_val, l, "scale", PC(g_param).urelax[l]);
    // Set lth-level parameters (solver)
    QOP_wilsonMgSet(ret_val, l, "cres", PC(g_param).cres[l]);
    QOP_wilsonMgSet(ret_val, l, "itmax", PC(g_param).cmaxiter[l]);
    QOP_wilsonMgSet(ret_val, l, "ngcr", PC(g_param).cngcr[l]);
  }

  timer += QDP_time(); QMP_max_double(&timer);
  printf0("QDP: multigrid structure setup time = %g secs\n", timer);
}

void MGP(destroy_subspace)(void *subspace)
{

  QOP_WilsonMg *qop_wilsmg = (QOP_WilsonMg *)subspace;  
  if (qop_wilsmg) {
    printf0("**** Destroying subspae at %lx\n", subspace);
    QOP_wilsonMgFree(qop_wilsmg);
  }
}

// -----------------------------------------------------------------
// Initialize the QDP environment
void MGP(initialize)( int *machsize, int *latsize, 
                     void (*peekpoke[4])(QLA(ColorMatrix) *dest, int coords[]) )
{
  double timer;
  int ndim = 4;
  // Initialize the QOP layout and basics
  if (!QDP_is_initialized()) {
    timer = -QDP_time();
    printf0("QDP: Initializing the QDP software stack!\n");
    
    int qmp_argc = ndim+2;
    char **qmp_argv = malloc(qmp_argc*sizeof(char*));
    for (int i=0; i<qmp_argc; i++) qmp_argv[i] = malloc(10*sizeof(char));
    sprintf(qmp_argv[0],"wilsonmg");
    sprintf(qmp_argv[1],"-qmp-geom");
    for (int i=0; i<ndim; i++) sprintf(qmp_argv[i+2],"%i",machsize[i]);
    // These pointers are usually &argc and &argv, passed down through QDP into
    // QMP_init_msg_passing; we try to reconstruct them from the machine layout.
    QDP_initialize(&qmp_argc, &qmp_argv);
    QDP_profcontrol(0);
    QDP_set_latsize(ndim, latsize);
    //QDP_set_default_layout(PC(QDP_layout_hyper_eo2));
    QDP_create_layout();
    
    QOP_layout_t qoplayout = QOP_LAYOUT_ZERO;
    qoplayout.latdim = ndim;
    qoplayout.latsize = (int*) malloc(ndim*sizeof(int));
    for(int i=0; i<ndim; i++) {
      qoplayout.latsize[i] = latsize[i];
    }
    qoplayout.machdim = -1;
    QOP_init(&qoplayout);
    timer += QDP_time(); QMP_max_double(&timer);
    printf0("QDP: stack initialization time = %g secs\n", timer);
  }
  
  // Initialize the gauge field
  timer = -QDP_time();
  printf0("QDP: Loading gauge field from Chroma\n");
  QDP_ColorMatrix **gauge =
    (QDP_ColorMatrix**)malloc(ndim*sizeof(QDP_ColorMatrix*));
  for (int d=0; d<ndim; d++) {
    gauge[d] = QDP_create_M();
    QDP_M_eq_func(gauge[d], peekpoke[d], QDP_all);
  }
    
  // Check the average plaquette of the configuration
  QLA_Real tsum, plaq, splaq=0.0, tplaq=0.0;
  QDP_ColorMatrix *temp1, *temp2, *temp3, *temp4;
  temp1 = QDP_create_M(); temp2 = QDP_create_M();
  temp3 = QDP_create_M(); temp4 = QDP_create_M();
  for (int mu=0; mu<ndim-1; ++mu) {
    for (int nu=mu+1; nu<ndim; ++nu) {
      QDP_M_eq_sM(temp1, gauge[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
      QDP_M_eq_sM(temp2, gauge[mu], QDP_neighbor[nu], QDP_forward, QDP_all);

      QDP_M_eq_Ma_times_M(temp3, gauge[nu], gauge[mu], QDP_all);

      QDP_M_eq_M_times_M(temp4, temp3, temp1, QDP_all);
      QDP_discard_M(temp1);

      QDP_r_eq_re_M_dot_M(&tsum, temp2, temp4, QDP_all);
      QDP_discard_M(temp2);
      if (nu==ndim-1) tplaq += tsum;
      else splaq += tsum;
    }
  }
  QDP_destroy_M(temp1); QDP_destroy_M(temp2);
  QDP_destroy_M(temp3); QDP_destroy_M(temp4);

  splaq /= (0.5*(ndim-1)*(ndim-2)*QDP_volume()*QLA_Nc);
  tplaq /= ((ndim-1)*QDP_volume()*QLA_Nc);
  plaq = ((ndim-2)*splaq + 2*tplaq)/ndim;
  
  printf0("QDP:   Average plaquette = %g\n", plaq);
  printf0("QDP:   Average spatial plaquette = %g\n", splaq);
  printf0("QDP:   Average temporal plaquette = %g\n", tplaq);

  // Set the boundary conditions -- All periodic -- Chroma 
  // sets boundaries now through ferm-state. Non gauge link bcs are 
  // not really allowed.
#if 0
  for (int d=0; d<ndim; d++) {
    if (PC(g_param).bc[d]!=1) {
      bcdir = d;
      bcphase = PC(g_param).bc[d];
      QDP_M_eq_func(gauge[d], change_bc, QDP_all);
    }
  }
#endif 

  timer += QDP_time(); QMP_max_double(&timer);
  printf0("QDP: load gauge field time = %g secs\n", timer);
  
  
  // Initialize the Wilson clover operator
  QOP(GaugeField) *qopgauge = QOP(create_G_from_qdp)(gauge);
  QOP_wilson_coeffs_t coeffs;
  coeffs.clov_s = PC(g_param).clov_s/PC(g_param).aniso_xi;
  coeffs.clov_t = PC(g_param).clov_t;
  // Flip it.
  coeffs.aniso = PC(g_param).aniso_nu/PC(g_param).aniso_xi;

  timer = -QDP_time();
  printf0("QDP: Creating Wilson-clover Dirac operator\n");
  wil = QOP(wilson_create_L_from_G)(&info, &coeffs, qopgauge);

  // After the gauge fields are safely inside the ops, we can delete them
  //for (int d=0; d<ndim; d++) QDP_destroy_M(gauge[d]);
  //free(gauge);

  // Initialize the single-precision version of the operator  
#if QDP_Precision == 1 || QDP_Precision == 'F'
  wilf = wil;
#else
  printf0("QDP:   Creating single-precision operator\n");
  wilf = QOP_FD3_wilson_create_L_from_L(wil);
#endif
  timer += QDP_time(); QMP_max_double(&timer);
  printf0("QDP: Dirac-op creation time = %g secs\n", timer);
 

  //  wilmg = create_subspace(latsize);




 
 
}

// -----------------------------------------------------------------
// Solve for a fermion source specified by the peekpoke function
int MGP(solve)( void peekpokesrc(QLA(DiracFermion) *dest, int coords[]),
		void peekpokeguess(QLA(DiracFermion) *dest, int coords[]),
		void peekpokesol(QLA(DiracFermion) *src,  int coords[]),
		void *subspace_in)
{
  QOP_WilsonMg* subspace = (QOP_WilsonMg *)subspace_in;
  if (subspace == NULL || subspace == 0x0 ) { 
    printf0("ERROR: MGP(solve) called with NULL subspace\n");
    return -1;
  }
  printf0("Using subspace with ptr: %lx\n", subspace);

  double timer;
  int its=0;
  QLA_Real bsq;
  // Initialize some fermion fields
  printf0("QDP: Copying fermion source from Chroma\n");
  timer = -QDP_time(); // Start timing source transfer
  QDP_DiracFermion *in = QDP_create_D();
  QDP_DiracFermion *out = QDP_create_D();
  //  QDP_D_eq_zero(out, QDP_all);
  QDP_D_eq_func(out,peekpokeguess, QDP_all);
  QDP_r_eq_norm2_D(&bsq, out, QDP_all);
  printf0("QDP:   Chroma initial guess norm2 = %g\n", bsq);
  
  QDP_D_eq_func(in, peekpokesrc, QDP_all);
  QDP_r_eq_norm2_D(&bsq, in, QDP_all);
  printf0("QDP:   Chroma in norm2 = %g\n", bsq);
  timer += QDP_time(); QMP_max_double(&timer); // End timing source transfer
  printf0("QDP: Source transfer took %g secs\n", timer);
 
#if 1 
  // Load the QOP structures with the appropriate parameters
  inv.max_iter = PC(g_param).maxiter;
  res.rsqmin = PC(g_param).res*PC(g_param).res;
  
  printf0("QDP: Apply the multigrid inverter!\n");
  timer = -QDP_time(); // Start timing inversion
  // Kludge: QOP_wilsonMgSolve does not seem to be properly registered for D/F
  #if QDP_Precision == 1 || QDP_Precision == 'F'
  QOP_F3_wilsonMgSolve(&info, subspace, wil, &inv, &res, PC(g_param).kappa, out, in);
  #else
  QOP_D3_wilsonMgSolve(&info, subspace, wil, &inv, &res, PC(g_param).kappa, out, in);
  #endif

  timer += QDP_time(); QMP_max_double(&timer); // End timing inversion
  its += res.final_iter;

  printf0("QDP: Inversion took %i iterations and %g secs\n", its, timer);
#else 
  printf0("Testing QDP Dirac operator\n");
  QOP_D3_wilson_dslash_qdp(&info, wil, PC(g_param).kappa, +1, out, in, QOP_EVENODD, QOP_EVENODD); 

#endif

  printf0("QDP: Copying fermion solution back to Chroma\n");
  timer = -QDP_time(); // Start timing source transfer
  QDP_r_eq_norm2_D(&bsq, out, QDP_all);
  printf0("QDP:   Chroma out norm2 = %g\n", bsq);
  QDP_D_eq_func(out, peekpokesol, QDP_all);
  timer += QDP_time(); QMP_max_double(&timer); // End timing solution transfer
  printf0("QDP: Solution transfer took %g secs\n", timer);

  QDP_destroy_D(in); QDP_destroy_D(out);
  return its;
}

// -----------------------------------------------------------------
// Close down the QDP environment
void MGP(finalize)() {
  printf0("QDP: Finalizing\n");
  //  if (wilmg) QOP_wilsonMgFree(wilmg);
  if (wil) QOP(wilson_destroy_L)(wil);
  if (wilf && wilf!=wil) QOP_F3_wilson_destroy_L(wilf);
  //QDP_finalize();
}

void MGP(teststuff)() {
#define test(PARAM) printf0("I'm QDP and I see that %s is",#PARAM); for (int d=0;d<4;d++) printf0(" %g",PC(g_param).PARAM[d]); printf0("!\n");
        test(bc);
#undef test
#define test(PARAM) printf0("I'm QDP and I see that %s is %g!\n",#PARAM,PC(g_param).PARAM);
        test(aniso_xi);test(aniso_nu);test(kappa);test(kappac);test(mass);test(massc);test(clov_s);test(clov_t);test(res);
#undef test
#define test(PARAM) printf0("I'm QDP and I see that %s is %i!\n",#PARAM,PC(g_param).PARAM);
        test(maxiter);test(verb);test(ngcr);test(levels);
#undef test
#define test(PARAM) printf0("I'm QDP and I see that %s is",#PARAM); for (int n=0; n<PC(g_param).levels; n++) {for (int d=0;d<4;d++) printf0(" %i",PC(g_param).PARAM[n][d]); if (n<PC(g_param).levels-1) printf0(",");} printf0("!\n");
        test(block);
#undef test
#define test(PARAM) printf0("I'm QDP and I see that %s is",#PARAM); for (int n=0; n<PC(g_param).levels; n++) printf0(" %i",PC(g_param).PARAM[n]); printf0("!\n");
        test(nNullVecs);test(nullMaxIter);test(nExtraVecs);test(npre);test(npost);test(cngcr);test(cmaxiter);
#undef test
#define test(PARAM) printf0("I'm QDP and I see that %s is",#PARAM); for (int n=0; n<PC(g_param).levels; n++) printf0(" %g",PC(g_param).PARAM[n]); printf0("!\n");
        test(nullRes);test(nullConv);test(urelax);test(cres);
#undef test
}

