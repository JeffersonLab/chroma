// =================================================================
// QDP Clover Multigrid Chroma Interface
//   Precision-dependent header
//     This header defines the master clover parameter object

//#ifndef __wilmg_h__
//#define __wilmg_h__

#ifndef __cplusplus
#include <qop.h>
#include <math.h>

QDP_Layout *PC(QDP_layout_hyper_eo2);

// Utility functions for changing QDP gauge field boundary conditions
static int bcdir;
static QLA_Real bcphase;
static void change_bc(QLA_ColorMatrix *cm, int coords[]) {
  if (coords[bcdir] == QDP_coord_size(bcdir)-1)
    QLA_M_eq_r_times_M(cm, &bcphase, cm);
}

#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MAXLEVELS 5
// Since we want the structure here to have a definite size, we want the arrays
// to be of known length. Given that the deepest multigrid in current use has
// 3 levels, a limit of 5 seems safe. The value can be changed on recompile.

struct MGP(Clover_Params) {
  // Lattice Action Parameters
    double bc[4];     // Boundary conditions (+1 periodic, -1 antiperiodic)
    QLA_Real aniso_xi;// Lattice bare anisotropy (xi_0)
    QLA_Real aniso_nu;// Lattice bare dispersion parameter (nu_s)
    QLA_Real kappa;   // Hopping parameter to solve (1/(2(m+1+3nu/xi))
    QLA_Real kappac;  // Critical hopping parameter (for null vectors)
    QLA_Real mass;    // Bare mass of fermion (sets kappa)
    QLA_Real massc;   // Bare critical mass (sets kappac)
    QLA_Real clov_s;  // Spatial clover parameter
    QLA_Real clov_t;  // Temporal clover parameter
  // Solver Parameters
    QLA_Real res;     // Stopping residual for solver
    int      maxiter; // Maximum number of iterations to allow in solver
    int      ngcr;    // Number of GCR vectors
  // Diagnostic Parameters
    int      verb;    // Level of diagnostic verbosity
  // Multigrid Parameters
    int      levels;  // Number of levels in multigrid
    int block[MAXLEVELS][4];  // Spacetime blocking of each multigrid level
    int nNullVecs[MAXLEVELS]; // Number of null vectors per multigrid level
  // Null-Vector Setup Parameters
    int nullMaxIter[MAXLEVELS]; // Maximum iterations for setup on each vector
    QLA_Real nullRes[MAXLEVELS]; // Residual to converge each vector
    QLA_Real nullConv[MAXLEVELS]; // Convergence criterion
      // This indicates the level at which a vector is considered to have
      // converged. That is, if it changes less than this amount during the
      // relaxation, a new random vector will be generated for further nullvecs.
    int nExtraVecs[MAXLEVELS]; // Number of extra vectors to generate and discard
  // Multigrid Solver Parameters
    QLA_Real urelax[MAXLEVELS]; // Underrelaxation for each V-cycle
    int npre[MAXLEVELS]; // Number of smoother pre-hits per V-cycle
    int npost[MAXLEVELS]; // Number of smoother post-hits per V-cycle
    int cngcr[MAXLEVELS]; // Number of GCR vectors in coarse-level solvers
    // Solver will run until one of the following two stopping criteria is met
    int cmaxiter[MAXLEVELS]; // Coarse-level maximum number of iterations
    QLA_Real cres[MAXLEVELS]; // Coarse-level relative stopping residual
  };
 

 
void MGP(initialize)( int *machsize, int *latsize,
                      void (*peekpoke[4])(QLA(ColorMatrix) *dest, int coords[]) );

  void*  MGP(create_subspace)( int *latsize );
  void MGP(reset_subspace)(int *latsize, void *subspace_in);
  void MGP(destroy_subspace)( void *subspace  );

int MGP(solve)( void peekpokesrc(QLA(DiracFermion) *dest, int coords[]),
		void peekpokeguess(QLA(DiracFermion) *dest, int coords[]),
                void peekpokesol(QLA(DiracFermion) *src,  int coords[]),
		void *subspace );

void MGP(finalize)();
void MGP(teststuff)();
 
#ifdef __cplusplus
}
#endif

//#endif 
