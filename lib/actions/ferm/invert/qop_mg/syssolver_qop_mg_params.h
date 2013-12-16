// -*- C++ -*-
// $Id: syssolver_qop_mg_params.h, v1.0 2012-07-06 16:00 sdcohen $
/*! \file
 *  \brief Parameters for the external QDP clover multigrid solver
 */

#ifndef __syssolver_qop_mg_params_h__
#define __syssolver_qop_mg_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Parameters for the external QDP multigrid inverter
  /*! \ingroup invert */
  struct SysSolverQOPMGParams
  {
    SysSolverQOPMGParams() {};
    SysSolverQOPMGParams(XMLReader& in, const std::string& path);
    
    /* These parameters are translated directly from the sample code in wilmg.c
     * in qcdalg. The comments give the command-line flag used to input the
     * parameters in compiled QDP, followed by the C variable name used,
     * follwed by a description of how the value is used. We will try to use
     * the default values from wilmg.c to initialize these parameters.
     */
    
    // Lattice Action Parameters
    Real AnisoXi;    /*!< Lattice bare anisotropy (xi_0) */
    Real AnisoNu;    /*!< Lattice bare dispersion parameter (nu_s) */
    Real Kappa;      /*!< Hopping parameter to solve */
    Real KappaCrit;  /*!< Critical hopping parameter (for null vectors) */
    Real Mass;       /*!< Bare mass of fermion (sets kappa) */
    Real MassCrit;   /*!< Bare critical mass (sets kappac) */
    Real Clover;     /*!< Spatial clover parameter */
    Real CloverT;    /*!< Temporal clover parameter */
  // Solver Parameters
    Real Residual;   /*!< Stopping residual for solver */
    int  MaxIter;    /*!< Maximum number of iterations to allow in solver */
    int  NumGCRVecs; /*!< Number of GCR vectors at top level */
  // Diagnostic Parameters
    int  Verbose;    /*!< Level of diagnostic verbosity */
  // Multigrid Parameters
    int  Levels;     /*!< Number of levels in multigrid 
      If Levels is specified to be nonpositive, the same multigrid setup
      previously created will be reused by this inversion.
      If Levels is negative, the multigrid structure will be deleted following
      the inversion. */
    multi1d< multi1d<int> > Blocking;/*!< Spacetime blocking of each multigrid level */
    multi1d<int>  NumNullVecs;     /*!< Number of null vectors per multigrid level */
  // Null-Vector Setup Parameters
    multi1d<int>  NullMaxIter;     /*!< Maximum iterations for setup on each vector */
    multi1d<Real> NullResidual;    /*!< Residual to converge each vector */
    multi1d<Real> NullConvergence; /*!< Convergence criterion
      This indicates the level at which a vector is considered to have
      converged. That is, if it changes less than this amount during the
      relaxation, a new random vector will be generated for further nullvecs. */
    multi1d<int>  NumExtraVecs;    /*!< Number of extra vectors to generate and discard */
  // Multigrid Solver Parameters
    multi1d<Real> Underrelax;      /*!< Underrelaxation for each V-cycle */
    multi1d<int>  NumPreHits;      /*!< Number of smoother pre-hits per V-cycle */
    multi1d<int>  NumPostHits;     /*!< Number of smoother post-hits per V-cycle */
    multi1d<int>  CoarseNumGCRVecs;/*!< Number of GCR vectors in coarse solver */
    // Solver will run until one of the following two stopping criteria is met
    multi1d<int>  CoarseMaxIter;   /*!< Coarse-level maximum number of iterations */
    multi1d<Real> CoarseResidual;  /*!< Coarse-level relative stopping residual */
  };


  // Reader and writer
  /*! \ingroup invert */
  void read(XMLReader& xml, const string& path, SysSolverQOPMGParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const string& path,
             const SysSolverQOPMGParams& param);

} // End namespace

#endif 

