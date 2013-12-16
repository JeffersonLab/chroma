// $Id: syssolver_qdp_mg_params.cc, v1.0 2012-07-06 16:03 sdcohen $
/*! \file
 *  \brief Parameters for the external QDP clover multigrid solver
 */

#include "actions/ferm/invert/qop_mg/syssolver_qop_mg_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SysSolverQOPMGParams& param)
  {
    XMLReader paramtop(xml, path);

#define defaultread(PARAM,DEFAULT) if (paramtop.count(#PARAM)) {read(paramtop, #PARAM, param.PARAM);} else {param.PARAM = DEFAULT;}

    defaultread( AnisoXi, 1.0 );
    defaultread( AnisoNu, 1.0 );

    if (paramtop.count("Mass") > 0) {
      read(paramtop, "Mass", param.Mass);
      param.Kappa = 0.5/(1.0 + 3.0*param.AnisoNu/param.AnisoXi + param.Mass);
    } else {
      read(paramtop, "Kappa", param.Kappa); // die if not input!
    }   
      
    if (paramtop.count("MassCrit") > 0) {
      read(paramtop, "MassCrit", param.MassCrit);
      param.KappaCrit = 0.5/(1.0 + 3.0*param.AnisoNu/param.AnisoXi + param.MassCrit);
    } else {
      if (paramtop.count("KappaCrit") > 0) {
        read(paramtop, "KappaCrit", param.KappaCrit);
      } else {
        param.KappaCrit = param.Kappa;
      }
    }
    
    defaultread( Clover,     0.0 );
    defaultread( CloverT,    param.Clover );   
    defaultread( Residual,   1e-6 );
    defaultread( MaxIter,    40 );
    defaultread( NumGCRVecs, 8 );
    defaultread( Verbose,    0 );
    defaultread( Levels,     1 );

    if (param.Levels > 0) {
      if (paramtop.count("Blocking")) {
        read(paramtop, "Blocking", param.Blocking);
      } else {
        param.Blocking.resize(param.Levels);
        for (int l=0; l<param.Levels; l++) {
          param.Blocking[l].resize(4);
          param.Blocking[l][0] = param.Blocking[l][1] = param.Blocking[l][2] =
            param.Blocking[l][3] = 4;
        }
      }
    
#define defaultreadvec(PARAM,DEFAULT) if (paramtop.count(#PARAM)) { \
      read(paramtop, #PARAM, param.PARAM); \
    } else { \
      param.PARAM.resize(param.Levels); \
      for (int l=0; l<param.Levels; l++) { \
        param.PARAM[l] = DEFAULT; \
      } \
    }

      defaultreadvec( NumNullVecs, 24 );
      defaultreadvec( NullMaxIter, 100 );
      defaultreadvec( NullResidual, 0.4 );
      defaultreadvec( NullConvergence, 0.5 );
      defaultreadvec( NumExtraVecs, 0 );
      defaultreadvec( Underrelax, 1.0 );
      defaultreadvec( NumPreHits, 0 );
      defaultreadvec( NumPostHits, 4 );
      defaultreadvec( CoarseNumGCRVecs, 8 );
      defaultreadvec( CoarseMaxIter, 12 );
      defaultreadvec( CoarseResidual, 0.2 );

#undef defaultread
#undef defaultreadvec
    }

  }

  // Writer parameters
  void write(XMLWriter& xml, const string& path, const SysSolverQOPMGParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
    write(xml, "invType", "QDP_WILSON_MULTIGRID");
    
#define writeparam(PARAM) write(xml, #PARAM, param.PARAM)
    writeparam(AnisoXi);
    writeparam(AnisoNu);
    writeparam(Kappa);
    writeparam(KappaCrit);
    writeparam(Clover);
    writeparam(CloverT);
    
    writeparam(Residual);
    writeparam(MaxIter);
    writeparam(NumGCRVecs);
    
    writeparam(Verbose);
    
    writeparam(Levels);
    
    if (param.Levels > 0) {
      writeparam(Blocking);
      writeparam(NumNullVecs);
      
      writeparam(NullMaxIter);
      writeparam(NullResidual);
      writeparam(NullConvergence);
      writeparam(NumExtraVecs);
      
      writeparam(Underrelax);
      writeparam(NumPreHits);
      writeparam(NumPostHits);
      writeparam(CoarseNumGCRVecs);
      writeparam(CoarseMaxIter);
      writeparam(CoarseResidual);
    }
#undef writeparam
    pop(xml);
  }

  //! Default constructor
  //SysSolverQOPMGParams::SysSolverQOPMGParams()
  //{
  //}

  //! Read parameters
  SysSolverQOPMGParams::SysSolverQOPMGParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
