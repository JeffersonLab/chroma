// -*- C++ -*-
// $Id: invtype.h,v 1.7 2004-09-22 17:24:59 bjoo Exp $

/*! @file
 * @brief Types of inverters
 */

#ifndef __invtype_h__
#define __invtype_h__


namespace Chroma
{

  //! Types of inverters
  enum InvType {
    CG_INVERTER = 21,
    MR_INVERTER = 22,
    BICG_INVERTER = 23,
    CR_INVERTER = 24,
    SUMR_INVERTER = 25,
    REL_CG_INVERTER = 26,
    REL_SUMR_INVERTER = 27,
    REL_GMRESR_SUMR_INVERTER = 28,
    REL_GMRESR_CG_INVERTER = 29,
    BICGSTAB_INVERTER = 30
  };


  //! Parameters for inverter
  struct InvertParam_t
  {
    InvType       invType;   // Inverter type
    Real          MROver;
    Real          RsdCG;
    Real          RsdCGPrec; // Precision for preconditioned iterations
    int           MaxCG;	   // Iteration parameters
    int           MaxCGPrec;
  };

  //! Parameters for MultiMass Inverter
  struct MultiInvertParam_t
  {
    InvType      invType;
    Real         MROver;
    int          MaxCG;
    multi1d<Real> RsdCG;
  };

}


using namespace Chroma;

#endif
