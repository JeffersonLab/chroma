// -*- C++ -*-
// $Id: dwf_qpropt_w.h,v 3.12 2008-05-05 19:11:58 bjoo Exp $
/*! \file
 * \brief Pick up possibly optimized DWF inverters.
 *
 * Typedefs to pick up possibly optimized DWF inverters
 */

#ifndef DWF_QPROPT_W_H
#define DWF_QPROPT_W_H

#include "chroma_config.h"

// The QDP naive class: PrecFermAct5DQprop
#include "eoprec_fermact_qprop_array.h"


// The following is an ifdef lis that switches in optimised
// Dslash-es. Currently only optimised inverters are the SSE and ALTICEC 
#if defined(BUILD_CG_DWF)
#include "actions/ferm/qprop/avp_inverter_interface.h"

// Internal configuration
#include <cg-dwf-config.h>

#ifdef CG_DWF_ARCH_SSE
// SSE architecture
#include "actions/ferm/qprop/avp_ssef_solver.h"
#include "actions/ferm/qprop/avp_ssed_solver.h"

// Enable both single and a double prec solver
#define SINGLE_PREC_SOLVER
#define DOUBLE_PREC_SOLVER


// Include the qprop file
#ifdef CHROMA_USE_CG_DWF_LOWMEM
#include "actions/ferm/qprop/eoprec_dwf_qprop_array_cg_dwf_lowmem_w.h"
#else
#include "actions/ferm/qprop/eoprec_dwf_qprop_array_cg_dwf_w.h"
#endif

// Define types 
namespace Chroma { 
  typedef Chroma::CGDWFQpropT< AVPSolver::SSEDWFSolverF, AVPSolver::SSEDWFSolverD>  DWFQpropT;
}
#elif defined CG_DWF_ARCH_BLUELIGHT // CG_DWF_ARCH_SSE
#include "actions/ferm/qprop/avp_bglf_solver.h"
#include "actions/ferm/qprop/avp_bgld_solver.h"
// Enable both single and a double prec solver
#define SINGLE_PREC_SOLVER
#define DOUBLE_PREC_SOLVER

// Include the qprop file
#ifdef CHROMA_USE_CG_DWF_LOWMEM
#include "actions/ferm/qprop/eoprec_dwf_qprop_array_cg_dwf_lowmem_w.h"
#else
#include "actions/ferm/qprop/eoprec_dwf_qprop_array_cg_dwf_w.h"
#endif

namespace Chroma { 
  typedef Chroma::CGDWFQpropT< AVPSolver::BGLDWFSolverF, AVPSolver::BGLDWFSolverD>  DWFQpropT;
}
#elif defined CG_DWF_ARCH_ALTIVEC
#include "actions/ferm/qprop/avp_altivecf_solver.h"
// Enable only single prec solver
#define SINGLE_PREC_SOLVER
#ifdef CHROMA_USE_CG_DWF_LOWMEM
#include "actions/ferm/qprop/eoprec_dwf_qprop_array_cg_dwf_lowmem_w.h"
#else
#include "actions/ferm/qprop/eoprec_dwf_qprop_array_cg_dwf_w.h"
#endif
// To satisfy type requirements pass the type of single prec 
// as the double prec solver. However this will be turned off by
// the fact that we have not defined DOUBLE_PREC_SOLVER
namespace Chroma {
  typedef Chroma::CGDWFQpropT< AVPSolver::AltiVecDWFSolverF, AVPSolver::AltiVecDWFSolverF>  DWFQpropT;
}
#endif

#elif defined(BUILD_SSE_DWF_CG)
// The file defines the SSE Dslash class
// The following typedef switches it in.
#include "eoprec_dwf_qprop_array_sse_w.h"
namespace Chroma {
typedef SSEDWFQpropT DWFQpropT;
}  // end namespace Chroma

#elif defined(BUILD_ALTIVEC_DWF_CG)
// The file defines the ALTIVEC Dslash class
// The following typedef switches it in.
#include "eoprec_dwf_qprop_array_altivec_w.h"
namespace Chroma {
typedef ALTIVECDWFQpropT DWFQpropT;
}  // end namespace Chroma

#elif defined(BUILD_MDWF)

// OK Set up Andrews new MDWF solver here.
#include "mdwf_solver.h"
namespace Chroma { 
  typedef MDWFQpropT DWFQpropT;
}


#else
// Bottom line, if no optimised DWF qpropT-s exist then the naive QDP qpropT
// becomes the DWFQpropT
namespace Chroma {
typedef PrecFermAct5DQprop<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > DWFQpropT;
}  // end namespace Chroma
#endif


#endif
