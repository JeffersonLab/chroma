// -*- C++ -*-
/*! \file
 *  \brief Include possibly optimized Clover terms
 */

#ifndef __exp_clover_term_w_h__
#define __exp_clover_term_w_h__

#include "chroma_config.h"
#include "qdp_config.h"

#include "actions/ferm/linop/exp_clover_term_qdp_w.h"

// The QDP naive clover term
//
namespace Chroma
{

  using ExpCloverTerm = QDPExpCloverTerm<>;
  using ExpCloverTermF = QDPExpCloverTermF<>;
  using ExpCloverTermD = QDPExpCloverTermD<>;

}
#endif
