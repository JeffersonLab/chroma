// -*- C++ -*-
#ifndef __wilson_flow_w_h__
#define __wilson_flow_w_h__

#include "chromabase.h"

namespace Chroma 
{


  //! Compute the Wilson flow
  /*!
   * \ingroup glue
   *
   * \param xml    wilson flow (Write)
   * \param u      gauge field      (Read)
   * \param nstep  number of steps  (Read)
   * \param wflow_eps  size of step (Read)
   * \param time direction (Read)

   */

  void wilson_flow(XMLWriter& xml,
		   multi1d<LatticeColorMatrix> & u, int nstep, 
		   Real  wflow_eps, int jomit)  ;


}  // end namespace Chroma

#endif
