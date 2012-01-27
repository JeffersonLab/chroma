// -*- C++ -*-
// $Id: wloop.h,v 1.1 2008-06-26 14:58:35 mcneile Exp $

#ifndef __wilson_flow_w_h__
#define __wilson_flow_w_h__

namespace Chroma 
{


  //! Print the value of the average plaquette normalized to 1
  /*!
   * \ingroup glue
   *
   * \param xml    plaquette average (Write)
   * \param u      gauge field      (Read)
   * \param nstep  number of steps  (Read)
   * \param wflow_eps  size of step (Read)
   * \param time direction (Read)

   */
  void wilson_flow(multi1d<LatticeColorMatrix> & u, int nstep, 
		   Real  wflow_eps, int jomit) ;

}  // end namespace Chroma

#endif
