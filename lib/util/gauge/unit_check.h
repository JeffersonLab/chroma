// -*- C++ -*-
// $Id: unit_check.h,v 3.1 2006-08-25 23:46:37 edwards Exp $

#ifndef UNIT_CHECK_INCLUDE
#define UNIT_CHECK_INCLUDE

namespace Chroma 
{
  //! Check the unitarity of color matrix in SU(N)
  /*!
   * \ingroup gauge
   *
   * \param  u  The multi1d<LatticeColorMatrix> to be tested
   */
  void unitarityCheck(const multi1d<LatticeColorMatrixF3>& u);
  void unitarityCheck(const multi1d<LatticeColorMatrixD3>& u);

}

#endif
