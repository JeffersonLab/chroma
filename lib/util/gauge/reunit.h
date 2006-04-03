// -*- C++ -*-
// $Id: reunit.h,v 3.0 2006-04-03 04:59:12 edwards Exp $

/*! \file
 *  \ingroup gauge
 *  \author Subsetting added by A. Hart
 *  \param[in,out] xa  The descriptor of matrices to be reunitarized.
 *            Must be of type LatticeColorMatrix
 *  \param[in] bad Descriptor of flags indicating sites violating unitarity.
 *            Only used if ruflag = REUNITARIZE_LABEL or
 *            REUNITARIZE_ERROR.
 *  \param[in] ruflag Can also be REUNITARIZE in which case the
 *            matrices are reunitarized but no complaints are made.
 *  \param[out] numbad Total number of matrices violating unitarity.
 *            ONLY USED IF ruflag is testing for ERROR or LABEL. 
 *  \param[in] mstag  An (un)ordered subset of sites
 *  \brief Reunitarize in place a color matrix to SU(N)
 *
 *  Reunitarize (to a SU(N)) inplace the matrix XA under some option
 */

#ifndef __reunit_h__
#define __reunit_h__

namespace Chroma {

  enum Reunitarize {REUNITARIZE, REUNITARIZE_ERROR, REUNITARIZE_LABEL};
  
  void reunit(LatticeColorMatrix& xa);
  
  void reunit(LatticeColorMatrix& xa,
	      const UnorderedSubset& mstag);
  
  void reunit(LatticeColorMatrix& xa,
	      const OrderedSubset& mstag);
  
  // With ruflag
  void reunit(LatticeColorMatrix& xa,
	      int& numbad, 
	      enum Reunitarize ruflag);
  
  void reunit(LatticeColorMatrix& xa,
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const UnorderedSubset& mstag);
  
  void reunit(LatticeColorMatrix& xa,
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const OrderedSubset& mstag);
  
  // With numbad and ruflag
  void reunit(LatticeColorMatrix& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag);
  
  void reunit(LatticeColorMatrix& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const UnorderedSubset& mstag);
  
  void reunit(LatticeColorMatrix& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const OrderedSubset& mstag);
  
}; // End namespace
#endif
