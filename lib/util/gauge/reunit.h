// -*- C++ -*-

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
  namespace ReunitEnv {
    extern double getTime();
  }

  enum Reunitarize {REUNITARIZE, REUNITARIZE_ERROR, REUNITARIZE_LABEL};
  
  void reunit(LatticeColorMatrixF& xa);
  void reunit(LatticeColorMatrixD& xa);
  
  void reunit(LatticeColorMatrixF& xa,
	      const Subset& mstag);

  void reunit(LatticeColorMatrixD& xa,
	      const Subset& mstag);
  
  // With ruflag
  void reunit(LatticeColorMatrixF& xa,
	      int& numbad, 
	      enum Reunitarize ruflag);

  void reunit(LatticeColorMatrixD& xa,
	      int& numbad, 
	      enum Reunitarize ruflag);
  
  
  void reunit(LatticeColorMatrixF& xa,
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const Subset& mstag);
  
  void reunit(LatticeColorMatrixD& xa,
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const Subset& mstag);
  
  // With numbad and ruflag
  void reunit(LatticeColorMatrixF& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag);
  
  void reunit(LatticeColorMatrixD& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag);
  
  void reunit(LatticeColorMatrixF& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const Subset& mstag);

  void reunit(LatticeColorMatrixD& xa, 
	      LatticeBoolean& bad, 
	      int& numbad, 
	      enum Reunitarize ruflag,
	      const Subset& mstag);
  
} // End namespace
#endif
