// -*- C++ -*-
// $Id: stout_utils.h,v 1.3 2008-01-25 22:22:39 edwards Exp $
/*! \file
 *  \brief Stout utilities
 */

#ifndef STOUT_UTILS_H
#define STOUT_UTILS_H

#include "chromabase.h"

#ifdef QDP_IS_QDPJIT
CUfunction function_get_fs_bs_build(const LatticeColorMatrix& Q,
                                    const LatticeColorMatrix& QQ,
                                    multi1d<LatticeComplex>& f,
                                    multi1d<LatticeComplex>& b1,
                                    multi1d<LatticeComplex>& b2,
                                    bool dobs);
CUfunction function_get_fs_bs_exec(CUfunction function,
                                   const LatticeColorMatrix& Q,
                                   const LatticeColorMatrix& QQ,
                                   multi1d<LatticeComplex>& f,
                                   multi1d<LatticeComplex>& b1,
                                   multi1d<LatticeComplex>& b2,
                                   bool dobs);
#endif


namespace Chroma 
{

  /*!
   * Stouting
   *
   * \ingroup gauge
   *
   * @{
   */

  /*! \ingroup gauge */
  namespace StoutLinkTimings { 
    double getForceTime();
    double getSmearingTime();
    double getFunctionsTime();
  }

  /*! \ingroup gauge */
  namespace Stouting 
  {
    //! Given field U, form Q and Q^2
    void getQs(const multi1d<LatticeColorMatrix>& u, 
	       LatticeColorMatrix& Q, 
	       LatticeColorMatrix& QQ,
	       int mu,
	       const multi1d<bool>& smear_in_this_dirP,
	       const multi2d<Real>& rho);

    //! Given field U, construct the staples into C, form Q and Q^2 and compute  c0 and c1
    void getQsandCs(const multi1d<LatticeColorMatrix>& u, 
		    LatticeColorMatrix& Q, 
		    LatticeColorMatrix& QQ,
		    LatticeColorMatrix& C, 
		    int mu,
		    const multi1d<bool>& smear_in_this_dirP,
		    const multi2d<Real>& rho);

    //! Given c0 and c1 compute the f-s and b-s
    /*! \ingroup gauge */
    void getFs(const LatticeColorMatrix& Q,
	       const LatticeColorMatrix& QQ,
	       multi1d<LatticeComplex>& f);

    //! Given c0 and c1 compute the f-s and b-s
    /*! Only compute b-s if do_bs is set to true (default) */
    void getFsAndBs(const LatticeColorMatrix& Q,
		    const LatticeColorMatrix& QQ,
		    multi1d<LatticeComplex>& f,
		    multi1d<LatticeComplex>& b1,
		    multi1d<LatticeComplex>& b2,
		    bool do_bs=true);
    
    //! Stout smear in a specific link direction
    void stout_smear(LatticeColorMatrix& next,
		     const multi1d<LatticeColorMatrix>& current, 
		     int mu,
		     const multi1d<bool>& smear_in_this_dirP,
		     const multi2d<Real>& rho);

    //! Do the smearing from level i to level i+1
    void smear_links(const multi1d<LatticeColorMatrix>& current,
		     multi1d<LatticeColorMatrix>& next, 
		     const multi1d<bool>& smear_in_this_dirP,
		     const multi2d<Real>& rho);
    
    //! Do the force recursion from level i+1, to level i
    void deriv_recurse(multi1d<LatticeColorMatrix>&  F,
		       const multi1d<bool>& smear_in_this_dirP,
		       const multi2d<Real>& rho,
		       const multi1d<LatticeColorMatrix>& u);

  }

  /*! @} */   // end of group gauge
}

#endif
