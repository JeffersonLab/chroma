// -*- C++ -*-
//  $Id: stout_smear_old.h,v 3.1 2008-01-20 03:07:51 edwards Exp $
/*! \file
 *  \brief Stout smear a gauge field using the old (non-gauge covariant) method
 */


#ifndef __stout_smear_old_h__
#define __stout_smear_old_h__

namespace Chroma 
{ 
  //! Construct the "stout-smeared" links using the non-gauge covariant method
  /*! 
   * \ingroup smear
   *
   *  Stout smear a gauge field using the old (non-gauge covariant) method
   *
   * Arguments:
   *
   * \param u_smear      smeared gauge field ( Write )
   * \param u            gauge field ( Read )
   * \param mu           direction of smeared gauge field ( Read )
   * \param sm_fact      smearing factor ( Read )
   * \param smear_dirs	 no staples in these directions (Read)
   */

  void stout_smear_old(LatticeColorMatrix& u_smear,
		       const multi1d<LatticeColorMatrix>& u,
		       int mu, 
		       const Real& sm_fact, 
		       multi1d<bool> smear_dirs);


}

#endif
