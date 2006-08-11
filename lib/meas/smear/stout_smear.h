// -*- C++ -*-
//  $Id: stout_smear.h,v 3.2 2006-08-11 18:12:00 edwards Exp $
/*! \file
 *  \brief Stout smear a gauge field
 */


#ifndef __stout_smear_h__
#define __stout_smear_h__

namespace Chroma 
{ 
  //! Construct the "stout-smeared" links
  /*! 
   * \ingroup smear
   *
   * Arguments:
   *
   * \param u_smear      smeared gauge field ( Write )
   * \param u            gauge field ( Read )
   * \param mu           direction of smeared gauge field ( Read )
   * \param sm_fact      smearing factor ( Read )
   * \param smear_dirs	 no staples in these directions (Read)
   */

  void stout_smear(LatticeColorMatrix& u_smear,
		   const multi1d<LatticeColorMatrix>& u,
		   int mu, 
		   const Real& sm_fact, 
		   multi1d<bool> smear_dirs);


}

#endif
