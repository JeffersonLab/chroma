// -*- C++ -*-
//  $Id: stout_smear.h,v 3.1 2006-08-11 16:13:30 edwards Exp $
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
   *  \param u		gauge field (Read)
   *  \param u_smr	"stout-smeared" gauge field (Write)
   */

  void stout_smear(LatticeColorMatrix& u_smear,
		   const multi1d<LatticeColorMatrix>& u,
		   int mu, 
		   const Real& sm_fact, int j_decay);


}

#endif
