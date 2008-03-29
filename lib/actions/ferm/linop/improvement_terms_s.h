// -*- C++ -*-
//  $Id: improvement_terms_s.h,v 3.3 2008-03-29 19:15:36 mcneile Exp $
/*! \file
 *  \brief Support for Asqtad
 */

#ifndef IMPROVEMENT_TERMS_S
#define IMPROVEMENT_TERMS_S

#include "chromabase.h"

namespace Chroma 
{ 
  //! FAT7_LINKS 
  /*!
   * \ingroup linop
   * Construct the "fat" links with staples up to 7 links long
   * used in the staggered "asqtad" action
   *
   * NOTE: the staggered phase factors are assumed to be included
   *       in the gauge fields u
   * Arguments:
   *
   *  \param u       gauge field (Read)
   *  \param u_fat   "fat-link" gauge field (Write) 
   *  \param u0      tapdole factor
   */
  void Fat7_Links(multi1d<LatticeColorMatrix>& u, multi1d<LatticeColorMatrix>& u_fat, Real u0);

  /*! \ingroup linop */
  void Triple_Links(multi1d<LatticeColorMatrix>& u, multi1d<LatticeColorMatrix>& u_triple, Real u0);

  void Triple_Links(multi1d<LatticeColorMatrix> & u,
		    multi1d<LatticeColorMatrix> & ut,
		    Real u0, Real c_3) ;


  //! Pass parameters to the fat link code
  /*! \ingroup linop */
  class fat7_param
  {
  public :
    Real c_1l ;
    Real c_3l ;
    Real c_5l ;
    Real c_7l ;
    Real c_Lepage ;
  };

  /*! \ingroup linop */
  void Fat7_Links(multi1d<LatticeColorMatrix> & u,
		  multi1d<LatticeColorMatrix> & uf,
		  fat7_param & pp);

} // End Namespace Chroma


#endif

