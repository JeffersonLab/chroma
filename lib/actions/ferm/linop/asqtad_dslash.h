// -*- C++ -*-
//  $Id: asqtad_dslash.h,v 3.1 2007-09-01 18:53:55 edwards Exp $
/*! \file
 *  \brief Include possibly optimized Asqtad dslash
 */

#ifndef DSLASH_S_H
#define DSLASH_S_H

#include "actions/ferm/linop/asq_dsl_s.h"

namespace Chroma 
{
  //! Generic QDP fersion of Asqtad dslash
  /*! \ingroup linop */ 
  typedef QDPStaggeredDslash AsqtadDslash; 

}  // end namespace Chroma

#endif
