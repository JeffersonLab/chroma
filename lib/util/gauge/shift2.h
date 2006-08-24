// -*- C++ -*-
// $Id: shift2.h,v 1.2 2006-08-24 03:14:39 edwards Exp $
/*! \file
 *  \brief Shift by a power of 2
 */

#ifndef __shift2_h__
#define __shift2_h__

namespace Chroma
{

  //! A simple not-fancy power of 2 shift
  /*! \ingroup gauge */
  LatticeColorMatrix shift2(const LatticeColorMatrix& s1, int isign, int dir, int level);

}

#endif
