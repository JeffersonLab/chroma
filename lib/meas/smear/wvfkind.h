// -*- C++ -*-
// $Id: wvfkind.h,v 2.0 2005-09-25 21:04:40 edwards Exp $
/*! \file
 *  \brief Wave-function types for smearing
 */

#ifndef __wvftype_h__
#define __wvftype_h__

namespace Chroma { 
//! Wave-function types for smearing
enum WvfKind {
  WVF_KIND_GAUSSIAN,
  WVF_KIND_EXPONENTIAL,
  WVF_KIND_GAUGE_INV_GAUSSIAN,
  WVF_KIND_WUPPERTAL,
  WVF_KIND_JACOBI
};

};

#endif
