// -*- C++ -*-
// $Id: wvfkind.h,v 1.2 2004-09-22 17:25:00 bjoo Exp $
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

using namespace Chroma;
#endif
