// -*- C++ -*-
// $Id: wvfkind.h,v 1.1 2004-01-13 03:58:32 edwards Exp $
/*! \file
 *  \brief Wave-function types for smearing
 */

#ifndef __wvftype_h__
#define __wvftype_h__

//! Wave-function types for smearing
enum WvfKind {
  WVF_KIND_GAUSSIAN,
  WVF_KIND_EXPONENTIAL,
  WVF_KIND_GAUGE_INV_GAUSSIAN,
  WVF_KIND_WUPPERTAL,
  WVF_KIND_JACOBI
};

#endif
