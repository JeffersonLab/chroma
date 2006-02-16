// -*- C++ -*-
// $Id: remez.h,v 2.1 2006-02-16 21:02:58 edwards Exp $
/*! \file
 *  \brief Redirector for Remez algorithm for finding nth roots
 */

#ifndef __remez_h__
#define __remez_h__

#include "chroma_config.h"

#ifdef BUILD_GMP_REMEZ        // If GMP is defined 

#include "remez_gmp.h"

// The following is an ifdef lis that switches in different
// versions of the Remez implementation

// GMP (Gnu Multi-Precision) version
namespace Chroma
{
  typedef RemezGMP   Remez;
}

#else

#include "remez_stub.h"

// Dummy version
namespace Chroma
{
  typedef RemezStub   Remez;
}

#endif  // ifdef GMP

#endif  // include guard



