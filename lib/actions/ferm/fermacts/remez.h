// -*- C++ -*-
// $Id: remez.h,v 2.0 2005-09-25 21:04:26 edwards Exp $
/*! \file
 *  \brief Redirector for Remez algorithm for finding nth roots
 */

#ifndef __remez_h__
#define __remez_h__

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



