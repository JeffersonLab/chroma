// -*- C++ -*-
/*! \file
 *  \brief Fermion action factories
 */

#ifndef __qphix_singleton_h__
#define __qphix_singleton_h__

#include "chroma_config.h"

// If we build QPhiX -- Define the singleton holding command line args
#if defined(BUILD_QPHIX)

// Define this for both regular QPhiX and for QDP-JIT QPhiX
#include "qphix/qphix_cli_args.h"
namespace Chroma {
// A Singleton to hold the command line args.
typedef Chroma::SingletonHolder<QPhiX::QPhiXCLIArgs> TheQPhiXParams;
}


#include "qdp_config.h"

#if defined(QDP_IS_QDPJIT) && defined(CHROMA_QPHIX_DSLASH_ENABLED)
#define CHROMA_BUILDING_QPHIX_DSLASH	1
#endif // QDP_IS_QDPJIT && CHROMA_QPHIX_DSLASH_ENABLED

#endif // BUILD_QPHIX
#endif // QPHIX_SINGLETON_H

