/*
 *  This include file supports various types of quark sources, including
 *  P-wave, D-wave etc
 */

#ifndef __sources_h__
#define __sources_h__

#include "srcsnktype.h"
#include "wavetype.h"

#ifdef CHROMA_BUILD_WILSON
#include "sources_w.h"
#elif CHROMA_BUILD_STAGGERED

#endif

#endif
