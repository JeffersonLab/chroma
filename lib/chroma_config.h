#ifndef CHROMA_CONFIG_H
#define CHROMA_CONFIG_H

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION


#include "chroma_config_internal.h"

#ifdef PACKAGE_STRING
static const char* const CHROMA_PACKAGE_STRING(PACKAGE_STRING);
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_VERSION
static const char* const CHROMA_PACKAGE_VERSION(PACKAGE_VERSION);
#undef PACKAGE_VERSION
#endif


/* Undef the unwanted from the environment -- eg the compiler command line */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#endif
