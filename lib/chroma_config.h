#ifndef CHROMA_CONFIG_H
#define CHROMA_CONFIG_H

/* Undef the unwanted from the environment -- eg the compiler command line */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#include "chroma_config_internal.h"

/* Prefix everything with CHROMA_QDP_ */
#define CHROMA_PACKAGE PACKAGE
#define CHROMA_PACKAGE_BUGREPORT PACKAGE_BUGREPORT
#define CHROMA_PACKAGE_NAME PACKAGE_NAME
#define CHROMA_PACKAGE_STRING PACKAGE_STRING
#define CHROMA_PACKAGE_TARNAME PACKAGE_TARNAME
#define CHROMA_QDP_PACKAGE_VERSION PACKAGE_VERSION
#define CHROMA_VERSION VERSION
                                                                                
                                                                                
/* Undef the unwanted */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION


#endif
