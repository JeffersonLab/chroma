#ifndef CPP_DSLASH_CONFIG_H
#define CPP_DSLASH_CONFIG_H

/* Undef the unwanted from the environment -- eg the compiler command line */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#include "dslash_config_internal.h"

/* Prefix everything with SSE */
#define CPP_DSLASH_PACKAGE PACKAGE
#define CPP_DSLASH_PACKAGE_BUGREPORT PACKAGE_BUGREPORT
#define CPP_DSLASH_PACKAGE_NAME PACKAGE_NAME
#define CPP_DSLASH_PACKAGE_STRING PACKAGE_STRING
#define CPP_DSLASH_PACKAGE_TARNAME PACKAGE_TARNAME
#define CPP_DSLASH_PACKAGE_VERSION PACKAGE_VERSION
#define CPP_DSLASH_VERSION VERSION
                                                                                
                                                                                
/* Undef the unwanted */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION


#endif
