#ifndef WILSONMG_CONFIG_H
#define WILSONMG_CONFIG_H


/* Undef the unwanted from the environment -- eg the compiler command line */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#include "wilsonmg_config_internal.h"

/* Prefix everything with WILSONMG_QDP_ */
static const char* const WILSONMG_PACKAGE(PACKAGE);
static const char* const WILSONMG_PACKAGE_BUGREPORT(PACKAGE_BUGREPORT);
static const char* const WILSONMG_PACKAGE_NAME(PACKAGE_NAME);
static const char* const WILSONMG_PACKAGE_STRING(PACKAGE_STRING);
static const char* const WILSONMG_PACKAGE_TARNAME(PACKAGE_TARNAME);
static const char* const WILSONMG_PACKAGE_VERSION(PACKAGE_VERSION);
static const char* const WILSONMG_VERSION(VERSION);

/* Undef the unwanted */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#endif
