// $Id: chromainc.h,v 1.12 2004-01-02 01:42:59 edwards Exp $
//
// Include file that includes all the include files.

#ifndef CHROMAINC_INCLUDE
#define CHROMAINC_INCLUDE

#include "qdp.h"
#include "qdp_iogauge.h"

#include "linearop.h"
//#include "primitives.h"
//#include "common_declarations.h"

#include "lib.h"

using namespace QDP;

enum Sources {POINT_SOURCE, WALL_SOURCE, SHELL_SOURCE, BNDST_SOURCE, POINT_AND_BNDST_SOURCE, 
	      SHELL_AND_BNDST_SOURCE, POINT_AND_SHELL_AND_BNDST_SOURCE};
enum Sinks {POINT_SINK, WALL_SINK, POINT_AND_WALL_SINK, SHELL_SINK, POINT_AND_SHELL_SINK, 
	    BNDST_SINK, POINT_AND_BNDST_SINK, SHELL_AND_BNDST_SINK, POINT_AND_SHELL_AND_BNDST_SINK};

#endif
