// $Id: chromainc.h,v 1.6 2003-03-30 17:14:36 edwards Exp $
//
// Include file that includes all the include files.

#include "qdp.h"

#include "linearop.h"
#include "primitives.h"
#include "common_declarations.h"

#include "actions/actions.h"
#include "util/util.h"
#include "meas/meas.h"
#include "update/update.h"
#include "io/io.h"


using namespace QDP;

enum Exponentiate {TWELWTH_ORDER, EXACT};
enum Sources {POINT_SOURCE, WALL_SOURCE, SHELL_SOURCE, BNDST_SOURCE, POINT_AND_BNDST_SOURCE, 
	      SHELL_AND_BNDST_SOURCE, POINT_AND_SHELL_AND_BNDST_SOURCE};
enum Sinks {POINT_SINK, WALL_SINK, POINT_AND_WALL_SINK, SHELL_SINK, POINT_AND_SHELL_SINK, 
	    BNDST_SINK, POINT_AND_BNDST_SINK, SHELL_AND_BNDST_SINK, POINT_AND_SHELL_AND_BNDST_SINK};

