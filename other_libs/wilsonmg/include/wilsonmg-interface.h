// =================================================================
// QDP Clover Multigrid Chroma Interface
//   Interface
//     This header provides the declarations that Chroma needs in order
//     to invoke functions in QDP without exposing parts of the QDP
//     software stack that would cause problems or conflicts.

#ifndef __wilmg_h__
#define __wilmg_h__

//#include <qcdlib.h>
#include <math.h>
#include <qla.h>

//extern QDP_Layout *QDP_layout_hyper_eo2;

// Define both single- and double-precision multigrids
#include "generic_undef.h"
#include "generic_f.h"
#include "wilsonmg_p.h"
#include "generic_undef.h"
#include "generic_d.h"
#include "wilsonmg_p.h"
#include "generic_undef.h"

// Restore whatever precision we're actually using
#if QDP_Precision == 'F' || QDP_Precision == 1
#include "generic_f.h"
#else
#include "generic_d.h"
#endif

#endif
