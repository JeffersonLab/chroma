// -*- C++ -*-
// $Id: lib.h,v 3.1 2006-10-19 16:01:26 edwards Exp $

/*! \file
 * \brief Chroma Lattice Field Theory library
 *
 * Central include file for all pieces of the Chroma lattice field theory library
 */

/*! \defgroup lib Chroma Lattice Field Theory library
 *
 * Central include file for all pieces of the Chroma lattice field theory library
 */


#ifndef __lib_h__
#define __lib_h__

#include "init/chroma_init.h"

#include "fermact.h"

#include "wilstype_fermact_w.h"
#include "unprec_wilstype_fermact_w.h"
#include "eoprec_wilstype_fermact_w.h"
#include "eoprec_constdet_wilstype_fermact_w.h"
#include "eoprec_logdet_wilstype_fermact_w.h"
#include "stagtype_fermact_s.h"

#include "fermbc.h"
#include "gaugeact.h"
#include "gaugebc.h"
#include "handle.h"
#include "state.h"

#include "linearop.h"
#include "eoprec_linop.h"
#include "eoprec_constdet_linop.h"
#include "eoprec_logdet_linop.h"

#include "tprec_linop.h"
#include "tprec_logdet_linop.h"
#include "teoprec_linop.h"
#include "teoprec_logdet_linop.h"

#include "actions/actions.h"
#include "update/update.h"
#include "util/util.h"
#include "meas/meas.h"

#include "update/update.h"
#include "io/io.h"

#endif


