// -*- C++ -*-
// $Id: gaugeact.h,v 1.2 2003-10-10 03:46:46 edwards Exp $

/*! @file
 * @brief Class structure for gauge actions
 */

#ifndef __gaugeactt_h__
#define __gaugeactt_h__

using namespace QDP;

#include "linearop.h"

//! Base class for gauge actions
/*! @ingroup actions
 *
 * Supports creation and application for gauge actions
 */

class GaugeAction
{
public:
  //! Virtual destructor to help with cleanup;
  virtual ~GaugeAction() {}
};


#endif
