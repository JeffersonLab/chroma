// -*- C++ -*-
// $Id: invtype.h,v 1.1 2003-11-13 04:09:33 edwards Exp $

/*! @file
 * @brief Types of inverters
 */

#ifndef __invtype_h__
#define __invtype_h__

enum InvType {
  CG_INVERTER = 21, 
  MR_INVERTER = 22,
  BICG_INVERTER = 23,
  CR_INVERTER = 24};

#endif
