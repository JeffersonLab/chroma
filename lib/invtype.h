// -*- C++ -*-
// $Id: invtype.h,v 1.1 2004-01-08 03:11:47 edwards Exp $

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
