// -*- C++ -*-
// $Id: printgeom.h,v 1.2 2004-01-29 20:46:29 edwards Exp $
/*! \file
 *  \brief Print out machine geometry and problem size info
 */

#ifndef __printgeom_h__
#define __printgeom_h__

//! Print out machine geometry and problem size info
/*!
 * \ingroup info
 *
 * Arguments:
 *
 *  \param xml          The xml stream to write the info
 */

void printgeom(XMLWriter& xml);

#endif
