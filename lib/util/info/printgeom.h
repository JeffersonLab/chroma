// -*- C++ -*-
// $Id: printgeom.h,v 1.1 2004-01-29 20:38:09 edwards Exp $
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

void print_geom(XMLWriter& xml);

#endif
