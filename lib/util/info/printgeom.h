// -*- C++ -*-
// $Id: printgeom.h,v 3.0 2006-04-03 04:59:13 edwards Exp $
/*! \file
 *  \brief Print out machine geometry and problem size info
 */

#ifndef __printgeom_h__
#define __printgeom_h__

namespace Chroma {

//! Print out machine geometry and problem size info
/*!
 * \ingroup info
 *
 * Arguments:
 *
 *  \param xml          The xml stream to write the info
 */

void printgeom(XMLWriter& xml);

}  // end namespace Chroma

#endif
