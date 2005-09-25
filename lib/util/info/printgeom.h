// -*- C++ -*-
// $Id: printgeom.h,v 2.0 2005-09-25 21:04:45 edwards Exp $
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
