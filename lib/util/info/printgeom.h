// -*- C++ -*-
// $Id: printgeom.h,v 1.3 2005-01-14 18:42:38 edwards Exp $
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
