// -*- C++ -*-
/*! \file
 *  \brief Print out basic info about this program
 */

#ifndef __proginfo_h__
#define __proginfo_h__

namespace Chroma {

//! Print out basic info about this program
/*!
 * \ingroup info
 *
 * Arguments:
 *
 *  \param xml          The xml stream to write the info
 */

void proginfo(XMLWriter& xml);

}  // end namespace Chroma

#endif
