// $Id: readszinferm_w.h,v 3.0 2006-04-03 04:58:55 edwards Exp $

#ifndef __readszinferm_h__
#define __readszinferm_h__

/*! \file
 *  \brief Read an old SZIN-style (checkerboarded) lattice Dirac fermion
 */

namespace Chroma {

//! Read an old SZIN-style (checkerboarded) lattice Dirac fermion
/*!
 * \param q          fermion ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readSzinFerm(LatticeFermion& q, const string& file);

}  // end namespace Chroma

#endif
