// $Id: readszinferm_w.h,v 2.0 2005-09-25 21:04:32 edwards Exp $

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
