// $Id: readszinferm_w.h,v 1.1 2003-05-09 20:27:11 edwards Exp $

#ifndef __readszinferm_h__
#define __readszinferm_h__

/*! \file
 *  \brief Read an old SZIN-style (checkerboarded) lattice Dirac fermion
 */

//! Read an old SZIN-style (checkerboarded) lattice Dirac fermion
/*!
 * \param q          fermion ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readSzinFerm(LatticeFermion& q, const string& file);

#endif
