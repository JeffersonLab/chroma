// -*- C++ -*-
// $Id: wilslp.h,v 1.1 2004-04-26 16:12:49 mcneile Exp $
/*! \file
 *  \brief Calculate Wilson loops
 */

#ifndef __wilslp_h__
#define __wilslp_h__

//! 
/*!
 * \ingroup glue
 *
*/

/* WILSLP -- calculates, depending on option, (1) "space-like" planar    */
/*           Wilson loops in the directions perpendicular to j_decay     */
/*           that have equal length, (2) "time-like" planar Wilson loops */
/*           with time direction j_decay and space directions the        */
/*           perpendicular ones that have equal length and (3) off-axis  */
/*           "time-like" Wilson loops along 3 paricular paths in the     */
/*           space directions that have equal length.                    */
/* Translated by ACI from the SZIN code originally by Urs Heller.        */
/*                                                                       */
/* u -- gauge field (Read)                                               */
/* j_decay -- time direction (Read)                                      */
/* kind -- binary-combined YES/NO [1/0] of the three options (Read)      */

void wilslp(const multi1d<LatticeColorMatrix>& u,
        int j_decay, int kind,
        XMLWriter& xml, const string& xml_group);

#endif
