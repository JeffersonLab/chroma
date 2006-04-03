// -*- C++ -*-
// $Id: wilslp.h,v 3.0 2006-04-03 04:58:58 edwards Exp $
/*! \file
 *  \brief Calculate Wilson loops
 */

#ifndef __wilslp_h__
#define __wilslp_h__

namespace Chroma {

//! Calculate Wilson loops
/*!
 * \ingroup glue
 *
 * Calculates, depending on option, (1) "space-like" planar   
 * Wilson loops in the directions perpendicular to j_decay    
 * that have equal length, (2) "time-like" planar Wilson loops
 * with time direction j_decay and space directions the       
 * perpendicular ones that have equal length and (3) off-axis 
 * "time-like" Wilson loops along 3 paricular paths in the    
 * space directions that have equal length.                   
 *
 * \param u          gauge field (Read)                                              
 * \param j_decay    decay direction (Read)                                     
 * \param t_dir      time direction (Read)                                     
 * \param kind       binary-combined YES/NO [1/0] of the three options (Read)      
 *                   e.g. kind = 2 gives planar t-like, kind=6 is 
 *                   planar + off-axis: sqrt(2), sqrt(5), sqrt(3)
 */

void wilslp(const multi1d<LatticeColorMatrix>& u,
	    int j_decay, int t_dir, int kind,
	    XMLWriter& xml, const string& xml_group);

}  // end namespace Chroma

#endif
