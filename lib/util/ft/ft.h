// -*- C++ -*-
// $Id: ft.h,v 3.0 2006-04-03 04:59:12 edwards Exp $
// $Log: ft.h,v $
// Revision 3.0  2006-04-03 04:59:12  edwards
// Major overhaul of fermion and gauge action interface. Basically,
// all fermacts and gaugeacts now carry out  <T,P,Q>  template parameters. These are
// the fermion type, the "P" - conjugate momenta, and "Q" - generalized coordinates
// in the sense of Hamilton's equations. The fermbc's have been rationalized to never
// be over multi1d<T>. The "createState" within the FermionAction is now fixed meaning
// the "u" fields are now from the coordinate type. There are now "ConnectState" that
// derive into FermState<T,P,Q> and GaugeState<P,Q>.
//
// Revision 2.0  2005/09/25 21:04:44  edwards
// Moved to version 2.0
//
// Revision 1.2  2003/04/01 02:38:26  edwards
// Added doxygen comments.
//
// Revision 1.1  2003/03/14 05:06:06  flemingg
// Initial version of SftMom class
//

/*! \file
 * \brief Fourier transform support
 *
 * Utility routines for support of fourier transorm
 */

/*! \defgroup ft Fourier transform support
 * \ingroup util
 *
 * Utility routines for support of fourier transorm
 */

#ifndef __ft_h__
#define __ft_h__

#include "sftmom.h"

#endif
