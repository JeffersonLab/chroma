// $Id: barcomp_io.h,v 1.3 2005-01-14 20:13:06 edwards Exp $
/*! \file
 * \brief Routines associated with barcomp generalized correlator IO
 */

#ifndef __barcomp_io_h__
#define __barcomp_io_h__

namespace Chroma {

/*
 *  Routines for reading and writing the QQQ correlators
 */
void convertBarcomp(multi1d<Complex>& barprop1d, const multiNd<Complex>& barprop, 
		    const int j_decay);

}  // end namespace Chroma

#endif
