// $Id: barcomp_io.h,v 1.2 2004-04-21 17:25:22 edwards Exp $
/*! \file
 * \brief Routines associated with barcomp generalized correlator IO
 */

#ifndef __barcomp_io_h__
#define __barcomp_io_h__

/*
 *  Routines for reading and writing the QQQ correlators
 */
void convertBarcomp(multi1d<Complex>& barprop1d, const multiNd<Complex>& barprop, 
		    const int j_decay);

#endif
