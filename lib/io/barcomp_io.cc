// $Id: barcomp_io.cc,v 1.2 2004-04-21 17:25:21 edwards Exp $
/*! \file
 * \brief Routines associated with QQQ generalized correlator IO
 */

#include "chromabase.h"
#include "io/barcomp_io.h"

//! Convert generalized correlator object
void convertBarcomp(multi1d<Complex>& barprop_1d, const multiNd<Complex>& barprop, 
		    const int j_decay)
{
  /*
   * Local variables
   */
  
  int length = Layout::lattSize()[j_decay]; // Temporal extent of lattice

  multi1d<int> ranks(7);

  barprop_1d.resize(length*Ns*Ns*Ns*Ns*Ns*Ns);

  int cnt = 0;
  for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])           // sf_3
    for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])         // sf_2
      for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])       // sf_1
	for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])     // si_3
	  for(ranks[4]=0; ranks[4] < Ns; ++ranks[4])   // si_2
	    for(ranks[5]=0; ranks[5] < Ns; ++ranks[5]) // si_1
	      for(ranks[6] = 0; ranks[6] < length; ++ranks[6])
		{
		  barprop_1d[cnt++] = barprop[ranks];
		}
}
