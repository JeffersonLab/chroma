/*
 *  $Id: follana_io_s.cc,v 3.1 2007-06-10 14:40:23 edwards Exp $
 *
 *  These are a few simple I/O routines that we can use until QIO makes its appearance
 *  I have tried to include a simple header by means of a structure.
 */

#include "chromabase.h"
#include "io/follana_io_s.h"
#include "qdp_util.h"
#include <string>

namespace Chroma {


void readQpropFollana(char file[], LatticePropagator& quark_prop, bool swap)
{
  /*
   *  Now the local variables
   */

  BinaryFileReader prop_in(file);
 
  int x, y, z, t, src_col, snk_col, src_spin, snk_spin, index;
  Propagator site_prop;
 
  Complex tmp_cmpx;
  const multi1d<int>& latt_size = Layout::lattSize();
  multi1d<int> coords(Nd);

  multi1d<Real64> buf( latt_size[0]*Nc*Nc*2 );
  /* Get size of the lattice */
  /* Read the file somehow */

  for( t = 0; t < latt_size[3]; t++) {  
    QDPIO::cout << "Reading timeslice: " << t << endl;

    for( z = 0; z < latt_size[2]; z++) { 
      for( y = 0; y < latt_size[1]; y++ ) {


        
	read(prop_in, buf, latt_size[0]*Nc*Nc*2);

	if( swap == true ) { 
	  QDPUtil::byte_swap((void *)&buf[0], sizeof(Real64), latt_size[0]*Nc*Nc*2 );
        } 
	index = 0;

	for( x = 0; x < latt_size[0]; x++) {

	  coords[0] = x;
	  coords[1] = y;
	  coords[2] = z;
	  coords[3] = t;
	 
	  for( snk_spin = 0; snk_spin < Ns; snk_spin++) { 
	    for( src_spin = 0; src_spin < Ns; src_spin++) { 

	      ColorMatrix tmp_col;

	      for( snk_col = 0; snk_col < Nc; snk_col++) { 
		for( src_col = 0; src_col < Nc; src_col++) { 

		  Real32 re, im;
		  
		  re = (Real32)(buf[index]);
		  index++;
		  im = (Real32)(buf[index]);
		  index++;

		  tmp_cmpx = cmplx(Real(re), Real(im));
		  
		  tmp_col = pokeColor(tmp_col,
			    tmp_cmpx,
			    snk_col,
			    src_col);
		}
	      }
		
	      site_prop = pokeSpin(site_prop,
		       tmp_col,
		       snk_spin,
		       src_spin);
	    }
	  }

	  pokeSite(quark_prop,
		   site_prop,
		   coords);

	}
      }
    }
  }
}


}  // end namespace Chroma
