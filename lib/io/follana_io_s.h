/*
 *  Routines associated with simple propagator IO
 */

/*
 *  First the simple propagator header
 */

#ifndef __follana_io_h__
#define __follana_io_h__

namespace Chroma {

/*
 *  Routines for reading and writing propagator
 */

void readQpropFollana(char file[], LatticePropagator& quark_prop, bool swap);

}  // end namespace Chroma

#endif
