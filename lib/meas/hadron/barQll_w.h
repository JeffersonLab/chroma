//$Id: barQll_w.h,v 1.8 2009-04-09 22:57:38 caubin Exp $
/*! \file
 *  \brief Heavy Baryon (Qll)  2-pt function : Orginos and Savage
 */

#ifndef __barQll_h__
#define __barQll_h__

#include "chromabase.h"
#include "util/ft/sftmom.h"

namespace Chroma {

//! Lambdaq and SigmaQ 2-pt functions
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * Construct baryon propagators for the LambdaQ and SigmaQ(*) with
 * degenerate "u" and "d" quarks.  In the heavy quark limit the Sigma
 * and Sigma* are degenerate.  The heavy quark is inserted in the 
 * infinitely heavy quark limit
 * by a Wilson-Line without spin indices.
 * We are effectively propagating a spin-0 diquark and a spin-1 diquark.
 *
 * \param u                gauge field (Read) 
 * \param quark_propagator quark propagator ( Read )
 * \param src_coord        cartesian coordinates of the source ( Read )
 * \param phases           holds list of momenta and Fourier phases ( Read )
 * \param xml              xml file object ( Read )
 * \param xml_group        group name for xml data ( Read )
 *
 */

void Qll(const multi1d<LatticeColorMatrix>& u, 
	 const LatticePropagator& quark_propagator,
	 const multi1d<int>& src_coord, 
	 const SftMom& phases,
       	 XMLWriter& xml,
	 const string& xml_group);

void Qll(const multi1d<LatticeColorMatrix>& u, 
	 const LatticePropagator& quark_prop1,
	 const LatticePropagator& quark_prop2,
	 const multi1d<int>& src_coord, 
	 const SftMom& phases,
       	 XMLWriter& xml,
	 const string& xml_group);

//! Heavy Quark Propagator
/*!
 * \ingroup hadron
 *
 * 
 * This constructs the propagator for a spinless Wilson-Line propagating from the 
 * point src_coord forward in time, and vanishing on previous time-slices.
 *
 * \param Qprop              Wilson-Line (write)
 * \param u                  Gauge Field (Read) 
 * \param src_coord          cartesian coordinates of the source ( Read )
 *
 */
void HeavyQuarkProp(LatticeColorMatrix& Qprop,
		    const multi1d<LatticeColorMatrix>& u,
		    const multi1d<int>& src_coord,
		    int length,
		    int bc = 0);

//! Backwards Heavy Quark Propagator
/*!
 * \ingroup hadron
 *
 * 
 * This constructs the propagator for a spinless Wilson-Line propagating from the 
 * point src_coord backward in time, and vanishing on later time-slices.
 *
 * \param Qprop              Wilson-Line (write)
 * \param u                  Gauge Field (Read) 
 * \param src_coord          cartesian coordinates of the source ( Read )
 *
 */
void HeavyQuarkPropBack(LatticeColorMatrix& Qprop,
		    const multi1d<LatticeColorMatrix>& u,
		    const multi1d<int>& src_coord,
		    int length,
		    int bc = 0);



 

//! Function object used for constructing the time-slice set
class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc(int dir): dir_decay(dir) {}

  int operator() (const multi1d<int>& coordinate) const
  {
    if ((dir_decay<0)||(dir_decay>=Nd)) {
      return 0 ;
    } else {
      return coordinate[dir_decay] ;
    }
  }

  int numSubsets() const
  {
    if ((dir_decay<0)||(dir_decay>=Nd)) {
      return 1 ;
    } else {
      return Layout::lattSize()[dir_decay] ;
    }
  }

private:
  TimeSliceFunc() {}  // hide default constructor

  int dir_decay;
};

}

#endif
