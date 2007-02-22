//  $Id: loops_w.cc,v 3.1 2007-02-22 21:11:49 bjoo Exp $
//

#include "chromabase.h"

namespace Chroma {


// I cant forward declare this for some reason
// Standard Time Slicery
class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc(int dir): dir_decay(dir) {}
                                                                                
  int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
  int numSubsets() const {return Layout::lattSize()[dir_decay];}
                                                                                
  int dir_decay;
                                                                                
private:
  TimeSliceFunc() {}  // hide default constructor
};




//! Fermion loop code
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * Compute fermion loops via noise source
 *
 * \param q_source -- noise source
 * \param psi      -- M^{-1} on source 
 * \param length   -- length of lattice in time direction  
 * \param xml -- namelist file object ( Read )
 * \param xml_tag -- string used for writing xml data ( Read )
 *
 *        ____
 *        \          dagger
 * m(t) =  >  < Source Gamma  quark_solution >
 *        /
 *        ----
 *          x
 *
 *    where Gamma is a Dirac gamma matrix.
 *
 */

void loops(const LatticeFermion &q_source,
	   const LatticeFermion &psi,
	   int length,
	   XMLWriter& xml_gamma,
	   const string& xml_tag)
{
  push(xml_gamma, xml_tag);

  // Machinery to do timeslice sums with 
  Set timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  multi1d<DComplex> corr_fn_t(length);
  LatticeReal corr_fn_re ;

  // Loop over gamma matrix insertions
  for (int gamma_value=0; gamma_value < (Ns*Ns); ++gamma_value)
  {
    push(xml_gamma,"loop_diagram");     // next array element
    write(xml_gamma,"gamma_value", gamma_value );

    // Construct the meson correlation function
    LatticeComplex corr_fn;
    corr_fn = innerProduct(q_source , Gamma(gamma_value) *psi ) ; 

    corr_fn_t = sumMulti(corr_fn, timeslice);


    multi1d<Real> mesprop(length);
    for (int t=0; t < length; ++t) 
      {
	mesprop[t] = real(corr_fn_t[t]);
      }

      write(xml_gamma, "mesprop", mesprop);

      pop(xml_gamma); // end of array element

  } // end for(gamma_value)

  pop(xml_gamma);

}

}  // end namespace Chroma
