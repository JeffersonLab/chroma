// $Id: fuzz_smear.cc,v 1.2 2004-02-08 14:35:07 mcneile Exp $

#include "chromabase.h"
#include "meas/smear/fuzz_smear.h"

using namespace QDP;

// mcneile: this routine is still being checked **

//! apply a fuzz_smear operator to a lattice field
/*!
 * Arguments:
 *
 *  \param ufuzz    fuzzed gauge field ( Read )
 *  \param psi      color vector field ( Read )
 *  \param psifuzz     color vector field ( Modify )
 *  \param length   length of fuzz_smear ( Read )
 *  \param j_decay      direction of time ( Read )
 *
 *
 * Discription:
 *
 *
 *  The sink fuzzed propagator is calculated by extending the path of the input
 *  propagator and is sum of the two contributions as follows:  
 *  In the positive mu direction
 *              ___3
 *     fuzz     \    +        +            +
 *  psi   (x) =  )  U (x-mu) U (x-2mu)... U (x-length mu) psi(x-length mu)
 *              /    mu       mu           mu
 *              ---    
 *             mu = 1
 *  In the negative mu direction
 *              ___3
 *     fuzz     \
 *  psi   (x) =  )  U (x) U (x+mu) U (x+2mu)... U (x+(length-1)mu) psi(x+length mu)
 *              /    mu    mu       mu           mu
 *              ---
 *             mu = 1
 *   
 *
 *
 *                                                  fuzz
 *  Where U  defined in the the above equations is U
 *         mu                                       mu
 *  length is an input parameter that describes the size of the
 *  smearing.
 * 
 *  This code assumes that the gauge field has been fuzzed
 *  externally.
 * 
 *  (Documentation from Peter Boyle's Fortran code). 
 * 
 * 
 * Reference:
 * 
 *  EFFICIENT HADRONIC OPERATORS IN LATTICE GAUGE THEORY.
 * By UKQCD Collaboration (P. Lacock et al.). 
 * Published in Phys.Rev.D51:6403-6410,1995, hep-lat/9412079 
 *
 */


template<typename T>
void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		  const T& psi, T& psifuzz, 
		  int length, int j_decay)
{
  // Initial ferm field
  T psi_dir ;


  // if staggered then direction must be even
  // or else an error
  if( Ns == 0 &&  length % 2 == 1)
    {
      cout << "fuzz_smear::Error fuzzing length =  " << length << endl ; 
      cout << "Fuzzing length must be even for staggered fermions" << endl ; 
      QDP_abort(1);
    }



  //
  // loop over directions x,y,z
  //
  bool is_initial = true ;

  for(int mu = 0 ; mu < Nd ; ++mu)
    { 
      if( mu != j_decay )
	{

	  // positive direction
	  psi_dir = psi ; 
	  for(int n = 0; n < length; ++n)
	    {
	      T tmp = shift(psi_dir, FORWARD, mu);
	      psi_dir = ufuzz[mu] * tmp;
	    }

	  if( is_initial )
	    {
	      psifuzz = psi_dir ;
	      is_initial = false ;
	    }
	  else
	    {
	      psifuzz = psifuzz + psi_dir ;
	    }


	  // negative direction
	  psi_dir = psi ; 
	  for(int n = 0; n < length; ++n)
	    {
	      T tmp = shift(psi_dir, BACKWARD, mu);
	      psi_dir = adj(ufuzz[mu]) * tmp;
	    }
	      psifuzz = psifuzz + psi_dir ;


	} // end of x,y,z direction

    } // loop over directions

}




void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		const LatticeColorVector & psi, 
		LatticeColorVector & psifuzz, 
		int length, int j_decay)
{
  fuzz_smear<LatticeColorVector>(ufuzz, psi, psifuzz, 
				 length,j_decay) ;
}



void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		const LatticePropagator  & psi, 
		LatticePropagator& psifuzz, 
		int length, int j_decay)
{
  fuzz_smear<LatticePropagator>(ufuzz, psi, psifuzz, 
				length,j_decay) ;
}



void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		const LatticeFermion  & psi, 
		LatticeFermion& psifuzz, 
		int length, int j_decay)
{
  fuzz_smear<LatticeFermion>(ufuzz, psi, psifuzz, 
				length,j_decay) ;
}


