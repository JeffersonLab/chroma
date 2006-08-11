// $Id: fuzz_smear.cc,v 3.1 2006-08-11 16:13:29 edwards Exp $
/*! \file
 *  \brief Fuzzed sources
 */

#include "chromabase.h"
#include "meas/smear/displacement.h"
#include "meas/smear/fuzz_smear.h"

namespace Chroma 
{

  // mcneile: this routine is still being checked **
  // mcneile: I have now checked that it is gauge invariant

  //! apply a fuzz_smear operator to a lattice field
  /*!
   * \ingroup smear
   *
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
   *  For smearing at the source the input is a local source.
   *
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
    if( Ns == 1 &&  length % 2 == 1)
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
	displacement(ufuzz,psi_dir,length,mu);
	if( is_initial )
	{
	  psifuzz = psi_dir ;
	  is_initial = false ;
	}
	else
	{
	  psifuzz +=  psi_dir ;
	}


	// negative direction
	psi_dir = psi ; 
	int neg_length = -length ; 
	displacement(ufuzz,psi_dir,neg_length,mu);
	psifuzz +=  psi_dir ;

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


  void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		  const LatticeStaggeredFermion  & psi, 
		  LatticeStaggeredFermion& psifuzz, 
		  int length, int j_decay)
  {
    fuzz_smear<LatticeStaggeredFermion>(ufuzz, psi, psifuzz, 
					length,j_decay) ;
  }

}  // end namespace Chroma

