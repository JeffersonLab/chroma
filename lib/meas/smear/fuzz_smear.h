// -*- C++ -*-
/*! \file
 *  \brief Fuzzed sources
 */

#ifndef __fuzz_smear__
#define __fuzz_smear__

namespace Chroma 
{

  //! Fuzzed source 
  /*! \ingroup smear */
  void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		  const LatticeColorVector & psi, 
		  LatticeColorVector & psifuzz, 
		  int length, int j_decay);


  //! Fuzzed source 
  /*! \ingroup smear */
  void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		  const LatticePropagator  & psi, 
		  LatticePropagator& psifuzz, 
		  int length, int j_decay);

  //! Fuzzed source 
  /*! \ingroup smear */
  void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		  const LatticeFermion  & psi, 
		  LatticeFermion& psifuzz, 
		  int length, int j_decay);

  //! Fuzzed source 
  /*! \ingroup smear */
  void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		  const LatticeStaggeredFermion  & psi, 
		  LatticeStaggeredFermion& psifuzz, 
		  int length, int j_decay); 

}  // end namespace Chroma

#endif
