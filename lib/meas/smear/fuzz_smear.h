// $Id: fuzz_smear.h,v 2.0 2005-09-25 21:04:39 edwards Exp $

#ifndef __fuzz_smear__
#define __fuzz_smear__

namespace Chroma {

//
//  Fuzzed source 
//
//

void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		const LatticeColorVector & psi, 
		LatticeColorVector & psifuzz, 
		int length, int j_decay) ;


void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		const LatticePropagator  & psi, 
		LatticePropagator& psifuzz, 
		int length, int j_decay) ;

void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		const LatticeFermion  & psi, 
		LatticeFermion& psifuzz, 
		int length, int j_decay) ; 


void fuzz_smear(const multi1d<LatticeColorMatrix>& ufuzz, 
		const LatticeStaggeredFermion  & psi, 
		LatticeStaggeredFermion& psifuzz, 
		int length, int j_decay) ; 

}  // end namespace Chroma

#endif
