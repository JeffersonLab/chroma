// $Id: fuzz_smear.h,v 1.3 2004-11-20 14:54:56 mcneile Exp $

#ifndef __fuzz_smear__
#define __fuzz_smear__

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

#endif
