// $Id: fuzz_smear.h,v 1.2 2004-02-08 14:35:07 mcneile Exp $

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


#endif
