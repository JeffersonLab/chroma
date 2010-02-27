#ifndef __POISSON_H__
#define __POISSON_H__

#include "chromabase.h"
#include "tower.h"

namespace Chroma {
 
  struct Poisson { 
    ComplexD sst;
    ComplexD tst;
    ComplexD tsst;
    ComplexD ssst;
    ComplexD ttst;
    ComplexD ttsst;
    ComplexD sttst;
    ComplexD tssst; 
    ComplexD stsst;
    ComplexD tttst;
    ComplexD sssst;

    Poisson() 
    {
      sst = tst = tsst = ssst = ttst = ttsst 
	= sttst = tssst = stsst = tttst = sssst = zero;
    }

    Poisson(const multi1d<Tower<LatticeColorMatrix> >& F,
	    const multi1d<Tower<LatticeColorMatrix> >& G,
	    const multi1d<LatticeColorMatrix>& P) {
      computePoisson(F,G,P);
    }

    void computePoisson(const multi1d<Tower<LatticeColorMatrix> >& F,
			const multi1d<Tower<LatticeColorMatrix> >& G,
			const multi1d<LatticeColorMatrix>& P)
    {

      sst = tst = tsst = ssst = ttst = ttsst 
	= sttst = tssst = stsst = tttst = sssst = zero;

      for(int mu=0; mu < Nd; mu++) {
	// {S,{S,T}}  = tr(F1^2) 
	sst       -= Real(2)*sum(trace( F[mu][0]*F[mu][0] ));

	// {T,{S,T}}  = - tr(F2 P)  
	tst       += Real(2)*sum(trace( F[mu][1]*P[mu] ));

	// {T,{S,{S,T}}} = 2 tr(F1 F2)
	tsst      -= Real(2)*sum(trace( F[mu][0]*F[mu][1] ));

	// {S,{S,{S,T}}} = 0 - so don't recompute

	// {T,{T,{S,T}}} = -tr( F3 P) 
	ttst      += Real(1)*sum(trace( F[mu][2]*P[mu]));

	// {T,{T,{S,{S,T}}}} = 2( tr(F1 F3} + tr(F2^2) )
	ttsst      -= Real(4)*sum(trace( F[mu][0]*F[mu][2]));  // F1 F3 term
	ttsst      -= Real(4)*sum(trace( F[mu][1]*F[mu][1]));  // F2^2 term


	// {{S,T},{T,{S,T}}} = -3 tr ( [F1,F2] P )
	//                     -  tr ( [F1,P] [F1,P]  )
	//                     +  tr ( F1,F3 )
	//                     -2 tr ( F2,F2 )
 
	// Commutator 1 = [F1,F2] = F1 F2  - F2 F1
	LatticeColorMatrix commut1 = F[mu][0]*F[mu][1];
	commut1 -= F[mu][1]*F[mu][0];

	// Commutator 2 = {F1, P] = F1 P - P F1 
	LatticeColorMatrix commut2 = F[mu][0]*P[mu];
	commut2 -= P[mu]*F[mu][0];


	sttst += Real(6)*sum(trace(commut1*P[mu]));  // -3 tr ( [F1,F2] P )
	sttst += Real(2)*sum(trace(commut2*commut2));    // -  tr ( [F1,P] [F1,P] )
	sttst -= Real(2)*sum(trace(F[mu][0]*F[mu][2]));  // +  tr ( F1,F3 )
	sttst += Real(4)*sum(trace(F[mu][1]*F[mu][1])); //  -2 tr ( F2,F2 )

	// {T,{S,{S,{S,T}}}} = 0 so dont compute

	// {{S,T}, {S,{S,T}}} = -2 tr (F1 G1)
	stsst -= Real(4)*sum(trace(F[mu][0]*G[mu][0]));
	
	// {T,{T,{T,{S,T}}}} = -tr( F4 P )
	tttst -= Real(4)*sum(trace(F[mu][3]*P[mu]));

	// {S,{S,{S,{S,T}}}} = 0 so don't compute
      }



      
    }
   

    void print(void) { 
      QDPIO::cout << "PoissonBrackets: " << endl;
      QDPIO::cout << "================="  << endl;
      QDPIO::cout << "         {S,{S,T}} = " << sst << endl;
      QDPIO::cout << "         {T,{S,T}} = " << tst << endl;
      QDPIO::cout << "     {T,{S,{S,T}}} = " << tsst << endl;
      QDPIO::cout << "     {S,{S,{S,T}}} = " << ssst << endl;
      QDPIO::cout << "     {T,{T,{S,T}}} = " << ttst << endl;
      QDPIO::cout << " {T,{T,{S,{S,T}}}} = " << ttsst << endl;
      QDPIO::cout << " {{S,T},{T,{S,T}}} = " << sttst << endl;
      QDPIO::cout << " {T,{S,{S,{S,T}}}} = " << tssst << endl;
      QDPIO::cout << " {{S,T},{S,{S,T}}} = " << stsst << endl;
      QDPIO::cout << " {T,{T,{T,{S,T}}}} = " << tttst << endl;
      QDPIO::cout << " {S,{S,{S,{S,T}}}} = " << sssst << endl;
      QDPIO::cout << "=================" << endl;
    }

  };
}

#endif
