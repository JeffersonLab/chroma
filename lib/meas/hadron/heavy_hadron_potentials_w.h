// -*- C++ -*-
// $id: heavy_hadron_potentials_w.h,v 1.1 2008/07/21 17:59:14 kostas Exp $ 
/*! \file
 *  \brief Potential between 2 heavy hadrons : Detmold
 *
 * Includes Lambda-b's etc
 */

#ifndef __HHPOT_h__
#define __HHPOT_h__

#include "chromabase.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{
 
  //! Heavy hadron potentials
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   *
   * \param u                  gauge field (Read) 
   * \param quark1             quark propagator ( Read )
   * \param quark2             quark propagator ( Read )
   * \param src1               cartesian coordinates of one source ( Read )
   * \param src2               cartesian coordinates of the other source ( Read )
   * \param phases             object holds list of momenta and Fourier phases ( Read )
   * \param xml                xml file object ( Read )
   * \param xml_group          group name for xml data ( Read )
   *
   */
  void QllQllPOT(const multi1d<LatticeColorMatrix>& u, 
		 const LatticePropagator& quark1,
		 const LatticePropagator& quark2,
		 const multi1d<int>& src1, 
		 const multi1d<int>& src2, 
		 const SftMom& phases,
		 XMLWriter& xml,
		 const string& xml_group);


  // Data classes
  class QllBlock{
  private:
    multi1d<DComplex> data;
    int mysize;
    int Nt;

  public:

    QllBlock(int len){
      Nt = len;
      mysize = Nt * Nc * Nc * Nc * Ns * Ns;
      data.resize(mysize);
      data = 0;    
    };
  
    QllBlock(int len,
	     multi1d<DPropagator> Q1,
	     multi1d<DPropagator> Q2,
	     const SpinMatrix& S,
	     multiNd<DComplex> HQ)
      {
	Nt = len;
	mysize = Nt * Nc * Nc * Nc * Ns * Ns;
	data.resize(mysize);
	data = 0;    

	data = makeQllBlock(Q1, Q2, S, HQ, mysize, Nt);
      };
  
    ~QllBlock(){};
  
    const DComplex& operator()(int t, int a, int b, int c, int alpha, int beta) const{
    
      return data[a + Nc*(b + Nc*(c + Nc* (alpha + Ns *(t + Nt* beta))))];
    };

    const DComplex& operator()(multi1d<int> ind) const {

      if (ind.size() != 6)
      {
	QDPIO::cerr << "Incorrect QllBlock index" << endl;
	QDP_abort(1);
      }
    
      int t=ind[0]; int a=ind[1]; int b=ind[2]; 
      int c=ind[3]; int alpha=ind[4]; int beta=ind[5];
    
      return data[a + Nc*(b + Nc*(c + Nc* (alpha + Ns *(t + Nt* beta))))];
    };

    int size() const {return mysize;};

    int length() const {return Nt;};

    multi1d<DComplex> makeQllBlock( multi1d<DPropagator> Q1,
				    multi1d<DPropagator> Q2,
				    const SpinMatrix& S,
				    multiNd<DComplex> HQ,
				    int size, 
				    int length) const
      {
	multi1d<DComplex> result(size);
	multi1d<int> Carray;
	Carray.resize(4);
	SpinMatrix Q1tmp, Q2tmp;
	DComplex thisres, Q1ttmp, Q2ttmp, HQtmp;
	for (int t=0; t<length; t++){
	  Carray[0]=t;
	  for (int a=0; a<Nc; a++){
	    for (int b=0; b<Nc; b++){
	      for (int c=0; c<Nc; c++){
		Carray[3]=c;
		for (int alpha=0; alpha<Nd; alpha++){
		  for (int beta=0; beta<Nd; beta++){
		    thisres = 0;
		    for (int aa=0; aa<Nc; aa++){
		      Carray[1]=aa;
		      for (int bb=0; bb<Nc; bb++){
			Carray[2]=bb;
			for (int rho=0; rho<Nd; rho++){
			  for (int sigma=0; sigma<Nd; sigma++){
			    Q1tmp = peekColor(Q1[t],aa,a);
			    Q2tmp = peekColor(Q2[t],bb,b);
			    HQtmp = HQ[Carray];
			    Q1ttmp = peekSpin(Q1tmp,rho,alpha);
			    Q2ttmp = peekSpin(Q2tmp,sigma,beta);
			    thisres += Q1ttmp * Q2ttmp * peekSpin(S,rho,sigma) * HQtmp;
			  }
			}
		      }
		    }
		    result[a + Nc*(b + Nc*(c + Nc* (alpha + Ns *(t + Nt* beta))))] = thisres;
		  }
		}
	      }
	    }
	  } 
	}
	return result;
      };

  };



  class HeavyMesonBlock{
    // This can be either Qbarl or Qlbar meson
  private:
    multi1d<DComplex> data;
    int mysize;
    int Nt;

  public:
    HeavyMesonBlock(int len){
      Nt = len;
      mysize = Nt * Nc * Nc * Nc * Ns * Ns;
      data.resize(mysize);
      data = 0;    
    };
  
    HeavyMesonBlock(int len,
		    multi1d<DPropagator> Q,
		    const SpinMatrix& S,
		    multi1d<ColorMatrix> HQ,
		    const SpinMatrix& HQspin)
      {
	Nt = len;
	mysize = Nt * Nc * Nc * Ns * Ns;
	data.resize(mysize);
	data = makeHeavyMesonBlock(Q, S, HQ, HQspin, mysize, Nt);
      };
  
    ~HeavyMesonBlock(){};
  
    const DComplex& operator()(int t, int a, int b, int alpha, int beta) const {
    
      return data[a + Nc*(b + Nc* (alpha + Ns *(t + Nt* beta)))];
    };

    const DComplex& operator()(multi1d<int> ind) const {

      if (ind.size() != 5)
      {
	QDPIO::cerr << "Incorrect HeavyMesonBlock index" << endl;
	QDP_abort(1);
      }
    
      int t=ind[0]; int a=ind[1]; int b=ind[2]; 
      int alpha=ind[3]; int beta=ind[4];
    
      return data[a + Nc*(b + Nc* (alpha + Ns *(t + Nt* beta)))];
    };

    int size() const {return mysize;};

    int length() const {return Nt;};

    multi1d<DComplex> makeHeavyMesonBlock(multi1d<DPropagator> Q,
					  const SpinMatrix& S,
					  multi1d<ColorMatrix> HQ,
					  const SpinMatrix& HQspin,
					  int size, 
					  int Nt) const
    // For Qbarl :  HQspin = 0.5 * (g_one - g_one * Gamma(8)); // This is the antiQ projector
    // For Qlbar :  HQspin = 0.5 * (g_one + g_one * Gamma(8)); // This is the Q projector
      {
	multi1d<DComplex> result(size);
	SpinMatrix Qtmp, U2tmp, opg4, g_one;
	g_one = 1.0;
	DComplex thisres, Qttmp, U2ttmp, HQtmp, Stmp;
	result=0;
	for (int t=0; t<Nt; t++){
	  for (int a=0; a<Nc; a++){  // output colour index
	    for (int b=0; b<Nc; b++){ // output colour index
	      for (int alpha=0; alpha<Nd; alpha++){ // output spin index
		for (int beta=0; beta<Nd; beta++){ // output spin index
		  thisres = 0;
		  for (int c=0; c<Nc; c++){               // colour index summed over
		    for (int rho=0; rho<Nd; rho++){          // spin index summed over
		      for (int sigma=0; sigma<Nd; sigma++){  // spin index summed over
			Qtmp = peekColor(Q[t],c,b);        
			HQtmp = peekColor(HQ[t],a,c) * peekSpin(HQspin,alpha,rho);
			Stmp =  peekSpin(S, sigma,rho);
			Qttmp = peekSpin(Qtmp,sigma,beta);
			thisres += Qttmp * Stmp * HQtmp;
		      }
		    }
		  }
		  result[a + Nc*(b + Nc* (alpha + Ns *(t + Nt* beta)))] = thisres; 
		}
	      }
	    }
	  }
	}
	return result;
      };
  
  };



  // Work
  multiNd<DComplex> HBQfunc(const multi1d<ColorMatrix>& HQ);

  multiNd<DComplex> antiHBQfunc(const multi1d<ColorMatrix>& HQ);

  multi1d<DComplex> c1contract(const QllBlock& B1, const QllBlock& B2,
				const QllBlock& B3, const QllBlock& B4,
				const QllBlock& B5,  const QllBlock& B6,
				const QllBlock& B7,  const QllBlock& B8,
				const SpinMatrix& S1, const SpinMatrix& S2);

  multi1d<DComplex> c4contract( const QllBlock& B1,  const QllBlock& B2,
				const QllBlock& B3,  const QllBlock& B4,
				const QllBlock& B5,  
				const SpinMatrix& S1, const SpinMatrix& S2);

  multi1d<DComplex> c5contract( const QllBlock& B1,  const QllBlock& B2,
				const QllBlock& B3,  const QllBlock& B4,
				const QllBlock& B5,  const QllBlock& B6,
				const QllBlock& B7,  const QllBlock& B8,
				const SpinMatrix& S1, const SpinMatrix& S2);

  multi1d<DComplex> c6contract( const QllBlock& B1,  const QllBlock& B2,
				const SpinMatrix& S1, const SpinMatrix& S2);

  multi1d<DComplex> c7contract(const QllBlock& BzU1zD3z0zCG5, const QllBlock& BzU2zU4zRzCGii, const QllBlock& BzU3zD3z0zCG5, 
			       const QllBlock& BzU4zU2zRzCGii, const QllBlock& BzU4zU4zRzCGii,
			       const SpinMatrix& S1, const SpinMatrix& S2);

  multi1d<DComplex> d1contract(const QllBlock& BzU1zD1z0zCG5, const QllBlock& BzU3zD1z0zCG5, 
			       const HeavyMesonBlock& HzU2zRzC, const HeavyMesonBlock& HzU4zRzC,
			       const SpinMatrix& mesonS1, const SpinMatrix& baryonS2);

  multi1d<DComplex> d2contract(const QllBlock& BzU1zU1z0zCGi, const QllBlock& BzU3zU1z0zCGi, const QllBlock& BzU1zU3z0zCGi, 
			       const HeavyMesonBlock& HzU2zRzG5, const HeavyMesonBlock& HzU4zRzG5,
			       const SpinMatrix& mesonS1, const SpinMatrix& baryonS2);

  multi1d<DComplex> d3contract(const QllBlock& BzU1zU1z0zCGi, const HeavyMesonBlock& HzU2zRzG5,
			       const SpinMatrix& mesonS1, const SpinMatrix& baryonS2);

  multi1d<DComplex> m1contract( const HeavyMesonBlock& H1,  const HeavyMesonBlock& H2,
				const HeavyMesonBlock& H3,  const HeavyMesonBlock& H4,
				const SpinMatrix& S1, const SpinMatrix& S2);

  multi1d<DComplex> m2contract( const HeavyMesonBlock& H1,  const HeavyMesonBlock& H2,
				const SpinMatrix& S1, const SpinMatrix& S2);

  multi1d<DComplex> bcontract(const HeavyMesonBlock& H1,
			      const SpinMatrix& S1);

  multi1d<DComplex> lambdabcontract(const QllBlock& H1,
				    const SpinMatrix& S1);

  multi1d<DComplex> sigmabpluscontract(const QllBlock& H1,
				       const SpinMatrix& S1);

  multi2d<DComplex> c1J2corr(const QllBlock& B1,const QllBlock& B2,const QllBlock& B3,
			     const QllBlock& B4,const QllBlock& B5,const QllBlock& B6,
			     const QllBlock& B7,const QllBlock& B8,const QllBlock& B9,
			     const QllBlock& B10,const QllBlock& B11,const QllBlock& B12,
			     const QllBlock& B13,const QllBlock& B14,const QllBlock& B15,
			     const QllBlock& B16,const QllBlock& B17,const QllBlock& B18,
			     const QllBlock& B19,const QllBlock& B20,const QllBlock& B21,
			     const QllBlock& B22,const QllBlock& B23,const QllBlock& B24);

  multi2d<DComplex> c4J2corr(const QllBlock& B1,const QllBlock& B2,const QllBlock& B3,
			     const QllBlock& B4,const QllBlock& B5,const QllBlock& B6,
			     const QllBlock& B7,const QllBlock& B8,const QllBlock& B9,
			     const QllBlock& B10,const QllBlock& B11,const QllBlock& B12,
			     const QllBlock& B13,const QllBlock& B14,const QllBlock& B15);

  multi2d<DComplex> c5J2corr(const QllBlock& B1,const QllBlock& B2,const QllBlock& B3,
			     const QllBlock& B4,const QllBlock& B5,const QllBlock& B6,
			     const QllBlock& B7,const QllBlock& B8,const QllBlock& B9,
			     const QllBlock& B10,const QllBlock& B11,const QllBlock& B12,
			     const QllBlock& B13,const QllBlock& B14,const QllBlock& B15,
			     const QllBlock& B16,const QllBlock& B17,const QllBlock& B18,
			     const QllBlock& B19,const QllBlock& B20,const QllBlock& B21,
			     const QllBlock& B22,const QllBlock& B23,const QllBlock& B24);

  multi2d<DComplex> c6J2corr(const QllBlock& B1,const QllBlock& B2,const QllBlock& B3,
			     const QllBlock& B4,const QllBlock& B5,const QllBlock& B6);


  multi2d<DComplex> d2J32corr(const QllBlock& BzU1zU1z0zCGplus, const QllBlock& BzU3zU1z0zCGplus, const QllBlock& BzU1zU3z0zCGplus, 
			      const HeavyMesonBlock& HzU2zRzGup, const HeavyMesonBlock& HzU4zRzGup,
			      const QllBlock& BzU1zU1z0zCG3, const QllBlock& BzU3zU1z0zCG3, const QllBlock& BzU1zU3z0zCG3, 
			      const QllBlock& BzU1zU1z0zCGminus, const QllBlock& BzU3zU1z0zCGminus, const QllBlock& BzU1zU3z0zCGminus, 
			      const HeavyMesonBlock& HzU2zRzGdown, const HeavyMesonBlock& HzU4zRzGdown);

  multi2d<DComplex> d3J32corr(const QllBlock& BzU1zU1z0zCGplus, 
			      const HeavyMesonBlock& HzU2zRzGup, 
			      const QllBlock& BzU1zU1z0zCG3, 
			      const QllBlock& BzU1zU1z0zCGminus,
			      const HeavyMesonBlock& HzU2zRzGdown);

}

#endif
