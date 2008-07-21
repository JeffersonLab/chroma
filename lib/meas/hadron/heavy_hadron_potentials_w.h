// $Id: heavy_hadron_potentials_w.h,v 1.1 2008-07-21 17:59:14 kostas Exp $ 
/*! \file
 *  \brief Potential between 2 heavy hadrons : Detmold
 *
 * Includes Lambda-b's etc
 */

#ifndef __HHPOT_h__
#define __HHPOT_h__

using namespace Chroma; 

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



  // pre-declaration of initialisation function
multi1d<DComplex> makeQllBlock(multi1d<DPropagator> Q1,
			       multi1d<DPropagator> Q2,
			       const SpinMatrix S,
			       multiNd<DComplex> HQ,
			       int size, 
			       int length);

multi1d<DComplex> makeHeavyMesonBlock(multi1d<DPropagator> Q,
				 const SpinMatrix S,
				 multi1d<ColorMatrix> HQ,
				 const SpinMatrix HQspin,
				 int size, 
				 int length);


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
	   const SpinMatrix S,
	   multiNd<DComplex> HQ)
    {
      Nt = len;
      mysize = Nt * Nc * Nc * Nc * Ns * Ns;
      data.resize(mysize);
      data = 0;    

      data = makeQllBlock(Q1, Q2, S, HQ, mysize, Nt);
    };
  
  ~QllBlock(){};
  
  DComplex& operator()(int t, int a, int b, int c, int alpha, int beta){
    
    return data[a + Nc*(b + Nc*(c + Nc* (alpha + Ns *(t + Nt* beta))))];
  };

  DComplex& operator()(multi1d<int> ind){

    if (ind.size() != 6)
    {
      QDPIO::cerr << "Incorrect QllBlock index" << endl;
      QDP_abort(1);
    }
    
    int t=ind[0]; int a=ind[1]; int b=ind[2]; 
    int c=ind[3]; int alpha=ind[4]; int beta=ind[5];
    
    return data[a + Nc*(b + Nc*(c + Nc* (alpha + Ns *(t + Nt* beta))))];
  };

  int size(){return mysize;};

  int length() {return Nt;};

  multi1d<DComplex> makeQllBlock( multi1d<DPropagator> Q1,
				  multi1d<DPropagator> Q2,
				  const SpinMatrix S,
				  multiNd<DComplex> HQ,
				  int size, 
				  int length)
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
	     const SpinMatrix S,
	     multi1d<ColorMatrix> HQ,
	     const SpinMatrix HQspin)
    {
      Nt = len;
      mysize = Nt * Nc * Nc * Ns * Ns;
      data.resize(mysize);
      data = makeHeavyMesonBlock(Q, S, HQ, HQspin, mysize, Nt);
    };
  
  ~HeavyMesonBlock(){};
  
  DComplex& operator()(int t, int a, int b, int alpha, int beta){
    
    return data[a + Nc*(b + Nc* (alpha + Ns *(t + Nt* beta)))];
  };

  DComplex& operator()(multi1d<int> ind){

    if (ind.size() != 5)
    {
      QDPIO::cerr << "Incorrect HeavyMesonBlock index" << endl;
      QDP_abort(1);
    }
    
    int t=ind[0]; int a=ind[1]; int b=ind[2]; 
    int alpha=ind[3]; int beta=ind[4];
    
    return data[a + Nc*(b + Nc* (alpha + Ns *(t + Nt* beta)))];
  };

  int size(){return mysize;};

  int length() {return Nt;};

  multi1d<DComplex> makeHeavyMesonBlock(multi1d<DPropagator> Q,
				   const SpinMatrix S,
				   multi1d<ColorMatrix> HQ,
				   const SpinMatrix HQspin,
				   int size, 
				   int Nt)
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




void QllQllPOT(const multi1d<LatticeColorMatrix>& u, 
	       const LatticePropagator& quark1,
	       const LatticePropagator& quark2,
	       const multi1d<int>& src1, 
	       const multi1d<int>& src2, 
	       const SftMom& phases,
	       XMLWriter& xml,
	       const string& xml_group);

multiNd<DComplex> HBQfunc(multi1d < ColorMatrix >  HQ);

multiNd<DComplex> antiHBQfunc(multi1d < ColorMatrix >  HQ);

multi1d<DComplex> c1contract( QllBlock B1,  QllBlock B2,
			      QllBlock B3,  QllBlock B4,
			      QllBlock B5,  QllBlock B6,
			      QllBlock B7,  QllBlock B8,
			     const SpinMatrix S1, const SpinMatrix S2);

multi1d<DComplex> c4contract( QllBlock B1,  QllBlock B2,
			      QllBlock B3,  QllBlock B4,
			      QllBlock B5,  
			     const SpinMatrix S1, const SpinMatrix S2);

multi1d<DComplex> c5contract( QllBlock B1,  QllBlock B2,
			      QllBlock B3,  QllBlock B4,
			      QllBlock B5,  QllBlock B6,
			      QllBlock B7,  QllBlock B8,
			     const SpinMatrix S1, const SpinMatrix S2);

multi1d<DComplex> c6contract( QllBlock B1,  QllBlock B2,
			     const SpinMatrix S1, const SpinMatrix S2);

multi1d<DComplex> c7contract(QllBlock BzU1zD3z0zCG5, QllBlock BzU2zU4zRzCGii, QllBlock BzU3zD3z0zCG5, 
			     QllBlock BzU4zU2zRzCGii, QllBlock BzU4zU4zRzCGii,
			     const SpinMatrix S1, const SpinMatrix S2);

multi1d<DComplex> d1contract(QllBlock BzU1zD1z0zCG5, QllBlock BzU3zD1z0zCG5, 
			     HeavyMesonBlock HzU2zRzC, HeavyMesonBlock HzU4zRzC,
			     const SpinMatrix mesonS1, const SpinMatrix baryonS2);

multi1d<DComplex> d2contract(QllBlock BzU1zU1z0zCGi, QllBlock BzU3zU1z0zCGi, QllBlock BzU1zU3z0zCGi, 
			     HeavyMesonBlock HzU2zRzG5, HeavyMesonBlock HzU4zRzG5,
			     const SpinMatrix mesonS1, const SpinMatrix baryonS2);

multi1d<DComplex> d3contract(QllBlock BzU1zU1z0zCGi, HeavyMesonBlock HzU2zRzG5,
			     const SpinMatrix mesonS1, const SpinMatrix baryonS2);

multi1d<DComplex> m1contract( HeavyMesonBlock H1,  HeavyMesonBlock H2,
			      HeavyMesonBlock H3,  HeavyMesonBlock H4,
			      const SpinMatrix S1, const SpinMatrix S2);

multi1d<DComplex> m2contract( HeavyMesonBlock H1,  HeavyMesonBlock H2,
			      const SpinMatrix S1, const SpinMatrix S2);

multi1d<DComplex> bcontract(HeavyMesonBlock H1,
			    const SpinMatrix S1);

multi1d<DComplex> lambdabcontract(QllBlock H1,
				  const SpinMatrix S1);

multi1d<DComplex> sigmabpluscontract(QllBlock H1,
				     const SpinMatrix S1);

multi2d<DComplex> c1J2corr(QllBlock B1,QllBlock B2,QllBlock B3,
			   QllBlock B4,QllBlock B5,QllBlock B6,
			   QllBlock B7,QllBlock B8,QllBlock B9,
			   QllBlock B10,QllBlock B11,QllBlock B12,
			   QllBlock B13,QllBlock B14,QllBlock B15,
			   QllBlock B16,QllBlock B17,QllBlock B18,
			   QllBlock B19,QllBlock B20,QllBlock B21,
			   QllBlock B22,QllBlock B23,QllBlock B24);

multi2d<DComplex> c4J2corr(QllBlock B1,QllBlock B2,QllBlock B3,
			   QllBlock B4,QllBlock B5,QllBlock B6,
			   QllBlock B7,QllBlock B8,QllBlock B9,
			   QllBlock B10,QllBlock B11,QllBlock B12,
			   QllBlock B13,QllBlock B14,QllBlock B15);

multi2d<DComplex> c5J2corr(QllBlock B1,QllBlock B2,QllBlock B3,
			   QllBlock B4,QllBlock B5,QllBlock B6,
			   QllBlock B7,QllBlock B8,QllBlock B9,
			   QllBlock B10,QllBlock B11,QllBlock B12,
			   QllBlock B13,QllBlock B14,QllBlock B15,
			   QllBlock B16,QllBlock B17,QllBlock B18,
			   QllBlock B19,QllBlock B20,QllBlock B21,
			   QllBlock B22,QllBlock B23,QllBlock B24);

multi2d<DComplex> c6J2corr(QllBlock B1,QllBlock B2,QllBlock B3,
			   QllBlock B4,QllBlock B5,QllBlock B6);


multi2d<DComplex> d2J32corr(QllBlock BzU1zU1z0zCGplus, QllBlock BzU3zU1z0zCGplus, QllBlock BzU1zU3z0zCGplus, 
			    HeavyMesonBlock HzU2zRzGup, HeavyMesonBlock HzU4zRzGup,
			    QllBlock BzU1zU1z0zCG3, QllBlock BzU3zU1z0zCG3, QllBlock BzU1zU3z0zCG3, 
			    QllBlock BzU1zU1z0zCGminus, QllBlock BzU3zU1z0zCGminus, QllBlock BzU1zU3z0zCGminus, 
			    HeavyMesonBlock HzU2zRzGdown, HeavyMesonBlock HzU4zRzGdown);

multi2d<DComplex> d3J32corr(QllBlock BzU1zU1z0zCGplus, 
			    HeavyMesonBlock HzU2zRzGup, 
			    QllBlock BzU1zU1z0zCG3, 
			    QllBlock BzU1zU1z0zCGminus,
			    HeavyMesonBlock HzU2zRzGdown);

#endif
