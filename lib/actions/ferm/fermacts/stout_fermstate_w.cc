// -*- C++ -*-
// $Id: stout_fermstate_w.cc,v 1.6 2006-08-15 13:17:24 bjoo Exp $
/*! @file 
 *  @brief Connection State for Stout state (.cpp file)
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/stout_fermstate_w.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"

#include "actions/ferm/fermacts/ferm_createstate_factory_w.h"
#include "actions/ferm/fermacts/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"


namespace Chroma 
{ 

  /*! \ingroup fermacts */
  namespace CreateStoutFermStateEnv 
  { 
    CreateFermState<LatticeFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createFerm(XMLReader& xml, 
							      const std::string& path) 
    {
      return new CreateStoutFermState(WilsonTypeFermBCEnv::reader(xml, path),
				      StoutFermStateParams(xml, path));
    }

    const std::string name = "STOUT_FERM_STATE";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheCreateFermStateFactory::Instance().registerObject(name, 
									  createFerm);
      return foo;
    }

    const bool registered = registerAll();
  }

  void StoutFermState::getFsAndBs(const LatticeDouble& c0,
				  const LatticeDouble& c1,
				  multi1d<LatticeDComplex>& f,
				  multi1d<LatticeDComplex>& b1,
				  multi1d<LatticeDComplex>& b2,
				  bool dobs) const
  {

    LatticeBoolean latboo_c0 = (c0 < Double(0));
    LatticeDouble c0abs = fabs(c0);
    LatticeDouble c0max = Double(2)*pow( c1/Double(3), Double(1.5));
    
    LatticeDouble theta;
    LatticeBoolean latboo_c0max =( c0max == 0 );
    
    theta = acos( c0abs/c0max );
  
    LatticeDouble u = sqrt(c1/Double(3))*cos(theta/Double(3));
    LatticeDouble w = sqrt(c1)*sin(theta/Double(3));
    
    LatticeDouble u_sq = u*u;
    LatticeDouble w_sq = w*w;
    
    LatticeDouble xi0,xi1;
    {
      LatticeBoolean latboo_w = (fabs(w) < Double(0.05));
      xi0 = where( latboo_w,
		   Double(1) - (Double(1)/Double(6))*w_sq*(Double(1) - (Double(1)/Double(20))*w_sq*(Double(1) - (Double(1)/Double(42))*w_sq)),
		   (sin(w)/w) );

      if( dobs==true) {
	xi1 = where( latboo_w,
		     Double(-1)*((Double(1)/Double(3)) - ( Double(1)/Double(30) )*w_sq*(Double(1) - (Double(1)/Double(28))*w_sq*(Double(1)- (Double(1)/Double(54))*w_sq))),
		     cos(w)/w_sq - sin(w)/(w_sq*w) );
      }
    }

    
    LatticeDouble cosu = cos(u);
    LatticeDouble sinu = sin(u);
    LatticeDouble cosw = cos(w);
    LatticeDouble sinw = sin(w);
    LatticeDouble sin2u = sin(2*u);
    LatticeDouble cos2u = cos(2*u);

    // exp(2iu) and exp(-iu)
    //LatticeDComplex exp2iu = cmplx(( 2*cosu*cosu - 1), 2*cosu*sinu);
    LatticeDComplex exp2iu = cmplx( cos2u, sin2u );
    LatticeDComplex expmiu = cmplx(cosu, -sinu);

    LatticeDouble denum = 9*u_sq - w_sq;
    f.resize(3);

    // f_i = f_i(c0, c1). Expand f_i by c1, if c1 is small.
    f[0] = ((u_sq - w_sq) * exp2iu + expmiu * cmplx(8*u_sq*cosw, 2*u*(3*u_sq+w_sq)*xi0))/denum;
    
	
    f[1] = (2*u*exp2iu - expmiu * cmplx(2*u*cosw, (w_sq-3*u_sq)*xi0))/denum;
    
    f[2] = (exp2iu - expmiu * cmplx(cosw, 3*u*xi0))/denum;

    if( dobs == true ) {
      
      multi1d<LatticeDComplex> r_1(3);
      multi1d<LatticeDComplex> r_2(3);
      
      r_1[0]=Double(2)*cmplx(u, u_sq-w_sq)*exp2iu
	+ 2.0*expmiu*( cmplx(8.0*u*cosw, -4.0*u_sq*cosw)
		       + cmplx(u*(3.0*u_sq+w_sq),9.0*u_sq+w_sq)*xi0 );
      
      r_1[1]=cmplx(2.0, 4.0*u)*exp2iu
	+ expmiu*cmplx(-2.0*cosw-(w_sq-3.0*u_sq)*xi0,
		       2.0*u*cosw+6.0*u*xi0);
      
      r_1[2]=2.0*timesI(exp2iu)
	+expmiu*cmplx(-3.0*u*xi0, cosw-3*xi0);
      
      
      r_2[0]=-2.0*exp2iu + 2*cmplx(0,u)*expmiu*cmplx(cosw+xi0+3*u_sq*xi1,
						     4*u*xi0);
      
      r_2[1]= expmiu*cmplx(cosw+xi0-3.0*u_sq*xi1, 2.0*u*xi0);
      r_2[1] = timesMinusI(r_2[1]);
      
      r_2[2]=expmiu*cmplx(xi0, -3.0*u*xi1);
    
      b1.resize(3);
      b2.resize(3);
      
      LatticeDouble b_denum=2.0*(9.0*u_sq -w_sq)*(9.0*u_sq-w_sq);
      
      for(int j=0; j < 3; j++) { 
	
	// This has to be a little more careful
	
	b1[j]=( 2.0*u*r_1[j]+(3.0*u_sq-w_sq)*r_2[j]-2.0*(15.0*u_sq+w_sq)*f[j] )/b_denum;
	b2[j]=( r_1[j]-3.0*u*r_2[j]-24.0*u*f[j] )/b_denum;
	
	
      }
      
      b1[0] = where(latboo_c0,
	 	    conj(b1[0]),
		    b1[0]);
      
      b1[1] = where(latboo_c0,
		    -1*conj(b1[1]),
		    b1[1]);
      
      b1[2] = where(latboo_c0,
		    conj(b1[2]),
		    b1[2]);
      
      b2[0] = where(latboo_c0,
		    -1*conj(b2[0]),
		    b2[0]);
      
      b2[1] = where(latboo_c0,
		    conj(b2[1]),
		    b2[1]);
    
      b2[2] = where(latboo_c0,
		    -1*conj(b2[2]),
		    b2[2]);

    }

    f[0] = where(latboo_c0,
		 conj(f[0]),
		 f[0]);
    
    f[1] = where(latboo_c0,
		 -1*conj(f[1]),
		 f[1]);
    
    f[2] = where(latboo_c0,
		 conj(f[2]),
		 f[2]);
    
  }


  void StoutFermState::getQsandCs(const multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& Q, 
				  LatticeColorMatrix& QQ,
				  LatticeColorMatrix& C, 
				  LatticeDouble& c0, 
				  LatticeDouble& c1, int mu) const
  {
   
    C = zero;

    // If rho is nonzero in this direction then accumulate the staples
    for(int nu=0; nu < Nd; nu++) { 
	
      // Accumulate mu-nu staple
      if( mu != nu ) {
	
	// Forward staple
	//             2
	//       ^ ---------> 
	//       |          |
	//    1  |          |  3
	//       |          |
	//       |          V
	//       x          x + mu
	//
	C += params.rho(mu, nu)*u[nu] * shift(u[mu], FORWARD, nu) * adj(shift(u[nu], FORWARD, mu));
	  
	  
	// Backward staple
	//             
	//       |          ^ 
	//       |          |
	//    1  |          |  3
	//       |     2    |
	//       V--------->|          
	//       x-nu        x - nu  + mu
	//
	//
	//  we construct it on x-nu and shift it up to x, 
	// (with a backward shift)
	
	// This is the staple on x-nu:
	// tmp_1(x) = u_dag(x,nu)*u(x,mu)*u(x+mu,nu)
	
	LatticeColorMatrix tmp = adj(u[nu])*u[mu] * shift(u[nu], FORWARD, mu);
	// and here we shift it 
	// u_staple(x) += shift(tmp_1_dag(x-nu), BACK, nu)
	//             += u_dag(x+mu-nu,nu)*u_dag(x-nu,mu)*u(x-nu,nu)
	C += params.rho(mu,nu)*shift(tmp, BACKWARD, nu);

	}
      
      }

      // Now I can form the Q
      LatticeColorMatrix Omega;
      Omega = C*adj(u[mu]); // Q_mu is Omega mu here (eq 2 part 2)
      
      LatticeColorMatrix tmp2 = adj(Omega) - Omega;
      LatticeColorMatrix tmp3 = trace(tmp2);
      tmp3 *= Real(1)/Real(Nc);
      tmp2 -= tmp3;
      tmp2 *= Real(0.5);
      Q = timesI(tmp2);
      QQ = Q*Q;

      LatticeColorMatrix QQQ=QQ*Q;

      c0    = (Double(1)/Double(3)) * real(trace(QQQ));  // eq 13
      c1    = (Double(1)/Double(2)) * real(trace(QQ));	 // eq 15 
      
  }

  // Do the force recursion from level i+1, to level i
  // The input fat_force F is modified.
  void StoutFermState::deriv_recurse(multi1d<LatticeColorMatrix>& F,
				     const int level ) const
  {
    START_CODE();

    // Things I need
    // C_{\mu} = staple multiplied appropriately by the rho
    // Lambda matrices asper eq(73) 
    multi1d<LatticeColorMatrix> F_plus(Nd);

    // Save the fat force
    F_plus = F;


    multi1d<LatticeColorMatrix> Lambda(Nd);
    multi1d<LatticeColorMatrix> C(Nd);

    // The links at this level (unprimed in the paper).
    const multi1d<LatticeColorMatrix>& u = smeared_links[level];

    for(int mu=0; mu < Nd; mu++) 
    {
      LatticeColorMatrix Q,QQ;   // This is the C U^{dag}_mu suitably antisymmetrized
      LatticeDouble c0, c1;

      // Get Q, Q^2, C, c0 and c1 -- this code is the same as used in stout_smear()
      getQsandCs(u, Q, QQ, C[mu], c0, c1, mu);

      // Now work the f-s and b-s
      multi1d<LatticeDComplex> f;
      multi1d<LatticeDComplex> b_1;
      multi1d<LatticeDComplex> b_2;

      // Get the fs and bs  -- does internal resize to make them arrays of length 3
      getFsAndBs(c0, c1, f, b_1, b_2, true);

 
      LatticeColorMatrix B_1 = b_1[0] + b_1[1]*Q + b_1[2]*QQ;
      LatticeColorMatrix B_2 = b_2[0] + b_2[1]*Q + b_2[2]*QQ;

      
      // Construct the Gamma ( eq 74 and 73 )
      LatticeColorMatrix USigma = u[mu]*F_plus[mu];
      LatticeColorMatrix Gamma = f[1]*USigma + f[2]*(USigma*Q + Q*USigma)
	+ trace(B_1*USigma)*Q
        + trace(B_2*USigma)*QQ;

      // Take the traceless hermitian part to form Lambda_mu (eq 72)
      Lambda[mu] = Gamma + adj(Gamma);    // Make it hermitian
      LatticeColorMatrix tmp3 = (Double(1)/Double(Nc))*trace(Lambda[mu]); // Subtract off the trace
      Lambda[mu] -= tmp3;
      Lambda[mu] *= Double(0.5);         // overall factor of 1/2

      // The first 3 terms of eq 75
      // Now the Fat force * the exp(iQ)
      F[mu]  = F_plus[mu]*(f[0] + f[1]*Q + f[2]*QQ);
           
    }
    
    // At this point we should have 
    //
    //  F[mu] = F_plus[mu]*exp(iQ) 
    //
    //  We need the 8 staple terms left in dOmega/dU (last 6 terms in eq 75 + the iC{+}Lambda
    //  term in eq 75 which in reality just covers 2 staples.
 
    //  We have to make this a separate loop from the above, because we need to know the 
    //  Lambda[mu] and [nu] for all the avaliable mu-nu combinations
    for(int mu = 0; mu < Nd; mu++) { 
      LatticeColorMatrix staple_sum = zero;
      LatticeColorMatrix staple_sum_dag = adj(C[mu])*Lambda[mu];
      for(int nu = 0; nu < Nd; nu++) { 
	if(mu != nu) { 
	  // Staple 1
	  //
	  // rho_nu_mu* U_nu(x+mu) * U^+_mu(x+nu) U^+_nu(x) Lambda_nu(x)
	  staple_sum += params.rho(nu,mu)  
	    *shift(u[nu],FORWARD,mu)              //   U_nu(x+mu)     
	    *shift(adj(u[mu]),FORWARD,nu)         //   U^+_mu(x+nu)
	    *adj(u[nu])                           //   U^+_nu(x)
	    *Lambda[nu];                                //   Lambda_nu(x)
	  
	
	  
	  // Staple 2
	  //
	  // rho_mu_nu * U^{+}_nu(x-nu+mu) U^+_mu(x-nu) Lambda_mu(x-nu) U_nu(x-nu)
	  staple_sum += params.rho(mu,nu)
	    *shift(shift(adj(u[nu]), BACKWARD,nu), FORWARD,mu) //  U^{+}_nu(x-nu+mu)
	    *shift(adj(u[mu]), BACKWARD, nu)                   //  U^{+}_mu(x-nu)
	    *shift(Lambda[mu],  BACKWARD, nu)                        //  Lambda_mu(x-nu)
	    *shift(u[nu], BACKWARD, nu);                       //  U_nu(x-nu)

	  
	  // Staple 3
	  // 
	  // rho_nu_mu * U^{+}_nu(x-nu+mu) Lambda_nu(x-nu+mu) U^{+}_mu(x-nu )U_{nu}(x-nu)
	  staple_sum += params.rho(nu,mu)
	    *shift( shift( adj(u[nu]), BACKWARD, nu), FORWARD, mu)   // U^{+}_nu(x-nu+mu)
	    *shift( shift( Lambda[nu], BACKWARD, nu), FORWARD, mu)         // Lambda_nu (x -nu + mu)
	    *shift( adj(u[mu]), BACKWARD, nu)                        // U^{+}_mu(x-nu)
	    *shift( u[nu], BACKWARD, nu);                            // U_nu(x-nu)
	  
	  
	  
	  // Staple 4
	  //
	  // - rho_nu_mu U^{+}_nu(x-nu+mu)U^{+}_mu(x-nu) Lambda_nu(x-nu) U_nu(x-nu)
	  staple_sum_dag += params.rho(nu,mu)
	    *shift( shift( adj(u[nu]), BACKWARD, nu), FORWARD, mu) // U^{+}_nu(x-nu+mu)
	    *shift( adj(u[mu]), BACKWARD, nu)                       // U^{+}_mu(x-nu)
	    *shift( Lambda[nu], BACKWARD, nu )                            // Lambda_nu(x-nu) 
	    *shift( u[nu], BACKWARD, nu);                           // U_nu(x-nu)
	  
	  
	  // Staple 5 
	  //
	  // - rho(nu,mu) * Lambda_nu(x+mu)*U_nu(x+mu)U^{+}_mu(x+nu)U^{+}_nu(x)
	  //
	  // 			       
	  //   Construct on x
	  //   Fuse 1 and 2
	  // 
	  staple_sum_dag += params.rho(nu,mu)
	    * shift( Lambda[nu], FORWARD, mu)                     // Lambda_nu( x + mu )
	    * shift( u[nu], FORWARD, mu)                    // U_nu(x + mu)
	    * shift(adj(u[mu]), FORWARD, nu)              // U^{+}_mu(x+nu)
	    * adj( u[nu] );                                 // U^{+}_nu(x)

	  // Finally staple 6
	  //
	  // rho(mu,nu) * U_nu(x + mu) U^{+}_mu(x+nu) Lambda_mu(x + nu) U^{+}_nu(x)
	  //
	  staple_sum += params.rho(mu,nu)
	    * shift( u[nu], FORWARD, mu)         //  U_nu(x + mu)
	    * shift( adj(u[mu]), FORWARD, nu)    //  U^{+}_mu(x+nu)
	    * shift( Lambda[mu], FORWARD, nu)          //  Lambda_mu(x + nu)
	    * adj( u[nu] );                      //  U^{+}_nu(x)
	}  
      } // end nu loop
	
      // The siz terms from staple sum
      F[mu] += timesI( staple_sum_dag );
      F[mu] -= timesI(staple_sum);
     
    }
  

    // Done
    END_CODE();
  }


  //! Drive the recursion for the force term
  //  Parameters: F is the force computed with the FAT Links
  void StoutFermState::fatForceToThin(const multi1d<LatticeColorMatrix>& F_fat, multi1d<LatticeColorMatrix>& F_thin) const
  {
    START_CODE();
    
    F_thin.resize(Nd);

    F_thin = F_fat;
    fbc->zero(F_thin);

    // Now if the state is smeared recurse down.

    for(int level=params.n_smear; level > 0; level--) {
      deriv_recurse(F_thin,level-1);
      fbc->zero(F_thin);
    }


    END_CODE();
  }


  void StoutFermState::smear_links(const multi1d<LatticeColorMatrix>& current, 
				   multi1d<LatticeColorMatrix>& next)
  {
    START_CODE();


    for(int mu = 0; mu < Nd; mu++) { 

      LatticeColorMatrix Q, QQ;
      LatticeColorMatrix C; // Dummy - not used here
      LatticeDouble c0, c1;

      // Q contains the staple term. C is a throwaway
      getQsandCs(current, Q, QQ, C, c0, c1, mu);

      // Now compute the f's -- use the same function as for computing the fs, bs etc in derivative
      // but don't compute the b-'s
      multi1d<LatticeDComplex> f;   // routine will resize these
      multi1d<LatticeDComplex> b_1; // Dummy - not used      -- throwaway -- won't even get resized
      multi1d<LatticeDComplex> b_2; // Dummy - not used here -- throwaway -- won't even get resized


      getFsAndBs(c0,c1,f,b_1,b_2,false);   // This routine computes the f-s
      
      // Assemble the stout links exp(iQ)U_{mu} 
      next[mu]=(f[0] + f[1]*Q + f[2]*QQ)*current[mu];      
    }
    
    END_CODE();
  }



  // Constructor only from a parameter structure
  StoutFermState::StoutFermState(Handle< FermBC<T,P,Q> > fbc_, 
				 const StoutFermStateParams& p_,
				 const multi1d<LatticeColorMatrix>& u_)
  {
    START_CODE();
    QDPIO::cout << __func__ << ": entering" << endl;
    create(fbc_, p_, u_);
    QDPIO::cout << __func__ << ": exiting" << endl;
    END_CODE();
  }


  // create function
  void StoutFermState::create(Handle< FermBC<T,P,Q> > fbc_, 
			      const StoutFermStateParams& p_,
			      const multi1d<LatticeColorMatrix>& u_)
  { 
    START_CODE();

    fbc = fbc_;
    params = p_;

    // Allocate smeared and thin links
    smeared_links.resize(params.n_smear + 1);
    for(int i=0; i <= params.n_smear; i++) { 
      smeared_links[i].resize(Nd);
    }
    

    // Copy thin links into smeared_links[0]
    for(int mu=0; mu < Nd; mu++) { 
      (smeared_links[0])[mu] = u_[mu];
    }
    fbc->modify( smeared_links[0] );    
    

    // Iterate up the smearings
    for(int i=1; i <= params.n_smear; i++) { 


      smear_links(smeared_links[i-1], smeared_links[i]);
      // If the fermbc-s are nontrivial 
      // ie the underlying gauge fields are nontrivial
      // such as SF BC-s then apply them at every level
      fbc->modify( smeared_links[i] );
    }

    // Always apply at the top level 
    // (essentially what this does is it takes care
    //  where modify only cares about fermion bc-s
    //  and assumes an underlying trivially periodic gauge bc)
 
    

    END_CODE();
  }


} // End namespace Chroma
