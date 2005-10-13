// -*- C++ -*-
// $Id: stout_state.cc,v 2.5 2005-10-13 22:07:06 bjoo Exp $
/*! @file 
 *  @brief Connection State for Stout state (.cpp file)
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/stout_state.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"



namespace Chroma { 

  // Do the force recursion from level i+1, to level i
  void StoutConnectState::deriv_recurse(const multi1d<LatticeColorMatrix>&  F_plus,
		     multi1d<LatticeColorMatrix>& F_minus,
		     const int level ) const
  {
    START_CODE();

    // Things I need
    // C_{\mu} = staple multiplied appropriately by the rho
    // Q_{\mu} = -i taproj( C_mu U_mu^{\dagger} ) 
    // f coefficients (as per exponentiation)
    // b_{ij} coefficients -temporary 
    // B matrices as per eq (69) - temporary 
    // Gamma_mu as per eq (74) - temporary
    // Lambda matrices asper eq(73) 
    multi1d<LatticeColorMatrix> C(Nd);
    multi1d<LatticeColorMatrix> Lambda(Nd);

    // The links at this level (unprimed in the paper).
    const multi1d<LatticeColorMatrix>& u_minus = smeared_links[level];


    // First we need the C-s from which we make the fattened links.
    // The C-s are basically the staples
    for(int mu=0; mu < Nd; mu++) {

      LatticeColorMatrix Q;

      // If we don't smear in a direction, the rho-s would be zero
      // so C would be zero -- to save work, we just set C 
      // to zero here.
      C[mu] = 0;

      // If rho is nonzero in this direction then accumulate the staples
      if( smear_in_this_dirP[mu] == true ) { 

	for(int nu=0; nu < Nd; nu++) { 

	  // Accumulate mu-nu staple
	  if( mu != nu && smear_in_this_dirP[nu] == true ) {

	    // Forward staple
	    //             2
	    //       ^ ---------> 
            //       |          |
            //    1  |          |  3
            //       |          |
            //       |          V
	    //       x          x + mu
	    //
	    C[mu] += rho(mu, nu)*u_minus[nu] * shift(u_minus[mu], FORWARD, nu) * adj(shift(u_minus[nu], FORWARD, mu));
	    

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

	    LatticeColorMatrix tmp_1 = adj(u_minus[nu]) 
	      * u_minus[mu] * shift(u_minus[nu], FORWARD, mu);
	    
	    // and here we shift it 
	    // u_staple(x) += shift(tmp_1_dag(x-nu), BACK, nu)
	    //             += u_dag(x+mu-nu,nu)*u_dag(x-nu,mu)*u(x-nu,nu)
	    C[mu] += rho(mu,nu)*shift(tmp_1, BACKWARD, nu);
	  
	  }  // end if(smear_in_this_dirP[nu]
	} // end for nu

	// C_mu is now created
	// Now I can form the Q[mu]. First I form iQ and multiply by
	// -i. iQ is formed with taproj()
	if(smear_in_this_dirP[mu]) {
	
	  LatticeColorMatrix Omega;
	  Omega = C[mu]*adj(u_minus[mu]); // Q_mu is Omega mu here (eq 2 part 2)

#if 0
	  taproj(Omega);           // This does eq(2 part 1)
	  Q = timesI(Omega);
#else
	  LatticeColorMatrix tmp = adj(Omega) - Omega;
	  LatticeColorMatrix tmp2 = traceColor(tmp);
	  tmp2 *= Real(1)/Real(Nc);
	  tmp -= tmp2;
	  tmp *= Real(0.5);
	  Q = timesI(tmp);

	  // Check Q is traceless hermitian.?
#endif
	}
	else { 
	  Q = 0;
	}

    
	// Now I need the c1, c0 etc etc from the exponentiator
	multi1d<LatticeComplex> f(3);
	
	LatticeColorMatrix QQ = Q*Q;



     
	// This is 
	LatticeReal c0    = real((1.0/3.0) * trace(Q*QQ));
	LatticeReal c1    = real((1.0/2.0) * trace(QQ));

	{
	  LatticeBoolean latboo_c1 = ( c1 < 1.0e-4 );
	  LatticeInteger c1trouble = 
	    where( latboo_c1,
		   1,
		   0);
	  
	  if ( toBool( sum(c1trouble) > 0 ) ) {
	    QDPIO::cout << "There are sites where c1 is very small" << endl;
	    QDPIO::cout << "You may experience numerical instability" << endl;
	  }

	}

	LatticeReal c0abs = fabs(c0);
	LatticeReal c0max = 2 * pow((c1 / 3.0), 1.5);
	LatticeReal theta = acos(c0abs/c0max);
	LatticeReal u     = sqrt(c1 / 3.0) * cos(theta / 3.0);
	LatticeReal w     = sqrt(c1) * sin(theta / 3.0);
	LatticeReal uu    = u*u;
	LatticeReal ww    = w*w;
	LatticeReal cosu  = cos(u);
	LatticeReal cosw  = cos(w);
	LatticeReal sinu  = sin(u);
	LatticeReal sinw  = sin(w);

	// exp(2iu) and exp(-iu)
	LatticeComplex exp2iu = cmplx((2*cosu*cosu - 1), 2*cosu*sinu);
	LatticeComplex expmiu = cmplx(cosu, -sinu);

	// c0 is negative
	LatticeBoolean latboo_c0 = (c0      <      0);
	
	// w is greater than 0.05
	LatticeBoolean latboo_w  = (fabs(w) >   0.05);
	
    
	LatticeReal denom = 9 * uu - ww;

	// xi0 = xi0(w).  Expand xi0 if w is small.
	LatticeReal xi0 = where(latboo_w, 
				sinw/w, 
				1 - (1.0/6.0)*ww*(1-(1.0/20.0)*ww*(1-(1.0/42.0)*ww)));
	
	// f_i = f_i(c0, c1). Expand f_i by c1, if c1 is small.
	f[0] = ((uu - ww) * exp2iu + expmiu 
		* cmplx(8*uu*cosw, 2*u*(3*uu+ww)*xi0))/denom;
		 
	
	f[1] = (2*u*exp2iu - expmiu * cmplx(2*u*cosw, (ww-3*uu)*xi0))/denom;
	
	f[2] = (exp2iu - expmiu * cmplx(cosw, 3*u*xi0))/denom;
	
	// apply everywhere were c0 is negative
	// f_j(-c0, c1) = (-1)^j f*_j(c0, c1)
	f[0] = where(latboo_c0,
		     adj(f[0]),
		     f[0]);
	
	f[1] = where(latboo_c0,
		     -1.0*adj(f[1]),
		     f[1]);
	
	f[2] = where(latboo_c0,
		     adj(f[2]),
		     f[2]);
	

	// xi1 = xi1(w) = cos(w)/w^2 - sin(w)/w^3  
	// Expand xi1 if w is small.
	// Not in paper but maple tells me 12th order expansion is
	// xi1= (-1/3)+ (w^2/30)*(1-(w^2/28)*(1-w^2/54));
	LatticeReal xi1;
	{
	  LatticeReal www = ww*w;
	  xi1 = where(latboo_w, 
		      cosw/ww-sinw/www, 
		      (-1.0/3.0) + (1.0/30.0)*ww*(1-(1.0/28.0)*ww*(1-(1.0/54.0)*ww)));    
	}

	LatticeReal b_denom=2.0*(9.0*uu -ww)*(9.0*uu-ww);
	
	multi1d<LatticeComplex> r_1(3);
	multi1d<LatticeComplex> r_2(3);

	r_1[0]=cmplx(2.0*u, 2.0*(uu-ww))*exp2iu
	  + 2.0*expmiu*( cmplx(8.0*u*cosw, -4.0*uu*cosw)
			 + cmplx(u*(3.0*uu+ww),9.0*uu+ww)*xi0 );
	
	r_1[1]=cmplx(2.0, 4.0*u)*exp2iu
	  + expmiu*cmplx(-2.0*cosw-(ww-3.0*uu)*xi0,
			 2.0*u*cosw+6.0*u*xi0);

	r_1[2]=2.0*timesI(exp2iu)
	  +expmiu*cmplx(-3.0*u*xi0, cosw-3*xi0);


	r_2[0]=-2.0*exp2iu + expmiu*cmplx(-8.0*uu*xi0,
					  2.0*(u*cosw+xi0)+6.0*u*uu*xi1);

	
	r_2[1]= expmiu*cmplx(cosw+xi0-3.0*uu*xi1, 2.0*u*xi0);
	r_2[1] = timesMinusI(r_2[1]);

	r_2[2]=expmiu*cmplx(xi0, 3.0*u*xi1);
	
	multi1d<LatticeComplex> b_1(3);
	multi1d<LatticeComplex> b_2(3);

	for(int j=0; j < 3; j++) { 

	  // This has to be a little more careful
	  b_1[j]=2.0*u*r_1[j]+(3.0*uu-ww)*r_2[j]-2.0*(15.0*uu+ww)*f[j];
	  b_1[j]/=b_denom;

	  b_2[j]=r_1[j]-3.0*u*r_2[j]-24.0*u*f[j];
	  b_2[j]/=b_denom;

	}
	 
	// Now take care of the fact that 
	// c_0 could be negative (a la eq 70)
	//
	// c_0 is negative where latboo_c0 is true
	//
	// where c_0 is negative we just do a 
	// 
	//    b_{i,j}(-c_0, c_1) = (-)^{i+j+1} adj( b_{i,j}(c0, c1) ) 
	
	//i=1, j=0, i+j+1 = 1+0+1=2  (-)^{i+j+1} = (-)^2 = +1
	b_1[0] = where(latboo_c0, 
		       adj(b_1[0]),
		       b_1[0]);

	//i=2, j=0, i+j+1 = 2+0+1=3, (-)^{i+j+1}= (-)^3 = -1
	b_2[0] = where(latboo_c0,
		       -1.0*adj(b_2[0]), 
		       b_2[0]);
	
	//i=1, j=1, i+j+1 = 1+1+1=3  (-)^{i+j+1} = (-)^3= -1
	b_1[1] = where(latboo_c0,
		       -1.0*adj(b_1[1]), 
		       b_1[1]);
	
	//i=2, j=1, i+j+1 = 2+1+1=4  (-)^{i+j+1} = (-)^4= +1
	b_2[1] = where(latboo_c0, 
		       adj(b_2[1]),
		       b_2[1]);
	
	//i=1, j=2, i+j+1 = 1+2+1=4  (-)^{i+j+1} = (-)^4 = +1
	b_1[2] = where(latboo_c0, 
		       adj(b_1[2]),
		       b_1[2]);

	//i=2, j=2, i+j+1 = 1+2+1=5  (-)^{i+j+1} = (-)^5 = -1
	b_2[2] = where(latboo_c0,
		       -1.0*adj(b_2[2]), 
		       b_2[2]);
   

	LatticeColorMatrix B_1 = b_1[0] + b_1[1]*Q + b_1[2]*QQ;
	LatticeColorMatrix B_2 = b_2[0] + b_2[1]*Q + b_2[2]*QQ;

	LatticeComplex g1=trace(F_plus[mu]*B_1*u_minus[mu]);
	LatticeComplex g2=trace(F_plus[mu]*B_2*u_minus[mu]);

	LatticeColorMatrix Gamma=g1*Q + g2*QQ + f[1]*u_minus[mu]*F_plus[mu]
	  + f[2]*Q*u_minus[mu]*F_plus[mu] + f[2]*u_minus[mu]*F_plus[mu]*Q;

 #if 0
	Lambda[mu]=Real(0.5)*(Gamma + adj(Gamma)) - (Real(0.5)/Real(Nc))*trace( Gamma + adj(Gamma));
#else
	Lambda[mu] = Gamma + adj(Gamma);
	LatticeColorMatrix tmp = traceColor(Lambda[mu]);
	tmp *= (Real(1)/Real(Nc));
	Lambda[mu] -= tmp;
	Lambda[mu] *= Real(0.5);
#endif
	F_minus[mu] = F_plus[mu]*(f[0] + f[1]*Q + f[2]*QQ);
	F_minus[mu] += timesI(adj(C[mu])*Lambda[mu]);
      }
      else {
	// If there is no smearing, copy the previous force.
	F_minus[mu] = F_plus[mu];
      }
    }


    for(int mu = 0; mu < Nd; mu++) { 

      if( smear_in_this_dirP[mu] ) { 
	// All the staple terms go here
	// There are six staple terms
	LatticeColorMatrix staple_sum = 0;

	for(int nu = 0; nu < Nd; nu++) { 
	  if( smear_in_this_dirP[nu] == true && nu != mu ) {


	    // Staple 1
	    //
	    // U_nu(x+mu) * U^+_mu(x+nu) U^+_nu(x) Lambda_nu(x)
	    //               2
	    //     ^    <---------^
	    //     ||   |         |                      ^nu
	    //   4 ||   |3        |1                     |
	    //     ||   |         |                       -> mu
	    //          V
	    staple_sum += rho(nu,mu)*shift(u_minus[nu],FORWARD,mu)*adj(shift(u_minus[mu],FORWARD,nu))*adj(u_minus[nu])*Lambda[nu];
	    
	

	    // Staple 2
	    //
	    // U^{+}_nu(x-nu+mu) U^+_mu(x-nu) Lambda_mu(x-nu) U_nu(x-nu)
	    //
	    //           x    ^           | 
	    //                |           | 1
	    //              4 |           |
	    //                |     2     |
	    //        x-nu    <---------- V x-nu+mu
	    //                 ==========>                  
	    //                      3
	    // 
	    //  Construct on x-nu and then shift

       
	    
	    LatticeColorMatrix st_tmp = adj(shift(u_minus[nu],FORWARD,mu))
	      * adj(u_minus[mu])*Lambda[mu]*u_minus[nu];
	    
	    staple_sum  += rho(mu,nu)*shift(st_tmp, BACKWARD, nu);
	    
	    
	    // Staple 3
	    //
	    // U^{+}_nu(x-nu+mu) Lamnda_nu(x-nu+mu) U^{+}_mu(x-nu)U_{nu}(x-nu)
	    //
	    // 
	    //
	    //           x    ^          ^   | 
	    //                |        2 ||  | 1
	    //              4 |          ||  |
	    //                |     3    ||  |
	    //        x-nu    <----------||  V x-nu+mu
	    //                                   
	    //  Construct on x-nu, but first do multiplies of 1 and 2 on 
	    //  every site
	    // 
	    LatticeColorMatrix st_tmp2 = adj(u_minus[nu])*Lambda[nu];
	    st_tmp  = shift(st_tmp2, FORWARD, mu)*adj(u_minus[mu])*u_minus[nu];
	    staple_sum  += rho(nu,mu)*shift(st_tmp, BACKWARD, nu);
	    
	    
	    
	    // Staple 4
	    //
	    // U^{+}_nu(x-nu+mu)U^{+}_mu(x-nu) Lambda_nu(x-nu)U_nu(x-nu)
	    //
	    //                   ^
	    //           x    ^  ||           | 
	    //                |  ||3          | 1
	    //              4 |  ||           |
	    //                |  ||    2      |
	    //             x-nu    <----------V x-nu+mu
	    //                                   
	    // Construct on x-nu
	    //

	    st_tmp = adj(shift(u_minus[nu], FORWARD, mu))
	      * adj(u_minus[mu])
	      * Lambda[nu]
	      * u_minus[nu];
	    
	    staple_sum  -= rho(nu,mu)*shift(st_tmp, BACKWARD, nu);
	    
	    
	    // Staple 5 
	    //
	    // Lambda_nu(x+mu)*U_nu(x+mu)U^{+}_mu(x+nu)U^{+}_nu(x)
	    //
	    //                           3
	    //             x+nu     <---------^  ^
	    //                      |         |  ||
	    //                   4  |         |  || 1
	    //                      |      2  |  ||
	    //              x       V         |  ||  x + mu
	    // 			       
	    //   Construct on x
	    //   Fuse 1 and 2
	    // 
	    st_tmp = Lambda[nu]*u_minus[nu];
	    staple_sum -= rho(nu,mu)*shift(st_tmp,FORWARD, mu)*adj(shift(u_minus[mu],FORWARD, nu))*adj(u_minus[nu]);
	    
	    // Finally staple 6
	    //
	    // Lambda_nu(x+mu)*U_nu(x+mu)U^{+}_mu(x+nu)U^{+}_nu(x)
	    //
	    //                           2
	    //             x+nu     <---------^  
	    //                      =========>|
	    //                   4  |    3    | 1
	    //                      |         |  
	    //              x       V         |   x + mu
	    // 			       
	    //   Construct on x
	    //   Fuse 2 and 3
	    // 
	    st_tmp  = adj(u_minus[mu])*Lambda[mu];
	    staple_sum += rho(mu,nu)*shift(u_minus[nu],FORWARD, mu)*shift(st_tmp, FORWARD, nu)*adj(u_minus[nu]);
	    
	  } // smear in nu dir
	} // end nu loop
    
	F_minus[mu] -= timesI(staple_sum);
      }
    }

    // Done
    END_CODE();
  }


  //! Drive the recursion for the force term
  //  Parameters: F is the force computed with the FAT Links
  void StoutConnectState::deriv(multi1d<LatticeColorMatrix>& F) const
  {
    START_CODE();
    multi1d<LatticeColorMatrix> F_minus(Nd); // Force at one level lower

    for(int level=n_smear-1; level >= 0; level--) {

      // Take the current force, and compute force one level down 
      // Level is index into smeared link array
      deriv_recurse(F,F_minus,level);   
      F = F_minus;
    }



    for(int mu=0; mu < Nd; mu++) { 
      F[mu] = (smeared_links[0])[mu]*F_minus[mu];
    }



    END_CODE();
  }


  void StoutConnectState::smear_links(const multi1d<LatticeColorMatrix>& current, 
				      multi1d<LatticeColorMatrix>& next)
  {
    START_CODE();

    // Construct and add the staples, where smear_this_dirP[mu] is true
    
    for(int mu = 0; mu < Nd; mu++) { 

      if( smear_in_this_dirP[mu] == true ) { 

	// Do smear
	LatticeColorMatrix u_staple = 0;
	for(int nu = 0; nu < Nd; ++nu) {
	  
	  if( nu != mu && smear_in_this_dirP[nu] == true ) {
	    
	    // Forward staple
	    u_staple += rho(mu,nu)* (current[nu] * shift(current[mu], FORWARD, nu) * adj(shift(current[nu], FORWARD, mu)));
	    
	    // Backward staple
	    // tmp_1(x) = u_dag(x,nu)*u(x,mu)*u(x+mu,nu)
	    LatticeColorMatrix tmp_1 = ( adj(current[nu]) * current[mu] * shift(current[nu], FORWARD, mu) );
	    
	    // u_staple(x) += tmp_1_dag(x-nu)
	    //             += u_dag(x+mu-nu,nu)*u_dag(x-nu,mu)*u(x-nu,nu)
	    u_staple += rho(mu,nu)*shift(tmp_1, BACKWARD, nu);
	  }
	}

	// The proto smeared link
	LatticeColorMatrix u_tmp = u_staple * adj(current[mu]);
	
	// Take the trace-less anti-hermitian projection of the staple
	taproj(u_tmp);
	
	// Exactly exponentiate the Lie Algebra matrix
	// Now u_tmp = exp(iQ)
	expmat(u_tmp,EXP_EXACT);
	
	next[mu]=u_tmp*current[mu];
      }
      else { 
	next[mu]=current[mu]; // No smearing in this dir. Just copyq
      }
    }

    END_CODE();
  }




  StoutConnectState::StoutConnectState(const multi1d<LatticeColorMatrix>& u_,
				       const multi2d<Real>& sm_fact_,
				       const int n_smear_, 
				       const multi1d<bool>& smear_in_this_dirP_)
  {
    START_CODE();
    create(u_, sm_fact_, n_smear_, smear_in_this_dirP_);
    END_CODE();
  }


  //! Explicitly specify smearing factor tensor
  StoutConnectState::StoutConnectState(const multi1d<LatticeColorMatrix>& u_,
				const multi2d<Real>& sm_fact_,
				const int n_smear_) 
  {
    START_CODE();
    // Specify all smearings but no mask. Assume 
    // smearing desired in all directions
    multi1d<bool> smear_in_this_dirP_aux(Nd);
    for(int mu=0; mu < Nd; mu++) { 
      smear_in_this_dirP_aux[mu] = true;
    }
    
    create(u_, sm_fact_, n_smear_, smear_in_this_dirP_aux);
    END_CODE();
  }

     //! Construct isotropic smearing in all 4 directions
  StoutConnectState::StoutConnectState(const multi1d<LatticeColorMatrix>& u_, 
				       const Real& sm_fact_, 
				       const int   n_smear_) 
  { 
    START_CODE();

    multi2d<Real> sm_fact_array(Nd, Nd);
    multi1d<bool> smear_in_this_dirP_aux(Nd);
    
    // For each (mu,nu) set sm_fact_array(mu,nu)=sm_fact
    // (Isotropy). Since mu != nu ever, we set those
    // to zero for safety
    for(int mu=0; mu < Nd; mu++) { 
      
      for(int nu=0; nu < Nd; nu++) { 
	
	if( mu==nu ) { 
	  sm_fact_array[mu][nu] = 0;
	}
	else { 
	  sm_fact_array[mu][nu] = sm_fact_;
	}
	
      }
      
      // smearing in all 4 directions
      smear_in_this_dirP_aux[mu]=true;
    }
    
    // call the create
    create(u_, sm_fact_array, n_smear_, smear_in_this_dirP_aux);
    END_CODE();
  }

  //! Construct isotopic smearing in directions orthogonal to orthog dir
  StoutConnectState::StoutConnectState(const multi1d<LatticeColorMatrix>& u_,
				       const Real& sm_fact_, 
				       const int   n_smear_,
				       const int   orthog_dir) 
  {
    START_CODE();
    multi2d<Real> sm_fact_array(Nd, Nd);
    multi1d<bool> smear_in_this_dirP_aux(Nd);
    
    // For each (mu,nu) set sm_fact_array(mu,nu)=sm_fact
    // (Isotropy). Since mu != nu ever, we set those
    // to zero for safety
    for(int mu=0; mu < Nd; mu++) { 
      
      for(int nu=0; nu < Nd; nu++) { 
	if( mu != nu ) {
	  sm_fact_array[mu][nu] = sm_fact_;
	}
	else {
	  // Set the rho to 0 if mu==nu
	  sm_fact_array[mu][nu] = 0;
	}
      }
      
      // Mask out the orthog dir
      if( mu == orthog_dir ) {  // Direction is same as orthog dir
	smear_in_this_dirP_aux[mu]=false;
      }
      else {                    // Direction orthogonal to orthog dir
	smear_in_this_dirP_aux[mu]=true;
      }

#if 0
      QDPIO::cout << "Smearing in direction " << mu << " ? ";
      if ( smear_in_this_dirP_aux[mu] == true ) { 
	QDPIO::cout << " yes" << endl;
      }
      else {
	QDPIO::cout << " no" << endl;
      }
#endif 
 
    }
      
    // call the create
    create(u_, sm_fact_array, n_smear_, smear_in_this_dirP_aux);
    END_CODE();
  }

  // create function
  void StoutConnectState::create(const multi1d<LatticeColorMatrix>& u_,
	      const multi2d<Real>& sm_fact_,
	      const int n_smear_, 
	      const multi1d<bool>& smear_in_this_dirP_) 
  { 
    START_CODE();

    // Copy smearing factors
    rho.resize(Nd, Nd);
    rho = sm_fact_;
    
    // set n_smear
    n_smear = n_smear_;

    // Copy the direction maske
    smear_in_this_dirP.resize(Nd);
    smear_in_this_dirP = smear_in_this_dirP_;

    // Allocate smeared links
    smeared_links.resize(n_smear + 1);
    for(int i=0; i <= n_smear; i++) { 
      smeared_links[i].resize(Nd);
    }
    

    // Copy thin links into smeared_links[0]
    for(int mu=0; mu < Nd; mu++) { 
      smeared_links[0][mu] = u_[mu];
    }

    // Iterate up the smearings
    for(int i=1; i <= n_smear; i++) { 
      smear_links(smeared_links[i-1], smeared_links[i]);
    }
    END_CODE();
  }

}; // End namespace Chroma
