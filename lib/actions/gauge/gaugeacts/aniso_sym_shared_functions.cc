#include "actions/gauge/gaugeacts/aniso_sym_shared_functions.h"

namespace Chroma { 

  namespace AnisoSym { 


    void deriv_part(const int mu, 
		    const int nu, 
		    const int t_dir,
		    const Real& c_plaq_munu, 
		    const Real& c_rect_munu, 
		    const bool noTemporal2Link,
		    multi1d<LatticeColorMatrix>& ds_u,
		    const multi1d<LatticeColorMatrix>& u) 
    {

          START_CODE();


	  LatticeColorMatrix ds_plaq_mu;
	  LatticeColorMatrix ds_plaq_nu;
	  LatticeColorMatrix ds_rect_mu;
	  LatticeColorMatrix ds_rect_nu;
	

	  LatticeColorMatrix u_mu_xplusnu = shift(u[mu], FORWARD, nu);
	  LatticeColorMatrix u_nu_xplusmu = shift(u[nu], FORWARD, mu);
	  LatticeColorMatrix t_1, t_2, t_3;



	  //  --<--
	  // |               U^+(x+nu, mu) U^+(x,nu)
	  // V           
	  // |
	  // x
	  t_1 = adj(u_mu_xplusnu)*adj(u[nu]);

	  //  --<--
	  //       |
	  //       ^
	  //       |
	           
	  t_3 = u_nu_xplusmu*adj(u_mu_xplusnu);



	  //  -->-- 
	  //       |
          //       V
	  //       |
	  //  --<--
	  LatticeColorMatrix right_staple = adj(t_3)*adj(u[mu]);

	  if( nu > mu ) {
	    //  --<--
	    // |      |
	    // V      ^
	    // |      |
	    
	    LatticeColorMatrix up_staple = u_nu_xplusmu*t_1;

	    ds_plaq_mu = up_staple;
	    ds_plaq_nu = right_staple;
	  }
	  //  --<--
	  // |
	  // V
	  // | 
	  //  -->--          U^+(x+nu, mu) U^+(x,nu) U(x,mu)
	  // x
	  t_2 = t_1 * u[mu];

	  //  --<--
	  // |
	  // V              U^+(x-mu+nu, mu) U^+(x-mu,nu) U(x-mu,mu)
	  // |              
	  //  -->--
	  //      x
	  LatticeColorMatrix left_staple = shift(t_2, BACKWARD, mu);
	  if( nu > mu ) { 
	    ds_plaq_nu += left_staple;
	  }
	  
	  //  
	  //  |
	  //  ^
	  //  |
	  //  --<--
	  LatticeColorMatrix l_left_corner = adj(u[mu])*u[nu];
	  
	  if( nu > mu ) { 
	    //  
	    //  |     |
	    //  ^     V
	    //  |     |
	    //   --<--
	    t_2 = adj(u_nu_xplusmu)*l_left_corner;
	    t_1 = shift(t_2, BACKWARD, nu);
	    ds_plaq_mu += t_1;

	    Real c =  c_plaq_munu/Real(-2*Nc);
	    ds_u[mu]  += c*ds_plaq_mu;

	    ds_u[nu]  += c*ds_plaq_nu;
	  }
	  // ------------- Plaquette done here, now rectangle ----
	  
	  // If mu is temporal then we are computing contributions from 
	  // rectangles with time extent = 2,  both to ds_u[mu] and ds_u[nu].
	  // 
	  // If contribs from such rectangles are unrequired, then set 
	  // the appropriate coefficient to zero, and/or set noTemporal2Link
	  // or do both
	  bool skip = ( mu == t_dir && noTemporal2Link );

	  if( !skip) { 


	  //        --<--
	  //        |          
	  //        V            U^+(x-mu+nu, mu) U^+(x-mu,nu) U(x-mu,mu) U(x,mu)
	  //        |               
	  //        -->-- -->--
	  //             x
	  t_1 = left_staple * u[mu];
	
	  //        --<-- --<--
	  //        |            
	  //        V           U^+(x+nu, mu) U^+(x-mu+nu, mu) U^+(x-mu,nu) U(x-mu,mu) U(x,mu)
	  //        |
	  //        -->-- -->--
	  //             x

	  t_2= adj(u_mu_xplusnu)*t_1;

	  //        --<-- --<--
	  //        |             
	  //        V           U^+(x-mu+nu, mu) U^+(x-2mu+nu, mu) U^+(x-2mu,nu) U(x-2mu,mu) U(x-mu,mu)
	  //        |
	  //        -->-- -->--
	  //                    x


	  //=U^dag(x+-mu+nu, mu)*U^dag(x-2mu+nu, mu)*U^dag(x-2mu, nu)*U(x-2mu,mu)*U(x-mu,mu)	  

	  ds_rect_nu = shift(t_2, BACKWARD, mu);	  
	  



	  


	  
	  //        -->--  
	  //             |
	  //             V   U(x+mu+nu, mu) U^+(x+2mu, nu) U^+(x+mu,mu)
	  //             |
	  // x      --<--
	  LatticeColorMatrix shift_right_staple = shift(right_staple, FORWARD, mu);

	  

	  // -->--  -->--  
	  //             |
	  //             V   U(x+nu, mu)U(x+mu+nu, mu) U^+(x+2mu, nu) U^+(x+mu,mu)
	  //             |
	  // x      --<--
	  t_2 = u_mu_xplusnu*shift_right_staple;

	  // -->--  -->--  
	  //             |
	  //             V   U(x+nu, mu)U(x+mu+nu, mu) U^+(x+2mu, nu) U^+(x+mu,mu) U^(x,mu)
	  //             |
	  // --<--  --<--
	  // x
	  
	  // = u(x+nu, mu)*u(x+mu+nu, mu)*u^dag(x+2mu, nu)*u^dag(x+mu, mu)*u^dag(x,mu)
	  
	  ds_rect_nu += t_2*adj(u[mu]);
	  

	  
	  
	  


	    //     --<-- --<--
	    //    |           |
	    //    V           ^
	    //    |           |
	    //           -->--
	    //    x 
  

	    ds_rect_mu = adj(t_2)*adj(u[nu]);

	  
	  
	    //     --<-- --<--
	    //    |           |
	    //    V           ^
	    //    |           |
	    //     -->--     
	    //         x 
	    
	    ds_rect_mu += t_3*left_staple;

	    //
	    //           -->--
	    //    |           |
	    //    ^           V
	    //    |           |
	    //    --<--  --<--
	    //    x 
	    
	    t_2 = shift_right_staple*l_left_corner;


	    //    
	    //         |
	    //         V
	    //         |
	    //    --<-- 
	    //    x 
	    
	    t_1 = adj(u_nu_xplusmu)*adj(u[mu]);

	 
	    //    -->--
	    //    |          |
	    //    ^          V
	    //    |          |
	    //    --<--  --<-- 
	    //         x 
	    t_3 = t_1*adj(left_staple);
	    t_2 += t_3;

	    ds_rect_mu += shift(t_2, BACKWARD, nu);
	    Real c = c_rect_munu/Real(-2*Nc);

	    ds_u[mu] += c * ds_rect_mu;
	    ds_u[nu] += c * ds_rect_nu;
	  }



	  END_CODE();

	  return;

    }


    // S_part
    void  S_part(int mu, int nu, int t_dir,
		 Real c_plaq_munu,
		 Real c_rect_munu, 
		 bool noTemporal2Link,
		 LatticeReal& lgimp,
		 const multi1d<LatticeColorMatrix>& u)
    {
      START_CODE();
      
      LatticeColorMatrix tmp1;
      LatticeColorMatrix tmp2;
      
      LatticeColorMatrix rectangle_2munu;
      LatticeColorMatrix rectangle_mu2nu;
      
      LatticeColorMatrix u_nu_xplus_mu = shift(u[nu], FORWARD, mu);
      LatticeColorMatrix u_mu_xplus_nu = shift(u[mu], FORWARD, nu);
      
      //                     ^
      // lr_corner =         |   = u(x, mu) * u(x + mu, nu)
      //                     |
      //                ----->
      //
      LatticeColorMatrix lr_corner = u[mu]*u_nu_xplus_mu;
      
      //              <----
      //                   ^
      // right_staple =    | = u(x, mu) * u(x + mu, nu) * u^{+}(x+nu,mu)
      //                   |
      //              ----->
      //
      
      LatticeColorMatrix right_staple = lr_corner*adj(u_mu_xplus_nu);
      
      {
	//              <----
	//            |      ^
	// plaq     = |      | = u(x, mu) * u(x + mu, nu) * u^{+}(x+nu,mu)
	//            V      |    *u^{+}(x,nu)
	//              ----->
	
	if( nu > mu) {
	  tmp1 = right_staple*adj(u[nu]);
	  lgimp += c_plaq_munu*real(trace(tmp1));
	}
	
      }
      
      // End of Plaquette bit...
      bool skip = ( noTemporal2Link && (t_dir == mu) );
      if( ! skip ) {
	//            <----  <----^
	// Make:                  |
	//            -----> ----->
	
	//           <----
	//                |
	//      ----> ---->
	
	
	tmp1 = shift(right_staple, FORWARD, mu);
	tmp2 = u[mu]*tmp1;
	
	
	//      <---- <----
	//                |
	//      ----> ---->
	right_staple = tmp2*adj(u_mu_xplus_nu);
	
	//   Loop 
	//        =  u(x, mu) u(x + mu, mu) * u(x + 2mu, nu)
	//         * u^dag(x + mu + nu, mu) * u^dag(x + nu, mu) * u^dag(nu)
	rectangle_2munu = right_staple*adj(u[nu]);
	
	
	lgimp += c_rect_munu * real(trace(rectangle_2munu));
      }

      END_CODE();
    }


  }

}
