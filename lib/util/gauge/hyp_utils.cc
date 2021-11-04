// -*- C++ -*-
/*! \file
 *  \brief Stout utilities
 */

#include "chroma_config.h"
#include "chromabase.h"
#include "util/gauge/hyp_utils.h"
#if QDP_NC != 3
#include "util/gauge/expm12.h"
#include "util/gauge/sun_proj.h"
#endif


namespace Chroma 
{ 

  //! Timings
  /*! \ingroup gauge */
  namespace HypLinkTimings { 
    static double smearing_secs = 0;
    double getSmearingTime() { 
      return smearing_secs;
    }

    static double force_secs = 0;
    double getForceTime() { 
      return force_secs;
    }

    static double functions_secs = 0;
    double getFunctionsTime() { 
      return functions_secs;
    }
  }

  //! Utilities
  /*! \ingroup gauge */
  namespace Hyping 
  {

    //! Get the level 1 staples

    //! Do the force recursion from level i+1, to level i
    void deriv_recurse(multi1d<LatticeColorMatrix>&  F,
		       const multi1d<bool>& smear_in_this_dirP,
                       const Real alpha1,
                       const Real alpha2,
                       const Real alpha3,
                       const int BlkMax,
                       const Real BlkAccu,
                       const multi1d<LatticeColorMatrix>& u)
    {
      START_CODE();

      END_CODE();
    }

    /*! \ingroup gauge */
    void smear_links(const multi1d<LatticeColorMatrix>& u, 
                     multi1d<LatticeColorMatrix>& u_hyp,
                     const multi1d<bool>& smear_in_this_dirP,
                     const Real alpha1,
                     const Real alpha2,
                     const Real alpha3,
                     const int BlkMax,
                     const Real BlkAccu)
    {
      multi1d<LatticeColorMatrix> u_lv1(Nd*(Nd-1));
      multi1d<LatticeColorMatrix> u_lv2(Nd*(Nd-1));
      LatticeColorMatrix u_tmp;
      LatticeColorMatrix tmp_1;
      Real ftmp1;
      Real ftmp2;
      int rho;
      int sigma;
      int ii;
      int jj;
      int kk;

      START_CODE();
  
      if (Nd > 4)
        QDP_error_exit("Hyp-smearing only implemented for Nd<=4",Nd);
      
      //QDPIO::cout << "HYP start " << std::endl;      
      /*
       * Construct "level 1" smeared links in mu-direction with
       * staples only in one orthogonal direction, nu
       */
      ftmp1 = 1.0 - alpha3;
      ftmp2 = alpha3 / 2;
      ii = -1;
      for(int mu = 0; mu < Nd; ++mu) {
        //QDPIO::cout << "HYP start mu " << mu << std::endl;      
        for(int nu = 0; nu < Nd; ++nu) {
          if(nu == mu) continue;
          //QDPIO::cout << "HYP start nu " << nu << std::endl;      
          ii++;
          /*
           * Forward staple
           *
           * u_tmp(x) = u(x,nu)*u(x+nu,mu)*u_dag(x+mu,nu)
           */
          u_tmp = u[nu] * shift(u[mu],FORWARD,nu) * adj(shift(u[nu],FORWARD,mu));
          //QDPIO::cout << "HYP 1 " << std::endl;      
          /*
           * Backward staple
           *
           * u_tmp(x) += u_dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)
           */
          LatticeColorMatrix u_nu_tmp = shift(u[nu],FORWARD,mu);
          u_tmp += shift(adj(u[nu]) * u[mu] * u_nu_tmp,BACKWARD,nu);
          //QDPIO::cout << "HYP 2 " << std::endl;                
          /*
           * Unprojected level 1 link
           */
          tmp_1 = ftmp1*u[mu] + ftmp2*u_tmp;
          u_tmp = adj(tmp_1);
          //QDPIO::cout << "HYP 3 " << std::endl;      
          /*
           * Project onto SU(Nc)
           */
          u_lv1[ii] = u[mu];
          //QDPIO::cout << "HYP 4 " << std::endl;      
          sun_proj(u_tmp, u_lv1[ii], BlkAccu, BlkMax);
          //QDPIO::cout << "HYP sun " << std::endl;      
        }
      }
      //QDPIO::cout << "HYP start lvl1 done" << std::endl;      

      /*
       * Construct "level 2" smeared links in mu-direction with
       * "level 1" staples not in the orthogonal direction, nu,
       * and the "level 1" links decorated in the 4-th orthogonal direction
       */
      ftmp1 = 1.0 - alpha2;
      ftmp2 = alpha2 / 4;
      ii = -1;
      for(int mu = 0; mu < Nd; ++mu) {
        for(int nu = 0; nu < Nd; ++nu) {
          if(nu == mu) continue;
          
          ii++;
          u_tmp = 0;
          for(rho = 0; rho < Nd; ++rho) {
            if(rho == mu || rho == nu) continue;
            
            /* 4-th orthogonal direction: sigma */
            for(jj = 0; jj < Nd; ++jj) {
              if(jj != mu && jj != nu && jj != rho) sigma = jj;
            }
            jj = (Nd-1)*mu + sigma;
            if(sigma > mu ) jj--;
            kk = (Nd-1)*rho + sigma;
            if(sigma > rho ) kk--;
            
            /*
             * Forward staple
             *
             * u_tmp(x) += u_lv1(x,kk)*u_lv1(x+rho,jj)*u_lv1_dag(x+mu,kk)
             */
            u_tmp += u_lv1[kk] * shift(u_lv1[jj],FORWARD,rho) * adj(shift(u_lv1[kk],FORWARD,mu));
            
            /*
             * Backward staple
             *
             * u_tmp(x) += u_lv1_dag(x-rho,kk)*u_lv1(x-rho,jj)*u_lv1(x-rho+mu,kk)
             */
            LatticeColorMatrix u_lv1_tmp = shift(u_lv1[kk],FORWARD,mu);
            u_tmp += shift(adj(u_lv1[kk]) * u_lv1[jj] * u_lv1_tmp , BACKWARD,rho);
          }
          
          /*
           * Unprojected level 2 link
           */
          tmp_1 = ftmp1*u[mu] + ftmp2*u_tmp;
          u_tmp = adj(tmp_1);
          
          /*
           * Project onto SU(Nc)
           */
          u_lv2[ii] = u[mu];
          sun_proj(u_tmp, u_lv2[ii], BlkAccu, BlkMax);
          //QDPIO::cout << "HYP start lvl1 sun done" << std::endl;      
        }
      }
      //QDPIO::cout << "HYP start lvl2 done" << std::endl;      
      /*
       * Construct hyp-smeared links in mu-direction with
       * "level 2" staples in the orthogonal direction, nu,
       * and the "level 2" links not decorated in the mu and nu directions
       */
      ftmp1 = 1.0 - alpha1;
      ftmp2 = alpha1 / 6;
      for(int mu = 0; mu < Nd; ++mu) {
        u_tmp = 0;
        for(int nu = 0; nu < Nd; ++nu) {
          if(nu == mu) continue;
          
          jj = (Nd-1)*mu + nu;
          if(nu > mu ) jj--;
          kk = (Nd-1)*nu + mu;
          if(mu > nu ) kk--;
          
          /*
           * Forward staple
           *
           * u_tmp(x) += u_lv2(x,kk)*u_lv2(x+nu,jj)*u_lv2_dag(x+mu,kk)
           */
          u_tmp += u_lv2[kk] * shift(u_lv2[jj],FORWARD,nu) * adj(shift(u_lv2[kk],FORWARD,mu));
          
          /*
           * Backward staple
           *
           * u_tmp(x) += u_lv2_dag(x-nu,kk)*u_lv2(x-nu,jj)*u_lv2(x-nu+mu,kk)
           */
          LatticeColorMatrix u_lv2_tmp = shift(u_lv2[kk],FORWARD,mu);
          u_tmp += shift(adj(u_lv2[kk]) * u_lv2[jj] * u_lv2_tmp , BACKWARD,nu);
        }
        
        /*
         * Unprojected hyp-smeared link
         */
        tmp_1 = ftmp1*u[mu] + ftmp2*u_tmp;
        u_tmp = adj(tmp_1);
        
        /*
         * Project onto SU(Nc)
         */
        u_hyp[mu] = u[mu];
        sun_proj(u_tmp, u_hyp[mu], BlkAccu, BlkMax);
      }
      //QDPIO::cout << "HYP start lvl3 done" << std::endl;      
      END_CODE();
    }
  } // End Namespace Hyping
} // End Namespace Chroma
	
    

