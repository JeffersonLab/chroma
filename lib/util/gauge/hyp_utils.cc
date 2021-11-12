// -*- C++ -*-
/*! \file
 *  \brief Hyp utilities
 */

#include "chroma_config.h"
#include "chromabase.h"
#include "util/gauge/hyp_utils.h"
#include "util/gauge/sun_proj.h"
#if QDP_NC != 3
#include "util/gauge/expm12.h"
#endif

#if defined(BUILD_JIT_CLOVER_TERM)
void function_hyp_smear_link_exec(JitFunction& function, 
                                  const LatticeColorMatrix& u, 
                                  LatticeColorMatrix& u_hyp,
                                  const Real alpha1,
                                  const Real alpha2,
                                  const Real alpha3,
                                  const int BlkMax,
                                  const Real BlkAccu);

void function_hyp_smear_link_build(JitFunction& function,
                                   const LatticeColorMatrix& u, 
                                   LatticeColorMatrix& u_hyp,
                                   const Real alpha1,
                                   const Real alpha2,
                                   const Real alpha3,
                                   const int BlkMax,
                                   const Real BlkAccu);

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
                       const int hyp_qr_max_iter,
                       const Real hyp_qr_tol,
                       const multi1d<LatticeColorMatrix>& U)
    {
      START_CODE();
      
      multi1d<LatticeColorMatrix> UH(Nd);
      multi1d<LatticeColorMatrix> P(Nd);
      
      for(int mu=0; mu<Nd; mu++) {        
        if(smear_in_this_dirP[mu]) {
          // Initialise P to identity
          P[mu] = 1;
          upper_hessenberg_link(U[mu], UH[mu], P[mu]);          
        }
      }

      END_CODE();
    }

    /* A namespace to hide the thread dispatcher in */
    namespace HypUtils { 
      
      struct HypUpperHessLinkArgs { 

        // Site link
        const LatticeColorMatrix &U;

        // Accmulated reflectors
        LatticeColorMatrix &P;

        // The Upper Hessenberg matrix
        LatticeColorMatrix &UH;
      };
      
      inline void hypUpperHessSiteLoop(int lo, int hi, int myId, 
                                       HypUpperHessLinkArgs* arg)
      {
#ifndef QDP_IS_QDPJIT
#if 0
        // We follow the computation described in arxiv:1606.01277
        //--------------------------------------------------------------
        // 1. Compute V = Omega * Q^(1/2) | Q = (Omega^dag * Omega)
        //    To compute the square root for all Nc we reduce the SU(Nc)
        //    matrix to UpperHessenberg form, then perform a QR 
        //    decomposition to further reduce it to UpperTriangular.
        //    The eigenvalues will then lie on the diagonal.
        
        const LatticeColorMatrix &U = arg->U;   // The original links
        LatticeColorMatrix &UH = arg->UH; // Upper Hessenberg reduction
        LatticeColorMatrix &P = arg->P;   // Accumulated reflectors
      
        Real tol = 1e-12;

        for(int site=lo; site < hi; site++) { 
      
          // Reflector at iteration `site`
          PColorMatrix<QDP::RComplex<REAL>, Nc> Pi;

          // Site local Upper Hessenberg reduction, intitialised to U.
          PColorMatrix<QDP::RComplex<REAL>, Nc> UH_site = arg->U.elem(site).elem();
          // Get accumulated reflecor.
          PColorMatrix<QDP::RComplex<REAL>, Nc> P_site = arg->P.elem(site).elem();
          
          for(int i = 0; i < Nc - 2; i++) {
            
            REAL col_norm, col_norm_inv;
            // get (partial) dot product of ith column vector
            col_norm = 0.0;
            for(int j = i+1; j < Nc; j++) col_norm += fabs(UH_site.elem().elem(i,j));
            col_norm = sqrt(col_norm);

            PColorVector<QDP::RComplex<REAL>, Nc> v;
            LatticeColorVector v; // vector of length (1 x N)
            QDP::RComplex<REAL> rho = (REAL)1.0;
            REAL abs_elem = fabs(UH_site.elem().elem(i+1,i));
            if(abs_elem > tol) {
              for(int k=0; k<Nc; k++) rho = -UH_site.elem().elem(i+1,k)/abs_elem;
            }
            v.elem().elem(i+1) = UH_site(i+1,i) - col_norm * rho;

            // reuse col_norm
            col_norm = fabs(v.elem().elem(i+1));

            // copy the rest of the column
            for(int j = i + 2; j < Nc; j++ ) {
              v.elem().elem(j) = UH_site.elem().elem(j,i);
              col_norm += fabs(v.elem().elem(j));
            }
            col_norm_inv = 1.0/std::max(sqrt(col_norm), 1e-30);

            // Normalise the column
            for(int j = i + 1; j < Nc; j++) v.elem().elem(j) *= col_norm_inv;
            
            // Construct the householder matrix P = I - 2 v * v^{\dag}
            Pi = 1;
            for(int j = i + 1; j < Nc; j++ ) {
              for(int k = i + 1; k < Nc; k++ ) {
                QDP::RComplex<REAL> P_elem;
                P_elem.elem().real()  = v.elem().elem(j).real() * v.elem().elem(k).real();
                P_elem.elem().real() += v.elem().elem.imag() * v.elem().elem.imag();
                P_elem.elem().imag()  = v.elem().elem.imag() * v.elem().elem.real();
                P_elem.elem().imag() -= v.elem().elem.real() * v.elem().elem.imag();
                Pi.elem().elem(j,k) -= 2.0 * P_elem;
              }
            }
            
            // Similarity transform
            UH = Pi * UH * Pi;
            // Store reflector
            P_site = P_site * Pi;
            
            // Load back into the lattice sized objects
            for(int j=0; j < Nc; j++) { 
              for(int k=0; k < Nc; k++) { 
                P.elem(site).elem().elem(i,k).real() = P_site.elem().elem(i,k).real();
                P.elem(site).elem().elem(i,k).imag() = P_site.elem().elem(i,k).imag();

                UH.elem(site).elem().elem(i,k).real() = UH_site.elem().elem(i,k).real();
                UH.elem(site).elem().elem(i,k).imag() = UH_site.elem().elem(i,k).imag();
                
              }
            }
          }
        }
#endif
        
#endif
      } // End Function
      
    }// End Namespace
        
    void upper_hessenberg_link(const LatticeColorMatrix &U, 
                               LatticeColorMatrix &UH,
                               LatticeColorMatrix &P)
      
    {
      START_CODE();
      QDP::StopWatch swatch;
      swatch.reset();
      swatch.start();
      
      int num_sites = Layout::sitesOnNode();
      HypUtils::HypUpperHessLinkArgs args={U, UH, P};
      
      
#if !defined(BUILD_JIT_CLOVER_TERM)

#ifdef QDP_IS_QDPJIT
#warning "Hyping disabled, as building with QDP-JIT but the JIT Clover term is disabled."
#else
#warning "Using QDP++ hyping"
      dispatch_to_threads(num_sites, args, HypUtils::hypUpperHessSiteLoop);
#endif
      
#else
#warning "Using QDP-JIT hyping"
      static JitFunction function;
      if (function.empty()) {
        //function_upper_hess_link_build(function, U, UH, P);
        //function_upper_hess_link_exec(function, U, UH, P);
      }
#endif

      swatch.stop();
      HypLinkTimings::functions_secs += swatch.getTimeInSeconds();
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
#if 0
      hyp_smear_link(u, u_hyp, alpha1, alpha2, alpha3, BlkMax, BlkAccu);
#else
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
      
      // Construct "level 1" smeared links in mu-direction with
      //staples only in one orthogonal direction, nu
      ftmp1 = 1.0 - alpha3;
      ftmp2 = alpha3 / 2;
      ii = -1;
      for(int mu = 0; mu < Nd; ++mu) {
        if( smear_in_this_dirP[mu] ) { 
          //QDPIO::cout << "HYP start mu " << mu << std::endl;      
          for(int nu = 0; nu < Nd; ++nu) {
            if(nu == mu) continue;
            //QDPIO::cout << "HYP start nu " << nu << std::endl;      
            ii++;
          
            // Forward staple
            // u_tmp(x) = u(x,nu)*u(x+nu,mu)*u_dag(x+mu,nu)
            u_tmp = u[nu] * shift(u[mu],FORWARD,nu) * adj(shift(u[nu],FORWARD,mu));
            
            //QDPIO::cout << "HYP 1 " << std::endl;      
            // Backward staple
            // u_tmp(x) += u_dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)
            LatticeColorMatrix u_nu_tmp = shift(u[nu],FORWARD,mu);
            u_tmp += shift(adj(u[nu]) * u[mu] * u_nu_tmp,BACKWARD,nu);
            //QDPIO::cout << "HYP 2 " << std::endl;                

            // Unprojected level 1 link
            tmp_1 = ftmp1*u[mu] + ftmp2*u_tmp;
            u_tmp = adj(tmp_1);
            //QDPIO::cout << "HYP 3 " << std::endl;      
          
            // Project onto SU(Nc)
            u_lv1[ii] = u[mu];
            //QDPIO::cout << "HYP 4 " << std::endl;      
            sun_proj(u_tmp, u_lv1[ii], BlkAccu, BlkMax);
            //QDPIO::cout << "HYP sun " << std::endl;      
          }
        }
      }
      //QDPIO::cout << "HYP start lvl1 done" << std::endl;      
      
      // Construct "level 2" smeared links in mu-direction with
      // "level 1" staples not in the orthogonal direction, nu,
      // and the "level 1" links decorated in the 4-th orthogonal direction
      ftmp1 = 1.0 - alpha2;
      ftmp2 = alpha2 / 4;
      ii = -1;
      for(int mu = 0; mu < Nd; ++mu) {
        if( smear_in_this_dirP[mu] ) { 
          for(int nu = 0; nu < Nd; ++nu) {
            if(nu == mu) continue;
            
            ii++;
            u_tmp = 0;
            for(rho = 0; rho < Nd; ++rho) {
              if(rho == mu || rho == nu) continue;
              
              // 4-th orthogonal direction: sigma 
              for(jj = 0; jj < Nd; ++jj) {
                if(jj != mu && jj != nu && jj != rho) sigma = jj;
              }
              jj = (Nd-1)*mu + sigma;
              if(sigma > mu ) jj--;
              kk = (Nd-1)*rho + sigma;
              if(sigma > rho ) kk--;
             
              // Forward staple
              // u_tmp(x) += u_lv1(x,kk)*u_lv1(x+rho,jj)*u_lv1_dag(x+mu,kk)
              u_tmp += u_lv1[kk] * shift(u_lv1[jj],FORWARD,rho) * adj(shift(u_lv1[kk],FORWARD,mu));
            
              // Backward staple
              // u_tmp(x) += u_lv1_dag(x-rho,kk)*u_lv1(x-rho,jj)*u_lv1(x-rho+mu,kk)
              LatticeColorMatrix u_lv1_tmp = shift(u_lv1[kk],FORWARD,mu);
              u_tmp += shift(adj(u_lv1[kk]) * u_lv1[jj] * u_lv1_tmp , BACKWARD,rho);
            }
          
            // Unprojected level 2 link
            tmp_1 = ftmp1*u[mu] + ftmp2*u_tmp;
            u_tmp = adj(tmp_1);
          
            // Project onto SU(Nc)
            u_lv2[ii] = u[mu];
            sun_proj(u_tmp, u_lv2[ii], BlkAccu, BlkMax);
            //QDPIO::cout << "HYP start lvl1 sun done" << std::endl;      
          }
        }
      }
      //QDPIO::cout << "HYP start lvl2 done" << std::endl;      
      
      // Construct hyp-smeared links in mu-direction with
      // "level 2" staples in the orthogonal direction, nu,
      // and the "level 2" links not decorated in the mu and nu directions
      ftmp1 = 1.0 - alpha1;
      ftmp2 = alpha1 / 6;
      for(int mu = 0; mu < Nd; ++mu) {
        if( smear_in_this_dirP[mu] ) { 
          u_tmp = 0;
          for(int nu = 0; nu < Nd; ++nu) {
            if(nu == mu) continue;
            
            jj = (Nd-1)*mu + nu;
            if(nu > mu ) jj--;
            kk = (Nd-1)*nu + mu;
            if(mu > nu ) kk--;
            
            
            // Forward staple
            // u_tmp(x) += u_lv2(x,kk)*u_lv2(x+nu,jj)*u_lv2_dag(x+mu,kk)
            u_tmp += u_lv2[kk] * shift(u_lv2[jj],FORWARD,nu) * adj(shift(u_lv2[kk],FORWARD,mu));
          
            // Backward staple
            // u_tmp(x) += u_lv2_dag(x-nu,kk)*u_lv2(x-nu,jj)*u_lv2(x-nu+mu,kk)
            LatticeColorMatrix u_lv2_tmp = shift(u_lv2[kk],FORWARD,mu);
            u_tmp += shift(adj(u_lv2[kk]) * u_lv2[jj] * u_lv2_tmp , BACKWARD,nu);
          }
          
          // Unprojected hyp-smeared link
          tmp_1 = ftmp1*u[mu] + ftmp2*u_tmp;
          u_tmp = adj(tmp_1);
        
          // Project onto SU(Nc)
          u_hyp[mu] = u[mu];
          sun_proj(u_tmp, u_hyp[mu], BlkAccu, BlkMax);
        }
      }
      //QDPIO::cout << "HYP start lvl3 done" << std::endl;      
      END_CODE();
#endif
    }
  } // End Namespace Hyping
} // End Namespace Chroma

    

