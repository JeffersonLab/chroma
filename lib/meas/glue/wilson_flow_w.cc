// $Id: wilson_flow_w.cc,v 1.12 2011/12/11 17:23:22 cmcneile Exp cmcneile $
/*! \file
 *  \brief Code for Wilson flow
 *   
 *   A collection of routines to compute the
 *   Wilson flow.
 *
 *   Essentially the code implements appendix C
 *   of
 *
 *    Properties and uses of the Wilson flow in lattice QCD.
 *    Martin Luscher 
 *    Published in JHEP 1008 (2010) 071 , arXiv:1006.4518 
 *
 *    The code calls the existing stout smearing routines
 *    -- which already exist in chroma.
 *
 *   The Wilson flow can be used to determine the lattice spacing.
 *
 *
 * See below for additional information about the Wilson flow.
 * 
 * Continuous smearing of Wilson Loops.
 * Robert Lohmayer, Herbert Neuberger. arXiv:1110.3522. 
 * Published in PoS LATTICE2011 (2011) 249 
 *
 *
 *
 *
 */

#include "meas/glue/wilson_flow_w.h"
#include "meas/glue/mesfield.h"
#include "util/gauge/stout_utils.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"

//using namespace Chroma;
namespace Chroma
{



  /**

   **/


  void measure_wilson_gauge(multi1d<LatticeColorMatrix> & u,
			    Real & gspace, Real & gtime,
			    int jomit)
  {

    multi1d<LatticeColorMatrix> field_st(10) ;
    LatticeColorMatrix tmp ; 

    mesField( field_st,u); 

    int offset = 0;

    gtime  = 0.0 ;
    gspace = 0.0 ;

    Double tr ; 

    //  for(int mu=0; mu < Nd-1; ++mu)

    for(int mu=0; mu < Nd; ++mu)
    {
      for(int nu=mu+1; nu < Nd; ++nu)
      {
	//	  tr = real(sum(trace(field_st[offset]))) ;
	//	  tr = imag(sum(trace(field_st[offset]))) ;
	tmp = field_st[offset] * field_st[offset] ;
	tr = real(sum(trace(tmp))) ; 

	// Real tt = 2.0 * tr ; 
	//	  cout << "DEBUG " << mu << " " << nu  << " " << tt << endl ;

	if (nu==jomit)
	{
	  gtime += 2.0*(tr);
	}
	else
	{
	  gspace += 2.0*(tr);
	}

	++offset  ; 
      }
    }

    gspace /= -2.0*Layout::vol() ;
    gtime /= -2.0*Layout::vol() ;

  }


  void wilson_flow_one_step(multi1d<LatticeColorMatrix> & u, Real rho)
  {
    int mu, dir;
    multi1d<LatticeColorMatrix> dest(Nd);
    multi1d<LatticeColorMatrix> next(Nd);


    // -------------------------------------

    multi1d<bool> smear_in_this_dirP(4) ;
    multi2d<Real> rho_a(4,4) ;
    multi2d<Real> rho_b1(4,4) ;
    multi2d<Real> rho_b2(4,4) ;
    multi2d<Real> rho_c(4,4) ;

    for (mu = 0; mu <= Nd-1; mu++)
    {
      smear_in_this_dirP(mu) = true ;
      for (dir = 0; dir <= Nd-1; dir++)
      {
	rho_a[mu][dir] = rho * 0.25 ;

	rho_b1[mu][dir] = rho * 8.0/9.0 ;
	rho_b2[mu][dir] = rho * 17.0/36.0 ;

	rho_c[mu][dir] = rho * 3.0/4.0 ;

      }
    }


    Stouting::smear_links(u, dest,smear_in_this_dirP, rho_a);

    LatticeColorMatrix  Q, QQ, C ;
    LatticeColorMatrix  Q2, QQ2  ;

    multi1d<LatticeColorMatrix> Q0(Nd);
    multi1d<LatticeColorMatrix> Q1(Nd);


    multi1d<LatticeComplex> f;   // routine will resize these


    for (mu = 0; mu <= Nd-1; mu++)
    {
      Stouting::getQsandCs(dest, Q1[mu],QQ,C,mu,smear_in_this_dirP,rho_b1) ;
      Stouting::getQsandCs(u   , Q0[mu],QQ,C,mu,smear_in_this_dirP,rho_b2) ;

      Q = Q1[mu] - Q0[mu] ;
      QQ = Q * Q ;
      Stouting::getFs(Q,QQ,f);   // This routine computes the f-s
          
      // Assemble the stout links exp(iQ)U_{mu} 
      next[mu]=(f[0] + f[1]*Q + f[2]*QQ)*dest[mu];      

    }

    for (mu = 0; mu <= Nd-1; mu++)
    {
      u[mu]    =  next[mu] ;
      dest[mu] =  next[mu] ;      
    }


    for (mu = 0; mu <= Nd-1; mu++)
    {
      Stouting::getQsandCs(dest, Q2,QQ,C,mu,smear_in_this_dirP,rho_c) ;

      Q = Q2 - Q1[mu] + Q0[mu] ;
      QQ = Q * Q ;
      Stouting::getFs(Q,QQ,f);   // This routine computes the f-s
          
      // Assemble the stout links exp(iQ)U_{mu} 
      next[mu]=(f[0] + f[1]*Q + f[2]*QQ)*dest[mu];      

    }


    for (mu = 0; mu <= Nd-1; mu++)
    {
      u[mu]    =  next[mu] ;
    }



  }


  void wilson_flow(XMLWriter& xml,
		   multi1d<LatticeColorMatrix> & u, int nstep, 
		   Real  wflow_eps, int jomit)
  {
    Real gact4i, gactij;
    int dim = nstep + 1 ;
    multi1d<Real> gact4i_vec(dim);
    multi1d<Real> gactij_vec(dim);
    multi1d<Real> step_vec(dim);



    measure_wilson_gauge(u,gactij,gact4i,jomit) ;
    gact4i_vec[0] = gact4i ;
    gactij_vec[0] = gactij ;
    step_vec[0] = 0.0 ;

    //  QDPIO::cout << "WFLOW " << 0.0 << " " << gact4i << " " << gactij <<  endl ; 

    QDPIO::cout << "START_ANALYZE_wflow" << endl ; 
    QDPIO::cout << "WFLOW time gact4i gactij" << endl ; 

    for(int i=0 ; i < nstep ; ++i)
    {
      wilson_flow_one_step(u,wflow_eps) ;

      measure_wilson_gauge(u,gactij,gact4i,jomit) ;
      gact4i_vec[i+1] = gact4i ;
      gactij_vec[i+1] = gactij ;


      Real xx = (i + 1) * wflow_eps ;
      QDPIO::cout << "WFLOW " << xx << " " << gact4i << " " << gactij <<  endl ; 

      step_vec[i+1] = xx ;

    }
    QDPIO::cout << "END_ANALYZE_wflow" << endl ; 

    push(xml, "wilson_flow_results");
    write(xml,"wflow_step",step_vec) ; 
    write(xml,"wflow_gact4i",gact4i_vec) ; 
    write(xml,"wflow_gactij",gactij_vec) ; 
    pop(xml);  // elem

  }


}  // end namespace Chroma


