// $Id: dwf_quarkprop4_w.cc,v 1.6 2004-01-29 16:18:52 kostas Exp $
// $Log: dwf_quarkprop4_w.cc,v $
// Revision 1.6  2004-01-29 16:18:52  kostas
// Removed the if 1 directive
//
// Revision 1.5  2004/01/29 16:17:50  kostas
// Fixed an ax-ps correlation function bug
//
// Revision 1.4  2004/01/29 14:59:46  kostas
// Fixed the bugs. Code compiles now
//
/*! \file
 * \brief Full quark propagator solver for domain wall fermions
 *
 * Given a complete propagator as a source, this does all the inversions needed
 */


#include "chromabase.h"
#include "fermact.h"
#include "actions/ferm/qprop/dwf_quarkprop4_w.h"
#include "actions/ferm/linop/dwffld_w.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"

using namespace QDP;

//! Given a complete propagator as a source, this does all the inversions needed
/*! \ingroup qprop
 *
 * This routine is actually generic to Domain Wall fermions (Array) fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */
template<typename T, template<class> class C>
static 
void dwf_quarkProp4_a(LatticePropagator& q_sol, 
		      XMLWriter& xml_out,
		      const LatticePropagator& q_src,
		      const C<T>& S_f,
		      Handle<const ConnectState> state,
		      enum InvType invType,
		      const Real& RsdCG, 
		      int MaxCG, int& ncg_had)
{
  START_CODE("dwf_quarkProp4");

  push(xml_out, "DWF_QuarkProp4");

  ncg_had = 0;

  multi1d<LatticePropagator> prop5d(S_f.size()) ;
  LatticePropagator q_mp  ;

  multi1d<LatticeFermion> psi(S_f.size()) ;
 
  psi = zero;  // note this is ``zero'' and not 0
   
  // This version loops over all color and spin indices
  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    for(int spin_source = 0; spin_source < Ns; ++spin_source)
    {
      multi1d<LatticeFermion> tmp(S_f.size()) ;
      tmp = zero ;
      PropToFerm(q_src, tmp[0], color_source, spin_source);

      // Use the last initial guess as the current initial guess
	   
      /* 
       * Normalize the source in case it is really huge or small - 
       * a trick to avoid overflows or underflows
       */

      Real fact = Real(1) / sqrt(norm2(tmp[0]));
      for(int i = 0; i < tmp.size(); ++i)
	tmp[i] *= fact; 
	   
      multi1d<LatticeFermion> chi(S_f.size()) ;
	   
      // Split the source to oposite walls according to chirality
      DwfFld(chi, tmp, PLUS);

      // now we are ready invert
	   
      // Compute the propagator for given source color/spin.
      int n_count;
	   
      S_f.qpropT(psi, state, chi, invType, RsdCG, MaxCG, n_count);
      ncg_had += n_count;
	   
      push(xml_out,"Qprop");
      write(xml_out, "RsdCG", RsdCG);
      Write(xml_out, n_count);
      pop(xml_out);

      // Unnormalize the source following the inverse 
      // of the normalization above
      fact = Real(1) / fact;
      for(int i = 0; i < tmp.size(); ++i)
	psi[i] *= fact; 

      /*
       * Move the solution to the appropriate components
       * of quark propagator.
       */
	   
      //First the 5D quark propagator
      for(int s(0);s<S_f.size();s++)
	FermToProp(psi[s], prop5d[s], color_source, spin_source);
      // Now get the 4D propagator too

      // move the light states to wall 0
      DwfFld(tmp, psi, MINUS); 
      // move solution to the appropriate components of the 4d
      // quark propagator
      FermToProp(tmp[0], q_sol, color_source, spin_source);

      // move solution to the appropriate components of the 4d
      // midpoint quark propagator which occures at N5/2 -1 
      FermToProp(tmp[S_f.size()/2 - 1], q_mp, color_source, spin_source);

	   
    }	/* end loop over spin_source */
  } /* end loop over color_source */

  // constuct the conserved axial current correlator
  LatticeComplex cfield ;
  dwf_conserved_axial_ps_corr(cfield,state->getLinks(),prop5d,3,S_f.size());
			       
  multi2d<DComplex> corr ;  
   
  SftMom trick(0,false,3) ;
   
  corr = trick.sft(cfield);
  // Length of lattice in decay direction
  int length = trick.numSubsets() ;
  multi1d<Real> mesprop(length);
  
  for(int t(0);t<length; t++)
    mesprop[t] = real(corr[0][t]) ; // only need the zero momentum

  //push(xml_out, "t_dir");
  //Write(xml_out, 3);
  //pop(xml_out);
  push(xml_out, "DWF_ConservedAxial");
  //write(xml_out, "corr", corr);
  Write(xml_out, mesprop); 
  pop(xml_out);

  //Now the midpoint Pseudoscalar correlator
  cfield = trace(q_mp * adj(q_mp));    

  corr = trick.sft(cfield);
  for(int t(0);t<length; t++)
    mesprop[t] = real(corr[0][t]) ; // only need the zero momentum

  push(xml_out, "DWF_MidPoint_Pseudo");
  Write(xml_out, mesprop);
  pop(xml_out);

  pop(xml_out);   // DWF_QuarkProp4

  END_CODE("dwf_quarkProp4");
}


/*!
  Corelation function:
\f[
  C(t) = \sum_x \sum_s \left[\bar{\Psi}(x+\hat{\mu},t,s)\frac{1+\gamma_{\mu}}{2} 
                                  \Psi (x,t,s) - 
                             \bar{\Psi}(x,t,s)\frac{1-\gamma_{\mu}}{2} 
                                  \Psi (x+\hat{\mu},t,s)\right] 
                       \bar{q}(0,0)\gamma_5 q(0,0) sign(s - \frac{L_s - 1}{2})
\f]
**/
void dwf_conserved_axial_ps_corr(LatticeComplex& corr,
				 const multi1d<LatticeColorMatrix>& u,
				 const multi1d<LatticePropagator>& p5d, 
				 const int mu,
				 const int N5)
{
  // gamma mapping G_0 --> Gamma(1)
  //               G_1 --> Gamma(2)
  //               G_2 --> Gamma(4)
  //               G_3 --> Gamma(8)
  int d(1);
  for(int c(0);c<mu;c++)
    d *= 2 ;

  int g5(Ns*Ns - 1) ;

  LatticeComplex  C(0.0) ;
 
  for(int s(0); s<N5;s++){
    // first the 1-gamma_mu term 
    C=.5*(trace(shift(p5d[s],FORWARD,mu)*adj(p5d[s])*Gamma(g5)*Gamma(d)*u[mu])-
	  trace(shift(p5d[s],FORWARD,mu)*adj(p5d[s])*Gamma(g5)         *u[mu])+
	  //now th 1+gamma_mu term
	  trace(p5d[s]*shift(adj(p5d[s]),FORWARD,mu)*Gamma(g5)*Gamma(d)*adj(u[mu]))+
	  trace(p5d[s]*shift(adj(p5d[s]),FORWARD,mu)*Gamma(g5)         *adj(u[mu]))
	  );
	  
      
    if(s<N5/2) 
      corr -= C ;
    else
      corr += C ;
  }  
}


//! Given a complete propagator as a source, this does all the inversions needed
/*! \ingroup qprop
 *
 * This routine is actually generic to Domain Wall fermions (Array) fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void dwf_quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    const EvenOddPrecDWFermActBaseArray<LatticeFermion>& S_f,
		    Handle<const ConnectState> state,
		    enum InvType invType,
		    const Real& RsdCG, 
		    int MaxCG, int& ncg_had)
{
  dwf_quarkProp4_a(q_sol, xml_out, q_src, S_f, state, invType, RsdCG, MaxCG, ncg_had);
}

