// $Id: dwf_quarkprop4_w.cc,v 1.1 2003-12-10 21:26:07 kostas Exp $
/*! \file
 * \brief Full quark propagator solver for domain wall fermions
 *
 * Given a complete propagator as a source, this does all the inversions needed
 */


#include "chromabase.h"
#include "util/util.h"
#include "fermact.h"
#include "actions/ferm/qprop/quarkprop4_w.h"

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
void dwf_quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    const C<T>& S_f,
		    const ConnectState& state,
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

	   tmp *= fact; 
	   
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
	   psi *= fact;

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

   pop(xml_out);

 
   // constuct the conserved axial current correlator
   LatticeComplex cfield ;
   dwf_conserved_axial_ps_corr(cfield,state.getLinks(),prop5d,3,S_f.size());
			       
   multi2d<DComplex> corr ;
   
   SftMom trick(0,false,3) ;
   
   corr = trick.sft(cfield);

   //push(xml_out, "t_dir");
   //Write(xml_out, 3);
   //pop(xml_out);
   push(xml_out, "DWF_ConservedAxial");
   Write(xml_out, corr);
   pop(xml_out);

   //Now the midpoint Pseudoscalar correlator
   cfield = trace(prop_mp * adj(prop_mp));

   corr = trick.sft(cfield);
   push(xml_out, "DWF_MidPoint_Pseudo");
   Write(xml_out, corr);
   pop(xml_out);


   END_CODE("dwf_quarkProp4");
}



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

  LatticeComplex  C(0.0) ;
  for(int s(0); s<N5;s++)
    {
      // first the 1-gamma_4 term 
      C=0.5*(trace(shift(p5d[s],FORWARD,mu) * Gamma(d)*u[mu]* adj(p5d[2]))-
	     trace(shift(p5d[s],FORWARD,mu) *          u[mu]* adj(p5d[2]))+
	     //now th 1+gamma_4 term
	     trace(p5d[s] *Gamma(d)*adj(u[mu])*shift(adj(p5d[s]),FORWARD,mu))+
	     trace(p5d[s] *         adj(u[mu])*shift(adj(p5d[s]),FORWARD,mu))
	     ) ;
      
      if(s<N5/2) 
	corr -= C ;
      else
	corr += C ;
    }  
}

