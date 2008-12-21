// $Id: mesQlPOT_w.cc,v 1.2 2008-12-21 21:22:36 edwards Exp $ 
/*! \file
 *  \brief Potential between 2 heavy mesons : Orginos and Savage
 *
 *   V 3 
 *
 * Revised by Savage Sept 8 2005 , S=0 states included.
 *
 * Revised by Savage March 20 2005 , included the case with vanishing separation.
 *
 * Revised by Savage March 18 2005 , d<Nd-2 -> d<Nd-1 for the quark shifting 
 * in  amplitude 2. And changed adj(spin..) to spin... in amp2
 *
 */

#include "barQll_w.h"
#include "mesQl_w.h"
#include "mesQlPOT_w.h"

namespace Chroma
{
  //! Heavy-light meson potential
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   * Construct propagators for two heavy mesons in all combinations of spin up and 
   * spin down light degrees of freedom..
   * In the heavy quark limit the D and the D* are degenerate.  
   * The heavy quark is inserted in the infinitely heavy quark limit
   * by a Wilson-Line without spin indices.
   * We are effectively propagating two  spin-1/2 light degrees of freedom.
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

  void QlQlPOT(const multi1d<LatticeColorMatrix>& u, 
	       const LatticePropagator& quark1,
	       const LatticePropagator& quark2,
	       const multi1d<int>& src1, 
	       const multi1d<int>& src2, 
	       const SftMom& phases,
	       XMLWriter& xml,
	       const string& xml_group)
  {
    START_CODE();
  
    if ( Ns != 4 )		/* Code is specific to Ns=4 */
      return;

    if (src1[Nd-1] != src2[Nd-1] ){
      QDPIO::cerr<<"Sources must be at the same time slice\n";
      return ;
    }

    int length = phases.numSubsets() ;
    int num_mom = phases.numMom();
  
    LatticeColorMatrix Qprop1;
    LatticeColorMatrix Qprop2;

    HeavyQuarkProp(Qprop1,u,src1,length);
    HeavyQuarkProp(Qprop2,u,src2,length);


    LatticeComplex   Hq1u , Hq2u  ,Hq1d , Hq2d ;

    // S_proj = (1 + \Sigma_3)*(1 + gamma_4) / 2 
    //        = (1 + Gamma(8) - i G(3) - i G(11)) / 2

    SpinMatrix g_one = 1.0;
    SpinMatrix S_proj_up = 
      0.5*((g_one + Gamma(8) * g_one) - timesI(Gamma(3) * g_one  +  Gamma(11) * g_one));
    SpinMatrix S_proj_down = 
      0.5*((g_one + Gamma(8) * g_one) + timesI(Gamma(3) * g_one  +  Gamma(11) * g_one));
    SpinMatrix S_proj_unpol = 
      0.5*((g_one + Gamma(8) * g_one));
  
    SpinMatrix cg5 = g_one * Gamma(5);
    SpinMatrix STcg5 = transpose(cg5);
  
    SpinMatrix  spinsinglet = g_one*S_proj_unpol*Gamma(5);
    SpinMatrix STspinsinglet = transpose(spinsinglet);


    // Single meson states, spin up and down at location 1 and 2 

    Hq1u = trace(Qprop1 * S_proj_up * adj(quark1) ) ;
    Hq2u = trace(Qprop2 * S_proj_up * adj(quark2) ) ;
    Hq1d = trace(Qprop1 * S_proj_down * adj(quark1) ) ;
    Hq2d = trace(Qprop2 * S_proj_down * adj(quark2) ) ;

    // fourier transform
    multi2d<DComplex> hsumHq1u ,  hsumHq2u , hsumHq1d ,  hsumHq2d ;
    hsumHq1u = phases.sft(Hq1u);
    hsumHq2u = phases.sft(Hq2u);
    hsumHq1d = phases.sft(Hq1d);
    hsumHq2d = phases.sft(Hq2d);
  
    // Make time slices for xml reader
    multi2d<DComplex> Hq1uprop(num_mom,length)    , Hq2uprop(num_mom,length)   , Hq1dprop(num_mom,length)     , Hq2dprop(num_mom,length) ;
    multi2d<DComplex> Hq1prop(num_mom,length)    , Hq2prop(num_mom,length)  ;
    for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
      for(int t = 0; t < length; ++t)
      {
	int t_eff = (t - src1[Nd-1] + length) % length;
	Hq1uprop[sink_mom_num][t_eff]  = hsumHq1u[sink_mom_num][t];
	Hq2uprop[sink_mom_num][t_eff]  = hsumHq2u[sink_mom_num][t];
	Hq1dprop[sink_mom_num][t_eff]  = hsumHq1d[sink_mom_num][t];
	Hq2dprop[sink_mom_num][t_eff]  = hsumHq2d[sink_mom_num][t];
	Hq1prop[sink_mom_num][t_eff]  = 0.5000000000 * (hsumHq1u[sink_mom_num][t] + hsumHq1d[sink_mom_num][t]);
	Hq2prop[sink_mom_num][t_eff]  = 0.5000000000 * (hsumHq2u[sink_mom_num][t] + hsumHq2d[sink_mom_num][t]);
      }

    QDPIO::cerr<<"Done with single meson correlators \n";





    LatticePropagator tmp , AdjQuark2Qprop2 ;
    tmp = adj(quark2)*Qprop2 ;
    AdjQuark2Qprop2=tmp;
    for (int d(0);d<Nd-1;d++){
      int r = src2[d] - src1[d] ;
      while(r>0){
	AdjQuark2Qprop2 = shift(tmp,FORWARD,d);
	tmp = AdjQuark2Qprop2 ;
	r--;
      }
      while(r<0){
	AdjQuark2Qprop2 = shift(tmp,BACKWARD,d);
	tmp = AdjQuark2Qprop2 ;
	r++;
      }
    }
    QDPIO::cerr<<" Done shifting q2 Q2\n";

    LatticePropagator tmp1,AdjQuark1Qprop2 ;
    tmp1 = adj(quark1)*Qprop2 ;
    AdjQuark1Qprop2=tmp1;
    for (int d(0);d<Nd-1;d++){
      int r = src2[d] - src1[d] ;
      while(r>0){
	AdjQuark1Qprop2 = shift(tmp1,FORWARD,d);
	tmp1 = AdjQuark1Qprop2 ;
	r--;
      }
      while(r<0){
	AdjQuark1Qprop2 = shift(tmp1,BACKWARD,d);
	tmp1 = AdjQuark1Qprop2 ;
	r++;
      }
    }
    QDPIO::cerr<<" Done shifting q1 Q2\n";



    // S=1 correlators

    LatticeComplex  A1S1Sp1 , A1S1Sm1 , A2S1Sp1 , A2S1Sm1;
    A1S1Sp1 = trace(S_proj_up * Qprop1 * adj(quark1)) * trace(S_proj_up * AdjQuark2Qprop2 );
    A1S1Sm1 = trace(S_proj_down * Qprop1 * adj(quark1)) * trace(S_proj_down * AdjQuark2Qprop2 );

    A2S1Sp1 = -trace(Qprop1 * AdjQuark1Qprop2 * S_proj_up* adj(quark2) * S_proj_up );
    A2S1Sm1 = -trace(Qprop1 * AdjQuark1Qprop2 * S_proj_down * adj(quark2) * S_proj_down );

    // Project S=1 correlators onto zero momentum
    multi2d<DComplex> hsumA1S1Sp1 ,  hsumA2S1Sp1 ,  hsumA1S1Sm1 , hsumA2S1Sm1 ;
    hsumA1S1Sp1 =  phases.sft(A1S1Sp1);
    hsumA1S1Sm1 =  phases.sft(A1S1Sm1);
    hsumA2S1Sp1 =  phases.sft(A2S1Sp1);
    hsumA2S1Sm1 =  phases.sft(A2S1Sm1);

    // Make time slices for xml reader
    multi2d<DComplex> A1S1Sp1prop(num_mom,length) , A1S1Sm1prop(num_mom,length) , A2S1Sp1prop(num_mom,length) , A2S1Sm1prop(num_mom,length);
    multi2d<DComplex> I1S1prop(num_mom,length) , I0S1prop(num_mom,length);

    for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
      for(int t = 0; t < length; ++t)
      {
	int t_eff = (t - src1[Nd-1] + length) % length;
	A1S1Sp1prop[sink_mom_num][t_eff]  = hsumA1S1Sp1[sink_mom_num][t];
	A1S1Sm1prop[sink_mom_num][t_eff]  = hsumA1S1Sm1[sink_mom_num][t];
	A2S1Sp1prop[sink_mom_num][t_eff]  = hsumA2S1Sp1[sink_mom_num][t];
	A2S1Sm1prop[sink_mom_num][t_eff]  = hsumA2S1Sm1[sink_mom_num][t];
	I1S1prop[sink_mom_num][t_eff]  = 0.50000000000000 * ( hsumA1S1Sp1[sink_mom_num][t] +  hsumA1S1Sm1[sink_mom_num][t] 
							      +  hsumA2S1Sp1[sink_mom_num][t] +  hsumA2S1Sm1[sink_mom_num][t] ); 
	I0S1prop[sink_mom_num][t_eff]  = 0.50000000000000 * ( hsumA1S1Sp1[sink_mom_num][t] +  hsumA1S1Sm1[sink_mom_num][t] 
							      -  hsumA2S1Sp1[sink_mom_num][t] -  hsumA2S1Sm1[sink_mom_num][t] ); 
      }
    QDPIO::cerr<<" Done with S=1 correlators\n";



    // S=0 correlators
    LatticeComplex  A1S0 , A2S0;


    LatticePropagator STL2Q2 , STL2 , adjL2;
    SpinTranspose(AdjQuark2Qprop2 , STL2Q2 ) ;
    adjL2 = adj(quark2);
    SpinTranspose(adjL2 , STL2 ) ;

    A1S0 = trace( traceColor(adj(quark1) * Qprop1) * spinsinglet *  traceColor(STL2Q2) * STspinsinglet );
    A2S0 = -trace( AdjQuark1Qprop2 * spinsinglet * STL2 * Qprop1 * spinsinglet );

    // Project S=0 correlators onto zero momentum
    multi2d<DComplex> hsumA1S0 ,  hsumA2S0;
    hsumA1S0 =  phases.sft(A1S0);
    hsumA2S0 =  phases.sft(A2S0);

    // Make time slices for xml reader
    multi2d<DComplex> A1S0prop(num_mom,length) , A2S0prop(num_mom,length) ;
    multi2d<DComplex> I1S0prop(num_mom,length) , I0S0prop(num_mom,length);
   
    for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
      for(int t = 0; t < length; ++t)
      {
	int t_eff = (t - src1[Nd-1] + length) % length;
	A1S0prop[sink_mom_num][t_eff]  = hsumA1S0[sink_mom_num][t];
	A2S0prop[sink_mom_num][t_eff]  = hsumA2S0[sink_mom_num][t];
	I1S0prop[sink_mom_num][t_eff]  = hsumA1S0[sink_mom_num][t] +  hsumA2S0[sink_mom_num][t]  ; 
	I0S0prop[sink_mom_num][t_eff]  = hsumA1S0[sink_mom_num][t] -  hsumA2S0[sink_mom_num][t]  ; 
      }
    QDPIO::cerr<<" Done with S=0 correlators\n";




    // XMLWriter xml_bar(xml);
    push(xml, xml_group);
    write(xml, "Meson1U",    Hq1uprop[0]);
    write(xml, "Meson1D",    Hq1dprop[0]);
    write(xml, "Meson2U",    Hq2uprop[0]);
    write(xml, "Meson2D",    Hq2dprop[0]);
    write(xml, "Meson1",    Hq1prop[0]);
    write(xml, "Meson2",    Hq2prop[0]);
    write(xml, "A1S1Sp1",  A1S1Sp1prop[0] );
    write(xml, "A1S1Sm1",  A1S1Sm1prop[0] );
    write(xml, "A2S1Sp1",  A2S1Sp1prop[0] );
    write(xml, "A2S1Sm1",  A2S1Sm1prop[0] );
    write(xml, "A1S0",  A1S0prop[0] );
    write(xml, "A2S0",  A2S0prop[0] );
    write(xml, "I1S1",  I1S1prop[0] );
    write(xml, "I0S1",  I0S1prop[0] );
    write(xml, "I1S0",  I1S0prop[0] );
    write(xml, "I0S0",  I0S0prop[0] );
    pop(xml);

    END_CODE();
  }



  //!  Spin Transpose Function
  /*!
   * \ingroup hadron
   *
   * This is a dumb way of taking the spin transpose of a propagator,
   * while leaving all other indices untouched...suggested to 
   * us by Bob Edwards
   * \param prop                input propagator
   * \param STprop              spin transposed propagator
   *
   */

  void SpinTranspose(const LatticePropagator& prop , LatticePropagator& STprop)
  {
    for(int j = 0; j < Ns; ++j){
      for(int i = 0; i < Ns; ++i){
	pokeSpin(STprop, peekSpin(prop, i, j), j, i);
      }
    }
  }

}
