// $Id: Ql_3pt_w.cc,v 1.5 2009-04-09 22:57:38 caubin Exp $ 
/*! \file
 *  \brief Heavy-Light 3pt function
 */

#include "Ql_3pt_w.h"
#include "barQll_w.h"

namespace Chroma {

//! Heavy-light meson 2-pt function
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * The heavy quark is inserted in the infinitely heavy quark limit
 * by a Wilson-Line without spin indices.
 * We are effectively propagating a spin-1/2 light quark
 * This generates with two quark propagators all three-point
 * functions, with all 16 gamma matrix insertions, of the form:
 *
 * Mat1 = Sc(z-y) [Gamma] Su(y-z) gamma5
 * Mat2 = Sb(x-z) gamma5 Ss(z-x) [Gamma]
 * 
 * And we output the 16 gamma's, all possible trace combinations (1 or 2 color traces and 1 or 
 * 2 spin traces
 *                                                                                             
 * \param u                  gauge field (Read)                                                
 * \param quark_propagator1   quark propagator1 ( Read )                                 
 * \param quark_propagator2   quark propagator2 ( Read )                             
 * \param heavy_quark_propagator1   heavy quark propagator1 ( Read; Opt, if not static  ) 
 * \param heavy_quark_propagator2   heavy quark propagator2 ( Read; Opt, if not static )
 * \param src_coord          cartesian coordinates of the source ( Read )        
 * \param snk_coord          cartesian coordinates of the sink ( Read )          
 * \param heavy_src          cartesian coordinates of the heavy quark if not static ( Read )
 *                           Not needed if both are static or both are not static.
 * \param phases             object holds list of momenta and Fourier phases ( Read ) 
 * \param xml                xml file object ( Read )  
 * \param xml_group          group name for xml data ( Read ) 
 *
 * \param bc                 if there, the bc's for the static quark props...
 */
  
void QlQl(const multi1d<LatticeColorMatrix>& u,
	  const LatticePropagator& quark_propagator1,
	  const LatticePropagator& quark_propagator2,
	  const multi1d<int>& src_coord,
	  const multi1d<int>& snk_coord,
	  const int& bc, 
	  const SftMom& phases,
	  XMLWriter& xml,
	  const string& xml_group)
{
    START_CODE();
  
  if ( Ns != 4 )		/* Code is specific to Ns=4 */
    return;

  int length = phases.numSubsets() ;
  int num_mom = phases.numMom();
  
  LatticeColorMatrix Qprop1; // b quark
  LatticeColorMatrix Qprop2; // c quark

  HeavyQuarkProp(Qprop1,u,src_coord,length,bc);
  HeavyQuarkPropBack(Qprop2,u,snk_coord,length,bc);

  // S_proj_unpol = (1/2)(1 + gamma_4)
  // This is the heavy quark dirac structure, but I am using the notation
  // that is currently used in Chroma!!!!!!!!

  SpinMatrix g_one = 1.0;
  SpinMatrix S_proj_unpol = 0.5 * (g_one + (g_one * Gamma(8)));
  multi1d<DComplex> Hq;
  
  for(int gind=0; gind<16;++gind){// Loop over gamma matrices...
    LatticePropagator Mat1 , Mat2;
    LatticeComplex Hq_prop11 , Hq_prop12, Hq_prop21, Hq_prop22 ;

    Mat1 = quark_propagator1 * adj(Qprop1*S_proj_unpol) * Gamma(15) * Gamma(gind) ;
    Mat2 = Qprop2*S_proj_unpol * adj(quark_propagator2) * Gamma(15) * Gamma(gind);
    //First is number of spin traces , second is number of color traces
    Hq_prop11 = trace(Mat1*Mat2);
    Hq_prop12 = traceSpin(traceColor(Mat1)*traceColor(Mat2));
    Hq_prop21 = traceColor(traceSpin(Mat1)*traceSpin(Mat2));
    Hq_prop22 = trace(Mat1)*trace(Mat2);
    
    // Project onto zero momentum
    multi2d<DComplex> hsumHq11, hsumHq12 , hsumHq21, hsumHq22 ;
    hsumHq11 = phases.sft(Hq_prop11);
    hsumHq12 = phases.sft(Hq_prop12);
    hsumHq21 = phases.sft(Hq_prop21);
    hsumHq22 = phases.sft(Hq_prop22);
    multi2d<DComplex> HQprop11(num_mom,length) ;
    multi2d<DComplex> HQprop12(num_mom,length) ;
    multi2d<DComplex> HQprop21(num_mom,length) ;
    multi2d<DComplex> HQprop22(num_mom,length) ;
    
    for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
      for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - src_coord[Nd-1] + length) % length;
	  HQprop11[sink_mom_num][t_eff] = hsumHq11[sink_mom_num][t];
	  HQprop12[sink_mom_num][t_eff] = hsumHq12[sink_mom_num][t];
	  HQprop21[sink_mom_num][t_eff] = hsumHq21[sink_mom_num][t];
	  HQprop22[sink_mom_num][t_eff] = hsumHq22[sink_mom_num][t];
	}
    
    //XMLWriter xml_bar(xml);
    push(xml, xml_group);
    write(xml, "gamma_value", gind);
    write(xml, "ThreePoint1Spin1ColorTrace", HQprop11[0]);
    write(xml, "ThreePoint1Spin2ColorTrace", HQprop12[0]);
    write(xml, "ThreePoint2Spin1ColorTrace", HQprop21[0]);
    write(xml, "ThreePoint2Spin2ColorTrace", HQprop22[0]);

    pop(xml);
    
    END_CODE();
  }// gind loop
}// QlQl routine
  
  
void QlQl(const multi1d<LatticeColorMatrix>& u,
	  const LatticePropagator& quark_propagator1,
	  const LatticePropagator& quark_propagator2,
	  const LatticePropagator& heavy_quark_propagator,
	  const multi1d<int>& src_coord,
	  const multi1d<int>& snk_coord,
	  const multi1d<int>& heavy_src, 
	  const int& bc,
	  const SftMom& phases,
	  XMLWriter& xml,
	  const string& xml_group)
{
    START_CODE();
  
  if ( Ns != 4 )		/* Code is specific to Ns=4 */
    return;

  int length = phases.numSubsets() ;
  int num_mom = phases.numMom();
  
  LatticeColorMatrix Qprop; // b quark

  if (heavy_src==src_coord)
    HeavyQuarkProp(Qprop,u,snk_coord,length,bc);
  else 
    HeavyQuarkProp(Qprop,u,src_coord,length,bc);

  // S_proj_unpol = (1/2)(1 + gamma_4)
  // This is the heavy quark dirac structure, but I am using the notation
  // that is currently used in Chroma!!!!!!!!

  SpinMatrix g_one = 1.0;
  SpinMatrix S_proj_unpol = 0.5 * (g_one + (g_one * Gamma(8)));
  multi1d<DComplex> Hq;
  
  for(int gind=0; gind<16;++gind){// Loop over gamma matrices...
    LatticePropagator Mat1 , Mat2;
    LatticeComplex Hq_prop11 , Hq_prop12, Hq_prop21, Hq_prop22 ;

	if (heavy_src==src_coord){
	  Mat1 = quark_propagator1 * adj(heavy_quark_propagator) * Gamma(15) * Gamma(gind) ;
	  Mat2 = Qprop*S_proj_unpol * adj(quark_propagator2) * Gamma(15) * Gamma(gind);
	}
	else{
	  Mat1 = quark_propagator1 * adj(Qprop*S_proj_unpol) * Gamma(15) * Gamma(gind) ;
	  Mat2 = heavy_quark_propagator * adj(quark_propagator2) * Gamma(15) * Gamma(gind);
    }

    Hq_prop11 = trace(Mat1*Mat2);
    Hq_prop12 = traceSpin(traceColor(Mat1)*traceColor(Mat2));
    Hq_prop21 = traceColor(traceSpin(Mat1)*traceSpin(Mat2));
    Hq_prop22 = trace(Mat1)*trace(Mat2);
    
    // Project onto zero momentum
    multi2d<DComplex> hsumHq11, hsumHq12 , hsumHq21, hsumHq22 ;
    hsumHq11 = phases.sft(Hq_prop11);
    hsumHq12 = phases.sft(Hq_prop12);
    hsumHq21 = phases.sft(Hq_prop21);
    hsumHq22 = phases.sft(Hq_prop22);
    multi2d<DComplex> HQprop11(num_mom,length) ;
    multi2d<DComplex> HQprop12(num_mom,length) ;
    multi2d<DComplex> HQprop21(num_mom,length) ;
    multi2d<DComplex> HQprop22(num_mom,length) ;
    
    for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
      for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - src_coord[Nd-1] + length) % length;
	  HQprop11[sink_mom_num][t_eff] = hsumHq11[sink_mom_num][t];
	  HQprop12[sink_mom_num][t_eff] = hsumHq12[sink_mom_num][t];
	  HQprop21[sink_mom_num][t_eff] = hsumHq21[sink_mom_num][t];
	  HQprop22[sink_mom_num][t_eff] = hsumHq22[sink_mom_num][t];
	}
    
    //XMLWriter xml_bar(xml);
    push(xml, xml_group);
    write(xml, "gamma_value", gind);
    write(xml, "ThreePoint1Spin1ColorTrace", HQprop11[0]);
    write(xml, "ThreePoint1Spin2ColorTrace", HQprop12[0]);
    write(xml, "ThreePoint2Spin1ColorTrace", HQprop21[0]);
    write(xml, "ThreePoint2Spin2ColorTrace", HQprop22[0]);
    pop(xml);
    
    END_CODE();
  }// gind loop
}// QlQl routine
  
  /***

  This version is for both heavy quarks being not static (thus, they could be light quarks).
   **/  
void QlQl(const multi1d<LatticeColorMatrix>& u,
	    const LatticePropagator& quark_propagator1,
	    const LatticePropagator& quark_propagator2,
	    const LatticePropagator& heavy_quark_propagator1,
	    const LatticePropagator& heavy_quark_propagator2,
	    const multi1d<int>& src_coord,
	    const multi1d<int>& snk_coord,
	    const SftMom& phases,
	    XMLWriter& xml,
	    const string& xml_group)
{
    START_CODE();
  
  if ( Ns != 4 )		/* Code is specific to Ns=4 */
    return;

  int length = phases.numSubsets() ;
  int num_mom = phases.numMom();
  
  multi1d<DComplex> Hq;
  
  for(int gind=0; gind<16;++gind){// Loop over gamma matrices...
    LatticePropagator Mat1 , Mat2;
    LatticeComplex Hq_prop11 , Hq_prop12, Hq_prop21, Hq_prop22 ;

    Mat1 = quark_propagator1 * adj(heavy_quark_propagator1) * Gamma(15) * Gamma(gind);
    Mat2 = heavy_quark_propagator2 * adj(quark_propagator2) * Gamma(15) * Gamma(gind);
    
    Hq_prop11 = trace(Mat1*Mat2);
    Hq_prop12 = traceSpin(traceColor(Mat1)*traceColor(Mat2));
    Hq_prop21 = traceColor(traceSpin(Mat1)*traceSpin(Mat2));
    Hq_prop22 = trace(Mat1)*trace(Mat2);
    
    // Project onto zero momentum
    multi2d<DComplex> hsumHq11, hsumHq12 , hsumHq21, hsumHq22 ;
    hsumHq11 = phases.sft(Hq_prop11);
    hsumHq12 = phases.sft(Hq_prop12);
    hsumHq21 = phases.sft(Hq_prop21);
    hsumHq22 = phases.sft(Hq_prop22);
    multi2d<DComplex> HQprop11(num_mom,length) ;
    multi2d<DComplex> HQprop12(num_mom,length) ;
    multi2d<DComplex> HQprop21(num_mom,length) ;
    multi2d<DComplex> HQprop22(num_mom,length) ;
    
    for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
      for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - src_coord[Nd-1] + length) % length;
	  HQprop11[sink_mom_num][t_eff] = hsumHq11[sink_mom_num][t];
	  HQprop12[sink_mom_num][t_eff] = hsumHq12[sink_mom_num][t];
	  HQprop21[sink_mom_num][t_eff] = hsumHq21[sink_mom_num][t];
	  HQprop22[sink_mom_num][t_eff] = hsumHq22[sink_mom_num][t];
	}
    
    //XMLWriter xml_bar(xml);
    push(xml, xml_group);
    write(xml, "gamma_value", gind);
    write(xml, "ThreePoint1Spin1ColorTrace", HQprop11[0]);
    write(xml, "ThreePoint1Spin2ColorTrace", HQprop12[0]);
    write(xml, "ThreePoint2Spin1ColorTrace", HQprop21[0]);
    write(xml, "ThreePoint2Spin2ColorTrace", HQprop22[0]);
    pop(xml);
    
    END_CODE();
  }// gind loop
}// QlQl routine
  
  
}// Namespace Chroma
