// $Id: barQll_w.cc,v 1.1 2004-09-08 14:00:13 kostas Exp $ 
/*! \file
 *  \brief Heavy Baryon (Qll)  2-pt function : Orginos and Savage
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/barQll_w.h"

using namespace QDP;
 

//! Lambdaq and SigmaQ 2-pt functions
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * Construct baryon propagators for the LambdaQ and SigmaQ(*) with
 * degenerate "u" and "d" quarks.  In the heavy quark limit the Sigma
 * and Sigma* are degenerate.  The heavy quark is inserted in the infinitely heavy quark limit
 * by a Wilson-Line without spin indices.
 * We are effectively propagating a spin-0 diquark and a spin-1 diquark.
 *
 * \param u                  gauge field (Read) 
 * \param quark_propagator   quark propagator ( Read )
 * \param src_coord          cartesian coordinates of the source ( Read )
 * \param phases             object holds list of momenta and Fourier phases ( Read )
 * \param xml                xml file object ( Read )
 * \param xml_group          group name for xml data ( Read )
 *
 */

void Qll(const multi1d<LatticeColorMatrix>& u, 
	 const LatticePropagator& quark_propagator,
	 const multi1d<int>& src_coord, 
	 const SftMom& phases,
       	 XMLWriter& xml,
	 const string& xml_group)
{
  START_CODE();

  if ( Ns != 4 || Nc != 3 )		/* Code is specific to Ns=4 and Nc=3. */
    return;

  LatticeColorMatrix Qprop;

  HeavyQuarkProp(Qprop,u,src_coord);

  multi1d<DComplex> barLamQ;
  multi1d<DComplex> barSigQ;
  LatticePropagator di_quark;
  LatticeComplex LamQ_prop, SigQx_prop, SigQy_prop, SigQz_prop;
 

  // LambdaQ  :  |LambdaQ> = (d C gamma_5 u) Q

  di_quark = quarkContract13(quark_propagator * Gamma(5),
			     Gamma(5) * quark_propagator); //spin-0 diquark

  LamQ_prop = traceColor(Qprop * traceSpin(di_quark));
 
 
 // SigmQ  :  |SigmaQ> = (d C gamma_mu u) Q  : The SigmaQ and SigmaQ* are degenerate!!!!

  di_quark = quarkContract13(quark_propagator * Gamma(11),
			     Gamma(11) * quark_propagator); //spin-1 diquark oriented in x-direction

  SigQx_prop = traceColor(Qprop * traceSpin(di_quark));
 

  di_quark = quarkContract13(quark_propagator * Gamma(8),
			     Gamma(8) * quark_propagator); //spin-1 diquark oriented in y-direction

  SigQy_prop = traceColor(Qprop * traceSpin(di_quark));
 

  di_quark = quarkContract13(quark_propagator * Gamma(14),
			     Gamma(14) * quark_propagator); //spin-1 diquark oriented in z-direction

  SigQz_prop = traceColor(Qprop * traceSpin(di_quark));

 
  // Project onto zero momentum
  multi2d<DComplex> hsumLamQ,hsumSigQx,hsumSigQy,hsumSigQz;
  hsumLamQ = phases.sft(LamQ_prop);
  hsumSigQx = phases.sft(SigQx_prop);
  hsumSigQy = phases.sft(SigQy_prop);
  hsumSigQz = phases.sft(SigQz_prop);

  for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
    for(int t = 0; t < Nt; ++t)
      {
	int t_eff = (t - t0 + Nt) % Nt;
	LQprop[sink_mom_num][t_eff]  = hsumLamQ[sink_mom_num][t];
	SQxprop[sink_mom_num][t_eff] = hsumSigQx[sink_mom_num][t];
	SQyprop[sink_mom_num][t_eff] = hsumSigQy[sink_mom_num][t];
	SQzprop[sink_mom_num][t_eff] = hsumSigQz[sink_mom_num][t];
      }


  XMLArrayWriter xml_bar(xml,4);
  push(xml_bar, xml_group);
  push(xml_bar, "LambdaQ");
  write(xml_bar, "prop", LQprop);
  pop(xml_bar);
  push(xml_bar, "SigmaQx");
  write(xml_bar, "prop", SQxprop);
  pop(xml_bar);
  push(xml_bar, "SigmaQy");
  write(xml_bar, "prop", SQyprop);
  pop(xml_bar);
  push(xml_bar, "SigmaQz");
  write(xml_bar, "prop", SQzprop);
  pop(xml_bar);
  pop(xml_bar);

  END_CODE();
}



//! Heavy Quark Propagator
/*!
 * \ingroup hadron
 *
 * 
 * This constructs the propagator for a spinless Wilson-Line propagating from the 
 * point src_coord forward in time, and vanishing on previous time-slices.
 *
 * \param Qprop              Wilson-Line (write)
 * \param u                  Gauge Field (Read) 
 * \param src_coord          cartesian coordinates of the source ( Read )
 *
 */

void HeavyQuarkProp(LatticeColorMatrix& Qprop,
		    const multi1d<LatticeColorMatrix>& u,
		    const multi1d<int>& src_coord)
{

  UnorderedSet slice ;
  slice.make(TimeSliceFunc(Nd-1));

  Qprop = 0.0 ;
  
  ColorMatrix one = 1.0 ;

  pokeSite(Qprop,one,src_coord);

  for(int t(t0+1);t<Nt;t++){
    Qprop[slice[t]] = shift(u,BACKWARD,3)*shift(Qprop,BACKWARD,3); 
  }
 
  
}

