//  $Id: qqbar_w.cc,v 3.3 2008-12-21 21:22:37 edwards Exp $
//  $Log: qqbar_w.cc,v $
//  Revision 3.3  2008-12-21 21:22:37  edwards
//  Some code cleanups. Put in chroma namespace. Moved around the chromabase.h
//  include.
//
//  Revision 3.2  2008/03/29 03:40:20  kostas
//  added vector mesons
//
//  Revision 3.1  2007/02/22 21:11:49  bjoo
//  Removed Ordered and Unordered Subsets and Sets. Now just have Subsets and Sets -not Unordered or Ordered - passes regressions with suitable QDP where QDP has SSE disabled
//
//  Revision 3.0  2006/04/03 04:59:00  edwards
//  Major overhaul of fermion and gauge action interface. Basically,
//  all fermacts and gaugeacts now carry out  <T,P,Q>  template parameters. These are
//  the fermion type, the "P" - conjugate momenta, and "Q" - generalized coordinates
//  in the sense of Hamilton's equations. The fermbc's have been rationalized to never
//  be over multi1d<T>. The "createState" within the FermionAction is now fixed meaning
//  the "u" fields are now from the coordinate type. There are now "ConnectState" that
//  derive into FermState<T,P,Q> and GaugeState<P,Q>.
//
//  Revision 2.1  2005/12/27 20:41:41  kostas
//  added NPLQCD code
//
//  constructs 2 quark propagators contracted at the sink
//

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "qqbar_w.h"

namespace Chroma
{

  //! Meson-Meson 4-pt functions
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * Construct meson-meson propagators
   * The two propagators can be identical or different.
   *
   * \param qqbar           -- the 2-quark propagator ( Write )
   * \param quark_prop_1 -- first quark propagator ( Read )
   * \param quark_prop_2 -- second (anti-) quark propagator ( Read )
   * \param t0 -- timeslice coordinate of the source ( Read )
   * \param phases -- object holds list of momenta and Fourier phases ( Read )
   *
   *          ____
   *          \                               +
   *qqbar(p,t)=> [g5 q2(t_src;t + t_src,x) g5] g5 q1(t+t_src,x;t_src) *  exp(ipx)
   *          /
   *          ----
   *            x
   */

/*** chroma bug
     void compute_qqbar( multi2d<DPropagator>& qqbar,
     const LatticePropagator& quark_prop_1,
     const LatticePropagator& quark_prop_2, 
     const SftMom& phases,
     int t0)
     {
     START_CODE();

     QDPIO::cout<<"Starting the qqbar code\n";

     // Length of lattice in decay direction
     Set sft_set(phases.getSet()) ;
     int length(sft_set.numSubsets());
     QDPIO::cout<<"Time length: "<<length<<endl ;

     // Construct the anti-quark propagator from quark_prop_2
     int G5 = Ns*Ns-1;
     LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

     LatticePropagator quarkloop;

     quarkloop = adj(anti_quark_prop)*Gamma(G5) *quark_prop_1 ;
  
     multi1d<DPropagator> foo(length);
     for (int mom_num(0); mom_num < phases.numMom(); ++mom_num){
     foo = sumMulti(phases[mom_num]*quarkloop, sft_set);
     for(int t = 0; t < length; ++t){
     int t_eff = (t - t0 + length) % length;
     qqbar[mom_num][t_eff]= foo[t] ;
     }
     }
    
     QDPIO::cout<<"Finished the qqbar code\n";

     END_CODE();
     }
***/

/*** bug work around ***/
  void compute_qqbar( multi2d<DPropagator>& qqbar,
		      const LatticePropagator& quark_prop_1,
		      const LatticePropagator& quark_prop_2, 
		      const SftMom& phases,
		      int t0)
  {
    START_CODE();
  
    QDPIO::cout<<"Starting the qqbar code\n";

    // Length of lattice in decay direction
    Set sft_set(phases.getSet()) ;
    int length(sft_set.numSubsets());
    //QDPIO::cout<<"Time length: "<<length<<endl ;

    // Construct the anti-quark propagator from quark_prop_2
    int G5 = Ns*Ns-1;
    LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

    LatticePropagator quarkloop;
    LatticeSpinMatrix sm ;
    LatticeColorMatrix cm ;

    quarkloop = adj(anti_quark_prop)*Gamma(G5) *quark_prop_1 ;
    LatticeComplex cc;
    multi2d<DPropagator> foo(phases.numMom(),length);
    multi2d<DColorMatrix> dcm(phases.numMom(),length);

    for(int s1(0);s1<Ns;s1++)
      for(int s2(0);s2<Ns;s2++){
	cm = peekSpin(quarkloop,s1,s2);
	for(int c1(0);c1<Nc;c1++)
	  for(int c2(0);c2<Nc;c2++){
	    cc = peekColor(cm,c1,c2);
	    multi2d<DComplex> fcc(phases.sft(cc));
	    for (int mom_num(0); mom_num < phases.numMom(); ++mom_num){
	      for(int t = 0; t < length; ++t){
		pokeColor(dcm[mom_num][t],fcc[mom_num][t],c1,c2);
	      }
	    }
	  }
	for (int mom_num(0); mom_num < phases.numMom(); ++mom_num)
	  for(int t = 0; t < length; ++t)
	    pokeSpin(foo[mom_num][t],dcm[mom_num][t],s1,s2);
      }

    for (int mom_num(0); mom_num < phases.numMom(); ++mom_num){
      for(int t = 0; t < length; ++t){
	// QDPIO::cout<<mom_num<<" "<<t<<" "<<trace(foo[mom_num][t]*Gamma(G5))<<endl;
	int t_eff = (t - t0 + length) % length;
	qqbar[mom_num][t_eff]= foo[mom_num][t] ;
      }
    }
    
    QDPIO::cout<<"Finished the qqbar code\n";
  
    END_CODE();
  }

  /*** bug work around ***/
  //New code that allows the vector mesons
  void compute_qqbar( multi2d<DPropagator>& qqbar,const int gg,
		      const LatticePropagator& quark_prop_1,
		      const LatticePropagator& quark_prop_2, 
		      const SftMom& phases,
		      int t0)
  {
    START_CODE();
  
    QDPIO::cout<<"Starting the qqbar code\n";

    // Length of lattice in decay direction
    Set sft_set(phases.getSet()) ;
    int length(sft_set.numSubsets());
    //QDPIO::cout<<"Time length: "<<length<<endl ;

    // Construct the anti-quark propagator from quark_prop_2
    int G5 = Ns*Ns-1;
    LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

    LatticePropagator quarkloop;
    LatticeSpinMatrix sm ;
    LatticeColorMatrix cm ;

    quarkloop = adj(anti_quark_prop)*Gamma(gg) *quark_prop_1 ;
    LatticeComplex cc;
    multi2d<DPropagator> foo(phases.numMom(),length);
    multi2d<DColorMatrix> dcm(phases.numMom(),length);

    for(int s1(0);s1<Ns;s1++)
      for(int s2(0);s2<Ns;s2++){
	cm = peekSpin(quarkloop,s1,s2);
	for(int c1(0);c1<Nc;c1++)
	  for(int c2(0);c2<Nc;c2++){
	    cc = peekColor(cm,c1,c2);
	    multi2d<DComplex> fcc(phases.sft(cc));
	    for (int mom_num(0); mom_num < phases.numMom(); ++mom_num){
	      for(int t = 0; t < length; ++t){
		pokeColor(dcm[mom_num][t],fcc[mom_num][t],c1,c2);
	      }
	    }
	  }
	for (int mom_num(0); mom_num < phases.numMom(); ++mom_num)
	  for(int t = 0; t < length; ++t)
	    pokeSpin(foo[mom_num][t],dcm[mom_num][t],s1,s2);
      }

    for (int mom_num(0); mom_num < phases.numMom(); ++mom_num){
      for(int t = 0; t < length; ++t){
	// QDPIO::cout<<mom_num<<" "<<t<<" "<<trace(foo[mom_num][t]*Gamma(G5))<<endl;
	int t_eff = (t - t0 + length) % length;
	qqbar[mom_num][t_eff]= foo[mom_num][t] ;
      }
    }
    
    QDPIO::cout<<"Finished the qqbar code with Gamma("<<gg<<")\n";

    END_CODE();
  }



  void write_qqbar(QDPFileWriter& to,
		   multi2d<DPropagator>& qqbar, 
		   const SftMom& phases,
		   string type,
		   string sink){
  
    for(int p(0) ; p<phases.numMom();p++){
      XMLBufferWriter record_xml;
      push(record_xml, "qqbar_desc");//write out the momemtum of each bit
      write(record_xml, "mom", phases.numToMom(p));
      write(record_xml, "type", type);
      write(record_xml, "sink", sink);
      pop(record_xml);

      write(to,record_xml,qqbar[p]);

    }
  
  }

} // namespace Chroma

