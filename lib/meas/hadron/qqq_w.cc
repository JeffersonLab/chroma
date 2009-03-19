//  $Id: qqq_w.cc,v 3.3 2009-03-19 17:17:20 mcneile Exp $
//  $Log: qqq_w.cc,v $
//  Revision 3.3  2009-03-19 17:17:20  mcneile
//  Added guard code so that chroma compiles with Nc <> 3
//
//  Revision 3.2  2008/12/21 21:22:37  edwards
//  Some code cleanups. Put in chroma namespace. Moved around the chromabase.h
//  include.
//
//  Revision 3.1  2008/09/09 20:30:42  kostas
//  added Decuplet to NPLQCD code
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
//  constructs 3 quark propagators contracted at the sink
//


#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "util/ferm/transf.h"
#include "barspinmat_w.h"
#include "qqq_w.h"

namespace Chroma 
{
  //! Baryon-Baryon 2-pt functions (C\gamma_5 diquark)
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   */
  void compute_qqq(multi2d<ThreeQuarks>& qqq, 
		   const LatticePropagator& q1,
		   const LatticePropagator& q2,
		   const LatticePropagator& q3,
		   const SftMom& phases,
		   int t0, int bc_spec)
  {
  
    START_CODE();
    if ( Nc != 3 ){    /* Code is specific to Ns=4 and Nc=3. */
      QDPIO::cerr<<" code only works for Nc=3 and Ns=4\n";
      QDP_abort(111) ;
    }
#if QDP_NC == 3

  
    QDPIO::cout<<"Starting the qqq code\n";

    // Length of lattice in decay direction
    int length  = phases.numSubsets();
    //int j_decay = phases.getDir();
    int num_mom = phases.numMom();
  
    multi2d<LatticeFermion> f1(Ns,Nc) ;

    { //save a little memory usage
      // C gamma_5 = Gamma(5)
      //LatticePropagator c5qc5 = Gamma(5)*q1*Gamma(5) ;
      //c5qc5=-c5qc5 ;
      //I may want to only multiply c5 only in the sink
      //this way the qqq object will be more versatile...
      //for(QuarkIndex i;i.NotEnd();++i)
      //  PropToFerm(c5qc5, f1[i.s][i.c] ,i.c,i.s) ;

      LatticePropagator c5q  = Gamma(5)*q1 ;
      for(QuarkIndex i;i.NotEnd();++i)
	PropToFerm(c5q, f1[i.s][i.c] ,i.c,i.s) ;
    }
  
    multi2d<LatticeFermion> f2(Ns,Nc) ;
    multi2d<LatticeFermion> f3(Ns,Nc) ;
    for(QuarkIndex i;i.NotEnd();++i){
      PropToFerm(q2, f2[i.s][i.c] ,i.c,i.s) ;
      PropToFerm(q3, f3[i.s][i.c] ,i.c,i.s) ;
    }
  
    LatticeComplex cc ;
    multi2d<DComplex> foo(num_mom,length) ;
    for(int s4(0);s4<Ns;s4++) //sink spin index of the 3rd quark 
      for(QuarkIndex i3;i3.NotEnd();++i3)
	for(QuarkIndex i2;i2.NotEnd();++i2)
	  for(QuarkIndex i1;i1.NotEnd();++i1){
	    //cout<<"QuarkIndex iterator: "<<i1.s<<" "<<i1.c<<endl;
	    cc=0.0;
	    // s is the sink spin index of the 1 and 2 quark
	    // that participate in the diquark
	    for(int s(0);s<Ns;s++) //contract the diquark on the sink
	      cc += colorContract(peekSpin(f1[i1.s][i1.c],s ),
				  peekSpin(f2[i2.s][i2.c],s ),
				  peekSpin(f3[i3.s][i3.c],s4));

	    foo = phases.sft(cc); 
	    for(int sink_mom_num(0); sink_mom_num < num_mom; sink_mom_num++) 
	      for(int t = 0; t < length; ++t){
		//shift source to 0 and take care of the antiperiodic BC
		//sign flip for the baryon
		int t_eff = (t - t0 + length) % length;
		qqq[sink_mom_num][t_eff](i1,i2,i3,s4) = 
		  (bc_spec < 0 && (t_eff+t0) >= length) ? -foo[sink_mom_num][t] :
		  foo[sink_mom_num][t] ;
		//cout<<t<<" "<<qqq[sink_mom_num][t](i1,i2,i3,s4)<<endl ;
	      }
	    //cout<<i1.s<<" "<<qqq[0][0](i1,i2,i3,s4)<<endl ;
	  }

    QDPIO::cout<<"Finished the qqq code\n";

#endif
    END_CODE();
    
  }

  void compute_qqq(multi2d<ThreeQuarks>& qqq, const int k, // k = 1 2 3 4(time)
		   const LatticePropagator& q1,
		   const LatticePropagator& q2,
		   const LatticePropagator& q3,
		   const SftMom& phases,
		   int t0, int bc_spec
    ){
  
    START_CODE();
    if ( Nc != 3 ){    /* Code is specific to Ns=4 and Nc=3. */
      QDPIO::cerr<<" code only works for Nc=3 and Ns=4\n";
      QDP_abort(111) ;
    }
#if QDP_NC == 3

  
    QDPIO::cout<<"Starting the qqq code\n";

    // Length of lattice in decay direction
    int length  = phases.numSubsets();
    //int j_decay = phases.getDir();
    int num_mom = phases.numMom();
  
    multi2d<LatticeFermion> f1(Ns,Nc) ;

    { 
      LatticePropagator c5q  = BaryonSpinMats::Cgmu(k)*q1 ;
      for(QuarkIndex i;i.NotEnd();++i)
	PropToFerm(c5q, f1[i.s][i.c] ,i.c,i.s) ;
    }
  
    multi2d<LatticeFermion> f2(Ns,Nc) ;
    multi2d<LatticeFermion> f3(Ns,Nc) ;
    for(QuarkIndex i;i.NotEnd();++i){
      PropToFerm(q2, f2[i.s][i.c] ,i.c,i.s) ;
      PropToFerm(q3, f3[i.s][i.c] ,i.c,i.s) ;
    }
  
    LatticeComplex cc ;
    multi2d<DComplex> foo(num_mom,length) ;
    for(int s4(0);s4<Ns;s4++) //sink spin index of the 3rd quark 
      for(QuarkIndex i3;i3.NotEnd();++i3)
	for(QuarkIndex i2;i2.NotEnd();++i2)
	  for(QuarkIndex i1;i1.NotEnd();++i1){
	    //cout<<"QuarkIndex iterator: "<<i1.s<<" "<<i1.c<<endl;
	    cc=0.0;
	    // s is the sink spin index of the 1 and 2 quark
	    // that participate in the diquark
	    for(int s(0);s<Ns;s++) //contract the diquark on the sink
	      cc += colorContract(peekSpin(f1[i1.s][i1.c],s ),
				  peekSpin(f2[i2.s][i2.c],s ),
				  peekSpin(f3[i3.s][i3.c],s4));

	    foo = phases.sft(cc); 
	    for(int sink_mom_num(0); sink_mom_num < num_mom; sink_mom_num++) 
	      for(int t = 0; t < length; ++t){
		//shift source to 0 and take care of the antiperiodic BC
		//sign flip for the baryon
		int t_eff = (t - t0 + length) % length;
		qqq[sink_mom_num][t_eff](i1,i2,i3,s4) = 
		  (bc_spec < 0 && (t_eff+t0) >= length) ? -foo[sink_mom_num][t] :
		  foo[sink_mom_num][t] ;
		//cout<<t<<" "<<qqq[sink_mom_num][t](i1,i2,i3,s4)<<endl ;
	      }
	    //cout<<i1.s<<" "<<qqq[0][0](i1,i2,i3,s4)<<endl ;
	  }

    QDPIO::cout<<"Finished the qqq code\n";

#endif
    END_CODE();
    
  }

  void write_qqq(QDPFileWriter& to,
		 multi2d<ThreeQuarks>& qqq, 
		 const SftMom& phases,
		 string type,
		 string sink){
    
    for(int p(0) ; p<phases.numMom();p++){
      XMLBufferWriter record_xml;
      push(record_xml, "qqq_desc");//write out the momemtum of each bit
      write(record_xml, "mom", phases.numToMom(p));
      write(record_xml, "type", type);
      write(record_xml, "sink", sink);
      pop(record_xml);

      multi1d<DComplex> data(phases.numSubsets()*qqq[p][0].handle().size());  
      int k(0);
      for(int t(0);t<phases.numSubsets();t++)
	for(int j(0);j<qqq[p][t].Size();j++)
	  data[k++] = qqq[p][t][j] ;
      
      write(to,record_xml,data);
    }

  }

} //Chroma namespace

