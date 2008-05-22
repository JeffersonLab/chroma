// $Id: delta_2pt_w.cc,v 3.9 2008-05-22 18:50:00 caubin Exp $
/*! \file
 *  \brief Construct meson 2pt correlators.
 */

#include "meas/hadron/delta_2pt_w.h"
#include "meas/hadron/barhqlq_w.h"
#include "meas/hadron/hadron_contract_factory.h"
#include "barspinmat_w.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, Delta2PtEnv::Params& param)
  {
    Delta2PtEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const Delta2PtEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Delta correlators
  /*! 
   * \ingroup hadron 
   *
   * @{
   */
  namespace Delta2PtEnv
  { 
    //! Anonymous namespace
    namespace
    {

      //-------------------- callback functions ---------------------------

      //! Construct pion correlator
      HadronContract* mesDeltaCorrs(XMLReader& xml_in,
					const std::string& path)
      {
	return new DeltaCorrs(Params(xml_in, path));   // all gammas
      }


      //! Local registration flag
      bool registered = false;

    } // end anonymous namespace


    //! Initialize
    Params::Params()
    {
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "mom2_max", mom2_max);
      read(paramtop, "avg_equiv_mom", avg_equiv_mom);
      read(paramtop, "mom_origin", mom_origin);
      if(paramtop.count("min_contractions") !=0 )
	read(paramtop, "min_contractions", min_contractions);
      else
	min_contractions = false ;
      if(paramtop.count("parity") !=0 ) 
	read(paramtop, "parity", parity);
      else
	parity = "all" ;//By default, do both pos and neg parity.
      read(paramtop, "first_id", first_id);
      read(paramtop, "second_id", second_id);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "mom2_max", mom2_max);
      write(xml, "avg_equiv_mom", avg_equiv_mom);
      write(xml, "mom_origin", mom_origin);
      write(xml, "min_contractions", min_contractions);
      write(xml, "parity", parity);
      write(xml, "first_id", first_id);
      write(xml, "second_id", second_id);

      pop(xml);
    }


    // Construct all the correlators
    std::list< Handle<HadronContractResult_t> >
    DeltaCorrs::operator()(const multi1d<LatticeColorMatrix>& u,
				    const std::string& xml_group,
				    const std::string& id_tag)
    {
      START_CODE();

      QDPIO::cout << "Hadron2Pt: Delta" << endl;

      multi1d<ForwardProp_t> forward_headers(2);
      forward_headers[0] = readForwardPropHeader(params.first_id);
      forward_headers[1] = readForwardPropHeader(params.second_id);
      
      multi1d<int> t_srce = getTSrce(forward_headers);
      int decay_dir       = getDecayDir(forward_headers);

      // Get references to the props
      const LatticePropagator& quark_prop1 = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.first_id);
      const LatticePropagator& quark_prop2 = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.second_id);

      // Parameters needed for the momentum projection
      SftMomParams_t sft_params;
      sft_params.mom2_max      = params.mom2_max;
      sft_params.origin_offset = t_srce;
      sft_params.mom_offset    = params.mom_origin;
      sft_params.avg_equiv_mom = params.avg_equiv_mom;
      sft_params.decay_dir     = decay_dir;

      std::list< Handle<Hadron2PtContract_t> > hadron;   // holds the contract lattice correlator
      std::map<std::string,SpinMatrix> Projector ;
      std::map<std::string,SpinMatrix> Parity ;
      multi1d<SpinMatrix> SrcDiQuark(Ns) ;
      multi1d<SpinMatrix> SnkDiQuark(Ns) ;
      Parity["PosPar"] = BaryonSpinMats::Tunpol();
      Parity["NegPar"] = BaryonSpinMats::TunpolNegPar();
      Projector["OnePlusSigma3"] = BaryonSpinMats::TspinUp() ;
      Projector["OneMinusSigma3"] = BaryonSpinMats::TspinDown() ;
      Projector["SigmaPlus"] = BaryonSpinMats::T_ig5XpiY();
      Projector["SigmaMinus"] = BaryonSpinMats::T_ig5XmiY();

	  for(int mu(0);mu<Ns ;mu++)
		SnkDiQuark[mu] = BaryonSpinMats::Cgmu(mu+1);
	  for(int mu(0);mu<Ns ;mu++)
		SrcDiQuark[mu] = BaryonSpinMats::CgmuTrans(mu+1);
	  
	  if(!params.min_contractions){
		/* 
		   If we do all contractions, we ignore the parity parameter.
		 */
		map<std::string,SpinMatrix>::iterator par;
		map<std::string,SpinMatrix>::iterator proj;
		for ( par=Parity.begin();par != Parity.end(); par++){
		  for ( proj=Projector.begin();proj != Projector.end(); proj++){
			SpinMatrix T = par->second * proj->second ;// the projector matrix
			QDPIO::cout<<" Parity: "<<par->first<<endl;
			QDPIO::cout<<" Projector: "<<proj->first<<endl;
			for( int src(0) ;src<Ns;src++)
			  for( int snk(0) ;snk<Ns;snk++){
				QDPIO::cout<<"   Computing C_"<<snk<<src<<endl;
				
				Handle<Hadron2PtContract_t> had(new Hadron2PtContract_t);
				had->corr = Baryon2PtContractions::sigmast2pt(quark_prop1, 
															  quark_prop2,
															  T,SrcDiQuark[src],
															  SnkDiQuark[snk]);
				push(had->xml, xml_group);
				write(had->xml, id_tag, "delta");
				write(had->xml, "SrcDiQuark", src);
				write(had->xml, "SnkDiQuark", snk);
				write(had->xml, "Parity", par->first);
				write(had->xml, "Projector", proj->first);	      
				write(had->xml, "PropHeaders", forward_headers);
				
				pop(had->xml);
				
				hadron.push_back(had); 
			  }
		  }
		}
		
		//projector for the spin averaged correlator
		//only works for zero momentum
		SpinMatrix g_one = 1.0 ;
		multi1d< multi1d<SpinMatrix> > ProjGmuGnu(Ns-1) ;
		for(int s1(0);s1<Ns-1;s1++){
		  ProjGmuGnu[s1].resize(Ns-1) ;
		  for(int s2(0);s2<Ns-1;s2++)
			ProjGmuGnu[s1][s2] = Gamma(1<<s1) * (Gamma(1<<s2)*g_one) ;
		}
		
		for ( par=Parity.begin();par != Parity.end(); par++){
		  Handle<Hadron2PtContract_t> had(new Hadron2PtContract_t);
		  had->corr = 0.0 ;
		  for( int src(0) ;src<Ns-1;src++)
			for( int snk(0) ;snk<Ns-1;snk++){
			  QDPIO::cout<<"   Computing C_"<<snk<<src<<endl;
			  SpinMatrix T =  (- 1.0/3.0) * ProjGmuGnu[src][snk] ;
			  if(src == snk ) T += g_one ;
			  T = par->second * T ;
			  had->corr += Baryon2PtContractions::sigmast2pt(quark_prop1, 
															 quark_prop2,
															 T,SrcDiQuark[src],
															 SnkDiQuark[snk]);
			}
		  push(had->xml, xml_group);
		  write(had->xml, id_tag, "delta");
		  write(had->xml, "Parity", par->first);
		  write(had->xml, "Projector", "SpinAveraged");
		  write(had->xml, "PropHeaders", forward_headers);
		  pop(had->xml);
	      
		  hadron.push_back(had); 
		}
	  }//Ends if(!params.min_contractions)

	  if(params.min_contractions && (params.parity!="all")){
		/**
		   Here we need all mu, nu combos for the 1 +/- Sigma3,and
		   just the 02, 20, 12, 21 for the SigmaMinus and SigmaPlus
		   This only does one of the parities, either positive or negative,
		   but default is positive
		 **/
		map<std::string,SpinMatrix>::iterator proj;
		SpinMatrix par;
		/**
		   if(params.parity=="Neg")
		   par = Parity["NegPar"];
		   else
		   par = Parity["PosPar"];
		**/
		par = Parity[params.parity];
		cout<<"Parity flag is"<<params.parity<<endl;
		for ( proj=Projector.begin();proj != Projector.end(); proj++){
		  SpinMatrix T = par * proj->second ;// the projector matrix
		  QDPIO::cout<<" Parity: "<<params.parity<<endl;
		  QDPIO::cout<<" Projector: "<<proj->first<<endl;
		  if((proj->first)=="SigmaPlus"||(proj->first)=="SigmaMinus"){
		  for( int src(0) ;src<Ns-1;src++)
			for( int snk(0) ;snk<Ns-1;snk++){
			  if((src!=snk) && ((src==2)||(snk==2)) ){
				QDPIO::cout<<"   Computing C_"<<snk<<src<<endl;
				
				Handle<Hadron2PtContract_t> had(new Hadron2PtContract_t);
				had->corr = Baryon2PtContractions::sigmast2pt(quark_prop1, 
															  quark_prop2,
															  T,SrcDiQuark[src],
															  SnkDiQuark[snk]);
				
				
				push(had->xml, xml_group);
				write(had->xml, id_tag, "delta");
				write(had->xml, "SrcDiQuark", src);
				write(had->xml, "SnkDiQuark", snk);
				write(had->xml, "Parity", "PosPar");
				write(had->xml, "Projector", proj->first);	      
				write(had->xml, "PropHeaders", forward_headers);
			  
				pop(had->xml);
			  
				hadron.push_back(had); 
			  }
			}
		  }
		  else{
			for( int src(0) ;src<Ns-1;src++)
			  for( int snk(0) ;snk<Ns-1;snk++){
				QDPIO::cout<<"   Computing C_"<<snk<<src<<endl;
				
				Handle<Hadron2PtContract_t> had(new Hadron2PtContract_t);
				had->corr = Baryon2PtContractions::sigmast2pt(quark_prop1, 
															  quark_prop2,
															  T,SrcDiQuark[src],
															  SnkDiQuark[snk]);
				
				
				push(had->xml, xml_group);
				write(had->xml, id_tag, "delta");
				write(had->xml, "SrcDiQuark", src);
				write(had->xml, "SnkDiQuark", snk);
				write(had->xml, "Parity", "PosPar");
				write(had->xml, "Projector", proj->first);	      
				write(had->xml, "PropHeaders", forward_headers);
				
				pop(had->xml);
				
				hadron.push_back(had); 
			  }
		  }
		}
		
	  }//Ends if(params.min_contractions)

	  if(params.min_contractions && (params.parity=="all")){
		/**
		   Here we need all mu, nu combos for the 1 +/- Sigma3,and
		   just the 02, 20, 12, 21 for the SigmaMinus and SigmaPlus
		   This does both positive and negative parities.
		 **/
		map<std::string,SpinMatrix>::iterator par;
		map<std::string,SpinMatrix>::iterator proj;
		for ( par=Parity.begin();par != Parity.end(); par++){		
		  for ( proj=Projector.begin();proj != Projector.end(); proj++){
			SpinMatrix T = par->second * proj->second ;// the projector matrix
			QDPIO::cout<<" Parity: "<<par->first<<endl;
			QDPIO::cout<<" Projector: "<<proj->first<<endl;
			if((proj->first)=="SigmaPlus"||(proj->first)=="SigmaMinus"){
			  for( int src(0) ;src<Ns-1;src++)
				for( int snk(0) ;snk<Ns-1;snk++){
				  if((src!=snk) && ((src==2)||(snk==2)) ){
					QDPIO::cout<<"   Computing C_"<<snk<<src<<endl;
					
					Handle<Hadron2PtContract_t> had(new Hadron2PtContract_t);
					had->corr = Baryon2PtContractions::sigmast2pt(quark_prop1, 
																  quark_prop2,
																  T,SrcDiQuark[src],
																  SnkDiQuark[snk]);
					
					
					push(had->xml, xml_group);
					write(had->xml, id_tag, "delta");
					write(had->xml, "SrcDiQuark", src);
					write(had->xml, "SnkDiQuark", snk);
					write(had->xml, "Parity", par->first);
					write(had->xml, "Projector", proj->first);	      
					write(had->xml, "PropHeaders", forward_headers);
					pop(had->xml);
					
					hadron.push_back(had); 
				  }
				}
			}
			else{
			  for( int src(0) ;src<Ns-1;src++)
				for( int snk(0) ;snk<Ns-1;snk++){
				  QDPIO::cout<<"   Computing C_"<<snk<<src<<endl;
				  
				  Handle<Hadron2PtContract_t> had(new Hadron2PtContract_t);
				  had->corr = Baryon2PtContractions::sigmast2pt(quark_prop1, 
																quark_prop2,
																T,SrcDiQuark[src],
																SnkDiQuark[snk]);
				  
				  
				  push(had->xml, xml_group);
				  write(had->xml, id_tag, "delta");
				  write(had->xml, "SrcDiQuark", src);
				  write(had->xml, "SnkDiQuark", snk);
				  write(had->xml, "Parity", par->first);
				  write(had->xml, "Projector", proj->first);	      
				  write(had->xml, "PropHeaders", forward_headers);
				  
				  pop(had->xml);
				  
				  hadron.push_back(had); 
				}
			}
		  }
		}
	  }//Ends if(params.min_contractions)
	  
      
      END_CODE();

      return this->project(hadron, sft_params);
    }


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//! Register all the factories
	success &= Chroma::TheHadronContractFactory::Instance().registerObject(string("Delta"), mesDeltaCorrs);

	registered = true;
      }
      return success;
    }

  }  // end namespace Delta2PtEnv

  /*! @} */   // end of group io

}  // end namespace Chroma


  
