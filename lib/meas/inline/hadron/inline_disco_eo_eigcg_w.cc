
// $Id: inline_disco_eo_eigcg_w.cc,v 3.3 2009-07-31 14:09:31 bjoo Exp $
/*! \file
 * \brief Inline measurement 3pt_prop
 *
 * This uses eo-preconditioning on the noise vector side as well as
 * with the eig-cg vectors. NOTE THIS CODE DOES NOT WORK.          
 */
#include <vector> 
#include <map> 
#include <qdp-lapack.h>
#include <qdp_config.h>

#include "handle.h"
#include "meas/inline/hadron/inline_disco_eo_eigcg_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"
#include "meas/smear/displacement.h"
#include "meas/sources/source_smearing_aggregate.h"
#include "meas/sources/source_smearing_factory.h"
#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "meas/sinks/sink_smearing_factory.h"
#include "meas/hadron/barspinmat_w.h"
#include "meas/hadron/baryon_operator_aggregate_w.h"
#include "meas/hadron/baryon_operator_factory_w.h"
#include "meas/hadron/barQll_w.h"
#include "meas/hadron/dilution_scheme_aggregate.h"
#include "meas/hadron/dilution_scheme_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "util/info/unique_id.h"
#include "util/ferm/transf.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "eoprec_logdet_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"

#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/fermacts/wilson_fermact_params_w.h"
#include "actions/ferm/linop/linop_w.h"

#include "util/ferm/key_val_db.h"

namespace Chroma{
  namespace InlineDiscoEoEigCGEnv{ 
    namespace{
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "DISCO_EO_EIGCG";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//success &= BaryonOperatorEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  
  // Reader for input parameters
  void read(XMLReader& xml, const string& path, Params::Param_t& param){
      XMLReader paramtop(xml, path);
      
      int version;
      read(paramtop, "version", version);
      
      switch (version) 
	{
	case 1:
	  /************************************************************/
	  read(paramtop,"max_path_length",param.max_path_length);
	  read(paramtop,"p2_max",param.p2_max);
	  read(paramtop,"mass_label",param.mass_label);
	  param.chi = readXMLArrayGroup(paramtop, "Quarks", "DilutionType");
	  param.action = readXMLGroup(paramtop, "FermionAction","FermAct");
	  
	  break;
	  
	default :
	  /**************************************************************/

	  QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	  QDP_abort(1);
	}
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params::Param_t& param){
      push(xml, path);
      
      int version = 1;
      
      write(xml, "version", version);

      write(xml,"max_path_length",param.max_path_length);
      write(xml,"p2_max",param.p2_max);
      write(xml,"mass_label",param.mass_label);
      push(xml,"FermionAction");
      xml<<param.action.xml ;
      pop(xml);
      push(xml,"Quarks");
      for( int t(0);t<param.chi.size();t++){
	push(xml,"elem");
	xml<<param.chi[t].xml;
	pop(xml);
      }
      pop(xml);
      
      pop(xml); // final pop
    }
    
    //! Gauge field parameters
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);
      
      read(inputtop, "gauge_id"   , input.gauge_id   ) ;
      // If no file is included, then we will turn off the eig_cg part (ie, use 0 vectors)
      if ( inputtop.count("evecs_file")!=0 )
	read(inputtop, "evecs_file" , input.evecs_file ) ;
      else
	input.evecs_file="";
	
      read(inputtop, "op_db_file" , input.op_db_file ) ;
    }
    
    //! Gauge field parameters
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input){
      push(xml, path);
      
      write(xml, "gauge_id"   , input.gauge_id   );
      write(xml, "evecs_file" , input.evecs_file );
      write(xml, "op_db_file" , input.op_db_file );
      pop(xml);
    }
    
    
    // Param stuff
    Params::Params(){ 
      frequency = 0;
    }
    
    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
	{
	  XMLReader paramtop(xml_in, path);
	  
	  if (paramtop.count("Frequency") == 1)
	    read(paramtop, "Frequency", frequency);
	  else
	    frequency = 1;
	  
	  // Read program parameters
	  read(paramtop, "Param", param);
	  
	  // Read in the output propagator/source configuration info
	  read(paramtop, "NamedObject", named_obj);
	  
	  // Possible alternate XML file pattern
	  if (paramtop.count("xml_file") != 0) 
	    {
	      read(paramtop, "xml_file", xml_file);
	    }
	}
      catch(const std::string& e) 
	{
	  QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
	  QDP_abort(1);
	}
    }


    void Params::write(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
      
      // Parameters for source construction
      InlineDiscoEoEigCGEnv::write(xml_out, "Param", param);
      
      // Write out the output propagator/source configuration info
      InlineDiscoEoEigCGEnv::write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }

    
    //! Meson operator     
    struct KeyOperator_t
    {
      multi1d<short int> mom     ; /*!< D-1 momentum of this operator */
      unsigned short int t_slice ; /*!< Meson operator time slice */
      multi1d<short int> disp    ; /*!< Displacement dirs of quark (right)*/
      
      KeyOperator_t(){
	mom.resize(Nd-1);
      }
    };
    
    bool operator<(const KeyOperator_t& a, const KeyOperator_t& b){
      return ((a.t_slice<b.t_slice)||(a.mom<b.mom)||(a.disp<b.disp));
    }

    std::ostream& operator<<(std::ostream& os, const KeyOperator_t& d)
    {
      os << "KeyOperator_t:"
	 << " t_slice = " << d.t_slice
	 << ", disp = ";
      for (int i=0; i<d.disp.size();i++)
	os << d.disp[i] << " " ;
      
      os << ", mom = ";
      for (int i=0; i<d.mom.size();i++)
	os << d.mom[i] << " " ;

      os << std::endl;

      return os;
    }

    class ValOperator_t{
    public:
      multi1d<ComplexD> op ;  
      ValOperator_t(){op.resize(Ns*Ns);} // Here go the 16 gamma matrices
      ~ValOperator_t(){}
    } ;

    //-------------------------------------------------------------------------                            
    //! stream IO                                                                                          
    std::ostream& operator<<(std::ostream& os, const ValOperator_t& d)
    {
      os << "ValOperator_t:\n";
      for (int i=0; i<d.op.size();i++)
        os <<"     gamma["<<i<<"] = "<< d.op[i] << endl ;
      
      return os;
    }
    
    struct KeyVal_t{
      SerialDBKey <KeyOperator_t> k ;
      SerialDBData<ValOperator_t> v ;
    };

    //! KeyOperator reader    
    void read(BinaryReader& bin, KeyOperator_t& d){
      read(bin,d.t_slice);
      unsigned short int n ;
      read(bin,n);
      d.disp.resize(n); 
      read(bin,d.disp);
      d.mom.resize(Nd-1) ;
      read(bin,d.mom);
    }
    //! KeyOperator writer
    void write(BinaryWriter& bin, const KeyOperator_t& d){
      write(bin,d.t_slice);
      unsigned short int n ;
      n = d.disp.size();
      write(bin,n);
      write(bin,d.disp);
      write(bin,d.mom);
    }

    //! ValOperator reader    
    void read(BinaryReader& bin, ValOperator_t& d){
      d.op.resize(Ns*Ns);
      read(bin,d.op);
    }
    //! ValOperator writer
    void write(BinaryWriter& bin, const ValOperator_t& d){
      write(bin,d.op);
    }

    namespace{
      StandardOutputStream& operator<<(StandardOutputStream& os, const multi1d<short int>& d){
	if (d.size() > 0){
	  os << d[0]; 
	  for(int i=1; i < d.size(); i++)
	    os << " " << d[i];
	}
	return os;
      }
    }

    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    Handle<EvenOddPrecLinearOperator<T,P,Q> > 
    createOddOdd_Op( const Params::Param_t& param, const P& u){     
      //
      // Initialize fermion action
      //
      std::istringstream  xml_s(param.action.xml);
      XMLReader  fermacttop(xml_s);
      QDPIO::cout << "FermAct = " << param.action.id << endl;
      //
      // Try the factories
      //
      Handle< FermState< T,P,Q> > state ;
      try{
	  QDPIO::cout << "Try the various Wilson fermion factories" << endl;
	  // Generic Wilson-Type stuff
	  Handle< FermionAction< T,P,Q > >
	    Sf(TheFermionActionFactory::Instance().createObject(param.action.id,
								fermacttop,
								param.action.path));
	  state = Sf->createState(u);//got the state
	  QDPIO::cout << "Suitable factory found: compute the trace quark prop"<<endl;
      }
      catch (const std::string& e){
	QDPIO::cout << name 
		    << ": caught exception instantiating the action: " << e << endl;
      }
      
      // Now need to construct the operator
      // this seems tricky and maybe there is a better way... (Robert-Balint help!)
      // We only work with Wilson and Clover Wilson fermions so check the following 2
      // For the moment we restrict to EO preconditioned only
      
      //this is the Odd-Odd piece
      Handle<EvenOddPrecLinearOperator<T,P,Q> > Doo ;
      if(param.action.id == "WILSON"){
	WilsonFermActParams wp(fermacttop,param.action.path);
	return new EvenOddPrecWilsonLinOp(state,wp.Mass,wp.anisoParam ) ;
      }
      else if(param.action.id == "CLOVER"){
	CloverFermActParams cp(fermacttop,param.action.path);
	return new EvenOddPrecCloverLinOp(state,cp) ;
      }
      else{
	QDPIO::cout<<name<<" : Tough luck dude! No code for you..."<<endl ;
	QDP_abort(1);
      }
      return NULL ;
    }

    struct CholeskyFactors{
      multi1d<Real>    evals ;
      multi1d<Complex> H     ;
      multi1d<Complex> HU    ;
      int              ldh   ;
      int              Nvec  ;
    } ;

    void do_disco(map< KeyOperator_t, ValOperator_t >& db,
		  const LatticeFermion& qbar,
		  const LatticeFermion& q,
		  const SftMom& p,
		  const int& t, 
		  const Subset& trb,
		  const multi1d<short int>& path,
	 	  const int& max_path_length ){
      QDPIO::cout<<" Computing Operator with path length "<<path.size()
		 <<" on timeslice "<<t<<".   Path: "<<path <<endl;
      /**
	 This do_disco routine should not be called, I don't think it's ever needed...
      
       **/
      ValOperator_t val ;
      KeyOperator_t key ;
      pair<KeyOperator_t, ValOperator_t> kv ; 
      kv.first.t_slice = t ;
      if(path.size()==0){
	kv.first.disp.resize(1);
	kv.first.disp[0] = 0 ;
      }
      else
	kv.first.disp = path ;
      
      multi1d< multi1d<ComplexD> > foo(p.numMom()) ;
      for (int m(0); m < p.numMom(); m++)
	foo[m].resize(Ns*Ns);
      for(int g(0);g<Ns*Ns;g++){
	LatticeComplex cc = localInnerProduct(qbar,Gamma(g)*q);
	for (int m(0); m < p.numMom(); m++){
	  // trb is the set of even/odd sites on timeslice t
	  foo[m][g] = sum(p[m]*cc,trb);
	}
      }
      for (int m(0); m < p.numMom(); m++){
	for(int i(0);i<(Nd-1);i++)
	  kv.first.mom[i] = p.numToMom(m)[i] ;
	
	kv.second.op = foo[m];
	pair<map< KeyOperator_t, ValOperator_t >::iterator, bool> itbo;
	
	itbo = db.insert(kv);
        if( itbo.second ){
	  QDPIO::cout<<"Inserting new entry in map\n";
        }
        else{ // if insert fails, key already exists, so add result
          for(int i(0);i<kv.second.op.size();i++)
	    itbo.first->second.op[i] += kv.second.op[i] ;
        }
      }//m 
      
      if(path.size()<max_path_length){
	QDPIO::cout<<" attempt to add new path. "
		   <<" current path length is : "<<path.size();
	multi1d<short int> new_path(path.size()+1);
	QDPIO::cout<<" new path length is : "<<new_path.size()<<endl;
	for(int i(0);i<path.size();i++)
	  new_path[i] = path[i] ;
	for(int sign(-1);sign<2;sign+=2)
	  for(int mu(0);mu<Nd;mu++){
	    new_path[path.size()]= sign*(mu+1) ;
	    //skip back tracking 
	    bool back_track=false ;
	    if(path.size()>0)
	      if(path[path.size()-1] == -new_path[path.size()])
		back_track=true;
	    if(!back_track){
	      QDPIO::cout<<" Added path: "<<new_path<<endl;
	      LatticeFermion q_mu ;
	      if(sign>0)
		q_mu = shift(q, FORWARD, mu);
	      else
		q_mu = shift(q, BACKWARD, mu);

	      do_disco(db, qbar, q_mu, p, t, trb,  new_path, max_path_length);
	    } // skip backtracking
	  } // mu
      }// path.size loop
      
    }// do_disco

    void do_disco(map< KeyOperator_t, ValOperator_t >& db,
		  const LatticeFermion& qbar,
		  const LatticeFermion& q,
		  const Handle<EvenOddPrecLinearOperator<T,P,Q> >& Doo,
		  const SftMom& p,
		  const int& t, 
		  const Subset& trb,
		  const multi1d<short int>& path,
	 	  const int& max_path_length ){
      QDPIO::cout<<" Computing Operator with path length "<<path.size()
		 <<" on timeslice "<<t<<".   Path: "<<path <<endl;
      
      /** This version of do_disco is the one we call 
	  
      Currently we have editted this so as to put the timeslice and phase factor
      on the correct part for the "even" sum, so right now this code is not terribly 
      efficient. Work on optimizing this all later.      

      **/
      
      ValOperator_t val ;
      KeyOperator_t key ;
      pair<KeyOperator_t, ValOperator_t> kv ; 
      kv.first.t_slice = t ;
      if(path.size()==0){
	kv.first.disp.resize(1);
	kv.first.disp[0] = 0 ;
      }
      else
	kv.first.disp = path ;

      multi1d< multi1d<ComplexD> > foo(p.numMom()) ;
      for (int m(0); m < p.numMom(); m++)
        foo[m].resize(Ns*Ns);

      // First do the odd sites, and put in foo.
      for(int g(0);g<Ns*Ns;g++){
        LatticeComplex cc = localInnerProduct(qbar,Gamma(g)*q);
        for (int m(0); m < p.numMom(); m++){
          // trb is the set of odd sites on timeslice t
          foo[m][g] = sum(p[m]*cc,trb);
        }
      }
      
      Set t0set;
      t0set.make(TimeSliceFunc(Nd-1));// 3 is decay dir....need to change this later.
      
      // Now we have to do the even sites
      LatticeFermion qtmp  = zero ;
      LatticeFermion qbar2 = zero ;
      Doo->evenOddLinOp(qtmp,q,PLUS);
      Doo->evenEvenInvLinOp(qbar2,qtmp,PLUS);                                            
      for(int g(0);g<Ns*Ns;g++){
        for (int m(0); m < p.numMom(); m++){
	  LatticeFermion q2 = zero ;
	  LatticeFermion chi = zero ;
	  qtmp = zero ;
	  q2[t0set[t]] = p[m]*(Gamma(g)*qbar2);//rhs should only set the t0 terms non-zero...
	  Doo->evenEvenInvLinOp(qtmp,q2,PLUS);
	  Doo->oddEvenLinOp(chi,qtmp,PLUS);                                                   
	  foo[m][g] += innerProduct(qbar,chi);
	}
      }
      
      for (int m(0); m < p.numMom(); m++){
	for(int i(0);i<(Nd-1);i++)
	  kv.first.mom[i] = p.numToMom(m)[i] ;
	
	kv.second.op = foo[m];
	pair<map< KeyOperator_t, ValOperator_t >::iterator, bool> itbo;
	
	itbo = db.insert(kv);
        if( itbo.second ){
	  QDPIO::cout<<"Inserting new entry in map\n";
        }
        else{ // if insert fails, key already exists, so add result
          for(int i(0);i<kv.second.op.size();i++)
	    itbo.first->second.op[i] += kv.second.op[i] ;
        }
      }//m 
      
      if(path.size()<max_path_length){
	QDPIO::cout<<" attempt to add new path. "
		   <<" current path length is : "<<path.size();
	multi1d<short int> new_path(path.size()+1);
	QDPIO::cout<<" new path length is : "<<new_path.size()<<endl;
	for(int i(0);i<path.size();i++)
	  new_path[i] = path[i] ;
	for(int sign(-1);sign<2;sign+=2)
	  for(int mu(0);mu<Nd;mu++){
	    new_path[path.size()]= sign*(mu+1) ;
	    //skip back tracking 
	    bool back_track=false ;
	    if(path.size()>0)
	      if(path[path.size()-1] == -new_path[path.size()])
		back_track=true;
	    if(!back_track){
	      QDPIO::cout<<" Added path: "<<new_path<<endl;
	      LatticeFermion q_mu ;
	      if(sign>0)
		q_mu = shift(q, FORWARD, mu);
	      else
		q_mu = shift(q, BACKWARD, mu);

	      do_disco(db, qbar, q_mu, Doo, p, t, trb,  new_path, max_path_length);
	    } // skip backtracking
	  } // mu
      }// path.size loop
      
    }// do_disco

    void do_disco(map< KeyOperator_t, ValOperator_t >& db,
                  const Params::Param_t& param,
                  const P& u,
		  const SftMom& p,
		  const int& t, 
		  const Subset& trb,
		  const multi1d<short int>& path){
      QDPIO::cout<<" Computing Operator with path length "<<path.size()
		 <<" on timeslice "<<t<<".   Path: "<<path <<endl;
      
      int max_path_length = param.max_path_length;
      ValOperator_t val ;
      KeyOperator_t key ;
      pair<KeyOperator_t, ValOperator_t> kv ; 
      kv.first.t_slice = t ;
      if(path.size()==0){
	kv.first.disp.resize(1);
	kv.first.disp[0] = 0 ;
      }
      else
	kv.first.disp = path ;
      
      Handle<EvenOddPrecLinearOperator<T,P,Q> > Doo = createOddOdd_Op(param,u);
      
      multi1d< multi1d<ComplexD> > foo(p.numMom()) ;
      for (int m(0); m < p.numMom(); m++){
	foo[m].resize(Ns*Ns);
	for (int g(0); g<Ns*Ns;g++){
	  foo[m][g] = 0.0;
	}
      }
      for(int g(0);g<Ns*Ns;g++){
	for(int col=0;col<3;col++){
	  for(int sp=0;sp<4;sp++){
	    Fermion tt = zero ;
	    ColorVector cv = zero ;
	    Complex z = cmplx(Real(1.0),0.0) ;
	    pokeColor(cv,z,col);
	    pokeSpin(tt,cv,sp);
	    LatticeFermion V = tt ;
	    for (int m(0); m < p.numMom(); m++){
	      //only on even sites
	      LatticeFermion DV = zero;
	      Doo->evenEvenInvLinOp(DV,V,PLUS);
	      foo[m][g] += sum(p[m]*localInnerProduct(V,Gamma(g)*DV),trb);
	    }//m
	  }//spin
	}//col
      }//g

      for (int m(0); m < p.numMom(); m++){
	for(int i(0);i<(Nd-1);i++)
	  kv.first.mom[i] = p.numToMom(m)[i] ;
	
	kv.second.op = foo[m];
	pair<map< KeyOperator_t, ValOperator_t >::iterator, bool> itbo;
	
	itbo = db.insert(kv);
        if( itbo.second ){
	  QDPIO::cout<<"Inserting new entry in map\n";
        }
        else{ // if insert fails, key already exists, so add result
          for(int i(0);i<kv.second.op.size();i++)
            itbo.first->second.op[i] += kv.second.op[i] ;
        }
      }//m 
      
    }// do_disco

    void do_disco(map< KeyOperator_t, ValOperator_t >& db,
		  CholeskyFactors Clsk , 
		  multi1d<LatticeFermion>& vec,
		  const Params::Param_t& param, 
                  const Handle<EvenOddPrecLinearOperator<T,P,Q> >& Doo,
		  const P& u,
                  const int& t,
		  const Subset& trb, 
		  const SftMom& p,
                  const multi1d<short int>& path){

      QDPIO::cout<<" Computing Operator with path length "<<path.size()
		 <<" on timeslice "<<t<<".   Path: "<<path <<endl;
      int max_path_length = param.max_path_length;
      
      ValOperator_t val ;
      KeyOperator_t key ;
      pair<KeyOperator_t, ValOperator_t> kv ; 

      kv.first.t_slice = t ;
      if(path.size()==0){
        kv.first.disp.resize(1);
        kv.first.disp[0] = 0 ;
      }
      else
        kv.first.disp = path ;
      
      int ldb = vec.size();
      int info;
      char U = 'U';
      int Nrhs = ldb; // Because we have no dilution vectors, but rhs's are made of EigCG vecs...
      multi2d<Complex> B(Nrhs,ldb);
      multi1d< multi1d<LatticeFermion> > DDvec(ldb);
      // Zero out the fermions
      for(int i(0); i<ldb;i++){
	DDvec[i].resize(Ns*Ns);
	DDvec[i] = zero;
      }
      cout<<"We've initialized DDvec!\n";

      for(int i(0); i<ldb;i++){
	LatticeFermion qtmp  = zero ;
	LatticeFermion qtmp2 = zero ;
	LatticeFermion qtmp3 = zero ;
	Doo->evenOddLinOp(qtmp,vec[i],PLUS);
	Doo->evenEvenInvLinOp(qtmp2,qtmp,PLUS);
	for(int g(0);g<Ns*Ns;g++){
	  Doo->evenEvenInvLinOp(qtmp3,Gamma(g)*qtmp2,PLUS);
	  qtmp = zero;
	  qtmp2 = zero;
	  Doo->oddEvenLinOp(qtmp,qtmp3,PLUS);
	  Doo->oddOddLinOp(qtmp2,qtmp,MINUS);
	  qtmp = zero;
	  Doo->oddOddLinOp(qtmp,Gamma(g)*vec[i],MINUS);
	  // I think this could run into memory problems with too many vectors
	  DDvec[i][g] = qtmp + qtmp2;
	}
      }

      multi1d< multi1d<ComplexD> > foo(p.numMom()) ;
      for (int m(0); m < p.numMom(); m++){
	foo[m].resize(Ns*Ns);
	foo[m] = zero;
      }
      // Okay, this is probably inefficient, b/c we are doing some things multiple times...
      for(int g(0);g<Ns*Ns;g++){
	for (int m(0); m < p.numMom(); m++){
	  for (int i(0); i<ldb;i++){
	    for (int j(0); j<ldb;j++){
	      //	      LatticeComplex cc1 = localInnerProduct(Svec[j],Gamma(g)*vec[i]);
	      LatticeComplex cc2 = localInnerProduct(vec[j],DDvec[i][g]);
	      B[i][j] = sum(p[m]*cc2,trb) ;
	    }//j
	  }//i 

#if BASE_PRECISION == 32
	  int r = QDPLapack::cpotrs(U, Clsk.Nvec, Nrhs, Clsk.HU, Clsk.ldh, B, ldb, info);
#else
	  int r = QDPLapack::zpotrs(U, Clsk.Nvec, Nrhs, Clsk.HU, Clsk.ldh, B, ldb, info);
#endif
	  // and at this point, B is H^-1 Vdag Sdag gamma V, so
	  // For debugging purposes, both r and info should be zero...
	  QDPIO::cout<<"do_disco cpotrs r = "<<r<<endl;
	  QDPIO::cout<<"do_disco cpotrs info = "<<info<<endl;

	  //	  foo[m][g] = 0.0;
	  for (int i(0); i<ldb;i++){
	    foo[m][g] += ComplexD(B[i][i]);// Trace over last set of indices
	  }
	}//m
      }//g
      for (int m(0); m < p.numMom(); m++){
	for(int i(0);i<(Nd-1);i++)
	  kv.first.mom[i] = p.numToMom(m)[i] ;
	
	kv.second.op = foo[m];

	pair<map< KeyOperator_t, ValOperator_t >::iterator, bool> itbo;
	
	itbo = db.insert(kv);
	if( !(itbo.second) ){  // if insert fails, key already exists, so add result
	  for(int i(0);i<kv.second.op.size();i++)
	    itbo.first->second.op[i] += kv.second.op[i] ;
	}
	else
	  QDPIO::cout<<"Inserting new entry in map\n";
      }

      if(path.size()<max_path_length){
	QDPIO::cout<<" attempt to add new path. "
		   <<" current path length is : "<<path.size();
	multi1d<short int> new_path(path.size()+1);
	QDPIO::cout<<" new path length is : "<<new_path.size()<<endl;
	for(int i(0);i<path.size();i++)
	  new_path[i] = path[i] ;
	for(int sign(-1);sign<2;sign+=2)
	  for(int mu(0);mu<Nd;mu++){
	    new_path[path.size()]= sign*(mu+1) ;
	    //skip back tracking 
	    bool back_track=false ;
	    if(path.size()>0)
	      if(path[path.size()-1] == -new_path[path.size()])
		back_track=true;
 	    if(!back_track){
	      QDPIO::cout<<" Added path: "<<new_path<<endl;
	      multi1d<LatticeFermion> vec_mu(vec.size()) ;
	      for(int j(0);j<vec.size();j++){
		if(sign>0)
		  vec_mu[j] = shift(vec[j], FORWARD, mu);
		else
		  vec_mu[j] = shift(vec[j], BACKWARD, mu);
	      }
 	      do_disco(db, Clsk , vec_mu, param, Doo, u, t, trb, p, new_path);
	    } // skip backtracking
	  } // mu
      }// path.size loop
    }// do_disco
    

    void ReadOPTEigCGVecs(multi1d<LatticeFermion>& vec,
			    CholeskyFactors& Clsk, 
			    const string& evecs_file)
    {
      QDPIO::cout<<name<<" : Reading vecs from "
		 << evecs_file <<endl ;
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      
      int Nvecs,ldh ;
      // File XML                                      
      XMLReader file_xml;
      // Open file
      QDPFileReader to(file_xml,evecs_file,QDPIO_SERIAL);
      read(file_xml, "/OptEigInfo/ncurEvals", Nvecs);
      read(file_xml, "/OptEigInfo/ldh", ldh);
      // Added Nvecs and ldh to the Cholesky structure for later...
      Clsk.Nvec = Nvecs;
      Clsk.ldh = ldh;
      vec.resize(Nvecs);
      Clsk.evals.resize(ldh);
      Clsk.H.resize(ldh*ldh);
      Clsk.HU.resize(ldh*ldh);

      for(int v(0);v<Nvecs;v++){
	XMLReader record_xml;
	read(to, record_xml, vec[v]);
      }
      
      XMLReader record_xml;
      read(to, record_xml, Clsk.evals);
      read(to, record_xml, Clsk.H);
      read(to, record_xml, Clsk.HU);
      swatch.stop();
      QDPIO::cout<<name<<" : Time to read vecs= "
		 << swatch.getTimeInSeconds() <<" secs "<<endl ;
    }

    // PRchi returns chitilde = (1 - V Hinv Vdag Sdag S)chi given
    // an input chi vector, and of course the vectors and H. 
    void PRchi(multi1d<multi1d< multi1d<LatticeFermion> > >& quarkstilde,
	       const multi1d< Handle< DilutionScheme<LatticeFermion> > >& quarks,
	       Handle<EvenOddPrecLinearOperator<T,P,Q> >& Doo,
	       CholeskyFactors Clsk , multi1d<LatticeFermion>& vec,
	       const Params::Param_t& param, const P& u){
      char U = 'U';
      int info;

      int ldb = vec.size();//This is the offset that will for now be the size of each vector

      Set timerb;
      timerb.make(TimeSliceRBFunc(3));

      // Now we have to create the Sdag * S * quarks object to put into B
      for(int n(0);n<quarks.size();n++){
	for (int it(0) ; it < quarks[n]->getNumTimeSlices() ; ++it){
	  int t = quarks[n]->getT0(it) ;
	  int Nrhs = quarks[n]->getDilSize(it) ;
	  multi2d<Complex> B(Nrhs, ldb);
	  for(int i(0); i<ldb;i++){
	    for(int j = 0 ; j <  Nrhs ; j++){
	      /**
		 Here, we are calculating 
		   Sdag S chi
		 where chi (the solution) is Sinv eta, and eta is the noise vector (the source)
		 Thus, we can save time by using the source, and calculate
		   Sdag eta,
		 and that's what is being done here
	       **/
	      LatticeFermion qsrc     = quarks[n]->dilutedSource(it,j);
	      LatticeFermion SdagSchi = zero ;
	      Doo->oddOddLinOp(SdagSchi,qsrc,MINUS); 
	      // Sum over (odd) lattice sites, so we return a complex number for B
	      // SHOULDN'T THIS BE ONLY SITES THAT INCLUDE THE TIMESLICE WE'RE ON?
	      //	      B[j][i] = sum(localInnerProduct(vec[i],SdagSchi),rb[1]);
	      B[j][i] = sum(localInnerProduct(vec[i],SdagSchi),timerb[2*t+1]);
	    }
	  }

	  // This could be better done.
#if BASE_PRECISION == 32
	  int r = QDPLapack::cpotrs(U, Clsk.Nvec, Nrhs, Clsk.HU, Clsk.ldh, B, ldb, info);
#else
	  int r = QDPLapack::zpotrs(U, Clsk.Nvec, Nrhs, Clsk.HU, Clsk.ldh, B, ldb, info);
#endif
	  QDPIO::cout<<"PRchi cpotrs r = "<<r<<endl;
	  QDPIO::cout<<"PRchi cpotrs info = "<<info<<endl;

	  for(int j = 0 ; j <  Nrhs ; j++){
	    LatticeFermion vB = zero; //B[j][0]*vec[0];
	    LatticeFermion q     = quarks[n]->dilutedSolution(it,j);
	    for(int i(0); i<ldb;i++)
	      vB += B[j][i]*vec[i];
	    quarkstilde[n][it][j] = q - vB;
	    cout<<"Norm of PRchi: "<<norm2(quarkstilde[n][it][j])<<endl;
	  }
	}
      }
    }// End of PRchi call...
    
    /**
       This version is called if we have no EigCG vectors, so PR = 1
     **/
    void PRchi(multi1d<multi1d< multi1d<LatticeFermion> > >& quarkstilde,
	       multi1d< Handle< DilutionScheme<LatticeFermion> > >& quarks){
      
      for(int n(0);n<quarks.size();n++){
	for (int it(0) ; it < quarks[n]->getNumTimeSlices() ; it++){
	  for (int j(0) ; j < quarks[n]->getDilSize(it) ; j++){
	    quarkstilde[n][it][j] = quarks[n]->dilutedSolution(it,j);
	  }
	}
      }
    }// End of PRchi call...
    
  //--------------------------------------------------------------
  // Function call
  //  void 
  //InlineDiscoEigCG::operator()(unsigned long update_no,
  //				XMLWriter& xml_out) 
    void InlineMeas::operator()(unsigned long update_no,
				XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != ""){
	string xml_file = makeXMLFileName(params.xml_file, update_no);
	
	push(xml_out, "discoEigCG");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);
	
	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else{
	func(update_no, xml_out);
      }
    }
    
    
    // Function call
    //void 
    //InlineDiscoEigCG::func(unsigned long update_no,
    //			  XMLWriter& xml_out) 
    void InlineMeas::func(unsigned long update_no,
			  XMLWriter& xml_out) 
    {
      START_CODE();
      
      StopWatch snoop;
      snoop.reset();
      snoop.start();
      
      /*** Stuff to set up everything! ***/
      // Test and grab a reference to the gauge field
      XMLBufferWriter gauge_xml;
      try
	{
	  TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	  TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
	}
      catch( std::bad_cast ) 
	{
	  QDPIO::cerr << InlineDiscoEoEigCGEnv::name << ": caught dynamic cast error" 
		      << endl;
	  QDP_abort(1);
	}
      catch (const string& e) 
	{
	  QDPIO::cerr << name << ": map call failed: " << e 
		      << endl;
	  QDP_abort(1);
	}
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      
      push(xml_out, "discoEigCG");
      write(xml_out, "update_no", update_no);
      
      QDPIO::cout << name << ": Disconnected diagrams with eigCG vectors" << endl;
      
      proginfo(xml_out);    // Print out basic program info
      
      // Write out the input
      params.write(xml_out, "Input");
      
      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);
      
      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);
      
      // First calculate some gauge invariant observables just for info.
      // This is really cheap.
      MesPlq(xml_out, "Observables", u);
      
      // Save current seed
      Seed ran_seed;
      QDP::RNG::savern(ran_seed);
      
      //
      // Read the source and solutions
      //
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      
      int N_quarks = params.param.chi.size() ;

      multi1d< Handle< DilutionScheme<LatticeFermion> > > quarks(N_quarks);

      try{
	// Loop over quark dilutions
	for(int n(0); n < params.param.chi.size(); n++){
	  const GroupXML_t& dil_xml = params.param.chi[n];
	  std::istringstream  xml_d(dil_xml.xml);
	  XMLReader  diltop(xml_d);
	  QDPIO::cout << "Dilution type = " << dil_xml.id << endl;
	  quarks[n] = 
	    TheFermDilutionSchemeFactory::Instance().createObject(dil_xml.id, 
								  diltop, 
								  dil_xml.path);
	}
      }
      catch(const std::string& e){
	QDPIO::cerr << name << ": Caught Exception constructing dilution scheme: " << e << endl;
	QDP_abort(1);
      }
      
      
      //-------------------------------------------------------------------
      //Sanity checks	

      //Another Sanity check, the three quarks must all be 
      //inverted on the same cfg
      for (int n = 1 ; n < N_quarks ; n++){
	if (quarks[0]->getCfgInfo() != quarks[n]->getCfgInfo()){
	  QDPIO::cerr << name << " : Quarks do not contain the same cfg info";
	  QDPIO::cerr << ", quark "<< n << endl;
	  QDP_abort(1);
	}		
      }

      //Also ensure that the cfg on which the inversions were performed 
      //is the same as the cfg that we are using
      {	
	std::string cfgInfo; 
	
	//Really ugly way of doing this (Robert Heeeelp!!)
	XMLBufferWriter top;
	write(top, "Config_info", gauge_xml);
	XMLReader from(top);
	XMLReader from2(from, "/Config_info");
	std::ostringstream os;
	from2.print(os);
	
	cfgInfo = os.str();
	
	if (cfgInfo != quarks[0]->getCfgInfo()){
	  QDPIO::cerr << name << " : Quarks do not contain the same";
	  QDPIO::cerr << " cfg info as the gauge field." ;
	  QDPIO::cerr << "gauge: XX"<<cfgInfo<<"XX quarks: XX" ;
	  QDPIO::cerr << quarks[0]->getCfgInfo()<<"XX"<<  endl;
	  QDP_abort(1);
	}
      }      

      int decay_dir = quarks[0]->getDecayDir();
      //
      // Initialize the slow Fourier transform phases
      //
      SftMom phases(params.param.p2_max, false, decay_dir);

      // Create even/odd subsets on each timeslice:
      Set timerb;
      timerb.make(TimeSliceRBFunc(decay_dir));
      
      // Another sanity check. The seeds of all the quarks must be different
      // and thier decay directions must be the same 
      for(int n = 1 ; n < quarks.size(); n++){
	if(toBool(quarks[n]->getSeed()==quarks[0]->getSeed())){
	  QDPIO::cerr << name << ": error, quark seeds are the same" << endl;
	  QDP_abort(1);
	}
	
	if(toBool(quarks[n]->getDecayDir()!=quarks[0]->getDecayDir())){
	  QDPIO::cerr<<name<< ": error, quark decay dirs do not match" <<endl;
	  QDP_abort(1);
	}
      }

      // Instantiate the Dirac operator:
      Handle<EvenOddPrecLinearOperator<T,P,Q> > Doo=createOddOdd_Op(params.param,u);

      //Now read evecs from disk, if file is specified.
      multi1d<LatticeFermion> vec; // the vectors
      CholeskyFactors Clsk; // the Cholesky Factors 
      if (params.named_obj.evecs_file!="")
	ReadOPTEigCGVecs(vec,Clsk,params.named_obj.evecs_file);
      
      map< KeyOperator_t, ValOperator_t > data ;
      
      multi1d<multi1d< multi1d<LatticeFermion> > > quarkstilde(quarks.size());
      for(int n(0);n<quarks.size();n++){
	quarkstilde[n].resize(quarks[n]->getNumTimeSlices());
	for (int it(0) ; it < quarks[n]->getNumTimeSlices() ; it++){
	  quarkstilde[n][it].resize(quarks[n]->getDilSize(it));
	}
      }
      
      /*** Everything is initialized, now to the real code ***/
      
      // chitilde = P_R S^-1 \eta_o
      // If there is no evecs file, then PR = 1
      if (params.named_obj.evecs_file!="")
	PRchi(quarkstilde, quarks, Doo, Clsk, vec, params.param, u);
      else
	PRchi(quarkstilde, quarks);

      for(int n(0);n<quarks.size();n++){
	for (int it(0) ; it < quarks[n]->getNumTimeSlices() ; it++){
	  int t = quarks[n]->getT0(it) ;
	  QDPIO::cout<<" Doing quark: "<< n <<endl ;
	  QDPIO::cout<<"   quark: "<< n <<" has "<<quarks[n]->getDilSize(it);
	  QDPIO::cout<<" dilutions on time slice "<< t <<endl ;
	  for(int i = 0 ; i <  quarks[n]->getDilSize(it) ; i++){
	    QDPIO::cout<<"   Doing dilution : "<<i<<endl ;
	    multi1d<short int> d ;
	    // Now, I want to do the two trace terms that have the chitilde. Thus
	    // The q chosen should be from quarkstilde = PR*quarks->dilutedSolution()
	    LatticeFermion qbar  = quarks[n]->dilutedSource(it,i);
	    LatticeFermion q     = quarkstilde[n][it][i];
	    QDPIO::cout<<"   Starting recursion "<<endl ;
	    // this does both terms with the noise vectors
	    // the 1 means to do sum over odd sites...
	    do_disco(data, qbar, q, Doo, phases, t, timerb[2*t+1], d, params.param.max_path_length);

	    QDPIO::cout<<" done with recursion! "
		       <<"  The length of the path is: "<<d.size()<<endl ;
	  }
	  QDPIO::cout<<" Done with dilutions for quark: "<<n <<endl ;
	}
      }

      // Okay, first we are going to normalize all the pieces above by the number of 
      // quarks...
      SerialDBKey <KeyOperator_t> key ;
      SerialDBData<ValOperator_t> val ;
      map< KeyOperator_t, ValOperator_t >::iterator it;

      for(it=data.begin();it!=data.end();it++){
	key.key()  = it->first  ;
	val.data().op.resize(it->second.op.size()) ;
	for(int i(0);i<it->second.op.size();i++){
          val.data().op[i] = it->second.op[i]/toDouble(quarks.size());
	}
      }
      
      // Now calculate the other pieces which need no normalization...
      // Note using just timeslices for quarks[0]...may be a problem
      // if quarks dilute on diff timeslices...
      
      for (int it(0) ; it < quarks[0]->getNumTimeSlices() ; it++){
	multi1d<short int> d ;
	int t = quarks[0]->getT0(it) ;
	// Now let's do the Tr[Dee^-1 gamma]
	do_disco(data,params.param, u, phases, t, timerb[2*t+0], d);
      }
      
      if (params.named_obj.evecs_file!=""){// Only if we have vectors...
	for (int it(0) ; it < quarks[0]->getNumTimeSlices() ; it++){
	  multi1d<short int> d ;
	  int t = quarks[0]->getT0(it) ;
	  QDPIO::cout<<"   Starting recursion again"<<endl ;
	  // Part with deflation..
	  do_disco(data, Clsk, vec, params.param, Doo, u, t, timerb[2*t+1], phases,d);
	  
	  QDPIO::cout<<" done with recursion! "
		     <<"  The length of the path is: "<<d.size()<<endl ;
	}
      }
            
      //After all the pieces are computed we write the final result to the database
      // DB storage          
      BinaryStoreDB<SerialDBKey<KeyOperator_t>,SerialDBData<ValOperator_t> > qdp_db;
      
      // Open the file, and write the meta-data and the binary for this operator
      {
	XMLBufferWriter file_xml;
	
	push(file_xml, "DBMetaData");
	write(file_xml, "id", string("discoEigElemOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", decay_dir);
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	std::string file_str(file_xml.str());
	qdp_db.setMaxUserInfoLen(file_str.size());

	qdp_db.open(params.named_obj.op_db_file, O_RDWR | O_CREAT, 0664);

	qdp_db.insertUserdata(file_str);
      }

      // Store all the data
      for(it=data.begin();it!=data.end();it++){
	key.key()  = it->first  ;
	val.data().op.resize(it->second.op.size()) ;
	// DON'T normalize to number of quarks here, because we only do it on a 
	// certain number of the terms in the trace...done above
	for(int i(0);i<it->second.op.size();i++)
          val.data().op[i] = it->second.op[i];
	qdp_db.insert(key,val) ;
      }

      // Close the namelist output file XMLDAT
      pop(xml_out);     // DiscoEigCG
      
      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;
      
      QDPIO::cout << name << ": ran successfully" << endl;
      
      END_CODE();
    } 
  }  // namespace InlineDiscoEoEigCGEnv
}// namespace chroma
