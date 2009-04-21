

#include "fermact.h"
#include "meas/inline/pbp/inline_psibarpsi_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "util/info/proginfo.h"
#include "meas/inline/io/named_objmap.h"

#include "meas/pbp/pbp.h"


namespace Chroma
{

  void read(XMLReader& xml, const string& path, InlinePsiBarPsiEnv::Params::Param_t& param)
  {
    XMLReader paramtop(xml, path);
	
	int version;
	read(paramtop, "version", version);
	
	switch (version)
	{
	case 1:
	  param.fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");
	  param.invParam = readXMLGroup(paramtop, "InvertParam", "invType");
	
	  read(paramtop, "Mass", param.mass);
	  if (paramtop.count("ichiral") == 1) {
		read(paramtop, "ichiral", param.ichiral);
	  } else {
		param.ichiral = 0;
	  }
	
	  break;
	
	default:
	  QDPIO::cerr << "InlinePsiBarPsiParams::Param_t: " << version
		<< " unsupported." << endl;
	  QDP_abort(1);
	}
  }
  
  void write(XMLWriter& xml, const string& path, InlinePsiBarPsiEnv::Params::Param_t& param)
  {
    push(xml, path);
	
	int version = 1;
	write(xml, "version", version);
	write(xml, "mass", param.mass);
	write(xml, "ichiral", param.ichiral);

	xml << param.fermact.xml;
	xml << param.invParam.xml;
	
	pop(xml);
  }
  
  void read(XMLReader& xml, const string& path, InlinePsiBarPsiEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);
	
	read(inputtop, "gauge_id", input.gauge_id);
  }
  
  void write(XMLWriter& xml, const string& path, const InlinePsiBarPsiEnv::Params::NamedObject_t& input)
  {
    push(xml, path);
	
	write(xml, "gauge_id", input.gauge_id);
	
	pop(xml);
  }


  namespace InlinePsiBarPsiEnv
  {
    namespace
	{
	  AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
						const std::string& path)
		{
	  return new InlineMeas(Params(xml_in, path));
	    }
		
		//! Local registration flag
		bool registered = false;
	  }
	  
	  const std::string name = "PSI_BAR_PSI";
	  
	  //! Register all the factories
	  bool registerAll()
	  {
	    bool success = true;
		if (! registered)
		{
	  success &= WilsonTypeFermActsEnv::registerAll();
	  success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	  registered = true;
	    }
		return success;
	  }
  
	Params::Params()
  {
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
		
	  read(paramtop, "Param", param);
	  
	  read(paramtop, "NamedObject", named_obj);
	}
	catch(const std::string& e)
	{
	  QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
	  QDP_abort(1);
	}
  }
  
  void InlineMeas::operator()(unsigned long update_no,
				XMLWriter& xml_out)
  {

	push(xml_out, "PsiBarPsi");
	write(xml_out, "update_no", update_no);
	pop(xml_out);
	
	func(update_no, xml_out);
  }

  void InlineMeas::func(unsigned long update_no,
				XMLWriter& xml_out)
  {
    START_CODE();
	
	Double pbp;
	int maxCG = 0;
	int n_congrd = 0; 
	XMLBufferWriter gauge_xml;
	
	multi1d<LatticeColorMatrix> u;
	try {
		u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
		TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
	} 
	catch ( std::bad_cast ) {
		QDPIO::cerr << name << ": caught dynamic cast error" << endl;
		QDP_abort(1);
	}
	catch (const string& e) {
		QDPIO::cerr << name << ": map call failed: " << e << endl;
		QDP_abort(1);
	}
	 
	proginfo(xml_out);
	
	write(xml_out, "Config_info", gauge_xml);
	
	// Initialize fermion action
	
	std::istringstream xml_s(params.param.fermact.xml);
	XMLReader	fermacttop(xml_s);
	QDPIO::cout << "FermAct = " << params.param.fermact.id << endl;
	
	
	try
	{
  StopWatch swatch;
  swatch.reset();
  QDPIO::cout << "Try the various factories" << endl;
  
  typedef LatticeFermion				T;
  typedef multi1d<LatticeColorMatrix>	P;
  typedef multi1d<LatticeColorMatrix>	Q;
  
  Handle< FermionAction<T,P,Q> >
    S_f(TheFermionActionFactory::Instance().createObject(params.param.fermact.id,
									fermacttop,
									params.param.fermact.path));
									
  Handle< FermState<T,P,Q> > state(S_f->createState(u));
  
  QDPIO::cout << "Suitable factory found: do the measurements" << endl;
  Handle< SystemSolver<T> > qprop(S_f->qprop(state, params.param.invParam));
  
  swatch.start();
  
  MesPbp(qprop, state, params.param.mass, params.param.ichiral, xml_out, "PsiBarPsi", params.param.fermact.id);
  
  swatch.stop();
  QDPIO::cout << "PsiBarPsi computed: time= "
			  << swatch.getTimeInSeconds()
			  << " secs" << endl;
	}
	catch ( std::bad_cast )
	{
  QDPIO::cerr << name << ": caught dynamic cast error"
		<< endl;
  QDP_abort(1);
    }
	catch (const std::string& e)
	{
  QDPIO::cout << name << ": caught exception with fermion action: " << e << endl;
    }
	
	pop(xml_out);
	
	QDPIO::cout << name << ": ran successfully" << endl;
	
	END_CODE();
  } // end func()
  
  } // end namespace InlinePsiBarPsiEnv

} // end namespace Chroma
 


  