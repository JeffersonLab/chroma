// $Id: inline_rotate_spin_w.cc,v 3.1 2008-05-16 21:50:49 edwards Exp $
/*! \file
 * \brief Inline task to spin rotate to a Dirac basis
 *
 * Spin rotate a named object
 */

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/hadron/inline_rotate_spin_w.h"
#include "util/ferm/diractodr.h"

namespace Chroma 
{ 
  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineRotateSpinEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "input_id", input.input_id);
    write(xml, "output_id", input.output_id);
    write(xml, "object_type", input.object_type);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineRotateSpinEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "input_id", input.input_id);
    read(inputtop, "output_id", input.output_id);
    read(inputtop, "object_type", input.object_type);
  }


  //! IO function map environment
  /*! \ingroup inlinehadron */
  namespace RotateSpinObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      struct DumbDisambiguator {};

      //! RotateSpin function map
      /*! \ingroup inlinehadron */
      typedef SingletonHolder< 
	FunctionMap<DumbDisambiguator,
		    void,
		    std::string,
		    TYPELIST_2(const string&, const string&),
		    void (*)(const string& output_id, const string& input_id),
		    StringFunctionMapError> >
      TheRotateSpinObjFuncMap;


      //! Transform a wilson-like fermion object. This works only for non-array objects.
      // This only works for a fermion-like object. A propagator needs 2 spin rotations
      void rotateDRtoDiracFerm(const string& output_id, const string& input_id)
      {
	// Save some typing
	typedef LatticeFermion   T;

	// Grab the source
	const T& input_obj = 
	  TheNamedObjMap::Instance().getData<T>(input_id);

	// The spin basis matrix to go from DR to Dirac
	SpinMatrix rotate_mat(adj(DiracToDRMat()));

	XMLReader input_file_xml, input_record_xml;
	TheNamedObjMap::Instance().get(input_id).getFileXML(input_file_xml);
	TheNamedObjMap::Instance().get(input_id).getRecordXML(input_record_xml);

	// Create space for the target
	TheNamedObjMap::Instance().create<T>(output_id);
	TheNamedObjMap::Instance().get(output_id).setFileXML(input_file_xml);
	TheNamedObjMap::Instance().get(output_id).setRecordXML(input_record_xml);

	// Do the actual rotation
	TheNamedObjMap::Instance().getData<T>(output_id) = rotate_mat * input_obj;
      }

      //! Transform a wilson-like fermion object. This works only for non-array objects.
      // This only works for a fermion-like object. A propagator needs 2 spin rotations
      void rotateDiractoDRFerm(const string& output_id, const string& input_id)
      {
	// Save some typing
	typedef LatticeFermion   T;

	// Grab the source
	const T& input_obj = 
	  TheNamedObjMap::Instance().getData<T>(input_id);

	// The spin basis matrix to go from Dirac to DR
	SpinMatrix rotate_mat(DiracToDRMat());

	XMLReader input_file_xml, input_record_xml;
	TheNamedObjMap::Instance().get(input_id).getFileXML(input_file_xml);
	TheNamedObjMap::Instance().get(input_id).getRecordXML(input_record_xml);

	// Create space for the target
	TheNamedObjMap::Instance().create<T>(output_id);
	TheNamedObjMap::Instance().get(output_id).setFileXML(input_file_xml);
	TheNamedObjMap::Instance().get(output_id).setRecordXML(input_record_xml);

	// Do the actual rotation
	TheNamedObjMap::Instance().getData<T>(output_id) = rotate_mat * input_obj;
      }


      //! Local registration flag
      bool registered = false;

    }  // end namespace RotateSpinCallMap

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
//	success &= TheRotateSpinObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
//									 rotateSpinObj<LatticePropagator>);
	success &= TheRotateSpinObjFuncMap::Instance().registerFunction(string("LatticeFermion:DR-to-Dirac"), 
									rotateDRtoDiracFerm);
	success &= TheRotateSpinObjFuncMap::Instance().registerFunction(string("LatticeFermion:Dirac-to-DR"), 
									rotateDRtoDiracFerm);
	registered = true;
      }
      return success;
    }
  }  // end CallMap namespace


  namespace InlineRotateSpinEnv 
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

    const std::string name = "ROTATE_SPIN_NAMED_OBJECT";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Rotation env
	success &= RotateSpinObjCallMapEnv::registerAll();

	// Inline measurement
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Ids
	read(paramtop, "NamedObject", named_obj);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Ids
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "rotate_spin_object");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": spin rotate an object of type "
		  << params.named_obj.object_type << endl;

      // Grab the input object
      try
      {
	// Spin rotate the object
	RotateSpinObjCallMapEnv::TheRotateSpinObjFuncMap::Instance().callFunction(params.named_obj.object_type,
										  params.named_obj.output_id,
										  params.named_obj.input_id);
      }
      catch (std::bad_cast) 
      {
	QDPIO::cerr << name << ": cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error message: " << e 
		    << endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << name << ": ran successfully" << endl;

      pop(xml_out);

      END_CODE();
    } 

  } // namespace InlineRotateSpinEnv 

} // namespace Chroma
