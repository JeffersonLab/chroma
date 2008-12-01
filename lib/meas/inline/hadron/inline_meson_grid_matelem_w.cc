// $Id: inline_meson_grid_matelem_w.cc,v 3.1 2008-12-01 03:09:17 kostas Exp $
/*! \file
 * \brief Inline measurement of meson operators via colorvector matrix elements
 */

#include "handle.h"
#include "meas/inline/hadron/inline_grid_prop_w.h"
#include "meas/inline/hadron/inline_meson_grid_matelem_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/displace.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_val_db.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/sources/diluteGrid_source_const.h"

#include "meas/inline/io/named_objmap.h"

#define MATELEM_TYPE_YALE       12
#define MATELEM_TYPE_DIAGONAL   11
#define MATELEM_TYPE_GENERIC    10

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineMesonGridMatElemEnv 
  { 

    //! Reader for Grid descriptor
    void read(XMLReader& xml, const string& path, InlineMesonGridMatElemEnv::Params::Param_t::Grid_t& g)
    {
      XMLReader inputtop(xml, path);
      read(inputtop, "spatial_mask_size",g.spatial_mask_size);
      read(inputtop, "spatial_masks", g.spatial_masks);
      read(inputtop, "decay_dir", g.decay_dir);
    }

    //! Writer for Grid descriptor
    void write(XMLWriter& xml,const string& path,const InlineMesonGridMatElemEnv::Params::Param_t::Grid_t& g)
    {
      push(xml, path);
      write(xml, "spatial_mask_size",g.spatial_mask_size);
      write(xml, "spatial_masks", g.spatial_masks);
      write(xml, "decay_dir", g.decay_dir);
      pop(xml);
    }

    // Reader for input parameters
    void read(XMLReader& xml, const string& path, InlineMesonGridMatElemEnv::Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);
    
      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	/**************************************************************************/
	break;

      default :
	/**************************************************************************/

	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "displacement_length", param.displacement_length);
      read(paramtop, "displacement_list", param.displacement_list);
      read(paramtop, "Grid", param.grid);

      param.link_smearing  = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
      
    }


    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const InlineMesonGridMatElemEnv::Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "displacement_length", param.displacement_length);
      write(xml, "displacement_list", param.displacement_list);
      write(xml, "Grid", param.grid);

      xml << param.link_smearing.xml;

      pop(xml);
    }

    //! Read named objects 
    void read(XMLReader& xml, const string& path, InlineMesonGridMatElemEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "meson_op_file", input.meson_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const InlineMesonGridMatElemEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "meson_op_file", input.meson_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const InlineMesonGridMatElemEnv::Params& param)
    {
      param.writeXML(xml, path);
    }
  }


  namespace InlineMesonGridMatElemEnv 
  { 
    // Anonymous namespace for registration
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

    const std::string name = "MESON_GRID_MATELEM";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= LinkSmearingEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    //! Anonymous namespace
    /*! Diagnostic stuff */
    namespace
    {
      StandardOutputStream& operator<<(StandardOutputStream& os, const multi1d<int>& d)
      {
	if (d.size() > 0)
	{
	  os << d[0];

	  for(int i=1; i < d.size(); ++i)
	    os << " " << d[i];
	}

	return os;
      }
    }


    //----------------------------------------------------------------------------
    // Param stuff
    Params::Params()
    { 
      frequency = 0; 
      param.mom2_max = 0;
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


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) const
    {
      push(xml_out, path);
    
      // Parameters for source construction
      write(xml_out, "Param", param);

      // Write out the output propagator/source configuration info
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }


    //----------------------------------------------------------------------------
    //! Meson operator
    struct KeyMesonOperator_t
    {
      int                t_slice;      /*!< Meson operator time slice */
      multi1d<int>       displacement; /*!< Displacement dirs of right colorvector */
      multi1d<int>       mom;          /*!< D-1 momentum of this operator */
    };

    //! Meson operator
    struct ValMesonOperator_t
    {
      int                type_of_data ;      /*!< Flag indicating type of data (maybe trivial) */
      //multi1d< multi2d< ComplexD> > y ; /*! yale sparse matrix format */
      //multi1d<int> row ; /*! grid row index for yale format*/
      //multi1d<int> col ; /*! grid column index for yale format */
      multi1d<ComplexD> diag ; /* diagnonal matrix */
      //generic dense matrix format: Ngrid*Nc x Ngrid*Nc not implemented
      // indexing i = c + Nc*g with c is color and g is the grid index
      multi2d<ComplexD>  gen;  /*!< grid source and sink with momentum projection */
    };


    //----------------------------------------------------------------------------
    //! Holds key and value as temporaries
    struct KeyValMesonOperator_t
    {
      SerialDBKey<KeyMesonOperator_t>  key;
      SerialDBData<ValMesonOperator_t> val;
    };


    //----------------------------------------------------------------------------
    //! KeyMesonElementalOperator reader
    void read(BinaryReader& bin, KeyMesonOperator_t& param)
    {
      read(bin, param.t_slice);
      read(bin, param.displacement);
      read(bin, param.mom);
    }

    //! MesonElementalOperator write
    void write(BinaryWriter& bin, const KeyMesonOperator_t& param)
    {
      write(bin, param.t_slice);
      write(bin, param.displacement);
      write(bin, param.mom);
    }

    //! MesonElementalOperator reader
    void read(XMLReader& xml, const std::string& path, KeyMesonOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "t_slice", param.t_slice);
      read(paramtop, "displacement", param.displacement);
      read(paramtop, "mom", param.mom);
    }

    //! MesonElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const KeyMesonOperator_t& p)
    {
      push(xml, path);

      write(xml, "t_slice", p.t_slice);
      write(xml, "displacement", p.displacement);
      write(xml, "mom", p.mom);

      pop(xml);
    }


    //----------------------------------------------------------------------------
    //! MesonElementalOperator reader
    void read(BinaryReader& bin, ValMesonOperator_t& param)
    {
      read(bin, param.type_of_data);

      if(param.type_of_data == MATELEM_TYPE_GENERIC){
	int n; // N=Ngrid*Nc
	read(bin, n);    // the size is always written, even if 0
	param.gen.resize(n,n);
	for(int i=0; i < n; ++i)
	  for(int j=0; j < n; ++j)
	    read(bin, param.gen(i,j));
      }
      if(param.type_of_data == MATELEM_TYPE_DIAGONAL){
	int n; // N=Ngrid*Nc
	read(bin, n);    // the size is always written, even if 0
	param.diag.resize(n); // always resize to Nc number of colors...
	for(int i=0; i < n; ++i)
	  read(bin, param.diag[i]);
      }
      if(param.type_of_data == MATELEM_TYPE_YALE){
      }
      
    }

    //! MesonElementalOperator write
    void write(BinaryWriter& bin, const ValMesonOperator_t& param)
    {
      write(bin, param.type_of_data);

      if(param.type_of_data == MATELEM_TYPE_GENERIC){
	int n = param.gen.size1();  // all sizes the same
	write(bin, n);
	for(int i=0; i < n; ++i)
	  for(int j=0; j < n; ++j)
	    write(bin, param.gen(i,j));
      }
      if(param.type_of_data == MATELEM_TYPE_DIAGONAL){
	int n = param.diag.size();  // all sizes the same
	write(bin, n);
	for(int i=0; i < n; ++i)
	  write(bin, param.diag[i]);
      }
    }


    //----------------------------------------------------------------------------
    //! Normalize just one displacement array
    multi1d<int> normDisp(const multi1d<int>& orig)
    {
      START_CODE();

      multi1d<int> disp;
      multi1d<int> empty; 
      multi1d<int> no_disp(1); no_disp[0] = 0;

      // NOTE: a no-displacement is recorded as a zero-length array
      // Convert a length one array with no displacement into a no-displacement array
      if (orig.size() == 1)
      {
	if (orig == no_disp)
	  disp = empty;
	else
	  disp = orig;
      }
      else
      {
	disp = orig;
      }

      END_CODE();

      return disp;
    } // void normDisp


    //------------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "MesonGridMatElem");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else
      {
	func(update_no, xml_out);
      }
    }


    // Function call
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      
      StopWatch swiss;
			
      // Test and grab a reference to the gauge field
      XMLBufferWriter gauge_xml;
      try
      {
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": map call failed: " << e << endl;
	QDP_abort(1);
      }
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "MesonGridMatElem");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Meson grid matrix element" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      //First calculate some gauge invariant observables just for info.
      //This is really cheap.
      MesPlq(xml_out, "Observables", u);

      //
      // Initialize the slow Fourier transform phases
      //
      SftMom phases(params.param.mom2_max, false, params.param.grid.decay_dir);

      //
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;

      try
      {
	std::istringstream  xml_l(params.param.link_smearing.xml);
	XMLReader  linktop(xml_l);
	QDPIO::cout << "Link smearing type = " << params.param.link_smearing.id << endl;
	
	
	Handle< LinkSmearing >
	  linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.param.link_smearing.id,
								       linktop, 
								       params.param.link_smearing.path));

	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception link smearing: " << e << endl;
	QDP_abort(1);
      }

      // Record the smeared observables
      MesPlq(xml_out, "Smeared_Observables", u_smr);

      // Keep track of no displacements and zero momentum
      multi1d<int> no_displacement;
      multi1d<int> zero_mom(3); zero_mom = 0;

      //prepare the grids
      StopWatch roloi;
      roloi.reset();
      roloi.start();
      DiluteGridQuarkSourceConstEnv::Params srcParams ;
      srcParams.smear = false ;
      srcParams.j_decay = params.param.grid.decay_dir ;
      srcParams.spatial_mask_size = params.param.grid.spatial_mask_size ;
      multi2d<LatticeFermion> vec(Nc,params.param.grid.spatial_masks.size()) ;
      srcParams.spin = -1 ;// do all spins at once                          
      srcParams.t_source = -1 ; // do all time slices at once             
      for(int c(0);c<Nc;c++){// loop over colors          
	srcParams.color = c;
	for(int g(0);g<params.param.grid.spatial_masks.size();g++){
	  srcParams.spatial_mask =params.param.grid.spatial_masks[g];
	  QDPIO::cout << name << ": Doing color " << c
		      << " and  grid " << g<< " ( of "
		      << params.param.grid.spatial_masks.size() << ")"<<endl;
	  DiluteGridQuarkSourceConstEnv::SourceConst<LatticeFermion> GridSrc(srcParams);
	  vec(c,g) = GridSrc(u);
	}//g                                                               
      }//c                                                          
      roloi.stop();
      QDPIO::cout << "Sources constructed: time= "
		  << roloi.getTimeInSeconds()
		  << " secs" << endl;

      //
      // Meson operators
      //
      // The creation and annihilation operators are the same without the
      // spin matrices.
      //
      QDPIO::cout << "Building meson operators" << endl;

      // DB storage
      BinaryFxStoreDB<SerialDBKey<KeyMesonOperator_t>, 
	              SerialDBData<ValMesonOperator_t> > 
	qdp_db(params.named_obj.meson_op_file, 10*1024*1024, 64*1024);

      push(xml_out, "MesonGridOps");

      // Loop over all time slices for the source. This is the same 
      // as the subsets for  phases

      int Ngrids = params.param.grid.spatial_masks.size() ;
      // Loop over each operator 
      int format_type = MATELEM_TYPE_DIAGONAL ;
      for(int l=0; l < params.param.displacement_list.size(); ++l){
	multi1d<int> disp = normDisp(params.param.displacement_list[l]);
	if(disp.size()!=0)
	  format_type = MATELEM_TYPE_GENERIC ;
      }

      for(int l=0; l < params.param.displacement_list.size(); ++l)
      {
	QDPIO::cout << "Elemental operator: op = " << l << endl;

	// Make sure displacement is something sensible
	multi1d<int> disp = normDisp(params.param.displacement_list[l]);

	QDPIO::cout << "displacement = " << disp << endl;

	// Build the operator
	swiss.reset();
	swiss.start();

	// Big loop over the momentum projection
	for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
	{
	  // The keys for the spin and displacements for this particular elemental operator
	  // No displacement for left colorvector, only displace right colorvector
	  // Invert the time - make it an independent key
	  multi1d<KeyValMesonOperator_t> buf(phases.numSubsets());
	  for(int t=0; t < phases.numSubsets(); ++t)
	  {
	    buf[t].key.key().t_slice       = t;
	    buf[t].key.key().mom           = phases.numToMom(mom_num);
	    buf[t].key.key().displacement  = disp; // only right vector
	    buf[t].val.data().type_of_data = format_type ;
	    if ( format_type == MATELEM_TYPE_DIAGONAL ){ // no displacement 
	      buf[t].val.data().diag.resize(Ngrids*Nc);
	    }
	    else{
	      buf[t].val.data().gen.resize(Ngrids*Nc,Ngrids*Nc);
	    }
	  }
	  int j(0);
	  for(int g(0) ; g<Ngrids; ++g){
	    for(int c(0) ; c<Nc ; c++,j++){
	      QDPIO::cout<<"source index: "<<j<<" ("<<g<<","<<c<<")"<<endl ;
	      // Displace the right vector and multiply by the momentum phase

	      LatticeColorVector q = peekSpin(vec(c,g),0);

	      LatticeColorVector shift_vec = 
		phases[mom_num]*displace(u_smr, q, 
					 params.param.displacement_length, disp);
	      if(format_type == MATELEM_TYPE_DIAGONAL){//diagonal data format
		LatticeComplex lop = localInnerProduct(q, shift_vec);
		
		// Slow fourier-transform                                          
		multi1d<ComplexD> op_sum = sumMulti(lop, phases.getSet());
		for(int t=0; t < op_sum.size(); ++t)
		  buf[t].val.data().diag[j] = op_sum[t];
	      }
	      else{
		int i(0);
		for(int gg(0) ; gg<Ngrids; ++gg){
		  for(int cc(0) ; cc<Nc ; cc++,i++){
		    QDPIO::cout<<"sink index: "<<i <<"("<<gg<<","<<cc<<")"<<endl ;
		    // Contract over color indices
		    // Do the relevant quark contraction
		    LatticeColorVector qbar = peekSpin(vec(cc,gg),0);
		    LatticeComplex lop = localInnerProduct(qbar, shift_vec);

		    // Slow fourier-transform
		    multi1d<ComplexD> op_sum = sumMulti(lop, phases.getSet());
		    
		    for(int t=0; t < op_sum.size(); ++t)
		      buf[t].val.data().gen(i,j) = op_sum[t];
		    //	      write(xml_out, "elem", key.key());  // debugging
		  } // end for cc
		}// end for gg
	      }//else
	    } // end for c
	  } // end for g

	  QDPIO::cout << "insert: mom= " << phases.numToMom(mom_num) << " displacement= " << disp << endl; 
	  for(int t=0; t < phases.numSubsets(); ++t)
	    qdp_db.insert(buf[t].key, buf[t].val);

	} // mom_num

	swiss.stop();

	QDPIO::cout << "Meson operator= " << l 
		    << "  time= "
		    << swiss.getTimeInSeconds() 
		    << " secs" << endl;

      } // for l

      pop(xml_out); // MesonGridOps

      // Write the meta-data and the binary for this operator
      swiss.reset();
      swiss.start();
      {
	XMLBufferWriter file_xml;

	push(file_xml, "MesonOperators");
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", params.param.grid.decay_dir);
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	write(file_xml, "Op_Info",params.param.displacement_list);
	pop(file_xml);

	qdp_db.insertUserdata(file_xml.str());
      }
      swiss.stop();

      QDPIO::cout << "Meson Operator written:"
		  << "  time= " << swiss.getTimeInSeconds() << " secs" << endl;

      // Close the namelist output file XMLDAT
      pop(xml_out);     // MesonGridMatElemtor

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } // func
  } // namespace InlineMesonGridMatElemEnv

  /*! @} */  // end of group hadron

} // namespace Chroma
