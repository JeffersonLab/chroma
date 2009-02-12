// $Id: inline_create_colorvecs.cc,v 3.6 2009-02-12 06:34:10 kostas Exp $
/*! \file
 * \brief make color vectors
 *
 * Propagator calculation on a colorvector
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_create_colorvecs.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/block_subset.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/key_grid_prop.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/smear/quark_smearing_aggregate.h"
#include "meas/smear/quark_smearing_factory.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/no_quark_displacement.h"
#include "meas/smear/no_link_smearing.h"

#include "meas/sources/diluteGrid_source_const.h"

#include "meas/inline/io/named_objmap.h"
#include "meas/smear/gaus_smear.h"

namespace Chroma 
{ 
  namespace InlineCreateColorVecsEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineCreateColorVecsEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineCreateColorVecsEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineCreateColorVecsEnv::Params::Param_t::Sources_t& input)
    {
      XMLReader inputtop(xml, path);
      read(inputtop, "spatial_mask_size",input.spatial_mask_size);
      read(inputtop, "spatial_masks", input.spatial_masks);

      read(inputtop, "decay_dir", input.decay_dir);
      
      input.smear = false ;
      if(inputtop.count("Smearing") !=0 ) {
	input.smr = readXMLGroup(inputtop, "Smearing", "wvf_kind");
	input.link_smear = readXMLGroup(inputtop, "LinkSmearing", "LinkSmearingType");
	input.smear = true ;
      }
      
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineCreateColorVecsEnv::Params::Param_t::Sources_t& out)
    {
      push(xml, path);

      write(xml, "spatial_mask_size", out.spatial_mask_size);
      write(xml, "spatial_masks", out.spatial_masks);

      write(xml, "decay_dir", out.decay_dir);

      if(out.smear){
	 xml << out.smr.xml;
	 xml << out.link_smear.xml ;
      }

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineCreateColorVecsEnv::Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Nhits", input.Nhits) ;
      read(inputtop, "OrthoNormal", input.OrthoNormal) ;
      read(inputtop, "Sources"   , input.src)  ;

      input.BlockOrthoNormal=false ;
      if(inputtop.count("BlockOrthoNormal")!=0){
	input.BlockOrthoNormal=true ;
	read(inputtop, "BlockOrthoNormal", input.block) ;
      }
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineCreateColorVecsEnv::Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "Nhits", input.Nhits) ;
      write(xml, "OrthoNormal", input.OrthoNormal) ;
      write(xml, "Sources"   , input.src)  ;

      if(input.BlockOrthoNormal){
	write(xml, "BlockOrthoNormal", input.block) ;
      }
      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineCreateColorVecsEnv::Params& input)
    {
      InlineCreateColorVecsEnv::Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineCreateColorVecsEnv::Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlineCreateColorVecsEnv 


  namespace InlineCreateColorVecsEnv 
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
      
    const std::string name = "CREATE_COLORVECS";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    //-------------------------------------------------------------------------
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

	// Parameters for source construction
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



    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "CreateColorVecs");
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


    // Real work done here
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Test and grab a reference to the gauge field
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;
      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
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

      push(xml_out, "CreateColorVecs");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Create color vectors" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      write(xml_out, "Input", params);

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      //
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;
      
      try
      {
	std::istringstream  xml_l(params.param.src.link_smear.xml);
	XMLReader  linktop(xml_l);
	QDPIO::cout << "Link smearing type = " 
		    << params.param.src.link_smear.id
		    << endl;
		
	Handle< LinkSmearing >
	  linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.param.src.link_smear.id, linktop,params.param.src.link_smear.path));
	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e){
	QDPIO::cerr << name << ": Caught Exception link smearing: "<<e<< endl;
	QDP_abort(1);
      }

      // Record the smeared observables
      MesPlq(xml_out, "Smeared_Observables", u_smr);


      //
      // Create the output files
      //
      try
	{
	  TheNamedObjMap::Instance().create< SubsetVectors<LatticeColorVector> >(params.named_obj.colorvec_id);
	}
      catch (std::bad_cast)
	{
	  QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	  QDP_abort(1);
	}
      catch (const string& e) 
	{
	  QDPIO::cerr << name << ": error creating prop: " << e << endl;
	  QDP_abort(1);
	}

      // Cast should be valid now
      SubsetVectors<LatticeColorVector>& color_vecs=TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(params.named_obj.colorvec_id);

      //QDPIO::cout << name << ": Created named buffer" << endl;
      
      DiluteGridQuarkSourceConstEnv::Params srcParams ;

      NoQuarkDisplacementEnv::Params noQuarkDisp ;
      NoLinkSmearingEnv::Params      noLinkSmear ;
      //create a no displacement xml
      {
	XMLBufferWriter tt ;
	push(tt,"foo");
	noQuarkDisp.writeXML(tt,"Displacement");
	noLinkSmear.writeXML(tt,"LinkSmearing");
	pop(tt);
	XMLReader from(tt);
	srcParams.displace = readXMLGroup(from, "/foo/Displacement", "DisplacementType");
	srcParams.link_smear = readXMLGroup(from, "/foo/LinkSmearing", "LinkSmearingType");
      }
      srcParams.j_decay = params.param.src.decay_dir ;
      srcParams.spatial_mask_size = params.param.src.spatial_mask_size ;
      if(params.param.src.smear){
	srcParams.smear = false ;
	srcParams.smr = params.param.src.smr ;
      }

      // The code goes here
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      /** 

       *Loop for Nhits
       *    apply smearing
       *    if orthonormalize is true to a GramSchmit on them
       *    Compute the v'Smearing v matrix elements for all color vectors
       *    call them evals and monitor convergence on one time slice
       *EndLoop

       **/

      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, params.param.src.decay_dir);
      
      int Nvecs(3*params.param.src.spatial_masks.size()) ;
      color_vecs.getEvectors().resize(Nvecs);
      color_vecs.getEvalues().resize(Nvecs);
      for(int i(0);i<Nvecs;i++){
	color_vecs.getEvalues()[i].weights.resize(phases.numSubsets());
      }

      color_vecs.getDecayDir() = params.param.src.decay_dir ;

      srcParams.spin = 0;
      int count(0);
      for(int c(0);c<Nc;c++){// loop over colors
	srcParams.color = c;
	for(int g(0);g<params.param.src.spatial_masks.size();g++){
	  srcParams.spatial_mask = params.param.src.spatial_masks[g];
	  QDPIO::cout << name << ": Doing  color " << c 
		      << " and  grid " << g<< " ( of " 
		      << params.param.src.spatial_masks.size() << ")"<<endl;
	  //Constuct the source
	  DiluteGridQuarkSourceConstEnv::SourceConst<LatticeFermion>  
	    GridSrc(srcParams);
	  LatticeFermion chi = GridSrc(u_smr);
	  color_vecs.getEvectors()[count] = peekSpin(chi,0) ;
	  if(params.param.OrthoNormal){
	    for(int k(0);k<count;k++){
	      multi1d<DComplex> cc = sumMulti(localInnerProduct(color_vecs.getEvectors()[k],color_vecs.getEvectors()[count]),phases.getSet());
	      for(int t(0);t<phases.numSubsets();t++)
		color_vecs.getEvectors()[count][phases.getSet()[t]] -= 
		  cc[t]*color_vecs.getEvectors()[k] ;
	    }
	    multi1d<Double> norm2 = 
	      sumMulti(localNorm2(color_vecs.getEvectors()[count]),phases.getSet());
	    for(int t(0);t<phases.numSubsets();t++)
	      color_vecs.getEvectors()[count][phases.getSet()[t]] /= sqrt(norm2[t]) ;
	  }
	      
	  count ++;
	}
      }

      Handle< QuarkSmearing<LatticeColorVector> >  Smearing ;
      try{
	std::istringstream  xml_l(params.param.src.smr.xml);
	XMLReader  smrtop(xml_l);
	QDPIO::cout << "ColorVector Smearing type = " <<params.param.src.smr.id ;
	QDPIO::cout << endl;

	Smearing =
	  TheColorVecSmearingFactory::Instance().createObject(params.param.src.smr.id,
							      smrtop,
							      params.param.src.smr.path);
      }
      catch(const std::string& e){
	QDPIO::cerr <<name
		    << ": Caught Exception creating ColorVector smearing object: " 
		    << e << endl;
	QDP_abort(1);
      }
      catch(...){
	QDPIO::cerr <<name<< ": Caught generic exception creating smearing object" 
		    <<   endl;
	QDP_abort(1);
      }
      
      // TheColorVecSmearingFactory::Instance().
      //GausQuarkSmearingEnv::Params QsmrParams(
      //GausQuarkSmearingEnv::QuarkSmear<LatticeColorVector> smearVecs(QsmrParams);
      //for the moment only gaussing smearing works...
      //but this is good enough for what we want to do
      for(int hit(1);hit<params.param.Nhits;hit++){
	for(int i(0);i<Nvecs;i++){
	  QDPIO::cout << name << ": Doing colorvec: "<<i
		      << " hit no: "<<hit<<endl ;
	  (*Smearing)(color_vecs.getEvectors()[i], u_smr);
	  if(params.param.OrthoNormal){
            for(int k(0);k<i;k++){
              multi1d<DComplex> cc = 
		sumMulti(localInnerProduct(color_vecs.getEvectors()[k],
					   color_vecs.getEvectors()[i]),
			 phases.getSet());
              for(int t(0);t<phases.numSubsets();t++)
                color_vecs.getEvectors()[i][phases.getSet()[t]] -=
                  cc[t]*color_vecs.getEvectors()[k] ;
            }
            multi1d<Double> norm2 =
	      sumMulti(localNorm2(color_vecs.getEvectors()[i]),phases.getSet());
            for(int t(0);t<phases.numSubsets();t++)
              color_vecs.getEvectors()[i][phases.getSet()[t]] /= sqrt(norm2[t]) ;
	  }

	}
	
      }
      if(params.param.BlockOrthoNormal){
	//create the block dims array
	Set blocks ;
	multi1d<int> block_dims(Nd);
	int k(0);
	QDPIO::cout<<"Block Orthonormalizing with block size: ";
	for(int d(0); d < Nd; d++){
	  if (d == params.param.src.decay_dir )
	    block_dims[d] = 1 ;
	  else{
	    if(k<params.param.block.size()){
	      block_dims[d] = params.param.block[k] ;
	      k++ ;
	    }
	    else
	      block_dims[d] = Layout::lattSize()[d] ;
	  }
	  QDPIO::cout<<	  block_dims[d]<< " ";
	}
	QDPIO::cout<<endl ;
	blocks.make(BlockFunc(block_dims));
	for(int i(0);i<Nvecs;i++){
	  for(int k(0);k<i;k++){
	      multi1d<DComplex> cc =
                sumMulti(localInnerProduct(color_vecs.getEvectors()[k],
                                           color_vecs.getEvectors()[i]),
                         blocks);
	      for(int b(0);b<blocks.numSubsets();b++)
                color_vecs.getEvectors()[i][blocks[b]] -=
                  cc[b]*color_vecs.getEvectors()[k] ;
	  }
	  multi1d<Double> norm2 =
	    sumMulti(localNorm2(color_vecs.getEvectors()[i]),blocks);
	  for(int b(0);b<blocks.numSubsets();b++)
	    color_vecs.getEvectors()[i][blocks[b]] /= sqrt(norm2[b]) ;
	}

      }//BlockOrthoNormalize

      {
	multi1d< multi1d<Double> > source_corrs(color_vecs.getNumVectors());
	for(int m=0; m < source_corrs.size(); ++m)
	  source_corrs[m] = sumMulti(localNorm2(color_vecs.getEvectors()[m]), 
				     phases.getSet());
	push(xml_out, "Source_correlators");
	write(xml_out, "source_corrs", source_corrs);
	pop(xml_out);
	//compute the "eigenvalues"
	push(xml_out,"SmearingEvals");
	for(int i(0);i<Nvecs;i++){
	  LatticeColorVector Svec = color_vecs.getEvectors()[i] ;
	  (*Smearing)(Svec, u_smr);
	  multi1d<DComplex> cc = 
	    sumMulti(localInnerProduct(color_vecs.getEvectors()[i],
				       Svec),  phases.getSet());

	  for(int t(0);t<phases.numSubsets();t++)
	    color_vecs.getEvalues()[i].weights[t] = real(cc[t]);

	  push(xml_out,"Vector") ;
	  write(xml_out, "VecNo",i);
	  write(xml_out, "Evals", color_vecs.getEvalues()[i].weights);
	  pop(xml_out);
	}
	pop(xml_out);
      }
      
      
      
      swatch.stop();
      QDPIO::cout << name << ": time for colorvec contstruction = "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;      

      pop(xml_out);  // CreateColorVecs

      /**/
      // Write the meta-data for this operator
      {
	XMLBufferWriter file_xml;

	push(file_xml, "GridColorVecs");
	write(file_xml, "Nvecs", Nvecs); 
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	XMLBufferWriter record_xml;
	push(record_xml, "EgenInfo");
	for(int i(0);i<Nvecs;i++){
	  push(record_xml, "EigenPair");
	  write(record_xml, "EigenPairNumber", i); 
	  write(record_xml, "EigenValues", color_vecs.getEvalues()[i].weights); 
	  pop(record_xml);
	}
	pop(record_xml);
	
	// Write the propagator xml info
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).setRecordXML(record_xml);
      }
      /**/

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

} // namespace Chroma
