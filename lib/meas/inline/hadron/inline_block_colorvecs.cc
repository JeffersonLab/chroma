// $Id: inline_block_colorvecs.cc,v 3.4 2009-02-23 19:50:51 edwards Exp $
/*! \file
 * \brief make color vectors
 *
 * Propagator calculation on a colorvector
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_block_colorvecs.h"
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

#include "meas/inline/io/named_objmap.h"
#include "meas/smear/gaus_smear.h"

#include "util/ferm/key_val_db.h"

namespace Chroma 
{ 
  namespace InlineBlockColorVecsEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineBlockColorVecsEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
      read(inputtop, "smearing_matrix_file", input.smearing_matrix_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineBlockColorVecsEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);
      write(xml, "smearing_matrix_file", input.smearing_matrix_file);
      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineBlockColorVecsEnv::Params::Param_t::Sources_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "decay_dir", input.decay_dir);
      
      input.smr = readXMLGroup(inputtop, "Smearing", "wvf_kind");
      input.link_smear = readXMLGroup(inputtop, "LinkSmearing", "LinkSmearingType");
      
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineBlockColorVecsEnv::Params::Param_t::Sources_t& out)
    {
      push(xml, path);
      write(xml, "decay_dir", out.decay_dir);

      xml << out.smr.xml;
      xml << out.link_smear.xml;

      pop(xml);
    }



    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineBlockColorVecsEnv::Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);
      
      read(inputtop, "Source", input.src);
      read(inputtop, "block", input.block);
      input.OrthoNormal=false;
      if(inputtop.count("OrthoNormal")!=0){
	read(inputtop, "OrthoNormal", input.OrthoNormal);
      }
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineBlockColorVecsEnv::Params::Param_t& input)
    {
      push(xml, path);
      
      write(xml, "Source", input.src);
      write(xml, "OrthoNormal", input.OrthoNormal);
      write(xml, "block", input.block);
      
      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineBlockColorVecsEnv::Params& input)
    {
      InlineBlockColorVecsEnv::Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineBlockColorVecsEnv::Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlineBlockColorVecsEnv 


  

    //-------------------------------------------------------------------------
    //! reweighting matrix operator
  struct KeySmearingMatrix_t
  {
    int                t_slice;      /*!< time slice */
  };

  //! Meson operator
  struct ValSmearingMatrix_t
  {
    multi2d< multi2d<ComplexD> > mat;  /*!< the matrix */
  };


  //----------------------------------------------------------------------------
  //! Holds key and value as temporaries
  struct KeyValSmearingMatrix_t
  {
    SerialDBKey<KeySmearingMatrix_t>  key;
    SerialDBData<ValSmearingMatrix_t> val;
  };


  //----------------------------------------------------------------------------
  //! KeySmearingMatrix reader
  void read(BinaryReader& bin, KeySmearingMatrix_t& param)
  {
    read(bin, param.t_slice);
  }

  //! SmearingMatrix write
  void write(BinaryWriter& bin, const KeySmearingMatrix_t& param)
  {
    write(bin, param.t_slice);
  }

  //! SmearingMatrix reader
  void read(XMLReader& xml, const std::string& path, KeySmearingMatrix_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "t_slice", param.t_slice);
  }

  //! SmearingMatrix writer
  void write(XMLWriter& xml, const std::string& path, const KeySmearingMatrix_t& param)
  {
    push(xml, path);

    write(xml, "t_slice", param.t_slice);

    pop(xml);
  }


  //----------------------------------------------------------------------------
  //! SmearingMatrix reader
  void read(BinaryReader& bin, ValSmearingMatrix_t& param)
  {

    int n1; //number of blocks
    int n2;
    read(bin, n2);    // the size is always written, even if 0
    read(bin, n1);    // the size is always written, even if 0
    param.mat.resize(n2,n1);
      
    int nv1; 
    int nv2; // number of color vectors
    read(bin, nv2);    // the size is always written, even if 0  
    read(bin, nv1);    // the size is always written, even if 0   
    //loops over blocks
    for(int bi=0; bi < param.mat.size1(); ++bi)
      for(int bj=0; bj < param.mat.size2(); ++bj){
	//loops over color vectors
	param.mat[bj][bi].resize(nv2,nv1);
	for(int ci=0; ci < param.mat[bj][bi].size1(); ++ci)
	  for(int cj=0; cj < param.mat[bj][bi].size2(); ++cj)
	    read(bin, param.mat[bj][bi][cj][ci]);
      }

  }

  //! SmearingMatrix write
  void write(BinaryWriter& bin, const ValSmearingMatrix_t& param)
  {
    write(bin, param.mat.size2());    // always write the size
    write(bin, param.mat.size1());    // always write the size
      
    write(bin, param.mat(0,0).size2());    // always write the size
    write(bin, param.mat(0,0).size1());    // always write the size
      
    //loops over blocks                                                   
    for(int bi=0; bi < param.mat.size1(); ++bi)
      for(int bj=0; bj < param.mat.size2(); ++bj)
	//loops over color vectors
	for(int ci=0; ci < param.mat[bj][bi].size1(); ++ci)
	  for(int cj=0; cj < param.mat[bj][bi].size2(); ++cj)
	    write(bin, param.mat[bj][bi][cj][ci]);
      
  }
  


  namespace InlineBlockColorVecsEnv 
  {
    namespace
    {
      AbsInlineMeasurement* blockMeasurement(XMLReader& xml_in, 
					     const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }
      
    const std::string name = "BLOCK_COLORVECS";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, blockMeasurement);
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

	push(xml_out, "BlockColorVecs");
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

      push(xml_out, "BlockColorVecs");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Block color vectors" << endl;

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
      

      Handle< QuarkSmearing<LatticeColorVector> >  Smearing;
      try{
	std::istringstream  xml_l(params.param.src.smr.xml);
	XMLReader  smrtop(xml_l);
	QDPIO::cout << "ColorVector Smearing type = " <<params.param.src.smr.id;
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

      

      /**
       //
       // create the output files
       //
       try
       {
       TheNamedObjMap::Instance().create< SubsetVectors<LatticeColorVector> >(params.named_obj.blocked_colorvec_id);
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
       SubsetVectors<LatticeColorVector>& blocked_color_vecs=TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(params.named_obj.blocked_colorvec_id);

      **/
      //
      // Read in the source along with relevant information.
      // 
      XMLReader source_file_xml, source_record_xml;

      QDPIO::cout << "Snarf the source from a named buffer" << endl;
      try
      {
	TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(params.named_obj.colorvec_id);

	// Snarf the source info. This is will throw if the colorvec_id is not there
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getFileXML(source_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getRecordXML(source_record_xml);

	// Write out the source header
	write(xml_out, "Source_file_info", source_file_xml);
	write(xml_out, "Source_record_info", source_record_xml);
      }    
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error extracting source_header: " << e << endl;
	QDP_abort(1);
      }
      SubsetVectors<LatticeColorVector>& color_vecs = 
	TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(params.named_obj.colorvec_id);
      
      QDPIO::cout << "Source successfully read and parsed" << endl;

      // The code goes here
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, params.param.src.decay_dir);
      
      int Nvecs = color_vecs.evectorsSize();

      color_vecs.openUpdate();

      if(params.param.OrthoNormal)
      {
	for(int i(0);i<Nvecs;i++) {
	  LatticeColorVector e_i; color_vecs.lookup(i,e_i);
	  for(int k(0);k<i;k++) {
	    LatticeColorVector e_k; color_vecs.lookup(k,e_k);
	    multi1d<DComplex> cc = 
	      sumMulti(localInnerProduct(e_k,
					 e_i),
		       phases.getSet());
	    for(int t(0);t<phases.numSubsets();t++) {
	      e_i[phases.getSet()[t]] -= cc[t]*e_k;
	    }
	  }

	  multi1d<Double> norm2 = sumMulti(localNorm2(e_i),phases.getSet());

	  for(int t(0);t<phases.numSubsets();t++) {
	    e_i[phases.getSet()[t]] /= sqrt(norm2[t]);
	  }
	  color_vecs.update(i, e_i);
	}

      }

      //block the block dims array
      Set blocks;
      multi1d<int> block_dims(Nd);
      int k(0);
      QDPIO::cout<<"Block Orthonormalizing with block size: ";
      for(int d(0); d < Nd; d++)
      {
	if (d == params.param.src.decay_dir )
	  block_dims[d] = 1;
	else{
	  if(k<params.param.block.size()){
	    block_dims[d] = params.param.block[k];
	    k++;
	  }
	  else
	    block_dims[d] = Layout::lattSize()[d];
	}
	QDPIO::cout<<	  block_dims[d]<< " ";
      }
      QDPIO::cout<<endl;

      blocks.make(BlockFunc(block_dims));
      for(int i(0);i<Nvecs;i++) {
	LatticeColorVector e_i; color_vecs.lookup(i, e_i);

	for(int k(0);k<i;k++) {
	  LatticeColorVector e_k; color_vecs.lookup(k, e_k);

	  multi1d<DComplex> cc =
	    sumMulti(localInnerProduct(e_k,e_i),blocks);

	  for(int b(0);b<blocks.numSubsets();b++) {
	    e_i[blocks[b]] -= cc[b]*e_k;
	  }
	}

	multi1d<Double> norm2 = sumMulti(localNorm2(e_i),blocks);
	for(int b(0);b<blocks.numSubsets();b++) {
	  e_i[blocks[b]] /= sqrt(norm2[b]);
	}
	color_vecs.update(i,e_i);
      }
      color_vecs.openRead(); // Finished update mode

      //#define PRINT_SMEARING_MATRIX
#ifdef PRINT_SMEARING_MATRIX
      // DB storage
      BinaryFxStoreDB< SerialDBKey<KeySmearingMatrix_t>, SerialDBData<ValSmearingMatrix_t> > 
	qdp_db(params.named_obj.smearing_matrix_file, DB_CREATE, db_cachesize, db_pagesize);

      //compute the smearing matrix
      QDPIO::cout<<"Computing the Smearing Matrix"<<endl;

      Set blks;
      blks.make(BlockFunc(params.param.src.decay_dir, params.param.block));
      LatticeColorVector vec;
      LatticeColorVector vecd;
      multi1d<DComplex> cc;

      for(int t(0);t<phases.numSubsets();t++)
      {
	KeyValSmearingMatrix_t kv;
	kv.key.key().t_slice = t;
	int Nblocks = blks.numSubsets();
	
	kv.val.data().mat.resize(Nblocks,Nblocks);
	for(int b(0);b<blks.numSubsets();b++)
	  for(int bb(0);bb<blks.numSubsets();bb++)
	    kv.val.data()(b,bb).resize(Nvecs,Nvecs);
	
	for(int b(0);b<blks.numSubsets();b++) {

	  for(int i(0);i<Nvecs;i++) {

	    vec = zero;
	    LatticeColorVector tmpvec; color_vecs.lookup(i,tmpvec);
	    vec[blks[b]] = tmpvec;
	    (*Smearing)(vec, u_smr);
	    for(int bb(b);bb<blks.numSubsets();bb++)  {

	      for(int j(i);j<Nvecs;j++) {
		LatticeColorVector tmpvec2; color_vecs.lookup(j,tmpvec2);
		vecd            = zero;
		vecd[blks[bb]] = tmpvec2;
		cc = sumMulti(localInnerProduct(vecd,vec),  phases.getSet());
		for(int t(0);t<phases.numSubsets();t++){
		  kv.val.data()(bb,b)(j,i) = cc[t];
		  kv.val.data()(b,bb)(i,j) = conj(cc[t]);
		}
	      }// j
	    }//bb
	  }//i
	}//b
	qdp_db.insert(kv.key, kv.val);
      }// loop over t
#endif

      {
	multi1d< multi1d<Double> > source_corrs(color_vecs.getNumVectors());
	for(int m=0; m < source_corrs.size(); ++m) {
	  LatticeColorVector tmpvec; color_vecs.lookup(m,tmpvec);
	  source_corrs[m] = sumMulti(localNorm2(tmpvec), phases.getSet());
	}
	push(xml_out, "Source_correlators");
	write(xml_out, "source_corrs", source_corrs);
	pop(xml_out);
	//compute the "eigenvalues"
	push(xml_out,"SmearingEvals");
	for(int i(0);i<Nvecs;i++)
	{
	  LatticeColorVector Svec; color_vecs.lookup(i,Svec);
	  LatticeColorVector tmpvec=Svec;

	  (*Smearing)(Svec, u_smr);
	  multi1d<DComplex> cc = 
	    sumMulti(localInnerProduct(tmpvec,Svec),  phases.getSet());

	  for(int t(0);t<phases.numSubsets();t++) {
	    color_vecs.getEvalue(i).weights[t] = real(cc[t]);
	  }
	  push(xml_out,"Vector");
	  write(xml_out, "VecNo",i);
	  write(xml_out, "Evals", color_vecs.getEvalue(i).weights);
	  pop(xml_out);
	}
	pop(xml_out);
      }
      
      swatch.stop();
      QDPIO::cout << name << ": time for colorvec contstruction = "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;      

      pop(xml_out);  // BlockColorVecs

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
	  write(record_xml, "EigenValues", color_vecs.getEvalue(i).weights); 
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
