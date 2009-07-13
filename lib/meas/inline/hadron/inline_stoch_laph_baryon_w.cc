// $Id: inline_stoch_laph_baryon_w.cc,v 3.2 2009-07-13 03:59:49 jbulava Exp $
/*! \file
 * \brief Inline measurement of stochastic group baryon operator
 *
 */

#include "handle.h"
#include "meas/inline/hadron/inline_stoch_laph_baryon_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/displacement.h"
#include "util/ferm/diractodr.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include <sstream> 

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineStochLapHBaryonEnv 
  { 
    //! Number of quarks to be used in this construction
    const int N_quarks = 3;
	
	
    //
    // The spin basis matrix to goto Dirac
    //
    SpinMatrix rotate_mat(adj(DiracToDRMat()));

    // Reader for input parameters
    void read(XMLReader& xml, const string& path, InlineStochLapHBaryonEnv::Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

			int version;
			read(paramtop, "version", version);

					multi1d< multi1d<int> > temp;
			switch (version) 
			{
				case 1:
					
					break;
				
				default :

					QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
					QDP_abort(1);
			}

			read(paramtop, "BaryonOperators", param.bops);

			param.link_smearing         = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");

			read(paramtop, "NumNoises", param.nnoises);

		}


		// Writer for input parameters
		void write(XMLWriter& xml, const string& path, const InlineStochLapHBaryonEnv::Params::Param_t& param)
		{
			push(xml, path);

			int version = 1;

			write(xml, "version", version);
			write(xml, "BaryonOperators", param.bops);
			write(xml, "NumNoises", param.nnoises);
      xml << param.link_smearing.xml;

      pop(xml);
    }


    //! Read named objects 
    void read(XMLReader& xml, const string& path, InlineStochLapHBaryonEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "BaryonOutfile", input.baryon_file);
      read(inputtop, "QuarkFiles", input.quark_files);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "BaryonOutfile", input.baryon_file);
      write(xml, "QuarkFiles", input.quark_files);

      pop(xml);
    }
  }


  namespace InlineStochLapHBaryonEnv 
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

    const std::string name = "STOCH_LAPH_BARYON";

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
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Parameters for source construction
      write(xml_out, "Param", param);

      // Write out the output propagator/source configuration info
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }

    

    //--------------------------------------------------------------------------
    //Support for the diquarks

    void makeDiquark( multi1d<LatticeComplex> & diquark, const multi1d<LatticeComplex> & q0,
		      const multi1d<LatticeComplex> & q1, const Subset & subset )
    {


      //The signs for the diquark are taken from
      //the colorContract function in qdp_primcolorvec.h
      diquark[0][subset] =  q0[0]*q1[1] - q0[1]*q1[0];
      diquark[1][subset] =  q0[1]*q1[2] - q0[2]*q1[1];
      diquark[2][subset] =  q0[2]*q1[0] - q0[0]*q1[2];


    }


    void makeColorSinglet (LatticeComplex & singlet, const multi1d<LatticeComplex> & diquark, 
			   const multi1d<LatticeComplex> & q2, const Subset & subset)
    {

      singlet[subset] = diquark[0] * q2[2];
      singlet[subset] += diquark[1] * q2[0];  
      singlet[subset] += diquark[2] * q2[1];  
    }
 

    //--------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "stoch_laph_baryon");
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
	QDPIO::cerr << InlineStochGroupBaryonEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineStochGroupBaryonEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "StochLapHBaryon");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << InlineStochGroupBaryonEnv::name << 
				": Stochastic LapH-Diluted Baryon Operators" << endl;

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
      
      // Sanity check - if this doesn't work we have serious problems
      if (phases.numSubsets() != QDP::Layout::lattSize()[decay_dir])
      {
	QDPIO::cerr << name << ": number of time slices not equal to that in the decay direction: " 
		    << QDP::Layout::lattSize()[decay_dir]
		    << endl;
	QDP_abort(1);
      }
 
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
								       linktop, params.param.link_smearing.path));

	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception link smearing: " << e << endl;
	QDP_abort(1);
      }

      MesPlq(xml_out, "Smeared_Observables", u_smr);

      

			//Initialize the quark database manager	
			


      //
      // Baryon operators
      //

      				
      QDPIO::cout << "Num Ordering = " << num_orderings << endl;

			

      for(int t0 = 0; t0 < participating_timeslices.size() ; ++t0)
      {
	
				//Loop over orderings
				
				//Loop over elementals

				//Initialize the Slow Fourier Transform Object for this elemental 
			  int decay_dir = diluted_quarks[0]->getDecayDir();
	
				SftMom phases(params.param.moms, decay_dir);

				//Loop over moms 
	
				
	      QDPIO::cout << "Ordering = " << ord << endl;

	  
				//Loop over d1
				//Loop over d2
		{

		  keySmearedDispColorVector[0].dil = i;
		  keySmearedDispColorVector[1].dil = j;

		  //Form the di-quark to save on recalculating 
		  multi1d<LatticeComplex> diquark(Nc);

		  const multi1d<LatticeComplex> &q0 = smrd_disp_srcs.getDispSource(n0, 
										   keySmearedDispColorVector[0]); 

		  const multi1d<LatticeComplex> &q1 = smrd_disp_srcs.getDispSource(n1, 
										   keySmearedDispColorVector[1]);

		 
		  watch.reset();
		  watch.start();
		  //For the source, restrict this operation to a subset
		  makeDiquark( diquark, q0 , q1, phases.getSet()[ participating_timeslices[t0] ] ); 
		  watch.stop();

		  /*QDPIO::cout<< " Made diquark : time = " << 
		    watch.getTimeInSeconds() << "secs" << endl;
		  */		

		  for(int k = 0 ; k < diluted_quarks[n2]->getDilSize(t0) ; ++k)	
		  {

		    keySmearedDispColorVector[2].dil = k;

		    // Contract over color indices with antisym tensor.
		    // NOTE: the creation operator only lives on a time slice, so restrict
		    // the operation to that time slice

		    LatticeComplex c_oper;

		    const multi1d<LatticeComplex> &q2 = smrd_disp_srcs.getDispSource(n2, 
										     keySmearedDispColorVector[2]);

		    watch.reset();
		    watch.start();

		    makeColorSinglet( c_oper, diquark, q2, phases.getSet()[ 
					participating_timeslices[t0] ] );

		    watch.stop();

		    /*QDPIO::cout<< "Made Color singlet : time =  " <<  
		      watch.getTimeInSeconds() << "secs" << endl;
		    */	
		    
		    // Slow fourier-transform
		    // We can restrict what the FT routine requires to a subset.
		    watch.reset();
		    watch.start();


		    
 				multi2d<DComplex> c_sum;
				int num_mom;

					c_sum = phases.sft(c_oper, participating_timeslices[t0]);
					num_mom = phases.numMom();
			
		    watch.stop();
 

		  } // end for k
		} // end for j
	      } // end for i
	    }//end ord 

	    swiss.stop();


	    QDPIO::cout << "Source operator construction: operator= " << l 
			<< "  time= "
			<< swiss.getTimeInSeconds() 
			<< " secs" << endl;

	    QDPIO::cout << "Source operator testval(t0 = " << 
	      participating_timeslices[t0] << ") = " << 
	      creat_oper.time_slices[0].orderings[0].dilutions(0,0,0).mom_projs[0].op[0];



	    //Hard code the elemental op name for now 
	    std::stringstream cnvrt;
	    cnvrt <<  creat_oper.id  << "_t" << participating_timeslices[t0] << "_src.lime";

	    std::string filename;

	    filename = cnvrt.str(); 

	    // Write the meta-data and the binary for this operator
	    swiss.reset();
	    swiss.start();
	    {
	      XMLBufferWriter     src_record_xml, file_xml;
	      BinaryBufferWriter  src_record_bin;

	      push(file_xml, "SourceBaryonOperator");
	      write(file_xml, "Params", params.param);
	      write(file_xml, "Config_info", gauge_xml);
	      write(file_xml, "Op_Info",qqq_oplist.ops[l]);

	      push(file_xml, "QuarkSources");

	      push(file_xml, "Quark_l");
	      push(file_xml, "TimeSlice");
	      push(file_xml, "Dilutions");
	      for (int dil = 0; dil < diluted_quarks[0]->getDilSize(t0) ; ++dil)
	      {
		write( file_xml, "elem", 
		       diluted_quarks[0]->getSourceHeader(t0, dil) );

		//	QDPIO::cout<< "t0 = " << t0 << " dil = "<< dil <<
		//	" srdhdr = XX"<<diluted_quarks[0]->getSourceHeader(t0,dil) << endl;
	      }
	      pop(file_xml); //dilutions 
	      pop(file_xml); //TimeSlice
	      pop(file_xml); //Quark_l

	      push(file_xml, "Quark_m");
	      push(file_xml, "TimeSlice");
	      push(file_xml, "Dilutions");
	      for (int dil = 0; dil < diluted_quarks[1]->getDilSize(t0) ; ++dil)
	      {
		write( file_xml, "elem", 
		       diluted_quarks[1]->getSourceHeader(t0, dil) );
	      }
	      pop(file_xml); //dilutions 
	      pop(file_xml); //TimeSlice
	      pop(file_xml); //Quark_m

	      push(file_xml, "Quark_r");
	      push(file_xml, "TimeSlice");
	      push(file_xml, "Dilutions");
	      for (int dil = 0; dil < diluted_quarks[2]->getDilSize(t0) ; ++dil)
	      {
		write( file_xml, "elem", 
		       diluted_quarks[2]->getSourceHeader(t0, dil) );
	      }
	      pop(file_xml); //dilutions 
	      pop(file_xml); //TimeSlice
	      pop(file_xml); //Quark_r

	      pop(file_xml);//QuarkSources
	      push(file_xml, "QuarkSinks");

	      push(file_xml, "Quark_l");
	      write(file_xml, "PropHeader", diluted_quarks[0]->getPropHeader(0,0) );
	      pop(file_xml);

	      push(file_xml, "Quark_m");
	      write(file_xml, "PropHeader", diluted_quarks[1]->getPropHeader(0,0) );
	      pop(file_xml);

	      push(file_xml, "Quark_r");
	      write(file_xml, "PropHeader", diluted_quarks[2]->getPropHeader(0,0) );
	      pop(file_xml);

	      pop(file_xml); //Quark Sinks  
	      pop(file_xml);//SourceBaryonOp

	      QDPFileWriter qdp_file(file_xml, filename,     // are there one or two files???
				     QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);


	      write(src_record_xml, "BaryonCreationOperator", creat_oper);
	      write(src_record_bin, creat_oper);

	      write(qdp_file, src_record_xml, src_record_bin);

	    }
	    swiss.stop();

	    QDPIO::cout << "Source Operator writing: operator = " << 
	      l	<< "  time= "
			<< swiss.getTimeInSeconds() 
			<< " secs" << endl;

	    pop(xml_out); // BaryonOperator 

	  } // end for l (operator )

	} //End Make creation operator



	//Make Annilation Operator
	{

	  // The object holding the smeared and displaced color vector maps  
	  SmearedDispObjects smrd_disp_snks(params.param.displacement_length,
					    diluted_quarks, quarkSmearing, u_smr );


	  // Annihilation operator
	  BaryonOperator_t  annih_oper;
	  annih_oper.mom2_max    = 0;
	  annih_oper.decay_dir   = decay_dir;
	  annih_oper.seed_l      = diluted_quarks[0]->getSeed();
	  annih_oper.seed_m      = diluted_quarks[1]->getSeed();
	  annih_oper.seed_r      = diluted_quarks[2]->getSeed();
	  annih_oper.dilution_l  = params.param.quark_dils[0];
	  annih_oper.dilution_m  = params.param.quark_dils[1];
	  annih_oper.dilution_r  = params.param.quark_dils[2];
	  annih_oper.smearing    = params.param.quark_smearing;
	  annih_oper.perms       = perms;
	  annih_oper.time_slices.resize( 1 );

	  // Construct annihilation operator
	  QDPIO::cout << "Building Sink operators" << endl;

	  // Loop over each operator 
	  for(int l=0; l < qqq_oplist.ops.size(); ++l)
	  {
	    QDPIO::cout << "Elemental operator: op = " << l << endl;

	    annih_oper.id = qqq_oplist.ops[l].name;

	    // Loop over all orderings and build the operator
	    swiss.reset();
	    swiss.start();

	    // The keys for the spin and displacements for this particular elemental operator
	    multi1d<KeySmearedDispColorVector_t> keySmearedDispColorVector(N_quarks);

	    for(int n = 0 ; n < N_quarks ; ++n)
	    {
	      keySmearedDispColorVector[n].displacement = qqq_oplist.ops[l].quark[n].displacement;
	      keySmearedDispColorVector[n].spin         = qqq_oplist.ops[l].quark[n].spin;
	    }

	    annih_oper.time_slices[0].t0 = participating_timeslices[t0];
	    //annih_oper.time_slices[0].orderings.resize(num_orderings);
	    annih_oper.time_slices[0].orderings.resize(1);


			int ord = 0;
	    //for(int ord = 0 ; ord < num_orderings ; ++ord)
	    {
	      QDPIO::cout << "Ordering = " << ord << endl;

	      annih_oper.time_slices[0].orderings[ord].perm = perms[ord];

	      const int n0 = perms[ord][0];
	      const int n1 = perms[ord][1];
	      const int n2 = perms[ord][2];

	      // The operator must hold all the dilutions
	      // We know that all time slices match. However, not all time slices of the
	      // lattice maybe used

	      // Creation operator
	      BaryonOperator_t::TimeSlices_t::Orderings_t& aop = annih_oper.time_slices[0].orderings[ord];

	      aop.dilutions.resize(diluted_quarks[n0]->getDilSize(t0), diluted_quarks[n1]->getDilSize(t0),
				   diluted_quarks[n2]->getDilSize(t0) );

	      for (int n = 0 ; n < N_quarks ; ++n)
	      {
		keySmearedDispColorVector[n].t0 = t0;
	      }

	      for(int i = 0 ; i <  diluted_quarks[n0]->getDilSize(t0) ; ++i)
	      {
		for(int j = 0 ; j < diluted_quarks[n1]->getDilSize(t0) ; ++j)	      
		{

		  keySmearedDispColorVector[0].dil = i;
		  keySmearedDispColorVector[1].dil = j;

		  //Form the di-quark to save on recalculating 
		  multi1d<LatticeComplex> diquark(Nc);

		  const multi1d<LatticeComplex> &q0 = smrd_disp_snks.getDispSolution(n0, 
										     keySmearedDispColorVector[0]); 

		  const multi1d<LatticeComplex> &q1 = smrd_disp_snks.getDispSolution(n1, 
										     keySmearedDispColorVector[1]);


		  //QDPIO::cout<<"q0[0] testval= "<< peekSite(q0[0], orig)
		  //	<< endl; 

		  //QDPIO::cout<<"q1[0] testval= "<< peekSite(q1[0], orig)
		  //	<< endl; 


		  watch.reset();
		  watch.start();

		  makeDiquark( diquark, q0 , q1, all ); 

		  watch.stop();
		  /*QDPIO::cout << "Made diquark: time = " << 
		    watch.getTimeInSeconds() << "secs " << endl;
		  */

		  for(int k = 0 ; k < diluted_quarks[n2]->getDilSize(t0) ; ++k)	
		  {

		    keySmearedDispColorVector[2].dil = k;

		    // Contract over color indices with antisym tensor.
		    // There is a potential optimization here - the colorcontract of
		    // the first two quarks could be pulled outside the innermost dilution
		    // loop.
		    // NOTE: the creation operator only lives on a time slice, so restrict
		    // the operation to that time slice

		    LatticeComplex a_oper;

		    const multi1d<LatticeComplex> &q2 = smrd_disp_snks.getDispSolution(n2, 
										       keySmearedDispColorVector[2]);

		    //QDPIO::cout<<"q2[0] testval= "<< peekSite(q2[0], orig)
		    //<< endl;

		    watch.reset();
		    watch.start();

		    makeColorSinglet( a_oper, diquark, q2, all);

		    watch.stop();

		    /*
		      QDPIO::cout <<	"Made Color Singlet: time = " <<
		      watch.getTimeInSeconds() << "secs" << endl;
		    */
		    /*QDPIO::cout << "testval = " << peekSite(a_oper, orig) 
		      << endl;
		    */

		    watch.reset();
		    watch.start();

		    // Slow fourier-transform
		    multi2d<DComplex> a_sum;
				int num_mom;

					a_sum = phases.sft(
							a_oper);
					num_mom = phases.numMom();

		    watch.stop();
		    /*
		      QDPIO::cout << "Spatial Sums completed: time " << 
		      watch.getTimeInSeconds() << "secs" << endl;
		    */		
		    // Unpack into separate momentum and correlator
		    aop.dilutions(i,j,k).mom_projs.resize(num_mom);

		    for(int mom_num = 0 ; mom_num < num_mom ; ++mom_num) 
		    {
		      aop.dilutions(i,j,k).mom_projs[mom_num].mom = params.param.moms[mom_num];

		      aop.dilutions(i,j,k).mom_projs[mom_num].op = a_sum[mom_num];

		    }

		  } // end for k
		} // end for j
	      } // end for i
	    }//end ord 
	    swiss.stop();


	    QDPIO::cout << "Sink operator construction: operator= " << l 
			<< "  time= "
			<< swiss.getTimeInSeconds() 
			<< " secs" << endl;

	    QDPIO::cout << "Sink op testval( t0 = " << 
	      participating_timeslices[t0] << ") = " << 
	      annih_oper.time_slices[0].orderings[0].dilutions(0,0,0).mom_projs[0].op[0] 
			<< endl;

	    //Hard code the elemental op name for now 
	    std::stringstream cnvrt;
	    cnvrt <<  annih_oper.id  << "_t" << participating_timeslices[t0] << "_snk.lime";

	    std::string filename;

	    filename = cnvrt.str(); 

	    // Write the meta-data and the binary for this operator
	    swiss.reset();
	    swiss.start();
	    {
	      XMLBufferWriter     src_record_xml, file_xml;
	      BinaryBufferWriter  src_record_bin;

	      push(file_xml, "SinkBaryonOperator");
	      write(file_xml, "Params", params.param);
	      write(file_xml, "Config_info", gauge_xml);
	      write(file_xml, "Op_Info",qqq_oplist.ops[l]);
	      push(file_xml, "QuarkSources");

	      push(file_xml, "Quark_l");
	      push(file_xml, "TimeSlice");
	      push(file_xml, "Dilutions");
	      for (int dil = 0; dil < diluted_quarks[0]->getDilSize(t0) ; ++dil)
	      {
		write( file_xml, "elem", 
		       diluted_quarks[0]->getSourceHeader(t0, dil) );
	      }
	      pop(file_xml); //dilutions 
	      pop(file_xml); //TimeSlice
	      pop(file_xml); //Quark_l

	      push(file_xml, "Quark_m");
	      push(file_xml, "TimeSlice");
	      push(file_xml, "Dilutions");
	      for (int dil = 0; dil < diluted_quarks[1]->getDilSize(t0) ; ++dil)
	      {
		write( file_xml, "elem", 
		       diluted_quarks[1]->getSourceHeader(t0, dil) );
	      }
	      pop(file_xml); //dilutions 
	      pop(file_xml); //TimeSlice
	      pop(file_xml); //Quark_m

	      push(file_xml, "Quark_r");
	      push(file_xml, "TimeSlice");
	      push(file_xml, "Dilutions");
	      for (int dil = 0; dil < diluted_quarks[2]->getDilSize(t0) ; ++dil)
	      {
		write( file_xml, "elem", 
		       diluted_quarks[2]->getSourceHeader(t0, dil) );
	      }
	      pop(file_xml); //dilutions 
	      pop(file_xml); //TimeSlice
	      pop(file_xml); //Quark_r

	      pop(file_xml);//QuarkSources
	      push(file_xml, "QuarkSinks");

	      push(file_xml, "Quark_l");
	      write(file_xml, "PropHeader", diluted_quarks[0]->getPropHeader(0,0) );
	      pop(file_xml);

	      push(file_xml, "Quark_m");
	      write(file_xml, "PropHeader", diluted_quarks[1]->getPropHeader(0,0) );
	      pop(file_xml);

	      push(file_xml, "Quark_r");
	      write(file_xml, "PropHeader", diluted_quarks[2]->getPropHeader(0,0) );
	      pop(file_xml);

	      pop(file_xml);//QuarkSinks 
	      pop(file_xml);//SinkBaryonOperator

	      QDPFileWriter qdp_file(file_xml, filename,     // are there one or two files???
				     QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);


	      write(src_record_xml, "BaryonAnnihilationOperator", annih_oper);
	      write(src_record_bin, annih_oper);

	      write(qdp_file, src_record_xml, src_record_bin);

	    }
	    swiss.stop();

	    QDPIO::cout << "Sink Operator writing: operator = " << l
			<< "  time= " << swiss.getTimeInSeconds() << " secs" << endl;

	  } // end for l (operator )

	} //End Make annihilation operator

      } //end t0 
      // Close the namelist output file XMLDAT
      pop(xml_out);     // StochBaryon

      snoop.stop();
      QDPIO::cout << InlineStochGroupBaryonEnv::name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << InlineStochGroupBaryonEnv::name << ": ran successfully" << endl;

      END_CODE();
    } // func

  } // namespace InlineStochGroupBaryonEnv

  /*! @} */  // end of group hadron

} // namespace Chroma
